/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_to_block.h"
#include "pdm_part_to_block_priv.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_hilbert.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  pdm_part_to_block_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int pdm_part_to_block_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*
 * Static can cause pb if we're call function in multiple contexte.
 *   For example python and other C++ program
 * No static : truly global
 *  https://stackoverflow.com/questions/1856599/when-to-use-static-keyword-before-global-variables
 */
double t_elaps[3] = {0., 0., 0.};
double t_cpu[3] = {0., 0., 0.};
PDM_timer_t *t_timer[3] = {NULL, NULL, NULL};

int min_exch_rank[2] = {INT_MAX, INT_MAX};
int max_exch_rank[2] = {-1, -1};

unsigned long long exch_data[2] = {0, 0};

int n_ptb = 0;


/*=============================================================================
 * Static function definitions
 *============================================================================*/

void
PDM_extents_conformize(int    dim,
                       double extents[],
                       double eps)
{
  double max_range = 1.e-12;
  for (int i = 0; i < dim; i++) {
    max_range = PDM_MAX (max_range, extents[i+dim] - extents[i]);
  }

  // eps = 1.e-3
  const double epsilon = eps * max_range; // Add eps in cas of only one point ...
  for (int i = 0; i < dim; i++) {
    extents[i    ] -= 1.1 * epsilon; // On casse la symetrie !
    extents[i+dim] +=       epsilon;
  }
}

static
void
_counting_sort_long
(
 PDM_part_to_block_t *ptb
)
{
  // On a un pb si la distrib est géant
  int n_elt_block_tot = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
  int* block_n = (int *) malloc(n_elt_block_tot * sizeof(int));

  PDM_g_num_t* block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = ptb->sorted_recv_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    block_n[lid] += 1;
    block_gnum[i] = ptb->sorted_recv_gnum[i];
  }

  int* block_idx = (int *) malloc( (n_elt_block_tot+1) * sizeof(int));
  block_idx[0] = 0;
  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_idx[i+1] = block_idx[i] + block_n[i];
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = block_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    int idx_write = block_idx[lid] + block_n[lid]++;
    ptb->sorted_recv_gnum[idx_write] = block_gnum[i];
    ptb->order[idx_write] = i;
  }

  free(block_n);
  free(block_idx);
  free(block_gnum);

}

/**
 *
 * \brief  Define active ranks
 *
 * \param [inout]   ptb          Part to block structure
 *
 */

static void
_active_ranks
(
 PDM_part_to_block_t  *ptb
)
{

  assert (ptb->active_ranks == NULL);

  if (ptb->s_comm == 1) {
    ptb->is_my_rank_active = 1;
    ptb->n_active_ranks = 1;
    ptb->active_ranks = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
    ptb->active_ranks[0] = ptb->i_rank;
  }

  else {

    switch (ptb->t_distrib) {

    case PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC : {
      ptb->is_my_rank_active = 1;
      ptb->n_active_ranks = ptb->s_comm;
      ptb->active_ranks   = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
      for (int i = 0; i < ptb->n_active_ranks; i++) {
        ptb->active_ranks[i] = i;
      }
      break;
    }

    case PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE :
    case PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE : {

      ptb->is_my_rank_active = 0;

      int rank_in_node = PDM_io_mpi_node_rank(ptb->comm);
      if (rank_in_node == 0) {
        ptb->is_my_rank_active = 1;
      }

      int *tag_active_ranks = (int *) malloc(sizeof(int) * ptb->s_comm);

      PDM_MPI_Allgather((void *) &ptb->is_my_rank_active, 1, PDM_MPI_INT,
                        (void *) tag_active_ranks, 1, PDM_MPI_INT,
                        ptb->comm);

      ptb->n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->n_active_ranks += 1;
        }
      }

      ptb->active_ranks   = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
      int n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->active_ranks[n_active_ranks++] = i;
        }
      }

      if (ptb->t_distrib == PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE) {
        break;
      }

      int n_node = ptb->n_active_ranks;
      double part_active_node = PDM_MIN (1, ptb->part_active_node);
      part_active_node = PDM_MAX (0, part_active_node);
      ptb->n_active_ranks = (int) floor (n_node * part_active_node);
      ptb->n_active_ranks = PDM_MAX (1, ptb->n_active_ranks);

      n_active_ranks = 0;
      int coeff = n_node / ptb->n_active_ranks;
      for (int i = 0; i < n_node; i++) {
        if (i % coeff == 0) {
          ptb->active_ranks[n_active_ranks++] = ptb->active_ranks[i];
        }
        if (n_active_ranks == ptb->n_active_ranks) {
          break;
        }
      }

      assert (n_active_ranks == ptb->n_active_ranks);

      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0,"Error cs_part_to_bloc : unknown distribute type\n");
      abort();
    }
    }

    /* Dump */

    if (0 == 1) {
      PDM_printf("active ranks : ");
      for(int i = 0; i < ptb->n_active_ranks; i++)
        PDM_printf("%i ", ptb->active_ranks[i]);
      PDM_printf("\n");
    }
  }
}


/**
 *
 * \brief Distrib data
 *
 * \param [inout] ptb              Part to block structure
 * \param [in]    n_totalData      Total number of data
 * \param [out]   data_distrib_index Element global number
 * \param [out]   s_block_min       Local number of elements
 * \param [out]   s_block_max       Number of partition
 *
 */

static void
_distrib_data
(
 PDM_part_to_block_t *ptb,
 int                  user_distrib
)
{
  PDM_g_num_t _id_max     = 0;
  PDM_g_num_t _id_max_max = 0;

  ptb->n_elt_proc= 0;
  for (int i = 0; i < ptb->n_part; i++) {
    ptb->n_elt_proc+= ptb->n_elt[i];
    if(user_distrib == 0) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        PDM_g_num_t gnum = PDM_ABS(ptb->gnum_elt[i][j]);
        _id_max = PDM_MAX (_id_max, gnum);
      }
    }
  }



  if(user_distrib == 0) {

    PDM_MPI_Allreduce (&_id_max,
                       &_id_max_max,
                       1,
                       PDM__PDM_MPI_G_NUM,
                       PDM_MPI_MAX,
                       ptb->comm);

    if (ptb->weight == NULL) {
      PDM_g_num_t _n_rank_data = _id_max_max / ptb->n_active_ranks;
      PDM_g_num_t _rest = _id_max_max % ptb->n_active_ranks;

      int n_rank_data = (int) (_n_rank_data);
      int rest       = (int) (_rest);

      ptb->s_block_max = n_rank_data;
      ptb->s_block_min = n_rank_data;

      if (rest != 0) {
        ptb->s_block_max += 1;
      }

      for (int i = 0; i < ptb->s_comm + 1; i++) {
        ptb->data_distrib_index[i] = 0;
      }

      int k = 0;
      int idx = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        ptb->data_distrib_index[i+1] +=  ptb->data_distrib_index[i];
        if (idx < ptb->n_active_ranks) {
          if (ptb->active_ranks[idx] == i) {
            ptb->data_distrib_index[i+1] += n_rank_data;
            if (k < rest)
              ptb->data_distrib_index[i+1] += 1;
            k += 1;
            idx++;
          }

        }
      }

    }

    else {
      const int dim = 2;
      const int  n_active_ranks = ptb->n_active_ranks;
      const int  sampling_factor = _sampling_factors[dim];

      double **weight = ptb->weight;
      int *_active_ranks = ptb->active_ranks;

      PDM_g_num_t *rank_index = NULL;
      PDM_distrib_weight(             sampling_factor,
                                      n_active_ranks,
                                      ptb->n_part,
                                      ptb->n_elt,
              (const PDM_g_num_t **)  ptb->gnum_elt,
              (const double      **)  weight,
                                      pdm_part_to_block_distrib_n_iter_max,
                                      pdm_part_to_block_distrib_tol,
                                      ptb->comm,
                                      &rank_index);

      if(0 == 1) {
        PDM_log_trace_array_long(rank_index, n_active_ranks, "rank_index :: ");
      }

      ptb->data_distrib_index[0] = 0;
      int k = 0;
      for (int i = 0; i < n_active_ranks; i++) {
        int i_activeRank = _active_ranks[i];
        while (k < i_activeRank) {
          ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
          k++;
        }
        ptb->data_distrib_index[k+1] = rank_index[i+1];
        k++;
      }
      while (k < ptb->s_comm) {
        ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
        k++;
      }

      free (rank_index);

    }
  } // If User

  /* Affichage */
  if (1 == 0) {
    if (ptb->i_rank == 0) {
      PDM_printf("data_distrib_index : ");
      for(int i = 0; i < ptb->s_comm + 1; i++)
        PDM_printf(PDM_FMT_G_NUM" ", ptb->data_distrib_index[i]);
      PDM_printf("\n");
    }
  }

  ptb->n_send_data = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_data = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Pour chaque donnee le proc ou elle va etre envoyee */

  ptb->dest_proc = (int *) malloc (sizeof(int) * ptb->n_elt_proc);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    for (int j = 0; j < ptb->n_elt[i]; j++) {

      PDM_g_num_t _gnum_elt = PDM_ABS(ptb->gnum_elt[i][j]) - 1;

      int iproc = PDM_binary_search_gap_long (_gnum_elt,
                                              ptb->data_distrib_index,
                                              ptb->s_comm + 1);

      ptb->dest_proc[++idx] = iproc;
      // assert (ptb->dest_proc[idx] >= 0);
      ptb->n_send_data[iproc] += 1;
    }
  }

  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  ptb->i_send_data = (int *) malloc(sizeof(int) * ptb->s_comm);
  ptb->i_recv_data = (int *) malloc(sizeof(int) * ptb->s_comm);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];

  PDM_g_num_t *send_gnum =
    (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * ptb->tn_send_data) ;

  for (int i = 0; i < ptb->s_comm; i++)
    ptb->n_send_data[i] = 0;

  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[idx];
      send_gnum[ptb->i_send_data[iproc] +
                ptb->n_send_data[iproc]] = PDM_ABS(ptb->gnum_elt[i][j]);
      idx++;
      ptb->n_send_data[iproc] += 1;
    }
  }

  ptb->sorted_recv_gnum =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * ptb->tn_recv_data);

  PDM_MPI_Alltoallv(send_gnum,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->sorted_recv_gnum,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->comm);

  free(send_gnum);

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->order = malloc (sizeof(int) * ptb->tn_recv_data);

    char *env_var = NULL;
    int counting_sort = 0;
    env_var = getenv ("PDM_PTB_COUNTING_SORT");
    if (env_var != NULL) {
      counting_sort = atoi(env_var);
    }

    if(counting_sort == 0) {
      for (int i = 0; i < ptb->tn_recv_data; i++) {
        ptb->order[i] = i;
      }

      PDM_sort_long (ptb->sorted_recv_gnum,
                     ptb->order,
                     ptb->tn_recv_data);
    } else {
      if(ptb->i_rank == 0) {
        printf("_counting_sort_long \n");
      }
      _counting_sort_long(ptb);
    }

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }

  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->n_elt_block = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
    // printf("ptb->n_elt_block::%d \n", ptb->n_elt_block);
    // printf(" ptb->data_distrib_index[ptb->i_rank+1]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank+1]);
    // printf(" ptb->data_distrib_index[ptb->i_rank]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank]);
  }

  /*
   * Cleanup
   */

  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    ptb->block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
      else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    ptb->block_gnum = realloc (ptb->block_gnum, sizeof(PDM_g_num_t) * ptb->n_elt_block);

    ptb->block_gnum_count = malloc (sizeof(int) * ptb->n_elt_block);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  }

  else {

    ptb->block_gnum = ptb->sorted_recv_gnum;

    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  /*
   * Create idx_partial for reverse
   */
  ptb->enable_reverse = 1;
  if ( (ptb->t_post         != PDM_PART_TO_BLOCK_POST_NOTHING)       &&
       (ptb->t_post         != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) &&
       (ptb->enable_reverse == 1)) {
    ptb->idx_partial = malloc(ptb->tn_recv_data * sizeof(int));
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      PDM_g_num_t gnum = ptb->sorted_recv_gnum[i];
      int idx_in_partial_block = PDM_binary_search_long(gnum, ptb->block_gnum, ptb->n_elt_block);
      assert(idx_in_partial_block != -1);
      ptb->idx_partial[i] = idx_in_partial_block;
    }
  }
}

static
void
_distrib_data_hilbert
(
  PDM_part_to_block_t   *ptb,
  double               **pvtx_coords,
  double               **weight
)
{
  ptb->n_elt_proc = 0;
  int *part_idx = malloc( (ptb->n_part + 1) * sizeof(int));
  part_idx[0] = 0;
  for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
    ptb->n_elt_proc += ptb->n_elt[i_part];
    part_idx[i_part+1] = part_idx[i_part] + ptb->n_elt[i_part];
  }

  double *concat_vtx_coord = NULL;
  double *concat_weight    = NULL;
  if(ptb->n_part == 1 ) {
    concat_vtx_coord = pvtx_coords[0];
    concat_weight    = weight     [0];
  } else if (ptb->n_part > 1) {
    concat_vtx_coord = malloc( 3 * ptb->n_elt_proc * sizeof(double));
    concat_weight    = malloc(     ptb->n_elt_proc * sizeof(double));

    int shift = 0;
    for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        int idx_write = 3*(shift + i_elt);
        concat_vtx_coord[idx_write    ] = pvtx_coords[i_part][3 * i_elt    ];
        concat_vtx_coord[idx_write + 1] = pvtx_coords[i_part][3 * i_elt + 1];
        concat_vtx_coord[idx_write + 2] = pvtx_coords[i_part][3 * i_elt + 2];
      }
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        concat_weight[shift+i_elt] = weight[i_part][i_elt];
      }
      shift += ptb->n_elt[i_part];
    }
  }

  /** Initialisation **/
  int dim = 3;
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_extents_conformize(dim, extents, 1e-3);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (ptb->n_elt_proc * sizeof(PDM_hilbert_code_t));

  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, ptb->n_elt_proc, concat_vtx_coord, hilbert_codes);

  if(ptb->n_part > 1 ) {
    free(concat_vtx_coord);
  }

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((ptb->s_comm+1) * sizeof(PDM_hilbert_code_t));

  PDM_hilbert_build_rank_index(3,
                               ptb->s_comm,  // Number of chunk
                               ptb->n_elt_proc,
                               hilbert_codes,
                               concat_weight,
                               NULL,
                               hilbert_codes_idx, // Is the distrib
                               ptb->comm);


  if(0 == 1) {
    PDM_log_trace_array_double(hilbert_codes, ptb->n_elt_proc, "hilbert_codes :: ");
    PDM_log_trace_array_double(hilbert_codes_idx, ptb->s_comm+1, "hilbert_codes_idx :: ");
  }

  if(ptb->n_part > 1 ) {
    free(concat_weight);
  }

  ptb->n_send_data = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_data = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Pour chaque donnee le proc ou elle va etre envoyee */
  ptb->dest_proc = (int *) malloc (sizeof(int) * ptb->n_elt_proc);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */
  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  /*
   * Caution :  We need to send in the ordering of partition !!!!
   */
  int idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = hilbert_order[i_concat_elt]; // Donc dans la frame de depart

      size_t t_rank = PDM_hilbert_quantile_search(ptb->s_comm,
                                                  hilbert_codes[i_concat_elt],
                                                  hilbert_codes_idx);
      ptb->dest_proc  [idx++ ]  = t_rank;
      ptb->n_send_data[t_rank] += 1;
    }
  }

  free(hilbert_codes_idx);
  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  ptb->i_send_data = (int *) malloc(sizeof(int) * ptb->s_comm);
  ptb->i_recv_data = (int *) malloc(sizeof(int) * ptb->s_comm);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];

  PDM_g_num_t        *send_gnum  = (PDM_g_num_t        *) malloc (sizeof(PDM_g_num_t       ) * ptb->tn_send_data);
  PDM_hilbert_code_t *send_codes = (PDM_hilbert_code_t *) malloc (sizeof(PDM_hilbert_code_t) * ptb->tn_send_data);

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  // Dans le cas géométrique on doit envoyer le gnum + le codes
  idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {

      int t_rank         = ptb->dest_proc[idx++];
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = hilbert_order[i_concat_elt]; // Donc dans la frame de depart

      send_gnum [ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = PDM_ABS(ptb->gnum_elt[i_part][j]);
      send_codes[ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = hilbert_codes[i_concat_elt]; // Not sure
      ptb->n_send_data[t_rank] += 1;
    }
  }
  // free(hilbert_order);
  free(hilbert_codes);

  ptb->sorted_recv_gnum                 = (PDM_g_num_t        *) malloc(sizeof(PDM_g_num_t       ) * ptb->tn_recv_data);
  PDM_hilbert_code_t *sorted_recv_codes = (PDM_hilbert_code_t *) malloc(sizeof(PDM_hilbert_code_t) * ptb->tn_recv_data);

  PDM_MPI_Alltoallv(send_gnum,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->sorted_recv_gnum,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->comm);

  PDM_MPI_Alltoallv(send_codes,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    PDM_MPI_DOUBLE,
                    sorted_recv_codes,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    PDM_MPI_DOUBLE,
                    ptb->comm);

  free(send_gnum);
  free(send_codes);
  free(part_idx);

  if(0 == 1) {
    PDM_log_trace_array_long  (ptb->sorted_recv_gnum, ptb->tn_recv_data, "ptb->sorted_recv_gnum :: ");
    PDM_log_trace_array_double(sorted_recv_codes    , ptb->tn_recv_data, "sorted_recv_codes :: ");
  }

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->order = malloc (sizeof(int) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      ptb->order[i] = i;
    }

    PDM_hilbert_local_order(ptb->tn_recv_data, sorted_recv_codes, ptb->order);
    PDM_order_array(ptb->tn_recv_data, sizeof(PDM_g_num_t), ptb->order, ptb->sorted_recv_gnum);

    // PDM_log_trace_array_double(sorted_recv_codes    , ptb->tn_recv_data, "sorted_recv_codes :: ");

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }

  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) {
    PDM_error(__FILE__, __LINE__, 0,"Error PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM : not implemented \n");
  }

  /*
   * Cleanup
   */
  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    ptb->block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      } else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    ptb->block_gnum = realloc (ptb->block_gnum, sizeof(PDM_g_num_t) * ptb->n_elt_block);

    ptb->block_gnum_count = malloc (sizeof(int) * ptb->n_elt_block);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  } else {
    ptb->block_gnum       = ptb->sorted_recv_gnum;
    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  // Generate distribution
  PDM_distrib_compute(ptb->n_elt_block, ptb->data_distrib_index, -1, ptb->comm);

  free(sorted_recv_codes);

  /*
   * To do : ptb->enable_reverse = 1;
   */
  ptb->enable_reverse = 0;


}


static
void
_distrib_data_morton
(
  PDM_part_to_block_t   *ptb,
  double               **pvtx_coords,
  double               **weight
)
{
  ptb->n_elt_proc = 0;
  int *part_idx = malloc( (ptb->n_part + 1) * sizeof(int));
  part_idx[0] = 0;
  for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
    ptb->n_elt_proc += ptb->n_elt[i_part];
    part_idx[i_part+1] = part_idx[i_part] + ptb->n_elt[i_part];
  }

  double *concat_vtx_coord = NULL;
  double *concat_weight    = NULL;
  if(ptb->n_part == 1 ) {
    concat_vtx_coord = pvtx_coords[0];
    concat_weight    = weight     [0];
  } else {
    concat_vtx_coord = malloc( 3 * ptb->n_elt_proc * sizeof(double));
    concat_weight    = malloc(     ptb->n_elt_proc * sizeof(double));

    int shift = 0;
    for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        int idx_write = 3*(shift + i_elt);
        concat_vtx_coord[idx_write    ] = pvtx_coords[i_part][3 * i_elt    ];
        concat_vtx_coord[idx_write + 1] = pvtx_coords[i_part][3 * i_elt + 1];
        concat_vtx_coord[idx_write + 2] = pvtx_coords[i_part][3 * i_elt + 2];
      }
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        concat_weight[shift+i_elt] = weight[i_part][i_elt];
      }
      shift += ptb->n_elt[i_part];
    }
  }


  /** Initialisation **/
  int dim = 3;
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  // PDM_morton_get_coord_extents_par(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_morton_get_coord_extents(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_extents_conformize(dim, extents, 1e-3);

  /** morton Coordinates Computation **/
  PDM_morton_code_t *morton_codes     = (PDM_morton_code_t *) malloc (ptb->n_elt_proc * sizeof(PDM_morton_code_t));

  const PDM_morton_int_t max_level = PDM_morton_max_level;
  double d[3];
  double s[3];
  PDM_morton_encode_coords(dim, max_level, extents, ptb->n_elt_proc, concat_vtx_coord, morton_codes, d, s);

  if(ptb->n_part != 1) {
    free(concat_vtx_coord);
  }

  PDM_morton_code_t *morton_codes_idx = (PDM_morton_code_t *) malloc ((ptb->s_comm+1) * sizeof(PDM_morton_code_t));
  PDM_morton_build_rank_index(3,
                              max_level, // Number of chunk
                              ptb->n_elt_proc,
                              morton_codes,
                              concat_weight,
                              NULL,
                              morton_codes_idx, // Is the distrib
                              ptb->comm);

  if(ptb->n_part != 1) {
    free(concat_weight);
  }

  ptb->n_send_data = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_data = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Pour chaque donnee le proc ou elle va etre envoyee */
  ptb->dest_proc = (int *) malloc (sizeof(int) * ptb->n_elt_proc);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */
  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  /*
   * Caution :  We need to send in the ordering of partition !!!!
   */
  int idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = morton_order[i_concat_elt]; // Donc dans la frame de depart

      size_t t_rank = PDM_morton_quantile_search(ptb->s_comm,
                                                 morton_codes[i_concat_elt],
                                                 morton_codes_idx);
      ptb->dest_proc  [idx++ ]  = t_rank;
      ptb->n_send_data[t_rank] += 1;
    }
  }

  free(morton_codes_idx);

  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  ptb->i_send_data = (int *) malloc(sizeof(int) * ptb->s_comm);
  ptb->i_recv_data = (int *) malloc(sizeof(int) * ptb->s_comm);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];


  PDM_g_num_t        *send_gnum  = (PDM_g_num_t        *) malloc (sizeof(PDM_g_num_t       ) * ptb->tn_send_data);
  PDM_morton_code_t  *send_codes = (PDM_morton_code_t *) malloc (sizeof(PDM_morton_code_t) * ptb->tn_send_data);

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  // Dans le cas géométrique on doit envoyer le gnum + le codes
  idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {

      int t_rank         = ptb->dest_proc[idx++];
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = morton_order[i_concat_elt]; // Donc dans la frame de depart

      send_gnum [ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = PDM_ABS(ptb->gnum_elt[i_part][j]);
      PDM_morton_copy(morton_codes[i_concat_elt], &send_codes[ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]]);

      ptb->n_send_data[t_rank] += 1;
    }
  }
  // free(morton_order);
  free(morton_codes);

  ptb->sorted_recv_gnum                 = (PDM_g_num_t        *) malloc(sizeof(PDM_g_num_t       ) * ptb->tn_recv_data);
  PDM_morton_code_t *sorted_recv_codes = (PDM_morton_code_t *) malloc(sizeof(PDM_morton_code_t) * ptb->tn_recv_data);

  PDM_MPI_Alltoallv(send_gnum,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->sorted_recv_gnum,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->comm);

  PDM_MPI_Datatype mpi_morton_type;
  PDM_MPI_Type_create_contiguous(4, PDM_MPI_INT, &mpi_morton_type);
  PDM_MPI_Type_commit(&mpi_morton_type);
  PDM_MPI_Alltoallv(send_codes,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    mpi_morton_type,
                    sorted_recv_codes,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    mpi_morton_type,
                    ptb->comm);
  PDM_MPI_Type_free(&mpi_morton_type);

  free(send_gnum);
  free(send_codes);
  free(part_idx);


  if(0 == 1) {
    PDM_log_trace_array_long  (ptb->sorted_recv_gnum, ptb->tn_recv_data, "ptb->sorted_recv_gnum :: ");
  }

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->order = malloc (sizeof(int) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      ptb->order[i] = i;
    }

    PDM_morton_local_order(ptb->tn_recv_data, sorted_recv_codes, ptb->order);
    PDM_order_array(ptb->tn_recv_data, sizeof(PDM_g_num_t      ), ptb->order, ptb->sorted_recv_gnum);
    // PDM_order_array(ptb->tn_recv_data, sizeof(PDM_morton_code_t), ptb->order, ptb->sorted_recv_gnum);

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }
  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) {
    PDM_error(__FILE__, __LINE__, 0,"Error PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM : not implemented \n");
  }

  /*
   * Cleanup
   */
  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    ptb->block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      } else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    ptb->block_gnum = realloc (ptb->block_gnum, sizeof(PDM_g_num_t) * ptb->n_elt_block);

    ptb->block_gnum_count = malloc (sizeof(int) * ptb->n_elt_block);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  } else {
    ptb->block_gnum       = ptb->sorted_recv_gnum;
    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  // Generate distribution
  PDM_distrib_compute(ptb->n_elt_block, ptb->data_distrib_index, -1, ptb->comm);

  free(sorted_recv_codes);

  /*
   * To do : ptb->enable_reverse = 1;
   */
  ptb->enable_reverse = 0;


}


static void
_compute_global_weights
(
 PDM_part_to_block_t *ptb
 )
{
  /* Send local weights */
  int *send_count = PDM_array_zeros_int (ptb->s_comm);
  double *part_weight = malloc (sizeof(double) * ptb->tn_send_data);
  int idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      part_weight[ptb->i_send_data[rank] + send_count[rank]++] = ptb->weight[i][j];
    }
  }

  double *recv_weight = malloc (sizeof(double) * ptb->tn_recv_data);
  PDM_MPI_Alltoallv (part_weight,
                     ptb->n_send_data,
                     ptb->i_send_data,
                     PDM_MPI_DOUBLE,
                     recv_weight,
                     ptb->n_recv_data,
                     ptb->i_recv_data,
                     PDM_MPI_DOUBLE,
                     ptb->comm);

  /* Sum received weights */
  double *block_weight = malloc (sizeof(double) * ptb->n_elt_block);
  for (int i = 0; i < ptb->n_elt_block; i++) {
    block_weight[i] = 0.;
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    block_weight[idx] += recv_weight[j];
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    recv_weight[j] = block_weight[idx];
  }
  free (block_weight);

  /* Send back global weights */
  PDM_MPI_Alltoallv (recv_weight,
                     ptb->n_recv_data,
                     ptb->i_recv_data,
                     PDM_MPI_DOUBLE,
                     part_weight,
                     ptb->n_send_data,
                     ptb->i_send_data,
                     PDM_MPI_DOUBLE,
                     ptb->comm);
  free (recv_weight);

  /* Store global weights */
  ptb->weight_g = malloc (sizeof(double *) * ptb->n_part);
  PDM_array_reset_int (send_count, ptb->s_comm, 0);
  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    ptb->weight_g[i] = malloc (sizeof(double) * ptb->n_elt[i]);
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      ptb->weight_g[i][j] = part_weight[ptb->i_send_data[rank] + send_count[rank]++];
    }
  }
  free (send_count);
  free (part_weight);
}

/**
 *
 * \brief  Define active ranks
 *
 * \param [inout]   ptb          Part to block structure
 *
 */

static
PDM_part_to_block_t *
_ptb_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  if (n_ptb == 0) {
    t_timer[0] = PDM_timer_create ();
    t_timer[1] = PDM_timer_create ();
    t_timer[2] = PDM_timer_create ();
  }
  n_ptb++;

  PDM_part_to_block_t *ptb = (PDM_part_to_block_t *) malloc (sizeof(PDM_part_to_block_t));

  ptb->t_distrib         = t_distrib;    /*!< Distribution type */
  ptb->t_post            = t_post;       /*!< Post processing type */
  ptb->n_active_ranks    = 0;            /*!< Number of active ranks */
  ptb->active_ranks      = NULL;         /*!< List of active ranks */
  ptb->comm              = comm;         /*!< MSG communicator */
  PDM_MPI_Comm_size (comm, &(ptb->s_comm));
  PDM_MPI_Comm_rank (comm, &(ptb->i_rank));
  ptb->is_my_rank_active   = 0;              /*!< Is active current rank */
  ptb->part_active_node    = part_active_node; /*!< Part of active nodes */

  ptb->n_part             = n_part;       /*!< Number of parts */
  ptb->n_elt              = n_elt;        /*!< Number of elements for any part */
  ptb->n_elt_proc         = 0;            /*!< Number of elements on the current rank */
  ptb->gnum_elt           = gnum_elt;     /*!< Global numbering of elements for any part */
  ptb->weight             = weight;
  ptb->dest_proc          = NULL;
  ptb->data_distrib_index = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (ptb->s_comm + 1));   /*!< Data distribution on ranks */

  ptb->s_block_min   = INT_MAX;
  ptb->s_block_max   = 0;

  ptb->i_send_data   = NULL;  /*!< Data to send to other processes index (size = s_comm) */
  ptb->i_recv_data   = NULL;  /*!< Received Data from other processes index (size = s_comm) */
  ptb->n_send_data   = NULL;  /*!< Number of data to send to other processes (size = s_comm) */
  ptb->n_recv_data   = NULL;  /*!< Number of received Data from other processes (size = s_comm) */

  ptb->tn_send_data      = 0;     /*!< Total number of sended data */
  ptb->tn_recv_data      = 0;     /*!< Total number of received data */
  ptb->sorted_recv_gnum  = NULL;  /*!< Sorted recv global num */
  ptb->order             = NULL;  /*!< Order */
  ptb->n_elt_block       = 0;
  ptb->block_gnum        = NULL;  /*!< Global number of reveived data (size = tn_recv_data) */
  ptb->block_gnum_count  = NULL;  /*!< Number of occurence of reveived data (size = tn_recv_data) */

  ptb->weight_g          = NULL; /*!< Global weights of elements for any part */
  ptb->enable_reverse    = 0;
  ptb->idx_partial       = NULL;

  /* Asynchone */
  ptb->max_exch_request = 10;
  ptb->next_request     = 0;
  ptb->s_data           = (size_t          * ) malloc ( ptb->max_exch_request * sizeof(size_t           ) );
  ptb->t_stride         = (PDM_stride_t    * ) malloc ( ptb->max_exch_request * sizeof(PDM_stride_t     ) );
  ptb->cst_stride       = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->wait_status      = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->request_mpi      = (PDM_MPI_Request * ) malloc ( ptb->max_exch_request * sizeof(PDM_MPI_Request  ) );

  ptb->send_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_stride      = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );

  ptb->block_stride     = (int           *** ) malloc ( ptb->max_exch_request * sizeof(int            **) );
  ptb->block_data       = (void          *** ) malloc ( ptb->max_exch_request * sizeof(void           **) );

  ptb->part_stride      = (int          **** ) malloc ( ptb->max_exch_request * sizeof(int           ***) );
  ptb->part_data        = (void         **** ) malloc ( ptb->max_exch_request * sizeof(void          ***) );

  ptb->comm_kind        = (PDM_mpi_comm_kind_t *) malloc ( ptb->max_exch_request * sizeof(PDM_mpi_comm_kind_t  ));
  ptb->win_send         = (PDM_MPI_Win         *) malloc ( ptb->max_exch_request * sizeof(PDM_MPI_Win         *));
  ptb->win_recv         = (PDM_MPI_Win         *) malloc ( ptb->max_exch_request * sizeof(PDM_MPI_Win         *));

  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    ptb->send_buffer  [i_req] = NULL;
    ptb->recv_buffer  [i_req] = NULL;
    ptb->recv_stride  [i_req] = NULL;
    ptb->n_send_buffer[i_req] = NULL;
    ptb->i_send_buffer[i_req] = NULL;
    ptb->n_recv_buffer[i_req] = NULL;
    ptb->i_recv_buffer[i_req] = NULL;
    ptb->block_stride [i_req] = NULL;
    ptb->block_data   [i_req] = NULL;
    ptb->part_stride  [i_req] = NULL;
    ptb->part_data    [i_req] = NULL;
    ptb->win_send     [i_req] = PDM_MPI_WIN_NULL;
    ptb->win_recv     [i_req] = PDM_MPI_WIN_NULL;
  }

  /*
   * Active ranks definition
   */
  _active_ranks (ptb);

  return ptb;
}


static
void
_prepare_exchange
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                  **part_stride,
  size_t                *i_send_buffer,
  size_t                *i_recv_buffer,
  int                   *n_send_buffer,
  int                   *n_recv_buffer,
  int                  **recv_stride
)
{
  /*
   * Exchange Stride and build buffer properties
   */
  int *_recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *send_stride = (int *) malloc (sizeof(int) * ptb->tn_send_data);

    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        send_stride[ptb->i_send_data[iproc] + n_send_buffer[iproc]] = part_stride[i][j];
        n_send_buffer[iproc] += 1;
      }
    }

    _recv_stride = (int *) malloc (sizeof(int) * ptb->tn_recv_data);
    assert (send_stride != NULL);
    assert (_recv_stride != NULL);

    *recv_stride = _recv_stride;

    PDM_MPI_Alltoallv (send_stride,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_INT,
                       _recv_stride,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_INT,
                       ptb->comm);

    /*
     * Build buffers
     */

    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_send_data[i];
      int iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_send_buffer[i] += send_stride[k];

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_recv_data[i];
      iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_recv_buffer[i] += _recv_stride[k];

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0)
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      else
        i_recv_buffer[i] = 0;

    }

    free(send_stride);

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_send_data[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = ptb->i_recv_data[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = ptb->n_send_data[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = ptb->n_recv_data[i] * cst_stride * (int) s_data;

    }
  }
}


static
void
_prepare_reverse_exchange
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                   *block_strid,
  size_t                *i_send_buffer,
  size_t                *i_recv_buffer,
  int                   *n_send_buffer,
  int                   *n_recv_buffer,
  int                  **send_stride,
  int                  **recv_stride
)
{
  /*
   * Exchange Stride and build buffer properties
   */
  int *_recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *_send_stride = (int *) malloc (sizeof(int) * ptb->tn_recv_data);
    *send_stride = _send_stride;

    /* Because it's reverse we have tn_recn_data for send -_-' */
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      _send_stride[ptb->order[i]] = block_strid[ptb->idx_partial[i]];
    }
    // PDM_log_trace_array_int(_send_stride, ptb->tn_recv_data, "_send_stride :: ");

    int s_recv = ptb->i_send_data[ptb->s_comm-1] + ptb->n_send_data[ptb->s_comm-1];
    _recv_stride = (int *) malloc (sizeof(int) * s_recv);
    assert (_send_stride != NULL);
    assert (_recv_stride != NULL);

    *recv_stride = _recv_stride;

    PDM_MPI_Alltoallv (_send_stride,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_INT,
                       _recv_stride,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_INT,
                       ptb->comm);

    if(0 == 1) {
      PDM_log_trace_array_int(_recv_stride, s_recv, "recv_stride :: ");
    }

    /*
     * Build buffers
     */
    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_recv_data[i];
      int iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++){
        n_send_buffer[i] += _send_stride[k];
      }

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_send_data[i];
      iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++){
        n_recv_buffer[i] += _recv_stride[k];
      }

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }

    }

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_recv_data[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = ptb->i_send_data[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = ptb->n_recv_data[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = ptb->n_send_data[i] * cst_stride * (int) s_data;

    }

  }

  if(0 == 1) {
    PDM_log_trace_array_size_t(i_send_buffer, ptb->s_comm, "i_send_buffer ::");
    PDM_log_trace_array_size_t(i_recv_buffer, ptb->s_comm, "i_recv_buffer ::");
    PDM_log_trace_array_int   (n_send_buffer, ptb->s_comm, "n_send_buffer ::");
    PDM_log_trace_array_int   (n_recv_buffer, ptb->s_comm, "n_recv_buffer ::");
  }
}








static
void
_prepare_send_buffer
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                  **part_stride,
  void                 **part_data,
  size_t                *i_send_buffer,
  int                   *n_send_buffer,
  unsigned char         *send_buffer
)
{
  unsigned char **_part_data = (unsigned char **) part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_send_buffer[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    size_t *i_part = NULL;
    if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
      i_part = (size_t *) malloc (sizeof(size_t) * (ptb->n_elt[i] + 1));
      assert (i_part != NULL);

    i_part[0] = 0;
      for (int j = 1; j < ptb->n_elt[i] + 1; j++)
        i_part[j] = i_part[j-1] + ((size_t) part_stride[i][j-1] * s_data);
    }

    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[++idx];
      size_t s_octet_elt = 0;
      size_t i_part_elt = 0;

      if (t_stride == PDM_STRIDE_CST_INTERLACED) {
        s_octet_elt = cst_stride * (int) s_data;
        i_part_elt  = cst_stride * (int) s_data * j;
      }

      else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
        s_octet_elt = i_part[j+1] - i_part[j];
        i_part_elt  = i_part[j];
      }

      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[i_send_buffer[iproc] + n_send_buffer[iproc]++] =
          _part_data[i][i_part_elt + k];
      }


    }

    if (i_part != NULL)
      free (i_part);
  }
}



static
void
_prepare_reverse_send_buffer
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                   *send_stride,
  int                   *block_stride,
  void                  *block_data,
  size_t                *i_send_buffer,
  int                   *n_send_buffer,
  unsigned char         *send_buffer
)
{
  // unsigned char **_part_data = (unsigned char **) part_data;
  unsigned char *_block_data = (unsigned char *) block_data;

  int s_send_buffer = i_send_buffer[ptb->s_comm-1] + n_send_buffer[ptb->s_comm-1];

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    /* Because it's reverse we have tn_recn_data for send -_-' */
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      send_stride[ptb->order[i]] = block_stride[ptb->idx_partial[i]];
    }

    int *send_idx = malloc((ptb->tn_recv_data+1) * sizeof(int));
    send_idx[0] = 0;
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      send_idx[i+1] = send_idx[i] + send_stride[i];
      // send_stride[i] = 0;
    }

    // Same for block_stride
    int *block_idx = malloc((ptb->n_elt_block+1) * sizeof(int));
    block_idx[0] = 0;
    for(int i = 0; i < ptb->n_elt_block; ++i) {
      block_idx[i+1] = block_idx[i] + block_stride[i];
    }


    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      size_t s_octet_elt = send_stride[ptb->order[i]] * s_data;
      int idx_write = send_idx [ptb->order      [i]] * s_data;
      int idx_read  = block_idx[ptb->idx_partial[i]] * s_data;
      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[idx_write+k] = _block_data[idx_read+k];
      }
    }

    free(send_stride);
    free(send_idx);
    free(block_idx);

  }
  else  {
    size_t s_octet_elt = cst_stride * (int) s_data;
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[s_octet_elt * ptb->order[i] + k] = _block_data[s_octet_elt*ptb->idx_partial[i]+k];
      }
    }
  }

  if(0 == 1) {
    int* tmp_send_buffer = (int *) send_buffer;
    PDM_log_trace_array_int(tmp_send_buffer, s_send_buffer/sizeof(int), "tmp_send_buffer ::");
  }

}



/**
 *
 * \brief  Post-treatment of the resulting buffer
 *
 * \param [inout]   ptb          Part to block structure
 *
 */
static
int
_post_treatment
(
  PDM_part_to_block_t  *ptb,
  size_t                s_data,
  PDM_stride_t          t_stride,
  int                   cst_stride,
  int                  *recv_stride,
  unsigned char        *recv_buffer,
  size_t                s_recv_buffer,
  int                 **block_stride,
  void                **block_data
)
{
  unsigned char *_block_data = malloc(sizeof(unsigned char) * s_recv_buffer);
  assert(_block_data != NULL);
  *block_data = _block_data;

  int *i_recv_stride = NULL;
  int *i_block_stride = NULL;
  int s_block_data = ((int) sizeof(unsigned char) * s_recv_buffer) / (int) s_data;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    *block_stride = NULL;
    int* _block_stride = NULL;
    if(ptb->tn_recv_data > 0){
      _block_stride = malloc(sizeof(int) * ptb->tn_recv_data);
    }
    *block_stride = _block_stride;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      _block_stride[i] = recv_stride[ptb->order[i]];
    }

    /*
     * Compute index in data
     */
    i_recv_stride  = malloc (sizeof(int) * (ptb->tn_recv_data + 1));
    i_block_stride = malloc (sizeof(int) * (ptb->tn_recv_data + 1));

    i_recv_stride [0] = 0;
    i_block_stride[0] = 0;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride [i+1] = i_recv_stride [i] + recv_stride[i];
      i_block_stride[i+1] = i_block_stride[i] + _block_stride[i];
    }

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride [i+1] *= (int) s_data;
      i_block_stride[i+1] *= (int) s_data;
    }

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int old    = ptb->order[i];
      int id_old = i_recv_stride[old];
      for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
        _block_data[k] = recv_buffer[id_old++];
      }
    }

    free (recv_stride);
    free (i_recv_stride);

    /*
     * post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

      int idx1 = 0;
      int idx2 = 0;

      if (ptb->tn_recv_data == 1) {
        idx2 = i_block_stride[1];
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        if (i == 1) {
          idx2 = i_block_stride[1];
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          _block_stride[idx1] = _block_stride[i];
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
              _block_data[idx2++] = _block_data[k];
            }
          }
        } else {
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
            _block_stride[idx1] += _block_stride[i];
          }
        }
      }

      /* Cleanup */

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;

        _block_stride = realloc (_block_stride, sizeof(int) * ptb->n_elt_block);

        *block_stride = _block_stride;
        s_block_data = idx2 / (int) s_data;
      }
    }
    free (i_block_stride);
  } else {

    /*
     * Sort Buffer
     */
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int n_octet = cst_stride * (int) s_data;
      int old = ptb->order[i];
      int id_old = old * n_octet;

      for (int k = i*n_octet; k < (i+1)*n_octet; k++) {
        _block_data[k] = recv_buffer[id_old++];
      }
    }

    /*
     * Post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {
      int idx2 = 0;
      int idx1 = 0;

      assert (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE);

      if (ptb->tn_recv_data == 1) {
        idx2 =  cst_stride * (int) s_data;
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        int n_octet = cst_stride * (int) s_data;
        if (i == 1) {
          idx2 = n_octet;
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            int idx3 = i * cst_stride * (int) s_data;
            for (int k = 0; k < n_octet; k++) {
              _block_data[idx2++] = _block_data[idx3++];
            }
          }
        }
      }

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;
        s_block_data = idx2 / (int) s_data;
      }

    }
  }

  return s_block_data;
}

/**
 *
 * \brief  Post-treatment of the resulting buffer
 *
 * \param [inout]   ptb          Part to block structure
 *
 */
static
void
_post_treatment_reverse
(
  PDM_part_to_block_t  *ptb,
  size_t                s_data,
  PDM_stride_t          t_stride,
  int                   cst_stride,
  int                  *recv_stride,
  unsigned char        *recv_buffer,
  int                  *n_recv_buffer,
  size_t               *i_recv_buffer,
  int                ***part_stride,
  void               ***part_data
)
{
  *part_data = malloc(ptb->n_part * sizeof(unsigned char *));
  unsigned char **_part_data = (unsigned char **) *part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_recv_buffer[i] = 0;
  }

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    *part_stride = malloc(ptb->n_part * sizeof(int *));
    int **_part_stride = *part_stride;
    int **_part_idx    = malloc(ptb->n_part * sizeof(int *));

    int *n_recv_strid = malloc(ptb->s_comm * sizeof(int));
    for(int i = 0; i < ptb->s_comm; ++i) {
      n_recv_strid[i] = 0;
    }

    /*
     * Post-treated recv_strid
     */
    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      _part_stride[i] = malloc( ptb->n_elt[i]    * sizeof(int));
      _part_idx   [i] = malloc((ptb->n_elt[i]+1) * sizeof(int));

      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        int idx_read = ptb->i_send_data[iproc] + n_recv_strid[iproc]++;
        _part_stride[i][j] = recv_stride[idx_read];
      }

      _part_idx[i][0] = 0;
      for(int j = 0; j < ptb->n_elt[i]; j++) {
        _part_idx[i][j+1] = _part_idx[i][j] + _part_stride[i][j];
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_part_idx[i], ptb->n_elt[i]+1, "_part_idx :: ");
      }

    }

    idx = -1;
    int n_octet = (int) s_data;
    for (int i = 0; i < ptb->n_part; i++) {
      _part_data[i] = malloc(_part_idx[i][ptb->n_elt[i]] * n_octet * sizeof(unsigned char));
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        int idx_write = n_octet * _part_idx[i][j];
        for (int k = 0; k < (int) n_octet * _part_stride[i][j]; k++) {
          int idx_read = i_recv_buffer[iproc] + n_recv_buffer[iproc]++;
          _part_data[i][idx_write + k] = recv_buffer[idx_read];
        }
      }
      free(_part_idx[i]);
    }
    free(n_recv_strid);
    free(_part_idx);
    free(recv_stride);

  } else { // PDM_STRIDE_CST_INTERLACED

    int n_octet = cst_stride * (int) s_data;
    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      _part_data[i] = malloc(ptb->n_elt[i] * n_octet * sizeof(unsigned char));
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        for (int k = 0; k < (int) n_octet; k++) {
          int idx_read = i_recv_buffer[iproc] + n_recv_buffer[iproc]++;
          _part_data[i][n_octet*j + k] = recv_buffer[idx_read];
        }
      }
    }
  }
}





/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Reset global statistic
 *
 */

void
PDM_part_to_block_global_statistic_reset
(
  void
)
{
  for (int i = 0; i < 3; i++) {
    t_elaps[i] = 0;
    t_cpu[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    min_exch_rank[i] = INT_MAX;
    max_exch_rank[i] = -1;
    exch_data[i] = 0;
  }
}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm                 MPI communicator
 * \param [out]  min_exch_rank_send   Global min part of ranks used to send
 * \param [out]  min_exch_rank_recv   Global min part of ranks used to receive
 * \param [out]  max_exch_rank_send   Global max part of ranks used to send
 * \param [out]  max_exch_rank_recv   Global max part of ranks used to receive
 * \param [out]  min_exch_data_send   Global min sent data for a rank
 * \param [out]  min_exch_data_recv   Global min received data for a rank
 * \param [out]  max_exch_data_send   Global max sent data for a rank
 * \param [out]  max_exch_data_recv   Global max received data for a rank
 *
 */

void
PDM_part_to_block_global_statistic_get
(
 PDM_MPI_Comm comm,
 int *min_exch_rank_send,
 int *min_exch_rank_recv,
 int *max_exch_rank_send,
 int *max_exch_rank_recv,
 unsigned long long *min_exch_data_send,
 unsigned long long *min_exch_data_recv,
 unsigned long long *max_exch_data_send,
 unsigned long long *max_exch_data_recv
)
{
  unsigned long long max_exch_data[2];
  unsigned long long min_exch_data[2];

  PDM_MPI_Allreduce (exch_data, min_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (exch_data, max_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MAX, comm);

  *min_exch_data_send = min_exch_data[0];
  *min_exch_data_recv = min_exch_data[1];
  *max_exch_data_send = max_exch_data[0];
  *max_exch_data_recv = max_exch_data[1];


  int max_max_exch_rank[2];
  int min_min_exch_rank[2];

  PDM_MPI_Allreduce (min_exch_rank, min_min_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (max_exch_rank, max_max_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MAX, comm);

  *min_exch_rank_send = min_min_exch_rank[0];
  *min_exch_rank_recv = min_min_exch_rank[1];
  *max_exch_rank_send = max_max_exch_rank[0];
  *max_exch_rank_recv = max_max_exch_rank[1];

}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm              MPI communicator
 * \param [out]  min_elaps         Min elapsed time
 * \param [out]  max_elaps         Max elapsed time
 * \param [out]  min_cpu           Min cpu time
 * \param [out]  max_cpu           Max cpu time
 * \param [out]  min_elaps_create  Global min elapsed for create function
 * \param [out]  max_elaps_create  Global max elapsed for create function
 * \param [out]  min_cpu_create    Global min cpu for create function
 * \param [out]  max_cpu_create    Global max cpu for create function
 * \param [out]  min_elaps_create2 Global min elapsed for create2 function
 * \param [out]  max_elaps_create2 Global max elapsed for create2 function
 * \param [out]  min_cpu_create2   Global min cpu for create2 function
 * \param [out]  max_cpu_create2   Global max cpu for create2 function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 *
 */

void
PDM_part_to_block_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_create2,
 double       *max_elaps_create2,
 double       *min_cpu_create2,
 double       *max_cpu_create2,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
)
{

  double min_elaps[3];
  double max_elaps[3];
  double min_cpu[3];
  double max_cpu[3];

  PDM_MPI_Allreduce (t_elaps, min_elaps, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_elaps, max_elaps, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  PDM_MPI_Allreduce (t_cpu, min_cpu, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_cpu, max_cpu, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  *min_elaps_create  = min_elaps[0];
  *max_elaps_create  = max_elaps[0];
  *min_cpu_create    = min_cpu[0];
  *max_cpu_create    = max_cpu[0];
  *min_elaps_create2 = min_elaps[1];
  *max_elaps_create2 = max_elaps[1];
  *min_cpu_create2   = min_cpu[1];
  *max_cpu_create2   = max_cpu[1];
  *min_elaps_exch    = min_elaps[2];
  *max_elaps_exch    = max_elaps[2];
  *min_cpu_exch      = min_cpu[2];
  *max_cpu_exch      = max_cpu[2];

}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */
  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         weight,
                                         n_elt,
                                         n_part,
                                         comm);

  double t1_elaps = PDM_timer_elapsed(t_timer[0]);
  double t1_cpu   = PDM_timer_cpu    (t_timer[0]);
  PDM_timer_resume(t_timer[0]);

  /*
   * Data distribution definition
   */

  _distrib_data (ptb, 0);


  /*
   * Compute global weight for each element
   */
  if (ptb->weight != NULL) {// && ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
    _compute_global_weights (ptb);
  }

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  PDM_timer_hang_on(t_timer[0]);
  double t2_elaps = PDM_timer_elapsed(t_timer[0]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[0]);

  t_elaps[0] += (t2_elaps - t1_elaps);
  t_cpu  [0] += (t2_cpu   - t1_cpu  );

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_from_distrib
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                         part_active_node,
 PDM_g_num_t                 **gnum_elt,
 const PDM_g_num_t            *data_distrib_index,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */
  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         NULL,
                                         n_elt,
                                         n_part,
                                         comm);

  for(int i_rank = 0; i_rank < ptb->s_comm+1; i_rank++){
    ptb->data_distrib_index[i_rank] = data_distrib_index[i_rank];
  }

  double t1_elaps = PDM_timer_elapsed(t_timer[1]);
  double t1_cpu = PDM_timer_cpu(t_timer[1]);
  PDM_timer_resume(t_timer[1]);

  /*
   * Data distribution definition
   */
  _distrib_data (ptb, 1);

  int n_rank_recv = 0;
  int n_rank_send = 0;


  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  PDM_timer_hang_on(t_timer[1]);
  double t2_elaps = PDM_timer_elapsed(t_timer[1]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[1]);

  t_elaps[1] += (t2_elaps - t1_elaps);
  t_cpu  [1] += (t2_cpu   - t1_cpu);

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_geom_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_part_geom_t               geom_kind,
 double                      **pvtx_coords,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */
  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         NULL,
                                         n_elt,
                                         n_part,
                                         comm);

  double t1_elaps = PDM_timer_elapsed(t_timer[0]);
  double t1_cpu = PDM_timer_cpu(t_timer[0]);
  PDM_timer_resume(t_timer[0]);

  if(geom_kind == PDM_PART_GEOM_HILBERT ) {
    _distrib_data_hilbert(ptb, pvtx_coords, weight);
  } else if (geom_kind == PDM_PART_GEOM_MORTON ){
    _distrib_data_morton(ptb, pvtx_coords, weight);
  }

  /*
   * Compute global weight for each element
   */
  if (ptb->weight != NULL) {// && ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
    _compute_global_weights (ptb);
  }

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  PDM_timer_hang_on(t_timer[0]);
  double t2_elaps = PDM_timer_elapsed(t_timer[0] );
  double t2_cpu = PDM_timer_cpu(t_timer[0]);

  t_elaps[0] += (t2_elaps - t1_elaps);
  t_cpu[0] += (t2_cpu - t1_cpu);

  return (PDM_part_to_block_t *) ptb;
}

/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of active ranks
 *
 */

int
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_active_ranks;
}


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  if current rank is active
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 )
{
  return ptb->is_my_rank_active;
}


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->active_ranks;
}


/**
 *
 * \brief Return number of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_elt_block;
}


/**
 *
 * \brief Return global numbers of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum;
}

/**
 *
 * \brief Return numbers of occurence of each gnum element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers counter
 *
 */

int *
PDM_part_to_block_block_gnum_count_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum_count;
}


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

int
PDM_part_to_block_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data,
 int                **block_stride,
 void               **block_data
)
{

  double t1_elaps = PDM_timer_elapsed(t_timer[2]);
  double t1_cpu = PDM_timer_cpu(t_timer[2]);
  PDM_timer_resume(t_timer[2]);

  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  int    *n_send_buffer = (int    *) malloc (sizeof(int   ) * ptb->s_comm);
  int    *n_recv_buffer = (int    *) malloc (sizeof(int   ) * ptb->s_comm);

  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    i_send_buffer,
                    i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);

  /*
   * Prepare buffer
   */
  size_t s_send_buffer = i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  size_t s_recv_buffer = i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  unsigned char *send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  unsigned char *recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Data exchange
   */
  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      ptb->comm);

  /*
   * Statistics
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i) {
      exch_data[1] += n_recv_buffer[i];
      exch_data[0] += n_send_buffer[i];
    }
  }

  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     block_stride,
                                     block_data);
  free (recv_buffer);

  PDM_timer_hang_on(t_timer[2]);
  double t2_elaps = PDM_timer_elapsed(t_timer[2]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[2]);

  t_elaps[2] += (t2_elaps - t1_elaps);
  t_cpu  [2] += (t2_cpu   - t1_cpu  );

  return s_block_data;
}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_reverse_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
)
{
  assert(ptb->enable_reverse == 1);

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  int    *n_send_buffer = (int    *) malloc (sizeof(int   ) * ptb->s_comm);
  int    *n_recv_buffer = (int    *) malloc (sizeof(int   ) * ptb->s_comm);

  /*
   * Exchange Stride and build buffer properties
   */
  int *send_stride = NULL; // Release in post-treatment
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_reverse_exchange(ptb,
                            s_data,
                            t_stride,
                            cst_stride,
                            block_stride,
                            i_send_buffer,
                            i_recv_buffer,
                            n_send_buffer,
                            n_recv_buffer,
                            &send_stride,
                            &recv_stride);

  /*
   * Prepare buffer
   */
  size_t s_send_buffer = i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  size_t s_recv_buffer = i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  unsigned char *send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  unsigned char *recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_reverse_send_buffer(ptb,
                               s_data,
                               t_stride,
                               cst_stride,
                               send_stride,
                               block_stride,
                               block_data,
                               i_send_buffer,
                               n_send_buffer,
                               send_buffer);
  // if(send_stride != NULL) {
  //   free(send_stride);
  // }

  /*
   * Data exchange
   */
  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      ptb->comm);



  _post_treatment_reverse(ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          recv_stride,
                          recv_buffer,
                          n_recv_buffer,
                          i_recv_buffer,
                          part_stride,
                          part_data);

  free (recv_buffer);
  free (send_buffer);

  free(i_send_buffer);
  free(i_recv_buffer);
  free(n_send_buffer);
  free(n_recv_buffer);
}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                 **part_stride,
       void                **part_data,
       int                 **block_stride,
       void                **block_data,
       int                  *request
)
{
  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;
  *request = request_id;

  ptb->i_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->i_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);

  ptb->block_stride[request_id] = block_stride;
  ptb->block_data  [request_id] = block_data;


  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;
  ptb->comm_kind [request_id] = k_comm;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *tmp_i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    tmp_i_send_buffer,
                    tmp_i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  int s_recv_buffer = tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  if(k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    PDM_MPI_Win_allocate(s_send_buffer, sizeof(unsigned char), ptb->comm, &ptb->send_buffer[request_id], &ptb->win_send[request_id]);
    PDM_MPI_Win_allocate(s_recv_buffer, sizeof(unsigned char), ptb->comm, &ptb->recv_buffer[request_id], &ptb->win_recv[request_id]);

  } else {

    ptb->send_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    ptb->recv_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);
  }
  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       tmp_i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  free(tmp_i_send_buffer);
  free(tmp_i_recv_buffer);

  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_P2P k_comm is not implemented yet\n");
    abort();
  } else if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    PDM_MPI_Ialltoallv(send_buffer,
                       n_send_buffer,
                       i_send_buffer,
                       PDM_MPI_BYTE,
                       recv_buffer,
                       n_recv_buffer,
                       i_recv_buffer,
                       PDM_MPI_BYTE,
                       ptb->comm,
                       &ptb->request_mpi[request_id]);
  } else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    // double t1 = PDM_MPI_Wtime();
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
    PDM_MPI_Get_ialltoallv(ptb->win_send[request_id],
                           ptb->win_recv[request_id],
                           send_buffer,
                           n_send_buffer,
                           i_send_buffer,
                           PDM_MPI_BYTE,
                           recv_buffer,
                           n_recv_buffer,
                           i_recv_buffer,
                           PDM_MPI_BYTE,
                           ptb->comm);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_MPI_Get_ialltoallv + fence dt = %12.5e \n", dt);

  }


  ptb->wait_status[request_id] = 0;
}



/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_reverse_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                  *block_stride,
       void                 *block_data,
       int                ***part_stride,
       void               ***part_data,
       int                  *request
)
{
  assert(ptb->enable_reverse == 1);

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;
  *request = request_id;

  ptb->i_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->i_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);

  ptb->part_stride[request_id] = part_stride;
  ptb->part_data  [request_id] = part_data;


  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;
  ptb->comm_kind [request_id] = k_comm;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *tmp_i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  /*
   * Exchange Stride and build buffer properties
   */
  int *send_stride = NULL; // Release in post-treatment
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_reverse_exchange(ptb,
                            s_data,
                            t_stride,
                            cst_stride,
                            block_stride,
                            tmp_i_send_buffer,
                            tmp_i_recv_buffer,
                            n_send_buffer,
                            n_recv_buffer,
                            &send_stride,
                            &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  int s_recv_buffer = tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  if(k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    PDM_MPI_Win_allocate(s_send_buffer, sizeof(unsigned char), ptb->comm, &ptb->send_buffer[request_id], &ptb->win_send[request_id]);
    PDM_MPI_Win_allocate(s_recv_buffer, sizeof(unsigned char), ptb->comm, &ptb->recv_buffer[request_id], &ptb->win_recv[request_id]);

  } else {

    ptb->send_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    ptb->recv_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);
  }
  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_reverse_send_buffer(ptb,
                               s_data,
                               t_stride,
                               cst_stride,
                               send_stride,
                               block_stride,
                               block_data,
                               tmp_i_send_buffer,
                               n_send_buffer,
                               send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  free(tmp_i_send_buffer);
  free(tmp_i_recv_buffer);

  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();
  } else if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    PDM_MPI_Ialltoallv(send_buffer,
                       n_send_buffer,
                       i_send_buffer,
                       PDM_MPI_BYTE,
                       recv_buffer,
                       n_recv_buffer,
                       i_recv_buffer,
                       PDM_MPI_BYTE,
                       ptb->comm,
                       &ptb->request_mpi[request_id]);
  } else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    // double t1 = PDM_MPI_Wtime();
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
    PDM_MPI_Get_ialltoallv(ptb->win_send[request_id],
                           ptb->win_recv[request_id],
                           send_buffer,
                           n_send_buffer,
                           i_send_buffer,
                           PDM_MPI_BYTE,
                           recv_buffer,
                           n_recv_buffer,
                           i_recv_buffer,
                           PDM_MPI_BYTE,
                           ptb->comm);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_MPI_Get_ialltoallv + fence dt = %12.5e \n", dt);

  }


  ptb->wait_status[request_id] = 0;
}




/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 * \return       Size of highest block
 */
int
PDM_part_to_block_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
  } else {
    int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
    assert(code == PDM_MPI_SUCCESS);
  }

  ptb->wait_status[request_id] = 1;

  /*
   *  Post-treatment
   */
  size_t s_recv_buffer = ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1];

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_send[request_id]);
    ptb->win_send[request_id] = PDM_MPI_WIN_NULL;
  } else {
    free(ptb->send_buffer  [request_id]);
  }
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];

  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  // recv_stride is free inside
  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     ptb->block_stride[request_id],
                                     ptb->block_data  [request_id]);

  /*
   * Free
   */

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_recv[request_id]);
    ptb->win_recv[request_id] = PDM_MPI_WIN_NULL;
  } else {
    free(ptb->recv_buffer  [request_id]);
  }
  ptb->recv_stride [request_id] = NULL;
  ptb->recv_buffer [request_id] = NULL;
  ptb->block_stride[request_id] = NULL;
  ptb->block_data  [request_id] = NULL;


  ptb->wait_status[request_id] = 2;

  return s_block_data;
}



/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 * \return       Size of highest block
 */
void
PDM_part_to_block_reverse_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
  } else {
    int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
    assert(code == PDM_MPI_SUCCESS);
  }

  ptb->wait_status[request_id] = 1;

  /*
   *  Post-treatment
   */
  // size_t s_recv_buffer = ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1];


  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_send[request_id]);
    ptb->win_send[request_id] = PDM_MPI_WIN_NULL;
  } else {
    free(ptb->send_buffer  [request_id]);
  }
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;

  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];

  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  size_t *tmp_i_recv_buffer = malloc((ptb->s_comm+1) * sizeof(size_t));
  tmp_i_recv_buffer[0] = 0;
  for(int i = 0; i < ptb->s_comm; ++i) {
    tmp_i_recv_buffer[i+1] = tmp_i_recv_buffer[i] + ptb->n_recv_buffer[request_id][i];
  }

  _post_treatment_reverse(ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          recv_stride,
                          recv_buffer,
                          ptb->n_recv_buffer[request_id],
                          tmp_i_recv_buffer,
                          ptb->part_stride  [request_id],
                          ptb->part_data    [request_id]);

  free(tmp_i_recv_buffer);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  /*
   * Free
   */

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_recv[request_id]);
    ptb->win_recv[request_id] = PDM_MPI_WIN_NULL;
  } else {
    free(ptb->recv_buffer  [request_id]);
  }
  ptb->recv_stride [request_id] = NULL;
  ptb->recv_buffer [request_id] = NULL;
  ptb->part_stride [request_id] = NULL;
  ptb->part_data   [request_id] = NULL;

  ptb->wait_status[request_id] = 2;

}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */
int
PDM_part_to_block_async_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data
)
{

  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;

  ptb->i_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->i_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *tmp_i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    tmp_i_send_buffer,
                    tmp_i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  int s_recv_buffer = tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  ptb->send_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  ptb->recv_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       tmp_i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  free(tmp_i_send_buffer);
  free(tmp_i_recv_buffer);

  PDM_MPI_Ialltoallv(send_buffer,
                     n_send_buffer,
                     i_send_buffer,
                     PDM_MPI_BYTE,
                     recv_buffer,
                     n_recv_buffer,
                     i_recv_buffer,
                     PDM_MPI_BYTE,
                     ptb->comm,
                     &ptb->request_mpi[request_id]);

  ptb->wait_status[request_id] = 0;

  return request_id;
}


/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 */
void
PDM_part_to_block_async_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
  assert(code == PDM_MPI_SUCCESS);

  ptb->wait_status[request_id] = 1;

}

/**
 *
 * \brief Get the raw exchange buffer and stride and deallocate memory
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_get_raw
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  // printf("PDM_part_to_block_asyn_get_raw::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 1);

  free(ptb->send_buffer  [request_id]);
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  /* Mv pointer */
  *block_stride = ptb->recv_stride[request_id];
  *block_data   = ptb->recv_buffer[request_id];

  /* Nulliffy - User is now owner of the excahnge data */
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  return ptb->tn_recv_data;

}


/**
 *
 * \brief Post-treatment of the recv buffer
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_post_treatment
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  assert(ptb->wait_status[request_id] == 1);

  // size_t s_send_buffer = ptb->i_send_buffer[request_id][ptb->s_comm - 1] + ptb->n_send_buffer[request_id][ptb->s_comm -1];
  size_t s_recv_buffer = ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1];

  free(ptb->send_buffer  [request_id]);
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];

  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  // recv_stride is free inside
  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     block_stride,
                                     block_data);

  /*
   * Free
   */
  free(ptb->recv_buffer  [request_id]);
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  return s_block_data;
}


/**
 *
 * \brief Free a part to block structure
 *
 * \param [inout] ptb         Part to block structure
 *
 * \return       NULL
 */

PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
)
{

  if (ptb->active_ranks != NULL) {
    free (ptb->active_ranks);
    ptb->active_ranks = NULL;
  }
  if (ptb->dest_proc != NULL) {
    free (ptb->dest_proc);
    ptb->dest_proc = NULL;
  }
  if (ptb->data_distrib_index != NULL) {
    free (ptb->data_distrib_index);
    ptb->data_distrib_index = NULL;
  }
  if (ptb->i_send_data != NULL) {
    free (ptb->i_send_data);
    ptb->i_send_data = NULL;
  }
  if (ptb->i_recv_data != NULL) {
    free (ptb->i_recv_data);
    ptb->i_recv_data = NULL;
  }
  if (ptb->n_send_data != NULL) {
    free (ptb->n_send_data);
    ptb->n_send_data = NULL;
  }
  if (ptb->n_recv_data != NULL) {
    free (ptb->n_recv_data);
    ptb->n_recv_data = NULL;
  }
  if (ptb->sorted_recv_gnum != NULL) {
    free (ptb->sorted_recv_gnum);
    ptb->sorted_recv_gnum = NULL;
  }

  if (ptb->enable_reverse == 1 && ptb->idx_partial != NULL) {
    free (ptb->idx_partial);
    ptb->idx_partial = NULL;
  }
  if (ptb->order != NULL) {
    free (ptb->order);
    ptb->order = NULL;
  }

  // if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) && (ptb->block_gnum != NULL)) {
  //   free (ptb->block_gnum);
  //   ptb->block_gnum = NULL;
  // }
  if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING      ) &&
      (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) && (ptb->block_gnum != NULL)) {
    free (ptb->block_gnum);
    ptb->block_gnum = NULL;
  }
  if (ptb->block_gnum_count != NULL) {
    free(ptb->block_gnum_count);
    ptb->block_gnum_count=NULL;
  }

  if (ptb->weight_g != NULL) {
    for (int i = 0; i < ptb->n_part; i++) {
      free (ptb->weight_g[i]);
    }
    free (ptb->weight_g);
  }

  /* This one check if all buffer has been correctly move or delete */
  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    assert(ptb->send_buffer  [i_req] == NULL);
    assert(ptb->recv_buffer  [i_req] == NULL);
    assert(ptb->recv_stride  [i_req] == NULL);
    assert(ptb->n_send_buffer[i_req] == NULL);
    assert(ptb->i_send_buffer[i_req] == NULL);
    assert(ptb->n_recv_buffer[i_req] == NULL);
    assert(ptb->i_recv_buffer[i_req] == NULL);
    assert(ptb->block_stride [i_req] == NULL);
    assert(ptb->block_data   [i_req] == NULL);
    assert(ptb->part_stride  [i_req] == NULL);
    assert(ptb->part_data    [i_req] == NULL);
    assert(ptb->win_send     [i_req] == PDM_MPI_WIN_NULL);
    assert(ptb->win_recv     [i_req] == PDM_MPI_WIN_NULL);
  }

  free(ptb->s_data       );
  free(ptb->t_stride     );
  free(ptb->cst_stride   );
  free(ptb->wait_status  );
  free(ptb->request_mpi  );
  free(ptb->send_buffer  );
  free(ptb->recv_buffer  );
  free(ptb->recv_stride  );
  free(ptb->n_send_buffer);
  free(ptb->i_send_buffer);
  free(ptb->n_recv_buffer);
  free(ptb->i_recv_buffer);
  free(ptb->block_stride );
  free(ptb->block_data   );
  free(ptb->part_stride  );
  free(ptb->part_data    );

  free(ptb->comm_kind);
  free(ptb->win_send );
  free(ptb->win_recv );

  free (ptb);

  n_ptb--;
  if (n_ptb == 0) {
    PDM_timer_free(t_timer[0]);
    PDM_timer_free(t_timer[1]);
    PDM_timer_free(t_timer[2]);
  }

  return NULL;
}


/**
 *
 * \brief Return block distribution index
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Distribution (size = communicator size + 1)
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->data_distrib_index;
}


/**
 *
 * \brief Return processus destination
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Destination (size = sum of partition elements)
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->dest_proc;
}

PDM_g_num_t*
PDM_part_to_block_adapt_partial_block_to_block
(
 PDM_part_to_block_t  *ptb,
 int                 **block_n,
 PDM_g_num_t           n_g_block
)
{
  PDM_g_num_t *_block_distrib_idx = malloc (sizeof(PDM_g_num_t) * (ptb->s_comm + 1));

  for (int i = 0; i < ptb->s_comm + 1; i++) {
    _block_distrib_idx[i] = ptb->data_distrib_index[i];
  }
  int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);

  int block_n_elt_tot = _block_distrib_idx[ptb->i_rank+1] - _block_distrib_idx[ptb->i_rank];
  int* block_n_tmp = PDM_array_zeros_int(block_n_elt_tot);

  int* _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt; i1++) {
    int i = (int) (ptb->block_gnum[i1] - _block_distrib_idx[ptb->i_rank] - 1);
    // printf(" ptb->block_gnum[%i] = %i --> %i (%i)\n", i1, ptb->block_gnum[i1], _block_distrib_idx[ptb->i_rank], i);
    block_n_tmp[i] = _block_n[i1];
  }

  *block_n = realloc(*block_n, sizeof(int) * block_n_elt_tot);
  _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt_tot; i1++) {
    _block_n[i1] = block_n_tmp[i1];
  }

  free(block_n_tmp);

  PDM_g_num_t old_max = _block_distrib_idx[ptb->s_comm];
  PDM_g_num_t new_max = n_g_block;
  int diff_last = (int) (new_max - old_max);

  _block_distrib_idx[ptb->s_comm] = new_max;

  if (ptb->i_rank == (ptb->s_comm - 1)) {
    int new_size = block_n_elt_tot + diff_last;
    *block_n = realloc(*block_n, sizeof(int) * new_size);
    _block_n = *block_n;
    for (int i = block_n_elt_tot; i < new_size; i++) {
      _block_n[i] = 0;
    }
  }

  return _block_distrib_idx;
}





/**
 *
 * \brief Return global weights of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global weights
 *
 */

double **
PDM_part_to_block_global_weight_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->weight_g;
}



/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_block_n_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->s_comm;
}


/**
 *
 * \brief Return total number of element in the current process (summed over all partitions)
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Total number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_proc_get
(
 PDM_part_to_block_t *ptb
 )
 {
  return ptb->n_elt_proc;
 }

#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
