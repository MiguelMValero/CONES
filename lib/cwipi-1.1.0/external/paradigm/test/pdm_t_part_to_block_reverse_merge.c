#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"

#include "pdm_part_to_block.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */
static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Total number of elements .\n\n"
     "  -f      <level>  Frequency of extract .\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_elmt,
           int           *freq)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_elmt = atol(argv[i]);
        *n_elmt = (PDM_g_num_t) _n_elmt;
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        int _freq = atoi(argv[i]);
        *freq = _freq;
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  PDM_g_num_t n_elmt = 10; // Number of elements in global configuration
  int         freq   = 1;
  _read_args(argc,
             argv,
             &n_elmt,
             &freq);
  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   *  Each proc have an equi repartition of data :
   *    - Define a distribution
   *    - Generate data among all pb by random
   */

  PDM_g_num_t *distrib_init_elmt = PDM_compute_uniform_entity_distribution(comm, n_elmt);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib_init_elmt, n_rank+1, "distrib_init_elmt : ");
  }
  int n_part  = 1;
  int pn_elmt = (distrib_init_elmt[i_rank+1] - distrib_init_elmt[i_rank]) / freq ;

  PDM_g_num_t *pln_to_to_gn = malloc(pn_elmt * sizeof(PDM_g_num_t));
  int         *pfield       = malloc(pn_elmt * sizeof(int        ));
  int         *pstrid       = malloc(pn_elmt * sizeof(int        ));
  for(int i = 0; i < pn_elmt; ++i) {
    unsigned int seed = (unsigned int) (distrib_init_elmt[i_rank] + i);
    srand(seed);
    pln_to_to_gn[i] = (rand() % n_elmt) + 1;
    pfield      [i] = i_rank;
    pstrid      [i] = 1;
  }

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn, pn_elmt, "pln_to_to_gn : ");
  }


  if(0 == 1) {
    PDM_log_trace_array_int(pfield, pn_elmt, "pfield : ");
  }


  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &pln_to_to_gn,
                                                      NULL,
                                                      &pn_elmt,
                                                      n_part,
                                                      comm);

  /*
   *  We can know the distribution
   */
  int n_elmt_in_block         = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t* distrib_elmt   = PDM_part_to_block_distrib_index_get(ptb);
  const PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get   (ptb);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib_elmt, n_rank+1       , "distrib_elmt : ");
    PDM_log_trace_array_long(blk_gnum    , n_elmt_in_block, "blk_gnum     : ");
  }

  int* dfield_strid = NULL;
  int* dfield       = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &pstrid,
              (void **)  &pfield,
                         &dfield_strid,
              (void **)  &dfield);

  if(0 == 1) {
    PDM_log_trace_array_int(dfield      , n_elmt_in_block, "dfield       : ");
    PDM_log_trace_array_int(dfield_strid, n_elmt_in_block, "dfield_strid : ");
  }

  PDM_g_num_t* dfield_post = malloc(n_elmt_in_block * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_elmt_in_block; ++i) {
    // dfield_post[i] = i;
    dfield_post[i] = blk_gnum[i];
  }

  PDM_g_num_t** tmp_pfield_post = NULL;
  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(PDM_g_num_t),
                                 PDM_STRIDE_CST_INTERLACED,
                                 1,
                                 NULL,
                      (void **)  dfield_post,
                                 NULL,
                      (void ***) &tmp_pfield_post);


  PDM_g_num_t *pfield_post = tmp_pfield_post[0];
  free(tmp_pfield_post);

  // PDM_log_trace_array_long(pfield_post, pn_elmt, "pfield_post : ");

  /*
   * Check
   */
  for(int i = 0; i < pn_elmt; ++i) {
    assert(pfield_post[i] == pln_to_to_gn[i]);
  }


  free(pfield_post);


  /*
   * Stride Var check
   */
  for(int i = 0; i < n_elmt_in_block; ++i) {
    dfield_strid[i] = 1; // TO DO --> Generation aléatoire de strid entre 1 et 6 par exemple
  }

  int** tmp_pfield_post_strid = NULL;
  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(PDM_g_num_t),
                                 PDM_STRIDE_VAR_INTERLACED,
                                 1,
                                 dfield_strid,
                      (void **)  dfield_post,
                                 &tmp_pfield_post_strid,
                      (void ***) &tmp_pfield_post);

  pfield_post = tmp_pfield_post[0];
  int *pfield_post_strid = tmp_pfield_post_strid[0];
  free(tmp_pfield_post);
  free(tmp_pfield_post_strid);

  /*
   * Check
   */
  for(int i = 0; i < pn_elmt; ++i) {
    assert(pfield_post      [i] == pln_to_to_gn[i]);
    assert(pfield_post_strid[i] == 1);
  }
  free(pfield_post);
  free(pfield_post_strid);


  /*
   * Stride Var check more complex
   */
  int dn_data = 0;
  for(int i = 0; i < n_elmt_in_block; ++i) {
    dfield_strid[i] = (int) rand() % 10; // TO DO --> Generation aléatoire de strid entre 1 et 6 par exemple
    dn_data += dfield_strid[i];
  }


  dfield_post = realloc(dfield_post, dn_data * sizeof(PDM_g_num_t));
  int idx_write = 0;
  for(int i = 0; i < n_elmt_in_block; ++i) {
    for(int k = 0; k < dfield_strid[i]; ++k) {
      dfield_post[idx_write++] = blk_gnum[i];
    }
  }
  // PDM_log_trace_array_int (dfield_strid, n_elmt_in_block, "dfield_strid ::");
  // PDM_log_trace_array_long(dfield_post , idx_write      , "dfield_post ::");

  PDM_part_to_block_reverse_exch(ptb,
                                 sizeof(PDM_g_num_t),
                                 PDM_STRIDE_VAR_INTERLACED,
                                 1,
                                 dfield_strid,
                      (void **)  dfield_post,
                                 &tmp_pfield_post_strid,
                      (void ***) &tmp_pfield_post);

  pfield_post       = tmp_pfield_post[0];
  pfield_post_strid = tmp_pfield_post_strid[0];
  free(tmp_pfield_post);
  free(tmp_pfield_post_strid);

  if(0 == 1) {
    int s_data = 0;
    for(int i = 0; i < pn_elmt; ++i) {
      for(int k = 0; k < pfield_post_strid[i]; ++k) {
        PDM_g_num_t check = pfield_post[s_data++];
        assert(check == pln_to_to_gn[i]);
      }
    }
    // PDM_log_trace_array_long(pfield_post, s_data, "pfield_post : ");
  }


  free(pfield_post);
  free(pfield_post_strid);


  free(dfield_post);

  PDM_part_to_block_free(ptb);

  free(pln_to_to_gn);
  free(distrib_init_elmt);
  free(pfield);
  free(dfield_strid);
  free(dfield);
  free(pstrid);

  PDM_MPI_Finalize ();
  return 0;
}
