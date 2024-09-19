
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_unique.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_deduce_combine_connectivity_impl
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_n,
        PDM_g_num_t   **dentity1_entity3
)
{
  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];
  int dn_entity2 = entity2_distrib[i_rank+1] - entity2_distrib[i_rank];

  /*
   * Some print to create unit_test
   */
  // PDM_log_trace_array_long(entity1_distrib  , n_rank+1               , "entity1_distrib::");
  // PDM_log_trace_array_long(entity2_distrib  , n_rank+1               , "entity2_distrib::");
  // PDM_log_trace_array_int (dentity1_entity2_idx, dn_entity1+1              , "dentity1_entity2_idx::");
  // PDM_log_trace_array_long(dentity1_entity2    , dentity1_entity2_idx[dn_entity1], "dentity1_entity2::");
  // PDM_log_trace_array_int (dentity2_entity3_idx , dn_entity2+1              , "dentity2_entity3_idx::");
  // PDM_log_trace_array_long(dentity2_entity3     , dentity2_entity3_idx[dn_entity2] , "dentity2_entity3::");

  int* dentity2_entity3_n = (int * ) malloc( dn_entity2 * sizeof(int));
  for(int i = 0; i < dn_entity2; ++i) {
    dentity2_entity3_n[i] = dentity2_entity3_idx[i+1] - dentity2_entity3_idx[i];
  }

  // PDM_log_trace_array_int(dentity1_entity2_idx, dn_entity1+1, "dentity1_entity2_idx::");
  // PDM_log_trace_array_long(dentity1_entity2_cur, dentity1_entity2_idx[dn_entity1], "dentity1_entity2::");

  /*
   * First compute the dentity1_entity3 connectivity
   *      -  We use ln_to_gn = dentity1_entity2 in order to have for each partition the entity2s
   *      -  We exhange the entity2_vtx, entity2_vtx_idx
   * So, for each entity2 describe by dentity1_entity2, we receive all vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity2_distrib,
                               (const PDM_g_num_t **) &dentity1_entity2,
                                                      &dentity1_entity2_idx[dn_entity1],
                                                      1,
                                                      comm);

  /*
   * Exchange
   */
  int**         pentity1_entity3_n;
  PDM_g_num_t** pentity1_entity3;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          dentity2_entity3_n,
             (void *  )   dentity2_entity3,
             (int  ***)  &pentity1_entity3_n,
             (void ***)  &pentity1_entity3);
  free(dentity2_entity3_n);

  if(is_signed) {
    int idx_read = 0;
    for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
      int sign = PDM_SIGN(dentity1_entity2[i]);
      for(int j = 0; j < pentity1_entity3_n[0][i]; ++j){
        pentity1_entity3[0][idx_read] = sign * pentity1_entity3[0][idx_read];
        idx_read++;
      }
    }
  }


  /*
   * Panic Verbose
   */
  // int* pentity1_entity3_idx = (int*) malloc( (dentity1_entity2_idx[dn_entity1] + 1) * sizeof(int));
  // pentity1_entity3_idx[0] = 0;
  // for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
  //   pentity1_entity3_idx[i+1] = pentity1_entity3_idx[i] + pentity1_entity3_n[0][i];
  // }

  // PDM_log_trace_array_int(pentity1_entity3_n[0], dentity1_entity2_idx[dn_entity1], "pentity1_entity3_n::");
  // PDM_log_trace_array_int(pentity1_entity3_idx, dentity1_entity2_idx[dn_entity1]+1, "pentity1_entity3_idx::");
  // PDM_log_trace_array_long(pentity1_entity3[0], pentity1_entity3_idx[dentity1_entity2_idx[dn_entity1]], "pentity1_entity3::");

  /*
   * Free
   */
  PDM_block_to_part_free(btp);

  /*
   * Assign pointer
   */
  *dentity1_entity3   = pentity1_entity3[0];
  *dentity1_entity3_n = pentity1_entity3_n[0];

  /*
   * Free first level of pointer - the second level is hold by dentity1_entity3/dentity1_entity3_n
   */
  free(pentity1_entity3  );
  free(pentity1_entity3_n);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compute the combine connectivty of entity1 with entity2 to entity3
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       Distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       Distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity2_entity3_idx
 * \param [in]   dentity2_entity3
 * \param [in]   is_signed             If connectivity is signed
 * \param [in]   dentity1_entity3_idx
 * \param [in]   dentity1_entity3
 */
void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
)
{
  int* pentity1_entity3_n;
  _deduce_combine_connectivity_impl(comm,
                                    entity1_distrib,
                                    entity2_distrib,
                                    dentity1_entity2_idx,
                                    dentity1_entity2,
                                    dentity2_entity3_idx,
                                    dentity2_entity3,
                                    is_signed,
                                    &pentity1_entity3_n,
                                    dentity1_entity3);

  PDM_g_num_t* _dentity1_entity3 = *dentity1_entity3;

  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];

  /*
   * Post-treat
   */
  *dentity1_entity3_idx = (int*) malloc( (dn_entity1 + 1) * sizeof(int));
  int* _dentity1_entity3_idx = *dentity1_entity3_idx;
  int*  dentity1_entity3_n   = (int*) malloc( (dn_entity1    ) * sizeof(int));

  int idx = 0;
  _dentity1_entity3_idx[0] = 0;
  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {

    dentity1_entity3_n   [i_entity1  ] = 0;
    _dentity1_entity3_idx[i_entity1+1] = _dentity1_entity3_idx[i_entity1];

    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    // printf("n_entity2_per_entity1::%i\n", n_entity2_per_entity1);

    for(int i_entity2 = 0; i_entity2 < n_entity2_per_entity1; ++i_entity2) {
      _dentity1_entity3_idx[i_entity1+1] += pentity1_entity3_n[idx];
      // if(is_signed) {
      //   int beg = dentity1_entity2_idx[i_entity1];
      //   int sgn = PDM_SIGN(dentity1_entity2[beg+i_entity2]);
      //   _dentity1_entity3[idx] = _dentity1_entity3[idx] * sgn;
      //   dentity1_entity3_n[i_entity1]      += pentity1_entity3_n[idx++];
      // } else {
        dentity1_entity3_n[i_entity1]      += pentity1_entity3_n[idx++];
      // }
    }
  }
  // printf("idx::%i\n", idx);
  // printf("dentity1_entity2_idx[dn_entity1]::%i\n", dentity1_entity2_idx[dn_entity1]);

  assert(idx == dentity1_entity2_idx[dn_entity1]);
  free(pentity1_entity3_n);

  // PDM_log_trace_array_int(_dentity1_entity3_idx, dn_entity1+1, "_dentity1_entity3_idx::");
  // PDM_log_trace_array_int(dentity1_entity3_n  , dn_entity1  , "dentity1_entity3_n::");

  PDM_para_graph_compress_connectivity(dn_entity1, _dentity1_entity3_idx, dentity1_entity3_n, _dentity1_entity3);

  /*
   * Realloc
   */
  *dentity1_entity3 = (PDM_g_num_t *) realloc(*dentity1_entity3, sizeof(PDM_g_num_t) * _dentity1_entity3_idx[dn_entity1] );
  _dentity1_entity3 = *dentity1_entity3;

  // PDM_log_trace_array_int (_dentity1_entity3_idx, dn_entity1+1              , "after -> _dentity1_entity3_idx::");
  // PDM_log_trace_array_int (dentity1_entity3_n   , dn_entity1                , "after -> dentity1_entity3_n::");
  // PDM_log_trace_array_long(_dentity1_entity3    , _dentity1_entity3_idx[dn_entity1], "after -> dentity1_entity3::");

  /*
   * Free
   */
  free(dentity1_entity3_n);
}


/**
 *
 * \brief Compute the combine connectivty of entity1 with entity2 to entity3
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       Distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       Distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity2_entity3_idx
 * \param [in]   dentity2_entity3
 * \param [in]   is_signed             If connectivity is signed
 * \param [in]   dentity1_entity3_idx
 * \param [in]   dentity1_entity3
 */
void
PDM_deduce_combine_connectivity_dual
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       PDM_g_num_t    **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
)
{
  int* pentity1_entity3_n;
  _deduce_combine_connectivity_impl(comm,
                                    entity1_distrib,
                                    entity2_distrib,
                                    dentity1_entity2_idx,
                                    dentity1_entity2,
                                    dentity2_entity3_idx,
                                    dentity2_entity3,
                                    is_signed,
                                    &pentity1_entity3_n,
                                    dentity1_entity3);

  PDM_g_num_t* _dentity1_entity3 = *dentity1_entity3;

  int i_rank, n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1 = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];

  /*
   * Post-treat
   */
  *dentity1_entity3_idx = (PDM_g_num_t*) malloc( (dn_entity1 + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dentity1_entity3_idx = *dentity1_entity3_idx;
  int*  dentity1_entity3_n   = (int*) malloc( (dn_entity1    ) * sizeof(int));

  int idx = 0;
  _dentity1_entity3_idx[0] = 0;
  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {

    dentity1_entity3_n   [i_entity1  ] = 0;
    _dentity1_entity3_idx[i_entity1+1] = _dentity1_entity3_idx[i_entity1];

    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    // printf("n_entity2_per_entity1::%i\n", n_entity2_per_entity1);

    for(int i_entity2 = 0; i_entity2 < n_entity2_per_entity1; ++i_entity2) {
      _dentity1_entity3_idx[i_entity1+1] += (PDM_g_num_t) pentity1_entity3_n[idx];
      dentity1_entity3_n[i_entity1]      += (PDM_g_num_t) pentity1_entity3_n[idx++];
    }
  }
  // printf("idx::%i\n", idx);
  // printf("dentity1_entity2_idx[dn_entity1]::%i\n", dentity1_entity2_idx[dn_entity1]);

  for(int i = 0; i < _dentity1_entity3_idx[dn_entity1]; ++i) {
    _dentity1_entity3[i] = PDM_ABS(_dentity1_entity3[i]);
  }

  assert(idx == dentity1_entity2_idx[dn_entity1]);
  free(pentity1_entity3_n);

  // PDM_log_trace_array_long(_dentity1_entity3_idx, dn_entity1+1, "_dentity1_entity3_idx::");
  // PDM_log_trace_array_int (dentity1_entity3_n   , dn_entity1  , "dentity1_entity3_n::");

  PDM_para_graph_compress_connectivity_dual(dn_entity1,
                                            entity1_distrib[i_rank],
                                            _dentity1_entity3_idx,
                                            dentity1_entity3_n,
                                            _dentity1_entity3);

  /*
   * Realloc
   */
  *dentity1_entity3 = (PDM_g_num_t *) realloc(*dentity1_entity3, sizeof(PDM_g_num_t) * _dentity1_entity3_idx[dn_entity1] );

  // PDM_log_trace_array_long(_dentity1_entity3_idx, dn_entity1+1              , "after -> _dentity1_entity3_idx::");
  // PDM_log_trace_array_int(dentity1_entity3_n   , dn_entity1                , "after -> dentity1_entity3_n::");
  // PDM_log_trace_array_long(dentity1_entity3[0]    , _dentity1_entity3_idx[dn_entity1], "after -> dentity1_entity3::");

  /*
   * Free
   */
  free(dentity1_entity3_n);
}

/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       Distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       Distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx  Index of dentitiy1->dentity2 connectivity
 * \param [in]   dentity1_entity2      Connectivity of dentitiy1->dentity2
 * \param [in]   is_signed             If connectivity is signed
 * \param [in]   dentity2_entity1_idx  Index of dentitiy2->dentity1 connectivity
 * \param [in]   dentity2_entity1      Connectivity of dentitiy2->dentity1
 */
void
PDM_dconnectivity_transpose
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
       PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
       int              is_signed,
       int            **dentity2_entity1_idx,
       PDM_g_num_t    **dentity2_entity1
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity1  = entity1_distrib[i_rank+1] - entity1_distrib[i_rank];
  // int dn_entity2  = entity2_distrib[i_rank+1] - entity2_distrib[i_rank];

  PDM_g_num_t* ln_to_gn = (PDM_g_num_t * ) dentity1_entity2;
  PDM_g_num_t* gnum     = (PDM_g_num_t * ) malloc( dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));

  PDM_g_num_t shift_g = 1 + entity1_distrib[i_rank]; // Entre 1 et N

  if(is_signed) {
    for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
      for(int j = dentity1_entity2_idx[i_entity]; j < dentity1_entity2_idx[i_entity+1]; ++j) {
        int g_sign = PDM_SIGN(dentity1_entity2[j]);
        gnum[j] = g_sign * (i_entity + shift_g);
      }
    }
  } else {
    ln_to_gn = (PDM_g_num_t * ) dentity1_entity2;
    for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
      for(int j = dentity1_entity2_idx[i_entity]; j < dentity1_entity2_idx[i_entity+1]; ++j) {
        gnum[j] = i_entity + shift_g;
      }
    }
  }

  // PDM_log_trace_array_long(ln_to_gn, dentity1_entity2_idx[dn_entity1], "ln_to_gn::");
  // PDM_log_trace_array_long(gnum,  dentity1_entity2_idx[dn_entity1], "gnum::");
  /*
   * In order to revert the conncectivty we use the global numbering property
   */
  PDM_part_to_block_t *ptb = NULL;

  int save_entity_distrib = 0;
  if(entity2_distrib != NULL && entity2_distrib[0] != -1) {
    ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                    &ln_to_gn,
                                    entity2_distrib,
                         (int *)    &dentity1_entity2_idx[dn_entity1],
                                     1,
                                     comm);
  } else {
    save_entity_distrib = 1;
    ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                   1.,
                                   &ln_to_gn,
                                   NULL,
                         (int *)   &dentity1_entity2_idx[dn_entity1],
                                    1,
                                    comm);
  }

  int* send_stri = PDM_array_const_int(dentity1_entity2_idx[dn_entity1], 1);

  int         *dentity2_entity1_n = NULL;
  PDM_g_num_t *recv_data          = NULL;

  int blk_size = PDM_part_to_block_exch (ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                         1,
                                         &send_stri,
                               (void **) &gnum,
                                         &dentity2_entity1_n,
                               (void **) &recv_data);
  PDM_UNUSED(blk_size);

  int dn_entity2_recv = PDM_part_to_block_n_elt_block_get(ptb);

  if(entity2_distrib != NULL && save_entity_distrib != 1) {
    PDM_g_num_t *distrib2_idx_full =
      PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                      &dentity2_entity1_n,
                                                      entity2_distrib[n_rank]);
    dn_entity2_recv = distrib2_idx_full[i_rank+1] - distrib2_idx_full[i_rank];
    free(distrib2_idx_full);
  }

  if(save_entity_distrib == 1) {
    // Update distrib
    PDM_g_num_t* ptb_distrib = PDM_part_to_block_distrib_index_get(ptb);
    for(int i = 0; i < n_rank+1; ++i) {
      entity2_distrib[i] = ptb_distrib[i];
    }
  }

  /*
   * Free
   */
  PDM_part_to_block_free(ptb);
  free(gnum);
  free(send_stri);

  /*
   * Allocate
   */
  *dentity2_entity1_idx = (int *) malloc( (dn_entity2_recv + 1) * sizeof(int));
  int* _dentity2_entity1_idx = *dentity2_entity1_idx;

  // printf("blk_size       ::%i\n", blk_size       );
  // printf("dn_entity2_recv::%i\n", dn_entity2_recv);
  // log_trace("blk_size = %i | dn_entity2_recv = %i \n", blk_size, dn_entity2_recv);

  // PDM_log_trace_array_long(dentity2_entity1_n, dn_entity2_recv, "Before : dentity2_entity1_n::");
  // PDM_log_trace_array_long(recv_data, blk_size, "Before : recv_data::");

  PDM_para_graph_compress_connectivity(dn_entity2_recv,
                                       _dentity2_entity1_idx,
                                       dentity2_entity1_n,
                                       recv_data);
  // printf("*dentity2_entity1_idx[dn_entity2_recv]       ::%i\n", _dentity2_entity1_idx[dn_entity2_recv]       );

  // *dentity2_entity1 = recv_data;

  /*
   * Realloc
   */
  *dentity2_entity1 = realloc(recv_data, _dentity2_entity1_idx[dn_entity2_recv] * sizeof(PDM_g_num_t));
  free(dentity2_entity1_n);
}

/**
 *
 * \brief Shortcut to PDM_dconnectivity_transpose for facecell like connectivity
 *
 * \param [in]   face_distri           distribution of faces over the procs (size=n_rank+1)
 * \param [in]   cell_distri           distribution of cells over the procs (size=n_rank+1)
 * \param [in]   dface_cell            Facecell array, including zeros (size=2*dn_face)
 * \param [out]  dcell_face_idx        Output cell_face_idx. Will be allocated (size=dn_cell+1)
 * \param [out]  dcell_face            Ouput cell_face. Will be allocated (size=cell_face_idx[dn_cell])
 * \param [in]   comm                  PDM_MPI communicator
 */
void PDM_dfacecell_to_dcellface
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const PDM_g_num_t* dface_cell,
  int**              dcell_face_idx,
  PDM_g_num_t**      dcell_face,
  PDM_MPI_Comm       comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int dn_face = face_distri[i_rank+1] - face_distri[i_rank];

  int*         _dface_cell_idx = (int *)         malloc((dn_face+1) * sizeof(int));
  PDM_g_num_t* _dface_cell     = (PDM_g_num_t *) malloc(2*dn_face   * sizeof(PDM_g_num_t));
  int w_idx = 0;
  _dface_cell_idx[0] = 0;
  for (int i_face = 0; i_face < dn_face; i_face++) {
    _dface_cell_idx[i_face+1] = _dface_cell_idx[i_face];
    if (dface_cell[2*i_face] != 0) {
      _dface_cell_idx[i_face+1]++;
      _dface_cell[w_idx++] = dface_cell[2*i_face];
    }
    if (dface_cell[2*i_face+1] != 0) {
      _dface_cell_idx[i_face+1]++;
      _dface_cell[w_idx++] = -dface_cell[2*i_face+1];
    }
  }

  /*for (int i_face = 0; i_face < dn_face; i_face++)*/
    /*log_trace("%d %d\n", dface_cell[2*i_face], dface_cell[2*i_face+1]);*/

  /*PDM_log_trace_array_int(_dface_cell_idx, dn_face+1, "dface_cell_idx");*/
  /*PDM_log_trace_array_long(_dface_cell, 2*dn_face, "dface_cell");*/
  
  PDM_dconnectivity_transpose(comm,
                              face_distri,
         (PDM_g_num_t *)      cell_distri,
                              _dface_cell_idx,
                              _dface_cell,
                              1,
                              dcell_face_idx,
                              dcell_face);
  /*int dn_cell = cell_distri[i_rank+1] - cell_distri[i_rank];*/
  /*PDM_log_trace_connectivity_long(dcell_face_idx[0], dcell_face[0], dn_cell, "dcell_face");*/

  free(_dface_cell_idx);
  free(_dface_cell);
}

/**
 *
 * \brief Shortcut to PDM_dconnectivity_transpose to obtain facecell like connectivity
 *
 * \param [in]   face_distri           distribution of faces over the procs (size=n_rank+1)
 * \param [in]   cell_distri           distribution of cells over the procs (size=n_rank+1)
 * \param [in]   dcell_face_idx        cell_face idx array (size = dn_cell)
 * \param [in]   dcell_face            cell_face array (size = cell_face_idx[dn_cell])
 * \param [out]  dface_cell            Output face_cell. Will be allocated (size = 2*dn_face)
 * \param [in]   comm                  PDM_MPI communicator
 */
void PDM_dcellface_to_dfacecell
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const int*         dcell_face_idx,
  const PDM_g_num_t* dcell_face,
  PDM_g_num_t**      dface_cell,
  PDM_MPI_Comm       comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int dn_face = (int) face_distri[i_rank+1] - face_distri[i_rank];

  int*         _dface_cell_idx = NULL;
  PDM_g_num_t* _dface_cell     = NULL;
  PDM_dconnectivity_transpose(comm,
                              cell_distri,
         (PDM_g_num_t *)      face_distri,
                              dcell_face_idx,
                              dcell_face,
                              1,
                             &_dface_cell_idx,
                             &_dface_cell);

  /*PDM_log_trace_connectivity_long(_dface_cell_idx, _dface_cell, dn_face, "dface_cell");*/

  *dface_cell = (PDM_g_num_t*) malloc(2*dn_face*sizeof(PDM_g_num_t));

  int r_idx = 0;
  for (int i = 0; i < dn_face; i++) {
    int n_conn = _dface_cell_idx[i+1] - _dface_cell_idx[i];
    if (n_conn == 2) { //Interior faces
      assert (PDM_SIGN(_dface_cell[r_idx]) * PDM_SIGN(_dface_cell[r_idx+1]) == -1);
      if (_dface_cell[r_idx] < 0) {
        (*dface_cell)[2*i]   = _dface_cell[r_idx+1];
        (*dface_cell)[2*i+1] = PDM_ABS(_dface_cell[r_idx]);
      }
      else {
        (*dface_cell)[2*i]   = _dface_cell[r_idx];
        (*dface_cell)[2*i+1] = PDM_ABS(_dface_cell[r_idx+1]);
      }
      r_idx += 2;
    }
    else { //Bnd faces
      assert (n_conn == 1);
      int sign = PDM_SIGN(_dface_cell[r_idx]);
      // PDM_g_num_t val = PDM_ABS(_dface_cell[r_idx]);
      if (sign > 0) {
        (*dface_cell)[2*i]   = _dface_cell[r_idx];
        (*dface_cell)[2*i+1] = 0;
      }
      else {
        (*dface_cell)[2*i]   = 0;
        (*dface_cell)[2*i+1] = _dface_cell[r_idx];
      }
      r_idx++;
    }
  }

  /*PDM_log_trace_array_long(*dface_cell, 2*dn_face, "dface_cell");*/
  free(_dface_cell_idx);
  free(_dface_cell);
}


/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       Distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2      Connectivity entity1->entity2
 * \param [in]   dentity2_entity1      Reversed connectivity of entity1->entity2
 */
void
PDM_dorder_reverse
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity_distrib,
 const PDM_g_num_t     *dentity1_entity2,
       PDM_g_num_t    **dentity2_entity1
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int dn_entity1  = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  PDM_g_num_t* gnum = (PDM_g_num_t * ) malloc( dn_entity1 * sizeof(PDM_g_num_t));

  for(int i_entity = 0; i_entity < dn_entity1; ++i_entity) {
    gnum[i_entity] = entity_distrib[i_rank] + i_entity + 1;
  }

  /*
   * In order to revert the conncectivty we use the global numbering property
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                      (PDM_g_num_t **) &dentity1_entity2,
                                                       entity_distrib,
                                            (int *)    &dn_entity1,
                                                       1,
                                                       comm);

  PDM_g_num_t *recv_data          = NULL;

  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &gnum,
                           NULL,
                 (void **) &recv_data);
  free(gnum);

  *dentity2_entity1 = recv_data;
  PDM_part_to_block_free(ptb);
}

void
PDM_dgroup_entity_transpose
(
 int            n_group,
 int           *dgroup_entity_idx,
 PDM_g_num_t   *dgroup_entity,
 PDM_g_num_t   *distrib_entity,
 int          **dentity_group_idx,
 int          **dentity_group,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  //printf(" dgroup_entity_idx[%i] = %i \n", n_group, dgroup_entity_idx[n_group]);
  PDM_part_to_block_t *ptb = PDM_part_to_block_create_from_distrib (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &dgroup_entity,
                                                        distrib_entity,
                                                        &dgroup_entity_idx[n_group],
                                                        1,
                                                        comm);


  int* pgroup_id_n = (int *) malloc(dgroup_entity_idx[n_group] * sizeof(int));
  int* pgroup_id   = (int *) malloc(dgroup_entity_idx[n_group] * sizeof(int));

  for(int i_group = 0; i_group < n_group; ++i_group) {
    for(int i = dgroup_entity_idx[i_group]; i < dgroup_entity_idx[i_group+1]; ++i){
      pgroup_id_n[i] = 1;
      pgroup_id  [i] = i_group;
    }
  }

  /*
   *  Exchange group id
   */

  int *tmp_dentity_group_n = NULL;
  int *tmp_dentity_group   = NULL;
  int s_block = PDM_part_to_block_exch (ptb,
                                        sizeof(int),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        1,
                                        &pgroup_id_n,
                              (void **) &pgroup_id,
                                        &tmp_dentity_group_n,
                              (void **) &tmp_dentity_group);
  free(pgroup_id_n);
  free(pgroup_id);


  int dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  PDM_g_num_t* tmp_distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb, &tmp_dentity_group_n, distrib_entity[n_rank]);
  free(tmp_distrib);
  PDM_part_to_block_free (ptb);

  /*
   * Post-treatment
   */
  *dentity_group     = malloc(s_block       * sizeof(int));
  *dentity_group_idx = malloc((dn_entity+1) * sizeof(int));
  int *_dentity_group     = *dentity_group;
  int *_dentity_group_idx = *dentity_group_idx;

  int idx_read  = 0;
  int idx_write = 0;
  _dentity_group_idx[0] = 0;
  for(int i = 0; i < dn_entity; ++i) {
    int n_id = tmp_dentity_group_n[i];
    int n_unique_id = 0;
    if(n_id > 0) {
      n_unique_id = PDM_inplace_unique(&tmp_dentity_group[idx_read], 0, n_id-1);
    }

    for(int j = 0; j < n_unique_id; ++j) {
      _dentity_group[idx_write++] = tmp_dentity_group[idx_read+j];
    }
    _dentity_group_idx[i+1] = _dentity_group_idx[i] + n_unique_id;

    idx_read += n_id;
  }

  free(tmp_dentity_group_n);
  free(tmp_dentity_group);

  if(0 == 1) {
    PDM_log_trace_connectivity_int(_dentity_group_idx, _dentity_group, dn_entity, "_dentity_group ::");
  }

  *dentity_group = realloc(*dentity_group, _dentity_group_idx[dn_entity] * sizeof(int));
  if (*dentity_group == NULL) {
    *dentity_group = malloc(sizeof(int) * _dentity_group_idx[dn_entity]);
    assert(*dentity_group != NULL);
  }
}



void
PDM_dentity_group_signed_transpose
(
 int            n_group,
 int           *dentity_group_idx,
 int           *dentity_group,
 PDM_g_num_t   *distrib_entity,
 int          **dgroup_entity_idx,
 PDM_g_num_t  **dgroup_entity,
 PDM_MPI_Comm   comm,
 const int      is_signed
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  int          *select_entity_n = malloc(n_group * sizeof(int          ));
  PDM_g_num_t **select_entity   = malloc(n_group * sizeof(PDM_g_num_t *));
  for(int i_group = 0; i_group < n_group; ++i_group) {
    select_entity_n[i_group] = 0;
  }

  for(int idx_entity = 0; idx_entity < dentity_group_idx[dn_entity]; ++idx_entity) {
    int i_group   = dentity_group[idx_entity];
      if (is_signed) {
        i_group = PDM_ABS(i_group) - 1;
      }
    select_entity_n[i_group]++;
  }

  if(0 == 1) {
    PDM_log_trace_array_int(select_entity_n, n_group, " select_entity_n : ");
  }

  for(int i_group = 0; i_group < n_group; ++i_group) {
    select_entity  [i_group] = malloc( select_entity_n[i_group] * sizeof(PDM_g_num_t));
    select_entity_n[i_group] = 0;
  }


  for(int i_entity = 0; i_entity < dn_entity; ++i_entity) {
    for(int idx_entity = dentity_group_idx[i_entity]; idx_entity <  dentity_group_idx[i_entity+1]; ++idx_entity) {
      int i_group   = dentity_group[idx_entity];
      if (is_signed) {
        i_group = PDM_ABS(i_group) - 1;
      }
      int idx_write = select_entity_n[i_group]++;
      select_entity[i_group][idx_write] = distrib_entity[i_rank] + i_entity + 1;
    }
  }

  if(0 == 1) {
    for(int i_group = 0; i_group < n_group; ++i_group) {
      PDM_log_trace_array_long(select_entity[i_group], select_entity_n[i_group], " select_entity : ");
    }
  }

  /*
   *  For each group we repart all information by block
   */
  PDM_part_to_block_t **ptb = malloc(n_group * sizeof(PDM_part_to_block_t *));

  int *_dgroup_entity_idx = (int * ) malloc( (n_group+1) * sizeof(int));
  _dgroup_entity_idx[0] = 0;
  for(int i_group = 0; i_group < n_group; ++i_group) {

    _dgroup_entity_idx[i_group+1] = _dgroup_entity_idx[i_group];
    ptb[i_group] = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                            PDM_PART_TO_BLOCK_POST_CLEANUP,
                                            1.,
                                            &select_entity[i_group],
                                            NULL,
                                            &select_entity_n[i_group],
                                            1,
                                            comm);
    _dgroup_entity_idx[i_group+1] += PDM_part_to_block_n_elt_block_get(ptb[i_group]);
  }

  // PDM_log_trace_array_long(_dgroup_entity_idx, n_group+1, " _dgroup_entity_idx : ");
  PDM_g_num_t *_dgroup_entity = malloc(_dgroup_entity_idx[n_group] * sizeof(PDM_g_num_t));
  for(int i_group = 0; i_group < n_group; ++i_group) {
    PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get(ptb[i_group]);


    // PDM_log_trace_array_long(blk_gnum, PDM_part_to_block_n_elt_block_get(ptb[i_group]), " blk_gnum : ");

    for(int i = _dgroup_entity_idx[i_group]; i < _dgroup_entity_idx[i_group+1]; ++i) {
      int j = i - _dgroup_entity_idx[i_group];
      _dgroup_entity[i] = blk_gnum[j];
    }
  }



  for(int i_group = 0; i_group < n_group; ++i_group) {
    free(select_entity[i_group]);
    PDM_part_to_block_free(ptb[i_group]);
  }
  free(ptb);
  free(select_entity_n);
  free(select_entity  );

  *dgroup_entity_idx = _dgroup_entity_idx;
  *dgroup_entity     = _dgroup_entity;

}



void
PDM_dentity_group_transpose
(
 int            n_group,
 int           *dentity_group_idx,
 int           *dentity_group,
 PDM_g_num_t   *distrib_entity,
 int          **dgroup_entity_idx,
 PDM_g_num_t  **dgroup_entity,
 PDM_MPI_Comm   comm
)
{
  PDM_dentity_group_signed_transpose(n_group,
                                     dentity_group_idx,
                                     dentity_group,
                                     distrib_entity,
                                     dgroup_entity_idx,
                                     dgroup_entity,
                                     comm,
                                     0);
}


void
PDM_dconnectivity_to_extract_dconnectivity_bis
(
 const PDM_MPI_Comm    comm,
       int             n_selected_entity1,
       PDM_g_num_t    *select_entity1,
       PDM_g_num_t    *entity1_distribution,
       int            *dentity1_entity2_idx,
       PDM_g_num_t    *dentity1_entity2,
       PDM_g_num_t   **extract_entity1_distribution,
       PDM_g_num_t   **extract_entity2_distribution,
       int           **dextract_entity1_entity2_idx,
       PDM_g_num_t   **dextract_entity1_entity2,
       PDM_g_num_t   **dparent_entity1_g_num,
       PDM_g_num_t   **dparent_entity2_g_num,
       PDM_g_num_t   **entity1_old_to_new
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Create implcite global numbering
   */
  int dn_entity1 = entity1_distribution[i_rank+1] - entity1_distribution[i_rank];

  /*
   *  Create global numbering from parent
   */
  PDM_gen_gnum_t* gen_gnum_entity1 = PDM_gnum_create(3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);


  PDM_gnum_set_from_parents (gen_gnum_entity1,
                             0,
                             n_selected_entity1,
                             select_entity1);

  PDM_gnum_compute (gen_gnum_entity1);


  PDM_g_num_t* extract_entity1_ln_to_gn = PDM_gnum_get (gen_gnum_entity1, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(select_entity1          , n_selected_entity1, "select_entity1:: ");
    PDM_log_trace_array_long(extract_entity1_ln_to_gn, n_selected_entity1, "extract_entity1_ln_to_gn:: ");
  }

  PDM_gnum_free(gen_gnum_entity1);

  /*
   * Caution we need a result independant of parallelism
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity1_distribution,
                               (const PDM_g_num_t **) &select_entity1,
                                                      &n_selected_entity1,
                                                      1,
                                                      comm);

  /*
   * Prepare data
   */
  int* dentity1_entity2_n = (int *) malloc( sizeof(int) * dn_entity1);
  for(int i_elmt = 0; i_elmt < dn_entity1; ++i_elmt){
    dentity1_entity2_n[i_elmt] = dentity1_entity2_idx[i_elmt+1] - dentity1_entity2_idx[i_elmt];
  }

  /*
   * Exchange
   */
  int**         tmp_pextract_entity1_entity2_n;
  PDM_g_num_t** tmp_pextract_entity1_entity2;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          dentity1_entity2_n,
             (void *  )   dentity1_entity2,
             (int  ***)  &tmp_pextract_entity1_entity2_n,
             (void ***)  &tmp_pextract_entity1_entity2);

  int*         pextract_entity1_entity2_n = tmp_pextract_entity1_entity2_n[0];
  PDM_g_num_t* pextract_entity1_entity2   = tmp_pextract_entity1_entity2[0];
  free(tmp_pextract_entity1_entity2_n);
  free(tmp_pextract_entity1_entity2);
  free(dentity1_entity2_n);

  PDM_block_to_part_free(btp);

  /*
   *  Remap inside true block to ensure parallelism independant
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      &extract_entity1_ln_to_gn,
                                                      NULL,
                                                      &n_selected_entity1,
                                                      1,
                                                      comm);

  int*         _dextract_entity1_entity2_n = NULL;
  PDM_g_num_t* _dextract_entity1_entity2   = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &pextract_entity1_entity2_n,
                (void **) &pextract_entity1_entity2,
                          &_dextract_entity1_entity2_n,
                (void **) &_dextract_entity1_entity2);

  PDM_g_num_t* _dparent_entity1_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &select_entity1,
                          NULL,
                (void **) &_dparent_entity1_g_num);

  free(pextract_entity1_entity2_n);
  free(pextract_entity1_entity2  );

  /*
   * Post-Treatment
   */
  PDM_g_num_t* distrib_extract_entity1 = PDM_part_to_block_distrib_index_get(ptb);
  int dn_extract_entity1 = distrib_extract_entity1[i_rank+1] - distrib_extract_entity1[i_rank];

  int* _dextract_entity1_entity2_idx = malloc( (dn_extract_entity1 + 1) * sizeof(int));
  _dextract_entity1_entity2_idx[0] = 0;
  for(int i = 0; i < dn_extract_entity1; ++i) {
    _dextract_entity1_entity2_idx[i+1] = _dextract_entity1_entity2_idx[i] + _dextract_entity1_entity2_n[i];
  }
  free(_dextract_entity1_entity2_n);

  /*
   * Revert information
   */
  PDM_g_num_t* entity1_init_distrib = PDM_compute_entity_distribution(comm, n_selected_entity1);

  PDM_dorder_reverse(comm,
                     entity1_init_distrib,
                     extract_entity1_ln_to_gn,
                     entity1_old_to_new);

  if(0 == 1) {
    PDM_log_trace_array_long(extract_entity1_ln_to_gn  , n_selected_entity1, "extract_entity1_ln_to_gn:: ");
    PDM_log_trace_array_long(*entity1_old_to_new, n_selected_entity1, "entity1_old_to_new:: ");
  }
  // free(old_to_new_entity_ln_to_gn);
  free(entity1_init_distrib);


  /*
   *  Create global numbering from parent entity2
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);


  PDM_gnum_set_from_parents (gen_gnum_entity2,
                             0,
                             _dextract_entity1_entity2_idx[dn_extract_entity1],
                             _dextract_entity1_entity2);

  PDM_gnum_compute (gen_gnum_entity2);


  PDM_g_num_t* extract_entity2_ln_to_gn = PDM_gnum_get (gen_gnum_entity2, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(_dextract_entity1_entity2, _dextract_entity1_entity2_idx[dn_extract_entity1], "select_entity2:: ");
    PDM_log_trace_array_long(extract_entity2_ln_to_gn, _dextract_entity1_entity2_idx[dn_extract_entity1], "extract_entity2_ln_to_gn:: ");
  }


  PDM_gnum_free(gen_gnum_entity2);


  PDM_part_to_block_t *ptb_entity2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              &extract_entity2_ln_to_gn,
                                                              NULL,
                                                              &_dextract_entity1_entity2_idx[dn_extract_entity1],
                                                              1,
                                                              comm);


  PDM_g_num_t* _dparent_entity2_g_num = NULL;
  PDM_part_to_block_exch (ptb_entity2,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &_dextract_entity1_entity2,
                          NULL,
                (void **) &_dparent_entity2_g_num);


  PDM_g_num_t* distrib_extract_entity2 = PDM_part_to_block_distrib_index_get(ptb_entity2);
  int dn_extract_entity2 = distrib_extract_entity2[i_rank+1] - distrib_extract_entity2[i_rank];

  if(0 == 1) {
    PDM_log_trace_array_long(_dparent_entity2_g_num, dn_extract_entity2, "_dparent_entity2_g_num:: ");
  }

  /*
   * Update connectivity with the new ones
   */
  for(int i = 0; i < _dextract_entity1_entity2_idx[dn_extract_entity1]; ++i) {
    _dextract_entity1_entity2[i] = extract_entity2_ln_to_gn[i];
  }

  free(extract_entity1_ln_to_gn);
  free(extract_entity2_ln_to_gn);
  *dextract_entity1_entity2_idx = _dextract_entity1_entity2_idx;
  *dextract_entity1_entity2     = _dextract_entity1_entity2;
  *dparent_entity1_g_num        = _dparent_entity1_g_num;
  *dparent_entity2_g_num        = _dparent_entity2_g_num;

  *extract_entity1_distribution = malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  *extract_entity2_distribution = malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t *_extract_entity1_distribution = *extract_entity1_distribution;
  PDM_g_num_t *_extract_entity2_distribution = *extract_entity2_distribution;

  for(int i = 0; i < n_rank+1; ++i) {
    _extract_entity1_distribution[i] = distrib_extract_entity1[i];
    _extract_entity2_distribution[i] = distrib_extract_entity2[i];
  }

  PDM_part_to_block_free(ptb);
  PDM_part_to_block_free(ptb_entity2);
}


void
PDM_dconnectivity_to_extract_dconnectivity_block
(
 const PDM_MPI_Comm          comm,
       int                   dn_extract_entity1,
       PDM_g_num_t          *dextract_gnum_entity1,
       PDM_g_num_t          *entity1_distribution,
       int                  *dentity1_entity2_idx,
       PDM_g_num_t          *dentity1_entity2,
       int                 **dextract_entity1_entity2_idx,
       PDM_g_num_t         **dextract_entity1_entity2,
       PDM_block_to_part_t **btp_entity1_to_extract_entity1,
       PDM_g_num_t         **extract_entity2_distribution,
       PDM_g_num_t         **dparent_entity2_g_num
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   * Caution we need a result independant of parallelism
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity1_distribution,
                               (const PDM_g_num_t **) &dextract_gnum_entity1,
                                                      &dn_extract_entity1,
                                                      1,
                                                      comm);
  *btp_entity1_to_extract_entity1 = btp;

  /*
   * Prepare data
   */
  int dn_entity1 = entity1_distribution[i_rank+1] - entity1_distribution[i_rank];
  int* dentity1_entity2_n = (int *) malloc( sizeof(int) * dn_entity1);
  for(int i_elmt = 0; i_elmt < dn_entity1; ++i_elmt){
    dentity1_entity2_n[i_elmt] = dentity1_entity2_idx[i_elmt+1] - dentity1_entity2_idx[i_elmt];
  }

  /*
   * Exchange
   */
  int**         tmp_pextract_entity1_entity2_n = NULL;
  PDM_g_num_t** tmp_pextract_entity1_entity2   = NULL;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          dentity1_entity2_n,
             (void *  )   dentity1_entity2,
             (int  ***)  &tmp_pextract_entity1_entity2_n,
             (void ***)  &tmp_pextract_entity1_entity2);

  int*         pextract_entity1_entity2_n = tmp_pextract_entity1_entity2_n[0];
  PDM_g_num_t* pextract_entity1_entity2   = tmp_pextract_entity1_entity2  [0];
  free(tmp_pextract_entity1_entity2_n);
  free(tmp_pextract_entity1_entity2);
  free(dentity1_entity2_n);

  int* _dextract_entity1_entity2_idx = malloc( (dn_extract_entity1 + 1) * sizeof(int));
  _dextract_entity1_entity2_idx[0] = 0;
  for(int i = 0; i < dn_extract_entity1; ++i) {
    _dextract_entity1_entity2_idx[i+1] = _dextract_entity1_entity2_idx[i] + pextract_entity1_entity2_n[i];
  }
  free(pextract_entity1_entity2_n);

  /*
   * Prepare next renumbering
   */
  PDM_part_to_block_t *ptb_entity2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              &pextract_entity1_entity2,
                                                              NULL,
                                                              &_dextract_entity1_entity2_idx[dn_extract_entity1],
                                                              1,
                                                              comm);

  int          dn_extract_entity2      = PDM_part_to_block_n_elt_block_get  (ptb_entity2);
  PDM_g_num_t *dextract_gnum_entity2   = PDM_part_to_block_block_gnum_get   (ptb_entity2);

  *extract_entity2_distribution = PDM_compute_entity_distribution(comm, dn_extract_entity2);

  *dparent_entity2_g_num = malloc(dn_extract_entity2 * sizeof(PDM_g_num_t));
  PDM_g_num_t *_dparent_entity2_g_num        = *dparent_entity2_g_num;
  for(int i = 0; i < dn_extract_entity2; ++i) {
    _dparent_entity2_g_num[i] = dextract_gnum_entity2[i];
  }
  PDM_part_to_block_free(ptb_entity2);


  /*
   * We use gen gnum to update numbering of descending connectivity
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);
  // PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 1);

  int *dextract_entity1_entity2_sgn = malloc(_dextract_entity1_entity2_idx[dn_extract_entity1] * sizeof(int));
  for(int i = 0; i < _dextract_entity1_entity2_idx[dn_extract_entity1]; ++i) {
    dextract_entity1_entity2_sgn[i] = PDM_SIGN(pextract_entity1_entity2[i]);
    pextract_entity1_entity2    [i] = PDM_ABS (pextract_entity1_entity2[i]);
  }

  PDM_gnum_set_from_parents (gen_gnum_entity2,
                             0,
                             _dextract_entity1_entity2_idx[dn_extract_entity1],
                             pextract_entity1_entity2);
  PDM_gnum_compute (gen_gnum_entity2);
  PDM_g_num_t* _dextract_entity1_entity2 = PDM_gnum_get(gen_gnum_entity2, 0);
  free(pextract_entity1_entity2);

  for(int i = 0; i < _dextract_entity1_entity2_idx[dn_extract_entity1]; ++i) {
    _dextract_entity1_entity2[i] = _dextract_entity1_entity2[i] * dextract_entity1_entity2_sgn[i];
  }
  free(dextract_entity1_entity2_sgn);


  PDM_gnum_free(gen_gnum_entity2);

  *dextract_entity1_entity2_idx = _dextract_entity1_entity2_idx;
  *dextract_entity1_entity2     = _dextract_entity1_entity2;
}


void
PDM_dconnectivity_to_extract_dconnectivity
(
 const PDM_MPI_Comm          comm,
       int                   n_selected_entity1,
       PDM_g_num_t          *select_entity1,
       PDM_g_num_t          *entity1_distribution,
       int                  *dentity1_entity2_idx,
       PDM_g_num_t          *dentity1_entity2,
       PDM_g_num_t         **extract_entity1_distribution,
       PDM_g_num_t         **dparent_entity1_g_num,
       int                 **dextract_entity1_entity2_idx,
       PDM_g_num_t         **dextract_entity1_entity2,
       PDM_block_to_part_t **btp_entity1_to_extract_entity1,
       PDM_g_num_t         **extract_entity2_distribution,
       PDM_g_num_t         **dparent_entity2_g_num
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   * Create implicit global numbering
   */
  double* weight = malloc(n_selected_entity1 * sizeof(double));
  for(int i = 0; i < n_selected_entity1; ++i) {
    weight[i] = 1.;
  }
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      &select_entity1,
                                                      &weight,
                                                      &n_selected_entity1,
                                                      1,
                                                      comm);
  free(weight);

  int          dn_extract_entity1      = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t *dextract_gnum_entity1   = PDM_part_to_block_block_gnum_get   (ptb);

  *extract_entity1_distribution = PDM_compute_entity_distribution(comm, dn_extract_entity1);

  *dparent_entity1_g_num        = malloc(dn_extract_entity1 * sizeof(PDM_g_num_t));
  PDM_g_num_t *_dparent_entity1_g_num        = *dparent_entity1_g_num;
  for(int i = 0; i < dn_extract_entity1; ++i) {
    _dparent_entity1_g_num[i] = dextract_gnum_entity1[i];
  }

  PDM_dconnectivity_to_extract_dconnectivity_block(comm,
                                                   dn_extract_entity1,
                                                   dextract_gnum_entity1,
                                                   entity1_distribution,
                                                   dentity1_entity2_idx,
                                                   dentity1_entity2,
                                                   dextract_entity1_entity2_idx,
                                                   dextract_entity1_entity2,
                                                   btp_entity1_to_extract_entity1,
                                                   extract_entity2_distribution,
                                                   dparent_entity2_g_num);

  PDM_part_to_block_free(ptb);
}



void
PDM_dconnectivity_dface_vtx_from_face_and_edge
(
 const PDM_MPI_Comm    comm,
       PDM_g_num_t    *distrib_face,
       PDM_g_num_t    *distrib_edge,
       int            *dface_edge_idx,
       PDM_g_num_t    *dface_edge,
       PDM_g_num_t    *dedge_vtx,
       PDM_g_num_t   **dface_vtx
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int dn_edge = distrib_edge[i_rank+1] - distrib_edge[i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];
  int* dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

  int pn_edge = dface_edge_idx[dn_face];
  PDM_g_num_t* edge_ln_to_gn = dface_edge;
  int pn_vtx = 0;
  PDM_g_num_t* pvtx_ln_to_gn = NULL;
  int* pedge_vtx_idx = NULL;
  int* pedge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           pn_edge,
                                                           edge_ln_to_gn,
                                                          &pn_vtx,
                                                          &pvtx_ln_to_gn,
                                                          &pedge_vtx_idx,
                                                          &pedge_vtx);
  free(dedge_vtx_idx);



  int *pface_edge = (int *) malloc(pn_edge*sizeof(int));
  for (int i=0; i < pn_edge; ++i) {
    pface_edge[i] = PDM_SIGN(dface_edge[i])*(i+1);
  }
  int *pface_vtx = NULL;
  PDM_compute_face_vtx_from_face_and_edge(dn_face,
                                          dface_edge_idx,
                                          pface_edge,
                                          pedge_vtx,
                                          &pface_vtx);

  *dface_vtx = (PDM_g_num_t*) malloc(dface_edge_idx[dn_face]*sizeof(PDM_g_num_t));
  for (int i=0; i < dface_edge_idx[dn_face]; ++i) {
    (*dface_vtx)[i] = pvtx_ln_to_gn[pface_vtx[i]-1];
  }
   
  free(pface_edge);
  free(pface_vtx);

  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
}


#ifdef  __cplusplus
}
#endif
