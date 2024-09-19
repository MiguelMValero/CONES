
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_multi_block_merge.h"
#include "pdm_multi_block_merge_priv.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_array.h"
#include "pdm_part_geom.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"

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

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

PDM_multi_block_merge_t*
PDM_multi_block_merge_create
(
       PDM_g_num_t  **block_distrib_idx,
 const int            n_block,
       int           *n_selected,
       PDM_g_num_t  **selected_g_num,
       int            graph_size,
       int           *dmerge_idx,
       int           *dmerge_block_id,
       PDM_g_num_t   *dmerge_g_num,
       PDM_MPI_Comm   comm
)
{
  PDM_multi_block_merge_t *mbm = (PDM_multi_block_merge_t *) malloc (sizeof(PDM_multi_block_merge_t));

  PDM_MPI_Comm_size (comm, &mbm->n_rank);
  PDM_MPI_Comm_rank (comm, &mbm->i_rank);

  mbm->comm        = comm;
  mbm->n_block     = n_block;

  mbm->blocks_ini_dn = (int *) malloc(n_block*sizeof(int));
  mbm->multi_block_distrib = malloc( (n_block + 1) * sizeof(PDM_g_num_t));
  mbm->multi_block_distrib[0] = 0;
  for(int i_block = 0; i_block < n_block; ++i_block) {
    mbm->multi_block_distrib[i_block+1] = mbm->multi_block_distrib[i_block] + block_distrib_idx[i_block][mbm->n_rank];
    mbm->blocks_ini_dn[i_block] = block_distrib_idx[i_block][mbm->i_rank+1] - block_distrib_idx[i_block][mbm->i_rank];
  }

  /*
   * Shift des données d'entrés
   */
  PDM_g_num_t **_selected_g_num  = malloc( (n_block + 1) * sizeof(PDM_g_num_t *));
  PDM_g_num_t **_send_orig_g_num = malloc( (n_block + 1) * sizeof(PDM_g_num_t *));
  int         **_send_stride     = malloc( (n_block + 1) * sizeof(int         *));
  int          *_n_selected      = malloc( (n_block + 1) * sizeof(int          ));
  for(int i_block = 0; i_block < n_block; ++i_block) {

    _n_selected     [i_block] = n_selected[i_block];
    _selected_g_num [i_block] = malloc(n_selected[i_block] * sizeof(PDM_g_num_t));
    _send_stride    [i_block] = malloc(n_selected[i_block] * sizeof(int        ));
    _send_orig_g_num[i_block] = _selected_g_num[i_block];

    for(int i = 0; i < n_selected[i_block]; ++i) {
      _selected_g_num[i_block][i] = mbm->multi_block_distrib[i_block] + selected_g_num[i_block][i];
      _send_stride   [i_block][i] = 1;
    }
  }

  // Add graph as partition
  _n_selected     [n_block] = dmerge_idx[graph_size];
  _selected_g_num [n_block] = malloc(dmerge_idx[graph_size] * sizeof(PDM_g_num_t));
  _send_orig_g_num[n_block] = malloc(2*(dmerge_idx[graph_size]-graph_size) * sizeof(PDM_g_num_t));
  _send_stride    [n_block] = malloc(dmerge_idx[graph_size] * sizeof(int        ));
  int w_idx = 0;
  for(int i = 0; i < graph_size; ++i) {
    // First vertex in the graph becomes the reference. We send to corresponding gnum
    //   the list of related gnums; and we send to the other vertices 0, which will be used
    //   to indicate that vertex must be removed

    int first = dmerge_idx[i];
    _send_stride[n_block][first] = dmerge_idx[i+1] - dmerge_idx[i] - 1;
    _selected_g_num[n_block][first] = mbm->multi_block_distrib[dmerge_block_id[first]] + dmerge_g_num[first];
    for(int j = first+1; j < dmerge_idx[i+1]; ++j) { //First (exlude itself ; it will be send by first part)
      _send_orig_g_num[n_block][w_idx++] = mbm->multi_block_distrib[dmerge_block_id[j]] + dmerge_g_num[j];
    }

    //Others : we send 0 to warn that entity will be removed
    for(int j = dmerge_idx[i]+1; j < dmerge_idx[i+1]; ++j) {
      _selected_g_num[n_block][j] = mbm->multi_block_distrib[dmerge_block_id[j]] + dmerge_g_num[j];
      _send_stride[n_block][j] = 1;
      _send_orig_g_num[n_block][w_idx++] = 0;
    }
  }
  assert (w_idx == 2*(dmerge_idx[graph_size] - graph_size));

  if (0 == 1) {
    PDM_log_trace_array_int(_n_selected, n_block+1, "n selected");
    for (int i_block = 0; i_block < n_block+1; i_block++) {
      log_trace("block %d\n", i_block);
      PDM_log_trace_array_long(_selected_g_num[i_block], _n_selected[i_block] ,"  selected gnum");
      PDM_log_trace_array_long(_send_orig_g_num[i_block], w_idx ,"  send gnum");
    }
  }

  /*
   * Use ptb to equilibrate AND to reforms a block and keep link with the concatenate global numbering
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      _selected_g_num,
                                                      NULL,
                                                      _n_selected,
                                                      n_block+1,
                                                      comm);

  int *recv_stride = NULL;
  PDM_g_num_t *blk_send_orig_g_num = NULL;
  int dn_blk_size = PDM_part_to_block_exch(ptb,
                                           sizeof(PDM_g_num_t),
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           _send_stride,
                                 (void **) _send_orig_g_num,
                                           &recv_stride,
                                 (void **) &blk_send_orig_g_num);

  int dn_merge  = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t* merge_gnum = PDM_part_to_block_block_gnum_get(ptb);

  for(int i_block = 0; i_block < n_block+1; ++i_block) {
    free(_send_stride [i_block]);
    free(_selected_g_num[i_block]);
    if (i_block == n_block)
      free(_send_orig_g_num[i_block]); //Only last was allocated
  }
  free(_send_stride);
  free(_selected_g_num);
  free(_send_orig_g_num);
  free(_n_selected);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_stride, dn_merge, "recv_n :: ");
    PDM_log_trace_array_long(blk_send_orig_g_num, dn_blk_size, "blk_send_orig_g_num :: ");
  }

  // Now we have to build the new global numbering and old to new array : 
  // for numbering, we just select gnums who did not receive a 0
  // for old to new, we will just invert the recv field which is in a fact the new to old order
  int* recv_idx = PDM_array_new_idx_from_sizes_int (recv_stride, dn_merge);
  int* dnew_to_old_idx = (int *) malloc( (dn_merge+1) * sizeof(int));
  dnew_to_old_idx[0] = 0;
  mbm->dn_new_block = 0;

  for(int i = 0; i < dn_merge; ++i ) {
    int dn_loc = recv_idx[i+1] - recv_idx[i];
    int to_remove = 0;
    for (int j = recv_idx[i]; j < recv_idx[i+1]; j++) {
      if (blk_send_orig_g_num[j] == 0) {
        to_remove = 1;
        break;
      }
    }
    if (to_remove == 0) {
      dnew_to_old_idx[mbm->dn_new_block+1] = dnew_to_old_idx[mbm->dn_new_block] + dn_loc;

      //Erase recv data instead (let stuff at end; array is already allocated. 
      //Be careful not to free PTB until merge_gnum is useless !
      for (int k = 0; k < dn_loc; k++) {
        blk_send_orig_g_num[dnew_to_old_idx[mbm->dn_new_block] + k] = blk_send_orig_g_num[recv_idx[i]+k];
      }
      merge_gnum[mbm->dn_new_block++] = merge_gnum[i];
    }
  }

  if (0 == 1) {
    PDM_log_trace_array_int(dnew_to_old_idx, mbm->dn_new_block+1, "dnew_to_old_idx");
    PDM_log_trace_array_long(blk_send_orig_g_num, dnew_to_old_idx[mbm->dn_new_block], "dnew_to_old_idx");
  }

  mbm->distrib_merge = PDM_compute_entity_distribution(comm, mbm->dn_new_block);
  mbm->old_distrib   = PDM_compute_uniform_entity_distribution(mbm->comm, mbm->multi_block_distrib[n_block]);

  PDM_dconnectivity_transpose(comm,
                              mbm->distrib_merge,
                              mbm->old_distrib,
                              dnew_to_old_idx,
                              blk_send_orig_g_num,
                              0,
                              &mbm->dold_to_new_idx,
                              &mbm->dold_to_new);

  if(0 == 1) {
    int dn_orig = mbm->old_distrib[mbm->i_rank+1] - mbm->old_distrib[mbm->i_rank];
    PDM_log_trace_array_int (mbm->dold_to_new_idx, dn_orig+1                    , "dold_to_new_idx :: ");
    PDM_log_trace_array_long(mbm->dold_to_new    , mbm->dold_to_new_idx[dn_orig], "dold_to_new     :: ");
    PDM_log_trace_array_long(merge_gnum   , mbm->dn_new_block    , "merge_gnum :: ");
  }

  free(recv_stride);
  free(recv_idx);
  free(dnew_to_old_idx);
  free(blk_send_orig_g_num);

  mbm->mbtp = PDM_multi_block_to_part_create(mbm->multi_block_distrib,
                                             mbm->n_block,
                       (const PDM_g_num_t**) block_distrib_idx,
                       (const PDM_g_num_t**)&merge_gnum,
                                            &mbm->dn_new_block,
                                             1,
                                             mbm->comm);
  PDM_part_to_block_free(ptb); //This free merge_gnum

  return mbm;
}

/*
 *
 */
void
PDM_multi_block_merge_exch
(
 PDM_multi_block_merge_t     *mbm,
 size_t                       s_data,
 PDM_stride_t                 t_stride,
 int                        **block_stride,
 void                       **block_data,
 int                        **merge_block_stride,
 void                       **merge_block_data
)
{
  int           **tmp_block_strid_out = NULL;
  unsigned char **tmp_block_data_out  = NULL;
  PDM_multi_block_to_part_exch2(mbm->mbtp,
                                s_data,
                                t_stride,
                                block_stride,
                      (void **) block_data,
                                &tmp_block_strid_out,
                    (void ***)  &tmp_block_data_out);

  if(t_stride == PDM_STRIDE_VAR_INTERLACED) {
    *merge_block_stride = tmp_block_strid_out[0];
    free(tmp_block_strid_out);
  }

  *merge_block_data = tmp_block_data_out[0];
  free(tmp_block_data_out);

}


void
PDM_multi_block_merge_get_old_to_new
(
 PDM_multi_block_merge_t  *mbm,
 PDM_g_num_t             **old_distrib,
 int                     **dold_to_new_idx,
 PDM_g_num_t             **dold_to_new
)
{
  *old_distrib     = mbm->old_distrib;
  *dold_to_new_idx = mbm->dold_to_new_idx;
  *dold_to_new     = mbm->dold_to_new;
}


// Update a distributed field using a distributed ordering
static void _dist_data_update
(
 const PDM_g_num_t   *old_to_new_distri,
 const int           *dold_to_new_idx,
 const PDM_g_num_t   *dold_to_new,
 const int            dn_data,
 int                 *stride,
 PDM_g_num_t        **data,
 PDM_MPI_Comm         comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  PDM_g_num_t *_data = *data;
  int n_to_update = 0; //This is the number of ids to update
  for (int i = 0; i < dn_data; i++)
    n_to_update += stride[i];

  int *data_sign = (int *) malloc (n_to_update * sizeof(int));
  for (int i = 0; i < n_to_update; i++) {
    data_sign[i] = PDM_SIGN(_data[i]);
    _data[i]      = PDM_ABS (_data[i]);
  }

  // Prepare stride for old_to_new to send. Should be 0 or 1 !
  int dn_old_to_new = old_to_new_distri[i_rank+1] - old_to_new_distri[i_rank];
  int *dold_to_new_n = (int * ) malloc( dn_old_to_new * sizeof(int));
  for(int i = 0; i < dn_old_to_new; ++i) {
    dold_to_new_n[i] = dold_to_new_idx[i+1] - dold_to_new_idx[i];
    assert(dold_to_new_n[i] <= 1);
  }

  PDM_block_to_part_t *btp_update = PDM_block_to_part_create (old_to_new_distri,
                                              (const PDM_g_num_t **) &_data,
                                                                     &n_to_update,
                                                                      1,
                                                                      comm);

  int *recv_stride = (int *) malloc(n_to_update*sizeof(int));
  PDM_block_to_part_exch_in_place (btp_update,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          dold_to_new_n,
                 (void *) dold_to_new,
                          &recv_stride,
               (void **)  &_data); //Recycle memory

  free(dold_to_new_n);
  PDM_block_to_part_free(btp_update);

  //Some ids may have been removed, so we have to update stride and to be carefull
  //when applying sign
  int r_idx = 0;
  int w_idx = 0;
  for (int i = 0; i < dn_data; i++) {
    int new_stride = 0;
    for (int j = 0; j < stride[i]; j++) {
      if (recv_stride[r_idx] > 0) {
        _data[w_idx++] *= data_sign[r_idx];
        new_stride ++;
      }
      r_idx++;
    }
    stride[i] = new_stride; //Update stride
  }
  assert (r_idx == n_to_update);

  //Should we realloc ?
  PDM_g_num_t *_data_realloc = realloc(_data, w_idx*sizeof(PDM_g_num_t));
  *data = _data_realloc;

  free(recv_stride);
  free(data_sign);
}

void
PDM_multi_block_merge_exch_and_update
(
 PDM_multi_block_merge_t  *mbm_cur,
 PDM_multi_block_merge_t  *mbm_for_update,
 PDM_stride_t              t_stride,
 int                     **block_stride,
 PDM_g_num_t             **block_data,
 int                     **update_domain,
 int                     **merge_block_stride,
 PDM_g_num_t             **merge_block_data
)
{
  
  // 1. Shift the connectivity to exchange to make it related to its domain
  int n_block = mbm_cur->n_block;
  PDM_g_num_t **shifted_block_data = (PDM_g_num_t **) malloc(n_block*sizeof(PDM_g_num_t*));
  for (int i_block = 0; i_block < n_block; i_block++) {
    int n_data_block = 0;
    if (t_stride == PDM_STRIDE_CST_INTERLACED) {
      n_data_block = mbm_cur->blocks_ini_dn[i_block] * block_stride[0][0];
    }
    else {
      for (int i = 0; i < mbm_cur->blocks_ini_dn[i_block]; i++) {
        n_data_block += block_stride[i_block][i];
      }
    }
    shifted_block_data[i_block] = (PDM_g_num_t *) malloc(n_data_block*sizeof(PDM_g_num_t));
    if (update_domain == NULL) {
      PDM_g_num_t shift = mbm_for_update->multi_block_distrib[i_block];
      for (int i = 0; i < n_data_block; i++) {
        shifted_block_data[i_block][i] = block_data[i_block][i] + shift;
      }
    }
    else {
      for (int i = 0; i < n_data_block; i++) {
        PDM_g_num_t shift = mbm_for_update->multi_block_distrib[update_domain[i_block][i]];
        shifted_block_data[i_block][i] = block_data[i_block][i] + shift;
      }
    }
  }
   

  // 2. Merge data using main mbm; 
  int         *dparent_child_strid = NULL;
  PDM_g_num_t *dparent_child_mbm = NULL;

  PDM_multi_block_merge_exch(mbm_cur,
                             sizeof(PDM_g_num_t),
                             t_stride,
                             block_stride,
                 (void **)   shifted_block_data,
                             &dparent_child_strid,
                 (void **)   &dparent_child_mbm);

  // 3. Update received (dist) connectivity to have ids in the new numbering
  int dn_merge = PDM_multi_block_merge_get_n_block(mbm_cur);
  if (t_stride == PDM_STRIDE_CST_INTERLACED) { //Next function does not allow cst stride
    int cst_stride = block_stride[0][0];
    dparent_child_strid = PDM_array_const_int(dn_merge, cst_stride);
  }

  // Get old to new ordering for entity to update
  PDM_g_num_t *old_distrib;
  int         *dold_to_new_idx;
  PDM_g_num_t *dold_to_new;
  PDM_multi_block_merge_get_old_to_new(mbm_for_update, &old_distrib, &dold_to_new_idx, &dold_to_new);

  _dist_data_update(old_distrib,
                    dold_to_new_idx,
                    dold_to_new,
                    dn_merge,
                    dparent_child_strid,
                   &dparent_child_mbm,
                    mbm_cur->comm);

  *merge_block_stride  = dparent_child_strid;
  *merge_block_data    = dparent_child_mbm;

  for (int i_block = 0; i_block < n_block; i_block++) {
    free(shifted_block_data[i_block]);
  }
  free(shifted_block_data);
}


/*
 *
 */
int
PDM_multi_block_merge_get_n_block
(
 PDM_multi_block_merge_t* mbm
)
{
  return mbm->dn_new_block;
}


/*
 *
 */
PDM_g_num_t*
PDM_multi_block_merge_get_distrib
(
 PDM_multi_block_merge_t* mbm
)
{
  return mbm->distrib_merge;
}

PDM_g_num_t*
PDM_multi_block_merge_get_multi_distrib
(
 PDM_multi_block_merge_t* mbm
)
{
  return mbm->multi_block_distrib;
}


/*
 *
 */
void
PDM_multi_block_merge_free
(
 PDM_multi_block_merge_t* mbm
)
{

  free(mbm->multi_block_distrib);
  free(mbm->blocks_ini_dn);
  free(mbm->dold_to_new_idx);
  free(mbm->dold_to_new);
  free(mbm->distrib_merge);
  free(mbm->old_distrib);

  PDM_multi_block_to_part_free(mbm->mbtp);
  free(mbm);
}


// PDM_multi_block_merge_renum_from

#ifdef __cplusplus
}
#endif /* __cplusplus */
