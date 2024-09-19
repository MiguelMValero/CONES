/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multi_block_to_part.h"
#include "pdm_multi_block_to_part_priv.h"
#include "pdm_binary_search.h"
#include "pdm_array.h"
#include "pdm_priv.h"
#include "pdm_logging.h"

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

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   multi_distrib_idx Multiple block distribution (size : \ref size of \ref nblock + 1)
 * \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt          Element global number (size : \ref n_part)
 * \param [in]   n_elt             Local number of elements (size : \ref n_part)
 * \param [in]   n_part            Number of partition
 * \param [in]   comm              MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */

PDM_multi_block_to_part_t *
PDM_multi_block_to_part_create
(
 const PDM_g_num_t   *multi_distrib_idx,
 const int            n_block,
 const PDM_g_num_t  **block_distrib_idx,
 const PDM_g_num_t  **gnum_elt,
 const int           *n_elt,
 const int            n_part,
 const PDM_MPI_Comm   comm
)
{
  PDM_multi_block_to_part_t *mbtp =
    (PDM_multi_block_to_part_t *) malloc (sizeof(PDM_multi_block_to_part_t));

  mbtp->comm        = comm;
  mbtp->pttopt_comm = 0;

  PDM_MPI_Comm_size (comm, &mbtp->n_rank);
  PDM_MPI_Comm_rank (comm, &mbtp->i_rank);

  if(0 == 1) {
    printf(" PDM_multi_block_to_part_create :: n_part = %i \n", n_part);
    printf(" PDM_multi_block_to_part_create :: n_block = %i \n", n_block);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      printf(" PDM_multi_block_to_part_create :: n_elt[%i] = %i \n", i_part, n_elt[i_part]);
    }
  }

  mbtp->n_block = n_block;

  mbtp->multi_distrib_idx = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbtp->n_block + 1));
  for(int i_block = 0; i_block < mbtp->n_block+1; ++i_block){
    mbtp->multi_distrib_idx[i_block] = multi_distrib_idx[i_block];
  }

  mbtp->block_distrib_idx = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t*) * (mbtp->n_block));

  PDM_g_num_t shift = 0;
  for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
    mbtp->block_distrib_idx[i_block] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mbtp->n_rank + 1));
    for(int i = 0; i < mbtp->n_rank + 1; ++i) {
      mbtp->block_distrib_idx[i_block][i] = shift + block_distrib_idx[i_block][i];
    }
    shift += block_distrib_idx[i_block][mbtp->n_rank];
    // PDM_log_trace_array_long(mbtp->block_distrib_idx[i_block], mbtp->n_rank + 1, "mbtp->block_distrib_idx:: ");
  }

  mbtp->n_part  = n_part;

  mbtp->n_data_block = mbtp->n_rank * mbtp->n_block; /* All ranks have the same number of blocks */

  mbtp->requested_block_idx = PDM_array_zeros_int(mbtp->n_data_block + 1);
  mbtp->requested_block_n   = PDM_array_zeros_int(mbtp->n_data_block);

  mbtp->n_elt = malloc (sizeof(int  ) * n_part);
  mbtp->ind   = malloc (sizeof(int *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    mbtp->n_elt[i_part] = n_elt[i_part];
    mbtp->ind[i_part]   = malloc (sizeof(int) * n_elt[i_part]);

    const PDM_g_num_t *_gnum_elt = gnum_elt[i_part];

    for (int j = 0; j < n_elt[i_part]; j++) {

      int idx_block = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                 multi_distrib_idx,
                                                 mbtp->n_block + 1);

      int idx_rank     = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                    mbtp->block_distrib_idx[idx_block],
                                                    mbtp->n_rank + 1);

      //printf("[i_part:%i | ielm : %i -> i_block : %i | ind : % i\n", i_part, j, idx_block, idx_rank);
      int idx_data_block = idx_block + idx_rank*mbtp->n_block;
      mbtp->requested_block_n[idx_data_block]++;

    }
  }

  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->requested_block_idx[i+1] = mbtp->requested_block_idx[i] +
                                     mbtp->requested_block_n  [i];
  }

  int s_requested_data = mbtp->requested_block_idx[mbtp->n_data_block];

  for (int i = 0; i < mbtp->n_data_block; i++) {
    mbtp->requested_block_n[i] = 0;
  }

  int *requested_data = malloc( sizeof(int) * s_requested_data);
  for (int i_part = 0; i_part < n_part; i_part++) {

    const PDM_g_num_t *_gnum_elt = gnum_elt[i_part];

    for (int j = 0; j < n_elt[i_part]; j++) {

      int idx_block = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                 multi_distrib_idx,
                                                 mbtp->n_block + 1);

      int idx_rank     = PDM_binary_search_gap_long(_gnum_elt[j] - 1,
                                                    mbtp->block_distrib_idx[idx_block],
                                                    mbtp->n_rank + 1);

      int idx_data_block = idx_block + idx_rank*mbtp->n_block;

      int idx = mbtp->requested_block_idx[idx_data_block] + mbtp->requested_block_n[idx_data_block]++;

      mbtp->ind[i_part][j] = idx;

      PDM_g_num_t _mshift         = mbtp->block_distrib_idx[idx_block][idx_rank]; // - mbtp->multi_distrib_idx[idx_block];
      PDM_g_num_t _requested_data = _gnum_elt[j] - 1 - _mshift;

      requested_data[idx] = (int) _requested_data;
    }
  }

  mbtp->distributed_block_n = malloc (sizeof(int) * mbtp->n_data_block);

  PDM_MPI_Alltoall (mbtp->requested_block_n,   mbtp->n_block, PDM_MPI_INT,
                    mbtp->distributed_block_n, mbtp->n_block, PDM_MPI_INT,
                    comm);


  mbtp->distributed_block_idx = PDM_array_new_idx_from_sizes_int(mbtp->distributed_block_n, mbtp->n_data_block);

  if( 0 == 1) {
    PDM_log_trace_array_int(mbtp->requested_block_n    , mbtp->n_data_block  , "mbtp->requested_block_n    :: ");
    PDM_log_trace_array_int(mbtp->requested_block_idx  , mbtp->n_data_block+1, "mbtp->requested_block_idx  :: ");
    PDM_log_trace_array_int(mbtp->distributed_block_n  , mbtp->n_data_block  , "mbtp->distributed_block_n  :: ");
    PDM_log_trace_array_int(mbtp->distributed_block_idx, mbtp->n_data_block+1, "mbtp->distributed_block_idx:: ");
  }


  // Les data se deduisent des blocks
  mbtp->requested_data_idx   = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->requested_data_n     = malloc (sizeof(int) * (mbtp->n_rank    ));
  mbtp->distributed_data_idx = malloc (sizeof(int) * (mbtp->n_rank + 1));
  mbtp->distributed_data_n   = malloc (sizeof(int) * (mbtp->n_rank    ));

  for(int i = 0; i < mbtp->n_rank; ++i) {
    int ind = i*mbtp->n_block;
    mbtp->requested_data_idx  [i] = mbtp->requested_block_idx  [ind];
    mbtp->distributed_data_idx[i] = mbtp->distributed_block_idx[ind];
    mbtp->requested_data_n    [i] = 0;
    mbtp->distributed_data_n  [i] = 0;
    for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
      mbtp->requested_data_n    [i] += mbtp->requested_block_n    [i_block+ind];
      mbtp->distributed_data_n  [i] += mbtp->distributed_block_n  [i_block+ind];
    }
  }
  mbtp->requested_data_idx  [mbtp->n_rank] = mbtp->requested_block_idx  [mbtp->n_data_block];
  mbtp->distributed_data_idx[mbtp->n_rank] = mbtp->distributed_block_idx[mbtp->n_data_block];

  if( 0 == 1) {
    PDM_log_trace_array_int(mbtp->requested_data_n    , mbtp->n_rank  , "mbtp->requested_data_n    :: ");
    PDM_log_trace_array_int(mbtp->requested_data_idx  , mbtp->n_rank+1, "mbtp->requested_data_idx  :: ");
    PDM_log_trace_array_int(mbtp->distributed_data_n  , mbtp->n_rank  , "mbtp->distributed_data_n  :: ");
    PDM_log_trace_array_int(mbtp->distributed_data_idx, mbtp->n_rank+1, "mbtp->distributed_data_idx:: ");
  }

  mbtp->distributed_data = malloc(sizeof(int) * mbtp->distributed_data_idx[mbtp->n_rank]);

  PDM_MPI_Alltoallv (requested_data,
                     mbtp->requested_data_n,
                     mbtp->requested_data_idx,
                     PDM_MPI_INT,
                     mbtp->distributed_data,
                     mbtp->distributed_data_n,
                     mbtp->distributed_data_idx,
                     PDM_MPI_INT,
                     comm);

  if( 0 == 1) {
    PDM_log_trace_array_int(mbtp->distributed_data    , mbtp->distributed_data_idx[mbtp->n_rank], "mbtp->distributed_data:: ");
  }

  free (requested_data);

  return (PDM_multi_block_to_part_t *) mbtp;
}

/**
 *
 * \brief Initialize an exchange
 * (part_stride and part_data are allocated in function)
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride or NULL
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_multi_block_to_part_exch2
(
 PDM_multi_block_to_part_t   *mbtp,
 size_t                       s_data,
 PDM_stride_t                 t_stride,
 int                        **block_stride,
 void                       **block_data,
 int                       ***part_stride,
 void                      ***part_data
)
{

  unsigned char **_part_data;

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * (mbtp->n_rank + 1));
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * (mbtp->n_rank + 1));

  int *n_send_buffer = (int *) malloc (sizeof(int) * mbtp->n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * mbtp->n_rank);

  for (int i = 0; i < mbtp->n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }
  i_send_buffer[mbtp->n_rank] = 0;
  i_recv_buffer[mbtp->n_rank] = 0;

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  int n_rank1 = mbtp->n_rank - 1;

  int  *recv_stride  = NULL;
  int **_part_stride = NULL;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_send_stride = mbtp->distributed_data_idx[mbtp->n_rank];
    int s_recv_stride = mbtp->requested_data_idx  [mbtp->n_rank];

    int *send_stride = (int *) malloc (sizeof(int) * s_send_stride);
    recv_stride      = (int *) malloc (sizeof(int) * s_recv_stride);

    int idxs = 0;
    for (int i = 0; i < mbtp->n_rank; i++) {
      for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
        int* _block_stride = block_stride[i_block];
        int idx = i_block + i*mbtp->n_block;
        for(int ielt = mbtp->distributed_block_idx[idx]; ielt < mbtp->distributed_block_idx[idx+1]; ++ielt) {
          int ind = mbtp->distributed_data[ielt];
          send_stride[idxs++] = _block_stride[ind];
        }
      }
    }

    PDM_MPI_Alltoallv (send_stride,
                       mbtp->distributed_data_n,
                       mbtp->distributed_data_idx,
                       PDM_MPI_INT,
                       recv_stride,
                       mbtp->requested_data_n,
                       mbtp->requested_data_idx,
                       PDM_MPI_INT,
                       mbtp->comm);

    *part_stride = (int **) malloc(sizeof(int *) * mbtp->n_part);
    _part_stride = *part_stride;

    for (int i_part = 0; i_part < mbtp->n_part; i_part++) {

      _part_stride[i_part] = malloc (sizeof(int) * mbtp->n_elt[i_part]);

      for (int j = 0; j < mbtp->n_elt[i_part]; j++) {

        int ielt = mbtp->ind[i_part][j];
        _part_stride[i_part][j] = recv_stride[ielt];

      }
    }

    /*
     * Build buffers
     */
    for (int i = 0; i < mbtp->n_rank; i++) {

      int ibeg = mbtp->distributed_data_idx[i];
      int iend = mbtp->distributed_data_idx[i] + mbtp->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = mbtp->requested_data_idx[i];
      iend = mbtp->requested_data_idx[i] + mbtp->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }

    }

    free(send_stride);

    s_send_buffer = i_send_buffer[n_rank1] + n_send_buffer[n_rank1];
    s_recv_buffer = i_recv_buffer[n_rank1] + n_recv_buffer[n_rank1];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    int **block_stride_idx = (int **) malloc(sizeof(int*) * (mbtp->n_block));
    for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
      int n_elt_block = mbtp->block_distrib_idx[i_block][mbtp->i_rank+1] - mbtp->block_distrib_idx[i_block][mbtp->i_rank];
      // printf(" n_elt_block :: %i \n", n_elt_block);
      block_stride_idx[i_block] = PDM_array_new_idx_from_sizes_int(block_stride[i_block], n_elt_block);
    }

    int idx1 = 0;
    for (int i = 0; i < mbtp->n_rank; i++) {

      for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {

        unsigned char* _block_data       = (unsigned char *) block_data[i_block];
        int*           _block_stride_idx = block_stride_idx[i_block];
        int*           _block_stride     = block_stride    [i_block];

        int idx = i_block + i*mbtp->n_block;

        for(int ielt = mbtp->distributed_block_idx[idx]; ielt < mbtp->distributed_block_idx[idx+1]; ++ielt) {
          int ind_blk      = mbtp->distributed_data[ielt];
          int shift        = _block_stride_idx[ind_blk] * (int) s_data;
          int s_block_unit = _block_stride    [ind_blk] * (int) s_data;
          // printf("[%i] with ind = %i \n", i_block, ind);
          unsigned char *_block_data_deb = _block_data + shift;
          for(int i_data = 0; i_data < s_block_unit; ++i_data) {
            send_buffer[idx1++] = _block_data_deb[i_data];
          }
        }
      }
    }
    for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
      free(block_stride_idx[i_block]);
    }
    free(block_stride_idx);

  } else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < mbtp->n_rank; i++) {
      i_send_buffer[i+1] = i_send_buffer[i];
      i_recv_buffer[i+1] = i_recv_buffer[i];
      for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
        int idx = i_block + i*mbtp->n_block;
        int cst_stride   = block_stride[i_block][0];
        int s_block_unit = cst_stride * (int) s_data;
        i_send_buffer[i+1] += mbtp->distributed_block_n[idx] * s_block_unit;
        i_recv_buffer[i+1] += mbtp->requested_block_n  [idx] * s_block_unit;
        n_send_buffer[i] += mbtp->distributed_block_n[idx] * s_block_unit;
        n_recv_buffer[i] += mbtp->requested_block_n  [idx] * s_block_unit;
      }
    }

    if( 0 == 1) {
      PDM_log_trace_array_size_t(i_send_buffer, mbtp->n_rank+1, "i_send_buffer:: ");
      PDM_log_trace_array_size_t(i_recv_buffer, mbtp->n_rank+1, "i_recv_buffer:: ");
      PDM_log_trace_array_int   (n_send_buffer, mbtp->n_rank  , "n_send_buffer:: ");
      PDM_log_trace_array_int   (n_recv_buffer, mbtp->n_rank  , "n_recv_buffer:: ");
    }

    s_send_buffer = i_send_buffer[mbtp->n_rank];
    s_recv_buffer = i_recv_buffer[mbtp->n_rank];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    int idx1 = 0;
    for (int i = 0; i < mbtp->n_rank; i++) {
      for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {

        unsigned char* _block_data = (unsigned char *) block_data[i_block];
        int idx = i_block + i*mbtp->n_block;

        int cst_stride   = block_stride[i_block][0];
        int s_block_unit = cst_stride * (int) s_data;

        for(int ielt = mbtp->distributed_block_idx[idx]; ielt < mbtp->distributed_block_idx[idx+1]; ++ielt) {
          int ind = mbtp->distributed_data[ielt]*s_block_unit;
          for(int i_data = 0; i_data < s_block_unit; ++i_data) {
            send_buffer[idx1++] = _block_data[ind+i_data];
          }
        }
      }
    }
  }

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
                      mbtp->comm);
  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  *part_data = malloc(sizeof(unsigned char *) * mbtp->n_part);
  _part_data = (*(unsigned char ***) part_data);

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    int s_recv_elt = mbtp->requested_data_idx[n_rank1] + mbtp->requested_data_n[n_rank1];

    int **part_idx = malloc (sizeof(int *) * (mbtp->n_part ));
    int  *recv_idx = PDM_array_new_idx_from_sizes_int(recv_stride, s_recv_elt);

    for (int i_part = 0; i_part < mbtp->n_part; i_part++) {
      part_idx[i_part] = PDM_array_new_idx_from_sizes_int(_part_stride[i_part], mbtp->n_elt[i_part]);
    }

    for (int i_part = 0; i_part < mbtp->n_part; i_part++) {

      int s_part =  part_idx[i_part][mbtp->n_elt[i_part]] * (int) s_data;

      _part_data[i_part] = malloc(sizeof(unsigned char) * s_part);

      for (int j = 0; j < mbtp->n_elt[i_part]; j++) {

        int idx1  =  part_idx   [i_part][j] * (int) s_data;
        int n_elt = _part_stride[i_part][j] * (int) s_data;

        int idx2 = recv_idx[mbtp->ind[i_part][j]] * (int) s_data;

        for (int k = 0; k < n_elt; k++) {
          _part_data[i_part][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }

    for (int i = 0; i < mbtp->n_part; i++) {
      free (part_idx[i]);
    }
    free(recv_idx);
    free(part_idx);
    free(recv_stride);

  } else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    int cst_stride = 0;
    for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
      assert(block_stride[i_block][0] == block_stride[0][0]);
      cst_stride = block_stride[i_block][0];
    }
    const int s_block_unit = cst_stride * (int) s_data;

    for (int i = 0; i < mbtp->n_part; i++) {

      _part_data[i] = malloc(sizeof(unsigned char) * s_block_unit * mbtp->n_elt[i]);

      for (int j = 0; j < mbtp->n_elt[i]; j++) {

        int idx1  = j * s_block_unit;
        int idx2 = mbtp->ind[i][j] * s_block_unit;

        for (int k = 0; k < s_block_unit; k++) {
           _part_data[i][idx1+k] = recv_buffer[idx2+k];
        }
      }
    }
  }

  free(recv_buffer);
}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btp  Block to part structure
 *
 * \return       NULL
 */

PDM_multi_block_to_part_t *
PDM_multi_block_to_part_free
(
 PDM_multi_block_to_part_t *mbtp
)
{

  for(int i_block = 0; i_block < mbtp->n_block; ++i_block) {
    free (mbtp->block_distrib_idx[i_block]);
  }
  free(mbtp->block_distrib_idx);
  free(mbtp->multi_distrib_idx);

  for (int i = 0; i < mbtp->n_part; i++) {
    free (mbtp->ind[i]);
  }
  free (mbtp->ind);

  free (mbtp->n_elt);
  free (mbtp->distributed_data);
  free (mbtp->distributed_data_idx);
  free (mbtp->distributed_data_n);
  free (mbtp->distributed_block_idx);
  free (mbtp->distributed_block_n);
  free (mbtp->requested_data_idx);
  free (mbtp->requested_data_n);
  free (mbtp->requested_block_idx);
  free (mbtp->requested_block_n);


  free (mbtp);

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
