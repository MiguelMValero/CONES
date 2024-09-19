/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_to_part.h"
#include "pdm_part_to_part_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_gnum_location.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
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

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/**
 *
 * \brief Free the properties of a asynchronous alltoall
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 *
 */

static void
_free_async_alltoall
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  ptp->async_alltoall_subrequest[3 * request]     = -1;
  ptp->async_alltoall_subrequest[3 * request + 1] = -1;
  ptp->async_alltoall_subrequest[3 * request + 2] = -1;
  ptp->async_alltoall_free[ptp->async_alltoall_n_free++] = request;
}


/**
 *
 * \brief Free the properties of a asynchronous data sending
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 *
 */

static void
_free_async_send
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  if (ptp->async_send_request[request] != NULL) {
    free (ptp->async_send_request[request]);
  }
  ptp->async_send_request[request] = NULL;

  ptp->async_send_s_data[request]     = -1;
  ptp->async_send_cst_stride[request] = -1;
  ptp->async_send_tag[request]        = -1;
  if (ptp->async_send_buffer[request] != NULL) {
    free (ptp->async_send_buffer[request]);
  }
  ptp->async_send_buffer[request] = NULL;

  if ((ptp->async_n_send_buffer[request] != NULL) &&
      (ptp->async_n_send_buffer[request] != ptp->default_n_send_buffer)) {
    free (ptp->async_n_send_buffer[request]);
  }
  ptp->async_n_send_buffer[request] = NULL;

  if ((ptp->async_i_send_buffer[request] != NULL) &&
      (ptp->async_i_send_buffer[request] != ptp->default_i_send_buffer)) {
    free (ptp->async_i_send_buffer[request]);
  }
  ptp->async_i_send_buffer[request] = NULL;

  ptp->async_send_free[ptp->async_send_n_free++] = request;
}


/**
 *
 * \brief Free the properties of a asynchronous data reception
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 *
 */

static void
_free_async_recv
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  if (ptp->async_recv_request[request] != NULL) {
    free (ptp->async_recv_request[request]);
  }
  ptp->async_recv_request[request] = NULL;

  ptp->async_recv_s_data[request]     = -1;
  ptp->async_recv_cst_stride[request] = -1;
  ptp->async_recv_tag[request]        = -1;

  if (ptp->async_recv_buffer[request]  != NULL) {
    free (ptp->async_recv_buffer[request]);
  }
  ptp->async_recv_buffer[request]   = NULL;

  if ((ptp->async_n_recv_buffer[request] != NULL) &&
      (ptp->async_n_recv_buffer[request] != ptp->default_n_recv_buffer)) {
    free (ptp->async_n_recv_buffer[request]);
  }
  ptp->async_n_recv_buffer[request] = NULL;

  if ((ptp->async_i_recv_buffer[request] != NULL) &&
      (ptp->async_i_recv_buffer[request] != ptp->default_i_recv_buffer)) {
    free (ptp->async_i_recv_buffer[request]);
  }
  ptp->async_i_recv_buffer[request]   = NULL;

  if (ptp->async_recv_part2_data[request] != NULL) {
    free (ptp->async_recv_part2_data[request]);
    ptp->async_recv_part2_data[request] = NULL;
  }

  ptp->async_recv_free[ptp->async_recv_n_free++] = request;
}


/**
 *
 * \brief Free the asynchronous properties of an exchange
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   request       Request
 *
 */

static void
_free_async_exch
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  for (int i = 0; i < ptp->async_exch_subrequest_s[request]; i++) {
    ptp->async_exch_subrequest[request][2*i]   = -1;
    ptp->async_exch_subrequest[request][2*i+1] = -1;
  }

  if (ptp->async_exch_recv_n[request] != NULL) {
    free (ptp->async_exch_recv_n[request]);
    ptp->async_exch_recv_n[request] = NULL;
  }

  if (ptp->async_exch_recv_idx[request] != NULL) {
    free (ptp->async_exch_recv_idx[request]);
    ptp->async_exch_recv_idx[request] = NULL;
  }
  ptp->async_exch_part2_stride[request]   = NULL;
  ptp->async_exch_t_stride[request]       = -1;
  ptp->async_exch_k_comm[request]         = -1;


  ptp->async_exch_free[ptp->async_exch_n_free++] = request;
}


/**
 *
 * \brief Check arrays of the asynchronous send
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_alltoall_alloc
(
 PDM_part_to_part_t *ptp
)
{
  if (ptp->async_alltoall_l_array == 0) {
    ptp->async_alltoall_l_array    = 10;
    ptp->async_alltoall_free       = malloc (sizeof(int) * ptp->async_alltoall_l_array);
    ptp->async_alltoall_subrequest = malloc (sizeof(int) * 3 * ptp->async_alltoall_l_array);

    for (int i = 0; i < ptp->async_alltoall_l_array; i++) {
      ptp->async_alltoall_free[ptp->async_alltoall_n_free++] = ptp->async_alltoall_l_array -1 - i;
      ptp->async_alltoall_subrequest[3*i]   = -1;
      ptp->async_alltoall_subrequest[3*i+1] = -1;
      ptp->async_alltoall_subrequest[3*i+2] = -1;
    }
  }

  if (ptp->async_alltoall_n_free == 0) {
    const int pre_val = ptp->async_alltoall_l_array;
    ptp->async_alltoall_l_array   *= 2;
    ptp->async_alltoall_free       = realloc (ptp->async_alltoall_free      , sizeof(int) * ptp->async_alltoall_l_array);
    ptp->async_alltoall_subrequest = realloc (ptp->async_alltoall_subrequest, sizeof(int) * 3 * ptp->async_alltoall_l_array);

    for (int i = pre_val; i < ptp->async_alltoall_l_array; i++) {
      ptp->async_alltoall_free[ptp->async_alltoall_n_free++] = i;
      ptp->async_alltoall_subrequest[3*i]   = -1;
      ptp->async_alltoall_subrequest[3*i+1] = -1;
      ptp->async_alltoall_subrequest[3*i+2] = -1;
    }
  }
}


/**
 *
 * \brief Check arrays of the asynchronous send
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_send_alloc
(
 PDM_part_to_part_t *ptp
)
{
  if (ptp->async_send_l_array == 0) {
    ptp->async_send_l_array    = 10;
    ptp->async_send_s_data     = malloc (sizeof(size_t) * ptp->async_send_l_array);
    ptp->async_send_cst_stride = malloc (sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_tag        = malloc (sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_request    = malloc (sizeof(PDM_MPI_Request *) * ptp->async_send_l_array);
    ptp->async_send_buffer     = malloc (sizeof(unsigned char *) * ptp->async_send_l_array);
    ptp->async_n_send_buffer   = malloc (sizeof(int *) * ptp->async_send_l_array);
    ptp->async_i_send_buffer   = malloc (sizeof(int *) * ptp->async_send_l_array);
    ptp->async_send_free       = malloc (sizeof(int) * ptp->async_send_l_array);

    for (int i = 0; i < ptp->async_send_l_array; i++) {
      ptp->async_send_free[ptp->async_send_n_free++] = ptp->async_send_l_array -1 - i;
      ptp->async_send_s_data[i]     = -1;
      ptp->async_send_cst_stride[i] = -1;
      ptp->async_send_tag[i]        = -1;
      ptp->async_send_request[i]    = NULL;
      ptp->async_send_buffer[i]     = NULL;
      ptp->async_n_send_buffer[i]   = NULL;
      ptp->async_i_send_buffer[i]   = NULL;
    }
  }

  if (ptp->async_send_n_free == 0) {
    const int pre_val = ptp->async_send_l_array;
    ptp->async_send_l_array *= 2;
    ptp->async_send_free       = realloc (ptp->async_send_free       , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_s_data     = realloc (ptp->async_send_s_data     , sizeof(size_t) * ptp->async_send_l_array);
    ptp->async_send_cst_stride = realloc (ptp->async_send_cst_stride , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_tag        = realloc (ptp->async_send_tag        , sizeof(int) * ptp->async_send_l_array);
    ptp->async_send_request    = realloc (ptp->async_send_request    , sizeof(PDM_MPI_Request *) * ptp->async_send_l_array);
    ptp->async_send_buffer     = realloc (ptp->async_send_buffer     , sizeof(unsigned char *) * ptp->async_send_l_array);
    ptp->async_n_send_buffer   = realloc (ptp->async_n_send_buffer   , sizeof(int *) * ptp->async_send_l_array);
    ptp->async_i_send_buffer   = realloc (ptp->async_i_send_buffer   , sizeof(int *) * ptp->async_send_l_array);

    for (int i = pre_val; i < ptp->async_send_l_array; i++) {
      ptp->async_send_free[ptp->async_send_n_free++] = i;
      ptp->async_send_s_data[i]     = -1;
      ptp->async_send_cst_stride[i] = -1;
      ptp->async_send_tag[i]        = -1;
      ptp->async_send_request[i]    = NULL;
      ptp->async_send_buffer[i]     = NULL;
      ptp->async_n_send_buffer[i]   = NULL;
      ptp->async_i_send_buffer[i]   = NULL;
    }
  }
}


/**
 *
 * \brief Check arrays of the asynchronous send
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_recv_alloc
(
 PDM_part_to_part_t *ptp
)
{
  if (ptp->async_recv_l_array == 0) {
    ptp->async_recv_l_array    = 10;
    ptp->async_recv_s_data     = malloc (sizeof(size_t) * ptp->async_recv_l_array);
    ptp->async_recv_cst_stride = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_tag        = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_request    = malloc (sizeof(PDM_MPI_Request *) * ptp->async_recv_l_array);
    ptp->async_recv_buffer     = malloc (sizeof(unsigned char *) * ptp->async_recv_l_array);
    ptp->async_n_recv_buffer   = malloc (sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_i_recv_buffer   = malloc (sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_recv_free       = malloc (sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_part2_data = malloc (sizeof(void *) * ptp->async_recv_l_array);

    for (int i = 0; i < ptp->async_recv_l_array; i++) {
      ptp->async_recv_free[ptp->async_recv_n_free++] = ptp->async_recv_l_array -1 - i;
      ptp->async_recv_s_data[i]     = -1;
      ptp->async_recv_cst_stride[i] = -1;
      ptp->async_recv_tag[i]        = -1;
      ptp->async_recv_request[i]    = NULL;
      ptp->async_recv_buffer[i]     = NULL;
      ptp->async_n_recv_buffer[i]   = NULL;
      ptp->async_i_recv_buffer[i]   = NULL;
      ptp->async_recv_part2_data[i] = NULL;
    }
  }

  if (ptp->async_recv_n_free == 0) {
    const int pre_val = ptp->async_recv_l_array;
    ptp->async_recv_l_array *= 2;
    ptp->async_recv_free       = realloc (ptp->async_recv_free       , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_s_data     = realloc (ptp->async_recv_s_data     , sizeof(size_t) * ptp->async_recv_l_array);
    ptp->async_recv_cst_stride = realloc (ptp->async_recv_cst_stride , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_tag        = realloc (ptp->async_recv_tag        , sizeof(int) * ptp->async_recv_l_array);
    ptp->async_recv_request    = realloc (ptp->async_recv_request    , sizeof(PDM_MPI_Request *) * ptp->async_recv_l_array);
    ptp->async_recv_buffer     = realloc (ptp->async_recv_buffer     , sizeof(unsigned char *) * ptp->async_recv_l_array);
    ptp->async_n_recv_buffer   = realloc (ptp->async_n_recv_buffer   , sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_i_recv_buffer   = realloc (ptp->async_i_recv_buffer   , sizeof(int *) * ptp->async_recv_l_array);
    ptp->async_recv_part2_data = realloc (ptp->async_recv_part2_data , sizeof(void *) * ptp->async_recv_l_array);

    for (int i = pre_val; i < ptp->async_recv_l_array; i++) {
      ptp->async_recv_free[ptp->async_recv_n_free++] = i;
      ptp->async_recv_s_data[i]     = -1;
      ptp->async_recv_cst_stride[i] = -1;
      ptp->async_recv_tag[i]        = -1;
      ptp->async_recv_request[i]    = NULL;
      ptp->async_recv_buffer[i]     = NULL;
      ptp->async_n_recv_buffer[i]   = NULL;
      ptp->async_i_recv_buffer[i]   = NULL;
      ptp->async_recv_part2_data[i] = NULL;
    }
  }
}

/**
 *
 * \brief Check arrays of the asynchronous data exchange
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static void
_check_async_exch_alloc
(
 PDM_part_to_part_t *ptp
)
{
  if (ptp->async_recv_l_array == 0) {
    ptp->async_exch_l_array       = 10;
    ptp->async_exch_free          = malloc (sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_subrequest_s  = malloc (sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_subrequest    = malloc (sizeof(int *) * ptp->async_exch_l_array);
    ptp->async_exch_t_stride      = malloc (sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_k_comm        = malloc (sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_recv_n        = malloc (sizeof(int *) * ptp->async_exch_l_array);
    ptp->async_exch_recv_idx      = malloc (sizeof(int *) * ptp->async_exch_l_array);
    ptp->async_exch_part2_stride  = malloc (sizeof(int **) * ptp->async_exch_l_array);

    for (int i = 0; i < ptp->async_recv_l_array; i++) {
      ptp->async_exch_free[ptp->async_exch_n_free++] = ptp->async_exch_l_array -1 - i;
      ptp->async_exch_recv_n[i]   = NULL;
      ptp->async_exch_recv_idx[i] = NULL;
      ptp->async_exch_part2_stride[i] = NULL;
      ptp->async_exch_subrequest_s[i] = 1;
      ptp->async_exch_subrequest[i] = malloc(sizeof(int) * 2 * ptp->async_exch_subrequest_s[i]);
      for (int j = 0; j < ptp->async_exch_subrequest_s[i]; j++) {
        ptp->async_exch_subrequest[i][j] = -1;
      }
      ptp->async_exch_t_stride[i] = -1;
      ptp->async_exch_k_comm[i]   = -1;
    }

  }

  if (ptp->async_recv_n_free == 0) {
    const int pre_val = ptp->async_recv_l_array;
    ptp->async_exch_l_array      *= 2;
    ptp->async_exch_free          = realloc (ptp->async_exch_free       , sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_subrequest_s  = realloc (ptp->async_exch_subrequest_s , sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_subrequest    = realloc (ptp->async_exch_subrequest , sizeof(int *) *  ptp->async_exch_l_array);
    ptp->async_exch_t_stride      = realloc (ptp->async_exch_t_stride,    sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_k_comm        = realloc (ptp->async_exch_k_comm,      sizeof(int) * ptp->async_exch_l_array);
    ptp->async_exch_recv_n        = realloc (ptp->async_exch_recv_n     , sizeof(int *) * ptp->async_exch_l_array);
    ptp->async_exch_recv_idx      = realloc (ptp->async_exch_recv_idx   , sizeof(int *) * ptp->async_exch_l_array);
    ptp->async_exch_part2_stride  = realloc (ptp->async_exch_part2_stride, sizeof(int **) * ptp->async_exch_l_array);

    for (int i = pre_val; i < ptp->async_exch_l_array; i++) {
      ptp->async_exch_free[ptp->async_exch_n_free++] = i;
      ptp->async_exch_recv_n[i]   = NULL;
      ptp->async_exch_recv_idx[i] = NULL;
      ptp->async_exch_part2_stride[i] = NULL;
      ptp->async_exch_t_stride[i] = -1;
      ptp->async_exch_k_comm[i]   = -1;
      ptp->async_exch_subrequest_s[i] = 1;
      ptp->async_exch_subrequest[i] = malloc(sizeof(int) * 2 * ptp->async_exch_subrequest_s[i]);
      for (int j = 0; j < ptp->async_exch_subrequest_s[i]; j++) {
        ptp->async_exch_subrequest[i][j] = -1;
      }
    }
  }
}


/**
 *
 * \brief Initialize an asynchronous data alltoall
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static int
_find_open_async_alltoall_exch
(
 PDM_part_to_part_t *ptp
)
{
  _check_async_alltoall_alloc (ptp);
  return ptp->async_alltoall_free[--ptp->async_alltoall_n_free];
}

/**
 *
 * \brief Initialize an asynchronous data sending
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static int
_find_open_async_send_exch
(
 PDM_part_to_part_t *ptp
)
{
  _check_async_send_alloc (ptp);

  return ptp->async_send_free[--ptp->async_send_n_free];
}

/**
 *
 * \brief Initialize an asynchronous data reception
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static int
_find_open_async_recv_exch
(
 PDM_part_to_part_t *ptp
)
{
  _check_async_recv_alloc (ptp);

  return ptp->async_recv_free[--ptp->async_recv_n_free];
}

/**
 *
 * \brief Initialize an asynchronous data exchange
 *
 * \param [in]   ptp           Block to part structure
 *
 */

static int
_find_open_async_exch
(
 PDM_part_to_part_t *ptp
)
{
  _check_async_exch_alloc (ptp);

  return ptp->async_exch_free[--ptp->async_exch_n_free];
}


/**
 *
 * \brief Wait a asynchronous stride reception
 *
 * \param [in]  ptp                  Part to part structure
 * \param [out] MPI_buffer_recv_n    Number of data received from each rank
 * \param [out] MPI_buffer_recv_idx  Index in receive MPI buffer used to receive data
 * \param [in]  request              Request
 *
 */

static void
_p2p_stride_var_irecv_stride_wait
(
 PDM_part_to_part_t  *ptp,
 int                **MPI_buffer_recv_n,
 int                **MPI_buffer_recv_idx,
 int                  request
)
{

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }

  size_t s_data  = ptp->async_recv_s_data[request];
  int cst_stride = ptp->async_recv_cst_stride[request];

  unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request];

  int n_blk_recv = ptp->async_i_recv_buffer[request][ptp->n_rank]/sizeof(int);

  int* blk_recv_stride = (int*) ptp->async_recv_buffer[request];
  // PDM_log_trace_array_int(blk_recv_stride,
  //                         n_blk_recv,
  //                         "blk_recv_stride : ");

  int* _MPI_buffer_recv_idx = malloc( (n_blk_recv + 1) * sizeof(int) );
  int* _MPI_buffer_recv_n   = (int * ) malloc(ptp->n_rank * sizeof(int));

  for(int i = 0; i < ptp->n_rank; ++i) {
    _MPI_buffer_recv_n[i] = 0;
  }
  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int dest = ptp->active_rank_recv[i];
    int beg =       ptp->async_i_recv_buffer[request][dest]/sizeof(int);
    int end = beg + ptp->async_n_recv_buffer[request][dest]/sizeof(int);
    for(int j = beg; j < end; ++j) {
      _MPI_buffer_recv_n[dest] += blk_recv_stride[j];
    }
  }

  _MPI_buffer_recv_idx[0] = 0;
  for(int i = 0; i < n_blk_recv; ++i) {
    _MPI_buffer_recv_idx[i+1] = _MPI_buffer_recv_idx[i] + blk_recv_stride[i];
  }

  *MPI_buffer_recv_n    = _MPI_buffer_recv_n;
  *MPI_buffer_recv_idx = _MPI_buffer_recv_idx;

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
        int idx = ptp->recv_buffer_to_ref_lnum2[i][k] * delta;
        int idx1 = k* delta;
        for (int k1 = 0; k1 < delta; k1++) {
          _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1];
        }
      }
    }
  }

  _free_async_recv (ptp, request);

}


/**
 *
 * \brief Wait a asynchronous stride reception in reverse exchange
 *
 * \param [in]  ptp                  Part to part structure
 * \param [out] MPI_buffer_recv_n    Number of data received from each rank
 * \param [out] MPI_buffer_recv_idx  Index in receive MPI buffer used to receive data
 * \param [in]  request              Request
 *
 */

static void
_p2p_stride_var_reverse_irecv_stride_wait
(
 PDM_part_to_part_t  *ptp,
 int                **MPI_buffer_recv_n,
 int                **MPI_buffer_recv_idx,
 int                  request
)
{

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }

  size_t s_data  = ptp->async_recv_s_data[request];
  int cst_stride = ptp->async_recv_cst_stride[request];

  int n_blk_recv = ptp->async_i_recv_buffer[request][ptp->n_rank]/sizeof(int);

  int* blk_recv_stride = (int*) ptp->async_recv_buffer[request];

  int* _MPI_buffer_recv_idx = malloc( (n_blk_recv + 1) * sizeof(int) );
  int* _MPI_buffer_recv_n   = (int * ) malloc(ptp->n_rank * sizeof(int));

  for(int i = 0; i < ptp->n_rank; ++i) {
    _MPI_buffer_recv_n[i] = 0;
  }
  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    int beg =       ptp->async_i_recv_buffer[request][dest]/sizeof(int);
    int end = beg + ptp->async_n_recv_buffer[request][dest]/sizeof(int);
    for(int j = beg; j < end; ++j) {
      _MPI_buffer_recv_n[dest] += blk_recv_stride[j];
    }
  }

  _MPI_buffer_recv_idx[0] = 0;
  for(int i = 0; i < n_blk_recv; ++i) {
    _MPI_buffer_recv_idx[i+1] = _MPI_buffer_recv_idx[i] + blk_recv_stride[i];
  }

  // PDM_log_trace_array_int(_MPI_buffer_recv_n, ptp->n_rank, "MPI_buffer_recv_n : ");
  // PDM_log_trace_array_int(_MPI_buffer_recv_idx, n_blk_recv+1, "MPI_buffer_recv_idx : ");


  *MPI_buffer_recv_n    = _MPI_buffer_recv_n;
  *MPI_buffer_recv_idx = _MPI_buffer_recv_idx;

  unsigned char ** _part1_data = (unsigned char **) ptp->async_recv_part2_data[request];

  int delta = (int) s_data * cst_stride;

  for (int i = 0; i < ptp->n_part1; i++) {
    for (int i1 = 0; i1 < ptp->n_elt1[i]; i1++) {
      for (int j = ptp->part1_to_part2_idx[i][i1]; j < ptp->part1_to_part2_idx[i][i1+1]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
                 k < ptp->gnum1_to_send_buffer_idx[i][j+1];
                 k++) {

          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
            int idx  = ptp->gnum1_to_send_buffer[i][k] * delta;
            int idx1 = j * delta;
            for (int k1 = 0; k1 < delta; k1++) {
              _part1_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1];
            }
          }
        }
      }
    }
  }

  _free_async_recv (ptp, request);

}


/**
 *
 * \brief Initialize a asynchronus data reception for an exchange with variable stride
 *
 * \param [in]  ptp               Part to part structure
 * \param [in]  s_data            Data size
 * \param [in]  MPI_buffer_recv_n Number of data received from each rank (stride is taking into acount)
 * \param [out] part2_data        Partition 2 data
 * \param [in]  tag               Tag of the exchange
 * \param [out] request           Request
 *
 */

static void
_p2p_stride_var_data_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int          *MPI_buffer_recv_n,
 void              **part2_data,
 int                 tag,
 int                *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;

  ptp->async_recv_s_data[_request]      = s_data;
  ptp->async_recv_cst_stride[_request]  = -1;
  ptp->async_recv_tag[_request]         = tag;

  ptp->async_recv_part2_data[_request]  = malloc(sizeof (void *) * ptp->n_part2);
  memcpy(ptp->async_recv_part2_data[_request], part2_data, sizeof (void *) * ptp->n_part2);

  ptp->async_recv_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i]   = MPI_buffer_recv_n[i] * (int) s_data;
    ptp->async_i_recv_buffer[_request][i+1] = ptp->async_i_recv_buffer[_request][i] + ptp->async_n_recv_buffer[_request][i];
  }
  ptp->async_recv_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[_request][ptp->n_rank]);

  // PDM_log_trace_array_int(ptp->async_i_recv_buffe0r[request], ptp->n_rank+1, "async_i_recv_buffer 1 : ");

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int source = ptp->active_rank_recv[i];
    unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source,
                    tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }

}

/**
 *
 * \brief Initialize a asynchronus data reception for a reverse exchange with variable stride
 *
 * \param [in]  ptp               Part to part structure
 * \param [in]  s_data            Data size
 * \param [in]  MPI_buffer_recv_n Number of data received from each rank (stride is taking into acount)
 * \param [out] part2_data        Partition 2 data
 * \param [in]  tag               Tag of the exchange
 * \param [out] request           Request
 *
 */

static void
_p2p_stride_var_data_reverse_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int          *MPI_buffer_recv_n,
 void              **part1_data,
 int                 tag,
 int                *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;

  ptp->async_recv_s_data[_request]      = s_data;
  ptp->async_recv_cst_stride[_request]  = -1;
  ptp->async_recv_tag[_request]         = tag;

  ptp->async_recv_part2_data[_request]  = malloc(sizeof (void *) * ptp->n_part1);
  memcpy(ptp->async_recv_part2_data[_request], part1_data, sizeof (void *) * ptp->n_part1);

  ptp->async_recv_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;

  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i]   = MPI_buffer_recv_n[i] * (int) s_data;
    ptp->async_i_recv_buffer[_request][i+1] = ptp->async_i_recv_buffer[_request][i] + ptp->async_n_recv_buffer[_request][i];
  }
  ptp->async_recv_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[_request][ptp->n_rank]);

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int source = ptp->active_rank_send[i];
    unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source,
                   tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }

}


/**
 *
 * \brief Initialize a data sending for an exchange with variable stride
 *
 * \param [in]  ptp                   Part to part structure
 * \param [in]  tag                   Tag for p2p exchange
 * \param [in]  s_data                Data size
 * \param [in]  MPI_buffer_send_n     Number of data received from each rank
 * \param [in]  MPI_buffer_send_idx   Index in receive MPI buffer used to receive data
 * \param [in]  part1_to_part2_stride Stride of partition 1 data
 * \param [in]  part1_to_part1_data   Partition 1 data
 * \param [out] request               Request
 *
 */


static void
_p2p_stride_var_data_issend
(
 PDM_part_to_part_t *ptp,
 const int           tag,
 const size_t        s_data,
 const int*          MPI_buffer_send_n,
 const int*          MPI_buffer_send_idx,
 int               **part1_to_part2_stride,
 void              **part1_to_part2_data,
 int                *request
)
{
  unsigned char ** _part1_data = (unsigned char **) part1_to_part2_data;

  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data    [_request]    = s_data;
  ptp->async_send_cst_stride[_request]    = 1;
  ptp->async_send_tag       [_request]    = tag;
  ptp->async_send_request   [_request]    = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_send_buffer  [_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer  [_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer  [_request][0] = 0;

  for (int i = 0; i < ptp->n_rank; i++) {
    // ptp->async_n_send_buffer[_request][i]   = send_n[i] * ptp->default_n_send_buffer[i  ] * (int) s_data;
    ptp->async_n_send_buffer[_request][i]   = MPI_buffer_send_n[i] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = ptp->async_i_send_buffer[_request][i] + ptp->async_n_send_buffer[_request][i];
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);

  /*
   * Compute idx
   */

  int **part1_to_part2_data_idx = malloc(ptp->n_part1 * sizeof(int * ));
  for (int i = 0; i < ptp->n_part1; i++) {
    part1_to_part2_data_idx[i] = malloc((ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]+1) * sizeof(int));
    part1_to_part2_data_idx[i][0] = 0;
    for(int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      part1_to_part2_data_idx[i][j+1] = part1_to_part2_data_idx[i][j] + part1_to_part2_stride[i][j];
    }
  }

  for (int i = 0; i < ptp->n_part1; i++) {
    for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
               k < ptp->gnum1_to_send_buffer_idx[i][j+1];
               k++) {

        if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
          int idx_elmt = ptp->gnum1_to_send_buffer[i][k];
          int idx      = MPI_buffer_send_idx[idx_elmt] * (int) s_data;
          int idx1     = part1_to_part2_data_idx[i][j] * (int) s_data;

          int delta    = part1_to_part2_stride[i][j] * s_data;

          // log_trace(" send at : (i=%i, j=%i / k=%i ) - idx_elmt = %i | idx = %i | idx1 = %i | delta = %i \n", i, j, k, idx_elmt, idx, idx1, delta);

          for (int k1 = 0; k1 < delta; k1++) {
            ptp->async_send_buffer[_request][idx+k1] = _part1_data[i][idx1+k1];
          }
        }
      }
    }
  }

  for (int i = 0; i < ptp->n_part1; i++) {
    free(part1_to_part2_data_idx[i]);
  }
  free(part1_to_part2_data_idx);

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    unsigned char *buf =  ptp->async_send_buffer[_request] + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
                    tag, ptp->comm, &(ptp->async_send_request[_request][i]));
  }
}




/**
 *
 * \brief Initialize a data sending for an exchange with variable stride
 *
 * \param [in]  ptp                   Part to part structure
 * \param [in]  tag                   Tag for p2p exchange
 * \param [in]  s_data                Data size
 * \param [in]  MPI_buffer_send_n     Number of data received from each rank
 * \param [in]  MPI_buffer_send_idx   Index in receive MPI buffer used to receive data
 * \param [in]  part1_to_part2_stride Stride of partition 1 data
 * \param [in]  part1_to_part1_data   Partition 1 data
 * \param [out] request               Request
 *
 */


static void
_p2p_stride_var_data_reverse_issend
(
 PDM_part_to_part_t *ptp,
 const int           tag,
 const size_t        s_data,
 const int*          MPI_buffer_send_n,
 const int*          MPI_buffer_send_idx,
 int               **part2_to_part1_stride,
 void              **part2_to_part1_data,
 int                *request
)
{
  unsigned char ** _part2_data = (unsigned char **) part2_to_part1_data;

  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data    [_request]    = s_data;
  ptp->async_send_cst_stride[_request]    = 1;
  ptp->async_send_tag       [_request]    = tag;
  ptp->async_send_request   [_request]    = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_send_buffer  [_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer  [_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer  [_request][0] = 0;

  for (int i = 0; i < ptp->n_rank; i++) {
    // ptp->async_n_send_buffer[_request][i]   = send_n[i] * ptp->default_n_send_buffer[i  ] * (int) s_data;
    ptp->async_n_send_buffer[_request][i]   = MPI_buffer_send_n[i] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = ptp->async_i_send_buffer[_request][i] + ptp->async_n_send_buffer[_request][i];
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);

  /*
   * Compute idx
   */

  int **part2_to_part1_data_idx = malloc(ptp->n_part2 * sizeof(int *));
  for (int i = 0; i < ptp->n_part2; i++) {

    part2_to_part1_data_idx[i]    = malloc((ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]+1) * sizeof(int));//
    part2_to_part1_data_idx[i][0] = 0;

    for (int j = 0; j < ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]; j++) {

      part2_to_part1_data_idx[i][j+1] = part2_to_part1_data_idx[i][j] + part2_to_part1_stride[i][j];

    }

  }

  int delta = s_data;

  for (int i = 0; i < ptp->n_part2; i++) {

    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {

        int idx = MPI_buffer_send_idx[ptp->recv_buffer_to_ref_lnum2[i][k]] * delta;
        int idx1 = part2_to_part1_data_idx[i][k] * delta;

        for (int k1 = 0; k1 < part2_to_part1_stride[i][k] * delta; k1++) {
          ptp->async_send_buffer[_request][idx+k1] = _part2_data[i][idx1+k1];
        }

      }

      for (int k = ptp->recv_buffer_to_duplicate_idx[i][j]; k < ptp->recv_buffer_to_duplicate_idx[i][j+1]; k++) {

        int idx      = MPI_buffer_send_idx[ptp->recv_buffer_to_duplicate[i][2*k  ]] * delta;
        int idx_data = part2_to_part1_data_idx[i][ptp->recv_buffer_to_duplicate[i][2*k+1]] * delta;

        for (int k1 = 0; k1 < part2_to_part1_stride[i][ptp->recv_buffer_to_duplicate[i][2*k+1]] * delta; k1++) {
          ptp->async_send_buffer[_request][idx+k1] = _part2_data[i][idx_data+k1];
        }

      }
    }
  }

  for (int i = 0; i < ptp->n_part2; i++) {
    free(part2_to_part1_data_idx[i]);
  }
  free(part2_to_part1_data_idx);

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int dest = ptp->active_rank_recv[i];
    unsigned char *buf =  ptp->async_send_buffer[_request] + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
                    tag, ptp->comm, &(ptp->async_send_request[_request][i]));
    // unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][dest];
    // int count = ptp->async_n_recv_buffer[_request][dest];
    // PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
    //                 tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }
}


/**
 *
 * \brief Initialize an asynchronus exchange with p2p communication kind
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   tag              Tag for p2p exchange
 * \param [in]   t_part1_data_def Kind of part1 data definition
 * \param [in]   s_data           Data size
 * \param [in]   part1_stride     Stride of partition 1 data
 * \param [in]   part1_data       Partition 1 data
 * \param [out]  part2_stride     Stride of partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  part2_data       Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  request          Request
 *
 */

static void
_alltotall_stride_var_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_part_to_part_data_def_t  t_part1_data_def,
 const size_t                       s_data,
 const int                        **part1_stride,
 const void                       **part1_data,
 int                             ***part2_stride,
 void                            ***part2_data,
 int                               *request
)
{

  *request         = _find_open_async_alltoall_exch(ptp);
  int request_send = _find_open_async_send_exch    (ptp);
  int request_recv = _find_open_async_recv_exch    (ptp);

  int _request = *request;
  ptp->async_alltoall_subrequest[3 * _request]     = request_send;
  ptp->async_alltoall_subrequest[3 * _request + 1] = request_recv;

  int            **_part1_to_part2_stride  = (int           **) part1_stride;
  unsigned char  **_part1_to_part2_data    = (unsigned char **) part1_data;
  int            **__part1_to_part2_stride = NULL;
  unsigned char  **__part1_to_part2_data   = NULL;

  /*
   *  Create __part1_to_part2_stride and __part1_to_part2_data if necessary
   */
  if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {
    __part1_to_part2_stride = (int           **) malloc (sizeof (int           *) * ptp->n_part1);
    __part1_to_part2_data   = (unsigned char **) malloc (sizeof (unsigned char *) * ptp->n_part1);

    _part1_to_part2_stride = __part1_to_part2_stride;
    _part1_to_part2_data   = __part1_to_part2_data;

    for (int i = 0; i < ptp->n_part1; i++) {
      _part1_to_part2_stride[i] = malloc (sizeof(int) * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);
      size_t k = 0;
      size_t s_part_data = 0;
      for (int j = 0; j < ptp->n_elt1[i]; j++) {
        for (int j1 = ptp->part1_to_part2_idx[i][j]; j1 < ptp->part1_to_part2_idx[i][j+1]; j1++) {
          _part1_to_part2_stride[i][k++] = part1_stride[i][j];
          s_part_data += part1_stride[i][j];
        }
      }

      _part1_to_part2_data[i] = malloc (s_data * s_part_data);
      unsigned char *map_part1_to_part2_data = (unsigned char*) _part1_to_part2_data[i];

      int beg_data = 0;
      k = 0;
      for (int j = 0; j < ptp->n_elt1[i]; j++) {
        unsigned char *tmp_part1_data = (unsigned char*) (part1_data[i]) + beg_data;
        for (int j1 = ptp->part1_to_part2_idx[i][j]; j1 < ptp->part1_to_part2_idx[i][j+1]; j1++) {
          for (int j2 = 0; j2 < (int) (part1_stride[i][j] * s_data); j2++) {
            map_part1_to_part2_data[k++] = tmp_part1_data[j2];
          }
        }
        beg_data += part1_stride[i][j] * s_data;
      }
    }
  }

  /*
   *  Stride exchange if necessary
   */
  int stride2_unknown = ((*part2_stride) == NULL);

  int n_blk_send = 0;
  int *blk_send_stride = NULL;
  int *blk_recv_stride = NULL;

  int **_part2_stride = NULL;
  if (stride2_unknown) {

    blk_send_stride = malloc (ptp->default_i_send_buffer[ptp->n_rank] * sizeof (int));
    blk_recv_stride = malloc (ptp->default_i_recv_buffer[ptp->n_rank] * sizeof (int));

    // Exchange stride
    for (int i = 0; i < ptp->n_part1; i++) {
      for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j]; k < ptp->gnum1_to_send_buffer_idx[i][j+1]; k++) {
          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
            int idx = ptp->gnum1_to_send_buffer[i][k];
            blk_send_stride[idx] = _part1_to_part2_stride[i][j];
          }
        }
      }
    }

    // PDM_log_trace_array_int(blk_send_strid, ptp->default_i_send_buffer[ptp->n_rank], "blk_send_strid :: ");

    PDM_MPI_Alltoallv(blk_send_stride,
                      ptp->default_n_send_buffer,
                      ptp->default_i_send_buffer,
                      PDM_MPI_INT,
                      blk_recv_stride,
                      ptp->default_n_recv_buffer,
                      ptp->default_i_recv_buffer,
                      PDM_MPI_INT,
                      ptp->comm);

    /*
     * Post-treatment stride
     */
    assert(_part2_stride == NULL);
    _part2_stride = malloc( ptp->n_part2 * sizeof(int*));
    for(int i = 0; i < ptp->n_part2; ++i) {
      _part2_stride[i] = malloc( ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]] * sizeof(int));
    }

    for (int i = 0; i < ptp->n_part2; i++) {
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          int idx = ptp->recv_buffer_to_ref_lnum2[i][k];
          _part2_stride[i][k] = blk_recv_stride[idx];
        }
      }
    }

    // PDM_log_trace_array_int(blk_recv_strid, ptp->default_i_recv_buffer[ptp->n_rank], "blk_recv_strid :: ");
  } else {

    n_blk_send = ptp->default_i_send_buffer[ptp->n_rank];
    blk_send_stride = malloc(sizeof(int) * n_blk_send);
    for (int i = 0; i < ptp->n_part1; i++) {
      for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j]; k < ptp->gnum1_to_send_buffer_idx[i][j+1]; k++) {
          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
            int idx = ptp->gnum1_to_send_buffer[i][k];
            blk_send_stride[idx] = _part1_to_part2_stride[i][k];
          }
        }
      }
    }

    _part2_stride = *part2_stride;
    int n_blk_recv = ptp->default_i_recv_buffer[ptp->n_rank];
    blk_recv_stride = malloc(sizeof(int) * n_blk_recv);
    for (int i = 0; i < ptp->n_part2; i++) {
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          int idx = ptp->recv_buffer_to_ref_lnum2[i][k];
          blk_recv_stride[idx] = _part2_stride[i][k];
        }
        for (int k = ptp->recv_buffer_to_duplicate_idx[i][j]; k < ptp->recv_buffer_to_duplicate_idx[i][j+1]; k++) {
          int idx      = ptp->recv_buffer_to_duplicate[i][2*k  ];
          int idx_data = ptp->recv_buffer_to_duplicate[i][2*k+1];
          blk_recv_stride[idx] = _part2_stride[i][idx_data];
        }
      }
    }
  }

  /*
   * Fill structure for asyncrhonous  exchange
   */
  ptp->async_send_s_data    [request_send] = s_data;
  ptp->async_send_cst_stride[request_send] = -1;
  ptp->async_send_tag       [_request    ] = -1;

  /*
   * Compute size of send / recv data
   */
  ptp->async_n_send_buffer[request_send]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer[request_send]    = malloc (sizeof(int) * (ptp->n_rank + 1));

  ptp->async_n_recv_buffer[request_recv]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[request_recv]    = malloc (sizeof(int) * (ptp->n_rank + 1));

  int* send_rank_n   = ptp->async_n_send_buffer[request_send];
  int* send_rank_idx = ptp->async_i_send_buffer[request_send];

  int* recv_rank_n   = ptp->async_n_recv_buffer[request_recv];
  int* recv_rank_idx = ptp->async_i_recv_buffer[request_recv];

  for (int i = 0; i < ptp->n_rank; i++) {
    send_rank_n[i] = 0;
    recv_rank_n[i] = 0;
  }

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    int beg =       ptp->default_i_send_buffer[dest];
    int end = beg + ptp->default_n_send_buffer[dest];
    for(int j = beg; j < end; ++j) {
      send_rank_n[dest] += blk_send_stride[j] * s_data;
    }
  }

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int dest = ptp->active_rank_recv[i];
    int beg =       ptp->default_i_recv_buffer[dest];
    int end = beg + ptp->default_n_recv_buffer[dest];
    for(int j = beg; j < end; ++j) {
      recv_rank_n[dest] += blk_recv_stride[j] * s_data;
    }
  }

  send_rank_idx[0] = 0;
  recv_rank_idx[0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    send_rank_idx[i+1] = send_rank_idx[i] + send_rank_n[i];
    recv_rank_idx[i+1] = recv_rank_idx[i] + recv_rank_n[i];
  }

  ptp->async_send_buffer[request_send] = malloc(sizeof(unsigned char) * send_rank_idx[ptp->n_rank] * s_data);
  ptp->async_recv_buffer[request_recv] = malloc(sizeof(unsigned char) * recv_rank_idx[ptp->n_rank] * s_data);

  unsigned char *send_buffer = ptp->async_send_buffer[request_send];
  unsigned char *recv_buffer = ptp->async_recv_buffer[request_recv];

  /*
   * Fill send buffer
   */
  int **part1_to_part2_data_idx = malloc(ptp->n_part1 * sizeof(int * ));
  for (int i = 0; i < ptp->n_part1; i++) {
    part1_to_part2_data_idx[i] = malloc((ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]+1) * sizeof(int));
    part1_to_part2_data_idx[i][0] = 0;
    for(int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      part1_to_part2_data_idx[i][j+1] = part1_to_part2_data_idx[i][j] + _part1_to_part2_stride[i][j];
    }
  }

  int *blk_send_idx = malloc ( (ptp->default_i_send_buffer[ptp->n_rank] + 1) * sizeof (int));
  blk_send_idx[0] = 0;
  for(int i = 0; i < ptp->default_i_send_buffer[ptp->n_rank]; ++i) {
    blk_send_idx[i+1] = blk_send_idx[i] + blk_send_stride[i];
  }

  for (int i = 0; i < ptp->n_part1; i++) {
    for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
               k < ptp->gnum1_to_send_buffer_idx[i][j+1];
               k++) {

        if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
          int idx_elmt = ptp->gnum1_to_send_buffer[i][k];
          int idx      = blk_send_idx      [idx_elmt] * (int) s_data;
          int idx1     = part1_to_part2_data_idx[i][j] * (int) s_data;

          int delta    = _part1_to_part2_stride[i][j] * s_data;

          // log_trace(" send at : (i=%i, j=%i / k=%i ) - idx_elmt = %i | idx = %i | idx1 = %i | delta = %i \n", i, j, k, idx_elmt, idx, idx1, delta);

          for (int k1 = 0; k1 < delta; k1++) {
            send_buffer[idx+k1] = _part1_to_part2_data[i][idx1+k1];
          }
        }
      }
    }
  }
  free(blk_send_idx);

  for (int i = 0; i < ptp->n_part1; i++) {
    free(part1_to_part2_data_idx[i]);
  }
  free(part1_to_part2_data_idx);

  /*
   *  Exchange data
   */
  // PDM_MPI_Ialltoallv(send_buffer,
  //                    send_rank_n,
  //                    send_rank_idx,
  //                    PDM_MPI_UNSIGNED_CHAR,
  //                    recv_buffer,
  //                    recv_rank_n,
  //                    recv_rank_idx,
  //                    PDM_MPI_UNSIGNED_CHAR,
  //                    ptp->comm,
  //                    &(ptp->async_alltoall_subrequest[3 * _request + 2]));
  PDM_MPI_Alltoallv(send_buffer,
                     send_rank_n,
                     send_rank_idx,
                     PDM_MPI_UNSIGNED_CHAR,
                     recv_buffer,
                     recv_rank_n,
                     recv_rank_idx,
                     PDM_MPI_UNSIGNED_CHAR,
                     ptp->comm);

  free(blk_send_stride);

  int *blk_recv_idx = malloc ( (ptp->default_i_recv_buffer[ptp->n_rank] + 1) * sizeof (int));
  blk_recv_idx[0] = 0;
  for(int i = 0; i < ptp->default_i_recv_buffer[ptp->n_rank]; ++i) {
    blk_recv_idx[i+1] = blk_recv_idx[i] + blk_recv_stride[i];
  }
  free(blk_recv_stride);

  // Keep recv stride for  post-treatment
  ptp->async_exch_part2_stride[_request] = _part2_stride;
  ptp->async_exch_recv_idx    [_request] = blk_recv_idx;

  ptp->async_recv_s_data    [request_recv] = s_data;
  ptp->async_recv_cst_stride[request_recv] = -1;
  ptp->async_recv_tag       [request_recv] = -1;

  /*
   * Prepare part2_stride/part2_data
   */
  if (stride2_unknown) {
    *part2_stride = _part2_stride;
  }

  unsigned char** _part2_data = malloc( ptp->n_part2 * sizeof(unsigned char *));
  for(int i = 0; i < ptp->n_part2; ++i) {
    int size = 0;
    for(int j = 0; j < ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]; ++j) {
      size += _part2_stride[i][j];
    }
    _part2_data[i] = malloc( size * s_data * sizeof(unsigned char));
  }

  ptp->async_recv_part2_data[request_recv]  = malloc(sizeof (void *) * ptp->n_part2);
  memcpy(ptp->async_recv_part2_data[request_recv], _part2_data, sizeof (void *) * ptp->n_part2);

  *part2_data = (void *) _part2_data;

  if (__part1_to_part2_stride != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {
      free (__part1_to_part2_stride[i]);
      free (__part1_to_part2_data[i]);
    }
    free (__part1_to_part2_stride);
    free (__part1_to_part2_data);
  }

}

static
void
_alltotall_stride_var_wait_and_post
(
 PDM_part_to_part_t                *ptp,
 int                 request
)
{

  // PDM_MPI_Wait (&(ptp->async_alltoall_subrequest[3 * request + 2]));

  int request_send = ptp->async_alltoall_subrequest[3 * request];
  int request_recv = ptp->async_alltoall_subrequest[3 * request + 1];

  size_t s_data  = ptp->async_recv_s_data[request_recv];

  int **part2_stri = ptp->async_exch_part2_stride[request];
  unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request_recv];

  int  *blk_recv_idx = ptp->async_exch_recv_idx[request];

  int **part2_idx = malloc(ptp->n_part2 * sizeof(int * ));
  for (int i = 0; i < ptp->n_part2; i++) {
    part2_idx[i] = malloc((ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]+1) * sizeof(int));
    part2_idx[i][0] = 0;
    for(int j = 0; j < ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]; j++) {
      part2_idx[i][j+1] = part2_idx[i][j] + part2_stri[i][j];
    }
  }

  for (int i = 0; i < ptp->n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {

        int idx_elmt = ptp->recv_buffer_to_ref_lnum2[i][k];
        int idx      = blk_recv_idx[idx_elmt] * (int) s_data;

        int idx1 = part2_idx[i][k] * (int) s_data;

        int delta = part2_stri[i][k] * s_data;

        // log_trace(" write at : (i=%i, j=%i / k=%i ) - idx_elmt = %i | idx = %i | idx1 = %i | delta = %i \n", i, j, k, idx_elmt, idx, idx1, delta);

        for (int k1 = 0; k1 < delta; k1++) {
          _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request_recv][idx+k1];
        }
      }
    }
  }

  free(ptp->async_exch_recv_idx[request]);
  ptp->async_exch_recv_idx[request] = NULL;

  _free_async_send (ptp, request_send);
  _free_async_recv (ptp, request_recv);
  _free_async_alltoall (ptp, request);

  for (int i = 0; i < ptp->n_part2; i++) {
    free(part2_idx[i]);
  }
  free(part2_idx);
}



static void
_p2p_stride_var_iexch
(
 PDM_part_to_part_t                *ptp,
 PDM_mpi_comm_kind_t                comm_kind,
 const int                          tag,
 const PDM_part_to_part_data_def_t  t_part1_data_def,
 const size_t                       s_data,
 const int                        **part1_stride,
 const void                       **part1_data,
 int                             ***part2_stride,
 void                            ***part2_data,
 int                                request
)
{
  PDM_UNUSED(comm_kind);

  int   **_part1_to_part2_stride  = (int  **) part1_stride;
  void  **_part1_to_part2_data    = (void **) part1_data;
  int   **__part1_to_part2_stride = NULL;
  void  **__part1_to_part2_data   = NULL;

  /*
   *  Create __part1_to_part2_stride and __part1_to_part2_data if necessary
   */

  if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {
    __part1_to_part2_stride = (int  **) malloc (sizeof (int  *) * ptp->n_part1);
    __part1_to_part2_data   = (void **) malloc (sizeof (void *) * ptp->n_part1);

    _part1_to_part2_stride = __part1_to_part2_stride;
    _part1_to_part2_data   = __part1_to_part2_data;

    for (int i = 0; i < ptp->n_part1; i++) {
      _part1_to_part2_stride[i] = malloc (sizeof(int) * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);
      size_t k = 0;
      size_t s_part_data = 0;
      for (int j = 0; j < ptp->n_elt1[i]; j++) {
        for (int j1 = ptp->part1_to_part2_idx[i][j]; j1 < ptp->part1_to_part2_idx[i][j+1]; j1++) {
          _part1_to_part2_stride[i][k++] = part1_stride[i][j];
          s_part_data += part1_stride[i][j];
        }
      }

      _part1_to_part2_data[i] = malloc (s_data * s_part_data);
      unsigned char *map_part1_to_part2_data = (unsigned char*) _part1_to_part2_data[i];

      int beg_data = 0;
      k = 0;
      for (int j = 0; j < ptp->n_elt1[i]; j++) {
        unsigned char *tmp_part1_data = (unsigned char*) (part1_data[i]) + beg_data;
        for (int j1 = ptp->part1_to_part2_idx[i][j]; j1 < ptp->part1_to_part2_idx[i][j+1]; j1++) {
          for (int j2 = 0; j2 < (int) (part1_stride[i][j] * s_data); j2++) {
            map_part1_to_part2_data[k++] = tmp_part1_data[j2];
          }
        }
        beg_data += part1_stride[i][j] * s_data;
      }
    }
  }


  /*
   *  Stride exchange if necessary
   */
  int stride2_unknown = ((*part2_stride) == NULL);

  int* send_n = malloc(ptp->n_rank * sizeof(int));
  for (int i = 0; i < ptp->n_rank; i++) {
    send_n[i] = 0;
  }
  int n_blk_send = 0;
  int *blk_send_stride = NULL;
  int *blk_send_idx    = NULL;

  int **_part2_stride = NULL;
  if (stride2_unknown) {
    int send_request_stri = -1;

    PDM_part_to_part_issend(ptp,
                            sizeof (int),
                            1, // Stride = 1
                            (const void **)  _part1_to_part2_stride,
                            tag,
                            &send_request_stri);


    n_blk_send = ptp->async_i_send_buffer[send_request_stri][ptp->n_rank]/sizeof(int);
    blk_send_stride = (int*) ptp->async_send_buffer[send_request_stri];
    blk_send_idx    = malloc( (n_blk_send + 1) * sizeof(int) );

    for (int i = 0; i < ptp->n_active_rank_send; i++) {
      int dest = ptp->active_rank_send[i];
      int beg =       ptp->async_i_send_buffer[send_request_stri][dest]/sizeof(int);
      int end = beg + ptp->async_n_send_buffer[send_request_stri][dest]/sizeof(int);
      for(int j = beg; j < end; ++j) {
        send_n[dest] += blk_send_stride[j];
      }
    }

    blk_send_idx[0] = 0;
    for(int i = 0; i < n_blk_send; ++i) {
      blk_send_idx[i+1] = blk_send_idx[i] + blk_send_stride[i];
    }

    _part2_stride = malloc( ptp->n_part2 * sizeof(int *));
    for(int i = 0; i < ptp->n_part2; ++i) {
      _part2_stride[i] = malloc( ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]] * sizeof(int));
    }

    int recv_request_stri = -1;
    PDM_part_to_part_irecv (ptp,
                            sizeof (int),
                            1,
                            (void **) _part2_stride,
                            tag,
                            &recv_request_stri);

    PDM_part_to_part_issend_wait(ptp, send_request_stri);

    _p2p_stride_var_irecv_stride_wait (ptp,
                                       &(ptp->async_exch_recv_n[request]),
                                       &(ptp->async_exch_recv_idx[request]),
                                       recv_request_stri);
  }
  else {
    n_blk_send = ptp->default_i_send_buffer[ptp->n_rank];

    blk_send_stride = malloc(sizeof(int) * n_blk_send);
    for (int i = 0; i < ptp->n_part1; i++) {
      for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
             k < ptp->gnum1_to_send_buffer_idx[i][j+1];
             k++) {
          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
            int idx = ptp->gnum1_to_send_buffer[i][k];
            blk_send_stride[idx] = _part1_to_part2_stride[i][k];
          }
        }
      }
    }

    for (int i = 0; i < ptp->n_active_rank_send; i++) {
      int dest = ptp->active_rank_send[i];
      int beg =       ptp->default_i_send_buffer[dest];
      int end = beg + ptp->default_n_send_buffer[dest];
      for (int j = beg; j < end; ++j) {
        send_n[dest] += blk_send_stride[j];
      }
    }

    blk_send_idx = PDM_array_new_idx_from_sizes_int(blk_send_stride,
                                                    n_blk_send);
    free(blk_send_stride);



    _part2_stride = *part2_stride;



    int n_blk_recv = ptp->default_i_recv_buffer[ptp->n_rank];
    int *blk_recv_stride = malloc(sizeof(int) * n_blk_recv);
    for (int i = 0; i < ptp->n_part2; i++) {
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          int idx = ptp->recv_buffer_to_ref_lnum2[i][k];
          blk_recv_stride[idx] = _part2_stride[i][k];
        }
        for (int k = ptp->recv_buffer_to_duplicate_idx[i][j]; k < ptp->recv_buffer_to_duplicate_idx[i][j+1]; k++) {
          int idx      = ptp->recv_buffer_to_duplicate[i][2*k  ];
          int idx_data = ptp->recv_buffer_to_duplicate[i][2*k+1];
          blk_recv_stride[idx] = _part2_stride[i][idx_data];
        }
      }
    }

    // PDM_log_trace_array_int(blk_recv_stride,
    //                         n_blk_recv,
    //                         "blk_recv_stride : ");

    ptp->async_exch_recv_n[request] = PDM_array_zeros_int(ptp->n_rank);
    for (int i = 0; i < ptp->n_active_rank_recv; i++) {
      int dest = ptp->active_rank_recv[i];
      int beg =       ptp->default_i_recv_buffer[dest];
      int end = beg + ptp->default_n_recv_buffer[dest];
      for (int j = beg; j < end; ++j) {
        ptp->async_exch_recv_n[request][dest] += blk_recv_stride[j];
      }
    }

    ptp->async_exch_recv_idx[request] = PDM_array_new_idx_from_sizes_int(blk_recv_stride,
                                                                         n_blk_recv);
    free(blk_recv_stride);
  }


  /*
   * Exchange data
   */
  int send_request_data = -1;

  _p2p_stride_var_data_issend (ptp,
                               tag,
                               s_data,
                               send_n,
                               blk_send_idx,
                               _part1_to_part2_stride,
                               _part1_to_part2_data,
                               &send_request_data);

  ptp->async_exch_subrequest[request][0] = send_request_data;

  free(send_n);
  free(blk_send_idx);


  unsigned char** _part2_data = malloc( ptp->n_part2 * sizeof(unsigned char *));
  for(int i = 0; i < ptp->n_part2; ++i) {
    int size = 0;
    for(int j = 0; j < ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]; ++j) {
      size += _part2_stride[i][j];
    }
    _part2_data[i] = malloc( size * s_data * sizeof(unsigned char));
  }

  if (stride2_unknown) {
    *part2_stride = _part2_stride;
  }
  *part2_data = (void *) _part2_data;


  ptp->async_exch_part2_stride[request] = _part2_stride;

  int recv_request_data = -1;

  _p2p_stride_var_data_irecv (ptp,
                              s_data,
                              ptp->async_exch_recv_n[request],
                              (void **) _part2_data,
                              tag,
                              &recv_request_data);

  ptp->async_exch_subrequest[request][1] = recv_request_data;


  if (__part1_to_part2_stride != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {
      free (__part1_to_part2_stride[i]);
      free (__part1_to_part2_data[i]);
    }
    free (__part1_to_part2_stride);
    free (__part1_to_part2_data);
  }
}



/**
 *
 * \brief Initialize an asynchronus exchange with p2p communication kind
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   tag              Tag for p2p exchange
 * \param [in]   t_part1_data_def Kind of part1 data definition
 * \param [in]   s_data           Data size
 * \param [in]   part1_stride     Stride of partition 1 data
 * \param [in]   part1_data       Partition 1 data
 * \param [out]  part2_stride     Stride of partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  part2_data       Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  request          Request
 *
 */

static void
_p2p_stride_var_reverse_iexch
(
 PDM_part_to_part_t                *ptp,
 const int                          tag,
 const PDM_part_to_part_data_def_t  t_part2_data_def,
 const size_t                       s_data,
 const int                        **part2_stride,
 const void                       **part2_data,
 int                             ***part1_stride,
 void                            ***part1_data,
 int                                request
)
{
  int   **_part2_to_part1_stride = (int **) part2_stride;
  void  **_part2_to_part1_data   = (void **) part2_data;
  int   **__part2_to_part1_stride = NULL;
  void  **__part2_to_part1_data   = NULL;

  /*
   *  Create __part2_to_part1_stride and __part2_to_part1_data if necessary
   */

  if (t_part2_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART2) {
    __part2_to_part1_stride = (int **) malloc (sizeof (int*) * ptp->n_part2);
    __part2_to_part1_data   = (void **) malloc (sizeof (void*) * ptp->n_part2);

    _part2_to_part1_stride = __part2_to_part1_stride;
    _part2_to_part1_data   = __part2_to_part1_data;

    for (int i = 0; i < ptp->n_part2; i++) {

      int *part2_idx = malloc(sizeof(int) * (ptp->n_elt2[i] + 1));
      part2_idx[0] = 0;

      for (int j = 0; j < ptp->n_elt2[i]; j++) {
        part2_idx[j+1] = part2_idx[j] + part2_stride[i][j] * s_data;
      }

      _part2_to_part1_stride[i] = malloc (sizeof(int) * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);

      int k = 0;
      int s = 0;

      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        int ielt2 = ptp->ref_lnum2[i][j] - 1;
        for (int k1 = ptp->gnum1_come_from_idx[i][j]; k1 < ptp->gnum1_come_from_idx[i][j+1]; k1++) {
          _part2_to_part1_stride[i][k++] = part2_stride[i][ielt2];
          s += part2_stride[i][ielt2];
        }
      }

      _part2_to_part1_data[i] = malloc (s_data * s);

      unsigned char *map_part2_to_part1_data = (unsigned char*) _part2_to_part1_data[i];
      unsigned char *map_part2_data = (unsigned char*) part2_data[i];

      k = 0;
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        int ielt2 = ptp->ref_lnum2[i][j] - 1;
        for (int k1 = ptp->gnum1_come_from_idx[i][j]; k1 < ptp->gnum1_come_from_idx[i][j+1]; k1++) {
          for (int k2 = part2_idx[ielt2]; k2 < part2_idx[ielt2+1]; k2++) {
            map_part2_to_part1_data[k++] = map_part2_data[k2];
          }
        }
      }

      free (part2_idx);
    }

  }

  /*
   *  Stride exchange if necessary
   */
  int stride1_unknown = ((*part1_stride) == NULL);

  int* send_n = malloc(ptp->n_rank * sizeof(int));
  for (int i = 0; i < ptp->n_rank; i++) {
    send_n[i] = 0;
  }
  int n_blk_send = 0;
  int *blk_send_stride = NULL;
  int *blk_send_idx    = NULL;

  int **_part1_stride = NULL;
  if (stride1_unknown) {
    int send_request_stri = -1;

    PDM_part_to_part_reverse_issend(ptp,
                                    sizeof (int),
                                    1, // Stride = 1
                                    (const void **)  _part2_to_part1_stride,
                                    tag,
                                    &send_request_stri);

    n_blk_send = ptp->async_i_send_buffer[send_request_stri][ptp->n_rank]/sizeof(int);
    blk_send_stride = (int*) ptp->async_send_buffer[send_request_stri];
    blk_send_idx    = malloc( (n_blk_send + 1) * sizeof(int) );

    for (int i = 0; i < ptp->n_active_rank_recv; i++) {
      int dest = ptp->active_rank_recv[i];
      int beg =       ptp->async_i_send_buffer[send_request_stri][dest]/sizeof(int);
      int end = beg + ptp->async_n_send_buffer[send_request_stri][dest]/sizeof(int);
      for(int j = beg; j < end; ++j) {
        send_n[dest] += blk_send_stride[j];
      }
    }

    blk_send_idx[0] = 0;
    for(int i = 0; i < n_blk_send; ++i) {
      blk_send_idx[i+1] = blk_send_idx[i] + blk_send_stride[i];
    }

    _part1_stride = malloc( ptp->n_part1 * sizeof(int*));
    for(int i = 0; i < ptp->n_part1;  ++i) {
      _part1_stride[i] = malloc(ptp->part1_to_part2_idx[i][ptp->n_elt1[i]] * sizeof(int));
    }

    int recv_request_stri = -1;
    PDM_part_to_part_reverse_irecv (ptp,
                                    sizeof (int),
                                    1,
                                    (void **) _part1_stride,
                                    tag,
                                    &recv_request_stri);

    PDM_part_to_part_reverse_issend_wait(ptp, send_request_stri);

    _p2p_stride_var_reverse_irecv_stride_wait (ptp,
                                               &(ptp->async_exch_recv_n[request]),
                                               &(ptp->async_exch_recv_idx[request]),
                                               recv_request_stri);
  }
  else {
    n_blk_send = ptp->default_i_recv_buffer[ptp->n_rank];

    blk_send_stride = malloc(sizeof(int) * n_blk_send);
    for (int i = 0; i < ptp->n_part2; i++) {
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          int idx = ptp->recv_buffer_to_ref_lnum2[i][k];
          blk_send_stride[idx] = _part2_to_part1_stride[i][k];
        }
        for (int k = ptp->recv_buffer_to_duplicate_idx[i][j]; k < ptp->recv_buffer_to_duplicate_idx[i][j+1]; k++) {
          int idx      = ptp->recv_buffer_to_duplicate[i][2*k  ];
          int idx_data = ptp->recv_buffer_to_duplicate[i][2*k+1];
          blk_send_stride[idx] = _part2_to_part1_stride[i][idx_data];
        }
      }
    }

    for (int i = 0; i < ptp->n_active_rank_recv; i++) {
      int dest = ptp->active_rank_recv[i];
      int beg =       ptp->default_i_recv_buffer[dest];
      int end = beg + ptp->default_n_recv_buffer[dest];
      for (int j = beg; j < end; ++j) {
        send_n[dest] += blk_send_stride[j];
      }
    }

    blk_send_idx = PDM_array_new_idx_from_sizes_int(blk_send_stride,
                                                    n_blk_send);
    free(blk_send_stride);



    _part1_stride = *part1_stride;



    int n_blk_recv = ptp->default_i_send_buffer[ptp->n_rank];
    int *blk_recv_stride = malloc(sizeof(int) * n_blk_recv);
    for (int i = 0; i < ptp->n_part1; i++) {
      for (int i1 = 0; i1 < ptp->n_elt1[i]; i1++) {
        for (int j = ptp->part1_to_part2_idx[i][i1]; j < ptp->part1_to_part2_idx[i][i1+1]; j++) {
          for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
               k < ptp->gnum1_to_send_buffer_idx[i][j+1];
               k++) {
            int idx = ptp->gnum1_to_send_buffer[i][k];
            blk_recv_stride[idx] = _part1_stride[i][k];
          }
        }
      }
    }


    ptp->async_exch_recv_n[request] = PDM_array_zeros_int(ptp->n_rank);
    for (int i = 0; i < ptp->n_active_rank_send; i++) {
      int dest = ptp->active_rank_send[i];
      int beg =       ptp->default_i_send_buffer[dest];
      int end = beg + ptp->default_n_send_buffer[dest];
      for (int j = beg; j < end; ++j) {
        ptp->async_exch_recv_n[request][dest] += blk_recv_stride[j];
      }
    }

    ptp->async_exch_recv_idx[request] = PDM_array_new_idx_from_sizes_int(blk_recv_stride,
                                                                         n_blk_recv);
    free(blk_recv_stride);
  }


  /*
   * Exchange data
   */

  int send_request_data = -1;

  _p2p_stride_var_data_reverse_issend (ptp,
                                       tag,
                                       s_data,
                                       send_n,
                                       blk_send_idx,
                                       _part2_to_part1_stride,
                                       _part2_to_part1_data,
                                       &send_request_data);

  ptp->async_exch_subrequest[request][0] = send_request_data;

  free(send_n);
  free(blk_send_idx);



  unsigned char** _part1_data = malloc( ptp->n_part1 * sizeof(unsigned char *));
  for(int i = 0; i < ptp->n_part1; ++i) {
    int size = 0;
    for(int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; ++j) {
      size += _part1_stride[i][j];
    }
    _part1_data[i] = malloc( size * s_data * sizeof(unsigned char));
  }

  *part1_stride =          _part1_stride;
  *part1_data   = (void *) _part1_data;

  ptp->async_exch_part2_stride[request] = _part1_stride;

  int recv_request_data = -1;

  _p2p_stride_var_data_reverse_irecv (ptp,
                                      s_data,
                                      ptp->async_exch_recv_n[request],
                                      (void **) _part1_data,
                                      tag,
                                      &recv_request_data);

  ptp->async_exch_subrequest[request][1] = recv_request_data;

  if (__part2_to_part1_stride != NULL) {
    for (int i = 0; i < ptp->n_part2; i++) {
      free (__part2_to_part1_stride[i]);
      free (__part2_to_part1_data[i]);
    }
    free (__part2_to_part1_stride);
    free (__part2_to_part1_data);

  }


}

/**
 *
 * \brief Wait the end of a data exchange with stride
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

static
void
_p2p_stride_var_iexch_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
)
{

  int  *blk_recv_idx = ptp->async_exch_recv_idx[request];
  int  request_irecv = ptp->async_exch_subrequest[request][1];

  int **part2_stri = ptp->async_exch_part2_stride[request];

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request_irecv][i]));
  }

  size_t s_data  = ptp->async_recv_s_data[request_irecv];

  //int cst_stride = ptp->async_recv_cst_stride[request];

  unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request_irecv];

  int **part2_idx = malloc(ptp->n_part2 * sizeof(int * ));
  for (int i = 0; i < ptp->n_part2; i++) {
    part2_idx[i] = malloc((ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]+1) * sizeof(int));
    part2_idx[i][0] = 0;
    for(int j = 0; j < ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]; j++) {
      part2_idx[i][j+1] = part2_idx[i][j] + part2_stri[i][j];
    }
  }

  // PDM_log_trace_array_int(ptp->async_i_recv_buffer[request], ptp->n_rank+1, "async_i_recv_buffer : ");
  // PDM_log_trace_array_int(ptp->async_n_recv_buffer[request], ptp->n_rank  , "async_n_recv_buffer : ");
  // double *buf = (double *) ptp->async_recv_buffer[request];
  // int size  = ptp->async_i_recv_buffer[request][ptp->n_rank]/(sizeof(double));
  // PDM_log_trace_array_double(buf, size, "buf : ");

  // int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {

        int idx_elmt = ptp->recv_buffer_to_ref_lnum2[i][k];
        int idx      = blk_recv_idx[idx_elmt] * (int) s_data;

        // int idx = ptp->recv_buffer_to_ref_lnum2[i][k] * delta;
        // int idx1 = k * delta;
        int idx1 = part2_idx[i][k] * (int) s_data;

        int delta = part2_stri[i][k] * s_data;

        // log_trace(" write at : (i=%i, j=%i / k=%i ) - idx_elmt = %i | idx = %i | idx1 = %i | delta = %i \n", i, j, k, idx_elmt, idx, idx1, delta);

        for (int k1 = 0; k1 < delta; k1++) {
          _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request_irecv][idx+k1];
        }
      }
    }
  }

  //_free_async_exch (ptp, request);
  free(ptp->async_recv_part2_data[request_irecv]);
  ptp->async_recv_part2_data[request_irecv] = NULL;

  for (int i = 0; i < ptp->n_part2; i++) {
    free(part2_idx[i]);
  }
  free(part2_idx);

}


/**
 *
 * \brief Wait the end of a data exchange with stride
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

static
void
_p2p_stride_var_reverse_iexch_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
)
{

  int  *blk_recv_idx = ptp->async_exch_recv_idx[request];
  int  request_irecv = ptp->async_exch_subrequest[request][1];

  int **part1_stri = ptp->async_exch_part2_stride[request];

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request_irecv][i]));
  }

  size_t s_data  = ptp->async_recv_s_data[request_irecv];

  //int cst_stride = ptp->async_recv_cst_stride[request];

  unsigned char ** _part1_data = (unsigned char **) ptp->async_recv_part2_data[request_irecv];

  int **part1_idx = malloc(ptp->n_part1 * sizeof(int * ));
  for (int i = 0; i < ptp->n_part1; i++) {
    part1_idx[i] = malloc((ptp->part1_to_part2_idx[i][ptp->n_elt1[i]] +1) * sizeof(int));
    part1_idx[i][0] = 0;
    for(int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      part1_idx[i][j+1] = part1_idx[i][j] + part1_stri[i][j];
    }
  }

  // PDM_log_trace_array_int(ptp->async_i_recv_buffer[request_irecv], ptp->n_rank+1, "async_i_recv_buffer : ");
  // PDM_log_trace_array_int(ptp->async_n_recv_buffer[request_irecv], ptp->n_rank  , "async_n_recv_buffer : ");
  // double *buf = (double *) ptp->async_recv_buffer[request_irecv];
  // int size  = ptp->async_i_recv_buffer[request_irecv][ptp->n_rank]/(sizeof(double));
  // PDM_log_trace_array_double(buf, size, "buf : ");

//  PDM_log_trace_array_double((double *) ptp->async_recv_buffer[request], ptp->async_i_recv_buffer[request][ptp->n_rank]/ (int) s_data, "ptp->async_recv_buffer[request] : ");

  for (int i = 0; i < ptp->n_part1; i++) {
    for (int i1 = 0; i1 < ptp->n_elt1[i]; i1++) {
      for (int j = ptp->part1_to_part2_idx[i][i1]; j < ptp->part1_to_part2_idx[i][i1+1]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
                 k < ptp->gnum1_to_send_buffer_idx[i][j+1];
                 k++) {

          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {

            int idx  = blk_recv_idx[ptp->gnum1_to_send_buffer[i][k]] * (int) s_data;
            int idx1 = part1_idx[i][j] * (int) s_data;

            int delta = part1_stri[i][j] * s_data;

            for (int k1 = 0; k1 < delta; k1++) {
              _part1_data[i][idx1+k1] = ptp->async_recv_buffer[request_irecv][idx+k1];
            }
          }
        }
      }
    }
  }

  // free(ptp->async_recv_part2_data[request_irecv]);
  // ptp->async_recv_part2_data[request_irecv] = NULL;

  _free_async_recv(ptp, request_irecv);
  //_free_async_exch (ptp, request);

  for (int i = 0; i < ptp->n_part1; i++) {
    free(part1_idx[i]);
  }
  free(part1_idx);

}

/**
 *
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1               Element global number (size : \ref n_part1)
 * \param [in]   n_elt1                  Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1                 Number of partition
 * \param [in]   gnum_elt2               Element global number (size : \ref n_part2)
 * \param [in]   n_elt2                  Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2                 Number of partition
 * \param [in]   part1_to_part2_idx      Index of data to send to gnum2 from gnum1
 *                                       (for each part size : \ref n_elt1+1)
 * \param [in]   part1_to_part2          Data to send to gnum2 from gnum1
 * \param [in]   part1_to_part2_triplet  Data to send to (irank2, ipart2, ielt2) from gnum1
 * \param [in]   comm                    MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

static
PDM_part_to_part_t *
_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const int           **part1_to_part2_triplet_idx,
 const PDM_g_num_t   **part1_to_part2,
 const int           **part1_to_part2_triplet,
 const int             from_triplet,
 const PDM_MPI_Comm    comm
)
{
  PDM_part_to_part_t *ptp = (PDM_part_to_part_t *) malloc (sizeof(PDM_part_to_part_t));

  char *env_var = NULL;
  env_var = getenv ("PDM_PART_TO_PART_USE_TAG");
  ptp->use_tag = 0;
  if (env_var != NULL) {
    ptp->use_tag = (int) atoi(env_var);
  }

  /* Init */
  if(ptp->use_tag == 1) {
    ptp->comm                     = comm;
  } else {
    PDM_MPI_Comm_dup(comm, &ptp->comm);
  }

  ptp->n_part1                  = n_part1;
  ptp->gnum_elt1                = gnum_elt1;

  ptp->n_part2                  = n_part2;
  ptp->gnum_elt2                = gnum_elt2;

  ptp->part1_to_part2           = part1_to_part2;

  ptp->part1_to_part2_triplet   = part1_to_part2_triplet;

  // assert( (gnum_elt2 == NULL) ? (part1_to_part2_triplet != NULL) : (part1_to_part2 != NULL));

  /* Copy */
  ptp->n_elt1 = malloc( ptp->n_part1 * sizeof(int));
  ptp->n_elt2 = malloc( ptp->n_part2 * sizeof(int));

  for(int i = 0; i < ptp->n_part1; ++i) {
    ptp->n_elt1[i] = n_elt1[i];
  }

  for(int i = 0; i < ptp->n_part2; ++i) {
    ptp->n_elt2[i] = n_elt2[i];
  }

  //
  // Pourquoi une copie ?

  ptp->part1_to_part2_idx = malloc( ptp->n_part1 * sizeof(int *));
  for(int i_part = 0; i_part < ptp->n_part1; ++i_part) {
    // log_trace("i_part1 %d, n_elt1 %d\n", i_part, ptp->n_elt1[i_part]);
    ptp->part1_to_part2_idx[i_part] = malloc((ptp->n_elt1[i_part] + 1) * sizeof(int));

    if (from_triplet != 1) {
      for(int i = 0; i < n_elt1[i_part]+1; ++i) {
        // log_trace("  i %d\n", i);
        ptp->part1_to_part2_idx[i_part][i] = part1_to_part2_idx[i_part][i];
      }
    }
    else {
      for(int i = 0; i < n_elt1[i_part]+1; ++i) {
        ptp->part1_to_part2_idx[i_part][i] = part1_to_part2_idx[i_part][i]/3;
      }
    }
  }

  if (0) {
    log_trace("--- Part1 ---\n");
    for (int i_part = 0; i_part < ptp->n_part1; i_part++) {
      log_trace("part %d: n_elt = %d\n", i_part, ptp->n_elt1[i_part]);
      PDM_log_trace_array_long(ptp->gnum_elt1[i_part], ptp->n_elt1[i_part], "gnum_elt : ");
      PDM_log_trace_connectivity_long(ptp->part1_to_part2_idx[i_part],
                                      ptp->part1_to_part2[i_part],
                                      ptp->n_elt1[i_part],
                                      "part1_to_part2 : ");
    }

    log_trace("\n\n--- Part2 ---\n");
    for (int i_part = 0; i_part < ptp->n_part2; i_part++) {
      log_trace("part %d: n_elt = %d\n", i_part, ptp->n_elt2[i_part]);
      PDM_log_trace_array_long(ptp->gnum_elt2[i_part], ptp->n_elt2[i_part], "gnum_elt : ");
    }
  }

  PDM_MPI_Comm_size (comm, &(ptp->n_rank));
  PDM_MPI_Comm_rank (comm, &(ptp->my_rank));
  int my_rank = ptp->my_rank;
  int n_rank = ptp->n_rank;

  int n_part2_max = 0;
  PDM_MPI_Allreduce ((int* )&n_part2, &n_part2_max, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  ptp->n_ref_lnum2                = NULL;
  ptp->ref_lnum2                  = NULL;
  ptp->n_unref_lnum2              = NULL;
  ptp->unref_lnum2                = NULL;
  ptp->gnum1_come_from_idx        = NULL;
  ptp->gnum1_come_from            = NULL;

  ptp->gnum1_to_send_buffer       = NULL;
  ptp->recv_buffer_to_ref_lnum2   = NULL;

  ptp->recv_buffer_to_duplicate_idx = NULL;
  ptp->recv_buffer_to_duplicate     = NULL;

  ptp->default_n_send_buffer      = malloc (sizeof(int) * n_rank);
  ptp->default_i_send_buffer      = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_send_buffer[0]   = 0;
  ptp->default_n_recv_buffer      = malloc (sizeof(int) * n_rank);
  ptp->default_i_recv_buffer      = malloc (sizeof(int) * (n_rank + 1));
  ptp->default_i_recv_buffer[0]   = 0;

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i] = 0;
    ptp->default_n_recv_buffer[i] = 0;
  }

  ptp->n_active_rank_send         = 1;
  ptp->active_rank_send           = NULL;

  ptp->n_active_rank_recv         = 1;
  ptp->active_rank_recv           = NULL;

  ptp->async_send_n_free          = 0;
  ptp->async_send_free            = NULL;
  ptp->async_send_l_array         = 0;
  ptp->async_send_s_data          = NULL;
  ptp->async_send_cst_stride      = NULL;
  ptp->async_send_tag             = NULL;
  ptp->async_send_request         = NULL;
  ptp->async_send_buffer          = NULL;
  ptp->async_n_send_buffer        = NULL;
  ptp->async_i_send_buffer        = NULL;

  ptp->async_recv_n_free          = 0;
  ptp->async_recv_free            = NULL;
  ptp->async_recv_l_array         = 0;
  ptp->async_recv_s_data          = NULL;
  ptp->async_recv_cst_stride      = NULL;
  ptp->async_recv_tag             = NULL;
  ptp->async_recv_request         = NULL;
  ptp->async_recv_buffer          = NULL;
  ptp->async_n_recv_buffer        = NULL;
  ptp->async_i_recv_buffer        = NULL;
  ptp->async_recv_part2_data      = NULL;

  ptp->async_alltoall_n_free  = 0;
  ptp->async_alltoall_l_array = 0;

  ptp->async_exch_n_free         = 0;
  ptp->async_exch_free           = NULL;
  ptp->async_exch_l_array        = 0;
  ptp->async_exch_subrequest     = NULL;
  ptp->async_exch_subrequest_s   = NULL;

  /* 1 - gnum_location in 2 1D array part1_to_part2_rank part1_to_part2_part   part1_to_part2_part elt*/

  PDM_gnum_location_t *gl = NULL;

  // if (gnum_elt2 != NULL || part1_to_part2_triplet == NULL ) {
  if (from_triplet == 0) {

    gl = PDM_gnum_location_create (n_part2, n_part1, comm, PDM_OWNERSHIP_KEEP);

    for (int i = 0; i < n_part2; i++) {
      // PDM_log_trace_array_long(gnum_elt2[i], n_elt2[i]  , "gnum_elt2::");
      PDM_gnum_location_elements_set (gl, i, n_elt2[i], gnum_elt2[i]);
    }

    for (int i = 0; i < n_part1; i++) {
      // PDM_log_trace_array_long(part1_to_part2[i], part1_to_part2_idx[i][n_elt1[i]]  , "part1_to_part2::");
      PDM_gnum_location_requested_elements_set (gl, i, part1_to_part2_idx[i][n_elt1[i]], part1_to_part2[i]);
    }

    PDM_gnum_location_compute(gl);
  }

  int n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {

    int *location_part1_to_part2_idx;
    int *location_part1_to_part2;

    if (from_triplet == 0) {
      PDM_gnum_location_get (gl,
                             i,
                             &location_part1_to_part2_idx,
                             &location_part1_to_part2);
      n_total_elt += location_part1_to_part2_idx[part1_to_part2_idx[i][n_elt1[i]]];
    }
    else if (from_triplet == 1) {
      location_part1_to_part2_idx = (int *) ptp->part1_to_part2_idx[i];
      location_part1_to_part2     = (int *)  part1_to_part2_triplet[i];

      n_total_elt += 3 * ptp->part1_to_part2_idx[i][n_elt1[i]];
      // n_total_elt += 3 * location_part1_to_part2_idx[ptp->part1_to_part2_idx[i][n_elt1[i]]];
    }
    else {
      n_total_elt += part1_to_part2_triplet_idx[i][part1_to_part2_idx[i][n_elt1[i]]];
    }

  }

  n_total_elt /= 3;

  int *merge_part1_to_part2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *merge_part1_to_part2_addr1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *order                      = (int *) malloc (sizeof(int) * n_total_elt);
  int *i_send_buffer              = (int *) malloc (sizeof(int) * n_total_elt);
  int *n_part1_to_part2_rank      = (int *) malloc (sizeof(int) * n_rank);
  int *idx_part1_to_part2_rank    = (int *) malloc (sizeof(int) * (n_rank + 1));

  for (int i = 0; i < n_rank; i++) {
    n_part1_to_part2_rank[i] = 0;
  }
  idx_part1_to_part2_rank[0] = 0;

  n_total_elt = 0;
  for (int i = 0; i < n_part1; i++) {

    int *location_part1_to_part2_idx;
    int *location_part1_to_part2;

    if (from_triplet == 0) {
      PDM_gnum_location_get (gl,
                             i,
                             &location_part1_to_part2_idx,
                             &location_part1_to_part2);
    }
    else if (from_triplet == 2) {
      location_part1_to_part2_idx = (int *) part1_to_part2_triplet_idx[i];
      location_part1_to_part2     = (int *) part1_to_part2_triplet[i];
    }

    if (from_triplet != 1) {
      for (int j = 0; j < n_elt1[i]; j++) {
        for (int k1 = part1_to_part2_idx[i][j]; k1 < part1_to_part2_idx[i][j+1]; k1++) {
          for (int k = location_part1_to_part2_idx[k1]/3;
                   k < location_part1_to_part2_idx[k1+1]/3; k++) {
            int i_rank2 = location_part1_to_part2[3*k];
            n_part1_to_part2_rank[i_rank2]++;
            merge_part1_to_part2_rank2[n_total_elt] = i_rank2;
            merge_part1_to_part2_part2[n_total_elt] = location_part1_to_part2[3*k+1];
            merge_part1_to_part2_lnum2[n_total_elt] = location_part1_to_part2[3*k+2];
            merge_part1_to_part2_rank1[n_total_elt] = my_rank;
            merge_part1_to_part2_part1[n_total_elt] = i;
            merge_part1_to_part2_lnum1[n_total_elt] = j;
            merge_part1_to_part2_addr1[n_total_elt] = k1;
            order[n_total_elt]                      = n_total_elt;
            n_total_elt++;
          }
        }
      }
    }

    else {

      location_part1_to_part2_idx = (int *) ptp->part1_to_part2_idx[i];
      location_part1_to_part2     = (int *) part1_to_part2_triplet[i];

      for (int j = 0; j < n_elt1[i]; j++) {
        for (int k = location_part1_to_part2_idx[j];
                 k < location_part1_to_part2_idx[j+1]; k++) {
          int i_rank2 = location_part1_to_part2[3*k];
          n_part1_to_part2_rank[i_rank2]++;
          merge_part1_to_part2_rank2[n_total_elt] = i_rank2;
          merge_part1_to_part2_part2[n_total_elt] = location_part1_to_part2[3*k+1];
          merge_part1_to_part2_lnum2[n_total_elt] = location_part1_to_part2[3*k+2];
          merge_part1_to_part2_rank1[n_total_elt] = my_rank;
          merge_part1_to_part2_part1[n_total_elt] = i;
          merge_part1_to_part2_lnum1[n_total_elt] = j;
          merge_part1_to_part2_addr1[n_total_elt] = k;
          order[n_total_elt]                      = n_total_elt;
          n_total_elt++;
        }
      }
    }
  }

  if (from_triplet == 0) {
    PDM_gnum_location_free (gl);
  }

  for (int i = 0; i < n_rank; i++) {
    idx_part1_to_part2_rank[i+1] = n_part1_to_part2_rank[i] +
                                   idx_part1_to_part2_rank[i];
  }

  /* 2 - Sort part1_to_part2_rank part1_to_part2_part according successively rank, part and elemt  */

  PDM_sort_int (merge_part1_to_part2_rank2, order, n_total_elt);

  int *_merge_part1_to_part2_rank2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_part2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_lnum2 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_rank1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_part1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_lnum1 = (int *) malloc (sizeof(int) * n_total_elt);
  int *_merge_part1_to_part2_addr1 = (int *) malloc (sizeof(int) * n_total_elt);

  for (int i = 0; i < n_total_elt; i++) {
    i_send_buffer[i] = -2;
    _merge_part1_to_part2_part2[i] = merge_part1_to_part2_part2[order[i]];
    _merge_part1_to_part2_lnum2[i] = merge_part1_to_part2_lnum2[order[i]];
    _merge_part1_to_part2_rank1[i] = merge_part1_to_part2_rank1[order[i]];
    _merge_part1_to_part2_part1[i] = merge_part1_to_part2_part1[order[i]];
    _merge_part1_to_part2_lnum1[i] = merge_part1_to_part2_lnum1[order[i]];
    _merge_part1_to_part2_addr1[i] = merge_part1_to_part2_addr1[order[i]];
  }

  int *_tmp_merge_part1_to_part2_rank2 = merge_part1_to_part2_rank2;
  int *_tmp_merge_part1_to_part2_part2 = merge_part1_to_part2_part2;
  int *_tmp_merge_part1_to_part2_lnum2 = merge_part1_to_part2_lnum2;
  int *_tmp_merge_part1_to_part2_rank1 = merge_part1_to_part2_rank1;
  int *_tmp_merge_part1_to_part2_part1 = merge_part1_to_part2_part1;
  int *_tmp_merge_part1_to_part2_lnum1 = merge_part1_to_part2_lnum1;
  int *_tmp_merge_part1_to_part2_addr1 = merge_part1_to_part2_addr1;

  merge_part1_to_part2_rank2 = _merge_part1_to_part2_rank2;
  merge_part1_to_part2_part2 = _merge_part1_to_part2_part2;
  merge_part1_to_part2_lnum2 = _merge_part1_to_part2_lnum2;
  merge_part1_to_part2_rank1 = _merge_part1_to_part2_rank1;
  merge_part1_to_part2_part1 = _merge_part1_to_part2_part1;
  merge_part1_to_part2_lnum1 = _merge_part1_to_part2_lnum1;
  merge_part1_to_part2_addr1 = _merge_part1_to_part2_addr1;

  _merge_part1_to_part2_rank2 = _tmp_merge_part1_to_part2_rank2;
  _merge_part1_to_part2_part2 = _tmp_merge_part1_to_part2_part2;
  _merge_part1_to_part2_lnum2 = _tmp_merge_part1_to_part2_lnum2;
  _merge_part1_to_part2_rank1 = _tmp_merge_part1_to_part2_rank1;
  _merge_part1_to_part2_part1 = _tmp_merge_part1_to_part2_part1;
  _merge_part1_to_part2_lnum1 = _tmp_merge_part1_to_part2_lnum1;
  _merge_part1_to_part2_addr1 = _tmp_merge_part1_to_part2_addr1;


  int *n_elt_part = malloc (sizeof(int) * n_part2_max);
  int *idx_elt_part = malloc (sizeof(int) * (n_part2_max + 1));
  idx_elt_part[0] = 0;

  int cpt_buff = -1;
  for (int i = 0; i < n_rank; i++) {
    int n_elt_rank = n_part1_to_part2_rank[i];
    for (int j = 0; j < n_elt_rank; j++) {
      order[j] = j;
    }

    for (int j = 0; j < n_part2_max; j++) {
      n_elt_part[j] = 0;
    }

    int *rank_i_send_buffer              = i_send_buffer              + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_part2 = merge_part1_to_part2_part2 + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_lnum2 = merge_part1_to_part2_lnum2 + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_rank1 = merge_part1_to_part2_rank1 + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_part1 = merge_part1_to_part2_part1 + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_lnum1 = merge_part1_to_part2_lnum1 + idx_part1_to_part2_rank[i];
    int *rank_merge_part1_to_part2_addr1 = merge_part1_to_part2_addr1 + idx_part1_to_part2_rank[i];

    int *_rank_merge_part1_to_part2_part2 = _merge_part1_to_part2_part2 + idx_part1_to_part2_rank[i];
    int *_rank_merge_part1_to_part2_lnum2 = _merge_part1_to_part2_lnum2 + idx_part1_to_part2_rank[i];
    int *_rank_merge_part1_to_part2_rank1 = _merge_part1_to_part2_rank1 + idx_part1_to_part2_rank[i];
    int *_rank_merge_part1_to_part2_part1 = _merge_part1_to_part2_part1 + idx_part1_to_part2_rank[i];
    int *_rank_merge_part1_to_part2_lnum1 = _merge_part1_to_part2_lnum1 + idx_part1_to_part2_rank[i];
    int *_rank_merge_part1_to_part2_addr1 = _merge_part1_to_part2_addr1 + idx_part1_to_part2_rank[i];

    PDM_sort_int (rank_merge_part1_to_part2_part2, order, n_elt_rank);

    int _max_part = 0;
    for (int k = 0; k < n_elt_rank; k++) {
      int i_part = rank_merge_part1_to_part2_part2[k];
      n_elt_part[i_part]++;
      _max_part = PDM_MAX (_max_part, i_part+1);
      _rank_merge_part1_to_part2_part2[k] = i_part;
      _rank_merge_part1_to_part2_lnum2[k] = rank_merge_part1_to_part2_lnum2[order[k]];
      _rank_merge_part1_to_part2_rank1[k] = rank_merge_part1_to_part2_rank1[order[k]];
      _rank_merge_part1_to_part2_part1[k] = rank_merge_part1_to_part2_part1[order[k]];
      _rank_merge_part1_to_part2_lnum1[k] = rank_merge_part1_to_part2_lnum1[order[k]];
      _rank_merge_part1_to_part2_addr1[k] = rank_merge_part1_to_part2_addr1[order[k]];
    }

    for (int k = 0; k < _max_part; k++) {
      idx_elt_part[k+1] = idx_elt_part[k] + n_elt_part[k];
    }

    for (int k1 = 0; k1 < _max_part; k1++) {

      int _n_elt_part = n_elt_part[k1];

      for (int j = 0; j < _n_elt_part; j++) {
        order[j] = j;
      }

      int *_part_rank_merge_part1_to_part2_part2 = _rank_merge_part1_to_part2_part2 + idx_elt_part[k1];
      int *_part_rank_merge_part1_to_part2_lnum2 = _rank_merge_part1_to_part2_lnum2 + idx_elt_part[k1];
      int *_part_rank_merge_part1_to_part2_rank1 = _rank_merge_part1_to_part2_rank1 + idx_elt_part[k1];
      int *_part_rank_merge_part1_to_part2_part1 = _rank_merge_part1_to_part2_part1 + idx_elt_part[k1];
      int *_part_rank_merge_part1_to_part2_lnum1 = _rank_merge_part1_to_part2_lnum1 + idx_elt_part[k1];
      int *_part_rank_merge_part1_to_part2_addr1 = _rank_merge_part1_to_part2_addr1 + idx_elt_part[k1];

      int *part_rank_merge_part1_to_part2_part2 = rank_merge_part1_to_part2_part2 + idx_elt_part[k1];
      int *part_rank_merge_part1_to_part2_lnum2 = rank_merge_part1_to_part2_lnum2 + idx_elt_part[k1];
      int *part_rank_merge_part1_to_part2_rank1 = rank_merge_part1_to_part2_rank1 + idx_elt_part[k1];
      int *part_rank_merge_part1_to_part2_part1 = rank_merge_part1_to_part2_part1 + idx_elt_part[k1];
      int *part_rank_merge_part1_to_part2_lnum1 = rank_merge_part1_to_part2_lnum1 + idx_elt_part[k1];
      int *part_rank_merge_part1_to_part2_addr1 = rank_merge_part1_to_part2_addr1 + idx_elt_part[k1];
      int *part_rank_i_send_buffer              = rank_i_send_buffer              + idx_elt_part[k1];

      PDM_sort_int (_part_rank_merge_part1_to_part2_lnum2, order, _n_elt_part);

      int pre_val = -1;
      for (int k2 = 0; k2 < _n_elt_part; k2++) {
        part_rank_merge_part1_to_part2_part2[k2] = _part_rank_merge_part1_to_part2_part2[k2];
        part_rank_merge_part1_to_part2_lnum2[k2] = _part_rank_merge_part1_to_part2_lnum2[k2];
        part_rank_merge_part1_to_part2_rank1[k2] = _part_rank_merge_part1_to_part2_rank1[order[k2]];
        part_rank_merge_part1_to_part2_part1[k2] = _part_rank_merge_part1_to_part2_part1[order[k2]];
        part_rank_merge_part1_to_part2_lnum1[k2] = _part_rank_merge_part1_to_part2_lnum1[order[k2]];
        part_rank_merge_part1_to_part2_addr1[k2] = _part_rank_merge_part1_to_part2_addr1[order[k2]];

        if (pre_val != part_rank_merge_part1_to_part2_lnum2[k2]) {
          cpt_buff++;
          ptp->default_n_send_buffer[i]++;
        }
        part_rank_i_send_buffer[k2] = cpt_buff;

      }

    }

  }

  free (idx_elt_part);

  /* 3 - Define Default_n_send_buffer and  Default_i_send_buffer */

  for (int i = 0; i < n_rank; i++) {
    ptp->default_i_send_buffer[i+1] = ptp->default_i_send_buffer[i] + ptp->default_n_send_buffer[i];
  }

  if (1 == 0) {
    printf ("ptp->default_i_send_buffer :");
    for (int i = 0; i < n_rank + 1; i++) {
      printf (" %d", ptp->default_i_send_buffer[i]);
    }
    printf("\n");

    printf ("ptp->default_n_send_buffer :");
    for (int i = 0; i < n_rank; i++) {
      printf (" %d", ptp->default_n_send_buffer[i]);
    }
    printf("\n");
  }

  /* 4 - Define gnum1_to_send_buffer */

  ptp->gnum1_to_send_buffer_idx = malloc (sizeof (int*) * n_part1);
  ptp->gnum1_to_send_buffer     = malloc (sizeof (int*) * n_part1);
  int **gnum1_to_send_buffer_n    = malloc (sizeof (int*) * n_part1);

  for (int i = 0; i < n_part1; i++) {
    ptp->gnum1_to_send_buffer_idx[i] = malloc (sizeof (int) * (ptp->part1_to_part2_idx[i][n_elt1[i]]+1));
    gnum1_to_send_buffer_n[i] = malloc (sizeof (int) * ptp->part1_to_part2_idx[i][n_elt1[i]]);

    ptp->gnum1_to_send_buffer_idx[i][0] = 0;
    for (int k = 0; k < ptp->part1_to_part2_idx[i][n_elt1[i]]; k++) {
      ptp->gnum1_to_send_buffer_idx[i][k+1] = 0;
      gnum1_to_send_buffer_n[i][k] = 0;
    }
  }

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1 = merge_part1_to_part2_part1[i];
    int addr1  = merge_part1_to_part2_addr1[i];
//    int lnum1  = merge_part1_to_part2_lnum1[i];
    ptp->gnum1_to_send_buffer_idx[ipart1][addr1+1]++;
  }


  for (int i = 0; i < n_part1; i++) {
    for (int k = 0; k < ptp->part1_to_part2_idx[i][n_elt1[i]]; k++) {
      ptp->gnum1_to_send_buffer_idx[i][k+1] += ptp->gnum1_to_send_buffer_idx[i][k] ;
    }
  }

  if (1 == 0) {
    for (int i = 0; i < n_part1; i++) {
      printf ("gnum1_to_send_buffer_idx 2 %d :", i);
      for (int k = 0; k < n_elt1[i] + 1; k++) {
        printf (" %d",ptp->gnum1_to_send_buffer_idx[i][k]);
      }
      printf("\n");
    }
  }

  for (int i = 0; i < n_part1; i++) {
    int size = ptp->gnum1_to_send_buffer_idx[i][ptp->part1_to_part2_idx[i][n_elt1[i]]];
    ptp->gnum1_to_send_buffer[i] = malloc (sizeof (int) * size);
    for (int k = 0; k < size; k++) {
      ptp->gnum1_to_send_buffer[i][k] = -1;
    }
  }

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1 = merge_part1_to_part2_part1[i];
    int addr1  = merge_part1_to_part2_addr1[i];
//    int lnum1  = merge_part1_to_part2_lnum1[i];
    int idx    = i_send_buffer[i];
    int idx2   = ptp->gnum1_to_send_buffer_idx[ipart1][addr1] +
                 gnum1_to_send_buffer_n[ipart1][addr1]++;
    ptp->gnum1_to_send_buffer[ipart1][idx2] = idx;
  }

  if (1 == 0) {
    for (int i = 0; i < n_part1; i++) {
      for (int j = 0; j < ptp->part1_to_part2_idx[i][n_elt1[i]]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
                 k < ptp->gnum1_to_send_buffer_idx[i][j+1];
                 k++) {
          printf ("2 - %d %d : %d\n", i, j, ptp->gnum1_to_send_buffer[i][k]);
        }
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    free (gnum1_to_send_buffer_n[i]);
  }
  free(gnum1_to_send_buffer_n);

  free (n_elt_part);
  free (order);

  free (_merge_part1_to_part2_rank2);
  free (_merge_part1_to_part2_part2);
  free (_merge_part1_to_part2_lnum2);
  free (_merge_part1_to_part2_rank1);
  free (_merge_part1_to_part2_part1);
  free (_merge_part1_to_part2_lnum1);
  free (_merge_part1_to_part2_addr1);

  free (idx_part1_to_part2_rank);
  free (n_part1_to_part2_rank);


  /* 5 - Define Default_n_recv_buffer and  Default_i_recv_buffer */

  PDM_MPI_Alltoall (ptp->default_n_send_buffer, 1, PDM_MPI_INT,
                    ptp->default_n_recv_buffer, 1, PDM_MPI_INT,
                    comm);

  for (int i = 0; i < n_rank; i++) {
    ptp->default_i_recv_buffer[i+1] = ptp->default_i_recv_buffer[i] + ptp->default_n_recv_buffer[i];
  }

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i  ] *= 4;
    ptp->default_i_send_buffer[i+1] *= 4;
    ptp->default_n_recv_buffer[i  ] *= 4;
    ptp->default_i_recv_buffer[i+1] *= 4;
  }

  int *int_s_buff = malloc (sizeof(int) * ptp->default_i_send_buffer[n_rank]);
  int *int_r_buff = malloc (sizeof(int) * ptp->default_i_recv_buffer[n_rank]);

  PDM_g_num_t *gnum_s_buff = malloc (sizeof(PDM_g_num_t) * ptp->default_i_send_buffer[n_rank]);
  PDM_g_num_t *gnum_r_buff = malloc (sizeof(PDM_g_num_t) * ptp->default_i_recv_buffer[n_rank]);

  for (int i = 0; i < n_total_elt; i++) {
    int ipart1      = merge_part1_to_part2_part1[i];
    int ielt1       = merge_part1_to_part2_lnum1[i];
    int ipart2      = merge_part1_to_part2_part2[i];
    int ielt2       = merge_part1_to_part2_lnum2[i];
    int idx         = i_send_buffer[i];
    if (idx >= 0) {
      int_s_buff[4 * idx    ] = ipart1;
      int_s_buff[4 * idx + 1] = ielt1;
      int_s_buff[4 * idx + 2] = ipart2;
      int_s_buff[4 * idx + 3] = ielt2;
      gnum_s_buff[idx]        = gnum_elt1[ipart1][ielt1];
    }
  }

  free (merge_part1_to_part2_rank2);
  free (merge_part1_to_part2_part2);
  free (merge_part1_to_part2_lnum2);
  free (merge_part1_to_part2_rank1);
  free (merge_part1_to_part2_part1);
  free (merge_part1_to_part2_lnum1);
  free (merge_part1_to_part2_addr1);
  free (i_send_buffer);

  PDM_MPI_Alltoallv (int_s_buff, ptp->default_n_send_buffer, ptp->default_i_send_buffer, PDM_MPI_INT,
                     int_r_buff, ptp->default_n_recv_buffer, ptp->default_i_recv_buffer, PDM_MPI_INT,
                     comm);

  for (int i = 0; i < n_rank; i++) {
    ptp->default_n_send_buffer[i  ] /= 4;
    ptp->default_i_send_buffer[i+1] /= 4;
    ptp->default_n_recv_buffer[i  ] /= 4;
    ptp->default_i_recv_buffer[i+1] /= 4;
  }

  PDM_MPI_Alltoallv (gnum_s_buff, ptp->default_n_send_buffer, ptp->default_i_send_buffer, PDM__PDM_MPI_G_NUM,
                     gnum_r_buff, ptp->default_n_recv_buffer, ptp->default_i_recv_buffer, PDM__PDM_MPI_G_NUM,
                     comm);

  /* 6 - Build the arrays for the reveived view */

  ptp->n_ref_lnum2                  = malloc (sizeof (int          ) * n_part2);
  ptp->ref_lnum2                    = malloc (sizeof (int         *) * n_part2);
  ptp->n_unref_lnum2                = malloc (sizeof (int          ) * n_part2);
  ptp->unref_lnum2                  = malloc (sizeof (int         *) * n_part2);
  ptp->gnum1_come_from_idx          = malloc (sizeof (int         *) * n_part2);
  ptp->gnum1_come_from              = malloc (sizeof (PDM_g_num_t *) * n_part2);
  ptp->recv_buffer_to_ref_lnum2     = malloc (sizeof (int         *) * n_part2);
  ptp->recv_buffer_to_duplicate_idx = malloc (sizeof (int         *) * n_part2);
  ptp->recv_buffer_to_duplicate     = malloc (sizeof (int         *) * n_part2);

  //ptp->gnum1_to_send_buffer     = NULL;
  //ptp->recv_buffer_to_ref_lnum2 = NULL;

  int **tag_elt2 = malloc (sizeof (int *) * n_part2);

  for (int i = 0; i < n_part2; i++) {
    tag_elt2[i] = malloc (sizeof (int) * n_elt2[i]);
    for (int j = 0; j < n_elt2[i]; j++) {
      tag_elt2[i][j] = 0;
    }
  }

  for (int i = 0; i < n_rank; i++) {
    //int iproc1 = i;
    for (int j = ptp->default_i_recv_buffer[i]; j < ptp->default_i_recv_buffer[i+1]; j++) {
      //int recv_ipart1 = int_r_buff[4 * j + 0];
      //int recv_ielt1  = int_r_buff[4 * j + 1];
      //int recv_gnum1  = gnum_r_buff[j];
      int recv_ipart2 = int_r_buff[4 * j + 2];
      int recv_ielt2  = int_r_buff[4 * j + 3];
      tag_elt2[recv_ipart2][recv_ielt2]++;
    }
  }

  int **ielt_to_ref = malloc (sizeof (int *) * n_part2);
  for (int i = 0; i < n_part2; i++) {
    ielt_to_ref[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->n_ref_lnum2[i]   = 0;
    ptp->n_unref_lnum2[i] = 0;
    ptp->ref_lnum2[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->unref_lnum2[i] = malloc (sizeof (int) * n_elt2[i]);
    ptp->gnum1_come_from_idx[i] = malloc (sizeof (int) * (n_elt2[i] + 1));
    ptp->gnum1_come_from_idx[i][0] = 0;

    for (int j = 0; j < n_elt2[i]; j++) {
      int _tag = tag_elt2[i][j];
      tag_elt2[i][j] = 0;
      if (_tag > 0) {
        ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]+1] = _tag;
        ielt_to_ref[i][j] = ptp->n_ref_lnum2[i];
        ptp->ref_lnum2[i][ptp->n_ref_lnum2[i]++] = j+1;
      }
      else {
        ptp->unref_lnum2[i][ptp->n_unref_lnum2[i]++] = j+1;
      }
    }

    ptp->ref_lnum2[i]           = realloc (ptp->ref_lnum2[i], sizeof (int) * ptp->n_ref_lnum2[i]);
    ptp->unref_lnum2[i]         = realloc (ptp->unref_lnum2[i], sizeof (int) * ptp->n_unref_lnum2[i]);
    ptp->gnum1_come_from_idx[i] = realloc (ptp->gnum1_come_from_idx[i], sizeof (int) * (ptp->n_ref_lnum2[i] + 1));

    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      ptp->gnum1_come_from_idx[i][j+1] += ptp->gnum1_come_from_idx[i][j];
    }

    ptp->gnum1_come_from[i]          = malloc (sizeof (PDM_g_num_t) * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);
    ptp->recv_buffer_to_ref_lnum2[i] = malloc (sizeof (int) * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);
  }

  cpt_buff = 0;
  for (int i = 0; i < n_rank; i++) {
    //int iproc1 = i;
    for (int j = ptp->default_i_recv_buffer[i]; j < ptp->default_i_recv_buffer[i+1]; j++) {
      //int recv_ipart1 = int_r_buff[4 * j + 0];
      //int recv_ielt1  = int_r_buff[4 * j + 1];
      int recv_gnum1  = gnum_r_buff[j];
      int recv_ipart2 = int_r_buff[4 * j + 2];
      int recv_ielt2  = int_r_buff[4 * j + 3];
      int iref = ielt_to_ref[recv_ipart2][recv_ielt2];
      int idx = ptp->gnum1_come_from_idx[recv_ipart2][iref] + tag_elt2[recv_ipart2][iref];
      ptp->gnum1_come_from[recv_ipart2][idx] = recv_gnum1;
      ptp->recv_buffer_to_ref_lnum2[recv_ipart2][idx] = cpt_buff++;
      tag_elt2[recv_ipart2][iref]++;
    }
  }

  /* 8 - Sort and remove duplicate */

  int max_n_gnum1 = 0;
  for (int i = 0; i < n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      int idx = ptp->gnum1_come_from_idx[i][j];
      int n_gnum1 = ptp->gnum1_come_from_idx[i][j+1] - idx;
      max_n_gnum1 = PDM_MAX(max_n_gnum1, n_gnum1);
    }
  }

  if (1 == 0) {
    printf ("ptp->gnum1_come_from 1\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          printf (" "PDM_FMT_G_NUM"", ptp->gnum1_come_from[i][k]);
        }
        printf ("\n");
      }
    }

    printf ("ptp->recv_buffer_to_ref_lnum2 1\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          printf (" %d", ptp->recv_buffer_to_ref_lnum2[i][k]);
        }
        printf ("\n");
      }
    }
  }

  order = (int *) malloc (sizeof(int) * max_n_gnum1);

  for (int i = 0; i < n_part2; i++) {
    int *_recv_buffer_to_ref_lnum2 = malloc(sizeof(int) * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);

    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      int idx = ptp->gnum1_come_from_idx[i][j];
      int n_gnum1 = ptp->gnum1_come_from_idx[i][j+1] - idx;
      PDM_g_num_t *_gnum1_come_from = ptp->gnum1_come_from[i] + idx;
      int *__recv_buffer_to_ref_lnum2 = _recv_buffer_to_ref_lnum2 + idx;
      int *___recv_buffer_to_ref_lnum2 = ptp->recv_buffer_to_ref_lnum2[i] + idx;

      for (int k = 0; k < n_gnum1; k++) {
        order[k] = k;
      }

      PDM_sort_long (_gnum1_come_from, order, n_gnum1);

      for (int k = 0; k < n_gnum1; k++) {
        __recv_buffer_to_ref_lnum2[k] = ___recv_buffer_to_ref_lnum2[order[k]];
      }
    }

    free (ptp->recv_buffer_to_ref_lnum2[i]);
    ptp->recv_buffer_to_ref_lnum2[i] = _recv_buffer_to_ref_lnum2;

    int *_old_gnum1_come_from_idx = malloc(sizeof(int) * (ptp->n_ref_lnum2[i] + 1));

    for (int j = 0; j < ptp->n_ref_lnum2[i] + 1; j++) {
      _old_gnum1_come_from_idx[j] = ptp->gnum1_come_from_idx[i][j];

      ptp->gnum1_come_from_idx[i][j] = 0;
    }

    ptp->recv_buffer_to_duplicate_idx[i] = malloc (sizeof (int) * (ptp->n_ref_lnum2[i]+1));
    ptp->recv_buffer_to_duplicate_idx[i][0] = 0;
    ptp->recv_buffer_to_duplicate[i] = malloc(sizeof(int) * 2 * _old_gnum1_come_from_idx[ptp->n_ref_lnum2[i]]);

    int cpt = 0;
    int cpt1 = 0;

    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {

      PDM_g_num_t current_val = -1;
      for (int k = _old_gnum1_come_from_idx[j]; k < _old_gnum1_come_from_idx[j+1]; k++) {
        if (ptp->gnum1_come_from[i][k] != current_val) {
          current_val = ptp->gnum1_come_from[i][k];
          ptp->gnum1_come_from[i][cpt] = ptp->gnum1_come_from[i][k];
          ptp->recv_buffer_to_ref_lnum2[i][cpt] = ptp->recv_buffer_to_ref_lnum2[i][k];
          cpt++;
        }
        else {
          ptp->recv_buffer_to_duplicate[i][2*cpt1] = ptp->recv_buffer_to_ref_lnum2[i][k];
          ptp->recv_buffer_to_duplicate[i][2*cpt1+1] = cpt-1;
          cpt1++;
        }
      }
      ptp->gnum1_come_from_idx[i][j+1] = cpt;
      ptp->recv_buffer_to_duplicate_idx[i][j+1] = cpt1;
    }

    ptp->gnum1_come_from[i]          = realloc (ptp->gnum1_come_from[i], cpt * sizeof(PDM_g_num_t));
    ptp->recv_buffer_to_ref_lnum2[i] = realloc (ptp->recv_buffer_to_ref_lnum2[i], cpt * sizeof(int));
    ptp->recv_buffer_to_duplicate[i] = realloc (ptp->recv_buffer_to_duplicate[i], sizeof(int) * 2 * cpt1);

    free (_old_gnum1_come_from_idx);

  }

  if (1 == 0) {
    printf ("ptp->gnum1_come_from 2\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          printf (" "PDM_FMT_G_NUM"", ptp->gnum1_come_from[i][k]);
        }
        printf ("\n");
      }
    }

    printf ("ptp->recv_buffer_to_ref_lnum2 2\n");
    for (int i = 0; i < n_part2; i++) {
      printf(" - ipart : %d\n", i);
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        printf ("%d :", j);
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          printf (" %d", ptp->recv_buffer_to_ref_lnum2[i][k]);
        }
        printf ("\n");
      }
    }
  }

  /* 7 - Look for the active ranks */

  ptp->n_active_rank_send = 0;
  ptp->active_rank_send = malloc (sizeof(int) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    if (ptp->default_n_send_buffer[i] > 0) {
      ptp->active_rank_send[ptp->n_active_rank_send++] = i;
    }
  }

  ptp->n_active_rank_recv = 0;
  ptp->active_rank_recv = malloc (sizeof(int) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    if (ptp->default_n_recv_buffer[i] > 0) {
      ptp->active_rank_recv[ptp->n_active_rank_recv++] = i;
    }
  }

  free (order);

  free (int_s_buff);
  free (int_r_buff);

  free (gnum_s_buff);
  free (gnum_r_buff);

  for (int i = 0; i < n_part2; i++) {
    free (ielt_to_ref[i]);
    free (tag_elt2[i]);
  }
  free (tag_elt2);
  free (ielt_to_ref);

  // Create tag for P2P
  void  *max_tag_tmp;
  int    flag = 0;

  if(ptp->use_tag == 1) {
    // Mandatory to call with PDM_MPI_COMM_WORLD becuase only this one keep attributes (openMPI implemntation for exemple)
    PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_COMM_WORLD, &max_tag_tmp, &flag);
    ptp->max_tag  = (long) (*((int *) max_tag_tmp));
    ptp->seed_tag = PDM_MPI_Rand_tag(comm);
    ptp->next_tag = 1;
  } else {
    // Mandatory to call with PDM_MPI_COMM_WORLD becuase only this one keep attributes (openMPI implemntation for exemple)
    PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_COMM_WORLD, &max_tag_tmp, &flag);
    ptp->max_tag  = (long) (*((int *) max_tag_tmp));
    ptp->seed_tag = 1;
    ptp->next_tag = 1;
  }

  return ptp;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1          Element global number (size : \ref n_part1)
 * \param [in]   n_elt1             Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1            Number of partition
 * \param [in]   gnum_elt2          Element global number (size : \ref n_part2)
 * \param [in]   n_elt2             Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2            Number of partition
 * \param [in]   part1_to_part2_idx Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 * \param [in]   part1_to_part2     Data to send to gnum2 from gnum1
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const PDM_g_num_t   **part1_to_part2,
 const PDM_MPI_Comm    comm
)
{
  int from_triplet = 0;
  return _create (gnum_elt1,
                  n_elt1,
                  n_part1,
                  gnum_elt2,
                  n_elt2,
                  n_part2,
                  part1_to_part2_idx,
                  NULL,
                  part1_to_part2,
                  NULL,
                  from_triplet,
                  comm);
}

/**
 *
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1                   Element global number (size : \ref n_part1)
 * \param [in]   n_elt1                      Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1                     Number of partition
 * \param [in]   n_elt2                      Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2                     Number of partition
 * \param [in]   part1_to_part2_idx          Index of data to send to gnum2 from gnum1
 *                                           (for each part size : \ref n_elt1+1)
 * \param [in]   part1_to_part2_triplet_idx  (for each part size : \ref part1_to_part2_idx[\ref n_elt] + 1)
 * \param [in]   part1_to_part2_triplet      Data to send to (irank2, ipart2, ielt2) from gnum1
 * \param [in]   comm                        MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create_from_num2_triplet
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const int            *n_elt2,
 const int             n_part2,
 const int           **part1_to_part2_idx,
 const int           **part1_to_part2_triplet_idx,
 const int           **part1_to_part2_triplet,
 const PDM_MPI_Comm    comm
)
{
  int from_triplet = 2;
  if (part1_to_part2_triplet_idx == NULL) {
    from_triplet = 1;
  }

  return _create (gnum_elt1,
                  n_elt1,
                  n_part1,
                  NULL,
                  n_elt2,
                  n_part2,
                  part1_to_part2_idx,
                  part1_to_part2_triplet_idx,
                  NULL,
                  part1_to_part2_triplet,
                  from_triplet,
                  comm);
}


/**
 *
 * \brief Get selected numbers of part2 index
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 */

void
PDM_part_to_part_part1_to_part2_idx_get
(
 PDM_part_to_part_t *ptp,
 int               **n_elt1,
 int              ***part1_to_part2_idx
)
{
  *n_elt1             = (int          *) ptp->n_elt1;
  *part1_to_part2_idx = (int         **) ptp->part1_to_part2_idx;
}


/**
 *
 * \brief Get selected numbers of part2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 * \param [out]  part1_to_part2      Data to send to gnum2 from gnum1 for each part
 *
 */

void
PDM_part_to_part_part1_to_part2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_elt1,
 int              ***part1_to_part2_idx,
 PDM_g_num_t      ***part1_to_part2
)
{
  *n_elt1             = (int          *) ptp->n_elt1;
  *part1_to_part2_idx = (int         **) ptp->part1_to_part2_idx;
  *part1_to_part2     = (PDM_g_num_t **) ptp->part1_to_part2;
}


/**
 *
 * \brief Initialize an exchange based on MPI_ialltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data in same order than part1_to_part2 array
 * \param [out]  ref_part2_data      Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ialltoall
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part1_to_part2_data,
 void              **ref_part2_data,
 int                *request
)
{

  unsigned char ** _part1_data = (unsigned char **) part1_to_part2_data;

  *request         = _find_open_async_alltoall_exch (ptp);
  int request_send = _find_open_async_send_exch (ptp);
  int request_recv = _find_open_async_recv_exch (ptp);

  int _request = *request;
  ptp->async_alltoall_subrequest[3 * _request]     = request_send;
  ptp->async_alltoall_subrequest[3 * _request + 1] = request_recv;

  ptp->async_send_s_data[request_send]      = s_data;
  ptp->async_send_cst_stride[request_send]  = cst_stride;
  ptp->async_send_tag[_request]             = -1;
  ptp->async_send_request[_request]         = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_send_buffer[_request]        = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer[_request]        = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer[_request][0]     = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_send_buffer[_request][i]   = cst_stride * ptp->default_n_send_buffer[i] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = cst_stride * ptp->default_i_send_buffer[i+1] * (int) s_data;
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part1; i++) {
    for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
               k < ptp->gnum1_to_send_buffer_idx[i][j+1];
               k++) {

        if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
          int idx = ptp->gnum1_to_send_buffer[i][k] * delta;
          int idx1 = j* delta;
          for (int k1 = 0; k1 < delta; k1++) {
            ptp->async_send_buffer[_request][idx+k1] = _part1_data[i][idx1+k1];
          }
        }
      }
    }
  }

  ptp->async_recv_s_data[request_recv]      = s_data;
  ptp->async_recv_cst_stride[request_recv]  = cst_stride;
  ptp->async_recv_tag[request_recv]         = -1;

  ptp->async_recv_part2_data[request_recv]  = malloc(sizeof (void *) * ptp->n_part2);
  memcpy(ptp->async_recv_part2_data[request_recv], ref_part2_data, sizeof (void *) * ptp->n_part2);

  ptp->async_recv_request[request_recv]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_recv_buffer[request_recv]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[request_recv]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[request_recv][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[request_recv][i]   = cst_stride * ptp->default_n_recv_buffer[i] * (int) s_data;
    ptp->async_i_recv_buffer[request_recv][i+1] = cst_stride * ptp->default_i_recv_buffer[i+1] * (int) s_data;
  }
  ptp->async_recv_buffer[request_recv]      = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[request_recv][ptp->n_rank]);

  PDM_MPI_Ialltoallv (ptp->async_send_buffer[_request], ptp->async_n_send_buffer[_request], ptp->async_i_send_buffer[_request], PDM_MPI_UNSIGNED_CHAR,
                      ptp->async_recv_buffer[_request], ptp->async_n_recv_buffer[_request], ptp->async_i_recv_buffer[_request], PDM_MPI_UNSIGNED_CHAR,
                      ptp->comm, &(ptp->async_alltoall_subrequest[3 * _request + 2]));

}


/**
 *
 * \brief Wait a asynchronus ialltoall
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_ialltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
)
{

  PDM_MPI_Wait (&(ptp->async_alltoall_subrequest[3 * request + 2]));

  int request_send = ptp->async_alltoall_subrequest[3 * request];
  int request_recv = ptp->async_alltoall_subrequest[3 * request + 1];

  size_t s_data  = ptp->async_recv_s_data[request_recv];
  int cst_stride = ptp->async_recv_cst_stride[request_recv];

  unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request_recv];

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part2; i++) {
    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
        int idx = ptp->recv_buffer_to_ref_lnum2[i][k] * delta;
        int idx1 = k* delta;
        for (int k1 = 0; k1 < delta; k1++) {
          _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request_recv][idx+k1];
        }
      }
    }
  }

  _free_async_send (ptp, request_send);
  _free_async_recv (ptp, request_recv);
  _free_async_alltoall (ptp, request);

}


/**
 *
 * \brief Initialize an exchange based on MPI_ineighbor_alltoall
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data in same order than part1_to_part2 array
 * \param [out]  ref_part2_data          Data to referenced part2 elements
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_ineighbor_alltoall
(
PDM_part_to_part_t *ptp,
 const size_t       s_data,
 const int          cst_stride,
 void             **part1_to_part2_data,
 void             **ref_part2_data,
 int               *request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (s_data);
  PDM_UNUSED (cst_stride);
  PDM_UNUSED (part1_to_part2_data);
  PDM_UNUSED (ref_part2_data);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part_to_part_ineighbor_alltoall not yet implemented\n");
}

/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_ineighbor_alltoall_wait
(
 PDM_part_to_part_t *ptp,
 int                 request
)
{
  PDM_UNUSED (ptp);
  PDM_UNUSED (request);
  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_part_to_part_ineighbor_alltoall_wait not yet implemented\n");
}


/**
 *
 * \brief Get referenced Part2 elements
 *
 * \param [in]   ptp           Part-to-Part structure
 * \param [out]  n_ref_lnum2   Number of referenced Part2 elements
 * \param [out]  ref_lnum2     Referenced Part2 elements (one-based local ids)
 *
 */

void
PDM_part_to_part_ref_lnum2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_ref_lnum2,
 int              ***ref_lnum2
)
{
  *n_ref_lnum2 = ptp->n_ref_lnum2;
  *ref_lnum2   = ptp->ref_lnum2;
}


/**
 *
 * \brief Get unreferenced Part2 elements
 *
 * \param [in]   ptp             Part-to-Part structure
 * \param [out]  n_unref_lnum2   Number of referenced Part2 elements
 * \param [out]  unref_lnum2     Unreferenced Part2 elements (one-based local ids)
 *
 */

void
PDM_part_to_part_unref_lnum2_get
(
 PDM_part_to_part_t *ptp,
 int               **n_unref_lnum2,
 int              ***unref_lnum2
)
{
  *n_unref_lnum2 = ptp->n_unref_lnum2;
  *unref_lnum2   = ptp->unref_lnum2;
}


/**
 *
 * \brief Get gnum come from gnum1 for each referenced gnum2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [out]  gnum1_come_from_idx Index for gnum1_come_from array (size = \ref n_part2)
 * \param [out]  gnum1_come_from     Gnum come from gnum1 for each referenced gnum2
 *
 */

void
PDM_part_to_part_gnum1_come_from_get
(
 PDM_part_to_part_t *ptp,
 int              ***gnum1_come_from_idx,
 PDM_g_num_t      ***gnum1_come_from
)
{
  *gnum1_come_from_idx = ptp->gnum1_come_from_idx;
  *gnum1_come_from     = ptp->gnum1_come_from;
}


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data (order given by part1_to_part2 array)
 * \param [in]   tag                 Tag of the exchange
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part1_to_part2_data,
 const int           tag,
 int                *request
)
{
  unsigned char ** _part1_data = (unsigned char **) part1_to_part2_data;

  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data[_request]      = s_data;
  ptp->async_send_cst_stride[_request]  = cst_stride;
  ptp->async_send_tag[_request]         = tag;
  ptp->async_send_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_send_buffer[_request]  = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer[_request]  = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_send_buffer[_request][i]   = cst_stride * ptp->default_n_send_buffer[i] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = cst_stride * ptp->default_i_send_buffer[i+1] * (int) s_data;
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);

  // copy part1_to_part2 to send_buffer
  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part1; i++) {
    for (int j = 0; j < ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]; j++) {
      for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
               k < ptp->gnum1_to_send_buffer_idx[i][j+1];
               k++) {

        if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
          int idx = ptp->gnum1_to_send_buffer[i][k] * delta;
          int idx1 = j* delta;
          for (int k1 = 0; k1 < delta; k1++) {
            ptp->async_send_buffer[_request][idx+k1] = _part1_data[i][idx1+k1];
          }
        }
      }
    }
  }

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    unsigned char *buf =  ptp->async_send_buffer[_request] + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
                    tag, ptp->comm, &(ptp->async_send_request[_request][i]));
  }
}


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part1_to_part2_data Data (order given by part1_to_part2 array)
 * \param [in]   tag                 Tag of the exchange
 * \param [out]  request             Request
 *
 */
void
PDM_part_to_part_issend_raw
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void         *raw_buffer,
 const int           tag,
 int                *request
)
{
  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data    [_request] = s_data;
  ptp->async_send_cst_stride[_request] = cst_stride;
  ptp->async_send_tag       [_request] = tag;
  ptp->async_send_request   [_request] = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_send_buffer  [_request] = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer  [_request] = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer  [_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_send_buffer[_request][i  ] = cst_stride * ptp->default_n_send_buffer[i  ] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = cst_stride * ptp->default_i_send_buffer[i+1] * (int) s_data;
  }
  ptp->async_send_buffer[_request] = NULL;

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int dest = ptp->active_rank_send[i];
    unsigned char *buf = (unsigned char *) raw_buffer + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
                    tag, ptp->comm, &(ptp->async_send_request[_request][i]));
  }
}

/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  // log_trace("PDM_part_to_part_issend_wait = %i \n", request);
  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    PDM_MPI_Wait (&(ptp->async_send_request[request][i]));
  }

  _free_async_send (ptp, request);

}


/**
 *
 * \brief Initialize an asynchronus reverse issend (part2 to part1)
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   s_data              Data size
 * \param [in]   cst_stride          Constant stride
 * \param [in]   part2_to_part1_data Data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [in]   tag                 Tag of the exchange
 * \param [out]  request             Request
 *
 */

void
PDM_part_to_part_reverse_issend
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 const void        **part2_to_part1_data,
 const int           tag,
 int                *request
)
{
  unsigned char ** _part2_data = (unsigned char **) part2_to_part1_data;

  *request = _find_open_async_send_exch (ptp);
  int _request = *request;

  ptp->async_send_s_data[_request]      = s_data;
  ptp->async_send_cst_stride[_request]  = cst_stride;
  ptp->async_send_tag[_request]         = tag;
  ptp->async_send_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_send_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_send_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_send_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_send_buffer[_request][i]   = cst_stride * ptp->default_n_recv_buffer[i] * (int) s_data;
    ptp->async_i_send_buffer[_request][i+1] = cst_stride * ptp->default_i_recv_buffer[i+1] * (int) s_data;
  }
  ptp->async_send_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_send_buffer[_request][ptp->n_rank]);

  int delta = (int) s_data * cst_stride;
  for (int i = 0; i < ptp->n_part2; i++) {

    for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
      for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
        int idx = ptp->recv_buffer_to_ref_lnum2[i][k] * delta;
        int idx1 = k* delta;
        for (int k1 = 0; k1 < delta; k1++) {
          ptp->async_send_buffer[_request][idx+k1] = _part2_data[i][idx1+k1];
        }
      }
      for (int k = ptp->recv_buffer_to_duplicate_idx[i][j]; k < ptp->recv_buffer_to_duplicate_idx[i][j+1]; k++) {
        int idx      = ptp->recv_buffer_to_duplicate[i][2*k  ] * delta;
        int idx_data = ptp->recv_buffer_to_duplicate[i][2*k+1] * delta;
        for (int k1 = 0; k1 < delta; k1++) {
          ptp->async_send_buffer[_request][idx+k1] = _part2_data[i][idx_data+k1];
        }
      }

    }
  }

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int dest = ptp->active_rank_recv[i];
    unsigned char *buf =  ptp->async_send_buffer[_request] + ptp->async_i_send_buffer[_request][dest];
    int count = ptp->async_n_send_buffer[_request][dest];
    PDM_MPI_Issend (buf, count, PDM_MPI_UNSIGNED_CHAR, dest,
                    tag, ptp->comm, &(ptp->async_send_request[_request][i]));
  }
}


/**
 *
 * \brief Wait an asynchronus reverse issend (part2 to part1)
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  tag           Tag of the exchange
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_issend_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_send_request[request][i]));
  }

  _free_async_send (ptp, request);

}


/**
 *
 * \brief Wait an asynchronus reverse issend (part2 to part1)
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  tag           Tag of the exchange
 * \param [in]  request       Request
 *
 */
int
PDM_part_to_part_reverse_issend_test
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int flag = -1;
    PDM_MPI_Test (&(ptp->async_send_request[request][i]), &flag);
    if(flag == 0) {
      return 0; // If one of message is not OK we return immediatly
    }
  }
  return 1;
}


void
PDM_part_to_part_reverse_issend_post
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  _free_async_send(ptp, request);
}

/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part2_data    Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [in]  tag           Tag of the exchange
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part2_data,
 const int           tag,
 int                *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;
  // log_trace("PDM_part_to_part_irecv = %i | tag = %i \n", _request, tag);

  ptp->async_recv_s_data[_request]      = s_data;
  ptp->async_recv_cst_stride[_request]  = cst_stride;
  ptp->async_recv_tag[_request]         = tag;

  ptp->async_recv_part2_data[_request]  = malloc(sizeof (void *) * ptp->n_part2);
  memcpy(ptp->async_recv_part2_data[_request], part2_data, sizeof (void *) * ptp->n_part2);

  ptp->async_recv_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i]   = cst_stride * ptp->default_n_recv_buffer[i] * (int) s_data;
    ptp->async_i_recv_buffer[_request][i+1] = cst_stride * ptp->default_i_recv_buffer[i+1] * (int) s_data;
  }
  ptp->async_recv_buffer[_request]      = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[_request][ptp->n_rank]);

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int source = ptp->active_rank_recv[i];
    unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source,
                    tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }

}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part2_data    Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [in]  tag           Tag of the exchange
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_irecv_raw
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void               *raw_buffer,
 const int           tag,
 int                *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;

  ptp->async_recv_s_data    [_request] = s_data;
  ptp->async_recv_cst_stride[_request] = cst_stride;
  ptp->async_recv_tag       [_request] = tag;

  ptp->async_recv_request [_request]    = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_recv);
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;
  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i  ] = cst_stride * ptp->default_n_recv_buffer[i  ] * (int) s_data;
    ptp->async_i_recv_buffer[_request][i+1] = cst_stride * ptp->default_i_recv_buffer[i+1] * (int) s_data;
  }
  ptp->async_recv_buffer[_request] = NULL;

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int source = ptp->active_rank_recv[i];
    unsigned char *buf = (unsigned char *) raw_buffer + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source,
                    tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }

}


/**
 *
 * \brief Test the reception/send completion
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_issend_test
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  // log_trace("PDM_part_to_part_issend_test = %i \n", request);
  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int flag = -1;
    PDM_MPI_Test (&(ptp->async_send_request[request][i]), &flag);
    if(flag == 0) {
      return 0; // If one of message is not OK we return immediatly
    }
  }
  return 1;
}

/**
 *
 * \brief Test the reception/send completion
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_irecv_test
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  // log_trace("PDM_part_to_part_irecv_test = %i \n", request);
  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    int flag = -1;
    PDM_MPI_Test (&(ptp->async_recv_request[request][i]), &flag);
    if(flag == 0) {
      return 0; // If one of message is not OK we return immediatly
    }
  }
  return 1;
}

void
PDM_part_to_part_issend_post
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  // log_trace("PDM_part_to_part_issend_post = %i \n", request);
  _free_async_send (ptp, request);
}

void
PDM_part_to_part_irecv_post
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  size_t s_data  = ptp->async_recv_s_data[request];
  int cst_stride = ptp->async_recv_cst_stride[request];

  if(ptp->async_recv_part2_data[request] != NULL) {

    unsigned char ** _part2_data = (unsigned char **) ptp->async_recv_part2_data[request];

    int delta = (int) s_data * cst_stride;
    for (int i = 0; i < ptp->n_part2; i++) {
      for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
        for (int k = ptp->gnum1_come_from_idx[i][j]; k < ptp->gnum1_come_from_idx[i][j+1]; k++) {
          int idx = ptp->recv_buffer_to_ref_lnum2[i][k] * delta;
          int idx1 = k* delta;
          for (int k1 = 0; k1 < delta; k1++) {
            _part2_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1];
          }
        }
      }
    }
    free(ptp->async_recv_part2_data[request]);
    ptp->async_recv_part2_data[request] = NULL;
  }

  _free_async_recv (ptp, request);
}

/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{


  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }

  PDM_part_to_part_irecv_post(ptp, request);

}


/**
 *
 * \brief Wait a asynchronus raw irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_irecv_wait_raw
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  for (int i = 0; i < ptp->n_active_rank_recv; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }

  _free_async_recv (ptp, request);

}


/**
 *
 * \brief Initialize a asynchronus reverse irecv (from part2)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part1_data    Partition 1 data (order given by part1_to_part2 array)
 * \param [in]  tag           Tag of the exchange
 * \param [out] request       Request
 *
 */

void
PDM_part_to_part_reverse_irecv
(
 PDM_part_to_part_t *ptp,
 const size_t        s_data,
 const int           cst_stride,
 void              **part1_data,
 const int           tag,
 int                *request
)
{

  *request = _find_open_async_recv_exch (ptp);
  int _request = *request;

  ptp->async_recv_s_data[_request]      = s_data;
  ptp->async_recv_cst_stride[_request]  = cst_stride;
  ptp->async_recv_tag[_request]         = tag;

  ptp->async_recv_part2_data[_request]  = malloc(sizeof (void *) * ptp->n_part1);
  memcpy(ptp->async_recv_part2_data[_request], part1_data, sizeof (void *) * ptp->n_part1);

  ptp->async_recv_request[_request]     = malloc (sizeof (PDM_MPI_Request) * ptp->n_active_rank_send);
  ptp->async_n_recv_buffer[_request]    = malloc (sizeof(int) * ptp->n_rank);
  ptp->async_i_recv_buffer[_request]    = malloc (sizeof(int) * (ptp->n_rank + 1));
  ptp->async_i_recv_buffer[_request][0] = 0;

  for (int i = 0; i < ptp->n_rank; i++) {
    ptp->async_n_recv_buffer[_request][i]   = cst_stride * ptp->default_n_send_buffer[i]   * (int) s_data;
    ptp->async_i_recv_buffer[_request][i+1] = cst_stride * ptp->default_i_send_buffer[i+1] * (int) s_data;
  }
  ptp->async_recv_buffer[_request] = malloc (sizeof (unsigned char) * ptp->async_i_recv_buffer[_request][ptp->n_rank]);

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int source = ptp->active_rank_send[i];
    unsigned char *buf =  ptp->async_recv_buffer[_request] + ptp->async_i_recv_buffer[_request][source];
    int count = ptp->async_n_recv_buffer[_request][source];
    PDM_MPI_Irecv (buf, count, PDM_MPI_UNSIGNED_CHAR, source,
                    tag, ptp->comm, &(ptp->async_recv_request[_request][i]));
  }
}


/**
 *
 * \brief Initialize a asynchronus reverse irecv (from part2)
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  tag           Tag of the exchange
 * \param [in]  request       Request
 *
 */

void
PDM_part_to_part_reverse_irecv_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    PDM_MPI_Wait (&(ptp->async_recv_request[request][i]));
  }
  PDM_part_to_part_reverse_irecv_post(ptp, request);

  // size_t s_data  = ptp->async_recv_s_data[request];
  // int cst_stride = ptp->async_recv_cst_stride[request];

  // unsigned char ** _part1_data = (unsigned char **) ptp->async_recv_part2_data[request];

  // int delta = (int) s_data * cst_stride;

  // for (int i = 0; i < ptp->n_part1; i++) {
  //   for (int i1 = 0; i1 < ptp->n_elt1[i]; i1++) {
  //     for (int j = ptp->part1_to_part2_idx[i][i1]; j < ptp->part1_to_part2_idx[i][i1+1]; j++) {
  //       for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
  //                k < ptp->gnum1_to_send_buffer_idx[i][j+1];
  //                k++) {

  //         if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
  //           int idx  = ptp->gnum1_to_send_buffer[i][k] * delta;
  //           int idx1 = j * delta;
  //           for (int k1 = 0; k1 < delta; k1++) {
  //             _part1_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1];
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // _free_async_recv (ptp, request);
  // free(ptp->async_recv_part2_data[request]);

}


/**
 *
 * \brief Test the reception/send completion
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

int
PDM_part_to_part_reverse_irecv_test
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  for (int i = 0; i < ptp->n_active_rank_send; i++) {
    int flag = -1;
    PDM_MPI_Test (&(ptp->async_send_request[request][i]), &flag);
    if(flag == 0) {
      return 0; // If one of message is not OK we return immediatly
    }
  }
  return 1;
}


void
PDM_part_to_part_reverse_irecv_post
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{
  size_t s_data  = ptp->async_recv_s_data    [request];
  int cst_stride = ptp->async_recv_cst_stride[request];

  unsigned char ** _part1_data = (unsigned char **) ptp->async_recv_part2_data[request];

  int delta = (int) s_data * cst_stride;

  for (int i = 0; i < ptp->n_part1; i++) {
    for (int i1 = 0; i1 < ptp->n_elt1[i]; i1++) {
      for (int j = ptp->part1_to_part2_idx[i][i1]; j < ptp->part1_to_part2_idx[i][i1+1]; j++) {
        for (int k = ptp->gnum1_to_send_buffer_idx[i][j];
                 k < ptp->gnum1_to_send_buffer_idx[i][j+1];
                 k++) {

          if (ptp->gnum1_to_send_buffer[i][k] >= 0) {
            int idx  = ptp->gnum1_to_send_buffer[i][k] * delta;
            int idx1 = j * delta;
            for (int k1 = 0; k1 < delta; k1++) {
              _part1_data[i][idx1+k1] = ptp->async_recv_buffer[request][idx+k1];
            }
          }
        }
      }
    }
  }

  _free_async_recv (ptp, request);
  free(ptp->async_recv_part2_data[request]);
}

/**
 *
 * \brief Initialize a partial asynchronus exchange
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   k_commm          Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part1_data_def Kind of part1 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part1_stride     Stride of partition 1 data
 * \param [in]   part1_data       Partition 1 data
 * \param [out]  part2_stride     Stride of partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  part2_data       Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm,
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part1_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part1_stride,
 const void                       **part1_data,
 int                             ***part2_stride,
 void                            ***part2_data,
 int                               *request
)
{
  int tag = -10000;
  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    tag  = ptp->seed_tag;
    tag += (ptp->next_tag++);
    tag %= ptp->max_tag;
  }

  PDM_UNUSED (cst_stride);

  *request = _find_open_async_exch (ptp);

  int _request = *request;

  ptp->async_exch_t_stride[_request] = t_stride;
  ptp->async_exch_k_comm[_request]   = k_comm;

  if (t_stride == PDM_STRIDE_CST_INTERLEAVED) {

    ptp->async_exch_subrequest_s[_request] = cst_stride;
    ptp->async_exch_subrequest[_request] = realloc (ptp->async_exch_subrequest[_request], sizeof(int) * 2 * cst_stride);
    for (int i = 0; i < 2*cst_stride; i++) {
      ptp->async_exch_subrequest[_request][i] = -1;
    }

    void  ** __part1_to_part2_data = (void **) malloc (sizeof (void*) * ptp->n_part1);
    void  ** _part1_to_part2_data  = __part1_to_part2_data;

    if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {

      for (int i = 0; i < ptp->n_part1; i++) {

        _part1_to_part2_data[i] = malloc (s_data * cst_stride * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);

      }
    }

    *part2_data = malloc(sizeof(void *) * ptp->n_part2);
    void **_part2_data = *part2_data;
    for (int i = 0; i < ptp->n_part2; i++) {
      _part2_data[i] = malloc(s_data * cst_stride * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);
    }

    unsigned char **___part2_data = malloc(sizeof(unsigned char*) * ptp->n_part2);
    for (int i = 0; i < ptp->n_part2; i++) {
      ___part2_data[i] = (unsigned char *) _part2_data[i];
    }

    for (int i = 0; i < cst_stride; i++) {

      if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {

        for (int ipart = 0; ipart < ptp->n_part1; ipart++) {

          unsigned char *map_part1_to_part2_data = (unsigned char*) _part1_to_part2_data[ipart];
          unsigned char *map_part1_data = (unsigned char*) part1_data[ipart] + i * (s_data * ptp->n_elt1[ipart]);

          int k = 0;
          for (int j = 0; j < ptp->n_elt1[ipart]; j++) {
            for (int k1 = ptp->part1_to_part2_idx[ipart][j]; k1 < ptp->part1_to_part2_idx[ipart][j+1]; k1++) {
              for (int k2 = 0; k2 <  (int) s_data; k2++) {
                map_part1_to_part2_data[k++] = map_part1_data[j * (int) s_data + k2];
              }
            }
          }
        }
      }

      else {
        for (int i1 = 0; i1 < ptp->n_part1; i1++) {
          _part1_to_part2_data[i1] = (void *) ((unsigned char *) part1_data[i1] + i * ptp->part1_to_part2_idx[i1][ptp->n_elt1[i1]] * (int) s_data);
        }
      }

      for (int i1 = 0; i1 < ptp->n_part2; i1++) {
        ___part2_data[i1] = ((unsigned char *) _part2_data[i1]) + s_data * i * ptp->gnum1_come_from_idx[i1][ptp->n_ref_lnum2[i1]];
      }

      if (k_comm == PDM_MPI_COMM_KIND_P2P) {
        PDM_part_to_part_issend (ptp,
                                 s_data,
                                 1,
                  (const void **)_part1_to_part2_data,
                                 tag,
                                 &(ptp->async_exch_subrequest[_request][2*i]));

        PDM_part_to_part_irecv (ptp,
                                s_data,
                                1,
                        (void **)___part2_data,
                                tag,
                                &(ptp->async_exch_subrequest[_request][2*i+1]));
      }

      else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

        PDM_part_to_part_ialltoall (ptp,
                                    s_data,
                                    1,
                                    _part1_to_part2_data,
                            (void **)___part2_data,
                                    &(ptp->async_exch_subrequest[_request][2*i]));

      }

      else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

        printf ("Error PDM_part_to_part_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

        printf ("Error PDM_part_to_part_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

        printf ("Error PDM_part_to_part_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

        printf ("Error PDM_part_to_part_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

        printf ("Error PDM_part_to_part_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
        abort();

      }



    }

    if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {
      for (int i = 0; i < ptp->n_part1; i++) {
        free (__part1_to_part2_data[i]);
      }
    }
    free(__part1_to_part2_data);
    free(___part2_data);
  } 

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    void  **_part1_to_part2_data    = (void **) part1_data;
    void  **__part1_to_part2_data   = NULL;

    /*
     *  Create __part1_to_part2_stride and __part1_to_part2_data if necessary
     */

    if (t_part1_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART1) {
      __part1_to_part2_data  = (void **) malloc (sizeof (void*) * ptp->n_part1);

      _part1_to_part2_data   = __part1_to_part2_data;

      for (int i = 0; i < ptp->n_part1; i++) {

        _part1_to_part2_data[i] = malloc (s_data * cst_stride * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);
        unsigned char *map_part1_to_part2_data = (unsigned char*) _part1_to_part2_data[i];

        int k        = 0;
        int beg_data = 0;
        for (int j = 0; j < ptp->n_elt1[i]; j++) {
          unsigned char *tmp_part1_data = (unsigned char*) (part1_data[i]) + beg_data;
          for (int j1 = ptp->part1_to_part2_idx[i][j]; j1 < ptp->part1_to_part2_idx[i][j+1]; j1++) {
            for (int j2 = 0; j2 < (int) (cst_stride * s_data); j2++) {
              map_part1_to_part2_data[k++] = tmp_part1_data[j2];
            }
          }
          beg_data += cst_stride * s_data;
        }
      }
    }

    *part2_data = malloc(sizeof(void *) * ptp->n_part2);
    void **_part2_data = *part2_data;
    for (int i = 0; i < ptp->n_part2; i++) {
      _part2_data[i] = malloc(s_data * cst_stride * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);
    }

    if (k_comm == PDM_MPI_COMM_KIND_P2P) {


      PDM_part_to_part_issend (ptp,
                               s_data,
                               cst_stride,
                (const void **)_part1_to_part2_data,
                               tag,
                               &(ptp->async_exch_subrequest[_request][0]));

      PDM_part_to_part_irecv (ptp,
                              s_data,
                              cst_stride,
                             *part2_data,
                              tag,
                              &(ptp->async_exch_subrequest[_request][1]));

    }

    else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

      PDM_part_to_part_ialltoall (ptp,
                                  s_data,
                                  cst_stride,
                                  _part1_to_part2_data,
                          (void **)*part2_data,
                                  &(ptp->async_exch_subrequest[_request][0]));

    }

    else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }


    else {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " unknown k_comm is not implemented yet\n");
      abort();

    }

    if (__part1_to_part2_data != NULL) {
      for (int i = 0; i < ptp->n_part1; i++) {
        free (__part1_to_part2_data[i]);
      }
      free (__part1_to_part2_data);
      __part1_to_part2_data = NULL;
    }
  }

  else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    if (k_comm == PDM_MPI_COMM_KIND_P2P) {

      _p2p_stride_var_iexch (ptp,
                             PDM_MPI_COMM_KIND_P2P,
                             tag,
                             t_part1_data_def,
                             s_data,
                             part1_stride,
                             part1_data,
                             part2_stride,
                             part2_data,
                             _request);

    }

    else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

      _alltotall_stride_var_iexch(ptp,
                                  t_part1_data_def,
                                  s_data,
                                  part1_stride,
                                  part1_data,
                                  part2_stride,
                                  part2_data,
                                  &(ptp->async_exch_subrequest[_request][0]));


    }

    else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

  }

}


/**
 *
 * \brief Wait a partial asynchronus exchange
 *
 * \param [in]  ptp      Part to part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_iexch_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  if (ptp->async_exch_t_stride[request] == PDM_STRIDE_CST_INTERLEAVED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      for (int i = 0; i < ptp->async_exch_subrequest_s[request]; i++) {

        PDM_part_to_part_irecv_wait (ptp, ptp->async_exch_subrequest[request][2*i+1]);

        PDM_part_to_part_issend_wait (ptp, ptp->async_exch_subrequest[request][2*i]);

      }

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      for (int i = 0; i < ptp->async_exch_subrequest_s[request]; i++) {

        PDM_part_to_part_ialltoall_wait(ptp, ptp->async_exch_subrequest[request][2*i]);

      }

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

  }

  else if (ptp->async_exch_t_stride[request] == PDM_STRIDE_CST_INTERLACED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      PDM_part_to_part_irecv_wait (ptp, ptp->async_exch_subrequest[request][1]);

      PDM_part_to_part_issend_wait (ptp, ptp->async_exch_subrequest[request][0]);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      PDM_part_to_part_ialltoall_wait(ptp, ptp->async_exch_subrequest[request][0]);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

    else {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " unknown k_comm is not implemented yet\n");
      abort();

    }

  }

  else if (ptp->async_exch_t_stride[request] == PDM_STRIDE_VAR_INTERLACED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      PDM_part_to_part_issend_wait(ptp, ptp->async_exch_subrequest[request][0]);

      _p2p_stride_var_iexch_wait (ptp, request);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      _alltotall_stride_var_wait_and_post(ptp, ptp->async_exch_subrequest[request][0]);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

  }

  _free_async_exch (ptp, request);
}


/**
 *
 * \brief Initialize a partial reverse asynchronus exchange
 *
 * \param [in]   ptp              Part to part structure
 * \param [in]   k_commm          Kind of MPI communication
 * \param [in]   t_stride         Kind of stride
 * \param [in]   t_part2_data_def Kind of part2 data definition
 * \param [in]   cst_stride       Constant stride
 * \param [in]   s_data           Data size
 * \param [in]   part2_stride     Stride of partition 1 data (Accordding to t_part2_data_def)
 * \param [in]   part2_data       Partition 1 data (Accordding to t_part2_data_def)
 * \param [out]  part1_stride     Stride of partition 2 data (order given by part1_to_part2)
 * \param [out]  part1_data       Partition 2 data (order given by part1_to_part2)
 * \param [out]  request          Request
 *
 */

void
PDM_part_to_part_reverse_iexch
(
 PDM_part_to_part_t                *ptp,
 const PDM_mpi_comm_kind_t          k_comm,
 const PDM_stride_t                 t_stride,
 const PDM_part_to_part_data_def_t  t_part2_data_def,
 const int                          cst_stride,
 const size_t                       s_data,
 const int                        **part2_stride,
 const void                       **part2_data,
 int                             ***part1_stride,
 void                            ***part1_data,
 int                               *request
)
{
  int tag = -10000;
  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    tag  = ptp->seed_tag;
    tag += (ptp->next_tag++);
    tag %= ptp->max_tag;
  }

  *request = _find_open_async_exch (ptp);

  int _request = *request;

  ptp->async_exch_t_stride[_request] = t_stride;
  ptp->async_exch_k_comm[_request]   = k_comm;

  if (t_stride == PDM_STRIDE_CST_INTERLEAVED) {

    ptp->async_exch_subrequest_s[_request] = cst_stride;
    ptp->async_exch_subrequest[_request] = realloc (ptp->async_exch_subrequest[_request], sizeof(int) * 2 * cst_stride);
    for (int i = 0; i < 2*cst_stride; i++) {
      ptp->async_exch_subrequest[_request][i] = -1;
    }

    void **__part2_to_part1_data = (void **) malloc (sizeof (void*) * ptp->n_part2);
    void **_part2_to_part1_data  = __part2_to_part1_data;

    if (t_part2_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART2) {

      for (int i = 0; i < ptp->n_part2; i++) {

        _part2_to_part1_data[i] = malloc (s_data * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);

      }
    }

    *part1_data = malloc(sizeof(void *) * ptp->n_part1);

    void **_part1_data = *part1_data;

    for (int i = 0; i < ptp->n_part1; i++) {
      _part1_data[i] = malloc(s_data * cst_stride * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);
    }

    unsigned char **___part1_data = malloc(sizeof(unsigned char*) * ptp->n_part1);
    for (int i = 0; i < ptp->n_part1; i++) {
      ___part1_data[i] = (unsigned char *) _part1_data[i];
    }

    for (int i = 0; i < cst_stride; i++) {

      if (t_part2_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART2) {

        for (int i1 = 0; i1 < ptp->n_part2; i1++) {

          unsigned char *map_part2_to_part1_data = (unsigned char*) _part2_to_part1_data[i];
          unsigned char *map_part2_data = (unsigned char*) part2_data[i] + i * (cst_stride * s_data * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]]);

          int k = 0;
          for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
            int ielt2 = ptp->ref_lnum2[i][j] - 1;
            for (int k1 = ptp->gnum1_come_from_idx[i][j]; k1 < ptp->gnum1_come_from_idx[i][j+1]; k1++) {
              for (int k2 = 0; k2 <  (int) s_data; k2++) {
                map_part2_to_part1_data[k++] = map_part2_data[ielt2 * (int) s_data + k2];
              }
            }
          }
        }
      }

      else {
        for (int i1 = 0; i1 < ptp->n_part2; i1++) {
          _part2_to_part1_data[i1] = (void *) ((unsigned char *) part2_data[i1] + i * ptp->gnum1_come_from_idx[i1][ptp->n_ref_lnum2[i1]] * (int) s_data);
        }
      }

      for (int i1 = 0; i1 < ptp->n_part1; i1++) {
        ___part1_data[i1] = ((unsigned char *) _part1_data[i1]) + s_data *i *  ptp->part1_to_part2_idx[i1][ptp->n_elt1[i1]];
      }

      if (k_comm == PDM_MPI_COMM_KIND_P2P) {

        PDM_part_to_part_reverse_issend (ptp,
                                         s_data,
                                         1,
                          (const void **)_part2_to_part1_data,
                                         tag,
                                         &(ptp->async_exch_subrequest[_request][2*i]));

        PDM_part_to_part_reverse_irecv (ptp,
                                        s_data,
                                        1,
                                (void **)___part1_data,
                                        tag,
                                        &(ptp->async_exch_subrequest[_request][2*i+1]));

      }

      else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
        abort();

      }

      else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

        printf ("Error PDM_part_to_part_reverse_iexch : "
                 "PDM_STRIDE_CST_INTERLEAVED stride with"
                 " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
        abort();

      }

    }

    if (t_part2_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART2) {
      for (int i = 0; i < ptp->n_part2; i++) {
        free (__part2_to_part1_data[i]);
      }
    }
    free (__part2_to_part1_data);
    __part2_to_part1_data = NULL;
    free (___part1_data);
    ___part1_data = NULL;
  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    void  **_part2_to_part1_data   = (void **) part2_data;
    void  **__part2_to_part1_data   = NULL;

    /*
     *  Create __part2_to_part1_stride and __part2_to_part1_data if necessary
     */

    if (t_part2_data_def == PDM_PART_TO_PART_DATA_DEF_ORDER_PART2) {
      __part2_to_part1_data   = (void **) malloc (sizeof (void*) * ptp->n_part2);

      _part2_to_part1_data   = __part2_to_part1_data;

      for (int i = 0; i < ptp->n_part2; i++) {

        _part2_to_part1_data[i] = malloc (s_data * ptp->gnum1_come_from_idx[i][ptp->n_ref_lnum2[i]] * cst_stride);

        unsigned char *map_part2_to_part1_data = (unsigned char*) _part2_to_part1_data[i];
        unsigned char *map_part2_data = (unsigned char*) part2_data[i];

        int k = 0;
        for (int j = 0; j < ptp->n_ref_lnum2[i]; j++) {
          int ielt2 = ptp->ref_lnum2[i][j] - 1;
          for (int k1 = ptp->gnum1_come_from_idx[i][j]; k1 < ptp->gnum1_come_from_idx[i][j+1]; k1++) {
            for (int k2 = 0; k2 < cst_stride * (int) s_data; k2++) {
              map_part2_to_part1_data[k++] = map_part2_data[ielt2 * cst_stride * (int) s_data + k2];
            }
          }
        }
      }
    }

    *part1_data = malloc(sizeof(void *) * ptp->n_part1);
    void **_part1_data = *part1_data;
    for (int i = 0; i < ptp->n_part1; i++) {
      _part1_data[i] = malloc(s_data * cst_stride * ptp->part1_to_part2_idx[i][ptp->n_elt1[i]]);
    }

    if (k_comm == PDM_MPI_COMM_KIND_P2P) {

      PDM_part_to_part_reverse_issend (ptp,
                                       s_data,
                                       cst_stride,
                        (const void **)_part2_to_part1_data,
                                       tag,
                                       &(ptp->async_exch_subrequest[_request][0]));

      PDM_part_to_part_reverse_irecv (ptp,
                                      s_data,
                                      cst_stride,
                                      _part1_data,
                                      tag,
                                      &(ptp->async_exch_subrequest[_request][1]));

    }

    else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

    else {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " unknown k_comm is not implemented yet\n");
      abort();

    }

    if (__part2_to_part1_data != NULL) {
      for (int i = 0; i < ptp->n_part2; i++) {
        free (__part2_to_part1_data[i]);
      }
      free (__part2_to_part1_data);
      __part2_to_part1_data = NULL;
    }

  }

  else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    if (k_comm == PDM_MPI_COMM_KIND_P2P) {

      _p2p_stride_var_reverse_iexch (ptp,
                                     tag,
                                     t_part2_data_def,
                                     s_data,
                                     part2_stride,
                                     part2_data,
                                     part1_stride,
                                     part1_data,
                                     _request);

    }

    else if (k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_reverse_iexch : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }
  }
}


/**
 *
 * \brief Wait a partial asynchronus exchange
 *
 * \param [in]  ptp      Part to part structure
 * \param [in]  request  Request
 *
 */

void
PDM_part_to_part_reverse_iexch_wait
(
 PDM_part_to_part_t *ptp,
 const int           request
)
{

  if (ptp->async_exch_t_stride[request] == PDM_STRIDE_CST_INTERLEAVED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      for (int i = 0; i < ptp->async_exch_subrequest_s[request]; i++) {

        PDM_part_to_part_reverse_irecv_wait (ptp, ptp->async_exch_subrequest[request][2*i+1]);

        PDM_part_to_part_reverse_issend_wait (ptp, ptp->async_exch_subrequest[request][2*i]);

      }

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLEAVED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

  }

  else if (ptp->async_exch_t_stride[request] == PDM_STRIDE_CST_INTERLACED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      PDM_part_to_part_reverse_irecv_wait (ptp, ptp->async_exch_subrequest[request][0]);

      PDM_part_to_part_reverse_issend_wait (ptp, ptp->async_exch_subrequest[request][1]);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

    else {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_CST_INTERLACED stride with"
               " unknown k_comm is not implemented yet\n");
      abort();

    }

  }

  else if (ptp->async_exch_t_stride[request] == PDM_STRIDE_VAR_INTERLACED) {

    if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_P2P) {

      PDM_part_to_part_reverse_issend_wait(ptp, ptp->async_exch_subrequest[request][0]);

      _p2p_stride_var_reverse_iexch_wait (ptp, request);

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
      abort();

    }

    else if (ptp->async_exch_k_comm[request] == PDM_MPI_COMM_KIND_WIN_RMA) {

      printf ("Error PDM_part_to_part_reverse_iexch_wait : "
               "PDM_STRIDE_VAR_INTERLACED stride with"
               " PDM_MPI_COMM_KIND_WIN_RMA k_comm is not implemented yet\n");
      abort();

    }

  }

  _free_async_exch (ptp, request);
}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_free
(
 PDM_part_to_part_t *ptp
)
{
  if (ptp == NULL) {
    return NULL;
  }

  free(ptp->n_elt1);
  free(ptp->n_elt2);
  for(int i_part = 0; i_part < ptp->n_part1; ++i_part) {
    free(ptp->part1_to_part2_idx[i_part]);
  }
  free(ptp->part1_to_part2_idx);

  if (ptp->gnum1_to_send_buffer != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {
      free (ptp->gnum1_to_send_buffer[i]);
    }
    free (ptp->gnum1_to_send_buffer);
  }

  if (ptp->gnum1_to_send_buffer_idx != NULL) {
    for (int i = 0; i < ptp->n_part1; i++) {
      free (ptp->gnum1_to_send_buffer_idx[i]);
    }
    free (ptp->gnum1_to_send_buffer_idx);
  }

  if (ptp->recv_buffer_to_ref_lnum2 != NULL) {
    for (int i = 0; i < ptp->n_part2; i++) {
      if (ptp->ref_lnum2[i] != NULL) {
        free (ptp->ref_lnum2[i]);
      }
      if (ptp->unref_lnum2[i] != NULL) {
        free (ptp->unref_lnum2[i]);
      }
      free (ptp->gnum1_come_from_idx[i]);
      free (ptp->gnum1_come_from[i]);
      free (ptp->recv_buffer_to_ref_lnum2[i]);
      free (ptp->recv_buffer_to_duplicate_idx[i]);
      free (ptp->recv_buffer_to_duplicate[i]);
    }
    free (ptp->recv_buffer_to_ref_lnum2);
    free (ptp->ref_lnum2);
    free (ptp->unref_lnum2);
    free (ptp->n_ref_lnum2);
    free (ptp->n_unref_lnum2);
    free (ptp->gnum1_come_from_idx);
    free (ptp->gnum1_come_from);
    free (ptp->recv_buffer_to_duplicate_idx);
    free (ptp->recv_buffer_to_duplicate);
  }

  free (ptp->active_rank_send);
  free (ptp->active_rank_recv);

  if (ptp->async_send_l_array != 0) {
    for (int i = 0; i < ptp->async_send_l_array; i++) {
      if (ptp->async_send_buffer[i] != NULL) {
        free (ptp->async_send_buffer[i]);
      }
      if (ptp->async_n_send_buffer[i] != NULL) {
        free (ptp->async_n_send_buffer[i]);
      }
      if (ptp->async_i_send_buffer[i] != NULL) {
        free (ptp->async_i_send_buffer[i]);
      }
      if (ptp->async_send_request[i] != NULL) {
        free (ptp->async_send_request[i]);
      }
    }
    free (ptp->async_send_free);
    free (ptp->async_send_s_data);
    free (ptp->async_send_cst_stride);
    free (ptp->async_send_tag);
    free (ptp->async_send_request);
    free (ptp->async_send_buffer);
    free (ptp->async_n_send_buffer);
    free (ptp->async_i_send_buffer);
  }

  if (ptp->async_recv_l_array != 0) {
    for (int i = 0; i < ptp->async_recv_l_array; i++) {
      if (ptp->async_recv_buffer[i] != NULL) {
        free (ptp->async_recv_buffer[i]);
      }
      if (ptp->async_n_recv_buffer[i] != NULL) {
        free (ptp->async_n_recv_buffer[i]);
      }
      if (ptp->async_i_recv_buffer[i] != NULL) {
        free (ptp->async_i_recv_buffer[i]);
      }
      if (ptp->async_recv_request[i] != NULL) {
        free (ptp->async_recv_request[i]);
      }
    }
    free (ptp->async_recv_free);
    free (ptp->async_recv_s_data);
    free (ptp->async_recv_cst_stride);
    free (ptp->async_recv_tag);
    free (ptp->async_recv_request);
    free (ptp->async_recv_buffer);
    free (ptp->async_n_recv_buffer);
    free (ptp->async_i_recv_buffer);
    free (ptp->async_recv_part2_data);
  }

  free (ptp->default_n_send_buffer);
  free (ptp->default_i_send_buffer);
  free (ptp->default_n_recv_buffer);
  free (ptp->default_i_recv_buffer);

  if (ptp->async_exch_l_array > 0) {
    for (int i = 0; i < ptp->async_exch_l_array; i++) {
      free (ptp->async_exch_subrequest[i]);
      if (ptp->async_exch_recv_n[i] != NULL) {
        free (ptp->async_exch_recv_n[i]);
      }
      if (ptp->async_exch_recv_idx[i] != NULL) {
        free (ptp->async_exch_recv_idx[i]);
      }
    }

    free (ptp->async_exch_free);
    free (ptp->async_exch_subrequest);
    free (ptp->async_exch_subrequest_s);
    free (ptp->async_exch_t_stride);
    free (ptp->async_exch_k_comm);
    free (ptp->async_exch_recv_n);
    free (ptp->async_exch_recv_idx);
    free (ptp->async_exch_part2_stride);

  }

  if (ptp->async_alltoall_l_array > 0) {
    free (ptp->async_alltoall_free);
    free (ptp->async_alltoall_subrequest);
  }

  ptp->async_exch_n_free  = 0;
  ptp->async_exch_l_array = 0;

  ptp->async_send_n_free   = 0;
  ptp->async_send_l_array  = 0;

  ptp->async_recv_n_free   = 0;
  ptp->async_recv_l_array  = 0;

  ptp->async_alltoall_n_free  = 0;
  ptp->async_alltoall_l_array = 0;

  free(ptp);
  return NULL;

}


/**
 *
 * \brief Get number of partitions
 *
 * \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
 * \param [out] n_part1   Number of partitions on side 1
 * \param [out] n_part2   Number of partitions on side 2
 *
 */

void
PDM_part_to_part_n_part_get
(
 PDM_part_to_part_t *ptp,
 int                *n_part1,
 int                *n_part2
 )
{
  assert(ptp != NULL);

  *n_part1 = ptp->n_part1;
  *n_part2 = ptp->n_part2;
}


/**
 *
 * \brief Get number of partitions and n_elt1 and n_elt2
 *
 * \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
 * \param [out] n_part1   Number of partitions on side 1
 * \param [out] n_part2   Number of partitions on side 2
 * \param [out] n_elt1    Number of gnum1 element
 * \param [out] n_elt2    Number of gnum2 element
 *
 */
void
PDM_part_to_part_n_part_and_n_elt_get
(
 PDM_part_to_part_t *ptp,
 int                *n_part1,
 int                *n_part2,
 int               **n_elt1,
 int               **n_elt2
 )
{
  assert(ptp != NULL);
  *n_part1 = ptp->n_part1;
  *n_part2 = ptp->n_part2;
  *n_elt1  = ptp->n_elt1;
  *n_elt2  = ptp->n_elt2;
}


/**
 *
 * \brief Get referenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   i_part        Id of partition
 * \param [out]  n_ref_lnum2   Number of referenced gnum2
 * \param [out]  ref_lnum2     Referenced gnum2
 *
 */

void
PDM_part_to_part_ref_lnum2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_ref_lnum2,
       int                **ref_lnum2
)
{
  assert(ptp != NULL);
  assert(i_part < ptp->n_part2);

  *n_ref_lnum2 = ptp->n_ref_lnum2[i_part];
  *ref_lnum2   = ptp->ref_lnum2[i_part];
}


/**
 *
 * \brief Get unreferenced gnum2 elements
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   i_part        Id of partition
 * \param [out]  n_unref_lnum2 Number of unreferenced gnum2
 * \param [out]  unref_lnum2   Unreferenced gnum2
 *
 */

void
PDM_part_to_part_unref_lnum2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_unref_lnum2,
       int                **unref_lnum2
)
{
  assert(ptp != NULL);
  assert(i_part < ptp->n_part2);

  *n_unref_lnum2 = ptp->n_unref_lnum2[i_part];
  *unref_lnum2   = ptp->unref_lnum2[i_part];
}


/**
 *
 * \brief Get gnum come from gnum1 for each referenced gnum2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   i_part        Id of partition
 * \param [out]  gnum1_come_from_idx Index for gnum1_come_from array
 * \param [out]  gnum1_come_from     Gnum come from gnum1 for each referenced gnum2
 *
 */

void
PDM_part_to_part_gnum1_come_from_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                **gnum1_come_from_idx,
       PDM_g_num_t        **gnum1_come_from
)
{
  assert(ptp != NULL);
  assert(i_part < ptp->n_part2);

  *gnum1_come_from_idx = ptp->gnum1_come_from_idx[i_part];
  *gnum1_come_from     = ptp->gnum1_come_from[i_part];
}


/**
 *
 * \brief Get selected numbers of part2 (only index)
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   i_part              Id of partition
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 *
 */

void
PDM_part_to_part_part1_to_part2_idx_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_elt1,
       int                **part1_to_part2_idx
)
{
  assert(ptp != NULL);
  assert(i_part < ptp->n_part1);

  *n_elt1             = (int          ) ptp->n_elt1[i_part];
  *part1_to_part2_idx = (int         *) ptp->part1_to_part2_idx[i_part];
}


/**
 *
 * \brief Get selected numbers of part2
 *
 * \param [in]   ptp                 Block to part structure
 * \param [in]   i_part              Id of partition
 * \param [out]  n_elt1              Number of gnum1 element
 * \param [out]  part1_to_part2_idx  Index of data to send to gnum2 from gnum1
 *                                  (for each part size : \ref n_elt1+1)
 * \param [out]  part1_to_part2      Data to send to gnum2 from gnum1 for each part
 *
 */

void
PDM_part_to_part_part1_to_part2_single_part_get
(
       PDM_part_to_part_t  *ptp,
 const int                  i_part,
       int                 *n_elt1,
       int                **part1_to_part2_idx,
       PDM_g_num_t        **part1_to_part2
)
{
  assert(ptp != NULL);
  assert(i_part < ptp->n_part1);

  *n_elt1             = (int          ) ptp->n_elt1[i_part];
  *part1_to_part2_idx = (int         *) ptp->part1_to_part2_idx[i_part];
  *part1_to_part2     = (PDM_g_num_t *) ptp->part1_to_part2[i_part];
}

/**
 *
 * \brief Get indirection from part1_to_part2 to buffer send (usefull to setup buffer outside ptp )
 *
 * \param [in]   ptp                       Block to part structure
 * \param [out]  gnum1_to_send_buffer_idx  Index of data to send to gnum2 from gnum1
 *                                           (for each part size : \ref n_elt1+1)
 * \param [out]  gnum1_to_send_buffer      For each gnum1 the position in send buffer
 *
 */
void
PDM_part_to_part_gnum1_to_send_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                 ***gnum1_to_send_buffer_idx,
 int                 ***gnum1_to_send_buffer
)
{
  *gnum1_to_send_buffer_idx = ptp->gnum1_to_send_buffer_idx;
  *gnum1_to_send_buffer     = ptp->gnum1_to_send_buffer;
}

/**
 *
 * \brief Get indirection from ref_lnum2 to buffer recv (usefull to setup buffer outside ptp )
 *
 * \param [in]   ptp                       Block to part structure
 * \param [out]  recv_buffer_to_ref_lnum2  For each gnum2 the position in recv buffer ( size = gnum1_come_from_idx[n_ref_lnum2])
 *
 */
void
PDM_part_to_part_recv_buffer_to_ref_lnum2_get
(
 PDM_part_to_part_t    *ptp,
 int                 ***recv_buffer_to_ref_lnum2
)
{
  *recv_buffer_to_ref_lnum2 = ptp->recv_buffer_to_ref_lnum2;
}


/**
 *
 * \brief Get buffer size and stride for send
 *
 * \param [in]   ptp                       Block to part structure
 * \param [out]  default_n_send_buffer     Number of entities to send (size = n_rank)
 * \param [out]  default_i_send_buffer     Index (size = n_rank + 1)
 *
 */
void
PDM_part_to_part_default_send_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                  **default_n_send_buffer,
 int                  **default_i_send_buffer
)
{
  *default_n_send_buffer = ptp->default_n_send_buffer;
  *default_i_send_buffer = ptp->default_i_send_buffer;
}


/**
 *
 * \brief Get buffer size and stride for recv
 *
 * \param [in]   ptp                       Block to part structure
 * \param [out]  default_n_recv_buffer     Number of entities to recv (size = n_rank)
 * \param [out]  default_i_recv_buffer     Index (size = n_rank + 1)
 *
 */
void
PDM_part_to_part_default_recv_buffer_get
(
 PDM_part_to_part_t    *ptp,
 int                  **default_n_recv_buffer,
 int                  **default_i_recv_buffer
)
{
  *default_n_recv_buffer = ptp->default_n_recv_buffer;
  *default_i_recv_buffer = ptp->default_i_recv_buffer;
}


/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptp          Part to part structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_part_n_ranks_get
(
 PDM_part_to_part_t    *ptp
)
{
  return ptp->n_rank;
}

#ifdef __cplusplus
}
#endif

