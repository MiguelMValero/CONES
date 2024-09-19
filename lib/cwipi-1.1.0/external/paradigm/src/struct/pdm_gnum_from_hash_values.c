/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_gnum_from_hash_values.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER_HASH_VALUES 9

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/**
 * \enum _gnum_from_hash_values_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  EQUILIBRATE_DISTIBUTION       = 1,
  FIRST_EXCHANGE_PREPARE        = 2,
  FIRST_EXCHANGE                = 3,
  SECOND_EXCHANGE_PREPARE       = 4,
  SECOND_EXCHANGE               = 5,
  BLOCK_SORT                    = 6,
  BLOCK_EQUAL                   = 7,
  REVERSE_EXCHANGE              = 8,
  END                           = 9,

} _gnum_from_hash_values_step_t;


/**
 * \struct _pdm_gnum_from_hv_t
 * \brief  Define a global numberring
 *
 */

typedef struct {
  int             n_part;          /*!< Number of partitions                     */
  PDM_bool_t      equilibrate;     /*!< Equilibrate the hash values distribution */
  PDM_MPI_Comm    comm;            /*!< MPI communicator                         */
  int             n_rank;          /*!< MPI communicator size                    */

  int            *n_elts;          /*!< Number of elements in partitions         */
  size_t        **part_hkeys;
  unsigned char **part_hdata;
  int           **part_hstri;
  size_t          s_data;

  gnum_from_hv_compare fcompare;
  gnum_from_hv_equal   fequal;

  PDM_g_num_t     n_g_elt;        /*!< Global number of elements                 */
  PDM_g_num_t   **g_nums;         /*!< Global numbering of elements              */

  size_t        *distribution;

  PDM_timer_t *timer;                        /*!< Timer */

  double times_elapsed[NTIMER_HASH_VALUES];  /*!< Elapsed time */

  double times_cpu[NTIMER_HASH_VALUES];      /*!< CPU time */

  double times_cpu_u[NTIMER_HASH_VALUES];    /*!< User CPU time */

  double times_cpu_s[NTIMER_HASH_VALUES];    /*!< System CPU time */

  PDM_ownership_t owner;       /*!< Ownership */
  int  tag_results_get ;       /*!< Tag call to PDM_dist_cellcenter_surf_get function */ 

} _pdm_gnum_from_hv_t;

/*============================================================================
 * Global variable
 *============================================================================*/

// static PDM_Handles_t *_gnums_from_hv   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

// static _pdm_gnum_from_hv_t *
// _get_from_id
// (
//  int  id
// )
// {

//   _pdm_gnum_from_hv_t *gnum_from_hv = (_pdm_gnum_from_hv_t *) PDM_Handles_get (_gnums_from_hv, id);

//   if (gnum_from_hv == NULL) {
//     PDM_error(__FILE__, __LINE__, 0, "PDM_gnum_from_hv error : Bad identifier\n");
//   }

//   return gnum_from_hv;
// }

/**
 *
 * \brief Compute with equilibrate algorithm
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_compute_distribution_equilibrate
(
 _pdm_gnum_from_hv_t *_gnum
)
{
  PDM_UNUSED(_gnum);
  printf("_gnum_from_hv_compute_equilibrate Not implemented \n");
  abort();

}

/**
 *
 * \brief Setup a naive distribution from min max of data
 */
static void
setup_distribution_from_min_max
(
 size_t       min_elt,
 size_t       max_elt,
 size_t*      distribution,
 int          n_dist
)
{
  size_t nelmt = max_elt - min_elt + 1;

  assert(nelmt > 0);

  size_t quotient  = nelmt/n_dist;
  int    remainder = nelmt%n_dist;

  // printf(PDM_FMT_G_NUM"\n", quotient);
  // printf(PDM_FMT_G_NUM"\n", remainder);

  distribution[0] = min_elt;
  for(int i = 1; i < n_dist+1; ++i) {
    distribution[i] = quotient;
    PDM_g_num_t i1 = i - 1;
    if(i1 < remainder){
      distribution[i] += 1;
    }
  }

  for(int i = 0; i < n_dist; ++i) {
    distribution[i+1] += distribution[i];
  }

  // PDM_log_trace_array_size_t(distribution, n_dist+1, "distribution:: ");
}

/**
 *
 * \brief Compute with equilibrate algorithm
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_compute_distribution
(
 _pdm_gnum_from_hv_t *_gnum_from_hv
)
{

  size_t max_key_loc = 0;
  size_t min_key_loc = SIZE_MAX;

  size_t max_key = 0;
  size_t min_key = SIZE_MAX;

  for(int i_part = 0; i_part < _gnum_from_hv->n_part; i_part++){
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){
      max_key_loc = PDM_MAX(max_key_loc, _gnum_from_hv->part_hkeys[i_part][ielt]);
      min_key_loc = PDM_MIN(min_key_loc, _gnum_from_hv->part_hkeys[i_part][ielt]);

    }
  }
  int ierr;
  ierr = PDM_MPI_Allreduce(&min_key_loc, &min_key, 1, PDM_MPI_UNSIGNED_LONG, PDM_MPI_MIN, _gnum_from_hv->comm);
  assert(ierr == 0);

  ierr = PDM_MPI_Allreduce(&max_key_loc, &max_key, 1, PDM_MPI_UNSIGNED_LONG, PDM_MPI_MAX, _gnum_from_hv->comm);
  assert(ierr == 0);

  // printf(" max_key:: %lu \n", max_key);
  // printf(" min_key:: %lu \n", min_key);

  /* Prepare distribution from min and max elements */
  setup_distribution_from_min_max(min_key, max_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank);

}

/**
 *
 * \brief Compute but without equilibrate
 *
 * \param [in]   _gnum          Current _pdm_gnum_from_hv_t structure
 *
 */
static void
_gnum_from_hv_compute
(
 _pdm_gnum_from_hv_t *_gnum_from_hv
)
{
  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  // printf("_gnum_from_hv_compute \n");

  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Equilibrate the block distribution
   */
  if(_gnum_from_hv->equilibrate) {
    _compute_distribution_equilibrate(_gnum_from_hv);
  } else {
    _compute_distribution(_gnum_from_hv);
  }

  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[EQUILIBRATE_DISTIBUTION] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[EQUILIBRATE_DISTIBUTION]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[EQUILIBRATE_DISTIBUTION]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[EQUILIBRATE_DISTIBUTION]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Remapping of partition data in block data according to the hash values distribution
   */
  int* n_key_send  = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_key_recv  = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_data_send = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));
  int* n_data_recv = (int *) malloc( _gnum_from_hv->n_rank * sizeof(int));

  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    n_key_send[i] = 0;
    n_data_send[i] = 0;
  }

  /*
   * Prepare send
   */
  int** elmt_to_proc = (int **) malloc(_gnum_from_hv->n_part * sizeof(int*));
  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){
    elmt_to_proc[i_part] = (int *) malloc( _gnum_from_hv->n_elts[i_part] * sizeof(int));
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){

      size_t g_key = _gnum_from_hv->part_hkeys[i_part][ielt];
      // log_trace(" Search for :: %lu\n", g_key);
      int t_rank = PDM_binary_search_gap_size_t(g_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank+1);
      elmt_to_proc[i_part][ielt] = t_rank;

      // log_trace(" Found in t_rank :: %d\n", t_rank);
      // n_data_send[t_rank] += _gnum_from_hv->s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      n_data_send[t_rank] += _gnum_from_hv->part_hstri[i_part][ielt];
      n_key_send[t_rank]++;

    }
  }

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[FIRST_EXCHANGE_PREPARE] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[FIRST_EXCHANGE_PREPARE]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[FIRST_EXCHANGE_PREPARE]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[FIRST_EXCHANGE_PREPARE]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Exchange
   */
  PDM_MPI_Alltoall(n_key_send , 1, PDM_MPI_INT, n_key_recv , 1, PDM_MPI_INT, _gnum_from_hv->comm);
  PDM_MPI_Alltoall(n_data_send, 1, PDM_MPI_INT, n_data_recv, 1, PDM_MPI_INT, _gnum_from_hv->comm);

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[FIRST_EXCHANGE] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[FIRST_EXCHANGE]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[FIRST_EXCHANGE]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[FIRST_EXCHANGE]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);


  /*
   * Prepare utility array to setup the second exchange
   */
  int* i_key_send  = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_key_recv  = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_data_send = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));
  int* i_data_recv = (int *) malloc( (_gnum_from_hv->n_rank+1) * sizeof(int));

  i_key_send[0] = 0;
  i_key_recv[0] = 0;
  i_data_send[0] = 0;
  i_data_recv[0] = 0;
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    i_key_send[i+1] = i_key_send[i] + n_key_send[i];
    i_key_recv[i+1] = i_key_recv[i] + n_key_recv[i];
    i_data_send[i+1] = i_data_send[i] + n_data_send[i];
    i_data_recv[i+1] = i_data_recv[i] + n_data_recv[i];
    n_key_send[i]  = 0;
    n_data_send[i] = 0;
  }

  int s_send_keys = i_key_send[_gnum_from_hv->n_rank];
  int s_recv_keys = i_key_recv[_gnum_from_hv->n_rank];
  int s_send_data = i_data_send[_gnum_from_hv->n_rank] * _gnum_from_hv->s_data;
  int s_recv_data = i_data_recv[_gnum_from_hv->n_rank] * _gnum_from_hv->s_data;

  if(0 == 1) {
    log_trace("s_send_keys::%d - %d \n", s_send_keys, s_send_keys/_gnum_from_hv->s_data);
    log_trace("s_recv_keys::%d - %d \n", s_recv_keys, s_recv_keys/_gnum_from_hv->s_data);
    log_trace("s_send_data::%d\n", s_send_data);
    log_trace("s_recv_data::%d\n", s_recv_data);
    log_trace("i_data_send[_gnum_from_hv->n_rank]::%d\n", i_data_send[_gnum_from_hv->n_rank]);
    log_trace("i_data_recv[_gnum_from_hv->n_rank]::%d\n", i_data_recv[_gnum_from_hv->n_rank]);
  }

  /*
   * Allocate
   */
  size_t *send_buffer_keys = (size_t *) malloc(sizeof(size_t) * s_send_keys );
  size_t *recv_buffer_keys = (size_t *) malloc(sizeof(size_t) * s_recv_keys );

  int *send_buffer_stri = (int *) malloc(sizeof(int) * ( s_send_keys + 1) );
  int *recv_buffer_stri = (int *) malloc(sizeof(int) * ( s_recv_keys + 1) );

  unsigned char *send_buffer_data = (unsigned char *) malloc(sizeof(unsigned char) * s_send_data);
  unsigned char *recv_buffer_data = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_data);

  /*
   * Setup send_buffer
   */
  int s_data = (int) _gnum_from_hv->s_data;
  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){

    unsigned char* _part_data = (unsigned char *) _gnum_from_hv->part_hdata[i_part];

    int idx = 0;
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){

      size_t g_key = _gnum_from_hv->part_hkeys[i_part][ielt];
      // log_trace(" Search for :: %lu\n", g_key);
      // int t_rank = PDM_binary_search_gap_size_t(g_key, _gnum_from_hv->distribution, _gnum_from_hv->n_rank+1);
      int t_rank = elmt_to_proc[i_part][ielt];

      // log_trace(" Found in t_rank :: %d\n", t_rank);

      /* Send key and stride */
      int idx_send = i_key_send[t_rank]+n_key_send[t_rank]++;
      send_buffer_keys[idx_send] = g_key;
      send_buffer_stri[idx_send] = _gnum_from_hv->part_hstri[i_part][ielt];

      /* Send data */
      int n_data = s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      int shift  = n_data_send[t_rank];
      // log_trace(" n_data = %d | shift = %d | idx = %d \n", n_data, shift, idx);
      for(int i_data = 0; i_data < n_data; ++i_data) {
        int shift_tot = ( i_data_send[t_rank] + shift)*s_data;
        send_buffer_data[shift_tot+i_data] = _part_data[idx*s_data+i_data];
      }
      idx += _gnum_from_hv->part_hstri[i_part][ielt];

      // n_data_send[t_rank] += s_data * _gnum_from_hv->part_hstri[i_part][ielt];
      n_data_send[t_rank] += _gnum_from_hv->part_hstri[i_part][ielt];

    }
  }
  // abort();

  // int* send_buffer_data_int = (int*) send_buffer_data;
  // PDM_log_trace_array_int(send_buffer_data_int, s_send_data/s_data, "send_buffer_data_int:: ");

  if(0 == 1){
    log_trace("s_data %d %d \n ", s_data, s_send_data/s_data);
    log_trace("send_buffer_data_char:: ");
    for(int i = 0; i < s_send_data/s_data; ++i){
      log_trace("%lu ", send_buffer_data[i]);
    }
    log_trace("\n");
  }

  for(int i = 0; i < _gnum_from_hv->n_rank+1; ++i){
    i_data_send[i] = i_data_send[i] * s_data;
    i_data_recv[i] = i_data_recv[i] * s_data;
  }
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    n_data_send[i] = n_data_send[i] * s_data;
    n_data_recv[i] = n_data_recv[i] * s_data;
  }


  if(0 == 1){
    PDM_log_trace_array_int(i_key_send, _gnum_from_hv->n_rank+1, "i_key_send:: ");
    PDM_log_trace_array_int(n_key_send, _gnum_from_hv->n_rank  , "n_key_send:: ");
    PDM_log_trace_array_int(i_key_recv, _gnum_from_hv->n_rank+1, "i_key_recv:: ");
    PDM_log_trace_array_int(n_key_recv, _gnum_from_hv->n_rank  , "n_key_recv:: ");

    PDM_log_trace_array_int(i_data_send, _gnum_from_hv->n_rank+1, "i_data_send:: ");
    PDM_log_trace_array_int(n_data_send, _gnum_from_hv->n_rank  , "n_data_send:: ");
    PDM_log_trace_array_int(i_data_recv, _gnum_from_hv->n_rank+1, "i_data_recv:: ");
    PDM_log_trace_array_int(n_data_recv, _gnum_from_hv->n_rank  , "n_data_recv:: ");
  }

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[SECOND_EXCHANGE_PREPARE] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[SECOND_EXCHANGE_PREPARE]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[SECOND_EXCHANGE_PREPARE]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[SECOND_EXCHANGE_PREPARE]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Exchange
   */
  PDM_MPI_Alltoallv(send_buffer_keys, n_key_send, i_key_send, PDM_MPI_UNSIGNED_LONG,
                    recv_buffer_keys, n_key_recv, i_key_recv, PDM_MPI_UNSIGNED_LONG, _gnum_from_hv->comm);

  PDM_MPI_Alltoallv(send_buffer_stri, n_key_send, i_key_send, PDM_MPI_INT,
                    recv_buffer_stri, n_key_recv, i_key_recv, PDM_MPI_INT, _gnum_from_hv->comm);

  PDM_MPI_Alltoallv(send_buffer_data, n_data_send, i_data_send, PDM_MPI_BYTE,
                    recv_buffer_data, n_data_recv, i_data_recv, PDM_MPI_BYTE, _gnum_from_hv->comm);

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[SECOND_EXCHANGE] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[SECOND_EXCHANGE]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[SECOND_EXCHANGE]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[SECOND_EXCHANGE]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Verbose
   */
  if(0 == 1){
    PDM_log_trace_array_size_t(send_buffer_keys, s_send_keys, "send_buffer_keys:: ");
    PDM_log_trace_array_size_t(recv_buffer_keys, s_recv_keys, "recv_buffer_keys:: ");
    PDM_log_trace_array_int(send_buffer_stri, s_send_keys, "send_buffer_stri:: ");
    PDM_log_trace_array_int(recv_buffer_stri, s_recv_keys, "recv_buffer_stri:: ");
  }

  // int* recv_buffer_data_int = (int*) recv_buffer_data;
  // PDM_log_trace_array_int(recv_buffer_data_int, s_recv_data/s_data, "recv_buffer_data_int:: ");
  // log_trace("s_data %d %d \n ", s_data, s_recv_data/s_data);
  // log_trace("recv_buffer_data_char:: ");
  // for(int i = 0; i < s_recv_data/s_data; ++i){
  //   log_trace("%lu ", recv_buffer_data[i]);
  // }
  // log_trace("\n");

  /*
   * Rebuild a total stride
   */
  int tmp1 = recv_buffer_stri[0];
  recv_buffer_stri[0] = 0;
  for(int i = 0; i < s_recv_keys; ++i){
    int tmp2 = recv_buffer_stri[i+1];
    recv_buffer_stri[i+1] = recv_buffer_stri[i] + tmp1;
    tmp1 = tmp2;
  }

  // if(0 == 0){
  //   // PDM_log_trace_array_int(recv_buffer_stri, s_recv_keys+1, "recv_buffer_stri:: ");
  //   char* t = (char*) malloc( sizeof(char) * (s_recv_data + 1));
  //   // char* t = (char* )recv_buffer_data;
  //   for(int k = 0; k < s_recv_data; ++k){
  //     t[k] = (char) recv_buffer_data[k];
  //   }
  //   t[s_recv_data] = '\0';
  //   log_trace(" -------------- \n ");
  //   log_trace("%s \n ", t);
  //   log_trace(" -------------- \n ");
  //   free(t);
  // }

  /*
   * Generate global numbering from the block_data
   */
  int* order = (int *) malloc( sizeof(int) * s_recv_keys);
  for(int i = 0; i < s_recv_keys; ++i){
    order[i] = i;
  }

  // _gnum_from_hv->fcompare = PDM_operator_compare_connectivity;
  // _gnum_from_hv->fequal   = PDM_operator_equal_connectivity;

  PDM_user_defined_sort* us = (PDM_user_defined_sort *) malloc( sizeof(PDM_user_defined_sort) );
  us->idx = recv_buffer_stri;
  us->arr = recv_buffer_data;
  us->key = recv_buffer_keys;

  // A faire avec l'autre structure pour trier les clés aussi ...

  PDM_sort_int_special(order, s_recv_keys, _gnum_from_hv->fcompare, (void*) us);

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[BLOCK_SORT] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[BLOCK_SORT]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[BLOCK_SORT]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[BLOCK_SORT]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);

  /*
   * Panic verbose
   */
  if(0 == 1){
    PDM_log_trace_array_int(order, s_recv_keys, "order:: ");
    log_trace("order = ");
    int* arr_tmp = (int*) us->arr;
    for(int i = 0; i < s_recv_keys; ++i){
      log_trace("%d --> ", (int)order[i]);
      int j   = order[i];
      for(int k = us->idx[j]; k < us->idx[j+1]; ++k ){
        log_trace(" %d ", arr_tmp[k]);
      }
      log_trace("\n");
    }
    log_trace("\n");
  }

  /*
   * Use operator == to have an global numbering
   */
  PDM_g_num_t* blk_ln_to_gn = (PDM_g_num_t*) malloc( sizeof(PDM_g_num_t*) * s_recv_keys);
  PDM_g_num_t next_id = 0;
  PDM_g_num_t n_id    = 0;
  PDM_g_num_t last_id = -1;
  for(int i = 0; i < s_recv_keys; ++i){
    // log_trace(" generate g_id :: %d \n", i);
    if(i != 0 && _gnum_from_hv->fequal(&order[i], &last_id, (void*) us)){
      // log_trace(" \t Cas 1 :: order[%d] = %d | next_id : %d\n", i, order[i], next_id);
      blk_ln_to_gn[order[i]] = next_id;
    } else {
      next_id++;
      n_id++;
      // log_trace(" \t Cas 2 :: order[%d] = %d | next_id : %d\n", i, order[i], next_id);
      blk_ln_to_gn[order[i]] = next_id;
      last_id = order[i];
    }
  }

  // log_trace("n_id   :: %d \n", n_id);
  // log_trace("next_id:: %d \n", next_id);

  /*
   *  Rebuild the global numbering
   */
  PDM_g_num_t shift_g;
  PDM_MPI_Scan(&n_id, &shift_g, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _gnum_from_hv->comm);
  shift_g -= n_id;

  PDM_g_num_t _max_loc = -1;
  for(int i = 0; i < s_recv_keys; ++i){
    blk_ln_to_gn[i] += shift_g;
    _max_loc = PDM_MAX (_max_loc, blk_ln_to_gn[i]);
  }

  /*
   *  Reduce the total number of elements
   */
  PDM_MPI_Allreduce (&_max_loc,
                     &_gnum_from_hv->n_g_elt,
                     1,
                     PDM__PDM_MPI_G_NUM,
                     PDM_MPI_MAX,
                     _gnum_from_hv->comm);

  // log_trace("_gnum_from_hv->n_g_elt:: %d \n", _gnum_from_hv->n_g_elt);

  /*
   * Panic verbose
   */
  if(0 == 1){
    PDM_log_trace_array_long(blk_ln_to_gn, s_recv_keys, "blk_ln_to_gn:: ");
  }

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[BLOCK_EQUAL] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[BLOCK_EQUAL]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[BLOCK_EQUAL]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[BLOCK_EQUAL]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);


  /*
   * Reverse all_to_all exchange in order to remap global id on current partition
   */
  PDM_g_num_t* part_ln_to_gn = (PDM_g_num_t*) malloc( sizeof(PDM_g_num_t*) * s_send_keys);

  PDM_MPI_Alltoallv(blk_ln_to_gn , n_key_recv, i_key_recv, PDM__PDM_MPI_G_NUM,
                    part_ln_to_gn, n_key_send, i_key_send, PDM__PDM_MPI_G_NUM, _gnum_from_hv->comm);

  /*
   *  Remise en place dans chaque partition
   */
  for(int i = 0; i < _gnum_from_hv->n_rank; ++i){
    n_key_send[i]  = 0;
  }

  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){
    for(int ielt = 0; ielt < _gnum_from_hv->n_elts[i_part]; ++ielt){
      int t_rank   = elmt_to_proc[i_part][ielt];
      int idx_send = i_key_send[t_rank] + n_key_send[t_rank]++;
      _gnum_from_hv->g_nums[i_part][ielt] = part_ln_to_gn[idx_send];
    }
  }

  free(n_key_send);
  free(n_key_recv);
  free(n_data_send);
  free(n_data_recv);
  free(i_key_send);
  free(i_key_recv);
  free(i_data_send);
  free(i_data_recv);
  free(send_buffer_keys);
  free(recv_buffer_keys);
  free(send_buffer_data);
  free(recv_buffer_data);
  free(send_buffer_stri);
  free(recv_buffer_stri);
  free(blk_ln_to_gn);
  free(part_ln_to_gn);
  for(int i_part = 0; i_part < _gnum_from_hv->n_part; ++i_part){
    free(elmt_to_proc[i_part]);
  }
  free(elmt_to_proc);
  free(order);
  free(us);

  /*
   * Timer
   */
  PDM_timer_hang_on(_gnum_from_hv->timer);
  e_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  e_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);

  _gnum_from_hv->times_elapsed[REVERSE_EXCHANGE] += e_t_elapsed - b_t_elapsed;
  _gnum_from_hv->times_cpu[REVERSE_EXCHANGE]     += e_t_cpu     - b_t_cpu;
  _gnum_from_hv->times_cpu_u[REVERSE_EXCHANGE]   += e_t_cpu_u   - b_t_cpu_u;
  _gnum_from_hv->times_cpu_s[REVERSE_EXCHANGE]   += e_t_cpu_s   - b_t_cpu_s;

  /* Reset for the next step */
  b_t_elapsed = PDM_timer_elapsed(_gnum_from_hv->timer);
  b_t_cpu     = PDM_timer_cpu(_gnum_from_hv->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(_gnum_from_hv->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(_gnum_from_hv->timer);
  PDM_timer_resume(_gnum_from_hv->timer);


}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   equilibrate  Use algorithm to equilibrate the block treatment (hash value is not a priori equi-reparti)
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 * \param [in]   owner        Owner
 *
 * \return     Pointer to \ref PDM_gnum_from_hv_t object
 */

PDM_gnum_from_hv_t *
PDM_gnum_from_hash_values_create
(
 const int            n_part,
 const PDM_bool_t     equilibrate,
 const size_t         s_data,
 gnum_from_hv_compare fcompare,
 gnum_from_hv_equal   fequal,
 const PDM_MPI_Comm   comm,
 const PDM_ownership_t owner
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);


  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) malloc(sizeof(_pdm_gnum_from_hv_t));

  _gnum_from_hv->fcompare    = fcompare;
  _gnum_from_hv->fequal      = fequal;

  _gnum_from_hv->n_part      = n_part;
  _gnum_from_hv->equilibrate = equilibrate;
  _gnum_from_hv->comm        = comm;
  _gnum_from_hv->n_rank      = n_rank;
  _gnum_from_hv->n_g_elt     = -1;
  _gnum_from_hv->g_nums      = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_part);

  _gnum_from_hv->s_data      = s_data;
  _gnum_from_hv->n_elts      = (int            *) malloc (sizeof(int            ) * n_part);
  _gnum_from_hv->part_hkeys  = (size_t        **) malloc (sizeof(size_t        *) * n_part);
  _gnum_from_hv->part_hstri  = (int           **) malloc (sizeof(int           *) * n_part);
  _gnum_from_hv->part_hdata  = (unsigned char **) malloc (sizeof(unsigned char *) * n_part);

  for (int i = 0; i < n_part; i++) {
    _gnum_from_hv->g_nums[i]     = NULL;
    _gnum_from_hv->part_hkeys[i] = NULL;
    _gnum_from_hv->part_hstri[i] = NULL;
    _gnum_from_hv->part_hdata[i] = NULL;
  }

  _gnum_from_hv->distribution = (size_t *) malloc (sizeof(size_t ) * ( n_rank + 1 ));

  _gnum_from_hv->timer = PDM_timer_create();

  for (int i = 0; i < NTIMER_HASH_VALUES; i++) {
    _gnum_from_hv->times_elapsed[i] = 0.;
    _gnum_from_hv->times_cpu[i] = 0.;
    _gnum_from_hv->times_cpu_u[i] = 0.;
    _gnum_from_hv->times_cpu_s[i] = 0.;
  }

  _gnum_from_hv->owner = owner;
  _gnum_from_hv->tag_results_get = 0;

  return (PDM_gnum_from_hv_t *) _gnum_from_hv;
}


/**
 *
 * \brief Set hash values for one partition
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   part_hkey    For each elements the hash value associated
 * \param [in]   part_strid   Stride between each data in part_hdata
 * \param [in]   part_hdata   Partition data which compute the hash value, we need it to setup in a block way
 *
 */

void
PDM_gnum_set_hash_values
(
 PDM_gnum_from_hv_t  *gnum_from_hv,
 const int            i_part,
 const int            n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) gnum_from_hv;

  assert(_gnum_from_hv->part_hkeys != NULL);
  assert(_gnum_from_hv->part_hstri != NULL);
  assert(_gnum_from_hv->part_hdata != NULL);

  assert(_gnum_from_hv->part_hkeys[i_part] == NULL);
  assert(_gnum_from_hv->part_hstri[i_part] == NULL);
  assert(_gnum_from_hv->part_hdata[i_part] == NULL);

  _gnum_from_hv->n_elts[i_part]      = n_elts;
  _gnum_from_hv->part_hkeys[i_part]  = (size_t        *) part_hkeys;
  _gnum_from_hv->part_hstri[i_part]  = (int           *) part_hstri;
  _gnum_from_hv->part_hdata[i_part]  = (unsigned char *) part_hdata;

  // log_trace("n_elts:: %d \n", n_elts);
  // int idx = 0;
  // for(int i = 0; i < n_elts; ++i){

  //   log_trace("part_hstri:: %d ",  _gnum_from_hv->part_hstri[i_part][i]);
  //   log_trace("part_hdata:: ");

  //   for(int i_data = 0; i_data < _gnum_from_hv->part_hstri[i_part][i]; ++i_data) {
  //     log_trace("%lu ", _gnum_from_hv->part_hdata[i_part][idx++]);
  //   }
  // }
  // log_trace("\n");
  // char* t = (char* )part_hdata;
  // log_trace(" -------------- \n ");
  // log_trace("%s \n ", t);
  // log_trace(" -------------- \n ");


  _gnum_from_hv->g_nums[i_part] = (PDM_g_num_t * ) malloc( n_elts * sizeof(PDM_g_num_t));

}


/**
 *
 * \brief Compute
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 *
 */

void
PDM_gnum_from_hv_compute
(
 PDM_gnum_from_hv_t *gnum_from_hv
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) gnum_from_hv;

  _gnum_from_hv_compute(_gnum_from_hv);

}


/**
 *
 * \brief Get the global ids for the current partition
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   i_part       Current partition
 *
 * \return  Array of global ids for the current partition
 *
 */

PDM_g_num_t *
PDM_gnum_from_hv_get
(
 PDM_gnum_from_hv_t *gnum_from_hv,
 const int           i_part
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) gnum_from_hv;

  _gnum_from_hv->tag_results_get = 1;
  return _gnum_from_hv->g_nums[i_part];
}


/**
 *
 * \brief Dump elapsed an CPU time
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 *
 */

void
PDM_gnum_from_hv_dump_times
(
 PDM_gnum_from_hv_t *gnum_from_hv
)
{

  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) gnum_from_hv;

  double t_elaps_max[NTIMER_HASH_VALUES];
  PDM_MPI_Allreduce (_gnum_from_hv->times_elapsed, t_elaps_max, NTIMER_HASH_VALUES,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, _gnum_from_hv->comm);

  double t_cpu_max[NTIMER_HASH_VALUES];
  PDM_MPI_Allreduce (_gnum_from_hv->times_cpu, t_cpu_max, NTIMER_HASH_VALUES,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, _gnum_from_hv->comm);

  int i_rank;
  PDM_MPI_Comm_rank (_gnum_from_hv->comm, &i_rank);

  if (i_rank == 0) {

    // PDM_printf( "hash_values timer : all (elapsed and cpu) : %12.5es %12.5es\n",
    //             t1max, t2max);
    PDM_printf( "PDM_gnum_from_hv timer : Equilibrate distribution (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[EQUILIBRATE_DISTIBUTION],
                t_cpu_max[EQUILIBRATE_DISTIBUTION]);
    PDM_printf( "PDM_gnum_from_hv timer : First exchange prepare (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[FIRST_EXCHANGE_PREPARE],
                t_cpu_max[FIRST_EXCHANGE_PREPARE]);
    PDM_printf( "PDM_gnum_from_hv timer : First exchange (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[FIRST_EXCHANGE],
                t_cpu_max[FIRST_EXCHANGE]);
    PDM_printf( "PDM_gnum_from_hv timer : Second exchange prepare (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[SECOND_EXCHANGE_PREPARE],
                t_cpu_max[SECOND_EXCHANGE_PREPARE]);
    PDM_printf( "PDM_gnum_from_hv timer : Second exchange (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[SECOND_EXCHANGE],
                t_cpu_max[SECOND_EXCHANGE]);
    PDM_printf( "PDM_gnum_from_hv timer : block sort (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BLOCK_SORT],
                t_cpu_max[BLOCK_SORT]);
    PDM_printf( "PDM_gnum_from_hv timer : block equal (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BLOCK_EQUAL],
                t_cpu_max[BLOCK_EQUAL]);
    PDM_printf( "PDM_gnum_from_hv timer : Reverse exchange (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[REVERSE_EXCHANGE],
                t_cpu_max[REVERSE_EXCHANGE]);
    PDM_printf_flush();
  }


}


/**
 *
 * \brief Free
 *
 * \param [in]   gnum_from_hv Pointer to \ref PDM_gnum_from_hv_t object
 * \param [in]   partial      1 to free partially, 0 else
 *
 */

void
PDM_gnum_from_hv_free
(
 PDM_gnum_from_hv_t *gnum_from_hv
)
{
  _pdm_gnum_from_hv_t *_gnum_from_hv = (_pdm_gnum_from_hv_t *) gnum_from_hv;

  if(( _gnum_from_hv->owner == PDM_OWNERSHIP_KEEP ) ||
     ( _gnum_from_hv->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !_gnum_from_hv->tag_results_get)) {
    for (int i = 0; i < _gnum_from_hv->n_part; i++) {
      free (_gnum_from_hv->g_nums[i]);
    }
  }

  free (_gnum_from_hv->g_nums);
  free (_gnum_from_hv->n_elts);
  free (_gnum_from_hv->part_hkeys);
  free (_gnum_from_hv->part_hstri);
  free (_gnum_from_hv->part_hdata);
  free (_gnum_from_hv->distribution);

  PDM_timer_free(_gnum_from_hv->timer);

  free (_gnum_from_hv);
  _gnum_from_hv = NULL;
}


// void
// PDM_generate_global_id_from
// (
//  const int              blk_size,
//  const unsigned char   *blk_data,
//  const int             *blk_stri,
//  gnum_from_hv_compare   fcompare,
//  gnum_from_hv_equal     fequal,
//  PDM_g_num_t          **gnum
// )
// {
//   printf(" TODO \n");
//   abort();
// }



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
