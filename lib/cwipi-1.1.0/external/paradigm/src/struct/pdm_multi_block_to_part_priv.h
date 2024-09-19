/*
 * File:   pdm_block_to_part_priv.h
 * Author: bmaugars
 *
 * Created on November 06, 2020, 6:43 AM
 */

#ifndef PDM_MULTI_BLOCK_TO_PART_PRIV_H
#define	PDM_MULTI_BLOCK_TO_PART_PRIV_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct _PDM_multi_block_to_part_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

struct _pdm_multi_block_to_part_t {
  PDM_g_num_t   *multi_distrib_idx;   /*!< Multiple distribution
                                       * (size : \ref size of \ref comm + 1) */
  int n_block;                        /*!< Number of blocks */
  PDM_g_num_t   **block_distrib_idx;  /*!< Block distribution
                                       * (size : \ref size of \ref comm + 1) */
  PDM_MPI_Comm   comm;                /*!< MSG communicator */
  int            n_rank;              /*!< Communicator size */
  int            i_rank;              /*!< Current rank in comm */
  int            n_data_block;        /*!< Total number of blocks */
  int            n_part;              /*!< Number of partitions */
  int           *n_elt;               /*!< Number of elements for each partition */
  int          **ind;                 /*!< Ind for each element partition in distributed_data */
  int           *requested_data_n;    /*!< Numer of requested data for each process index
                                       * (size : n_rank) */
  int           *requested_data_idx;  /*!< Requested data for each process index
                                       * (size : n_rank) */
  int           *distributed_data_n;  /*!< Numer of distributed data for each process index
                                       * (size : n_rank) */
  int           *distributed_data_idx;/*!< Distributed data for each process index
                                       * (size : n_rank) */
  int           *distributed_data;    /*!< Distributed data for each process
                                       * (size : requestd_data_idx[n_rank - 1]
                                             +requestd_data_n[n_rank - 1] ) */
  int           *requested_block_n;    /*!< Numer of requested block for each process index
                                       * (size : n_block*n_rank) */
  int           *requested_block_idx;  /*!< Requested block for each process index
                                       * (size : n_block*n_rank) */
  int           *distributed_block_n;  /*!< Numer of distributed block for each process index
                                       * (size : n_block*n_rank) */
  int           *distributed_block_idx;/*!< Distributed block for each process index
                                       * (size : n_block*n_rank) */
  int           pttopt_comm;          /*!< Use point to point communication if pttopt_comm == 1 */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_MULTI_BLOCK_TO_PART_PRIV_H */
