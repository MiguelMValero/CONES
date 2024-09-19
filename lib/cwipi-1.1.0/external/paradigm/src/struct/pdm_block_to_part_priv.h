/*
 * File:   pdm_block_to_part_priv.h
 * Author: equemera
 *
 * Created on April 14, 2016, 8:16 AM
 */

#ifndef PDM_BLOCK_TO_PART_PRIV_H
#define	PDM_BLOCK_TO_PART_PRIV_H

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
 * \struct _PDM_block_to_part_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

struct _pdm_block_to_part_t {
  PDM_g_num_t   *block_distrib_idx;   /*!< Block distribution
                                       * (size : \ref size of \ref comm + 1) */
  PDM_MPI_Comm   comm;                /*!< MSG communicator */
  int            n_rank;              /*!< Communicator size */
  int            i_rank;              /*!< Current rank in comm */
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
  int            n_elt_partial_block;
  int           *idx_partial;         /*!< In case of partial block - contains index of non void idx or -1 if no data in current block
                                       * (size : requestd_data_idx[n_rank - 1]
                                             +requestd_data_n[n_rank - 1] ) */
  int           pttopt_comm;          /*!< Use point to point communication if pttopt_comm == 1 */

} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_PART_PRIV_H */
