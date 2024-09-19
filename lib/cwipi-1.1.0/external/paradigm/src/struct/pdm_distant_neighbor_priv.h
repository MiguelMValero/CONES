#ifndef __PDM_PART_DISTANT_NEIGHBOR_PRIV_H__
#define __PDM_PART_DISTANT_NEIGHBOR_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _distant_neighbor_t
 * \brief  Define a point merge structures
 *
 */
struct _pdm_distant_neighbor_t {
  PDM_MPI_Comm   comm;                  /*!< MPI communicator */
  int            n_part;                /*!< Number of partitions */
  const int     *n_entity;              /*!< Number of entities for each partition */
  int          **neighbor_idx;          /*!< Indexes of candidate for each current part point
                                         *   (size = number of entities in the current part + 1) */
  int          **neighbor_desc;         /*!< Candidates description (process,
                                         *                           part in the process,
                                         *                           entitiy number in the part) */
  int**          order;
  int**          order_unique;
  int           *requested_data_n;      /*!< Numer of requested data for each process index
                                         * (size : s_comm) */
  int           *requested_data_idx;    /*!< Requested data for each process index
                                         * (size : s_comm) */
  int           *distributed_data_n;    /*!< Numer of distributed data for each process index
                                         * (size : s_comm) */
  int           *distributed_data_idx;  /*!< Distributed data for each process index
                                         * (size : s_comm) */
  int           *distributed_data;      /*!< Distributed data for each process
                                         * (size : requestd_data_idx[s_comm - 1] */
  int           *distributed_part_n;
  int           *distributed_part_idx; /*!< For each part the shift to apply on recv buffer
                                       * (size : n_partà )*/

  int         **ind;

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_DISTANT_NEIGHBOR_PRIV_H__ */
