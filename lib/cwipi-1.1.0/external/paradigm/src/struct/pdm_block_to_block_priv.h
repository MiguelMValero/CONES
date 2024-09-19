/*
 * File:   pdm_block_to_part_priv.h
 * Author: equemera
 *
 * Created on April 14, 2016, 8:16 AM
 */

#ifndef PDM_BLOCK_TO_BLOCK_PRIV_H
#define	PDM_BLOCK_TO_BLOCK_PRIV_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
// #include "pdm_block_to_block_priv.h"

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
 * \struct _PDM_block_to_block_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

typedef struct {

  PDM_MPI_Comm   comm;                    /*!< MSG communicator     */
  int            n_rank;                  /*!< Communicator size    */
  int            i_rank;                  /*!< Current rank in comm */

  PDM_g_num_t    *block_distrib_ini_idx;  /*!< Block distribution initiale
                                           * (size : \ref size of \ref comm + 1) */
  PDM_g_num_t    *block_distrib_end_idx;  /*!< Block distribution final
                                           * (size : \ref size of \ref comm + 1) */

} _pdm_block_to_block_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_BLOCK_PRIV_H */

