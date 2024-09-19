#ifndef __PDM_GNUM_LOCATION_PRIV_H__
#define __PDM_GNUM_LOCATION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_timer.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _pdm_gnum_location_t
 * \brief  Define a global numbering location structure
 *
 */

struct _pdm_gnum_location_t {
  PDM_MPI_Comm        comm;         /*!< Communicator */
  int                 n_part_in;    /*!< Number of local partitions  */
  int                 n_part_out;   /*!< Number of local partitions
                                        for requested locations */
  int                *n_elts_in;    /*!< Number of elements of each partition */
  const PDM_g_num_t **g_nums_in;    /*!< Global numbering  */
  int                *n_elts_out;   /*!< Number of elements requesting location */
  const PDM_g_num_t **g_nums_out;   /*!< Global numbering of elements requesting location */
  int               **location_idx; /*!< Location index of elements requesting location */
  int               **location;     /*!< Location of elements requesting location */

  PDM_ownership_t owner;       /*!< Ownership */
  int  tag_results_get ;       /*!< Tag call to PDM_dist_cellcenter_surf_get function */ 

} ;



#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_GNUM_LOCATION_PRIV_H__ */
