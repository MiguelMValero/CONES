#ifndef __PDM_GLOBAL_MEAN_PRIV_H__
#define __PDM_GLOBAL_MEAN_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_surf_mesh.h"
#include "pdm_mesh_nodal.h"
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_global_point_mean_t
 * \brief  Define a global point mean
 *
 */

struct _pdm_global_point_mean_t {

  int          n_part;            /*!< Number of partitions */
  PDM_MPI_Comm comm;              /*!< MPI communicator */
  int          *n_elts;           /*!< Number of elements in partitions */
  PDM_g_num_t **g_nums;           /*!< Global numbering of elements */
  PDM_part_to_block_t *ptb;       /*!< Part to block structure */
  PDM_block_to_part_t *btp;       /*!< Block to part structure */
  int           stride;           /*!< Current Field stride */
  int         **strides;          /*!< Strides array storage
                                   *   (In ptb strides are variable) */
  double      **local_field;      /*!< Local field */
  double       *s_weight;         /*!< Sume of weights */
  double      **local_weight;     /*!< Weight */
  double      **global_mean_field;/*!< Global mean field */

} ;


/*=============================================================================
 * Static global variables
 *============================================================================*/
#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_GLOBAL_MEAN_PRIV_H__ */
