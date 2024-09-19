#ifndef PDM_MESH_INTERSECTION_SURF_SURF_ATOMIC_H
#define PDM_MESH_INTERSECTION_SURF_SURF_ATOMIC_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"


/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/


double PDM_mesh_intersection_surf_surf_atomic_compute
(
 double triaA_coord[9],
 double edgeB_coord[6],
 double edgeB_normal[6]
 );

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_MESH_INTERSECTION_H */
