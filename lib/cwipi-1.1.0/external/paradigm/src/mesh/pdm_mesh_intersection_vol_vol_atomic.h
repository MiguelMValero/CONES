#ifndef PDM_MESH_INTERSECTION_VOL_VOL_ATOMIC_H
#define PDM_MESH_INTERSECTION_VOL_VOL_ATOMIC_H

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

double
PDM_mesh_intersection_vol_vol_atomic_compute
(
 double    triaB_coord[9],
 double  **vtx_coordA,
 int      *n_vtxA,
 int     **face_vtxA,
 int      *n_faceA,
 double  **vtx_coordB,
 int      *n_vtxB,
 int     **face_vtxB,
 int      *n_faceB
);


double
PDM_mesh_intersection_vol_vol_atomic_compute2
(
 double triaB_coord[9]
);

double
PDM_mesh_intersection_vol_vol_atomic_compute3
(
 double triaB_coord[9]
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_MESH_INTERSECTION_H */
