#ifndef __PDM_DMESH_TO_DMESH_NODAL_H__
#define __PDM_DMESH_TO_DMESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_part_mesh.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_dmesh_to_dmesh_nodal_t PDM_dmesh_to_dmesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/


PDM_dmesh_to_dmesh_nodal_t*
PDM_dmesh_to_dmesh_nodal_create
(
 const int             n_mesh,
 const PDM_MPI_Comm    comm
);

void
PDM_dmesh_to_dmesh_nodal_compute
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
);

void
PDM_dmesh_to_dmesh_nodal_set_dmesh
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_t                *dm
);


void
PDM_dmesh_to_dmesh_nodal_distribution_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_g_num_t                *distrib_cell,
        PDM_g_num_t                *distrib_face,
        PDM_g_num_t                *distrib_edge,
        PDM_g_num_t                *distrib_vtx
);

void
PDM_dmesh_to_dmesh_nodal_connectivity_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        int                        *dcell_face_idx,
        PDM_g_num_t                *dcell_face,
        int                        *dface_edge_idx,
        PDM_g_num_t                *dface_edge,
        PDM_g_num_t                *dedge_vtx,
        int                        *dface_vtx_idx,
        PDM_g_num_t                *dface_vtx,
        double                     *dvtx_coords
);

void
PDM_dmesh_to_dmesh_nodal_group_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_bound_type_t            bound_type,
        int                         n_bound,
        int                        *dbound_idx,
        PDM_g_num_t                *dbound
);

void
PDM_dmesh_to_dmesh_nodal_dmesh_nodal_get
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_nodal_t         **dmn,
        PDM_ownership_t             ownership
);

void
PDM_dmesh_to_dmesh_nodal_free
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
);

void
PDM_dmesh_to_dmesh_nodal_update_group
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
 int                         i_mesh,
 PDM_bound_type_t            bound_type,
 int                         dn_elmt_bound,
 PDM_g_num_t                *dentity_bound,
 PDM_g_num_t               **delmt_bound
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_TO_DMESH_NODAL_H__ */
