/*
 * \file
 */

#ifndef __PDM_PARTITIONING_NODAL_ALGORITHM_H__
#define __PDM_PARTITIONING_NODAL_ALGORITHM_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_part_mesh_nodal_elmts_t*
PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           n_part,
 int                          *pn_vtx,
 PDM_g_num_t                 **vtx_ln_to_gn,
 int                          *pn_elmt,
 int                         **pelmt_to_entity,
 PDM_g_num_t                 **elmt_ln_to_gn,
 PDM_g_num_t                 **pparent_entitity_ln_to_gn
);


PDM_dmesh_nodal_elmts_t*
PDM_dmesh_nodal_elmts_to_extract_dmesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           dn_elmt,
 PDM_g_num_t                  *delmt_selected,
 PDM_g_num_t                 **extract_vtx_distribution,
 PDM_g_num_t                 **extract_parent_vtx_gnum,
 PDM_g_num_t                 **extract_entity_distribution,
 PDM_g_num_t                 **extract_parent_elt_gnum
);

void
PDM_reverse_dparent_gnum
(
       PDM_g_num_t    *dparent_gnum,
       int            *dparent_sign,
       PDM_g_num_t    *delmt_child_distrib,
       int             n_part,
       int            *pn_parent,
       PDM_g_num_t   **pparent_gnum,
       int           **pn_child,
       int          ***pelmt_to_entity,
       PDM_g_num_t  ***pchild_gnum,
       PDM_g_num_t  ***pchild_parent_gnum,
       int          ***pchild_parent_sign,
 const PDM_MPI_Comm    comm
);

void
PDM_generate_ho_vtx_ln_to_gn
(
 PDM_dmesh_nodal_t      *dmn,
 int                     n_part,
 int                    *pn_cell,
 PDM_g_num_t           **pcell_ln_to_gn,
 int                    *pn_face,
 PDM_g_num_t           **pface_ln_to_gn,
 int                    *pn_edge,
 PDM_g_num_t           **pedge_ln_to_gn,
 int                    *pn_vtx,
 PDM_g_num_t           **pvtx_ln_to_gn,
 int                   **pn_vtx_all,
 PDM_g_num_t          ***vtx_all_ln_to_gn
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_NODAL_ALGORITHM_H__ */
