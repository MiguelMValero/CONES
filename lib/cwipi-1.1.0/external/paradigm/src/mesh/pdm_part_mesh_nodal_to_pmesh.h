#ifndef __PDM_PART_MESH_NODAL_ELMTS_TO_PART_MESH_H__
#define __PDM_PART_MESH_NODAL_ELMTS_TO_PART_MESH_H__

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
#include "pdm_mesh_nodal_priv.h"

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

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

PDM_part_mesh_t*
PDM_part_mesh_nodal_to_part_mesh
(
        PDM_part_mesh_nodal_t                      *pmn,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_TO_PART_MESH_H__ */
