#ifndef __PDM_PART_MESH_NODAL_PRIV_H__
#define __PDM_PART_MESH_NODAL_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"

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
 * Type
 *============================================================================*/


/**
 * \struct  PDM_Mesh_nodal_geom_prepa_blocks_t
 *
 * \brief   Used to build blocks from cell to face face to edge connectivity
 *
 */

struct _pdm_part_mesh_nodal_t {

  PDM_MPI_Comm                        comm;                      /*!< MPI Communicator            */
  int                                 mesh_dimension;

  int                                 n_part;           /*!< Number of partitions */

  PDM_Mesh_nodal_vtx_t              **vtx;              /*!< Description des sommmets de chaque partition */

  PDM_l_num_t                        *n_vol;
  PDM_l_num_t                        *n_surf;
  PDM_l_num_t                        *n_ridge;
  PDM_l_num_t                        *n_corner;

  PDM_part_mesh_nodal_elmts_t       *volumic;
  PDM_part_mesh_nodal_elmts_t       *surfacic;
  PDM_part_mesh_nodal_elmts_t       *ridge;
  PDM_part_mesh_nodal_elmts_t       *corner;

  int                                is_vtx_def_from_parent; /*<! Are the points defined from parents */

  int s_section;
  int n_section;
  PDM_geometry_kind_t *section_kind;
  int                 *section_id;
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_PRIV_H__ */
