#ifndef __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__
#define __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"

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

// typedef struct _pdm_part_mesh_nodal_elmts_t _pdm_part_mesh_nodal_elmts_t;
struct _pdm_part_mesh_nodal_elmts_t {

  PDM_MPI_Comm                         comm;                      /*!< MPI Communicator            */
  int                                  mesh_dimension;
  PDM_Mesh_nodal_prepa_blocks_t       *prepa_blocks;              /*!< Blocks preparation          */

  int                                  n_part;
  PDM_l_num_t                         *n_elmts;                   /*!< Nombre de blocs d'elements  */

  int                                  n_section;                 /*!< Total number of sections */
  int                                  n_section_std;             /*!< Total number of standard sections   */
  int                                  n_section_poly2d;          /*!< Total number of polyhedron sections */
  int                                  n_section_poly3d;          /*!< Total number of polyhedron sections */
  int                                 *sections_id;               /*!< sections identifier           */
  PDM_Mesh_nodal_block_std_t         **sections_std;              /*!< Standard sections             */
  PDM_Mesh_nodal_block_poly2d_t      **sections_poly2d;           /*!< Polygon sections              */
  PDM_Mesh_nodal_block_poly3d_t      **sections_poly3d;           /*!< Polyhedron sections           */

  /* Group */
  int                                  n_group;
  int                                **n_group_elmt;              // (i_part, i_group)
  int                               ***group_elmt;
  PDM_g_num_t                       ***group_ln_to_gn;

  PDM_l_num_t                        **num_elmt_parent_to_local;  /*!< Initial local numbering to local numbering
                                                                   *   imposed by blocks */
  PDM_g_num_t                        **numabs;                    /*!< Global numbering per elmts per partition */

  /* Ownerships */
  PDM_ownership_t                      ownership_group;          /*!< group_elmt, group_ln_to_gn */
  PDM_ownership_t                      ownership_numabs;         /*!< numabs */
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__ */
