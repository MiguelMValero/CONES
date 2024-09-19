#ifndef __PDM_DMESH_NODAL_PRIV_H__
#define __PDM_DMESH_NODAL_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_dmesh_nodal_elmts_priv.h"

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
 * \struct PDM_Mesh_nodal_som_t
 * \brief  Vertices of a mesh partition
 *
 */

typedef struct PDM_DMesh_nodal_vtx_t PDM_DMesh_nodal_vtx_t;

struct PDM_DMesh_nodal_vtx_t {
  PDM_l_num_t       n_vtx;          /*!< Number of vertices */
  PDM_real_t       *_coords;        /*!< Coordinates
                                       * (Memory mapping) (size = 3 * \ref n_vtx) */
  PDM_g_num_t      *distrib;        /*!< Distribution on the processes
                                       * (size = \ref n_rank + 1) */

  int              *dvtx_tag;
  PDM_g_num_t      *dvtx_parent_g_num;

  PDM_ownership_t   owner;
};

/**
 * \struct  PDM_Mesh_nodal_geom_prepa_sections_t
 *
 * \brief   Used to build sections from cell to face face to edge connectivity
 *
 */
struct _pdm_dmesh_nodal_t {

  PDM_MPI_Comm           comm;                     /*!< MPI Communicator */
  int                    n_rank;                   /*!< Number of processes */
  int                    i_rank;                   /*!< Number of processes */
  int mesh_dimension;                              /*! Principal dimension of meshes */

  PDM_g_num_t            n_cell_abs;               /*!< Global number of elements */
  PDM_g_num_t            n_face_abs;               /*!< Global number of faces    */
  PDM_g_num_t            n_edge_abs;               /*!< Global number of edges    */
  PDM_g_num_t            n_vtx_abs;                /*!< Global number of vertices */

  PDM_DMesh_nodal_vtx_t *vtx;                      /*!< Description des sommmets de chaque partition */

  PDM_dmesh_nodal_elmts_t* volumic;
  PDM_dmesh_nodal_elmts_t* surfacic;
  PDM_dmesh_nodal_elmts_t* ridge;
  PDM_dmesh_nodal_elmts_t* corner;

  PDM_bool_t   is_computed_g_extents;
  double       g_extents[6];
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_PRIV_H__ */
