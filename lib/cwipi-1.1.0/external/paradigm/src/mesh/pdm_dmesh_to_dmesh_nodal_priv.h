#ifndef __PDM_DMESH_TO_DMESH_NODAL_PRIV_H__
#define __PDM_DMESH_TO_DMESH_NODAL_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_priv.h"

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
 * \struct _pdm_dmesh_to_dmesh_nodal_t
 * \brief  Transform a dmesh into dmesh_nodal
 *
 */

struct _pdm_dmesh_to_dmesh_nodal_t {

  PDM_MPI_Comm         comm;                    /*!< MPI communicator */
  PDM_ownership_t      owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t           results_is_getted;       /*!< Flags to indicate if result is getted      */
  int                  n_mesh;                  /*!< Number of meshes to manages                */

  /* Distrib */
  PDM_g_num_t         **distrib_cell;
  PDM_g_num_t         **distrib_face;
  PDM_g_num_t         **distrib_edge;
  PDM_g_num_t         **distrib_vtx ;


  /* Stored connectivity */
  PDM_g_num_t         **dcell_face;
  int                 **dcell_face_idx;
  PDM_g_num_t         **dface_edge;
  int                 **dface_edge_idx;
  PDM_g_num_t         **dface_vtx;
  int                 **dface_vtx_idx;
  PDM_g_num_t         **dedge_vtx;
  double              **dvtx_coords;
  int                 **dparent_elmt_position;

  int                  **n_bound;
  int                 ***dbound_idx;
  PDM_g_num_t         ***dbound;

  /* Result */
  PDM_dmesh_nodal_t    **dmn;
  int                  **n_blk_gnum;
  PDM_g_num_t         ***blk_entity_gnum;
  PDM_g_num_t         ***blk_elmt_gnum;

};

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_DMESH_TO_DMESH_NODAL_PRIV_H__ */
