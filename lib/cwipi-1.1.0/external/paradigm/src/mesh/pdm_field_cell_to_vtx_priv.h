#ifndef __PDM_FIELD_CELL_TO_VTX_PRIV_H__
#define __PDM_FIELD_CELL_TO_VTX_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

// #include "pdm_multipart.h"
#include "pdm_part_priv.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_mesh_nodal.h"

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
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_field_cell_to_vtx_t
 * \brief  Distributed cube
 *
 * _pdm_field_cell_to_vtx_t define a distributed mesh of a cube
 *
 */

struct _pdm_field_cell_to_vtx_t {
  PDM_MPI_Comm      comm;            /*!< MPI communicator                          */

  int                           n_depth;
  PDM_cell_to_vtx_interp_kind_t interp_kind;
  int                           idw_p;

  int               n_domain;
  int              *n_group;
  int              *n_part;
  int              *n_part_idx;
  int              *n_part_g_idx;


  PDM_part_domain_interface_t  *pdi;
  int                           n_interface;

  _part_t               **parts;
  PDM_part_mesh_nodal_t **pmn;
  int is_nodal;

  double **cell_center;

  int   n_part_loc_all_domain;
  int  *pn_vtx;
  int **neighbor_idx;
  int **neighbor_desc;
  int **neighbor_interface;

  int **pvtx_cell_n;
  int **pvtx_cell_idx;
  int **pvtx_cell;

  int    **vtx_face_bound_idx;
  int    **vtx_face_bound_n;
  int    **vtx_face_bound;
  int    **vtx_face_bound_group;
  double **vtx_face_bound_coords;

  PDM_distant_neighbor_t* dn;


  /* Graphe communication inside domain */
  int ****entity_part_bound_proc_idx;
  int ****entity_part_bound_part_idx;
  int ****entity_part_bound;
  int    *graph_comm_is_defined;

  int *****group_entity;     // (i_kind, i_domain, i_part, i_group)
  int  ****n_group_entity;
  int     *group_is_defined;

  /* Save graphe of  border */
  int    **pvtx_cell_coords_opp_n;
  double **pvtx_cell_coords_opp;
  int    **pvtx_face_bound_coords_opp_n;
  double **pvtx_face_bound_coords_opp;

  double ***user_cell_center;
  double ***user_vtx_center;

};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_field_cell_to_vtx_PRIV_H__ */
