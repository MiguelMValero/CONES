/*
 * \file
 */

#ifndef __PDM_DMESH_NODAL_ELMTS_H__
#define __PDM_DMESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

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

typedef struct _pdm_dmesh_nodal_elts_t PDM_dmesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/


PDM_dmesh_nodal_elmts_t*
PDM_DMesh_nodal_elmts_create
(
const PDM_MPI_Comm comm,
      int          mesh_dimension,
      PDM_g_num_t  n_elmts
);

int
PDM_DMesh_nodal_elmts_section_add
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const PDM_Mesh_nodal_elt_t     t_elt
);


int
PDM_DMesh_nodal_elmts_section_ho_add
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const PDM_Mesh_nodal_elt_t     t_elt,
const int                      order,
const char                    *ho_ordering
);

void
PDM_DMesh_nodal_elmts_update_ownership
(
  PDM_dmesh_nodal_elmts_t *dmn_elts,
  PDM_ownership_t          owner
);

void
PDM_DMesh_nodal_elmts_section_std_set
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                id_section,
const int                n_elt,
      PDM_g_num_t       *connec,
      PDM_ownership_t    owner
);

void
PDM_DMesh_nodal_elmts_group_set
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                n_group_elmt,
      int               *dgroup_elmt_idx,
      PDM_g_num_t       *dgroup_elmt,
      PDM_ownership_t    owner
);

void
PDM_DMesh_nodal_elmts_group_get
(
 PDM_dmesh_nodal_elmts_t  *dmn_elts,
 int                      *n_group_elmt,
 int                     **dgroup_elmt_idx,
 PDM_g_num_t             **dgroup_elmt
);

void
PDM_DMesh_nodal_elmts_free
(
PDM_dmesh_nodal_elmts_t* dmn_elts
);


const PDM_g_num_t *
PDM_DMesh_nodal_elmts_distrib_section_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
);

void
PDM_dmesh_nodal_elmts_generate_distribution
(
 PDM_dmesh_nodal_elmts_t *dmn_elts
);


PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_elmts_section_type_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
);

PDM_g_num_t *
PDM_DMesh_nodal_elmts_section_std_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
);

PDM_g_num_t *
PDM_DMesh_nodal_elmts_section_std_ho_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      int                      *order,
const char                    **ho_ordering
);

int
PDM_DMesh_nodal_elmts_section_n_elt_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
);

/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_elmts_section_poly2d_set
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section,
const PDM_l_num_t              n_elt,
      PDM_l_num_t             *connec_idx,
      PDM_g_num_t             *connec,
      PDM_ownership_t          owner
);

void
PDM_DMesh_nodal_elmts_section_poly2d_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      PDM_l_num_t             **connec_idx,
      PDM_g_num_t             **connec
);

/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */
void
PDM_DMesh_nodal_elmts_section_poly3d_set
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
const PDM_l_num_t               n_elt,
const PDM_l_num_t               n_face,
      PDM_l_num_t              *facvtx_idx,
      PDM_g_num_t              *facvtx,
      PDM_l_num_t              *cellfac_idx,
      PDM_g_num_t              *cellfac,
      PDM_ownership_t           owner
);

/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_elmts_section_poly3d_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      PDM_l_num_t              *n_face,
      PDM_l_num_t             **facvtx_idx,
      PDM_g_num_t             **facvtx,
      PDM_l_num_t             **cellfac_idx,
      PDM_g_num_t             **cellfac
);

void
PDM_dmesh_nodal_elmts_generate_distribution
(
 PDM_dmesh_nodal_elmts_t *dmn_elts
);

PDM_g_num_t
PDM_DMesh_nodal_elmts_total_n_elmt_get
(
 PDM_dmesh_nodal_elmts_t  *dmn_elts
);


void
PDM_dmesh_nodal_elmts_decompose_faces_get_size
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
int                     *n_edge_elt_tot,
int                     *n_sum_vtx_edge_tot
);

void
PDM_dmesh_nodal_elmts_decompose_edges_get_size
(
PDM_dmesh_nodal_elmts_t *dmn_elts,
int                     *n_edge_elt_tot,
int                     *n_sum_vtx_edge_tot
);


void
PDM_DMesh_nodal_elmts_section_std_ho_reorder
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section,
const char                    *ho_ordering
);


int
PDM_dmesh_nodal_elmts_have_ho
(
 PDM_dmesh_nodal_elmts_t *dmn_elts
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_NODAL_ELMTS_H__ */
