/*
 * \file
 */

#ifndef __PDM_DMESH_NODAL_H__
#define __PDM_DMESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_elmts.h"
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

typedef struct _pdm_dmesh_nodal_t      PDM_dmesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_dmesh_nodal_t*
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm,
      int          mesh_dimension,
      PDM_g_num_t  n_vtx,
      PDM_g_num_t  n_cell,
      PDM_g_num_t  n_face,
      PDM_g_num_t  n_edge
);


void
PDM_DMesh_nodal_free
(
 PDM_dmesh_nodal_t* dmesh_nodal
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
       PDM_dmesh_nodal_t *dmesh_nodal,
 const int                n_vtx,
       PDM_real_t        *coords,
       PDM_ownership_t    owner
);

void
PDM_DMesh_nodal_vtx_tag_set
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 int               *dvtx_tag
);

void
PDM_DMesh_nodal_vtx_parent_gnum_set
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 PDM_g_num_t       *dvtx_parent_g_num
);

int*
PDM_DMesh_nodal_vtx_tag_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
);


PDM_g_num_t *
PDM_DMesh_nodal_vtx_parent_gnum_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
);


void
PDM_DMesh_nodal_section_g_dims_get
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  PDM_g_num_t       *n_cell_abs,
  PDM_g_num_t       *n_face_abs,
  PDM_g_num_t       *n_edge_abs,
  PDM_g_num_t       *n_vtx_abs
);


/**
 * \brief  Return vertices distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
(
PDM_dmesh_nodal_t *dmesh_nodal
);


/**
 * \brief  Return section distribution
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
);


/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
  PDM_dmesh_nodal_t *dmesh_nodal
);


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

double *
PDM_DMesh_nodal_vtx_get
(
  PDM_dmesh_nodal_t *dmesh_nodal
);


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_section_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
);


/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_DMesh_nodal_sections_id_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
);

/**
 * \brief  Return type of element of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_elt_type_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
);


/**
 * \brief  Return type of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
);

/**
 * \brief  Return distri of section
 *
 * \param [in] dmesh_nodal
 * \param [in]  id_section   Block identifier
 *
 * \return  distri
 *
 */
PDM_g_num_t*
PDM_DMesh_nodal_section_distri_std_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
);

/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory
 *                             when the section is destroyed
 * \param [in]  id_section       Block identifier
 *
 * \return Block identifier
 *
 */

int
PDM_DMesh_nodal_section_add
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const PDM_Mesh_nodal_elt_t  t_elt
);

void
PDM_DMesh_nodal_update_ownership
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_ownership_t      owner
);


/**
 * \brief Define a standard section
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */

void
PDM_DMesh_nodal_section_std_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const int                  n_elt,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
);

void
PDM_DMesh_nodal_section_group_elmt_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  n_group_elmt,
      int                 *dgroup_elmt_idx,
      PDM_g_num_t         *dgroup_elmt,
      PDM_ownership_t      owner
);

void
PDM_DMesh_nodal_section_group_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind,
 int                 *n_group_elmt,
 int                 **dgroup_elmt_idx,
 PDM_g_num_t         **dgroup_elmt
);


/**
 * \brief Return standard section description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRISM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_section_std_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
);

PDM_g_num_t *
PDM_DMesh_nodal_section_std_ho_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
      int                 *order,
const char               **ho_ordering
);


/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_DMesh_nodal_section_n_elt_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
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
PDM_DMesh_nodal_section_poly2d_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
      PDM_l_num_t         *connec_idx,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
);


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
      PDM_l_num_t        **connec_idx,
      PDM_g_num_t        **connec
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
PDM_DMesh_nodal_section_poly3d_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
const PDM_l_num_t          n_face,
      PDM_l_num_t         *facvtx_idx,
      PDM_g_num_t         *facvtx,
      PDM_l_num_t         *cellfac_idx,
      PDM_g_num_t         *cellfac,
      PDM_ownership_t      owner
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
PDM_DMesh_nodal_section_poly3d_get
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const int                   id_section,
      PDM_l_num_t          *n_face,
      PDM_l_num_t         **facvtx_idx,
      PDM_g_num_t         **facvtx,
      PDM_l_num_t         **cellfac_idx,
      PDM_g_num_t         **cellfac
);


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 *
 */

PDM_g_num_t
PDM_dmesh_nodal_total_n_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
);

/**
 * \brief  Return vtx distribution of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return vtx distribution
 *
 */
PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
);


PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_copy_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
);

/**
 * \brief  Return total number of vertices of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return total number of vertices
 *
 */
PDM_g_num_t
PDM_dmesh_nodal_total_n_vtx_get
(
PDM_dmesh_nodal_t *dmesh_nodal
);


/**
 *
 * \brief Setup global distribution of all elements register in current structure
 *
 * \param [inout]  mesh
 *
 * \return         Null
 *
 */
void
PDM_dmesh_nodal_generate_distribution
(
PDM_dmesh_nodal_t *dmesh_nodal
);

/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     hdl                Distributed nodal mesh handle
* \param [inout]  n_face_elt_tot     Number of faces
* \param [inout]  n_sum_vtx_face_tot Number of vtx for all faces (cumulative)
*
*/
void
PDM_dmesh_nodal_decompose_faces_get_size
(
PDM_dmesh_nodal_t *dmesh_nodal,
int               *n_face_elt_tot,
int               *n_sum_vtx_face_tot
);

/**
*
* \brief Concatenates all element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_elt_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **section_idx,
  int               **cat_delt_vtx_idx,
  PDM_g_num_t       **cat_dcell_vtx
);

/**ss
*
* \brief Concatenates 3D element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_cell_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **section_idx,
  int               **cat_delt_vtx_idx,
  PDM_g_num_t       **cat_delt_vtx
);

/**
 * \brief  Compute cell->cell connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 *
 */
void
PDM_dmesh_nodal_dual_graph
(
  PDM_g_num_t*   vtx_dist,
  PDM_g_num_t*   elt_dist,
  int           *delt_vtx_idx,
  PDM_g_num_t   *delt_vtx,
  PDM_g_num_t  **dual_graph_idx,
  PDM_g_num_t  **dual_graph,
  PDM_MPI_Comm   comm
);

void
PDM_dmesh_nodal_transfer_to_new_dmesh_nodal
(
 PDM_dmesh_nodal_t   *dmn_in,
 PDM_dmesh_nodal_t   *dmn_out,
 PDM_geometry_kind_t  geom_kind,
 PDM_g_num_t         *dparent_vtx_distrib,
 PDM_g_num_t         *blk_parent_to_new_vtx_gnum
);

void
PDM_dmesh_nodal_transfer_to_new_dmesh_nodal_gen
(
 PDM_dmesh_nodal_t   *dmn_in,
 PDM_dmesh_nodal_t   *dmn_out,
 PDM_geometry_kind_t  geom_kind,
 PDM_g_num_t         *dparent_vtx_distrib,
 int                 *blk_parent_to_new_vtx_gnum_idx,
 PDM_g_num_t         *blk_parent_to_new_vtx_gnum
);

void
PDM_dmesh_nodal_dump_vtk
(
       PDM_dmesh_nodal_t   *dmn,
       PDM_geometry_kind_t  geom_kind,
 const char                *filename_patter
);

void
PDM_dmesh_nodal_reorder
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 const char        *ordering_name
 );



PDM_part_mesh_nodal_elmts_t*
PDM_dmesh_nodal_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_t            *dmn,
 PDM_geometry_kind_t           geom_kind,
 int                           n_part,
 int                          *pn_vtx,
 PDM_g_num_t                 **vtx_ln_to_gn,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn,
 PDM_g_num_t                 **pparent_entitity_ln_to_gn
);


const double *
PDM_dmesh_nodal_global_extents_get
(
 PDM_dmesh_nodal_t         *dmn
);


int
PDM_dmesh_nodal_have_ho
(
 PDM_dmesh_nodal_t         *dmn
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
