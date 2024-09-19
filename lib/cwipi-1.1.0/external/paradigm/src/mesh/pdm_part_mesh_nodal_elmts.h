/*
 * \file
 */

#ifndef __PDM_PART_MESH_NODAL_ELMTS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
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

typedef struct _pdm_part_mesh_nodal_elmts_t PDM_part_mesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);


int
PDM_part_mesh_nodal_elmts_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt
);

void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

void
PDM_part_mesh_nodal_elmts_std_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
      PDM_ownership_t              owner
);

void
PDM_part_mesh_nodal_elmts_std_ho_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
const int                          order,
const char                        *ho_ordering,
      PDM_ownership_t              owner
);

/**
 * \brief Define a polygon block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_section_poly2d_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
const int                          n_elt,
const int                         *connec_idx,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
      PDM_ownership_t              owner
);

/**
 * \brief Define a polyhedra block
 *
 * \param [in]  mesh           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_section_poly3d_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
const int                          n_elt,
const int                          n_face,
const int                         *facvtx_idx,
const int                         *facvtx,
const PDM_g_num_t                 *face_ln_to_gn,
const int                         *cellfac_idx,
const int                         *cellfac,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
      PDM_ownership_t              owner
);



void
PDM_part_mesh_nodal_elmts_section_std_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_section,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      PDM_ownership_t               ownership
);

void
PDM_part_mesh_nodal_elmts_section_std_ho_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_section,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      int                          *order,
const char                        **ho_ordering,
      PDM_ownership_t               ownership
);

/**
 * \brief Return a polygon block description
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  ownership      Who owns the getted arrays?
 */

void
PDM_part_mesh_nodal_elmts_section_poly2d_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_section,
 const int                           id_part,
       int                         **connec_idx,
       int                         **connec,
       PDM_ownership_t               ownership
);

/**
 * \brief Get the cell-vertex connectivity of a polyhedra block
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] cell_vtx_idx   Index of cell vertex connectivity
 * \param [out] cell_vtx       Cell vertex connectivity
 * \param [in]  ownership      Who owns the getted arrays?
 *
 */

void
PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_section,
 const int                           id_part,
       int                         **cell_vtx_idx,
       int                         **cell_vtx,
       PDM_ownership_t               ownership
);


void
PDM_part_mesh_nodal_elmts_section_poly3d_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_section,
const int                           id_part,
      int                          *n_face,
      PDM_g_num_t                 **face_ln_to_gn,
      int                         **face_vtx_idx,
      int                         **face_vtx,
      PDM_g_num_t                 **numabs,
      int                         **cell_face_idx,
      int                         **cell_face,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      PDM_ownership_t               ownership
);

int
PDM_part_mesh_nodal_elmts_section_n_elt_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part
);

int
PDM_part_mesh_nodal_elmts_n_section_get
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

int *
PDM_part_mesh_nodal_elmts_sections_id_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
);

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_elmts_section_type_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section
);

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create_from_part3d
(
  const int                n_part,
  const int               *n_cell,
  const int               *n_face,
  const int              **face_vtx_idx,
  const int              **face_vtx,
  const PDM_g_num_t      **face_ln_to_gn,
  const int              **cell_face_idx,
  const int              **cell_face,
  const double           **vtx_coord,
  const PDM_g_num_t      **numabs,
        PDM_MPI_Comm       comm
);

int *
PDM_part_mesh_nodal_elmts_parent_num_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
      PDM_ownership_t              ownership
);

PDM_g_num_t *
PDM_part_mesh_nodal_elmts_g_num_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_section,
const int                          id_part,
      PDM_ownership_t              ownership
);

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create_from_part2d
(
  const int                n_part,
  const int               *n_face,
  const int               *n_edge,
  const int               *n_vtx,
  const int              **edge_vtx_idx,
  const int              **edge_vtx,
  const int              **face_edge_idx,
  const int              **face_edge,
  const PDM_g_num_t      **numabs,
        PDM_MPI_Comm       comm
);

/**
 * \brief Compute element extents of a part of a section
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  tolerance      Expansion tolerance for bounding boxes
 * \param [in]  vtx_coord      Coordinates of vertices
 * \param [out] extents        Extents of mesh elements in current part of current block
 *
 */

void
PDM_part_mesh_nodal_elmts_elt_extents_compute
(
       PDM_part_mesh_nodal_elmts_t *pmne,
 const int                          id_section,
 const int                          id_part,
 const double                       tolerance,
       double                      *vtx_coord,
       double                      *extents
 );

/**
 * \brief Compute element centers of a part of a section
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  vtx_coord      Coordinates of vertices
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_nodal_elmts_elt_center_compute
(
       PDM_part_mesh_nodal_elmts_t *pmne,
 const int                          id_section,
 const int                          id_part,
 const int                          n_vtx,
       double                      *vtx_coord,
 const PDM_ownership_t              ownership
 );

/**
 * \brief Compute element centers of a part of a section
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  vtx_coord      Coordinates of vertices
 * \param [in]  ownership      Who owns the getted arrays?
 *
 */

const double *
PDM_part_mesh_nodal_elmts_elt_center_get
(
       PDM_part_mesh_nodal_elmts_t *pmne,
 const int                          id_section,
 const int                          id_part,
       PDM_ownership_t              ownership
);

/**
 * \brief Reset element centers of a part of a section
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 *
 */

void
PDM_part_mesh_nodal_elmts_elt_center_reset
(
       PDM_part_mesh_nodal_elmts_t *pmne,
 const int                          id_section,
 const int                          id_part
 );

/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 *
 * \return      NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_reset
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

/**
 * \brief  Compute a global numbering in a section
 *
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_section   Section identifier
 * \param [in]  ownership    Ownership
 *
 */

void
PDM_part_mesh_nodal_elmts_g_num_in_section_compute
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_section,
const PDM_ownership_t               ownership
);

/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  pmne      Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Return number elements of a partition
 *
 */

int
PDM_part_mesh_nodal_elmts_n_elmts_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part
);

/**
 * \brief Get the element global numbering taking into account parent_num
 *
 * \param [in]  pmne      Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_part   Partition identifier
 * \param [in]  ownership Who owns the getted arrays?
 *
 * \return  Global ids of element in current partition
 *
 */

PDM_g_num_t *
PDM_part_mesh_nodal_elmts_g_num_get_from_part
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part,
      PDM_ownership_t               ownership
);

/**
 * \brief Free partially a part_mesh_nodal_elmts structure
 *
 * \param [in]  pmne      Pointer to \ref PDM_part_mesh_nodal_elmts object
 *
 * \return      NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_partial_free
(
PDM_part_mesh_nodal_elmts_t *pmne
);

/**
 * \brief Get global element numbering of block elements inside the block
 *
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_section   Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Who owns the getted arrays?
 *
 * \return      Return global numbering of block elements inside the block
 *
 */

PDM_g_num_t *
PDM_part_mesh_nodal_elmts_section_g_num_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_section,
const int                           id_part,
      PDM_ownership_t               ownership
);

/**
 * \brief  Return parent element number to local number
 *
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent element number to local number
 *
 */

int *
PDM_part_mesh_nodal_elmts_num_elmt_parent_to_local_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part
);

/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  face_ln_to_gn  Face global numbering
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  cell_ln_to_gn  Global numbering
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_elmts_nodal_cell3d_cellface_add
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part,
const int                           n_cell,
const int                           n_face,
const int                          *face_vtx_idx,
const int                          *face_vtx,
const PDM_g_num_t                  *face_ln_to_gn,
const int                          *cell_face_idx,
const int                          *cell_face,
const PDM_g_num_t                  *cell_ln_to_gn,
      PDM_Mesh_nodal_vtx_t        **vtx,
const PDM_ownership_t               ownership
);

/**
 * \brief  Add some 2D faces from face edge conectivity.
 *
 * For each face, this function searchs the type of the face (triangles, quandrangles, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polyhedra
 * \param [in]  n_edge         Number of edges used to describe polyhedra
 * \param [in]  edge_vtx       edge vertex connectivity
 * \param [in]  face_edge_idx  Index of face edge connectivity
 * \param [in]  face_edge      face edge connectivity
 * \param [in]  face_ln_to_gn  Global numbering
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_nodal_elmts_face2d_faceedge_add
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part,
const int                           n_face,
const int                           n_edge,
const int                          *edge_vtx,
const int                          *face_edge_idx,
const int                          *face_edge,
const PDM_g_num_t                  *face_ln_to_gn,
const int                           n_vtx,
const PDM_ownership_t               ownership
);

/**
 * \brief  Add some standard 3D cells from cell vertex conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_cell         Number of cells
 * \param [in]  cell_vtx_idx   Index of cell vertex connectivity
 * \param [in]  cell_vtx       Cell vertex connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_nodal_elmts_cells_cellvtx_add
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part,
const int                           n_cell,
const int                          *cell_vtx_idx,
const int                          *cell_vtx,
const PDM_g_num_t                  *numabs,
const PDM_ownership_t               ownership
);

/**
 * \brief  Add some 2D faces from face vertex connectivity.
 *
 * For each face, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygon
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each edge
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_nodal_elmts_faces_facevtx_add
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_part,
const int                           n_face,
const int                          *face_vtx_idx,
const int                          *face_vtx,
const PDM_g_num_t                  *numabs,
const PDM_ownership_t               ownership
);

/**
 * \brief  Adapt pmne to fit operation communicator.
 *
 * pmne has its own communicator which might be encompassed into
 * the operation communicator.
 *
 * \param [in]  comm           Operation communicator
 * \param [in]  n_part         Number of partitions
 * \param [in]  pmne           Part mesh nodal element
 *
 */

void
PDM_part_mesh_nodal_elmts_extend_to_encompassing_comm
(
 const PDM_MPI_Comm                  comm,
 const int                           n_part,
       PDM_part_mesh_nodal_elmts_t **pmne
);

void
PDM_part_mesh_nodal_elmts_group_get
(
       PDM_part_mesh_nodal_elmts_t   *pmne,
 const int                            i_part,
 const int                            i_group,
       int                           *n_group_elmt,
       int                          **group_elmt,
       PDM_g_num_t                  **group_ln_to_gn,
       PDM_ownership_t                ownership_group
);

int*
PDM_part_mesh_nodal_elmts_compute_sections_idx
(
 PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                     id_part
);


void
PDM_part_mesh_nodal_elmts_group_set
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           i_part,
 const int                           i_group,
       int                           n_group_elmt,
       int                          *group_elmt,
       PDM_g_num_t                  *group_ln_to_gn
);

void
PDM_part_mesh_nodal_elmts_n_group_set
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           n_group,
       PDM_ownership_t               ownership_group
);


int
PDM_part_mesh_nodal_elmts_n_group_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_H__ */
