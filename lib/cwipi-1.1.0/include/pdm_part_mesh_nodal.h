/*
 * \file
 */

#ifndef __PDM_PART_MESH_NODAL_H__
#define __PDM_PART_MESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal_elmts.h"
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

typedef struct _pdm_part_mesh_nodal_t PDM_part_mesh_nodal_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a \ref PDM_part_mesh_nodal_t structure
 *
 * \param [in]   mesh_dimension   Mesh dimension
 * \param [in]   n_part           Number of partition on the current process
 * \param [in]   comm             MPI communicator
 *
 * \return       Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */
//---> PDM_Mesh_nodal_create

PDM_part_mesh_nodal_t*
PDM_part_mesh_nodal_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 * \param [in]  owner      Vertices ownship
 *
 */
//---> PDM_Mesh_nodal_coord_set

void
PDM_part_mesh_nodal_coord_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const double                *coords,
 const PDM_g_num_t           *numabs,
       PDM_ownership_t        owner
);

/**
 * \brief Define partition vertices from parents
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part       Partition identifier
 * \param [in]  n_vtx         Number of vertices
 * \param [in]  n_vtx_parent  Number of parent vertices
 * \param [in]  numabs        Global numbering (size = \ref n_vtx)
 * \param [in]  numabs        Global numbering of parent vertices (size = \ref n_vtx_parent)
 * \param [in]  coords        Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  coords        Interlaced coordinates of parent vertices (size = 3 * \ref n_vtx_parent)
 * \param [in]  owner         Vertices ownship
 *
 */
//---> PDM_Mesh_nodal_coord_from_parent_set

void
PDM_part_mesh_nodal_coord_from_parent_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const int                    n_vtx_parent,
 const PDM_g_num_t           *numabs,
 const int                   *num_parent,
 const PDM_real_t            *coords_parent,
 const PDM_g_num_t           *numabs_parent,
 const PDM_ownership_t        ownership
);

/**
 * \brief  Return number of partitions
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return  Number of partitions
 *
 */
//---> PDM_Mesh_nodal_n_part_get

int
PDM_part_mesh_nodal_n_part_get
(
       PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return the mesh dimension
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return  mesh dimension
 *
 */

int
PDM_part_mesh_nodal_mesh_dimension_get
(
       PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return number of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */
//---> PDM_Mesh_nodal_n_vertices_get

int
PDM_part_mesh_nodal_n_vtx_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
);

/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Coordinates of vertices
 *
 */
//---> PDM_Mesh_nodal_vertices_get

double*
PDM_part_mesh_nodal_vtx_coord_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
);

/**
 * \brief  Return global ids of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Golbal ids of vertices
 *
 */
//---> PDM_Mesh_nodal_vertices_g_num_get
PDM_g_num_t*
PDM_part_mesh_nodal_vtx_g_num_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
);

/**
 * \brief  Return number of sections in a specific geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Number of sections
 *
 */

int
PDM_part_mesh_nodal_n_section_in_geom_kind_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
);

/**
 * \brief  Return ids of sections in a specific geometry kind
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Ids of sections
 *
 */

int *
PDM_part_mesh_nodal_sections_id_in_geom_kind_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
);

/**
 * \brief  Return type of section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 *
 * \return  Type of section
 *
 */
//---> PDM_Mesh_nodal_block_type_get

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section
);

/**
 * \brief  Return type of section in a specific geometry kind
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_section
);

/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  t_elt        Section type
 *
 * \return Section identifier
 *
 */
//---> PDM_Mesh_nodal_block_add

int
PDM_part_mesh_nodal_section_add
(
      PDM_part_mesh_nodal_t *pmn,
const PDM_Mesh_nodal_elt_t   t_elt
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
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  owner                   Ownership
 *
 */
//---> PDM_Mesh_nodal_block_std_set
void
PDM_part_mesh_nodal_section_std_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
      PDM_ownership_t        owner
);

/**
 * \brief Define a standard high-order section
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  order                   Element order
 * \param [in]  ho_ordering             HO ordering
 * \param [in]  owner                   Ownership
 *
 */

void
PDM_part_mesh_nodal_section_std_ho_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
const int                    order,
const char                  *ho_ordering,
      PDM_ownership_t        owner
);

/**
 * \brief Get number of section elements
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section  Section identifier
 * \param [in]  id_part    Partition identifier
 *
 * \return      Number of elements
 *
 */
//---> PDM_Mesh_nodal_block_n_elt_get
int
PDM_part_mesh_nodal_section_n_elt_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part
);

/**
 * \brief Return standard section description
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global ids
 * \param [out] parent_num              Parent local ids or NULL
 * \param [out] parent_entity_g_num     Parent global ids or NULL
 * \param [in]  ownership               Data ownership
 *
 */
//---> PDM_Mesh_nodal_block_std_get
void
PDM_part_mesh_nodal_section_std_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      PDM_ownership_t         ownership
);

/**
 * \brief Return standard high-order section description
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section               Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global ids
 * \param [out] parent_num              Parent local ids or NULL
 * \param [out] parent_entity_g_num     Parent global ids or NULL
 * \param [out] order                   Element order
 * \param [out] ho_ordering             HO ordering
 * \param [in]  ownership               Data ownership
 *
 */

void
PDM_part_mesh_nodal_section_std_ho_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      int                    *order,
const char                  **ho_ordering,
      PDM_ownership_t         ownership
);

/**
 * \brief Get parent numbering of block elements
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return parent numbering of block elements
 *
 */
//---> PDM_Mesh_nodal_block_parent_num_get
int *
PDM_part_mesh_nodal_section_parent_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
);

/**
 * \brief Get global element numbering of section elements
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return global element numbering of section elements
 *
 */
//---> PDM_Mesh_nodal_g_num_get
PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
);

/**
 * \brief Add a \ref PDM_part_mesh_nodal_elmts_t to a \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 *
 */

void
PDM_part_mesh_nodal_add_part_mesh_nodal_elmts
(
 PDM_part_mesh_nodal_t       *pmn,
 PDM_part_mesh_nodal_elmts_t *pmne
);

/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */
//---> PDM_Mesh_nodal_free
void
PDM_part_mesh_nodal_free
(
 PDM_part_mesh_nodal_t* pmn
);

void
PDM_part_mesh_nodal_dump_vtk
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind,
 const char            *filename_patter
);

/**
 * \brief Compute element extents of a part of a section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  tolerance      Expansion tolerance for bounding boxes
 * \param [out] extents        Extents of mesh elements in current part of current block
 *
 */
//---> PDM_Mesh_nodal_compute_cell_extents
void
PDM_part_mesh_nodal_section_elt_extents_compute
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    i_section,
 const int                    i_part,
 const double                 tolerance,
       double                *extents
);


/**
 * \brief Compute cell centers of a part of section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_Mesh_nodal_cell_centers_compute
void
PDM_part_mesh_nodal_section_elt_center_compute
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part,
const PDM_ownership_t        ownership
);

/**
 * \brief  Return cell centers
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Data ownership
 *
 * \return  Return cell centers
 *
 */
//---> PDM_Mesh_cell_centers_get
const double *
PDM_part_mesh_nodal_section_elt_center_get
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part,
      PDM_ownership_t        ownership
);


/**
 * \brief Reset cell center computation
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 *
 */
//---> PDM_Mesh_nodal_cell_centers_reset
void
PDM_part_mesh_nodal_section_elt_center_reset
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part
);


/**
 * \brief Define a polygon section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connec_idx     Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connec         Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */
//---> PDM_Mesh_nodal_block_poly2d_set
void
PDM_part_mesh_nodal_section_poly2d_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec_idx,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
      PDM_ownership_t        owner
);


/**
 * \brief Return a polygon section description
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connec_idx     Connectivity index
 * \param [out] connec         Connectivity
 * \param [in]  ownership      Data ownership
 *
 */
//---> PDM_Mesh_nodal_block_poly2d_get
void
PDM_part_mesh_nodal_section_poly2d_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec_idx,
      int                   **connec,
      PDM_ownership_t         ownership
);


/**
 * \brief Define a polyhedron section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section      Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  n_face         Number of faces
 * \param [in]  facvtx_idx     Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [in]  facvtx         Face->vertex connectivity
 * \param [in]  face_ln_to_gn  Face global ids
 * \param [in]  cellfac_idx    Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [in]  cellfac        Cell->face connectivity
 * \param [in]  numabs         Cell global ids
 * \param [in]  parent_num     Parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */
//---> PDM_Mesh_nodal_block_poly3d_set
void
PDM_part_mesh_nodal_section_poly3d_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                    n_face,
const int                   *facvtx_idx,
const int                   *facvtx,
const PDM_g_num_t           *face_ln_to_gn,
const int                   *cellfac_idx,
const int                   *cellfac,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
      PDM_ownership_t        owner
);

/**
 * \brief Return a polyhedron section
 *
 * \param [in]  pmn                  Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section            Section identifier
 * \param [in]  id_part              Partition identifier
 * \param [out] n_face               Number of faces
 * \param [out] face_ln_to_gn        Face global ids
 * \param [out] facvtx_idx           Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [out] facvtx               Face->vertex connectivity
 * \param [out] numabs               Cell global ids
 * \param [out] cell_face_idx        Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [out] cell_face            Cell->face connectivity
 * \param [out] parent_num           Parent numbering or NULL
 * \param [out] parent_entity_g_num  Parent global ids or NULL
 * \param [in]  ownership            Data ownership
 *
 */
//---> PDM_Mesh_nodal_block_poly3d_get
void
PDM_part_mesh_nodal_section_poly3d_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                    *n_face,
      PDM_g_num_t           **face_ln_to_gn,
      int                   **face_vtx_idx,
      int                   **face_vtx,
      PDM_g_num_t           **numabs,
      int                   **cell_face_idx,
      int                   **cell_face,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      PDM_ownership_t         ownership
);

/**
 * \brief Get the cell-vertex connectivity of a polyhedron section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [out] cellvtx_idx  Index of cell vertex connectivity
 * \param [out] cellvtx      Cell vertex connectivity
 * \param [in]  ownership    Data ownership
 *
 */

void
PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **cellvtx_idx,
      int                   **cellvtx,
      PDM_ownership_t         ownership
);

/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return      NULL
 *
 */
//---> PDM_Mesh_nodal_reset
void
PDM_part_mesh_nodal_reset
(
 PDM_part_mesh_nodal_t *pmn
);


/**
 * \brief  Compute a global numbering in a section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 * \param [in]  ownership    Ownership
 *
 */
//---> PDM_Mesh_nodal_g_num_in_block_compute

void
PDM_part_mesh_nodal_g_num_in_section_compute
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const PDM_ownership_t         ownership
);

/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Return number elements of a partition
 *
 */
// ---> PDM_Mesh_nodal_n_cell_get

int
PDM_part_mesh_nodal_n_elmts_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part
);

/**
 * \brief Get the element global numbering taking into account parent_num
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership Data ownership
 *
 * \return  Global ids of element in current partition
 *
 */
// ---> PDM_Mesh_nodal_g_num_get_from_part

PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get_from_part
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part,
      PDM_ownership_t         ownership
);

/**
 * \brief Free partially a part_mesh_nodal structure
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return      NULL
 *
 */
//---> PDM_Mesh_nodal_partial_free

void
PDM_part_mesh_nodal_partial_free
(
 PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return true if the vertices are defined from parents
 */
//---> PDM_Mesh_nodal_is_set_coord_from_parent

int
PDM_part_mesh_nodal_is_set_coord_from_parent
(
 PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief Get global element numbering of block elements inside the block
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section    Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [in]  ownership    Data ownership
 *
 * \return      Return global numbering of block elements inside the block
 *
 */
//---> PDM_Mesh_nodal_block_g_num_get

PDM_g_num_t *
PDM_part_mesh_nodal_section_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
);

/**
 * \brief  Return parent element number to local number
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent element number to local number
 *
 */
//---> PDM_Mesh_nodal_num_cell_parent_to_local_get

int *
PDM_part_mesh_nodal_num_elmt_parent_to_local_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part
);


/**
 * \brief  Return parent num of vertices
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent of vertices
 *
 */
//---> PDM_Mesh_nodal_vertices_parent_get
const int *
PDM_part_mesh_nodal_vertices_parent_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part
 );

/**
 * \brief  Return parent  absolute number
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent of vertices
 *
 */
//---> PDM_Mesh_nodal_vertices_g_num_parent_get
const PDM_g_num_t *
PDM_part_mesh_nodal_vertices_g_num_parent_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part
 );

/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_cell         Number of polyhedra
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
//---> PDM_Mesh_nodal_cell3d_cellface_add
void
PDM_part_mesh_nodal_cell3d_cellface_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_cell,
const int                     n_face,
const int                    *face_vtx_idx,
const int                    *face_vtx,
const PDM_g_num_t            *face_ln_to_gn,
const int                    *cell_face_idx,
const int                    *cell_face,
const PDM_g_num_t            *cell_ln_to_gn,
const PDM_ownership_t         ownership
);

/**
 * \brief  Add some 2D faces from face edge conectivity.
 *
 * For each face, this function searchs the type of the face (triangles, quandrangles, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
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
//---> PDM_Mesh_nodal_cell2d_celledge_add
void
PDM_part_mesh_nodal_face2d_faceedge_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_face,
const int                     n_edge,
const int                    *edge_vtx,
const int                    *face_edge_idx,
const int                    *face_edge,
const PDM_g_num_t            *face_ln_to_gn,
const PDM_ownership_t         ownership
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
//---> PDM_Mesh_nodal_cells_cellvtx_add
void
PDM_part_mesh_nodal_cells_cellvtx_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_cell,
const int                    *cell_vtx_idx,
const int                    *cell_vtx,
const PDM_g_num_t            *numabs,
const PDM_ownership_t         ownership
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
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_Mesh_nodal_faces_facevtx_add
void
PDM_part_mesh_nodal_faces_facevtx_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_face,
const int                    *face_vtx_idx,
const int                    *face_vtx,
const PDM_g_num_t            *numabs,
const PDM_ownership_t         ownership
);


/**
 * \brief  Return geom_kind and id (local to this geom_kind) of a section
 *
 * \param [in]  pmn                      Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  i_section                Unique section identifier
 * \param [out] geom_kind                Geometry kind (corner, ridge, surface or volume)
 * \param [out] id_section_in_geom_kind  Section identifier local to the geometry kind
 *
 */

void
PDM_part_mesh_nodal_section_id_and_geom_kind_get
(
       PDM_part_mesh_nodal_t  *pmn,
 const int                     i_section,
       PDM_geometry_kind_t    *geom_kind,
       int                    *id_section_in_geom_kind
 );

/**
 * \brief  Return unique identifier of a section
 *
 * \param [in]  pmn                      Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind                Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section_in_geom_kind  Section identifier local to the geometry kind
 *
 * \return   Unique section identifier
 *
 */

int
PDM_part_mesh_nodal_section_id_from_geom_kind_get
(
       PDM_part_mesh_nodal_t  *pmn,
 const PDM_geometry_kind_t     geom_kind,
 const int                     id_section_in_geom_kind
 );

/**
 * \brief  Return number of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return  Number of sections
 *
 */
//---> PDM_Mesh_nodal_n_blocks_get

int
PDM_part_mesh_nodal_n_section_get
(
 PDM_part_mesh_nodal_t *pmn
);

/**
 * \brief  Return ids of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return  Ids of sections
 *
 */
//---> PDM_Mesh_nodal_blocks_id_get

int *
PDM_part_mesh_nodal_sections_id_get
(
 PDM_part_mesh_nodal_t *pmn
);


//********************************************************************************//
// /**
//  * \brief Create a new Mesh nodal from elements selected in a parent Mesh nodal
//  *
//  * \param [in]   parent_mesh       Parent Mesh nodal structure
//  * \param [in]   n_select_elt      Number of selected element for each partition of each nodal block
//  * \param [in]   select_elt_l_num  Local numbers of selected elements (for each partition of each nodal block)
//  *
//  * \return       Pointer to new \ref PDM_Mesh_nodal object
//  *
//  */
//---> not used in PDM/CWIPI

// PDM_Mesh_nodal_t *
// PDM_Mesh_nodal_extract_selection
// (
//  PDM_Mesh_nodal_t  *parent_mesh,
//  const int        **n_select_elt,
//  const int       ***select_elt_l_num
//  );

void
PDM_part_mesh_nodal_group_get
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind,
 const int                     i_part,
 const int                     i_group,
       int                    *n_group_elmt,
       int                   **group_elmt,
       PDM_g_num_t           **group_ln_to_gn,
       PDM_ownership_t         ownership
);

int*
PDM_part_mesh_nodal_compute_sections_idx
(
 PDM_part_mesh_nodal_t  *pmn,
 PDM_geometry_kind_t     geom_kind,
 const int               id_part
);

int
PDM_part_mesh_nodal_n_group_get
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind
);


PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_part_mesh_nodal_elmts_get
(
 PDM_part_mesh_nodal_t  *pmn,
 PDM_geometry_kind_t     geom_kind
);


/**
 * \brief Return the geometry kind of highest dimension
 * for a given \ref PDM_part_mesh_nodal_t object
 *
 * \param [in] pmn    Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return Geometry kind of highest dimension
 *
 */

PDM_geometry_kind_t
PDM_part_mesh_nodal_principal_geom_kind_get
(
 PDM_part_mesh_nodal_t  *pmn
 );


/**
 * \brief Return the cell->vertex connectivity
 * The output pointers are owned by the user.
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind     Geometry kind (corner, ridge, surface or volume)
 * \param [in]  i_part        Partition identifier
 * \param [out] cell_vtx_idx  Index for the cell->vertex connectivity
 * \param [out] cell_vtx      Cell->vertex connectivity
 *
 * \return Number of cells in current partition
 *
 */

int
PDM_part_mesh_nodal_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_t  *pmn,
        PDM_geometry_kind_t     geom_kind,
  const int                     i_part,
        int                   **cell_vtx_idx,
        int                   **cell_vtx
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_H__ */
