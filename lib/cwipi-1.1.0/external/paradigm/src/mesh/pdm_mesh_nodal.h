/*
 * \file
 */

#ifndef __PDM_MESH_NODAL_H__
#define __PDM_MESH_NODAL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
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

typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

} PDM_block_id_block_t;

/*----------------------------------------------------------------------------
 * Geometric element type
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_MESH_NODAL_POINT,
  PDM_MESH_NODAL_BAR2,
  PDM_MESH_NODAL_TRIA3,
  PDM_MESH_NODAL_QUAD4,
  PDM_MESH_NODAL_POLY_2D,
  PDM_MESH_NODAL_TETRA4,
  PDM_MESH_NODAL_PYRAMID5,
  PDM_MESH_NODAL_PRISM6,
  PDM_MESH_NODAL_HEXA8,
  PDM_MESH_NODAL_POLY_3D,
  PDM_MESH_NODAL_BARHO,
  PDM_MESH_NODAL_TRIAHO,
  PDM_MESH_NODAL_QUADHO,
  PDM_MESH_NODAL_TETRAHO,
  PDM_MESH_NODAL_PYRAMIDHO,
  PDM_MESH_NODAL_PRISMHO,
  PDM_MESH_NODAL_HEXAHO,
  PDM_MESH_NODAL_BARHO_BEZIER, // temporary add to visualize Bezier curves
  PDM_MESH_NODAL_TRIAHO_BEZIER, // temporary add to visualize Bezier triangles
  PDM_MESH_NODAL_N_ELEMENT_TYPES
} PDM_Mesh_nodal_elt_t;


typedef struct _PDM_Mesh_nodal_t PDM_Mesh_nodal_t;
typedef struct PDM_Mesh_nodal_vtx_t PDM_Mesh_nodal_vtx_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 *
 * \brief Get the dimension of an element
 *
 * \param[in]  type    Element type
 *
 * \return     Dimension of the element
 *
 */

int
PDM_Mesh_nodal_elt_dim_get
(
 PDM_Mesh_nodal_elt_t type
 );

int
PDM_Mesh_nodal_is_2D_element
(
  PDM_Mesh_nodal_elt_t type
);

int
PDM_Mesh_nodal_is_3D_element
(
  PDM_Mesh_nodal_elt_t type
);


/**
 * \brief Get the number of vertices of an element type
 *
 * \param [in]   type     Element type
 * \param [in]   comm     Element order
 *
 * \return       Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vtx_elt_get
(
  PDM_Mesh_nodal_elt_t type,
  const int            order
);

/**
 * \brief Get if the element is consider as HO or not
 *
 * \param [in]   type     Element type
 *
 * \return       0 if element is not ho else 1
 *
 */
int
PDM_Mesh_nodal_elmt_is_ho
(
  PDM_Mesh_nodal_elt_t type
);


/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       Pointer to \ref PDM_Mesh_nodal object
 *
 */
//---> PDM_part_mesh_nodal_create
PDM_Mesh_nodal_t*
PDM_Mesh_nodal_create
(
const int     n_part,
const PDM_MPI_Comm comm
);

/**
 * \brief Free partially a nodal mesh structure
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return      NULL
 *
 */
//---> PDM_part_mesh_nodal_partial_free
void
PDM_Mesh_nodal_partial_free
(
PDM_Mesh_nodal_t *mesh
);

/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return      NULL
 *
 */
//---> PDM_part_mesh_nodal_free
void
PDM_Mesh_nodal_free
(
PDM_Mesh_nodal_t *mesh
);

/**
 * \brief Define partition vertices
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 * \param [in]  ownership Vertices ownship
 *
 */
//---> PDM_part_mesh_nodal_coord_set
void
PDM_Mesh_nodal_coord_set
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part,
 const int               n_vtx,
 const PDM_real_t       *coords,
 const PDM_g_num_t      *numabs,
 PDM_ownership_t         ownership
);


/**
 * \brief  Return number of vertices
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */
//---> PDM_part_mesh_nodal_n_vtx_get
int
PDM_Mesh_nodal_n_vertices_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
);


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Coordinates of vertices
 *
 */
//---> PDM_part_mesh_nodal_vtx_coord_get
const double *
PDM_Mesh_nodal_vertices_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
);

/**
 * \brief  Return cell centers
 *
 * \param [in]  mesh       Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part    Partition identifier
 * \param [in]  id_block   Block identifier
 *
 * \return  Return cell centers
 *
 */
//---> PDM_part_mesh_nodal_section_elt_center_get
const double *
PDM_Mesh_cell_centers_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_block,
 const int               id_part
);


/**
 * \brief  Return parent num of vertices
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return  Parent of vertices
 *
 */
//---> PDM_part_mesh_nodal_vertices_parent_get
const int *
PDM_Mesh_nodal_vertices_parent_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
 );


/**
 * \brief  Return global numbering of vertices
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Global numbering of vertices
 *
 */
//---> PDM_part_mesh_nodal_vtx_g_num_get
const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
);


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return true if the vertices are defined from parents
 */
//---> PDM_part_mesh_nodal_is_set_coord_from_parent
int
PDM_Mesh_nodal_is_set_coord_from_parent
(
 PDM_Mesh_nodal_t *mesh
);


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  n_vtx_parent   Number of parent vertices
 * \param [in]  numabs         Global numbering (size = \ref n_vtx)
 * \param [in]  num_parent     Numbering in the parent numbering (size = \ref n_vtx)
 * \param [in]  coords_parent  Parent interlaced coordinates (size = 3 * \ref n_vtx_parent)
 * \param [in]  numabs_parent  Parent global numbering (size = \ref n_vtx_parent)
 *
 */
//---> PDM_part_mesh_nodal_coord_from_parent_set
void
PDM_Mesh_nodal_coord_from_parent_set
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part,
 const int               n_vtx,
 const int               n_vtx_parent,
 const PDM_g_num_t      *numabs,
 const int              *num_parent,
 const PDM_real_t       *coords_parent,
 const PDM_g_num_t      *numabs_parent,
 const PDM_ownership_t  ownership
);


/**
 * \brief  Return number of blocks
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return  Number of blocks
 *
 */
//---> PDM_part_mesh_nodal_n_section_get
int
PDM_Mesh_nodal_n_blocks_get
(
 PDM_Mesh_nodal_t *mesh
);


/**
 * \brief  Return blocks identifier
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return  Blocks identifier
 *
 */
//---> PDM_part_mesh_nodal_sections_id_in_geom_kind_get
int *
PDM_Mesh_nodal_blocks_id_get
(
 PDM_Mesh_nodal_t *mesh
);


/**
 * \brief  Return number of partitions
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return  Number of partitions
 *
 */
//---> PDM_part_mesh_nodal_n_part_get
int
PDM_Mesh_nodal_n_part_get
(
 PDM_Mesh_nodal_t *mesh
);


/**
 * \brief  Return type of block
 *
 * \param [in]  mesh       Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block   Block identifier
 *
 * \return  Type of block
 *
 */
//---> PDM_part_mesh_nodal_section_elt_type_get
PDM_Mesh_nodal_elt_t
PDM_Mesh_nodal_block_type_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block
);


/**
 * \brief  Add a new block to the current mesh
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  ownership Vertices ownship
 *
 * \return Block identifier
 *
 */
//---> PDM_part_mesh_nodal_section_add
int
PDM_Mesh_nodal_block_add
(
      PDM_Mesh_nodal_t     *mesh,
const PDM_Mesh_nodal_elt_t  t_elt,
const PDM_ownership_t       ownership
);


/**
 * \brief Define a standard block
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
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */
//---> PDM_part_mesh_nodal_section_std_set
void
PDM_Mesh_nodal_block_std_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const int               n_elt,
const PDM_l_num_t      *connec,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
);


/**
 * \brief Return standard block description
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
 * \param [in]   mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]   id_block       Block identifier
 * \param [in]   id_part        Partition identifier
 * \param [out]  connect        Connectivity
 *
 */
//---> PDM_part_mesh_nodal_section_std_get
void
PDM_Mesh_nodal_block_std_get
(
      PDM_Mesh_nodal_t  *mesh,
const int                id_block,
const int                id_part,
      PDM_l_num_t      **connec
);


/**
 * \brief Get number of block elements
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *
 */
//---> PDM_part_mesh_nodal_section_n_elt_get
int
PDM_Mesh_nodal_block_n_elt_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
);

/**
 * \brief Get global element numbering of block elements inside the block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global numbering of block elements inside the block
 *
 */
//---> PDM_part_mesh_nodal_section_g_num_get
PDM_g_num_t *
PDM_Mesh_nodal_block_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
);

/**
 * \brief Get global element numbering of block elements
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global element numbering of block elements
 *
 */
//---> PDM_part_mesh_nodal_g_num_get
PDM_g_num_t *
PDM_Mesh_nodal_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
);


/**
 * \brief Compute cell centers of a part of a block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_section_elt_center_compute
void
PDM_Mesh_nodal_cell_centers_compute
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               i_part,
const PDM_ownership_t  ownership
);


/**
 * \brief Get parent numbering of block elements
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return parent numbering of block elements
 *
 */
//---> PDM_part_mesh_nodal_section_parent_num_get
int *
PDM_Mesh_nodal_block_parent_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
);


/**
 * \brief Define a polygon block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */
//---> PDM_part_mesh_nodal_section_poly2d_set
void
PDM_Mesh_nodal_block_poly2d_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const PDM_l_num_t       n_elt,
const PDM_l_num_t      *connec_idx,
const PDM_l_num_t      *connec,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
);


/**
 * \brief Return a polygon block description
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] n_elt          Number of elements
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [out] numabs         Global numbering in the mesh
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
 *
 */
//---> PDM_part_mesh_nodal_section_poly2d_get
void
PDM_Mesh_nodal_block_poly2d_get
(
      PDM_Mesh_nodal_t  *mesh,
 const int               id_block,
 const int               id_part,
       PDM_l_num_t     **connec_idx,
       PDM_l_num_t     **connec
);



/**
 * \brief Get global inside numbering of block elements
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global inside numbering of block elements
 *
 */

// IDENTIQUE à PDM_Mesh_nodal_block_g_num_get

PDM_g_num_t *
PDM_Mesh_nodal_block_inside_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
);



/**
 * \brief Define a polyhedra block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  face_ln_to_gn  Face global numbering
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */
//---> PDM_part_mesh_nodal_section_poly3d_set
void
PDM_Mesh_nodal_block_poly3d_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const PDM_l_num_t       n_elt,
const PDM_l_num_t       n_face,
const PDM_l_num_t      *facvtx_idx,
const PDM_l_num_t      *facvtx,
const PDM_g_num_t      *face_ln_to_gn,
const PDM_l_num_t      *cellfac_idx,
const PDM_l_num_t      *cellfac,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
);

/**
 * \brief Define a polyhedra block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] n_elt          Number of polyhedra
 * \param [out] n_face         Number of faces used to describe polyhedra
 * \param [out] facvtx_idx     Index of face vertex connectivity
 * \param [out] facvtx         Face vertex connectivity
 * \param [out] cellfac_idx    Index of cell face connectivity
 * \param [out] cellfac        Cell face connectivity
 * \param [out] numabs         Global numbering
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
 *
 */
//---> PDM_part_mesh_nodal_section_poly3d_get
void
PDM_Mesh_nodal_block_poly3d_get
(
      PDM_Mesh_nodal_t  *mesh,
const int                id_block,
const int                id_part,
      PDM_l_num_t       *n_face,
      PDM_l_num_t      **facvtx_idx,
      PDM_l_num_t      **facvtx,
      PDM_l_num_t      **cellfac_idx,
      PDM_l_num_t      **cellfac
);

/**
 * \brief Get the cell-vertex connectivity of a polyhedra block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] cellvtx_idx    Index of cell vertex connectivity
 * \param [out] cellvtx        Cell vertex connectivity
 *
 */
//---> PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
void
PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get
(
      PDM_Mesh_nodal_t  *mesh,
 const int               id_block,
 const int               id_part,
 PDM_l_num_t           **cellvtx_idx,
 PDM_l_num_t           **cellvtx
);


/**
 *
 * \brief Get cell-vertex connectivity
 *
 * \param [in]   mesh           Pointer to \ref PDM_mesh_nodal object
 * \param [in]   id_part        Partition identifier
 * \param [out]  cellvtx_idx    Index of cell vertex connectivity
 * \param [out]  cellvtx        Cell vertex connectivity
 *
 */

void
PDM_Mesh_nodal_cell_vtx_connectivity_get
(
 PDM_Mesh_nodal_t       *mesh,
 const int               id_part,
 PDM_l_num_t           **cellvtx_idx,
 PDM_l_num_t           **cellvtx
);


/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each face
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  face_ln_to_gn  Face global numbering
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face_nb   Number of faces for each cell
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_cell3d_cellface_add
void
PDM_Mesh_nodal_cell3d_cellface_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_cell,
const int               n_face,
const PDM_l_num_t      *face_vtx_idx,
const PDM_l_num_t      *face_vtx_nb,
const PDM_l_num_t      *face_vtx,
const PDM_g_num_t      *face_ln_to_gn,
const PDM_l_num_t      *cell_face_idx,
const PDM_l_num_t      *cell_face_nb,
const PDM_l_num_t      *cell_face,
const PDM_g_num_t      *numabs,
const PDM_ownership_t  ownership
);

/**
 * \brief  Add some standard 3D cells from cell vertex conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_cell         Number of cells
 * \param [in]  cell_face_idx  Index of cell vertex connectivity
 * \param [in]  cell_face_nb   Number of vertices for each cell
 * \param [in]  cell_face      Cell vertex connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_cells_cellvtx_add
void
PDM_Mesh_nodal_cells_cellvtx_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_cell,
const PDM_l_num_t      *cell_vtx_idx,
const PDM_l_num_t      *cell_vtx_nb,
const PDM_l_num_t      *cell_vtx,
const PDM_g_num_t      *numabs,
const PDM_ownership_t  ownership
);

/**
 * \brief  Add some 2D cells from cell edge conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_edge         Number of edges used to describe polyhedra
 * \param [in]  edge_vtx_idx   Index of edge vertex connectivity
 * \param [in]  edge_vtx_nb    Number of vertices for each edge
 * \param [in]  edge_vtx       Edge vertex connectivity
 * \param [in]  cell_edge_idx  Index of cell edge connectivity
 * \param [in]  cell_edge_nb   Number of edges for each cell
 * \param [in]  cell_edge      Cell edge connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_face2d_faceedge_add
void
PDM_Mesh_nodal_cell2d_celledge_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_elt,
const int               n_edge,
const PDM_l_num_t      *edge_vtx_idx,
const PDM_l_num_t      *edge_vtx_nb,
const PDM_l_num_t      *edge_vtx,
const PDM_l_num_t      *cell_edge_idx,
const PDM_l_num_t      *cell_edge_nb,
const PDM_l_num_t      *cell_edge,
const PDM_g_num_t      *numabs,
const PDM_ownership_t  ownership
);


/**
 * \brief  Add some 2D cells from cell vertex connectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygon
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each edge
 * \param [in]  face_vtx       Edge vertex connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_faces_facevtx_add
void
PDM_Mesh_nodal_faces_facevtx_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_face,
const PDM_l_num_t      *face_vtx_idx,
const PDM_l_num_t      *face_vtx_nb,
const PDM_l_num_t      *face_vtx,
const PDM_g_num_t      *numabs,
const PDM_ownership_t  ownership
);


/**
 * \brief  Compute a global numbering in a block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  ownership      Ownership
 *
 */
//---> PDM_part_mesh_nodal_g_num_in_section_compute
void
PDM_Mesh_nodal_g_num_in_block_compute
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const PDM_ownership_t  ownership
);


/**
 * \brief  Return parent cell number to local number
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Parent cell number to local number
 *
 */
//---> PDM_part_mesh_nodal_num_elmt_parent_to_local_get
int *
PDM_Mesh_nodal_num_cell_parent_to_local_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part
);


/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  mesh      Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Return number elements of a partition
 *
 */
//---> PDM_part_mesh_nodal_n_elmts_get
int
PDM_Mesh_nodal_n_cell_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part
);


/**
 * \brief  Return parent  absolute number
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return  Parent of vertices
 *
 */
//---> PDM_part_mesh_nodal_vertices_g_num_parent_get
const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_parent_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
);

/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return      NULL
 *
 */
//---> PDM_part_mesh_nodal_reset
void
PDM_Mesh_nodal_reset
(
 PDM_Mesh_nodal_t *mesh
);

/**
 * \brief Returns the number of vertices in an element
 *
 * \param [in]  elt_type   Element type
 * \param [in]  order      Element order
 *
 * \return      Number of vertices in element
 *
 */
// A REMPLACER PAR PDM_Mesh_nodal_n_vtx_elt_get?

int
PDM_Mesh_nodal_n_vertices_element
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
);



/**
 * \brief Compute cell extents of a part of block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  tolerance      Expansion tolerance for bounding boxes
 * \param [out] extents        Extents of mesh elements in current part of current block
 *
 */
//---> PDM_part_mesh_nodal_section_elt_extents_compute
void
PDM_Mesh_nodal_compute_cell_extents
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_block,
 const int               i_part,
 const double            tolerance,
       double           *extents
);



/**
 * \brief Get the cell global numbering
 *
 * \param [in]  mesh     Pointer to \ref PDM_Mesh_nodal object
 *
 * \return      NULL
 *
 */
//---> PDM_part_mesh_nodal_g_num_get_from_part
PDM_g_num_t *
PDM_Mesh_nodal_g_num_get_from_part
(
 PDM_Mesh_nodal_t *mesh,
 const int i_part
);


/**
 * \brief Create a new Mesh nodal from elements selected in a parent Mesh nodal
 *
 * \param [in]   parent_mesh       Parent Mesh nodal structure
 * \param [in]   n_select_elt      Number of selected element for each partition of each nodal block
 * \param [in]   select_elt_l_num  Local numbers of selected elements (for each partition of each nodal block)
 *
 * \return       Pointer to new \ref PDM_Mesh_nodal object
 *
 */
// NOT USED
PDM_Mesh_nodal_t *
PDM_Mesh_nodal_extract_selection
(
 PDM_Mesh_nodal_t  *parent_mesh,
 const int        **n_select_elt,
 const int       ***select_elt_l_num
 );



/**
 * \brief Reset cell center computation
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 */
//---> PDM_part_mesh_nodal_section_elt_center_reset
void
PDM_Mesh_nodal_cell_centers_reset
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               i_part
);

void
PDM_Mesh_nodal_ho_parent_node
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering,
       int                  *parent_node
 );


void
PDM_Mesh_nodal_reorder_elt_vtx
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const char                 *ho_ordering_in,
 const char                 *ho_ordering_out,
 const int                   n_elt,
       int                  *elt_vtx_in,
       int                  *elt_vtx_out
 );



PDM_geometry_kind_t
PDM_Mesh_nodal_geom_kind_from_elt_type
(
 PDM_Mesh_nodal_elt_t t_elt
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
