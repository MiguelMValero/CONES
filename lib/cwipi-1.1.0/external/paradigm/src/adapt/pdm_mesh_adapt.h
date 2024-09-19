/*
 * \file
 */

#ifndef PDM_MESH_ADAPT_H
#define PDM_MESH_ADAPT_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2021       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/**

This file defines a parallel and genric workflow for the mesh adaptation.
This workflow can be called in a scientfic software without stopping its execution.

Entries are  a partitioned mesh with fields defined on it. The outpouts are the
partitioned adapted mesh with interpolated fields defined on it. To ensure a good
load balancing, the paritioning of the adapted mesh does not come from the source
mesh. It is recomputed.

The workflow is divided into 9 steps :
  1. Initialization with \ref PDM_Mesh_adapt_create functions. This step defines
     the global parameters (adaptation method, adapatation tool,
     graph partitioner, interpolation method, ...)
  2. Definition of source mesh with \p PDM_Mesh_adapt_src*
  3. Definition of boundary surface geometry with \p PDM_Mesh_adapt_geom_repr*
  4. Definition of specific options of the adaptation tool with
     \ref PDM_Mesh_adapt_tool_param_set
  5. Definition of field families with \ref PDM_Mesh_adapt_field_family_add.
     A family has 3 properties : the location of degrees of freedom, the
     interpolation method and the storage of components
  6. Computation of the adapted mesh and its partitioning with
     \ref PDM_Mesh_adapt_geom_compute.
     This function also computes the geometric algorithms necessary for the field
     interpolation. This function groups together all the blocking MPI calls which
     makes it possible to transfer the interpolated fields without blocking MPI
     exchange
  7. Getting of the target mesh properties with \p PDM_Mesh_adapt_tgt*
  8. Send source fields with \ref PDM_Mesh_adapt_src_field_issend
  9. Recevive target fields with \ref PDM_Mesh_adapt_tgt_field_irecv

This describes a simple workflow. it is possible to particularize a workflow
using some callbacks to take into the numerical methods of softwares. User can
define some user locations of the dregrees of freedom and he can define some user
interpolation method from 3 type of geometric results :
  1.  k-nearest neighbors into the source mesh (parallel algorithm)
  2.  Location of target of degrees of freedon into the source mesh
  3.  Intersection of source and target meshes (parallel algorithm)

*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Enum definitions
 *============================================================================*/

/**
 * \enum  PDM_Mesh_adapt_method_t
 * \brief Mesh adaptation method
 *
 */

typedef enum {

  PDM_MESH_ADAPT_REMESHING,   /*!< Remeshing */
  PDM_MESH_ADAPT_REFINMENT,   /*!< Refinment */
  PDM_MESH_ADAPT_METHOD_N     /*!< Number of methods */

} PDM_Mesh_adapt_method_t;


/**
 * \enum  PDM_Mesh_adapt_periodicity_t
 * \brief Periodicity
 *
 */

typedef enum {

  PDM_MESH_ADAPT_PERIO_TRANS,  /*!< Translation periodicity */
  PDM_MESH_ADAPT_PERIO_ROT,    /*!< Rotation periodicity */
  PDM_MESH_ADAPT_PERIO_N       /*!< Number of type of periodicity */

} PDM_Mesh_adapt_perio_t;


/**
 * \enum  PDM_Mesh_adapt_method_t
 * \brief Mesh adaptation method
 *
 */

typedef enum {

  PDM_MESH_ADAPT_ADAPTCELLS, /*!< Nuga - Adapt cells */
  PDM_MESH_ADAPT_FEFLO,      /*!< Feflo.a */
  PDM_MESH_ADAPT_PARMMG,     /*!< PARMMG */
  PDM_MESH_ADAPT_TREEPART,   /*!< Treepart */
  PDM_MESH_ADAPT_TOOL_N      /*!< Number of tools */

} PDM_Mesh_adapt_tool_t;

/**
 * \enum  PDM_Mesh_adapt_geom_repr_t
 * \brief Geometric representation of the boundary
 *
 */

typedef enum {

  PDM_MESH_ADAPT_GEOM_REPR_FROM_SRC_MESH,  /*!< From Boundary source mesh (default) */
  PDM_MESH_ADAPT_GEOM_REPR_DEDICATED_MESH, /*!< Dedicated */
  PDM_MESH_ADAPT_GEOM_REPR_STL,            /*!< STL file */
  PDM_MESH_ADAPT_GEOM_REPR_IGES,           /*!< IGES file  */
  PDM_MESH_ADAPT_GEOM_REPR_STEP,           /*!< STEP file */
  PDM_MESH_ADAPT_GEOM_REPR_EGADS_LITE,     /*!< Egads lite file */
  PDM_MESH_ADAPT_GEOM_REPR_N               /*!< Number of type of geometric
                                                representations */
} PDM_Mesh_adapt_geom_repr_t;


/**
 * \enum  PDM_Mesh_adapt_method_t
 * \brief Geometric adaptation criterion
 *
 */

typedef enum {

  PDM_MESH_ADAPT_METRIC_VTX,            /*!< Metric defined to vertices */
  PDM_MESH_ADAPT_METRIC_CELL,           /*!< Metric defined to elements (not recommended) */
  PDM_MESH_ADAPT_SUBDIVSION_LEVEL_VTX,  /*!< Subdvision level defined to vertices */
  PDM_MESH_ADAPT_SUBDIVSION_LEVEL_CELL, /*!< Subdvision level defined to elements */
  PDM_MESH_ADAPT_CRITERION_N            /*!< Number of cirteria */

} PDM_Mesh_adapt_criterion_t;


/**
 * \enum  PDM_Mesh_adapt_part_t
 * \brief Partition method for the target mesh
 *
 */

typedef enum {

  PDM_MESH_ADAPT_PART_PTSCOTCH, /*!< PT-Scotch */
  PDM_MESH_ADAPT_PART_PARMETIS, /*!< ParMETIS */
  PDM_MESH_ADAPT_PART_TREEPART, /*!< TreePart */
  PDM_MESH_ADAPT_PART_HILBERT,  /*!< Hilbert space filling curve*/
  PDM_MESH_ADAPT_PART_MORTON,   /*!< Morton space filling curve */
  PDM_MESH_ADAPT_PART_N         /*!< Number of partition method */

} PDM_Mesh_adapt_part_tool_t;


/**
 * \enum  PDM_Mesh_adapt_inter_part_graph_t
 * \brief MPI graph communication between partitions
 *
 */

typedef enum {

  PDM_MESH_ADAPT_INTER_PART_GRAPH_NONE,            /*!< No communication */
  PDM_MESH_ADAPT_INTER_PART_GRAPH_ELT_THROUGH_FACE,/*!< Communication between
                                                        elements through faces */
  PDM_MESH_ADAPT_INTER_PART_GRAPH_ELT_THROUGH_NODE,/*!< Communication between
                                                        elements through nodes */
  PDM_MESH_ADAPT_INTER_PART_GRAPH_FACE_THROUGH_FACE,/*!< Communication between
                                                        faces through faces */
  PDM_MESH_ADAPT_INTER_PART_GRAPH_NODE_THROUGH_NODE,/*!< Communication  between
                                                         nodes through nodes*/
  PDM_MESH_ADAPT_INTER_PART_GRAPH_N

} PDM_Mesh_adapt_graph_t;


/**
 * \enum  PDM_Mesh_adapt_dof_location_t
 * \brief Degrees of freedom location
 *
 */

typedef enum {

  PDM_MESH_ADAPT_DOF_LOC_NODE, /*!< degrees of freedom defined at node  */
  PDM_MESH_ADAPT_DOF_LOC_ELTS, /*!< degrees of freedom defined at element */
  PDM_MESH_ADAPT_DOF_LOC_USER, /*!< degrees of freedom defined by user */
  PDM_MESH_ADAPT_DOF_LOC_N     /*!< Number of  degrres of freedom location  */

} PDM_Mesh_adapt_dof_location_t;


/**
 * \enum  PDM_Mesh_adapt_field_interp_t
 * \brief Interpolation methods
 *
 */

typedef enum {

  PDM_MESH_ADAPT_INTERP_DEFAULT_FROM_KNN,          /*!< Default interpolation
                                                        from the 'k' closest vertices
                                                        of the source mesh */
  PDM_MESH_ADAPT_INTERP_DEFAULT_FROM_LOCATION,     /*!< Default interpolation
                                                        from the location of
                                                        target degrees of freedom
                                                        into the source mesh  */
  PDM_MESH_ADAPT_INTERP_DEFAULT_FROM_INTERSECTION, /*!< Default interpolation from
                                                        the intersection between
                                                        source and target meshes  */
  PDM_MESH_ADAPT_INTERP_USER_FROM_KNN,             /*!< user interpolation (callback)
                                                        from the 'k' closest vertices
                                                        of the source mesh */
  PDM_MESH_ADAPT_INTERP_USER_FROM_LOCATION,        /*!< User interpolation (callback)
                                                        from the location of
                                                        target degrees of freedom
                                                        into the source mesh  */
  PDM_MESH_ADAPT_INTERP_USER_FROM_INTERSECTION,    /*!< Field defined at node  */
  PDM_MESH_ADAPT_INTERP_INTERPOL,                  /*!< Tool : Interpol */
  PDM_MESH_ADAPT_INTERP_TREEPART,                  /*!< Tool : Treepart */
  PDM_MESH_ADAPT_INTERP_ADAPTCELLS,                /*!< Tool : Adaptcells */
  PDM_MESH_ADAPT_INTERP_N                          /*!< Number of interpolation
                                                        methods */

} PDM_Mesh_adapt_field_interp_t;

/**
 * \enum PDM_Mesh_adapt_group_entity_t
 * \brief Group entities
 *
 */

typedef enum {
  PDM_MESH_ADAPT_GROUP_CELL,  /*!< Cell group */
  PDM_MESH_ADAPT_GROUP_FACE,  /*!< Face group */
  PDM_MESH_ADAPT_GROUP_EDGE,  /*!< Edge group */
  PDM_MESH_ADAPT_GROUP_NODE,  /*!< Node group */
  PDM_MESH_ADAPT_GROUP_N      /*!< Number of group types */
} PDM_Mesh_adapt_group_entity_t;

/*============================================================================
 * Callback function prototypes
 *============================================================================*/

/**
 * \typedef void (*CWP_Interp_from_location_t)
 * \brief User interpolation function interface from location into a mesh.
 *
 * void (*CWP_Interp_from_location_t) defines the user interpolation
 * interface to take into account an user interpolation from location of target
 * points into the source mesh. Use \ref CWP_Interp_from_location_set to activate
 * the function.
 *
 * \param [in]  user_struct_about_src       Generic pointer to user structure about target
 *                                          setted from
 *                                          \ref PDM_adapt_mesh_user_src_struct_set function
 *                                          or NULL
 * \param [in]  n_tgt_vtx                   Local number of target vertices
 * \param [in]  tgt_pts_coords              Target points coordinates
 *                                          (size = 3 * n_tgt_pts)
 * \param [in]  tgt_pts_location            target points location
 *                                          (size = n_tgt_pts)
 * \param [in]  tgt_pts_dist                target points distance to location element
 *                                          (size = n_tgt_pts)
 * \param [in]  tgt_pts_bary_coords_idx     Index of Barycentric coordinates target points
 *                                          in location element
 *                                          (tgt_pts_bary_coords_idx[0] = 0 and
 *                                          size = n_tgt_pts + 1)
 * \param [in]  tgt_pts_bary_coords         Barycentric coordinates target points
 *                                          in location element
 *                                          (size = tgt_pts_bary_coords_idx[n_tgt_pts])
 * \param [in]  dof_location                Degrees of freedom location
 * \param [in]  field_stride                Field stride
 * \param [in]  src_field                   Source field
 *                                          (size depends on field type and stride)
 * \param [out] tgt_field                   Target field
 *                                          (size = stride * n_tgt_pts)
 *
 */

typedef void (*PDM_Mesh_adapt_interp_from_location_t)
(
 // TODO: Add some properties about intersected mesh. To be continued
 void                                *user_struct_about_src,
 const int                            n_tgt_pts,
 const double                         tgt_pts_coords[],
 const int                            tgt_pts_location[],
 const double                         tgt_pts_dist[],
 const int                            tgt_pts_bary_coords_idx[],
 const double                         tgt_pts_bary_coords[],
 const PDM_Mesh_adapt_dof_location_t  dof_location,
 const int                            field_stride,
 const double                         src_field[],
 double                               tgt_field[]
);


/**
 * \brief User interpolation function interface from intersection between meshes. <b>(Not implemented yet)</b>
 *
 * void (*PDM_Mesh_adapt_interp_from_intersect_t) defines the user interpolation
 * interface to take into account an user interpolation from intersection
 * between source and target meshes. Use \ref CWP_Interp_from_intersect_set to activate
 * the function.
 *
 * TODO: Add some properties about intersected mesh. To be continued
 *
 *  \param [in]  user_struct_about_src    Generic pointer to user structure about target
 *                                        setted from
 *                                        \ref PDM_adapt_mesh_user_src_struct_set function
 *                                        or NULL
 *  \param [in]  dof_location             Degrees of freedom location
 *  \param [in]  field_stride             Field stride
 *  \param [in]  src_field                Source field
 *  \param [out] tgt_field                Target field
 *
 */

typedef void (*PDM_Mesh_adapt_interp_from_intersect_t)
(
 // TODO: Add some properties about intersected mesh. To be continued
 void                                *user_struct_about_src,
 const PDM_Mesh_adapt_dof_location_t  dof_location,
 const int                            field_stride,
 const double                         src_field[],
 double                               tgt_field[]
);


/**
 * \brief User interpolation function from closest points. <b>(Not implemented yet)</b>
 *
 * void (*PDM_mesh_adapt_interp_from_closest_pts_t) defines the user interpolation
 * interface to take into account an user interpolation from <i>k</i> closest
 * points. For this method the source field is sent to the target
 * before the interpolation. Interpolation is done on the target side.
 * To transfer more informations about source, you can stored them in additional
 * components in the source fields.
 *
 * TODO: Add some properties about intersected mesh. To be continued
 *  \param [in]  user_struct_about_tgt    Generic pointer to an user structure about target
 *                                        setted from \ref PDM_adapta_mesh_user_tgt_struct_set
 *                                        function or NULL
 *  \param [in]  k_closest                value of 'k'
 *  \param [in]  n_tgt_vtx                Local number of target vertices
 *  \param [in]  tgt_vtx_coords           Coordinates of the target vertices
 *  \param [in]  tgt_src_closest_g_num    For each target vertex the global numbers of
 *                                        the 'k' closest source vertices
 *  \param [in]  tgt_src_closest_dist     For each target vertex the distance to
 *                                        the 'k' closest source vertices
 *  \param [in]  tgt_src_closest_coords   For each target vertex the coordinates of
 *                                        the 'k' closest source vertices
 *  \param [in]  dof_location             Degrees of freedom location
 *  \param [in]  field_stride             Field stride
 *  \param [in]  src_field_at_closest_vtx For each target vertex the field defined at
 *                                        the 'k' closest source vertices
 *  \param [out] tgt_field                Target field
 *
 */

typedef void (*PDM_Mesh_adapt_interp_from_closest_pts_t)
(
 // TODO: Add some properties about intersected mesh. To be continued
 void                                *user_struct_about_tgt,
 const int                           k_closest,
 const int                           n_tgt_vtx,
 const double                        tgt_vtx_coords[],
 const PDM_g_num_t                   tgt_src_closest_g_num[],
 const double                        tgt_src_closest_dist[],
 const double                        tgt_src_closest_coords[],
 const PDM_Mesh_adapt_dof_location_t dof_location,
 const int                           field_stride,
 const double                        src_field_at_closest_vtx[],
       double                        tgt_field[]
);


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct PDM_Mesh_adapt_t
 * \brief  Mesh adaptation workflow
 *
 */

typedef struct _PDM_Mesh_adapt_t PDM_Mesh_adapt_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 *
 * General functions
 *   - PDM_Mesh_adapt_create
 *   - PDM_Mesh_adapt_free
 *
 * The intial mesh is called : "source mesh"
 * The adapted mesh is called : "target mesh"
 *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a parallel mesh adaptation workflow
 *
 * \param [in] comm   PDM_MPI communicator
 * \param [in] method Mesh adaptation method
 *                       - PDM_MESH_ADAPT_REMESHING
 *                       - PDM_MESH_ADAPT_REFINMENT
 * \param [in] tool   Mesh adaptation tool
 *                       - PDM_MESH_ADAPT_ADAPTCELLS
 *                       - PDM_MESH_ADAPT_FEFLO
 *                       - PDM_MESH_ADAPT_PARMMG
 *                       - PDM_MESH_ADAPT_TREEADAPT
 * \param [in] criterion   Geometric adaptation criterion
 *                       - PDM_MESH_ADAPT_METRIC
 *                       - PDM_MESH_ADAPT_SUBDIVSION_LEVEL
 * \param [in] geom_repr  Geometric representation of the boundary
 *                       - PDM_MESH_ADAPT_GEOM_REPR_BOUNDARY_SRC_MESH
 *                       - PDM_MESH_ADAPT_GEOM_REPR_DEDICATED_MESH
 *                       - PDM_MESH_ADAPT_GEOM_REPR_STL
 *                       - PDM_MESH_ADAPT_GEOM_REPR_IGES
 *                       - PDM_MESH_ADAPT_GEOM_REPR_STEP
 * \param [in] part_tool Target mesh partitionning
 *                       - PDM_MESH_ADAPT_PART_PTSCOTCH
 *                       - PDM_MESH_ADAPT_PART_PARMETIS
 *                       - PDM_MESH_ADAPT_PART_TREEPART
 *                       - PDM_MESH_ADAPT_PART_HILBERT
 *                       - PDM_MESH_ADAPT_PART_MORTON
 * \param [in] mesh_order  Mesh order
 * \param [in] n_dom  Number of meshes
 * \param [in] n_part Number of local mesh partition for each domain (size = \p n_dom)
 *                    (same value for source and target meshes)
 * \param [in] owner  Results owner
 *                       - PDM_OWNERSHIP_KEEP                 : Paradigm will
 *                       - PDM_OWNERSHIP_USER                 : Ownership is given
 *                       - PDM_OWNERSHIP_UNGET_RESULT_IS_FREE : Free all memory that
 *                                                              not be getted by user
 *
 * \return     Identifier
 *
 */

PDM_Mesh_adapt_t *
PDM_Mesh_adapt_create
(
 const PDM_MPI_Comm                      comm,
 const PDM_Mesh_adapt_method_t           method,
 const PDM_Mesh_adapt_tool_t             tool,
 const PDM_Mesh_adapt_criterion_t        criterion,
 const PDM_Mesh_adapt_geom_repr_t        geom_repr,
 const PDM_Mesh_adapt_part_tool_t        part_tool,
 const int                               mesh_order,
 const int                               n_dom,
 const int                               n_part[],
 const PDM_ownership_t                   owner
);


/**
 *
 * \brief Free a parallel mesh adaptation workflow according to the selected
 *        ownership property
 *
 * \param [in]  ma    Mesh adaptation workflow
 *
 */

void
PDM_Mesh_adapt_free
(
 PDM_Mesh_adapt_t *ma
);


/*----------------------------------------------------------------------------*
 *
 * Functions about source mesh definition
 *   - Volume mesh    : call PDM_Mesh_adapt_src_block* to define it
 *   - Boundary mesh  : call PDM_Mesh_adapt_boundary_src_block* to define it
 *   - Entity groups (To describe boundary conditions) : call PDM_Mesh_adapt_entity_group*
 *                      to define them
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set source mesh vertices.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Current domain
 * \param [in]  i_part       Current partition
 * \param [in]  n_nodes      Number of vertices
 * \param [in]  coord        Coordinates (size = 3 * \p n_vtx)
 * \param [in]  g_num        Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_vtx_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         n_nodes,
 double            coord[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Add a connectivity block to the source mesh.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom        Current domain
 * \param [in]  block_type       Block type
 *                                 -  PDM_MESH_NODAL_POINT
 *                                 -  PDM_MESH_NODAL_BAR2
 *                                 -  PDM_MESH_NODAL_BARHO
 *                                 -  PDM_MESH_NODAL_TRIA3
 *                                 -  PDM_MESH_NODAL_TRIAHO
 *                                 -  PDM_MESH_NODAL_QUAD4
 *                                 -  PDM_MESH_NODAL_QUADHO
 *                                 -  PDM_MESH_NODAL_POLY_2D
 *                                 -  PDM_MESH_NODAL_TETRA4
 *                                 -  PDM_MESH_NODAL_TETRAHO
 *                                 -  PDM_MESH_NODAL_PYRAMID5
 *                                 -  PDM_MESH_NODAL_PYRAMIDHO
 *                                 -  PDM_MESH_NODAL_PRISM6
 *                                 -  PDM_MESH_NODAL_PRISMHO
 *                                 -  PDM_MESH_NODAL_HEXA8
 *                                 -  PDM_MESH_NODAL_HEXAHO
 *                                 -  PDM_MESH_NODAL_POLY_3D
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_src_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const int                   i_dom,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the source mesh.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
 *
 *   \code
 *             x 4
 *            /|\
 *           / | \
 *          /  |  \
 *       1 x- -|- -x 3
 *          \  |  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
 *
 *   \code
 *              5 x
 *               /|\
 *              //| \
 *             // |  \
 *          4 x/--|---x 3
 *           //   |  /
 *          //    | /
 *       1 x-------x 2
 *   \endcode
 *
 *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
 *
 *   \code
 *       4 x-------x 6
 *         |\     /|
 *         | \   / |
 *       1 x- \-/ -x 3
 *          \ 5x  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
 *
 *   \code
 *          8 x-------x 7
 *           /|      /|
 *          / |     / |
 *       5 x-------x6 |
 *         | 4x----|--x 3
 *         | /     | /
 *         |/      |/
 *       1 x-------x 2
 *   \endcode
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  l_num        Local element number in the partition (size = n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_src_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workflow
 * \param [in]  i_dom                Domain identifier
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  l_num                Local element number
 *                                   in the partition (size = n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the source mesh.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  l_num            Local element number in the partition
 *                               (size = n_elts)
 * \param [in]  g_num            Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);

/**
 *
 * \brief Set the connectivity of a polyhedron block of the source mesh.
 *
 * Connectivity is supposed to be oriented. Connectivity can be oriented by calling
 * PDM_cellface_orient
 *
 * \param [in]  ma                Mesh adaptation workflow
 * \param [in]  i_dom             Domain identifier
 * \param [in]  i_part            Partition identifier
 * \param [in]  i_block           Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [in]  face_vtx          Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  face_l_num        Face local element number in the partition
 *                                (size = n_elts)
 * \param [in]  face_g_num        Face global element number (or NULL)
 * \param [in]  cell_face_idx     Polyhedron to face index (or NULL)
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Polyhedron to face connectivity (or NULL)
 *                                The connectivity is oriented :
 *                                  - > 0 if outgoing normal,
 *                                  - < 0 otherwise
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  face_cell         Face to polyhedron connectivity (or NULL)
 *                                  - left value  : outgoing normal,
 *                                  - right value : incoming normal
 *                                (size = 2 * \p n_faces)
 * \param [in]  cell_l_num        Cell local element number in the partition
 *                                (size = n_elts)
 * \param [in]  cell_g_num        Cell global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_block_c_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 const int         n_faces,
 int               face_vtx_idx[],
 int               face_vtx[],
 int               face_l_num[],
 PDM_g_num_t       face_g_num[],
 int               cell_face_idx[],
 int               cell_face[],
 int               face_cell[],
 int               cell_l_num[],
 PDM_g_num_t       cell_g_num[]
);


/**
 * \brief Add a connectivity block to the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom             Domain identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_boundary_src_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const int                   i_dom,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the boundary source mesh.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  l_num        Local element number in the partition (size = n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_boundary_src_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the boundary source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workflow
 * \param [in]  i_dom                Domain identifier
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  l_num                Local element number in the partition (size = n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_boundary_src_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         Face to vertex connectivity
 *                               (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  l_num            Local element number in the partition
 *                               (size = \p n_elts)
 * \param [in]  g_num            Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_boundary_src_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Add an intra-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  t_graph     Graph type
 *
 * \return Graph identifier
 */

int
PDM_Mesh_adapt_src_intra_dom_graph_add
(
 PDM_Mesh_adapt_t             *ma,
 const PDM_Mesh_adapt_graph_t t_graph
);


/**
 * \brief Set an intra-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_graph     Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  n_elt_graph Number of elements in the
 * \param [in]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [in]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - Lists of 3 values :
 *                                  + Connected process
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_src_intra_dom_graph_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_graph,
 const int         i_part,
 const int         n_elt_graph,
 int               graph_idx[],
 int               graph[]
);


/**
 * \brief Set the number of entity groups
 *
 * Face groups are used to define boundary conditions.
 *
 * \param [in]  ma        Mesh adaptation workflow
 * \param [in]  g_entity  Group entity
 * \param [in]  i_dom     Domain identifier
 * \param [in]  n_group   Number of face groups
 *
 */

void
PDM_Mesh_adapt_entity_group_n_set
(
 PDM_Mesh_adapt_t              *ma,
 PDM_Mesh_adapt_group_entity_t g_entity,
 const int                     i_dom,
 const int                     n_group
);


/**
 * \brief Set a entity group
 *
 * An entity can be contained in several groups, the notion of group  is more general than
 * the notion of boudary condition
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  g_entity    Enity group
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  i_group     Group identifier
 * \param [in]  n_entities  Number of entities in the group
 * \param [in]  entities    List of entities
 * \param [in]  g_num       Global element number in the group (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_entity_group_set
(
 PDM_Mesh_adapt_t              *ma,
 PDM_Mesh_adapt_group_entity_t g_entity,
 const int                     i_dom,
 const int                     i_part,
 const int                     i_group,
 const int                     n_entities,
 int                           entities[],
 PDM_g_num_t                   g_num[]
);


/**
 * \brief Add an inter-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  t_graph     Graph type
 *
 * \return Graph identifier
 */

int
PDM_Mesh_adapt_src_inter_dom_graph_add
(
 PDM_Mesh_adapt_t             *ma,
 const PDM_Mesh_adapt_graph_t t_graph
);


/**
 * \brief Set an inter-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_graph     Graph identifier
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  n_elt_graph Number of elements in the
 * \param [in]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [in]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - Lists of 4 values :
 *                                  + Connected process
 *                                  + Connected domain
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_src_inter_dom_graph_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_graph,
 const int         i_dom,
 const int         i_part,
 const int         n_elt_graph,
 int               graph_idx[],
 int               graph[]
);


/**
 * \brief Add a periodicity relation
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  t_perio     Periodicity type
 *
 * \return Periodicity identifier
 */

int
PDM_Mesh_adapt_src_perio_add
(
 PDM_Mesh_adapt_t             *ma,
 const PDM_Mesh_adapt_perio_t t_perio
);


/**
 * \brief Set properties of a translation periodicity
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_perio     Periodicity identifier
 * \param [in]  vect        Translation vector
 *
 * \return Periodicity identifier
 */

int
PDM_Mesh_adapt_src_perio_translation_set
(
 PDM_Mesh_adapt_t  *ma,
 const int          i_perio,
 double             vect[]
);


/**
 * \brief Set properties of a rotation periodicity
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_perio     Periodicity identifier
 * \param [in]  axis        Axis
 * \param [in]  center      Center
 * \param [in]  angle       Angle
 *
 * \return Periodicity identifier
 */

int
PDM_Mesh_adapt_src_perio_rotation_set
(
 PDM_Mesh_adapt_t  *ma,
 const int          i_perio,
 double             axis[],
 double             center[],
 double             angle
);


/**
 * \brief Set a graph communication between 2 periodicity surfaces
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_perio     Periodicity identifier
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  n_elt_graph Number of elements in the
 * \param [in]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [in]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - List of 4 values :
 *                                  + Connected process
 *                                  + Connected domain
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_src_perio_graph_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_perio,
 const int         i_dom,
 const int         i_part,
 const int         n_elt_graph,
 int               graph_idx[],
 int               graph[]
);


/**
 * \brief Finalize the source mesh definition.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided (MPI collective communications)
 *
 * \param [in]  ma                Mesh adaptation workflow
 *
 */

void
PDM_Mesh_adapt_src_finalize
(
 PDM_Mesh_adapt_t *ma
);


/*----------------------------------------------------------------------------*
 *
 * Functions about the geometric representation of the boundary
 *     - Set the geometry reprsentation of the boundary
 *         * Dedicated boundary mesh :
 *             + Default : Boudary source mesh
 *                         (defined by PDM_Mesh_adapt_boundary_src*)
 *             + Another boundary mesh :
 *                         call PDM_Mesh_adapt_geom_repr* to define it
 *                         define it
 *             + Define ridges and corners : Entities that must not be destroyed
 *               by remeshing (PDM_Mesh_adapt_geom_repr_rigde_set,
 *                             PDM_Mesh_adapt_geom_repr_corners_set)
 *         * From STL file
 *         * From iges file
 *         * From step file
 *         * ...
 *
 *----------------------------------------------------------------------------*/


/**
 * \brief Add a connectivity block to the mesh that describes the geometric
 *        representation of the boundary.
 *
 *  This function is called only if \ref PDM_MESH_ADAPT_GEOM_REPR_ANOTHER_MESH
 *  is selected.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  n_dom        Number of domains
 * \param [in]  n_part       Number of local partitions for each domain
 * \param [in]  order        Mesh order
 *
 */

void
PDM_Mesh_adapt_geom_repr_mesh_init
(
 PDM_Mesh_adapt_t   *ma,
 const int           n_dom,
 const int           n_part[],
 const int           order
);


/**
 * \brief Set vertices of the geometric representation of the boundary
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part       Current partition
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coord        Coordinates (size = 3 * \p n_vtx)
 * \param [in]  g_num        Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_geom_repr_mesh_vtx_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         n_vtx,
 double            coord[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Add a connectivity block to the mesh that describes the geometric
 *        representation of the boundary.
 *
 *  This function is called only if \ref PDM_MESH_ADAPT_GEOM_REPR_ANOOTHER_MESH
 *  is selected.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  block_type       Block type
 *                                 -  PDM_MESH_NODAL_POINT
 *                                 -  PDM_MESH_NODAL_BAR2
 *                                 -  PDM_MESH_NODAL_BARHO
 *                                 -  PDM_MESH_NODAL_TRIA3
 *                                 -  PDM_MESH_NODAL_TRIAHO
 *                                 -  PDM_MESH_NODAL_QUAD4
 *                                 -  PDM_MESH_NODAL_QUADHO
 *                                 -  PDM_MESH_NODAL_POLY_2D
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_geom_repr_mesh_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const int                   i_dom,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the the mesh that describes the geometric
 *        representation of the boundary.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref PDM_MESH_NODAL_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref PDM_MESH_NODAL_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref PDM_MESH_NODAL_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  l_num        Local element number in the partition (size = n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_geom_repr_mesh_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the the mesh that describes
 *        the geometric representation of the boundary.the boundary source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workflow
 * \param [in]  i_dom                Domain identifier
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  l_num                Local element number in the partition
 *                                   (size = n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_geom_repr_mesh_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the geometric representation
 *        of the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         Connectitivty face to vertex (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  l_num            Local element number in the partition (size = \p n_elts)
 * \param [in]  g_num            Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_geom_repr_mesh_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 int               l_num[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Finalize the geometric representation mesh
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided (MPI collective communications)
 *
 * \param [in]  ma                Mesh adaptation workflow
 *
 */

void
PDM_Mesh_adapt_geom_mesh_repr_finalize
(
 PDM_Mesh_adapt_t *ma
);


/**
 * \brief Set ridges
 *
 * This functions defines the ridges in the mesh of the geometric representation
 * Under the form of the list of couples of points (vtx1/vtx2)
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  n_ridge          Number of ridges
 * \param [in]  ridges           List of ridges
 *
 */

void
PDM_Mesh_adapt_geom_repr_ridge_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         n_ridge,
 int               ridges
);


/**
 * \brief Set corners
 *
 * This functions defines the corners in the mesh of the geometric representation
 * under the form of the list of points
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  i_dom            Domain identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  n_corner         Number of corners
 * \param [in]  corners          List of corners
 *
 */

void
PDM_Mesh_adapt_geom_repr_corner_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         n_corner,
 int               corners
);

/**
 * \brief Load CAD
 *
 * This functions loads the CAD file
 *
 * \param [in]  ma               Mesh adaptation workflow
 * \param [in]  cad_file         Computer Aided Desgin file
 *
 */

void
PDM_Mesh_adapt_geom_cad_load
(
 PDM_Mesh_adapt_t *ma,
 const char*       cad_file
);


/*----------------------------------------------------------------------------*
 *
 * Functions about geometric criteria
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set the geometric adaptation criterion
 *
 *  Depends on the type of selected criterion.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  criterion    Geometric adaptation criterion (size = n_elt or n_vtx)
 *                              - type int * : local refinment criterion (size = n_elt)
 *                              - type double * : metric (size = 6 * n_nodes)
 */

void
PDM_Mesh_adapt_geom_criterion_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 void             *criterion
);


/*----------------------------------------------------------------------------*
 *
 * Functions about specific tool parameters
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set a tool parameter
 *
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  name         Parameter name
 * \param [in]  value        Parameter value

 */

void
PDM_Mesh_adapt_tool_param_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const  char      *name,
 void             *value
);


/*----------------------------------------------------------------------------*
 *
 * About field family
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set
 *
 * \param [in]  ma             Mesh adaptation workflow
 * \param [in]  dof_location   Location of the Degrees of freedom
 * \param [in]  field_interp   Intepolation method
 *
 * \return  i_field_family  field family identifier (used by transfer
 *                                                   field functions)
 */

int
PDM_Mesh_adapt_field_family_add
(
 PDM_Mesh_adapt_t                 *ma,
 PDM_Mesh_adapt_dof_location_t    *dof_location,
 PDM_Mesh_adapt_field_interp_t    *field_interp
);


// TODO : maque fonction PDM_Mesh_adapt_field_user_loc_set

/*----------------------------------------------------------------------------*
 *
 * Compute taget mesh
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Computes target mesh, splits it, redistributes it  and computes
 *        MPI communication graphs between source and target meshes according to
 *        the properties of field families.
 *
 * \param [in]  ma           Mesh adaptation workflow
 *
 */

void
PDM_Mesh_adapt_compute
(
 PDM_Mesh_adapt_t *ma
);


/*----------------------------------------------------------------------------*
 *
 * Functions about target mesh
 *
 *----------------------------------------------------------------------------*/


/**
 * \brief Get target mesh global size
 *
 * \param [in]   ma                Mesh adaptation workflow
 * \param [in]   i_dom             Domain identifier
 * \param [out]  g_n_vtx           Global number of nodes
 * \param [out]  g_n_elt           Global number of elements
 * \param [out]  g_n_face          Global number of faces
 * \param [out]  g_n_boundary_elt  Global number of elements on the boundary
 *
 */

void
PDM_Mesh_adapt_tgt_global_size
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 PDM_g_num_t      *g_n_vtx,
 PDM_g_num_t      *g_n_elt,
 PDM_g_num_t      *g_n_face,
 PDM_g_num_t      *g_n_boundary_elt
);


/**
 * \brief Get target mesh local size
 *
 * \param [in]   ma              Mesh adaptation workflow
 * \param [in]   i_dom           Domain identifier
 * \param [in]   i_part          Current partition
 * \param [out]  n_vtx           Number of vertices
 * \param [out]  n_elt           Number of elements
 * \param [out]  n_fac           Number of faces
 * \param [out]  n_boundary_elt  Number of elements on the boundary
 *
 */

void
PDM_Mesh_adapt_tgt_size
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 int              *n_vtx,
 int              *n_elt,
 int              *n_fac,
 int              *n_boundary_elt
);


/**
 * \brief Get target mesh vertices.
 *
 * \param [in]   ma           Mesh adaptation workflow
 * \param [in]   i_dom        Domain identifier
 * \param [in]   i_part       Current partition
 * \param [out]  coord        Coordinates (size = 3 * \p n_vtx)
 * \param [out]  g_num        Global element number
 *
 */

void
PDM_Mesh_adapt_tgt_vtx_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 double           *coord[],
 PDM_g_num_t      *g_num[]
);


/**
 * \brief Get the number of blocks in the target mesh
 *
 * \param [in]  ma      Mesh adaptation workflow
 * \param [in]  i_dom   Domain identifier
 *
 * \return      Number of blocks
 */

int
PDM_Mesh_adapt_tgt_n_block_get
(
 PDM_Mesh_adapt_t   *ma,
 const int         i_dom
);


/**
 * \brief Get the type of block
 *
 * \param [in]  ma            Mesh adaptation workflow
 * \param [in]  i_dom           Domain identifier
 * \param [in]  i_block       Block identifier
 *
 * \return   Type of the current block
 */

PDM_Mesh_nodal_elt_t
PDM_Mesh_adapt_tgt_block_type_get
(
 PDM_Mesh_adapt_t  *ma,
 const int         i_dom,
 const int         i_block
);


/**
 * \brief Get a standard block of the target mesh.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [out] n_elts       Number of elements
 * \param [out] elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [out] l_num        Local element number in the partition (size = n_elts)
 * \param [out] g_num        Global element number (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_tgt_block_std_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              *elt_vtx[],
 int              *l_num[],
 PDM_g_num_t      *g_num[]
);


/**
 * \brief Get a generic high order block of the target mesh.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [out] n_elts       Number of elements
 * \param [out] elt_node     Connectivity (size = n_nodes_elt * n_elts)
 * \param [out] l_num        Local element number in the partition (size = n_elts)
 * \param [out] g_num        Global element number (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_tgt_block_ho_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              *elt_node[],
 int              *l_num[],
 PDM_g_num_t      *g_num[]
);


/**
 * \brief Get a polygon block of the target mesh.
 *
 * \param [in]  ma             Mesh adaptation workflow
 * \param [in]  i_dom          Domain identifier
 * \param [in]  i_part         Partition identifier
 * \param [in]  i_block        Block identifier
 * \param [out] n_elts         Number of elements
 * \param [out] face_vtx_idx   Connectivity index (\p connec_id[0] = 0 and
 *                             size = \p n_elts + 1)
 * \param [out] face_vtx       face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [out] l_num          Local element number in the partition (size = n_elts)
 * \param [out] g_num          Global element number (or NULL)
 *
 *
 */

void
PDM_Mesh_adapt_tgt_block_f_poly_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              *face_vtx_idx[],
 int              *face_vtx[],
 int              *l_num[],
 PDM_g_num_t      *g_num[]
);


/**
 *
 * \brief Set the connectivity of a polyhedron block of the source mesh.
 *
 * Connectivity is supposed to be oriented. Connectivity can be oriented by calling
 * PDM_cellface_orient
 *
 * \param [in]  ma                Mesh adaptation workflow
 * \param [in]  i_dom             Domain identifier
 * \param [in]  i_part            Partition identifier
 * \param [in]  i_block           Block identifier
 * \param [out]  n_elts            Number of elements
 * \param [out]  n_faces           Number of faces
 * \param [out]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [out]  face_vtx          Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [out]  face_l_num        Face local element number
 * \param [out]  face_g_num        Face global element number
 * \param [out]  cell_face_idx     Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [out]  cell_face         Polyhedron to face connectivity
 *                                The connectivity is oriented :
 *                                  - > 0 if outgoing normal,
 *                                  - < 0 otherwise
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [out]  face_cell         Face to polyhedron connectivity
 *                                  - left value  : outgoing normal,
 *                                  - right value : incoming normal
 *                                (size = 2 * \p n_faces)
 * \param [out]  cell_l_num        Cell local element number
 * \param [out]  cell_g_num        Cell global element number
 *
 *
 */

void
PDM_Mesh_adapt_tgt_block_c_poly_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              *n_faces,
 int              *face_vtx_idx[],
 int              *face_vtx[],
 int              *face_l_num[],
 PDM_g_num_t      *face_g_num[],
 int              *cell_face_idx[],
 int              *cell_face[],
 int              *face_cell[],
 int              *cell_l_num[],
 PDM_g_num_t      *cell_g_num[]
);


/**
 * \brief Get the number of blocks in the target boundary mesh
 *
 * \param [in]  ma      Mesh adaptation workflow
 * \param [in]  i_dom   Domain identifier
 *
 * \return      Number of boundary blocks
 */

int
PDM_Mesh_adapt_tgt_n_boundary_block_get
(
 PDM_Mesh_adapt_t   *ma,
 const int         i_dom
);


/**
 * \brief Get the type of a boundary block block
 *
 * \param [in]  ma            Mesh adaptation workflow
 * \param [in]  i_dom         Domain identifier
 * \param [in]  i_block       Block identifier
 *
 * \return   Type of the current block
 */

PDM_Mesh_nodal_elt_t
PDM_Mesh_adapt_tgt_boundary_block_type_get
(
 PDM_Mesh_adapt_t  *ma,
 const int         i_dom,
 const int         i_block
);


/**
 * \brief Get a standard block of the target boundary mesh.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [out] n_elts       Number of elements
 * \param [out] elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [out] l_num        Local element number in the partition (size = n_elts)
 * \param [out] g_num        Global element number (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_tgt_boundary_block_std_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              (*elt_vtx)[],
 int              (*l_num)[],
 PDM_g_num_t      (*g_num)[]
);


/**
 * \brief Get a generic high order block of the target boundary mesh.
 *
 * \param [in]  ma           Mesh adaptation workflow
 * \param [in]  i_dom        Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [out] n_elts       Number of elements
 * \param [out] elt_node     Connectivity (size = n_nodes_elt * n_elts)
 * \param [out] l_num        Local element number in the partition (size = n_elts)
 * \param [out] g_num        Global element number (size = n_elts)
 *
 */

void
PDM_Mesh_adapt_tgt_boundary_block_ho_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              (*elt_node)[],
 int              (*l_num)[],
 PDM_g_num_t      (*g_num)[]
);


/**
 * \brief Get a polygon block of the target boundary mesh.
 *
 * \param [in]  ma             Mesh adaptation workflow
 * \param [in]  i_dom          Domain identifier
 * \param [in]  i_part         Partition identifier
 * \param [in]  i_block        Block identifier
 * \param [out] n_elts         Number of elements
 * \param [out] face_vtx_idx   Connectivity index (\p connec_id[0] = 0 and
 *                             size = \p n_elts + 1)
 * \param [out] face_vtx       face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [out] l_num          Local element number in the partition (size = n_elts)
 * \param [out] g_num          Global element number (or NULL)
 *
 *
 */

void
PDM_Mesh_adapt_tgt_boundary_block_f_poly_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 const int         i_block,
 int              *n_elts,
 int              (*face_vtx_idx)[],
 int              (*face_vtx)[],
 int              (*l_num)[],
 PDM_g_num_t      (*g_num)[]
);


/**
 * \brief Get ancestor in the source mesh (only for \ref PDM_MESH_ADAPT_REFINMENT method)
 *
 * \param [in]  ma                           Mesh adaptation workflow
 * \param [in]  i_dom                        Domain identifier
 * \param [in]  i_part                       Partition identifier
 * \param [out] g_num_vtx_ancestor           Ancestor global vertex number
 *                                           (size = \p n_elts)
 * \param [out] g_num_elt_ancestor           Ancestor global element number
 *                                           (size = \p n_elts)
 * \param [out] g_num_boundary_elt_ancestor  Ancestor global boundary element number
 *                                           (size = \p n_boundary_elts)
 *
 */

void
PDM_Mesh_adapt_tgt_ancestor_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_dom,
 const int         i_part,
 PDM_g_num_t      (*g_num_vtx_ancestor)[],
 PDM_g_num_t      (*g_num_elt_ancestor)[],
 PDM_g_num_t      (*g_num_boundary_elt_ancestor)[]
);

/**
 * \brief Get a entity group
 *
 * An entity can be contained in several groups, the notion of group  is more general than
 * the notion of boudary condition
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  g_entity    Entity group
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  i_group     Group identifier
 * \param [out]  n_entities Number of entities in the group
 * \param [out]  entities   List of entities
 * \param [out]  g_num      Global element number in the group (or NULL)
 *
 */

void
PDM_Mesh_adapt_tgt_entity_group_get
(
 PDM_Mesh_adapt_t              *ma,
 PDM_Mesh_adapt_group_entity_t g_entity,
 const int                     i_dom,
 const int                     i_part,
 const int                     i_group,
 const int                    *n_entities,
 int                         (*entities)[],
 PDM_g_num_t                 (*g_num)[]
);


/**
 * \brief Get an intra-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_graph     Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [out]  n_elt_graph Number of elements in the
 * \param [out]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [out]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - Lists of 3 values :
 *                                  + Connected process
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_tgt_intra_dom_graph_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_graph,
 const int         i_part,
 int              *n_elt_graph,
 int              (*graph_idx)[],
 int              (*graph)[]
);


/**
 * \brief Get an inter-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_graph     Graph identifier
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [out]  n_elt_graph Number of elements in the
 * \param [out]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [out]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - Lists of 4 values :
 *                                  + Connected process
 *                                  + Connected domain
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_tgt_inter_dom_graph_get
(
 PDM_Mesh_adapt_t  *ma,
 const int          i_graph,
 const int          i_dom,
 const int          i_part,
 int               *n_elt_graph,
 int              (*graph_idx)[],
 int              (*graph)[]
);


/**
 * \brief Set an intra-domain graph communication between partitions
 *
 * \param [in]  ma          Mesh adaptation workflow
 * \param [in]  i_perio     Periodicity identifier
 * \param [in]  i_dom       Domain identifier
 * \param [in]  i_part      Partition identifier
 * \param [in]  n_elt_graph Number of elements in the
 * \param [in]  graph_idx   Element index in \p graph
 *                          (size = \a n_elt_graph)
 * \param [in]  graph       For each element graph :
 *                             - Local entity number in the partition
 *                             - List of 4 values :
 *                                  + Connected process
 *                                  + Connected domain
 *                                  + Local partition number in the
 *                                    connected process
 *                                  + Local element number in the local
 *                                    partition number
 *
 */

void
PDM_Mesh_adapt_tgt_perio_graph_get
(
 PDM_Mesh_adapt_t *ma,
 const int         i_perio,
 const int         i_dom,
 const int         i_part,
 int              *n_elt_graph,
 int             (*graph_idx)[],
 int             (*graph)[]
);

/*----------------------------------------------------------------------------*
 *
 * Functions about interpolation callback
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set a user structure that can be used in the callback called from
 * the source mesh
 *
 * \param [in]  ma              Mesh adaptation workflow
 * \param [in]  i_dom           Domain identifier
 * \param [in]  user_struct     user_struct
 *
 */

void
PDM_Mesh_adapt_callback_src_user_struct_set
(
 PDM_Mesh_adapt_t   *ma,
 const int         i_dom,
 void               *user_struct
);


/**
 * \brief Set a user structure that can be used in the callback called from
 * the target mesh
 *
 * \param [in]  ma              Mesh adaptation workflow
 * \param [in]  i_dom           Domain identifier
 * \param [in]  user_struct     user_struct
 *
 */

void
PDM_Mesh_adapt_callback_tgt_user_struct_set
(
 PDM_Mesh_adapt_t   *ma,
 const int         i_dom,
 void               *user_struct
);


/**
 * \brief Set the function used to interpolate fields from the location of
 * the degrees of freedom into the src mesh
 *
 * \param [in]  ma                           Mesh adaptation workflow
 * \param [in]  callback                     Partition identifier
 *
 */

void
PDM_Mesh_adapt_interp_from_location_set
(
 PDM_Mesh_adapt_t                     *ma,
 PDM_Mesh_adapt_interp_from_location_t callback
);


/**
 * \brief Set the function used to interpolate fields from the intersectio of
 * source ans target meshes
 *
 * \param [in]  ma                           Mesh adaptation workflow
 * \param [in]  callback                     Partition identifier
 *
 */

void
PDM_Mesh_adapt_interp_from_intersect_set
(
 PDM_Mesh_adapt_t                     *ma,
 PDM_Mesh_adapt_interp_from_intersect_t callback
);


/**
 * \brief Set the function used to interpolate fields from the 'k' closest
 * degrees of freedom in the source mesh
 *
 * \param [in]  ma                           Mesh adaptation workflow
 * \param [in]  callback                     Partition identifier
 *
 */

void
PDM_Mesh_adapt_interp_from_closest_points_set
(
 PDM_Mesh_adapt_t                        *ma,
 PDM_Mesh_adapt_interp_from_closest_pts_t callback
);

/*----------------------------------------------------------------------------*
 *
 * Functions about field transfer
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Send a field to the target mesh (non-blocking communications)
 *
 * \param [in]  ma                  Mesh adaptation workflow
 * \param [in]  i_dom               Domain identifier
 * \param [in]  i_field_family      Family of the field
 *                                  (from \ref PDM_Mesh_adapt_field_family_add)
 * \param [in]  stride              Stride of the field
 * \param [in]  src_field           Array of pointers to field data (size = \n_part)
 * \param [out] request             Request is used by \ref PDM_Mesh_adapt_src_wait_field
 *
 */

void
PDM_Mesh_adapt_src_field_issend
(
 PDM_Mesh_adapt_t  *ma,
 const int         i_dom,
 const int         i_field_family,
 const int         stride,
 double           *src_field[],
 int              *request
);


/**
 * \brief Receive a field from the source mesh (non-blocking communications)
 *
 * \param [in]  ma                  Mesh adaptation workflow
 * \param [in]  i_dom               Domain identifier
 * \param [in]  i_field_family      Family of the field
 *                                  (from \ref PDM_Mesh_adapt_field_family_add)
 * \param [in]  stride              Stride of the field
 * \param [in]  tgt_field           Array of pointers to field data (size = \n_part)
 * \param [out] request             Request is used by \ref PDM_Mesh_adapt_src_wait_field
 *
 */

void
PDM_Mesh_adapt_tgt_field_irecv
(
 PDM_Mesh_adapt_t  *ma,
 const int         i_dom,
 const int         i_field_family,
 const int         stride,
 double           *tgt_field[],
 int              *request
);


/**
 * \brief Wait for the field exchange
 *
 * \param [in] ma           Mesh adaptation workflow
 * \param [in] request      Request from \ref PDM_Mesh_adapt_tgt_field_irecv or
 *                          \ref PDM_Mesh_adapt_tgt_field_issend
 *
 */

void
PDM_Mesh_adapt_wait
(
 PDM_Mesh_adapt_t  *ma,
 int              request
);

#ifdef __cplusplus
}
#endif

#endif // PDM_CLOSEST_POINTS_H
