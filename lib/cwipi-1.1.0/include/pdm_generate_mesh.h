#ifndef __PDM_GENERATE_MESH_H__
#define __PDM_GENERATE_MESH_H__

#include "pdm.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// TO DO: deformation methods on partitioned data(warning for randomization, same seed)

/**
 *
 * \brief Create a simple partitioned sphere mesh (2D).
 *
 * \param [in]   comm        MPI communicator
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_sphere_simplified
(
 const PDM_MPI_Comm   comm,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitioned rectangle mesh (2D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   n_vtx_seg   Number of vertices along each side of the rectangle
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_rectangle_simplified
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t    n_vtx_seg,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitioned ball mesh (3D).
 *
 * \param [in]   comm        MPI communicator
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_ball_simplified
(
 const PDM_MPI_Comm   comm,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitioned parallelepiped mesh (3D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   n_vtx_seg   Number of vertices along each side of the parallelepiped
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_parallelepiped_simplified
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t    n_vtx_seg,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a partitioned sphere mesh (2D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering High order nodes ordering type
 * \param [in]  radius      Radius of the sphere
 * \param [in]  center_x    x-coordinate of the sphere center
 * \param [in]  center_y    y-coordinate of the sphere center
 * \param [in]  center_z    z-coordinate of the sphere center
 * \param [in]  n_u         Number of points in longitude
 * \param [in]  n_v         Number of points in latitude
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitioning method
 *
 * \return PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_sphere
(
 const PDM_MPI_Comm           comm,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    order,
 const char                  *ho_ordering,
 const double                 radius,
 const double                 center_x,
 const double                 center_y,
 const double                 center_z,
 const PDM_g_num_t            n_u,
 const PDM_g_num_t            n_v,
 const int                    n_part,
 const PDM_split_dual_t       part_method
);

/**
 *
 * \brief Create a partitioned rectangle mesh (2D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering High order nodes ordering type
 * \param [in]  xmin        x-coordinate of the rectangle minimum corner
 * \param [in]  ymin        y-coordinate of the rectangle minimum corner
 * \param [in]  zmin        z-coordinate of the rectangle minimum corner
 * \param [in]  lengthx     Length of the rectangle in the x-direction
 * \param [in]  lengthy     Length of the rectangle in the y-direction
 * \param [in]  n_x         Number of points in the x-direction
 * \param [in]  n_y         Number of points in the y-direction
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitioning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_rectangle
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char             *ho_ordering,
 double                  xmin,
 double                  ymin,
 double                  zmin,
 double                  lengthx,
 double                  lengthy,
 PDM_g_num_t             n_x,
 PDM_g_num_t             n_y,
 const int               n_part,
 const PDM_split_dual_t  part_method
);

/**
 *
 * \brief Create a partitioned ball mesh (3D).
 *
 * \param [in]  comm            MPI communicator
 * \param [in]  elt_type        Mesh element type
 * \param [in]  order           Mesh element order
 * \param [in]  ho_ordering     High order nodes ordering type
 * \param [in]  radius          Radius of the ball
 * \param [in]  hole_radius     Radius of the hole of the ball
 * \param [in]  center_x        x-coordinate of the ball center
 * \param [in]  center_y        y-coordinate of the ball center
 * \param [in]  center_z        z-coordinate of the ball center
 * \param [in]  n_x             Number of vertices on segments in x-direction
 * \param [in]  n_y             Number of vertices on segments in y-direction
 * \param [in]  n_z             Number of vertices on segments in z-direction
 * \param [in]  n_layer         Number of extrusion layers
 * \param [in]  geometric_ratio Geometric ratio for layer thickness
 * \param [in]  n_part          Number of mesh partitions
 * \param [in]  part_method     Mesh partitioning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_ball
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char             *ho_ordering,
 const double            radius,
 const double            hole_radius,
 const double            center_x,
 const double            center_y,
 const double            center_z,
 const PDM_g_num_t       n_x,
 const PDM_g_num_t       n_y,
 const PDM_g_num_t       n_z,
 const PDM_g_num_t       n_layer,
 const double            geometric_ratio,
 const int               n_part,
 const PDM_split_dual_t  part_method
);

/**
 *
 * \brief Create a partitioned parallelepiped mesh (3D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering High order nodes ordering type
 * \param [in]  xmin        x-coordinate of the rectangle minimum corner
 * \param [in]  ymin        y-coordinate of the rectangle minimum corner
 * \param [in]  zmin        z-coordinate of the rectangle minimum corner
 * \param [in]  lengthx     Length of the rectangle in the x-direction
 * \param [in]  lengthy     Length of the rectangle in the y-direction
 * \param [in]  lengthz     Length of the rectangle in the z-direction
 * \param [in]  n_x         Number of points in the x-direction
 * \param [in]  n_y         Number of points in the y-direction
 * \param [in]  n_z         Number of points in the z-direction
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitioning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_parallelepiped
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char             *ho_ordering,
 double                  xmin,
 double                  ymin,
 double                  zmin,
 double                  lengthx,
 double                  lengthy,
 double                  lengthz,
 PDM_g_num_t             n_x,
 PDM_g_num_t             n_y,
 PDM_g_num_t             n_z,
 const int               n_part,
 const PDM_split_dual_t  part_method
);


/**
 *
 * \brief Create a partitioned rectangular mesh (2D) with descending connectivities
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   elt_type       Element type
 * \param [in]   xmin           Minimal x-coordinate
 * \param [in]   ymin           Minimal y-coordinate
 * \param [in]   zmin           Minimal z-coordinate
 * \param [in]   lengthx        Length of the rectangle in the x-direction
 * \param [in]   lengthy        Length of the rectangle in the y-direction
 * \param [in]   n_x            Number of points in the x-direction
 * \param [in]   n_y            Number of points in the y-direction
 * \param [in]   n_part         Number of partitions
 * \param [in]   part_method    Partitioning method
 * \param [in]   random_factor  Randomization factor (between 0 and 1)
 * \param [out]  pn_vtx         Number of vertices (size = \p n_part)
 * \param [out]  pn_edge        Number of edges (size = \p n_part)
 * \param [out]  pn_face        Number of faces (size = \p n_part)
 * \param [out]  pvtx_coord     Vertex coordinates (for each part, size = \p pn_vtx)
 * \param [out]  pedge_vtx      Edge->vertex connectivity (for each part, size = 2 * \p pn_edge)
 * \param [out]  pface_edge_idx Index of face->edge connectivity (for each part, size = \p pn_face + 1)
 * \param [out]  pface_edge     Face->edge connectivity (for each part, size = \p face_edge_idx[\p pn_face])
 * \param [out]  pface_vtx      Face->vertex connectivity (for each part, size = \p face_edge_idx[\p pn_face])
 * \param [out]  pvtx_ln_to_gn  Vertex global ids (for each part, size = \p pn_vtx)
 * \param [out]  pedge_ln_to_gn Edge global ids (for each part, size = \p pn_edge)
 * \param [out]  pface_ln_to_gn Face global ids (for each part, size = \p pn_face)
 *
 * \note Admissible values for \p elt_type:
 *   - \p PDM_MESH_NODAL_TRIA3   : triangles
 *   - \p PDM_MESH_NODAL_QUAD4   : quadrangles
 *   - \p PDM_MESH_NODAL_POLY_2D : mixed polygons (triangles, quadrangles and octagons)
 *
 */

void
PDM_generate_mesh_rectangle_ngon
(
 const PDM_MPI_Comm            comm,
 const PDM_Mesh_nodal_elt_t    elt_type,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  lengthx,
 const double                  lengthy,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const double                  random_factor,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn
);

/**
 *
 * \brief Create a partitioned sphere mesh (2D) with descending connectivities.
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   elt_type       Element type
 * \param [in]   order          Element order
 * \param [in]   ho_ordering    Ordering of nodes of the HO element
 * \param [in]   radius         Radius of the sphere
 * \param [in]   center_x       x-coordinate of the sphere center
 * \param [in]   center_y       y-coordinate of the sphere center
 * \param [in]   center_z       z-coordinate of the sphere center
 * \param [in]   n_u            Number of vertices in the u-direction
 * \param [in]   n_v            Number of vertices in the v-direction
 * \param [in]   n_part         Number of partitions
 * \param [in]   part_method    Partitioning method
 * \param [in]   pn_vtx         Number of vertices
 * \param [in]   pn_edge        Number of edges
 * \param [in]   pn_face        Number of faces
 * \param [in]   pvtx_coord     Vertex coordinates
 * \param [in]   pedge_vtx      edge->vertex connectivity
 * \param [in]   pface_edge_idx Index of face->edge connectivity
 * \param [in]   pface_edge     face->edge connectivity
 * \param [in]   pface_vtx      face->vertex connectivity
 * \param [in]   pvtx_ln_to_gn  Vertex global number
 * \param [in]   pedge_ln_to_gn Edge global number
 * \param [in]   pface_ln_to_gn Face global number
 *
 */

void
PDM_generate_mesh_sphere_ngon
(
 const PDM_MPI_Comm           comm,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    order,
 const char                  *ho_ordering,
 const double                 radius,
 const double                 center_x,
 const double                 center_y,
 const double                 center_z,
 const PDM_g_num_t            n_u,
 const PDM_g_num_t            n_v,
 const int                    n_part,
 const PDM_split_dual_t       part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn
);

/**
 *
 * \brief Create a partitioned ball mesh (3D) with descending connectivities.
 *
 * \param [in]  comm                      MPI communicator
 * \param [in]  elt_type                  Mesh element type
 * \param [in]  order                     Mesh element order
 * \param [in]  ho_ordering               High order nodes ordering type
 * \param [in]  radius                    Radius of the ball
 * \param [in]  hole_radius               Radius of the hole of the ball
 * \param [in]  center_x                  x-coordinate of the ball center
 * \param [in]  center_y                  y-coordinate of the ball center
 * \param [in]  center_z                  z-coordinate of the ball center
 * \param [in]  n_x                       Number of vertices on segments in x-direction
 * \param [in]  n_y                       Number of vertices on segments in y-direction
 * \param [in]  n_z                       Number of vertices on segments in z-direction
 * \param [in]  n_layer                   Number of extrusion layers
 * \param [in]  geometric_ratio           Geometric ratio for layer thickness
 * \param [in]  n_part                    Number of mesh partitions
 * \param [in]  part_method               Mesh partitioning method
 * \param [out] pn_vtx                    Number of vertices
 * \param [out] pn_edge                   Number of edges
 * \param [out] pn_face                   Number of faces
 * \param [out] pvtx_coord                Vertex coordinates
 * \param [out] pedge_vtx                 Edge->vertex connectivity
 * \param [out] pface_edge_idx            Index of face->edge connectivity
 * \param [out] pface_edge                Face->edge connectivity
 * \param [out] pface_vtx                 Face->vertex connectivity
 * \param [out] pvtx_ln_to_gn             Vertex global number
 * \param [out] pedge_ln_to_gn            Edge global number
 * \param [out] pface_ln_to_gn            Face global number
 * \param [out] pn_surface                Number of surfaces
 * \param [out] psurface_face_idx         Surface->face connectivity index
 * \param [out] psurface_face             Surface->face connectivity
 * \param [out] psurface_face_ln_to_gn    Surface->face connectivity with global numbers
 *
 */

void
PDM_generate_mesh_ball_ngon
(
 const PDM_MPI_Comm            comm,
 PDM_Mesh_nodal_elt_t          elt_type,
 int                           order,
 const char                   *ho_ordering,
 const double                  radius,
 const double                  hole_radius,
 const double                  center_x,
 const double                  center_y,
 const double                  center_z,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const PDM_g_num_t             n_z,
 const PDM_g_num_t             n_layer,
 const double                  geometric_ratio,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 int                         **pn_cell,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pcell_ln_to_gn,
 int                         **pn_surface,
 int                        ***psurface_face_idx,
 int                        ***psurface_face,
 PDM_g_num_t                ***psurface_face_ln_to_gn
 );

/**
 *
 * \brief Create a partitioned parallelepiped mesh (3D) with descending connectivities.
 *
 * \param [in]  comm                      MPI communicator
 * \param [in]  elt_type                  Mesh element type
 * \param [in]  order                     Mesh element order
 * \param [in]  ho_ordering               High order nodes ordering type
 * \param [in]  radius                    Radius of the ball
 * \param [in]  hole_radius               Radius of the hole of the ball
 * \param [in]  center_x                  x-coordinate of the ball center
 * \param [in]  center_y                  y-coordinate of the ball center
 * \param [in]  center_z                  z-coordinate of the ball center
 * \param [in]  n_x                       Number of vertices on segments in x-direction
 * \param [in]  n_y                       Number of vertices on segments in y-direction
 * \param [in]  n_z                       Number of vertices on segments in z-direction
 * \param [in]  n_layer                   Number of extrusion layers
 * \param [in]  geometric_ratio           Geometric ratio for layer thickness
 * \param [in]  n_part                    Number of mesh partitions
 * \param [in]  part_method               Mesh partitioning method
 * \param [out] pn_vtx                    Number of vertices
 * \param [out] pn_edge                   Number of edges
 * \param [out] pn_face                   Number of faces
 * \param [out] pvtx_coord                Vertex coordinates
 * \param [out] pedge_vtx                 Edge->vertex connectivity
 * \param [out] pface_edge_idx            Index of face->edge connectivity
 * \param [out] pface_edge                Face->edge connectivity
 * \param [out] pface_vtx                 Face->vertex connectivity
 * \param [out] pvtx_ln_to_gn             Vertex global number
 * \param [out] pedge_ln_to_gn            Edge global number
 * \param [out] pface_ln_to_gn            Face global number
 * \param [out] pn_surface                Number of surfaces
 * \param [out] psurface_face_idx         Surface->face connectivity index
 * \param [out] psurface_face             Surface->face connectivity
 * \param [out] psurface_face_ln_to_gn    Surface->face connectivity with global numbers
 * \param [out] pn_ridge                  Number of ridges
 * \param [out] pridge_edge_idx           Ridge->edge connectivity index
 * \param [out] pridge_edge               Ridge->edge connectivity
 * \param [out] pridge_edge_ln_to_gn      Ridge->edge connectivity with global numbers
 *
 */

void
PDM_generate_mesh_parallelepiped_ngon
(
 const PDM_MPI_Comm            comm,
 PDM_Mesh_nodal_elt_t          elt_type,
 int                           order,
 const char                   *ho_ordering,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  lengthx,
 const double                  lengthy,
 const double                  lengthz,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const PDM_g_num_t             n_z,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 int                         **pn_cell,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pcell_ln_to_gn,
 int                         **pn_surface,
 int                        ***psurface_face_idx,
 int                        ***psurface_face,
 PDM_g_num_t                ***psurface_face_ln_to_gn,
 int                         **pn_ridge,
 int                        ***pridge_edge_idx,
 int                        ***pridge_edge,
 PDM_g_num_t                ***pridge_edge_ln_to_gn
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_GENERATE_MESH__ */
