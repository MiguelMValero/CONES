/*
 * \file
 */

#ifndef PDM_MESH_LOCATION_H
#define PDM_MESH_LOCATION_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_to_part.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_mesh_location_t PDM_mesh_location_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute the location of point clouds inside a mesh
 *
 * \param [in]   n_point_cloud  Number of point clouds
 * \param [in]   comm           MPI communicator
 * \param [in]   owner          Ownership
 *
 * \return     Pointer to \ref PDM_mesh_location_t instance
 *
 */

PDM_mesh_location_t*
PDM_mesh_location_create
(
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  n_part
);



/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates (size : 3 * \p n_points)
 * \param [in]   gnum            Point global ids (size : \p n_points)
 *
 */

void
PDM_mesh_location_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);


/**
 *
 * \brief Get a point cloud
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 * \param [out]  n_points        Number of points
 * \param [out]  coords          Point coordinates
 * \param [out]  gnum            Point global ids
 *
 */

void
PDM_mesh_location_cloud_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       int                  *n_points,
       double              **coords,
       PDM_g_num_t         **gnum
);


/**
 *
 * \brief Get the number of located points
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 *
 * \return     Number of located points
 *
 */

int
PDM_mesh_location_n_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the number of unlocated points
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 *
 * \return     Number of unlocated points
 *
 */

int
PDM_mesh_location_n_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the list of unlocated points
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 *
 * \return     List of unlocated points (1-based ids)
 *
 */

int *
PDM_mesh_location_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the list of located points
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 *
 * \return     List of located points (1-based ids)
 *
 */

int *
PDM_mesh_location_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);



/**
 *
 * \brief Set the nodal mesh
 *
 * \param [in]   ml             Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   mesh_nodal     Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 PDM_mesh_location_t   *ml,
 PDM_part_mesh_nodal_t *mesh_nodal
);


/**
 *
 * \brief Set the number of partitions of the mesh
 *
 * \param [in]   ml             Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   n_part         Number of partitions
 *
 */

void
PDM_mesh_location_mesh_n_part_set
(
       PDM_mesh_location_t *ml,
 const int                  n_part
);


/**
 *
 * \brief Set a *volume* mesh partition
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_part          Partition identifier
 * \param [in]   n_cell          Number of cells
 * \param [in]   cell_face_idx   Index for cell -> face connectivity (size : \p n_cell + 1)
 * \param [in]   cell_face       Cell -> face connectivity (size : \p cell_face_idx[\p n_cell])
 * \param [in]   cell_ln_to_gn   Cell global ids (size : \p n_cell)
 * \param [in]   n_face          Number of faces
 * \param [in]   face_vtx_idx    Index for face -> vertex connectivity (size : \p n_face + 1)
 * \param [in]   face_vtx        Face -> vertex connectivity (size : \p face_vtx_idx[\p n_face])
 * \param [in]   face_ln_to_gn   Face global ids (size : \p n_face)
 * \param [in]   n_vtx           Number of vertices
 * \param [in]   coords          Vertex coordinates (size : 3 * \p n_vtx)
 * \param [in]   vtx_ln_to_gn    Vertex global ids (size : \p n_vtx)
 *
 */

void
PDM_mesh_location_part_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_face_idx,
 const int                 *cell_face,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_face,
 const int                 *face_vtx_idx,
 const int                 *face_vtx,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Set a *volume* mesh partition defined by nodal connectivity
 *
 * The mesh is assumed to contain only standard elements
 * (tetrahedra, pyramids, prisms, hexahedra).
 *
 * \param [in]   ml            Pointer to \ref PDM_mesh_location instance
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_vtx_idx  Index in the cell -> vertex connectivity
 * \param [in]   cell_vtx      Cell -> vertex connectivity
 * \param [in]   cell_ln_to_gn Global cell ids
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Vertex global ids
 *
 */
void
PDM_mesh_location_nodal_part_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_vtx_idx,
 const int                 *cell_vtx,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Select cells to extract before computing location
 *
 * \param [in]   ml                     Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_part                 Partition identifier
 * \param [in]   is_elmt_select_by_user Flag to determine which cells must be extracted (size : *n_cell*)
 *
 */
void
PDM_mesh_location_user_extract_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                 *is_elmt_select_by_user
);



/**
 *
 * \brief Set a *surface* mesh partition
 *
 * \param [in]   ml             Pointer to \ref PDM_mesh_location instance
 * \param [in]   i_part         Partition to define
 * \param [in]   n_face         Number of faces
 * \param [in]   face_edge_idx  Index for face -> edge connectivity (size : \p n_face + 1)
 * \param [in]   face_edge      Face -> edge connectivity (size : \p face_edge_idx[\p n_face])
 * \param [in]   face_ln_to_gn  Face global ids (size : \p n_face)
 * \param [in]   n_edge         Number of edges
 * \param [in]   edge_vtx       Edge -> vertex connectivity
 * \param [in]   n_vtx          Number of vertices
 * \param [in]   coords         Coordinates
 * \param [in]   vtx_ln_to_gn   Vertex global ids
 *
 */

void
PDM_mesh_location_part_set_2d
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_face,
 const int                 *face_edge_idx,
 const int                 *face_edge,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_edge,
 const int                 *edge_vtx,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Set a *surface* mesh partition with nodal connectivity
 *
 * \param [in]   ml             Pointer to \ref PDM_mesh_location instance
 * \param [in]   i_part         Partition to define
 * \param [in]   n_face         Number of faces
 * \param [in]   face_vtx_idx   Index for face -> vertex connectivity
 * \param [in]   face_vtx       Face -> vertex connectivity
 * \param [in]   face_ln_to_gn  Face global ids
 * \param [in]   n_vtx          Number of vertices
 * \param [in]   coords         Coordinates
 * \param [in]   vtx_ln_to_gn   Vertex global ids
 *
 */

void
PDM_mesh_location_nodal_part_set_2d
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_face,
 const int                 *face_vtx_idx,
 const int                 *face_vtx,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Set the relative tolerance for bounding boxes
 *
 * \note This is an optional setting. By default a relative tolerance equal to 0 is used.
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   tol             Tolerance
 *
 */

void
PDM_mesh_location_tolerance_set
(
       PDM_mesh_location_t *ml,
 const double               tol
);


/**
 *
 * \brief Set the method for computing location (preconditioning stage)
 *
 * \note This is an optional setting.
 *
 * Admissible values are :
 *    - \p PDM_MESH_LOCATION_OCTREE         : Use point octree (default method)
 *    - \p PDM_MESH_LOCATION_DBBTREE        : Use bounding-box tree
 *    - \p PDM_MESH_LOCATION_LOCATE_ALL_TGT : All target points are guaranteed to be located
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   method          Preconditioning method
 *
 */

void
PDM_mesh_location_method_set
(
       PDM_mesh_location_t        *ml,
 const PDM_mesh_location_method_t  method
);


/**
 *
 * \brief Compute point location
 *
 * \param [in]   ml              Pointer to \ref PDM_mesh_location_t instance
 *
 */

void
PDM_mesh_location_compute
(
 PDM_mesh_location_t        *ml
);


/**
 *
 * \brief Get point location
 *
 * \note The results are related to located points only
 *
 * \param [in]   ml                    Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud         Point cloud identifier
 * \param [in]   i_part                Partition identifier
 * \param [out]  location              Global id of nearest mesh element if the point is located, -1 otherwise
 * \param [out]  dist2                 Signed squared distance from nearest element (negative if the point is located inside that element)
 * \param [out]  projected_coord       Cartesian coordinates of projection onto the nearest element (identity if the point is located inside that element)
 *
 */

void
PDM_mesh_location_point_location_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       PDM_g_num_t         **location,
       double              **dist2,
       double              **projected_coord
);


/**
 *
 * \brief Get the cell->vertex connectivity used for internal computations
 *
 * \note For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate
 *       the `points_weights` array (returned by \ref PDM_mesh_location_points_in_elt_get)
 *       to the appropriate mesh vertices.
 *
 * \param [in]   ml                    Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_part                Partition identifier
 * \param [out]  cell_vtx_idx          Index for cell -> vertex connectivity (size = *n_elt* + 1)
 * \param [out]  cell_vtx              Cell -> vertex connectivity
 *
 */

void
PDM_mesh_location_cell_vertex_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
       int                 **cell_vtx_idx,
       int                 **cell_vtx
);


/**
 *
 * \brief Get location data for points located in elements
 *
 * \param [in]   ml                      Pointer to \ref PDM_mesh_location_t instance
 * \param [in]   i_point_cloud           Point cloud identifier
 * \param [in]   i_part                  Partition identifier
 * \param [out]  elt_pts_inside_idx      Index for element -> points mapping (size = *n_elt* + 1)
 * \param [out]  points_gnum             Located points global ids (size : \p elt_pts_inside_idx[ *n_elt* ])
 * \param [out]  points_coords           Located points cartesian coordinates (size : 3 * \p elt_pts_inside_idx[ *n_elt* ])
 * \param [out]  points_uvw              Located points parametric coordinates (size : 3 * \p elt_pts_inside_idx[ *n_elt* ])
 * \param [out]  points_weights_idx      Index for interpolation weights (size : \p elt_pts_inside_idx[ *n_elt* ] + 1)
 * \param [out]  points_weights          Interpolation weights (size : \p points_weights_idx[\p elt_pts_inside_idx[ *n_elt* ]])
 * \param [out]  points_dist2            Signed squared distance element-points (< 0 if the point is inside) (size : \p elt_pts_inside_idx[ *n_elt* ])
 * \param [out]  points_projected_coords Cartesian coordinates of projection on element (identity if the point is inside) (size : 3 * \p elt_pts_inside_idx[ *n_elt* ])
 *
 */

void
PDM_mesh_location_points_in_elt_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       int                 **elt_pts_inside_idx,
       PDM_g_num_t         **points_gnum,
       double              **points_coords,
       double              **points_uvw,
       int                 **points_weights_idx,
       double              **points_weights,
       double              **points_dist2,
       double              **points_projected_coords
);


/**
 *
 * \brief Free a mesh location structure
 *
 * \param [in]  ml       Pointer to \ref PDM_mesh_location_t instance
 *
 */

void
PDM_mesh_location_free
(
 PDM_mesh_location_t  *ml
);


/**
 *
 * \brief Get the number of cells
 *
 * \param [in]  ml       Pointer to \ref PDM_mesh_location_t instance
 * \param [in]  i_part   Partition identifier
 *
 * \return Number of cells
 */

int
PDM_mesh_location_n_cell_get
(
       PDM_mesh_location_t *ml,
 const int                  i_part
);


/**
 *
 * \brief  Dump elapsed and CPU times
 *
 * \param [in]  ml       Pointer to \ref PDM_mesh_location_t instance
 *
 */

void
PDM_mesh_location_dump_times
(
PDM_mesh_location_t *ml
);


/**
 *
 * \brief Get the nodal mesh
 *
 * \param [in]  ml  Pointer to \ref PDM_mesh_location_t instance
 *
 * \return          Pointer to \ref PDM_part_mesh_nodal_t instance
 *
 */

PDM_part_mesh_nodal_t*
PDM_mesh_location_mesh_nodal_get
(
PDM_mesh_location_t *ml
);


/**
 * \brief Get \ref PDM_part_to_part_t instance to exchange data between
 * the source mesh and a target point cloud (both in user frame)
 *
 * \param [in ] ml            Pointer to \ref PDM_mesh_location_t instance
 * \param [in ] i_point_cloud Point cloud identifier
 * \param [out] ptp           Pointer to \ref PDM_part_to_part_t instance
 * \param [in ] ownership     Ownership for \p ptp
 *
 */

void
PDM_mesh_location_part_to_part_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership
 );

#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_LOCATION_H
