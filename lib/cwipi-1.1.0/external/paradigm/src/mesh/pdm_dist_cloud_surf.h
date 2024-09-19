/*
 * \file
 */

#ifndef PDM_DIST_CLOUD_SURF_H
#define PDM_DIST_CLOUD_SURF_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_surf_mesh.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dist_cloud_surf_t PDM_dist_cloud_surf_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute distance between point clouds and a surface mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point clouds
 * \param [in]   comm           MPI communicator
 * \param [in]   owner          Ownership of \ref PDM_dist_cloud_surf_t
 *
 * \return     Pointer to \ref PDM_dist_cloud_surf_t object
 */


PDM_dist_cloud_surf_t*
PDM_dist_cloud_surf_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_dist_cloud_surf_n_part_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    n_part
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates (size : 3 * \p n_points)
 * \param [in]   gnum            Point global ids (size : 3 * \p n_points)
 *
 */

void
PDM_dist_cloud_surf_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
 const int                    n_points,
       double                *coords,
       PDM_g_num_t           *gnum
);



/**
 *
 * \brief Set the nodal mesh
 *
 * \param [in]   dist        Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   mesh_nodal  Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */

void
PDM_dist_cloud_surf_nodal_mesh_set
(
 PDM_dist_cloud_surf_t *dist,
 PDM_part_mesh_nodal_t *mesh_nodal
);

/**
 *
 * \brief Map a surface mesh
 *
 * \param [in]   dist       Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   surf_mesh  Pointer to \ref PDM_surf_mesh_t object
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_map
(
 PDM_dist_cloud_surf_t *dist,
 PDM_surf_mesh_t       *surf_mesh
);


/**
 *
 * \brief Set the number of partitions of the mesh
 *
 * \param [in]   dist           Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   n_part         Number of partitions
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_global_data_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    n_part
);


/**
 *
 * \brief Set a surface mesh partition
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_part          Partition identifier
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
PDM_dist_cloud_surf_surf_mesh_part_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_part,
 const int                    n_face,
 const int                   *face_vtx_idx,
 const int                   *face_vtx,
 const PDM_g_num_t           *face_ln_to_gn,
 const int                    n_vtx,
 const double                *coords,
 const PDM_g_num_t           *vtx_ln_to_gn
);


/**
 *
 * \brief Compute distance
 *
 * \param [in]   dist  Pointer to \ref PDM_dist_cloud_surf_t object
 *
 */

void
PDM_dist_cloud_surf_compute
(
 PDM_dist_cloud_surf_t *dist
);


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   dist                  Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud         Point cloud identifier
 * \param [in]   i_part                Partition identifier
 * \param [out]  closest_elt_distance  Squared distance from nearest elements
 * \param [out]  closest_elt_projected Coordinates of projection onto the nearest elements
 * \param [out]  closest_elt_gnum      Global ids of the nearest elements
 *
 */

void
PDM_dist_cloud_surf_get
(
       PDM_dist_cloud_surf_t  *dist,
 const int                     i_point_cloud,
 const int                     i_part,
       double                **distance,
       double                **projected,
       PDM_g_num_t           **closest_elt_gnum
);


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf_t object
 *
 */

void
PDM_dist_cloud_surf_free
(
 PDM_dist_cloud_surf_t  *dist
);


/**
 *
 * \brief  Dump elapsed and CPU times
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf_t object
 *
 */

void
PDM_dist_cloud_surf_dump_times
(
 PDM_dist_cloud_surf_t  *dist
);



/**
 *
 * \brief Get the dimensions of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud   Point cloud identifier
 * \param [in]   i_part          Partition identifier
 * \param [out]  n_points        Number of points in current part of current cloud
 *
 */

void
PDM_dist_cloud_surf_cloud_dim_get
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
       int                   *n_points
);


/**
 *
 * \brief Get the number of partitions in a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud   Point cloud identifier
 *
 * \return   Number of partitions
 *
 */

int
PDM_dist_cloud_surf_cloud_n_part_get
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud
);



/**
 *
 * \brief Distribute data from the surface mesh to a point cloud
 *
 * \param [in]   dist              Pointer to \ref PDM_dist_cloud_surf_t object
 * \param [in]   i_point_cloud     Point cloud identifier
 * \param [in]   stride            Stride
 * \param [in]   surf_data         Surface mesh data (send)
 * \param [out]  cloud_data        Point cloud data  (recv)
 *
 */

void
PDM_dist_cloud_surf_distri_data
(
       PDM_dist_cloud_surf_t  *dist,
 const int                     i_point_cloud,
 const int                     stride,
 const void                  **surf_data,
       void                 ***cloud_data
);


#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_DIST_H
