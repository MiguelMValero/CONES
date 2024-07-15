/*
 * \file
 */

#ifndef __PDM_INTERPOLATE_FROM_MESH_LOCATION_H__
#define __PDM_INTERPOLATE_FROM_MESH_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_interpolate_from_mesh_location_t PDM_interpolate_from_mesh_location_t;

typedef enum {

  PDM_INTERPOLATE_KIND_FROM_CENTER = 0,
  PDM_INTERPOLATE_KIND_FROM_NODE   = 1

} PDM_interpolate_kind_t;



/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute interpolation from mesh_location information
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

PDM_interpolate_from_mesh_location_t*
PDM_interpolate_from_mesh_location_create
(
 const int                    n_cloud_target,
       PDM_interpolate_kind_t interp_kind,
 const PDM_MPI_Comm           comm
);

/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_interpolate_from_mesh_location_n_part_cloud_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_point_cloud,
 const int                                   n_part
);


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_interpolate_from_mesh_location_mesh_global_data_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   n_part
);

/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_interpolate_from_mesh_location_cloud_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_point_cloud,
 const int                                   i_part,
 const int                                   n_points,
       double                               *coords,
       PDM_g_num_t                          *gnum
);

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */
void
PDM_interpolate_from_mesh_location_part_set
(
       PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 const int                                    i_part,
 const int                                    n_cell,
 const int                                   *cell_face_idx,
 const int                                   *cell_face,
 const PDM_g_num_t                           *cell_ln_to_gn,
 const int                                    n_face,
 const int                                   *face_vtx_idx,
 const int                                   *face_vtx,
 const PDM_g_num_t                           *face_ln_to_gn,
 const int                                    n_vtx,
 const double                                *coords,
 const PDM_g_num_t                           *vtx_ln_to_gn
);

/**
 *
 * \brief Set point list located in elements
 *
 * \param [in]   id                      Identifier
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [in]   elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [in]   points_gnum             Points global number
 * \param [in]   points_coords           Points coordinates
 * \param [in]   points_uvw              Points parametric coordinates in elements
 * \param [in]   points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [in]   points_weights          Interpolation weights
 * \param [in]   points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [in]   points_projected_coords Point projection on element if the point is outside
 *
 */
void
PDM_interpolate_from_mesh_location_points_in_elt_set
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 const int                              i_part,
 const int                              i_point_cloud,
 const int                              n_elts,
 int                                   *elt_pts_inside_idx,
 PDM_g_num_t                           *points_gnum,
 double                                *points_coords,
 double                                *points_uvw,
 int                                   *points_weights_idx,
 double                                *points_weights,
 double                                *points_dist2,
 double                                *points_projected_coords
);

void
PDM_interpolate_from_mesh_location_compute
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
);

void
PDM_interpolate_from_mesh_location_exch
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 int                                     i_point_cloud,
 size_t                                  s_data,
 double                                **part_data_in,
 double                               ***cloud_data_out
);

void
PDM_interpolate_from_mesh_location_exch_inplace
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 int                                     i_point_cloud,
 size_t                                  s_data,
 double                                **part_data_in,
 double                                **cloud_data_out
);


void
PDM_interpolate_from_mesh_location_send
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 size_t                                  s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
);

void
PDM_interpolate_from_mesh_location_recv
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 size_t                                  s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
);

void
PDM_interpolate_from_mesh_location_free
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_INTERPOLATE_FROM_MESH_LOCATION_H__ */

