/*
 * \file
 */

#ifndef __PDM_MEAN_VALUES_H__
#define __PDM_MEAN_VALUES_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Compute mean value coordinates of points in a polygon in 2d
 *
 * See "Mean value coordinates for arbitrary planar polygons", Kai Hormann, and Michael S. Floater. (2006).
 *
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    vtx_coord        xyz-coordinates of polygon vertices (size = 2 * \ref n_vtx)
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts_coord        xyz-coordinates of points to locate (size = 2 * \ref n_pts)
 * \param [out]   mean_value_coord Mean value coordinates of points to locate (size = \ref n_vtx * \ref n_pts)
 *
 */

void
PDM_mean_values_polygon_2d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );


/**
 * \brief Compute mean value coordinates of points in a polygon in 3d
 *
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    vtx_coord        xyz-coordinates of polygon vertices (size = 3 * \ref n_vtx)
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts_coord        xyz-coordinates of points to locate (size = 3 * \ref n_pts)
 * \param [out]   mean_value_coord Mean value coordinates of points to locate (size = \ref n_vtx * \ref n_pts)
 *
 */

void
PDM_mean_values_polygon_3d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );



/**
 * \brief Compute mean value coordinates of a point in a polyhedron
 *
 * See "Mean value coordinates for closed triangular meshes", T. Ju et al. (2005).
 *
 * \param [in]    n_vtx            Number of polyhedron vertices
 * \param [in]    vtx_coord        xyz-coordinates of polyhedron vertices (size = 3 * \ref n_vtx)
 * \param [in]    n_face           Number of polyhedron faces
 * \param [in]    face_vtx_idx     Index for face-vertex connectivity (size = \ref n_face + 1)
 * \param [in]    face_vtx         Face-vertex connectivity (size = \ref face_vtx_idx[\ref n_face])
 * \param [in]    face_orientation Face orientation (size = \ref n_face)
 * \param [in]    pt_coord         xyz-coordinates of point to locate (size = 3)
 * \param [out]   weights          Mean value coordinates of point to locate (size = \ref n_vtx)
 *
 */

void
PDM_mean_values_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            weights[]
 );





/*  UNUSED FUNCTIONS  */

void
PDM_mean_values_polygon
(
 const int         n_vtx,
 const PDM_l_num_t face_vtx[],
 const double      vtx_coord[],
 const double      pt_coord[],
 double            weights[]
 );

void
PDM_mean_value_coordinates_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEAN_VALUES_H__ */
