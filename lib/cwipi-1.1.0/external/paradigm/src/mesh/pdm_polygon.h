/*
 * \file
 */

#ifndef __PDM_POLYGON_H__
#define __PDM_POLYGON_H__

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

/**
 * \enum PDM_polygon_status_t
 * \brief Polygon status type
 *
 */

typedef enum {

  PDM_POLYGON_INSIDE      = 0,  /*!< Inside  */
  PDM_POLYGON_OUTSIDE     = 1,  /*!< Outside */
  PDM_POLYGON_DEGENERATED = 2,  /*!< Degenerated */

} PDM_polygon_status_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Get bounds
 *
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 *
 * \return      Bounds
 *
 */


double *
PDM_polygon_bounds_get
(
 const int     numPts,
 const double *pts
);


/**
 * \brief Evaluates the position in a polygon
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 * \param [out] closest  Closest Point in Polygon or NULL
 * \param [out] minDist2 Square of the distance
 *
 * \return      \ref PDM_POLYGON_INSIDE or \ref PDM_POLYGON_OUTSIDE
 *              if the projected is in the polygon or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_evaluate_position
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double        closestPoint[3],
 double        *minDist2
);


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref PDM_TRUE except for a triangle
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_bool_t
PDM_polygon_parameterize
(
 const int     numPts,
 const double *pts,
 double       *p0,
 double       *p10,
 double       *l10,
 double       *p20,
 double       *l20,
 double       *n
);


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_point_in
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double       *bounds,
 double       *n
);

PDM_polygon_status_t
PDM_polygon_point_in_new
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 double        bounds[6],
 double        normal[3]
);
/**
 * \brief Computes polygon barycenter
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  bary    Barycenter
 *
 *
 */

void
PDM_polygon_compute_barycenter
(
 const int numPts,
 const double *pts,
 double bary[3]
);


/**
 * \brief Test if a point is inside a 2d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xy            Point (x,y)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xy        Polygon vertices (x,y)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */


PDM_polygon_status_t PDM_polygon_point_in_2d
(
 const double  xy[2],
 const int     n_vtx,
 const double *vtx_xy,
 // const double  char_length,
 double        bounds[4]
 );

/**
 * \brief Test if a point is inside a 3d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xyz           Point (x,y,z)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xyz       Polygon vertices (x,y,z)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax, zmin, zmax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

PDM_polygon_status_t PDM_polygon_point_in_3d
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 // const double  char_length,
 double        bounds[6],
 double        normal[3]
 );


/**
 * \brief Compute parametric coordinates of a polygon's vertices
 * and a set of point inside the polygon's median plane.
 *
 * The uv-coordinate system is defined as follows:
 *   - origin at first polygon vertex ;
 *   - u-axis oriented from first to second polygon vertex.
 * Aspect-ratio and scale is preserved by the projection (no normalization).
 *
 * \param [in]  n_vtx    Number of polygon vertices
 * \param [in]  vtx_xyz  xyz-coordinates of polygon vertices (size = 3 * \ref n_vtx)
 * \param [out] vtx_uv   uv-coordinates of polygon vertices (size = 2 * \ref n_vtx)
 * \param [in]  n_pts    Number of points
 * \param [in]  pts_xyz  xyz-coordinates of points (size = 3 * \ref n_pts)
 * \param [out] pts_uv   uv-coordinates of points (size = 2 * \ref n_pts)
 * \param [in]  normal   Polygon's normal (size = 3 or NULL)
 *
 * \return      0 if the polygon is degenerate, 1 else.
 *
 */

int PDM_polygon_3d_to_2d
(
 const int    n_vtx,
 const double vtx_xyz[],
 double       vtx_uv[],
 const int    n_pts,
 const double pts_xyz[],
 double       pts_uv[],
 double       normal[3]
 );


/**
 * \brief Compute intersection point between a polygon and a semi-infinite ray
 *
 * \param[in]  ray_origin        Ray origin
 * \param[in]  ray_direction     Ray direction (need not be normalized)
 * \param[in]  n_vtx             Number of vertices
 * \param[in]  vtx_coord         Coordinates of the polygon's vertices
 * \param[in]  poly_center       Coordinates of the polygon's centroid (or NULL)
 * \param[in]  poly_normal       Normal to polygon's median plane (or NULL)
 * \param[in]  poly_bound        Polygon's bounds ([xmin, xmax, ymin, ymax, zmin, zmax] or NULL)
 * \param[out] intersection      Coordinates of the intersection point
 * \param[out] t                 Ray-parameter of the intersection point
 * \param[out] weight            Barycentric coordinates in polygon of intersection point (or NULL)
 *
 * \return Intersection status
 *
 */

PDM_polygon_status_t
PDM_polygon_ray_intersection
(
 const double  ray_origin[3],
 const double  ray_direction[3],
 const int     n_vtx,
 const double *vtx_coord,
       double *poly_center,
       double *poly_normal,
       double *poly_bound,
       double  intersection[3],
       double *t,
       double *weight
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
