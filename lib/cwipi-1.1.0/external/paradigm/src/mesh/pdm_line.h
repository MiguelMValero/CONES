/*
 * \file
 */

#ifndef __PDM_LINE_H__
#define __PDM_LINE_H__

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
 * \enum PDM_line_status_t
 * \brief Polygon status type
 *
 */

typedef enum {

  PDM_LINE_INTERSECT_UNDEF   = -1,  /*!< No intersection */
  PDM_LINE_INTERSECT_NO      = 0,  /*!< No intersection */
  PDM_LINE_INTERSECT_YES     = 1,  /*!< Intersection */
  PDM_LINE_INTERSECT_ON_LINE = 2,  /*!< On line  */

} PDM_line_intersect_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Performs intersection of two finite 3D lines
 *
 *  An intersection is found if the projection of the two lines onto the plane
 *  perpendicular to the cross product of the two lines intersect.
 *  The parameters (u,v) are the parametric coordinates of the lines at the
 *  position of closest approach
 *  This function is more robust than
 *  PDM_line_intersection_mean_square.
 *
 * \param [in]  a1 Coordinates of the first line vertex of 'a'
 * \param [in]  a2 Coordinates of the second line vertex of 'a'
 * \param [in]  b1 Coordinates of the first line vertex of 'b'
 * \param [in]  b2 Coordinates of the second line vertex of 'b'
 * \param [out] u  Parameter of the intersection in line 'a' parametric coordinates
 * \param [out] v  Parameter of the intersection in line 'b' parametric coordinates
 *
 * \return      \ref PDM_TRUE or \ref PDM_FALSE
 *
 */

PDM_line_intersect_t
PDM_line_intersection
(
 const double a1[3],
 const double a2[3],
 const double b1[3],
 const double b2[3],
 double *u,
 double *v
 );


/**
 * \brief Performs intersection of two finite lines
 *
 *  Performs the intersection of two lines with mean square method
 *
 * \param [in]  a1 Coordinates of the first line vertex of 'a'
 * \param [in]  a2 Coordinates of the second line vertex of 'a'
 * \param [in]  b1 Coordinates of the first line vertex of 'b'
 * \param [in]  b2 Coordinates of the second line vertex of 'b'
 * \param [out] u  Parameter of the intersection in line 'a' parametric coordinates
 * \param [out] v  Parameter of the intersection in line 'b' parametric coordinates
 *
 * \return      \ref PDM_TRUE or \ref PDM_FALSE
 *
 */


PDM_line_intersect_t
PDM_line_intersection_mean_square
(
 const double a1[3],
 const double a2[3],
 const double b1[3],
 const double b2[3],
 double *u,
 double *v
 );

/**
 * \brief Computes point-line distance
 *
 * \param [in]  x             Point coordinates
 * \param [in]  a1            First line vertex coordinates
 * \param [in]  a2            Second line vertex coordinates
 * \param [out] t             Parameter of the intersection in line parametric coordinates
 * \param [out] closest_point Closest point
 *
 * \return   The square of the distance
 *
 */

double
PDM_line_distance
(
 const double x[3],
 const double p1[3],
 const double p2[3],
 double *t,
 double closestPoint[3]
 );


PDM_line_intersect_t
PDM_line_intersection_2d
(
 const double a1[2],
 const double a2[2],
 const double b1[2],
 const double b2[2],
 double *u,
 double *v
 );

double
PDM_line_distance_2d
(
 const double uv[2],
 const double p1[2],
 const double p2[2],
 double *t,
 double closest_point[2]
 );


/**
 * \brief Evaluates the position on an line segment
 *
 * \param [in]  x               Point coordinates to evaluate position
 * \param [in]  vtx_coord       Line segment vertices coordinates
 * \param [out] closest_point   Closest Point on Line segment or NULL
 * \param [out] min_dist2       Square of the distance
 * \param [out] weights         Vertices weights or NULL
 *
 * \return      -1 if the line segment is degenerate, 0 else
 *
 */

int PDM_line_evaluate_position
(
 const double  x[3],
 const double *vtx_coord,
 double       *closest_point,
 double       *dist2,
 double        weights[2]
 );




PDM_line_intersect_t
PDM_line_intersection_2drobust
(
 double  a[3],
 double  b[3],
 double  c[3],
 double  d[3],
 double *u,
 double *v
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
