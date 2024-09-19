/*
 * \file
 */

#ifndef __PDM_HO_BEZIER_H__
#define __PDM_HO_BEZIER_H__

/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/


/**
 *
 * \brief De Casteljau algorithm for Bézier curves
 *
 * Evaluates the point with parameter t on the Bézier curve
 * and (optionally) builds the control points of the two
 * subcurves sharing the evaluated point as a common vertex.
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   t       Parametric coordinate
 * \param[in]   b       Bézier control points
 * \param[out]  val     Evaluated point
 * \param[out]  l       Control points of the 1st subcurve
 * \param[out]  r       Control points of the 2nd subcurve
 *
 */

void
PDM_ho_bezier_de_casteljau_curve
(
 const int     dim,
 const int     order,
 const double  t,
 double       *b,
 double       *val,
 double       *l,
 double       *r
 );

/**
 *
 * \brief De Casteljau algorithm for Bézier triangles
 *
 * Evaluates the point (u,v) on the Bézier triangle
 * and (optionally) builds the control points of the three
 * subtriangles sharing the evaluated point as a common vertex.
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   u       Parametric coordinate u
 * \param[in]   v       Parametric coordinate v
 * \param[in]   b       Bézier control points
 * \param[out]  val     Evaluated point
 * \param[out]  atr     Control points of the 1st subtriangle
 * \param[out]  ars     Control points of the 2nd subtriangle
 * \param[out]  ast     Control points of the 3rd subtriangle
 *
 */

void
PDM_ho_bezier_de_casteljau_triangle
(
 const int     dim,
 const int     order,
 const double  u,
 const double  v,
 double       *b,
 double       *val,
 double       *atr,
 double       *ars,
 double       *ast
 );


/**
 *
 * \brief Build control points for derivative of a Bézier curve
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   b       Bézier control points
 * \param[out]  db_dt   Bézier control points of derivative
 *
 */

void
PDM_ho_bezier_curve_derivative
(
 const int     dim,
 const int     order,
       double *b,
       double *db_dt
 );


/**
 *
 * \brief Build control points for partial derivatives of a Bézier triangle
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   b       Bézier control points
 * \param[out]  bu      Bézier control points of 1st partial derivative
 * \param[out]  bv      Bézier control points of 2nd partial derivative
 *
 */

void
PDM_ho_bezier_triangle_derivatives
(
 const int  dim,
 const int  order,
 double    *b,
 double    *bu,
 double    *bv
 );


/**
 *
 * \brief Point location in a high-order Bézier curve
 *
 * \param[in]   order             Order
 * \param[in]   n_node            Number of nodes
 * \param[in]   node_coord        Coordinates of the Bézier control points (size = 3 * \ref n_node)
 * \param[in]   point_coord       Coordinates of the point to locate (size = 3)
 * \param[out]  projected_coords  Coordinates of the projection on the Bézier curve (size = 3)
 * \param[out]  uvw               Parametric coordinates of the projection on the Bézier curve
 *
 */

double
PDM_ho_bezier_curve_location
(
 const int     order,
 const int     n_node,
       double *node_coord,
       double *point_coord,
       double *projected_coord,
       double *u
 );


/**
 *
 * \brief Point location in a high-order Bézier triangle
 *
 * \param[in]   order             Order
 * \param[in]   n_node            Number of nodes
 * \param[in]   node_coord        Coordinates of the Bézier control points (size = 3 * \ref n_node)
 * \param[in]   point_coord       Coordinates of the point to locate (size = 3)
 * \param[out]  projected_coords  Coordinates of the projection on the Bézier triangle (size = 3)
 * \param[out]  uvw               Parametric coordinates of the projection on the Bézier triangle
 *
 */

double
PDM_ho_bezier_triangle_location
(
 const int     order,
 const int     n_node,
       double *node_coord,
       double *point_coord,
       double *projected_coord,
       double *uv
 );

#endif /* __PDM_HO_BEZIER_H__ */
