#ifndef __FVMC_POINT_LOCATION_H__
#define __FVMC_POINT_LOCATION_H__

/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011-2018  ONERA

  Copyright (C) 2007-2009  EDF

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

/*----------------------------------------------------------------------------*/

#include "fvmc_config.h"
#include "fvmc_config_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define FVMC_POLYGON_FAILURE -1
#define FVMC_POLYGON_OUTSIDE 0
#define FVMC_POLYGON_INSIDE 1
#define FVMC_POLYGON_INTERSECTION 2
#define FVMC_POLYGON_ON_LINE 3

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_nodal        <-- pointer to nodal mesh representation structure
 *   tolerance         <-- associated tolerance
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   projected_coords  <-> coordinates of projected points in location elements
 *                         point (size: n_points * dim)
 *   uvw               <-> parametric coordinates of the point if inside the element
 *                         parametric coordinates of the projected point if outside the element
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_location_nodal(const fvmc_nodal_t  *this_nodal,
                          double                tolerance,
                          _Bool                 locate_on_parents,
                          fvmc_lnum_t           n_points,
                          const fvmc_coord_t    point_coords[],
                          fvmc_coord_t         *projected_coords,
                          double               *uvw,
                          fvmc_lnum_t           location[],
                          float                 distance[]);

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * This function currently only handles elements of lower dimension than
 * the spatial dimension.
 *
 * parameters:
 *   this_nodal        <-- pointer to nodal mesh representation structure
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_location_closest_nodal(const fvmc_nodal_t  *this_nodal,
                                 _Bool               locate_on_parents,
                                 fvmc_lnum_t          n_points,
                                 const fvmc_coord_t   point_coords[],
                                 fvmc_lnum_t          location[],
                                 float               distance[]);


/*----------------------------------------------------------------------------
 * Compute distance to polygons
 *
 * parameters:
 *   dim               <-- dimension
 *   n_poly            <-- number of polygon
 *   connectivity_idx  <-- polygon connectivity index
 *   connectivity      <-- polygon connectivity
 *   vertex_coords     <-- polygon connectivity
 *   n_points          <-- polygon connectivity
 *   point_coords      <-- polygon connectivity
 *   distance          --> 3d surf : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *                         2d or 3d : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_dist_closest_polygon(const int            dim,
                                const fvmc_lnum_t    n_poly,
                                const fvmc_lnum_t    connectivity_idx[],
                                const fvmc_lnum_t    connectivity[],
                                const fvmc_coord_t   vertex_coords[],
                                const fvmc_lnum_t    n_points,
                                const fvmc_lnum_t    point_ids[],
                                const fvmc_coord_t   point_coords[],
                                fvmc_lnum_t          location[],
                                float                distance[]);


int fvmc_point_in_polygon (double x[3],
                           int numPts,
                           double *pts,
                           double *bounds,
                           double *n);

double fvmc_distance_to_line(double x[3], double p1[3], double p2[3],
                             double *t, double closestPoint[3]);

double fvmc_distant_to_polygon (double x[3], int numPts, double *pts,
                                double bounds[6], double closest[3]);

int fvmc_parameterize_polygon(int numPts, double *pts, double *p0, double *p10, double *l10,
                              double *p20,double *l20, double *n);

int  fvmc_triangle_evaluate_Position (double x[3], double *pts, double* closestPoint,
                                      double closestPointpcoords[2],
                                      double *dist2,
                                      double closestPointweights[3]);

int fvmc_polygon_evaluate_Position(double x[3], int numPts, double *pts, double* closestPoint,
                                   double pcoords[3], double* minDist2);

int fvmc_edge_evaluate_Position (double x[3], double *pts, double* closestPoint,
                                 double closestPointpcoords[1], double* dist2,
                                 double closestPointweights[2]);

int  fvmc_tetrahedron_evaluate_Position (double x[3], double *pts,
                                         double* closestPoint,
                                         double closestPointpcoords[3],
                                         double *dist2,
                                         double closestPointweights[4]);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_POINT_LOCATION_H__ */
