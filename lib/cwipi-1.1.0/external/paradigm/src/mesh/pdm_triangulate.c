/*============================================================================
 * Triangulation of a polygon
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2005-2007  EDF

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_triangulate.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * State of current triangulation (working structure)
 *----------------------------------------------------------------------------*/

struct _PDM_triangulate_state_t {

  int          *triangle_vertices; /* current triangle vertices list */
  double       *coords;            /* vertex coordinates */
  int          *list_previous;     /* indices of previous vertices in polygon;
                                      size:n_vertices; */
  int          *list_next;         /* indices of next vertices in polygon;
                                      size:n_vertices; */
  int          *edge_vertices;     /* edges connectivity; size: n_edges * 2 */
  int          *edge_neighbors;    /* triangles sharing a given edge
                                      (2 values per edge);
                                      - (-1,-1): edge does not exist
                                      - (x ,-1): boundary edge, triangle x
                                      - (x ,-1): internal edge, triangles x
                                      and y */
  PDM_bool_t   *edge_is_delaunay;  /* Delaunay edge indicator */
  PDM_bool_t   *concave;           /* Concave vertex indicator */

  int           n_vertices_max;    /* Maximum number vertices; a larger
                                      size requires resizing work arrays */

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Project a face in 3d on a plane parallel to this face, which is then
 * treated as the (Oxy) plane.
 *
 * parameters:
 *   n_vertices  <-- number of vertices defining the polygon.
 *   coords      <-> coordinates of the polygon's vertices (3d in, 2d out)
 *----------------------------------------------------------------------------*/

static void
_polygon_plane_3d(const int    n_vertices,
                  double       coords[])
{
  int i, j;

  double  face_center[3], face_normal[3];
  double  v1[3], v2[3], cross_v1_v2[3];
  double  dot_v1_v2, module_v2;
  double  tmp_coord;

  double   cost;
  double   sint;

#define _N_VERTICES_AUTO_MAX   20 /* Size of local temporary coordinates
                                     buffer; above this size, allocation
                                     is necessary */

  double  _tmp_coords[_N_VERTICES_AUTO_MAX * 3];
  double  *_tmp_coords_p = NULL;
  double  *tmp_coords = _tmp_coords;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Estimate position of polygon center */

  for (i = 0; i < 3; i++)
    face_center[i] = 0.;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < n_vertices; j++)
      face_center[i] += coords[j*3 + i];
    face_center[i] /= n_vertices;
  }

  /* Estimate face normal */

  for (i = 0; i < 3; i++)
    face_normal[i] = 0.;

  for (j = 0; j < n_vertices; j++) {

    for (i = 0; i < 3; i++) {
      v1[i] = coords[j*3 + i] - face_center[i];
      if (j < n_vertices - 1)
        v2[i] = coords[(j+1)*3 + i] - face_center[i];
      else
        v2[i] = coords[          i] - face_center[i];
    }

    face_normal[0] += v1[1]*v2[2] - v1[2]*v2[1];
    face_normal[1] += v1[2]*v2[0] - v1[0]*v2[2];
    face_normal[2] += v1[0]*v2[1] - v1[1]*v2[0];

  }

  /* Project coordinates in a plane parallel to face */
  /*-------------------------------------------------*/

  /* We place the coordinate system origin at the estimated face center */

  for (j = 0; j < n_vertices; j++)
    for (i = 0; i < 3; i++)
      coords[j*3 + i] -= face_center[i];

  if (PDM_ABS(face_normal[0]) > 1.e-12 || PDM_ABS(face_normal[1]) > 1.e-12) {

    /* First rotation of axis (Oz) and angle (Ox, normal proj. on Oxy) */

    if (n_vertices > _N_VERTICES_AUTO_MAX) {
      _tmp_coords_p  = malloc (sizeof(double) * n_vertices*3);
      tmp_coords = _tmp_coords_p;
    }

    v1[0] = 1.;
    v1[1] = 0.;
    v1[2] = 0.;

    v2[0] = face_normal[0];
    v2[1] = face_normal[1];
    v2[2] = 0.;

    PDM_CROSS_PRODUCT(cross_v1_v2, v1, v2);

    dot_v1_v2 = PDM_DOT_PRODUCT(v1, v2);
    module_v2 = PDM_MODULE(v2);
    cost = dot_v1_v2 / module_v2;

    if (cross_v1_v2[2] > 0.)
      sint =   PDM_MODULE(cross_v1_v2) / module_v2;
    else
      sint = - PDM_MODULE(cross_v1_v2) / module_v2;

    for (j = 0; j < n_vertices; j++) {

      tmp_coords[j*3    ] =  cost*coords[j*3] + sint*coords[j*3 + 1];
      tmp_coords[j*3 + 1] = -sint*coords[j*3] + cost*coords[j*3 + 1];
      tmp_coords[j*3 + 2] =  coords[j*3 +2];

    }

    /* Second rotation with axis (Oy) and angle (Oz', normal proj. on  Ox'z) */

    v1[0] =  0.;
    v1[1] =  0.;
    v1[2] =  1.;

    v2[0] = sqrt(face_normal[0]*face_normal[0] + face_normal[1]*face_normal[1]);
    v2[1] = 0.;
    v2[2] = face_normal[2];

    PDM_CROSS_PRODUCT(cross_v1_v2, v1, v2);

    cost = PDM_DOT_PRODUCT(v1, v2) / PDM_MODULE(v2);

    if (cross_v1_v2[2] > 0.)
      sint =  PDM_MODULE(cross_v1_v2) / PDM_MODULE(v2);
    else
      sint = -PDM_MODULE(cross_v1_v2) / PDM_MODULE(v2);

    for (j = 0; j < n_vertices; j++) {

      coords[j*3    ] = cost*tmp_coords[j*3] + sint*tmp_coords[j*3 + 2];
      coords[j*3 + 1] = tmp_coords[j*3 + 1];
      coords[j*3 + 2] = 0.;

    }

    if (_tmp_coords_p != NULL) {
      free (_tmp_coords_p);
      tmp_coords = NULL;
    }

  }
  else {

    /* We only need to set the vertices z coordinate to 0, possibly
       swapping the coordinates in the (Oxy) projection plane.  */

    if (face_normal[2] > 0.)
      for (j = 0; j < n_vertices; j++)
        coords[j*3 + 2] = 0.;

    else
      for (j = 0; j < n_vertices; j++) {
        tmp_coord = coords[j*3];
        coords[j*3    ] = coords[j*3 + 1];
        coords[j*3 + 1] = tmp_coord;
        coords[j*3 + 2] = 0.;
      }

  }

  /* Convert coords to 2d */

  for (j = 0; j < n_vertices; j++) {
    coords[j*2    ] = coords[j*3    ];
    coords[j*2 + 1] = coords[j*3 + 1];
  }

#undef _N_VERTICES_AUTO_MAX
}

/*----------------------------------------------------------------------------
 * Check if a (2D projection of a) polygon vertex is convex.
 *
 * parameters:
 *   previous    <-- index of previous vertex in polygon.
 *   current     <-- index of current vertex in polygon.
 *   next        <-- index of following vertex in polygon.
 *   coords      <-- coordinates of the polygon's vertices (2d).
 *----------------------------------------------------------------------------*/

static PDM_bool_t
_polygon_vertex_is_convex(const int     previous,
                          const int     current,
                          const int     next,
                          const double  coords[])
{
  /* sin(theta) = (v1 x v2) / (|v1| x |v2|), so the sign of the sine
     is that of the cross product's z component (as v1 and v2 lie in
     the same plane) */

  if (  (  (coords[current*2]     - coords[previous*2])
           * (coords[next*2 + 1]    - coords[current*2 + 1]))
        - (  (coords[next*2]        - coords[current*2])
             * (coords[current*2 + 1] - coords[previous*2 + 1])) > 0.0)
    return PDM_TRUE;

  else
    return PDM_FALSE;
}

/*----------------------------------------------------------------------------
 * Check if a (2D projection of a) polygon vertex is an ear.
 *
 * parameters:
 *   n_vertices    <-- number of vertices defining the polygon.
 *   current       <-- index of current vertex in polygon.
 *   list_previous <-- indices of previous vertices in polygon.
 *   list_next     <-- index of previous vertex in polygon.
 *   concave       <-- flag concave vertices in polygon.
 *   coords        <-- coordinates of the polygon's vertices (2d).
 *   epsilon       <-- associated relative tolerance.
 *----------------------------------------------------------------------------*/

static PDM_bool_t
_polygon_vertex_is_ear(const int           n_vertices,
                       const int           current,
                       const int           list_previous[],
                       const int           list_next[],
                       const PDM_bool_t    concave[],
                       const double        coords[],
                       const double        epsilon)

{
  int i, previous, next;
  double surf_2, x_iso, y_iso;
  double vect1[2], vect2[2], vect3[2];

  /* If no vertex is concave, we have an ear */

  for (i = 0; i < n_vertices && concave[i] == PDM_FALSE; i++);

  if (i == n_vertices)
    return PDM_TRUE;

  /* If current vertex is convex */

  else {

    if (concave[current] == PDM_FALSE) {

      /* We check if the triangle formed by
         list_previous[current], current, and list_next[current]
         contains a concave vertex */

      previous = list_previous[current];
      next     = list_next[current];

      vect2[0] = coords[current*2    ] - coords[previous*2    ];
      vect2[1] = coords[current*2 + 1] - coords[previous*2 + 1];
      vect3[0] = coords[next*2    ]    - coords[previous*2    ];
      vect3[1] = coords[next*2 + 1]    - coords[previous*2 + 1];

      surf_2 = vect2[0]*vect3[1] - vect3[0]*vect2[1];

      /*
       * FIXME : Fonction _polygon_vertex_is_ear, add a temporary constant epsilon (1e-32)
       *         to protect jacobian computation. It must be variable
       */

      const double eps_surf2 = 1e-32;

      if (surf_2 <= eps_surf2) {
        return PDM_TRUE;
      }

      for (i = list_next[next]; i != previous; i = list_next[i]) {

        if (concave[i] == PDM_TRUE) {

          vect1[0] = coords[i*2    ] - coords[previous*2    ];
          vect1[1] = coords[i*2 + 1] - coords[previous*2 + 1];

          x_iso = (vect1[0]*vect3[1] - vect1[1]*vect3[0]) / surf_2;
          y_iso = (vect2[0]*vect1[1] - vect2[1]*vect1[0]) / surf_2;

          if (   (1.0 - x_iso - y_iso > - epsilon)
                 && (      x_iso         > - epsilon)
                 && (              y_iso > - epsilon))
            return PDM_FALSE;

        }

      }

      return PDM_TRUE;

    }
    else
      return PDM_FALSE;

  }

}

/*----------------------------------------------------------------------------
 * Check if an edge (between two triangles) is locally Delaunay.
 *
 * We compute the power of a point compared to the circle circumscribed
 * to the triangle formed by the three others (of which two define
 * the diagonal considered).
 *
 * parameters:
 *   edge_vertex_0     <-- index of first edge vertex in coords.
 *   edge_vertex_1     <-- index of second edge vertex in coords.
 *   flip_vertex_0     <-- index of first flip edge vertex in coords.
 *   flip_vertex_1     <-- index of second flip edge vertex in coords.
 *   coords            <-- coordinates of the triangulation's vertices (2d).
 *----------------------------------------------------------------------------*/

static PDM_bool_t
_edge_is_locally_delaunay(const int          edge_vertex_0,
                          const int          edge_vertex_1,
                          const int          flip_vertex_0,
                          const int          flip_vertex_1,
                          const double       coords[])
{
  double   a, b, delta;
  double   lambda[4];
  double   x_center, y_center, radius;
  double   point_power;

  lambda[0] = 2*(coords[edge_vertex_1*2    ] - coords[edge_vertex_0*2    ]);
  lambda[1] = 2*(coords[edge_vertex_1*2 + 1] - coords[edge_vertex_0*2 + 1]);

  lambda[2] = 2*(coords[flip_vertex_0*2    ] - coords[edge_vertex_0*2    ]);
  lambda[3] = 2*(coords[flip_vertex_0*2 + 1] - coords[edge_vertex_0*2 + 1]);

  delta = lambda[1]*lambda[2] - lambda[0]*lambda[3];

  /* If the triangle is flat, we automatically switch diagonals to
     avoid a division by zero. */

  if (PDM_ABS(delta) < 1.e-12)
    return PDM_FALSE;

  a =   coords[edge_vertex_1*2    ] * coords[edge_vertex_1*2    ]
    - coords[edge_vertex_0*2    ] * coords[edge_vertex_0*2    ]
    + coords[edge_vertex_1*2 + 1] * coords[edge_vertex_1*2 + 1]
    - coords[edge_vertex_0*2 + 1] * coords[edge_vertex_0*2 + 1];

  b =   coords[flip_vertex_0*2    ] * coords[flip_vertex_0*2    ]
    - coords[edge_vertex_0*2    ] * coords[edge_vertex_0*2    ]
    + coords[flip_vertex_0*2 + 1] * coords[flip_vertex_0*2 + 1]
    - coords[edge_vertex_0*2 + 1] * coords[edge_vertex_0*2 + 1];

  /* Center and radius of the circle passing through the vertices of
     the diagonal and through a third vertex. */

  x_center = (lambda[1]*b - lambda[3]*a)/delta;
  y_center = (lambda[2]*a - lambda[0]*b)/delta;

  radius = sqrt( (  (x_center - coords[edge_vertex_0*2    ])
                    * (x_center - coords[edge_vertex_0*2    ]))
                 + (  (y_center - coords[edge_vertex_0*2 + 1])
                      * (y_center - coords[edge_vertex_0*2 + 1])));

  /* Compute the power of the quadrilateral's fourth vertex compared
     to the circle. */

  point_power =   (  (coords[flip_vertex_1*2    ] - x_center)
                     * (coords[flip_vertex_1*2    ] - x_center))
    + (  (coords[flip_vertex_1*2 + 1] - y_center)
         * (coords[flip_vertex_1*2 + 1] - y_center))
    - radius*radius;

  /* Keep a slight margin in case there is no perceptible gain
     in switching diagonals */

  if (point_power > -1.e-12)
    return PDM_TRUE;
  else
    return PDM_FALSE;
}

/*----------------------------------------------------------------------------
 * Sort vertex indexes defining a triangle.
 *
 * parameters:
 *   triangle_vertices <-> triangles connectivity.
 *----------------------------------------------------------------------------*/

static void
_triangle_by_sorted_vertices(int triangle_vertices[])
{
  int  i_min, i_max, i_mid;

  /* First step */

  if (triangle_vertices[0] < triangle_vertices[1]) {
    i_min = triangle_vertices[0];
    i_mid = triangle_vertices[1];
    i_max = i_mid;
  }
  else {
    i_min = triangle_vertices[1];
    i_mid = triangle_vertices[0];
    i_max = i_mid;
  }

  /* Second step */

  if (triangle_vertices[2] < i_min) {
    i_mid = i_min;
    i_min = triangle_vertices[2];
  }
  else if (triangle_vertices[2] > i_max) {
    i_max = triangle_vertices[2];
  }
  else {
    i_mid = triangle_vertices[2];
  }

  /* Reordering */

  triangle_vertices[0] = i_min;
  triangle_vertices[1] = i_mid;
  triangle_vertices[2] = i_max;
}

/*----------------------------------------------------------------------------
 * Convert a given triangulation to a Delaunay triangulation using
 * edge flips.
 *
 * parameters:
 *   n_vertices        <-- number of vertices defining the triangulation.
 *   triangle_vertices <-> triangles connectivity.
 *   edge_vertices     --> edges connectivity.
 *   edge_neigbors     --> triangles sharing a given edge (2 values per edge);
 *                          - (-1,-1): edge does not exist
 *                          - (x ,-1): boundary edge, triangle x
 *                          - (x ,-1): internal edge, triangles x and y
 *   edge_is_delaunay  --> delaunay edge indicator.
 *   coords            <-- coordinates of the triangulation's vertices (2d).
 *----------------------------------------------------------------------------*/

static int
_polygon_delaunay_flip(const int           n_vertices,
                       int                 triangle_vertices[],
                       int                 edge_vertices[],
                       int                 edge_neighbors[],
                       PDM_bool_t          edge_is_delaunay[],
                       const double        coords[])
{
  int is_delaunay = 1;

  int    triangle_id, edge_id, vertex_id;
  int    triangle_id_0, triangle_id_1;

  int    current_edge_id, flip_edge_id;

  int    vertex_flip[2];

  int    i, i_0, i_1, i_min, i_max, j;

  PDM_bool_t   face_is_delaunay, edge_locally_delaunay;
  PDM_bool_t   restart, convex_quad;

  const int  n_edges = n_vertices*(n_vertices - 1)/2;
  const int  n_triangles = n_vertices - 2;
  int        i_previous[2] = {-1, -1}, i_next[2] = {-1, -1};

  /*
    There are (n_vertices*(n_vertices - 1)/2) possible edge combinations.
    An edge's number can be given by:
    -> SUM(k = 0, k < i_min) {n_vertices - k} + (i_max - i_min)
    i.e. n_vertices*i_min - i_min*(i_min+1)/2 + (i_max - i_min)
    An edge's index equals it's number - 1.

    Be careful with the order of arguments! arg[0]=i_min, arg[1]=i_max
  */

#undef _EDGE_INDEX
#define _EDGE_INDEX(i_min, i_max)                           \
  (n_vertices*i_min - i_min*(i_min+1)/2 + i_max-i_min - 1)

  /* Initialization */
  /*----------------*/

  for (i_0 = 0; i_0 < n_vertices; i_0++) {
    for (i_1 = i_0 + 1; i_1 < n_vertices; i_1++) {

      edge_id = _EDGE_INDEX(i_0, i_1);

      /* Define edges */

      edge_vertices[2*edge_id    ] = i_0;
      edge_vertices[2*edge_id + 1] = i_1;

      /*
        Liste of triangles sharing an edge:
        - (-1,-1): edge does not exist
        - (x ,-1): boundary edge, triangle x
        - (x ,-1): internal edge, triangles x and y
      */

      edge_neighbors[2*edge_id    ] = -1;
      edge_neighbors[2*edge_id + 1] = -1;

      /* Initialize an array indicating if an edge is locally Delaunay */

      edge_is_delaunay[edge_id] = PDM_TRUE;

    }
  }

  /* First traversal of triangles to determine initial neighborhood,
     as well as the list of edges which are not locally Delaunay */

  for (j = 0; j < n_triangles; j++) {
    for (i = 0; i < 3; i++) {

      i_0 = triangle_vertices[(j*3) +   i      ];
      i_1 = triangle_vertices[(j*3) + ((i+1)%3)];

      i_min = PDM_MIN(i_0, i_1);
      i_max = PDM_MAX(i_0, i_1);

      edge_id = _EDGE_INDEX(i_min, i_max);

      /* Update edge neighbors */

      if (edge_neighbors[2*edge_id] == -1)
        edge_neighbors[2*edge_id] = j;
      else
        edge_neighbors[2*edge_id + 1] = j;

      /* If edge is not on the boundary: */

      if (   !(i_max == i_min + 1)
             && !(i_min == 0 && i_max == n_vertices - 1))
        edge_is_delaunay[edge_id] = PDM_FALSE;

    }
  }

  /* Main flip algorithm */
  /*---------------------*/

  edge_id = 0;

  restart = PDM_FALSE;
  face_is_delaunay = PDM_FALSE;

  int n_restart = 0;
  const int max_n_restart = 1000;

  while(face_is_delaunay == PDM_FALSE) {

    if (edge_is_delaunay[edge_id] == PDM_FALSE) {

      edge_is_delaunay[edge_id] = PDM_TRUE;

      i_0 = edge_vertices[2*edge_id];
      i_1 = edge_vertices[2*edge_id + 1];

      for (j = 0; j < 2; j++) { /* Loop on triangles on each side of edge */

        triangle_id = edge_neighbors[2*edge_id + j];

        for (i = 0; i < 3; i++) { /* Loop on triangle's vertices */

          vertex_id = triangle_vertices[3*triangle_id + i];

          /* Seek opposite vertices */

          if (vertex_id != i_0 && vertex_id != i_1)
            vertex_flip[j] = vertex_id;

          /* Seek preceding and following vertices */

          if (   vertex_id == i_0
                 && triangle_vertices[3*triangle_id + ((i + 1)%3)] == i_1) {
            i_previous[0] = triangle_vertices[3*triangle_id + ((i + 2)%3)];
            i_next[1] = i_previous[0];
          }
          else if (   vertex_id == i_1
                      && triangle_vertices[3*triangle_id + ((i + 1)%3)] == i_0) {
            i_next[0] = triangle_vertices[3*triangle_id + ((i + 2)%3)];
            i_previous[1] = i_next[0];
          }

        } /* End of loop on triangle's vertices */

      } /* End of loop on triangles on each side of edge */

      /* Test quadrilateral's convexity */

      if (   _polygon_vertex_is_convex(i_previous[0],
                                       i_0,
                                       i_next[0],
                                       coords) == PDM_TRUE
             && _polygon_vertex_is_convex(i_previous[1],
                                          i_1,
                                          i_next[1],
                                          coords) == PDM_TRUE)
        convex_quad = PDM_TRUE;

      else
        convex_quad = PDM_FALSE;

      /* Check if edge is locally Delaunay */

      edge_locally_delaunay =
        _edge_is_locally_delaunay(i_0,
                                  i_1,
                                  vertex_flip[0],
                                  vertex_flip[1],
                                  coords);

      /* If edge is not locally Delaunay */
      /*---------------------------------*/

      if (   edge_locally_delaunay == PDM_FALSE
             && convex_quad == PDM_TRUE) {

        i_min = PDM_MIN(vertex_flip[0], vertex_flip[1]);
        i_max = PDM_MAX(vertex_flip[0], vertex_flip[1]);

        flip_edge_id = _EDGE_INDEX(i_min, i_max);

        for (j = 0; j < 2; j++) { /* Loop on triangles on each side of edge */

          triangle_id_0 = edge_neighbors[2*edge_id + j];
          triangle_id_1 = edge_neighbors[2*edge_id + ((j + 1)%2)];

          /* Redefine adjacent triangles */

          for (i = 0; i < 3; i++) {

            vertex_id = triangle_vertices[3*triangle_id_0 + i];

            if (vertex_id == edge_vertices[2*edge_id + ((j + 1)%2)])
              triangle_vertices[3*triangle_id_0 + i] = vertex_flip[(j + 1)%2];

          }

          /* Redefine triangle so that vertices appear in increasing order */

          _triangle_by_sorted_vertices(triangle_vertices + 3*triangle_id_0);

          /* Update neighborhood and non locally Delaunay edges */

          for (i = 0; i < 3; i++) { /* Loop on triangle's vertices */

            i_0 = triangle_vertices[3*triangle_id_0 + i];
            i_1 = triangle_vertices[3*triangle_id_0 + ((i + 1)%3)];

            i_min = PDM_MIN(i_0, i_1);
            i_max = PDM_MAX(i_0, i_1);

            current_edge_id = _EDGE_INDEX(i_min, i_max);

            if (current_edge_id < edge_id)
              restart = PDM_TRUE;

            /* If current edge is not the new (flip) egde ... */

            if (current_edge_id != flip_edge_id) {

              if (edge_neighbors[2*current_edge_id] == triangle_id_1)
                edge_neighbors[2*current_edge_id] = triangle_id_0;
              else if (edge_neighbors[2*current_edge_id + 1] == triangle_id_1)
                edge_neighbors[2*current_edge_id + 1] = triangle_id_0;

              /* ... and is not either a boundary edge */

              if (edge_neighbors[2*current_edge_id + 1] != -1)
                edge_is_delaunay[current_edge_id] = PDM_FALSE;

            }

          } /* End of loop on triangle's vertices */

        } /* End of loop on triangles on each side of edge */

        triangle_id_0 = edge_neighbors[2*edge_id];
        triangle_id_1 = edge_neighbors[2*edge_id + 1];

        edge_neighbors[2*flip_edge_id] = triangle_id_0;
        edge_neighbors[2*flip_edge_id + 1] = triangle_id_1;

        edge_neighbors[2*edge_id] = -1;
        edge_neighbors[2*edge_id + 1] = -1;

      } /* End if edge was not locally Delaunay */

    } /* End if edge was initially supposed non locally Delaunay */

    if (edge_id == n_edges - 1 && restart == PDM_TRUE) {
      restart = PDM_FALSE;
      edge_id = 0;
      n_restart += 1;

      //FIXME: Interruption de la boucle infinie du flip delaunay en ajoutant un nb de restart max (Reecrire l'algo de triangulation)

      if (n_restart > max_n_restart) {
        is_delaunay = 0;
        break;
      }
    }
    else if (edge_id == n_edges - 1 && restart == PDM_FALSE)
      face_is_delaunay = PDM_TRUE;
    else
      edge_id++;


  } /* End of flip algorithm */

  return is_delaunay;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a structure necessary to the polygon triangulation algorithm
 *
 * \param [in]   n_vertices_max  Maximum expected number of vertices per polygon
 *
 * \return  Pointer to polygon triangulation state structure
 *
 */

PDM_triangulate_state_t *
PDM_triangulate_state_create(const int  n_vertices_max)
{
  PDM_triangulate_state_t  *this_state = NULL;

  int n_edges_max = (2*n_vertices_max) - 3;
  int n_edges_tot_max = n_edges_max * (n_edges_max - 1) / 2;

  this_state = malloc (sizeof(PDM_triangulate_state_t));

  if (n_vertices_max > 3) {
    this_state->triangle_vertices = malloc (sizeof(int) * (n_vertices_max - 2) * 3);
    this_state->coords            = malloc (sizeof(double) * n_vertices_max*3);
    this_state->list_previous     = malloc (sizeof(int) * n_vertices_max);
    this_state->list_next         = malloc (sizeof(int) * n_vertices_max);
    this_state->edge_vertices     = malloc (sizeof(int) * n_edges_tot_max*2);
    this_state->edge_neighbors    = malloc (sizeof(int) * n_edges_tot_max*2);
    this_state->edge_is_delaunay  = malloc (sizeof(PDM_bool_t) * n_edges_tot_max);
    this_state->concave           = malloc (sizeof(PDM_bool_t) * n_vertices_max);
  }
  else {
    this_state->triangle_vertices = NULL;
    this_state->coords            = NULL;
    this_state->list_previous     = NULL;
    this_state->list_next         = NULL;
    this_state->edge_vertices     = NULL;
    this_state->edge_neighbors    = NULL;
    this_state->edge_is_delaunay  = NULL;
    this_state->concave           = NULL;
  }

  this_state->n_vertices_max = n_vertices_max;

  return this_state;
}



/**
 * \brief Destroy a structure necessary to the polygon triangulation algorithm
 *
 * \param [in]   this_state  Pointer to structure that should be destroyed
 *
 * \return  NULL pointer
 *
 */

PDM_triangulate_state_t *
PDM_triangulate_state_destroy(PDM_triangulate_state_t  *this_state)
{
  if (this_state != NULL) {
    if (this_state->triangle_vertices != NULL) {
      free (this_state->triangle_vertices);
      free (this_state->coords);
      free (this_state->list_previous);
      free (this_state->list_next);
      free (this_state->edge_vertices);
      free (this_state->edge_neighbors);
      free (this_state->edge_is_delaunay);
      free (this_state->concave);
    }
    free (this_state);
  }

  return NULL;
}



/**
 * \brief Triangulate a polygonal face
 *
 * For a polygon with n vertices, we should obtain a triangluation with
 * (n-2) triangles and (2n-3) edges. If the polygon_vertices argument
 * is NULL, 1, 2, ...,n local numbering is implied.
 *
 * \param [in]   dim                Spatial dimension (2 or 3)
 * \param [in]   n_vertices         Number of vertices defining the polygon
 * \param [in]   coords             Coordinates of the triangulation's vertices
 * \param [in]   parent_vertex_num  Optional indirection to vertex coordinates (1 to n)
 * \param [in]   polygon_vertices   Polygon connectivity; size: n_vertices or empty
 * \param [in]   mode               Triangles connectivity by vertex number or
 *                                  polygon vertex index (1 to n)
 * \param [out]   triangle_vertices Triangles connectivity
 *                                  (size = (\ref n_vertices - 2) * 3)
 * \param [in]   state              Associated triangulation state structure
 *
 * \return  Number of resulting triangles
 *
 */

int
PDM_triangulate_polygon(int                             dim,
                        int                             n_vertices,
                        const double                    coords[],
                        const PDM_l_num_t               parent_vertex_num[],
                        const PDM_l_num_t               polygon_vertices[],
                        PDM_triangulate_def_t           mode,
                        PDM_l_num_t                     triangle_vertices[],
                        PDM_triangulate_state_t  *const state)
{
  int i, j;

  int n_triangles = 0;
  int n_tries = 0;
  double epsilon[] = {1.0e-1, 1.0e-2, 0.0, -1.0e-2, -1.0e-1};

  int  *const list_previous = state->list_previous;
  int  *const list_next = state->list_next;
  PDM_bool_t  *const concave = state->concave;

  /* Initialize state structure */

  if (n_vertices > state->n_vertices_max) {

    int n_vertices_max = n_vertices*2;
    int n_edges_max = (2*n_vertices_max) - 3;
    int n_edges_tot_max = n_edges_max * (n_edges_max - 1) / 2;

    state->n_vertices_max = n_vertices_max;
    state->triangle_vertices = realloc (state->triangle_vertices, sizeof(int) * (n_vertices_max - 2) * 3);
    state->coords            = realloc (state->coords,            sizeof(double) * n_vertices_max*3);
    state->list_previous     = realloc (state->list_previous,     sizeof(int) * n_vertices_max);
    state->list_next         = realloc (state->list_next,         sizeof(int) * n_vertices_max);
    state->edge_vertices     = realloc (state->edge_vertices,     sizeof(int) * n_edges_tot_max*2);
    state->edge_neighbors    = realloc (state->edge_neighbors,    sizeof(int) * n_edges_tot_max*2);
    state->edge_is_delaunay  = realloc (state->edge_is_delaunay,  sizeof(PDM_bool_t) * n_edges_tot_max);
    state->concave           = realloc (state->concave,           sizeof(PDM_bool_t) * n_vertices_max);
  }

  if (parent_vertex_num != NULL) {
    if (polygon_vertices != NULL) {
      for (i = 0; i < n_vertices; i++) {
        int vertex_id = parent_vertex_num[polygon_vertices[i]-1] - 1;
        for (j = 0; j < dim; j++)
          state->coords[i*dim + j] = coords[vertex_id*dim + j];
      }
    }
    else {
      for (i = 0; i < (n_vertices * dim); i++) {
        int vertex_id = parent_vertex_num[i] - 1;
        state->coords[i] = coords[vertex_id];
      }
    }

  }
  else { /* (if parent_vertex_num == NULL) */

    if (polygon_vertices != NULL) {
      for (i = 0; i < n_vertices; i++) {
        for (j = 0; j < dim; j++)
          state->coords[i*dim + j] = coords[(polygon_vertices[i]-1)*dim + j];
      }
    }
    else {
      for (i = 0; i < (n_vertices * dim); i++)
        state->coords[i] = coords[i];
    }

  }

  /* Determine the work plane (3d coords are overwritten) */

  if (dim == 3)
    _polygon_plane_3d(n_vertices, state->coords);

  /* Initialization */

  while ((n_triangles != (n_vertices - 2)) &&
         n_tries < 5) {

    n_triangles = 0; /* Reset if previous try was a failure */

    for (i = 0; i < n_vertices; i++) {
      list_previous[i] = i - 1;
      list_next[i] = i + 1;
    }
    list_previous[0] = n_vertices - 1;
    list_next[n_vertices - 1] = 0;

    for (i = 0; i < n_vertices; i++) {
      if (_polygon_vertex_is_convex(list_previous[i],
                                    i,
                                    list_next[i],
                                    state->coords) == PDM_TRUE)
        concave[i] = PDM_FALSE;
      else
        concave[i] = PDM_TRUE;
    }

    i = 2;

    while (i != 0 && i != n_vertices) {

      if (_polygon_vertex_is_ear(n_vertices,
                                 list_previous[i],
                                 list_previous,
                                 list_next,
                                 state->concave,
                                 state->coords,
                                 epsilon[n_tries]) == PDM_TRUE) {

        /* Add a triangle with vertices list_previous[list_previous[i]],
           list_previous[i], i */

        state->triangle_vertices[n_triangles*3    ] = list_previous
          [list_previous[i]];
        state->triangle_vertices[n_triangles*3 + 1] = list_previous[i];
        state->triangle_vertices[n_triangles*3 + 2] = i;

        n_triangles += 1;

        /* Cut the ear corresponding to list_previous[i] */

        list_previous[i] = list_previous[list_previous[i]];
        list_next[list_previous[i]] = i;

        if (   (concave[i] == PDM_TRUE)
               && (_polygon_vertex_is_convex(list_previous[i],
                                             i,
                                             list_next[i],
                                             state->coords) == PDM_TRUE))
          concave[i] = PDM_FALSE;

        if (   (concave[list_previous[i]] == PDM_TRUE)
               && (_polygon_vertex_is_convex(list_previous[list_previous[i]],
                                             list_previous[i],
                                             i,
                                             state->coords) == PDM_TRUE))
          concave[list_previous[i]] = PDM_FALSE;

        if (list_previous[i] == 0)
          i = list_next[i];

      }
      else /* ! _polygon_vertex_is_ear(...) */

        i = list_next[i];

    }

    n_tries++;

  }


  /* Now that we have an initial triangulation, apply flip algorithm
     to obtain a Delaunay triangulation */

  int *_triangle_vertices;
  _triangle_vertices = malloc (sizeof(int) * (state->n_vertices_max - 2) * 3);

  for (int iii = 0; iii <  (state->n_vertices_max - 2) * 3; iii++) {
    _triangle_vertices[iii] =  state->triangle_vertices[iii];
  }

  int *__triangle_vertices = state->triangle_vertices;
  if (n_triangles == n_vertices - 2) {
    int is_delaunay = _polygon_delaunay_flip(n_vertices,
                                             state->triangle_vertices,
                                             state->edge_vertices,
                                             state->edge_neighbors,
                                             state->edge_is_delaunay,
                                             state->coords);
    if (!is_delaunay) {

      printf("Warning PDM_triangulate_polygon : Impossible to obtain a delaunay triangulation, keep initial triangulation from the 'ear' algortihm\n");
      printf("Coord des sommets : ");
      if (parent_vertex_num != NULL) {
        if (polygon_vertices != NULL) {
          for (i = 0; i < n_vertices; i++) {
            int vertex_id = parent_vertex_num[polygon_vertices[i]-1] - 1;
            printf("%16.9e %16.9e %16.9e ;", coords[3* vertex_id],
                   coords[3* vertex_id + 1],
                   coords[3* vertex_id + 2]);
          }
        }
        else {
          for (i = 0; i < n_vertices; i++) {
            int vertex_id = parent_vertex_num[i] - 1;
            printf("%16.9e %16.9e %16.9e ;", coords[3* vertex_id],
                   coords[3* vertex_id + 1],
                   coords[3* vertex_id + 2]);
          }
        }
      }
      else { /* (if parent_vertex_num == NULL) */

        if (polygon_vertices != NULL) {
          for (i = 0; i < n_vertices; i++) {
            int vertex_id = polygon_vertices[i]-1;
            printf("%16.9e %16.9e %16.9e ;", coords[3* vertex_id],
                   coords[3* vertex_id + 1],
                   coords[3* vertex_id + 2]);
          }
        }
        else {
          for (i = 0; i < n_vertices; i++) {
            printf("%16.9e %16.9e %16.9e ;", coords[3* i],
                   coords[3*i + 1],
                   coords[3*i + 2]);
          }
        }

      }

      printf("\n");
      __triangle_vertices = _triangle_vertices;
    }

  }

  /* Update triangle_vertices argument */

  if (polygon_vertices != NULL && mode == PDM_TRIANGULATE_MESH_DEF) {
    for (i = 0; i < n_triangles * 3; i++)
      triangle_vertices[i] = polygon_vertices[__triangle_vertices[i]];
  }
  else {
    for (i = 0; i < n_triangles * 3; i++)
      triangle_vertices[i] = __triangle_vertices[i] + 1;
  }

  free (_triangle_vertices);
  return n_triangles;
}



/**
 * \brief Triangulate a quadrangle
 *
 * A convex quadrangle is divided into two triangles along its shortest
 * diagonal. A non-convex quadrangle may only be divided along the diagonal
 * which lies inside the quadrangle.
 *
 * \param [in]   dim                  Spatial dimension (2 or 3)
 * \param [in]   coords               Coordinates of the triangulation's vertices
 * \param [in]   parent_vertex_num    Optional indirection to vertex coordinates (1 to 4)
 * \param [in]   quadrangle_vertices  Quadrangle connectivity (size = 4 or empty)
 * \param [out]  triangle_vertices    Triangles connectivity (size = 2 * 3)
 * \param [in]   state                Associated triangulation state structure
 *
 * \return  Number of resulting triangles
 *
 */

int
PDM_triangulate_quadrangle(int                 dim,
                           const double        coords[],
                           const PDM_l_num_t   parent_vertex_num[],
                           const PDM_l_num_t   quadrangle_vertices[],
                           PDM_l_num_t         triangle_vertices[])
{
  int i, j;
  double d2_02, d2_13;
  int o_count = 0, o_id = 0;
  PDM_l_num_t  vertex_id[4] = {0, 1, 2, 3};
  double v1[3] = {0.0, 0.0, 0.0}, v2[3] = {0.0, 0.0, 0.0};
  double n0[3] = {0.0, 0.0, 0.0}, ni[3] = {0.0, 0.0, 0.0};
  int    min_angle_idx;
  double min_cos;

  if (quadrangle_vertices != NULL) {
    for (i = 0; i < 4 ; i++)
      vertex_id[i] = quadrangle_vertices[i] - 1;
  }

  if (parent_vertex_num != NULL) {
    for (i = 0; i < 4 ; i++)
      vertex_id[i] = parent_vertex_num[i] - 1;
  }

  /* Check for an obtuse angle */

  for (i = 0; i < dim; i++) {
    v1[i] = coords[vertex_id[1]*dim + i] - coords[vertex_id[0]*dim + i];
    v2[i] = coords[vertex_id[3]*dim + i] - coords[vertex_id[0]*dim + i];
  }

  PDM_CROSS_PRODUCT(n0, v1, v2);
  min_cos = PDM_DOT_PRODUCT(v1, v2);
  min_angle_idx = 0;

  for (j = 1; j < 4; j++) {
    for (i = 0; i < dim; i++) {
      v1[i] = coords[vertex_id[(j+1)%4]*dim + i] - coords[vertex_id[j]*dim + i];
      v2[i] = coords[vertex_id[ j-1   ]*dim + i] - coords[vertex_id[j]*dim + i];
    }

    PDM_CROSS_PRODUCT(ni, v1, v2);
    double curr_cos = PDM_DOT_PRODUCT(v1, v2);
    if (curr_cos < min_cos) {
      min_cos = curr_cos;
      min_angle_idx = j;
    }

    if (PDM_DOT_PRODUCT(n0, ni) < 0) {
      o_count++;
      o_id = j;
    }
  }

  /* With an obtuse angle, only one diagonal lies inside the quadrangle;
     we define it as "shorter" */

  if (o_count > 0) {

    if (o_count > 1) {
      o_id = 0;
    }

    if (o_id%2 == 0) {
      d2_02 = 0.;
      d2_13 = 1.;
    }
    else {
      d2_02 = 1.;
      d2_13 = 0.;
    }

  }

  /* With no obtuse angle, we choose the largest angle */

  else {

    if ((min_angle_idx == 0) || (min_angle_idx == 2)) {
      d2_02 = 0.;
      d2_13 = 1.;
    }
    else {
      d2_02 = 1.;
      d2_13 = 0.;
    }

  }

  /* Now define first triangulation */

  if (quadrangle_vertices != NULL) {
    if (d2_02 < d2_13) {
      triangle_vertices[0] = quadrangle_vertices[0]; /* 1st triangle */
      triangle_vertices[1] = quadrangle_vertices[1];
      triangle_vertices[2] = quadrangle_vertices[2];
      triangle_vertices[3] = quadrangle_vertices[2]; /* 2nd triangle */
      triangle_vertices[4] = quadrangle_vertices[3];
      triangle_vertices[5] = quadrangle_vertices[0];
    }
    else {
      triangle_vertices[0] = quadrangle_vertices[0]; /* 1st triangle */
      triangle_vertices[1] = quadrangle_vertices[1];
      triangle_vertices[2] = quadrangle_vertices[3];
      triangle_vertices[3] = quadrangle_vertices[2]; /* 2nd triangle */
      triangle_vertices[4] = quadrangle_vertices[3];
      triangle_vertices[5] = quadrangle_vertices[1];
    }
  }
  else { /* if (quadrangle_vertices == NULL) */
    if (d2_02 < d2_13) {
      triangle_vertices[0] = 1; // 1st triangle
      triangle_vertices[1] = 2;
      triangle_vertices[2] = 3;
      triangle_vertices[3] = 3; // 2nd triangle
      triangle_vertices[4] = 4;
      triangle_vertices[5] = 1;
    }
    else {
      triangle_vertices[0] = 1; // 1st triangle
      triangle_vertices[1] = 2;
      triangle_vertices[2] = 4;
      triangle_vertices[3] = 3; // 2nd triangle
      triangle_vertices[4] = 4;
      triangle_vertices[5] = 2;
    }
  }

  /* Return number of triangles (for consistency with polygon triangulation) */

  return 2;
}



/**
 * \brief Tetrahedralize a pyramid
 *
 * \param [in]   dim               Spatial dimension (3)
 * \param [in]   coords            Coordinates of the vertices
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 5)
 * \param [in]   pyramid_vertices  Pyramid connectivity (size = 5)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 2 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_pyramid (int               dim,
                         const double      coords[],
                         const PDM_l_num_t parent_vertex_num[],
                         const PDM_l_num_t pyramid_vertices[],
                         PDM_l_num_t       tetra_vertices[])
{
  /* Triangulate base (quad) */
  PDM_l_num_t base_tri_vtx[6];
  PDM_triangulate_quadrangle(dim,
                             coords,
                             parent_vertex_num,
                             pyramid_vertices,
                             base_tri_vtx);

  /* 1st tetrahedron */
  tetra_vertices[0] = base_tri_vtx[0];
  tetra_vertices[1] = base_tri_vtx[1];
  tetra_vertices[2] = base_tri_vtx[2];
  tetra_vertices[3] = pyramid_vertices[4];

  /* 2nd tetrahedron */
  tetra_vertices[4] = base_tri_vtx[3];
  tetra_vertices[5] = base_tri_vtx[4];
  tetra_vertices[6] = base_tri_vtx[5];
  tetra_vertices[7] = pyramid_vertices[4];

  return 2;
}



/**
 * \brief Tetrahedralize a prism
 *
 * A simple look-up table is currently used,
 * thus no validity check is performed on the tetrahedra.
 *
 * \param [in]   dim               Spatial dimension (3) (unused)
 * \param [in]   coords            Coordinates of the vertices (unused)
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 6) (unused)
 * \param [in]   prism_vertices    Prism connectivity (size = 6)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 3 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_prism (int               dim,
                       const double      coords[],
                       const PDM_l_num_t parent_vertex_num[],
                       const PDM_l_num_t prism_vertices[],
                       PDM_l_num_t       tetrahedron_vertices[])
{

  PDM_UNUSED (dim);
  PDM_UNUSED (coords);
  PDM_UNUSED (parent_vertex_num);


  tetrahedron_vertices[ 0] = prism_vertices[0]; // 1st tetrahedron
  tetrahedron_vertices[ 1] = prism_vertices[1];
  tetrahedron_vertices[ 2] = prism_vertices[2];
  tetrahedron_vertices[ 3] = prism_vertices[3];
  tetrahedron_vertices[ 4] = prism_vertices[1]; // 2nd tetrahedron
  tetrahedron_vertices[ 5] = prism_vertices[2];
  tetrahedron_vertices[ 6] = prism_vertices[3];
  tetrahedron_vertices[ 7] = prism_vertices[4];
  tetrahedron_vertices[ 8] = prism_vertices[2]; // 3rd tetrahedron
  tetrahedron_vertices[ 9] = prism_vertices[3];
  tetrahedron_vertices[10] = prism_vertices[4];
  tetrahedron_vertices[11] = prism_vertices[5];

  return 3;


}



/**
 * \brief Tetrahedralize a hexahedron
 *
 * A simple look-up table is currently used,
 * thus no validity check is performed on the tetrahedra.
 *
 * \param [in]   dim               Spatial dimension (3) (unused)
 * \param [in]   coords            Coordinates of the vertices (unused)
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 8) (unused)
 * \param [in]   hexa_vertices     Hexahedron connectivity (size = 8)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 5 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_hexahedron (int               dim,
                            const double      coords[],
                            const PDM_l_num_t parent_vertex_num[],
                            const PDM_l_num_t hexa_vertices[],
                            PDM_l_num_t       tetrahedron_vertices[])
{
  PDM_UNUSED (dim);
  PDM_UNUSED (coords);
  PDM_UNUSED (parent_vertex_num);

  tetrahedron_vertices[ 0] = hexa_vertices[0]; // 1st tetrahedron
  tetrahedron_vertices[ 1] = hexa_vertices[1];
  tetrahedron_vertices[ 2] = hexa_vertices[2];
  tetrahedron_vertices[ 3] = hexa_vertices[4];
  tetrahedron_vertices[ 4] = hexa_vertices[1]; // 2nd tetrahedron
  tetrahedron_vertices[ 5] = hexa_vertices[4];
  tetrahedron_vertices[ 6] = hexa_vertices[5];
  tetrahedron_vertices[ 7] = hexa_vertices[7];
  tetrahedron_vertices[ 8] = hexa_vertices[1]; // 3rd tetrahedron
  tetrahedron_vertices[ 9] = hexa_vertices[2];
  tetrahedron_vertices[10] = hexa_vertices[3];
  tetrahedron_vertices[11] = hexa_vertices[7];
  tetrahedron_vertices[12] = hexa_vertices[2]; // 4th tetrahedron
  tetrahedron_vertices[13] = hexa_vertices[4];
  tetrahedron_vertices[14] = hexa_vertices[6];
  tetrahedron_vertices[15] = hexa_vertices[7];
  tetrahedron_vertices[16] = hexa_vertices[1]; // 5th tetrahedron
  tetrahedron_vertices[17] = hexa_vertices[2];
  tetrahedron_vertices[18] = hexa_vertices[4];
  tetrahedron_vertices[19] = hexa_vertices[7];

  return 5;


}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
