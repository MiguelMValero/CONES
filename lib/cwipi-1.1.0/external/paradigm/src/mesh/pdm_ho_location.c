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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * PDM library headers
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_nodal.h"
#include "pdm_ho_basis.h"
#include "pdm_geom_elem.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_triangulate.h"
#include "pdm_tetrahedron.h"

#include "pdm_vtk.h"
#include "pdm_ho_ordering.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_ho_location.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define S_HEAP 50

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _PDM_ho_location_user_elt_t {
  PDM_ho_location_fct_t location_in_elt;
} PDM_ho_location_user_elt_t;

/*============================================================================
 * LINEAR
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Sorted heap for sub-edge storage
 *----------------------------------------------------------------------------*/

typedef struct {

  int idx;

  double vtx_edge[3*2*S_HEAP];
  double uInPn_edge[1*2*S_HEAP];

  double closest_pt[3*S_HEAP];
  double closest_pt_uP1[1*S_HEAP];
  double closest_pt_uInPn[1*S_HEAP];
  double dist2[S_HEAP];

  int child[S_HEAP];

  int free_idx[S_HEAP];
  int sorted_idx[S_HEAP];

} _heap_l_t;


/*----------------------------------------------------------------------------
 * Function pointer to define a initial tesselation and push it in the heap
 *
 * parameters:
 *   heap              <-> Heap
 *   order             <-- element order
 *   n_nodes           <-- number of nodes
 *   local_to_user     <-- local to user ordering (for type)
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point coordinates
 *
 *----------------------------------------------------------------------------*/

typedef void (*_heap_fill_init_sub_edge_e)
(
 _heap_l_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 );

/*============================================================================
 * SURFACE
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Sorted heap for sub-triangle storage
 *----------------------------------------------------------------------------*/

typedef struct {

  int idx;

  double vtx_tria[3*3*S_HEAP];
  double uvInPn_tria[2*3*S_HEAP];

  double closest_pt[3*S_HEAP];
  double closest_pt_uvP1[2*S_HEAP];
  double closest_pt_uvInPn[2*S_HEAP];
  double dist2[S_HEAP];

  int child[S_HEAP];

  int free_idx[S_HEAP];
  int sorted_idx[S_HEAP];

} _heap_s_t;


/*----------------------------------------------------------------------------
 * Function pointer to define a initial tesselation and push it in the heap
 *
 * parameters:
 *   heap              <-> Heap
 *   order             <-- element order
 *   n_nodes           <-- number of nodes
 *   local_to_user     <-- local to user ordering (for type)
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point coordinates
 *
 *----------------------------------------------------------------------------*/

typedef void (*_heap_fill_init_sub_tria_t)
(
 _heap_s_t     *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 );

/*============================================================================
 * VOLUME
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Sorted heap for sub-tetrahedron storage
 *----------------------------------------------------------------------------*/

typedef struct {

  int idx;

  double vtx_tetra[3*4*S_HEAP];
  double uvwInPn_tetra[3*4*S_HEAP];

  double closest_pt[3*S_HEAP];
  double closest_pt_uvwP1[3*S_HEAP];
  double closest_pt_uvwInPn[3*S_HEAP];
  double dist2[S_HEAP];

  int child[S_HEAP];

  int free_idx[S_HEAP];
  int sorted_idx[S_HEAP];

} _heap_v_t;

/*----------------------------------------------------------------------------
 * Function pointer to define a initial tesselation and push it in the heap
 *
 * parameters:
 *   heap              <-> Heap
 *   order             <-- element order
 *   n_nodes           <-- number of nodes
 *   local_to_user     <-- local to user ordering (for type)
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point coordinates
 *
 *----------------------------------------------------------------------------*/

typedef void (*_heap_fill_init_sub_tetra_te)
(
 _heap_v_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 );


/*============================================================================
 * Static global variables
 *============================================================================*/
static PDM_ho_location_user_elt_t *_user_edge = NULL;

static PDM_ho_location_user_elt_t *_user_tria = NULL;
static PDM_ho_location_user_elt_t *_user_quad = NULL;

static PDM_ho_location_user_elt_t *_user_tetra = NULL;
static PDM_ho_location_user_elt_t *_user_pyra  = NULL;
static PDM_ho_location_user_elt_t *_user_prism = NULL;
static PDM_ho_location_user_elt_t *_user_hexa  = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Get u coordinate of ho edge nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   uv               --> u (size = n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_u_ho_edge_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 double       *u
 )
{
  int k = 0;
  int _order = order+1;


  double ustep = (umax - umin) / order;

  for (int i = 0; i < _order; i++) {
    double _u = umin + i * ustep;
    u[k++] = _u;
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uv coordinates of ho triangle nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   uv               --> uv (size = 2*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uv_ho_tria_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 double       *uv
 )
{
  int k = 0;
  int _order = order+1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;

  for (int j = 0; j < _order; j++) {
    double v = vmin + j * vstep;
    for (int i = 0; i < _order - j; i++) {
      double u = umin + i * ustep;

      uv[k++] = u;
      uv[k++] = v;

    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uv coordinates of ho quadrangle nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   uv               --> uv (size = 2*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uv_ho_quad_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 double       *uv
 )
{
  int k = 0;
  int _order = order+1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;

  for (int j = 0; j < _order; j++) {
    double v = vmin + j * vstep;
    for (int i = 0; i < _order; i++) {
      double u = umin + i * ustep;

      uv[k++] = u;
      uv[k++] = v;

    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uvw coordinates of ho tetrahedron nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   wmin             <-- w min
 *   wmax             <-- w max
 *   uvw              --> uvw (size = 3*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uvw_ho_tetra_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 const double  wmin,
 const double  wmax,
 double       *uvw
 )
{
  int p = 0;
  int _order = order + 1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  double wstep = (wmax - wmin) / order;


  for (int k = 0; k < _order; k++){
    double w = wmin + k * wstep;
    for (int j = 0; j < _order - k; j++) {
      double v = vmin + j * vstep;
      for (int i = 0; i < _order - j - k; i++) {
        double u = umin + i * ustep;

        uvw[p++] = u;
        uvw[p++] = v;
        uvw[p++] = w;

      }
    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uvw coordinates of ho pyramid nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   wmin             <-- w min
 *   wmax             <-- w max
 *   uvw              --> uvw (size = 3*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uvw_ho_pyra_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 const double  wmin,
 const double  wmax,
 double       *uvw
 )
{
  int p = 0;
  int _order = order + 1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  double wstep = (wmax - wmin) / order;


  for (int k = 0; k < _order; k++){
    double w = wmin + k * wstep;
    for (int j = 0; j < _order - k; j++) {
      double v = vmin + j * vstep;
      for (int i = 0; i < _order - k; i++) {
        double u = umin + i * ustep;

        uvw[p++] = u;
        uvw[p++] = v;
        uvw[p++] = w;

      }
    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uvw coordinates of ho prism nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   wmin             <-- w min
 *   wmax             <-- w max
 *   uvw              --> uvw (size = 3*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uvw_ho_prism_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 const double  wmin,
 const double  wmax,
 double       *uvw
 )
{
  int p = 0;
  int _order = order + 1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  double wstep = (wmax - wmin) / order;


  for (int k = 0; k < _order; k++){
    double w = wmin + k * wstep;
    for (int j = 0; j < _order; j++) {
      double v = vmin + j * vstep;
      for (int i = 0; i < _order - j; i++) {
        double u = umin + i * ustep;

        uvw[p++] = u;
        uvw[p++] = v;
        uvw[p++] = w;

      }
    }
  }
}


/*----------------------------------------------------------------------------
 *
 * Get uv coordinates of ho hexahedron nodes
 *
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   wmin             <-- w min
 *   wmax             <-- w max
 *   uvw              --> uvw (size = 3*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uvw_ho_hexa_nodes
(
 const int     order,
 const double  umin,
 const double  umax,
 const double  vmin,
 const double  vmax,
 const double  wmin,
 const double  wmax,
 double       *uvw
 )
{
  int p = 0;
  int _order = order+1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  double wstep = (wmax - wmin) / order;


  for (int k = 0; k < _order; k++){
    double w = wmin + k * wstep;
    for (int j = 0; j < _order; j++) {
      double v = vmin + j * vstep;
      for (int i = 0; i < _order; i++) {
        double u = umin + i * ustep;

        uvw[p++] = u;
        uvw[p++] = v;
        uvw[p++] = w;

      }
    }
  }
}


/*============================================================================
 * LINEAR
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Init heap 1D
 *
 * parameters:
 *   heap             <-- heap to initialize
 *
 *----------------------------------------------------------------------------*/

static void
_heap_l_init
(
 _heap_l_t *heap
 )
{
  heap->idx = -1;

  for (int i = 0; i < S_HEAP; i++) {
    heap->free_idx[i] = i;
  }
}


/*----------------------------------------------------------------------------
 *
 * Get top of the heap 1D
 *
 * parameters:
 *   heap             <-> heap to initialize
 *   order            <-> element order
 *   n_node           <-> number of nodes
 *   ho_vertex_num    <-> high order vertex num (internal ordering)
 *   vertex_coords    <-> vertex coordinates
 *   point_coords     <-> point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static int
_heap_l_top_get
(
 _heap_l_t *heap,
 double    *vtx_edge_current,
 double    *uInPn_edge_current,
 double    *closest_pt_current,
 double    *closest_pt_uP1_current,
 double    *closest_pt_uInPn_current,
 double    *dist2_current,
 int       *child
 )
{

  if (heap->idx < 0) {
    return 1;
  }

  int idx = heap->sorted_idx[heap->idx];
  heap->free_idx[heap->idx--]= idx;

  double *_vtx_edge_current = heap->vtx_edge + 6 * idx;
  for (int i = 0; i < 6; i++) {
    vtx_edge_current[i] = _vtx_edge_current[i];
  }

  double *_uInPn_edge_current = heap->uInPn_edge + 2 *idx;
  for (int i = 0; i < 2; i++) {
    uInPn_edge_current[i] = _uInPn_edge_current[i];
  }

  double *_closest_pt_current = heap->closest_pt + 3 *idx;
  for (int i = 0; i < 3; i++) {
    closest_pt_current[i] = _closest_pt_current[i];
  }

  double *_closest_pt_uP1_current = heap->closest_pt_uP1 + idx;
  double *_closest_pt_uInPn_current = heap->closest_pt_uInPn + idx;
  closest_pt_uP1_current[0] = _closest_pt_uP1_current[0];
  closest_pt_uInPn_current[0] = _closest_pt_uInPn_current[0];


  *child = heap->child[idx];
  *dist2_current = heap->dist2[idx];

  return 0;
}


/*----------------------------------------------------------------------------
 *
 * Add sub-edged of a pn-edge in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_l_insert
(
 _heap_l_t *heap,
 double    *vtx_edge,
 double    *uInPn_edge,
 double    *closest_pt,
 double    *closest_pt_uP1,
 double    *closest_pt_uInPn,
 double     dist2,
 int        child
 )
{
  /* Look for index (binary search) */

  int curr_idx = heap->idx;
  int *sorted_idx = heap->sorted_idx;
  double *sorted_dist2 = heap->dist2;

  int beg = 0;
  int end = curr_idx;

  while (beg <= end) {
    double dist2_beg = sorted_dist2[sorted_idx[beg]];
    double dist2_end = sorted_dist2[sorted_idx[end]];

    if (dist2 >= dist2_beg) {
      end = beg - 1;
    }

    else if (dist2 <= dist2_end) {
      beg = end + 1;
    }

    else {

      const int middle = (end + beg) / 2;
      if (beg == middle) {
        beg = beg + 1;
        end = beg - 1;
      }
      else {
        const double dist2_middle = sorted_dist2[sorted_idx[middle]];
        if (dist2 > dist2_middle) {
          end = middle;
        }
        else {
          beg = middle;
        }
      }
    }
  }


  /* If the heap is full remove the most distant */

  if (curr_idx >= (S_HEAP - 1)) {
    if (beg == 0) {
      return;
    }
    else {
      const int free_idx = sorted_idx[0];

      for (int i = 1; i < S_HEAP; i++) {
        sorted_idx[i-1] = sorted_idx[i];
      }

      heap->free_idx[heap->idx] = free_idx;
      heap->idx--;

      beg = beg - 1;
      end = end - 1;

    }
  }

  /* Add the element to the heap */

  heap->idx++;
  assert (heap->free_idx[heap->idx] != -1);

  for (int j = heap->idx; j > beg; j--) {
    sorted_idx[j] = sorted_idx[j-1];
  }

  sorted_idx[beg] = heap->free_idx[heap->idx];

  heap->free_idx[heap->idx] = -1;

  int _idx = sorted_idx[beg];

  for (int j = 0; j < 6; j++) {
    heap->vtx_edge[6*_idx+j] = vtx_edge[j];
  }

  for (int j = 0; j < 2; j++) {
    heap->uInPn_edge[2*_idx+j] = uInPn_edge[j];
  }

  for (int j = 0; j < 3; j++) {
    heap->closest_pt[3*_idx+j]  = closest_pt[j];
  }

  for (int j = 0; j < 1; j++) {
    heap->closest_pt_uP1[_idx+j] = closest_pt_uP1[j];
    heap->closest_pt_uInPn[_idx+j] = closest_pt_uInPn[j];
  }

  heap->dist2[_idx] = dist2;
  heap->child[_idx] = child;
}


/*----------------------------------------------------------------------------
 *
 * Add sub-edges of a pn-edge in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_sub_edge
(
 _heap_l_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{
  double *uNodes = malloc (sizeof(double) * n_nodes);

  _u_ho_edge_nodes (order,
                    0., 1.,
                    uNodes);

  int child = 0;
  for (int i = 0; i < n_nodes-1; i++) {
    int idx1 = i;
    int idx2 = i+1;

    double x1 = nodes_coords[3*idx1    ];
    double y1 = nodes_coords[3*idx1 + 1];
    double z1 = nodes_coords[3*idx1 + 2];


    double x2 = nodes_coords[3*idx2    ];
    double y2 = nodes_coords[3*idx2 + 1];
    double z2 = nodes_coords[3*idx2 + 2];

    double __vertex_coords[6] = {x1, y1, z1,
                                 x2, y2, z2};
    double _closest_pointP1[3];
    double _uClosestPointPn[1];
    double _weightsClosestPointP1[2];
    double _dist2;


    int is_degenerate = PDM_line_evaluate_position (point_coords,
                                                    __vertex_coords,
                                                    _closest_pointP1,
                                                    &_dist2,
                                                    _weightsClosestPointP1);

    double _uPn_sub_edge[2];
    double aux1 = uNodes[idx1];
    double aux2 = uNodes[idx2];
    _uPn_sub_edge[0] = aux1;
    _uPn_sub_edge[1] = aux2;


    if (is_degenerate != -1) {

      double *_uClosestPointP1 = _weightsClosestPointP1 + 1;

      _uClosestPointPn[0] = 0;
      for (int j = 0; j < 2; j++) {
        _uClosestPointPn[0] += _weightsClosestPointP1[j] * _uPn_sub_edge[j];
      }

      _heap_l_insert (heap,
                      __vertex_coords,
                      _uPn_sub_edge,
                      _closest_pointP1,
                      _uClosestPointP1,
                      _uClosestPointPn,
                      _dist2,
                      child++);
    }



  }
  free (uNodes);
}


/*-----------------------------------------------------------------------------
 *
 * Add children to a heap
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   vtx_tria_current <-- current triangle
 *   uvPn_tria_current<-- uv of current triangle vertices in the on element
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static void
_insert_subedge
(
 _heap_l_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      vtx_edge_current[],
 double                      uPn_edge_current[]
 )
{
  double _vtx_edge_children[9];
  double _uPn_edge_children[3];


  const int idx_sub_edge[4] = {0, 2,
                               2, 1};

  /* Compute middle vertices */

  for (int i = 0; i < 6; i++) {
    _vtx_edge_children[i] = vtx_edge_current[i];
  }

  for (int i = 0; i < 2; i++) {
    _uPn_edge_children[i] = uPn_edge_current[i];
  }

  _uPn_edge_children[2] = 0.5 * (uPn_edge_current[0] + uPn_edge_current[1]);


  PDM_ho_basis (type,
                order,
                n_nodes,
                1,
                _uPn_edge_children + 2,
                weightsPn);


  for (int j = 0; j < 3; j++) {
    _vtx_edge_children[6+j] = 0.;
  }
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < n_nodes; j++) {
      _vtx_edge_children[6 + k] += weightsPn[j] * nodes_coords[k + 3*j];
    }
  }

  int child = 0;
  for (int i = 0; i < 2; i++) {

    double _vtx_edge_child[6];
    double _uPn_edge_child[2];

    for (int j = 0; j < 2; j++) {
      int _j = idx_sub_edge[2 * i + j];
      for (int k = 0; k < 3; k++) {
        _vtx_edge_child[3*j + k] = _vtx_edge_children[3*_j + k];
      }
      _uPn_edge_child[j] = _uPn_edge_children[_j];

    }


    double _closest_pt_child[3];
    double _closest_pt_uPn_child[1];
    double _dist2_child = 0.;
    double _closest_pt_weights_child[2];

    int is_degenerate = PDM_line_evaluate_position (point_coords,
                                                    _vtx_edge_child,
                                                    _closest_pt_child,
                                                    &_dist2_child,
                                                    _closest_pt_weights_child);


    if (is_degenerate == -1) {
      continue;
    }

    double *_closest_pt_uP1_child = _closest_pt_weights_child + 1;

    _closest_pt_uPn_child[0] = 0.;
    for (int k = 0; k < 2; k++) {
      _closest_pt_uPn_child[0] += _closest_pt_weights_child[k] * _uPn_edge_child[k];
    }


    _heap_l_insert (heap,
                    _vtx_edge_child,
                    _uPn_edge_child,
                    _closest_pt_child,
                    _closest_pt_uP1_child,
                    _closest_pt_uPn_child,
                    _dist2_child,
                    child++);
  }
}


/*-----------------------------------------------------------------------------
 *
 * compute distance from closest edge subdivision
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current edge
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   uncertain_result --> 1 if the result is uncertain
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_closest_edge_subdivision
(
 _heap_l_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj,
 int                        *uncertain_result
 )
{
  *uncertain_result = 0;
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;
  double dist2_pre = HUGE_VAL;
  int distance_extension = 0;

  while (1) {

    double _vtx_edge_current[6];
    double _uPn_edge_current[2];

    double _closest_pt_current[3];
    double _closest_pt_uP1_current[1];
    double _closest_pt_uPn_current[1];
    double _dist2_current;

    /* Get closest edge stored in the heap */

    int _child;
    int is_empty = _heap_l_top_get (heap,
                                    _vtx_edge_current,
                                    _uPn_edge_current,
                                    _closest_pt_current,
                                    _closest_pt_uP1_current,
                                    _closest_pt_uPn_current,
                                    &_dist2_current,
                                    &_child);

    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty %s\n");
      abort();
    }

    if ((distance_extension == 0) && (_dist2_current > dist2_pre)) {
      distance_extension = 1;
    }

    else if (distance_extension == 1) {
      if (_dist2_current <= dist2_min_min) {
        distance_extension = 0;
      }
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);
    dist2_pre = _dist2_current;

    /* Compute projected from current P1 edge */
    double weightsP1[2];
    PDM_ho_basis (PDM_MESH_NODAL_BARHO,
                  1,
                  2,
                  1,
                  _closest_pt_uP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_edge_current[3*j+k];
      }
    }

    /* Compute projected from current Pn edge */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }



    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */
    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uPn_current[0];

      break;
    }

    /*
     * Insert sub-edge in the heap
     */

    _insert_subedge (heap,
                     order,
                     type,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_edge_current,
                     _uPn_edge_current);

  }

  if (distance_extension) {
    *uncertain_result = 1;
  }

  return dist2;
}


/*-----------------------------------------------------------------------------
 *
 * compute distance from uniform edge subdivision
 *
 * parameters:
 *   heap1            <-> heap
 *   heap2            <-> work heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_uniform_edge_subdivision
(
 _heap_l_t                  *heap1,
 _heap_l_t                  *heap2,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj
 )
{
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;

  _heap_l_t *heap      = heap1;
  _heap_l_t *next_heap = heap2;

  while (1) {

    double _vtx_edge_current[6];
    double _uPn_edge_current[2];

    double _closest_pt_current[3];
    double _closest_pt_uP1_current[1];
    double _closest_pt_uPn_current[1];
    double _dist2_current;

    /* Get closest edge stored in the heap */

    int _child;
    int is_empty = _heap_l_top_get (heap,
                                    _vtx_edge_current,
                                    _uPn_edge_current,
                                    _closest_pt_current,
                                    _closest_pt_uP1_current,
                                    _closest_pt_uPn_current,
                                    &_dist2_current,
                                    &_child);


    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty %s\n");
      abort();
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);

    /* Compute projected from current P1 edge */
    double weightsP1[2];
    PDM_ho_basis (PDM_MESH_NODAL_BARHO,
                  1,
                  2,
                  1,
                  _closest_pt_uP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_edge_current[3*j+k];
      }
    }

    /* Compute projected from current Pn edge */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }

    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */
    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uPn_current[0];

      break;
    }

    /*
     * Insert sub-edge in the next heap
     */

    _heap_l_init (next_heap);

    _insert_subedge (next_heap,
                     order,
                     type,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_edge_current,
                     _uPn_edge_current);

    double _vtx_edge_current2[6];
    double _uPn_edge_current2[2];

    double _closest_pt_current2[3];
    double _closest_pt_uP1_current2[1];
    double _dist2_current2;
    int _child_current2;

    while ( !_heap_l_top_get (heap,
                              _vtx_edge_current2,
                              _uPn_edge_current2,
                              _closest_pt_current2,
                              _closest_pt_uP1_current2,
                              _closest_pt_uPn_current,
                              &_dist2_current2,
                              &_child_current2)) {


      _insert_subedge (next_heap,
                       order,
                       type,
                       n_nodes,
                       nodes_coords,
                       point_coords,
                       weightsPn,
                       _vtx_edge_current2,
                       _uPn_edge_current2);

    }

    _heap_l_t *heap_tmp = heap;
    heap = next_heap;
    next_heap = heap_tmp;

  }


  return dist2;
}


/*-----------------------------------------------------------------------------
 *
 * Default point location on a high order segment
 *
 * parameters:
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   uvw              --> parametric coordinates in the element
 *
 * return:
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double
_default_location_generic_1d
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const double                char_size,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw,
 _heap_fill_init_sub_edge_e  fill_init_fct
 )
{
  const int n_it_max = 100;
  double err_max = PDM_MAX (char_size * 1e-6, 1e-15);

  double dist2 = HUGE_VAL;

  _heap_l_t heap;
  _heap_l_t heap2;

  double *weightsPn = malloc (sizeof(double) * n_nodes);

  /* Initialize heap */

  _heap_l_init (&heap);

  /* Build initial sub-edge and store them in the heap */

  (fill_init_fct) (&heap,
                   order,
                   n_nodes,
                   nodes_coords,
                   point_coords);

  /*
   *  While error > error_max
   *    - Get closest edge in the heap
   *    - Cut it in sub-edge
   *    - Store them in the heap
   */


  const int method = 0;

  int n_it;
  double err_proj = HUGE_VAL;
  int uncertain_result = 0;;

  if (method == 0) {
    dist2 = _compute_dist2_from_closest_edge_subdivision (&heap,
                                                          order,
                                                          type,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj,
                                                          &uncertain_result);

    if (uncertain_result) {

      /* Initialize heap */

      _heap_l_init (&heap);
      _heap_l_init (&heap2);

      /* Build initial sub-edge and store them in the heap */

      (fill_init_fct) (&heap,
                       order,
                       n_nodes,
                       nodes_coords,
                       point_coords);


      dist2 = _compute_dist2_from_uniform_edge_subdivision (&heap,
                                                            &heap2,
                                                            order,
                                                            type,
                                                            n_nodes,
                                                            n_it_max,
                                                            err_max,
                                                            nodes_coords,
                                                            point_coords,
                                                            weightsPn,
                                                            projected_coords,
                                                            uvw,
                                                            &n_it,
                                                            &err_proj);

    }
  }

  else {

    _heap_l_init (&heap2);
    dist2 = _compute_dist2_from_uniform_edge_subdivision (&heap,
                                                          &heap2,
                                                          order,
                                                          type,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj);

  }

  free (weightsPn);

  return dist2;

}


/*============================================================================
 * SURFACE
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Init heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *
 *----------------------------------------------------------------------------*/

static void
_heap_s_init
(
 _heap_s_t *heap
 )
{
  heap->idx = -1;

  for (int i = 0; i < S_HEAP; i++) {
    heap->free_idx[i] = i;
  }
}


/*----------------------------------------------------------------------------
 *
 * Get top of the heap
 *
 * parameters:
 *   heap             <-> heap to initialize
 *   order            <-> element order
 *   n_node           <-> number of nodes
 *   ho_vertex_num    <-> high order vertex num (internal ordering)
 *   vertex_coords    <-> vertex coordinates
 *   point_coords     <-> point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static int
_heap_s_top_get
(
 _heap_s_t *heap,
 double    *vtx_tria_current,
 double    *uvInPn_tria_current,
 double    *closest_pt_current,
 double    *closest_pt_uvP1_current,
 double    *closest_pt_uvInPn_current,
 double    *dist2_current,
 int       *child
 )
{
  if (heap->idx < 0) {
    return 1;
  }

  int idx = heap->sorted_idx[heap->idx];
  heap->free_idx[heap->idx--]= idx;

  double *_vtx_tria_current = heap->vtx_tria + 9 * idx;
  for (int i = 0; i < 9; i++) {
    vtx_tria_current[i] = _vtx_tria_current[i];
  }

  double *_uvInPn_tria_current = heap->uvInPn_tria + 6 *idx;
  for (int i = 0; i < 6; i++) {
    uvInPn_tria_current[i] = _uvInPn_tria_current[i];
  }

  double *_closest_pt_current = heap->closest_pt + 3 *idx;
  for (int i = 0; i < 3; i++) {
    closest_pt_current[i] = _closest_pt_current[i];
  }

  double *_closest_pt_uvP1_current = heap->closest_pt_uvP1 + 2 *idx;
  double *_closest_pt_uvInPn_current = heap->closest_pt_uvInPn + 2 *idx;
  for (int i = 0; i < 2; i++) {
    closest_pt_uvP1_current[i] = _closest_pt_uvP1_current[i];
    closest_pt_uvInPn_current[i] = _closest_pt_uvInPn_current[i];
  }

  *child = heap->child[idx];
  *dist2_current = heap->dist2[idx];

  return 0;
}


/*----------------------------------------------------------------------------
 *
 * Add sub-triangles of a pn-triangle in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_s_insert
(
 _heap_s_t *heap,
 double    *vtx_tria,
 double    *uvInPn_tria,
 double    *closest_pt,
 double    *closest_pt_uvP1,
 double    *closest_pt_uvInPn,
 double     dist2,
 int        child
 )
{
  /* Look for index (binary search) */

  int curr_idx = heap->idx;
  int *sorted_idx = heap->sorted_idx;
  double *sorted_dist2 = heap->dist2;

  int beg = 0;
  int end = curr_idx;

  while (beg <= end) {
    double dist2_beg = sorted_dist2[sorted_idx[beg]];
    double dist2_end = sorted_dist2[sorted_idx[end]];

    if (dist2 >= dist2_beg) {
      end = beg - 1;
    }

    else if (dist2 <= dist2_end) {
      beg = end + 1;
    }

    else {

      const int middle = (end + beg) / 2;
      if (beg == middle) {
        beg = beg + 1;
        end = beg - 1;
      }
      else {
        const double dist2_middle = sorted_dist2[sorted_idx[middle]];
        if (dist2 > dist2_middle) {
          end = middle;
        }
        else {
          beg = middle;
        }
      }
    }
  }


  /* If the heap is full remove the most distant */

  if (curr_idx >= (S_HEAP - 1)) {
    if (beg == 0) {
      return;
    }
    else {
      const int free_idx = sorted_idx[0];

      for (int i = 1; i < S_HEAP; i++) {
        sorted_idx[i-1] = sorted_idx[i];
      }

      heap->free_idx[heap->idx] = free_idx;
      heap->idx--;

      beg = beg - 1;
      end = end - 1;

    }
  }

  /* Add the element to the heap */

  heap->idx++;
  assert (heap->free_idx[heap->idx] != -1);

  for (int j = heap->idx; j > beg; j--) {
    sorted_idx[j] = sorted_idx[j-1];
  }

  sorted_idx[beg] = heap->free_idx[heap->idx];

  heap->free_idx[heap->idx] = -1;

  int _idx = sorted_idx[beg];

  for (int j = 0; j < 9; j++) {
    heap->vtx_tria[9*_idx+j] = vtx_tria[j];
  }

  for (int j = 0; j < 6; j++) {
    heap->uvInPn_tria[6*_idx+j] = uvInPn_tria[j];
  }

  for (int j = 0; j < 3; j++) {
    heap->closest_pt[3*_idx+j]  = closest_pt[j];
  }

  for (int j = 0; j < 2; j++) {
    heap->closest_pt_uvP1[2*_idx+j] = closest_pt_uvP1[j];
    heap->closest_pt_uvInPn[2*_idx+j] = closest_pt_uvInPn[j];
  }

  heap->dist2[_idx] = dist2;
  heap->child[_idx] = child;
}


/*----------------------------------------------------------------------------
 *
 * Add sub-triangles of a pn-triangle in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_sub_tria
(
 _heap_s_t     *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{
  int ibeg = 0;
  int iend = order;

  double *uv_nodes = malloc (sizeof(double) * 2 * n_nodes);

  _uv_ho_tria_nodes (order,
                     0., 1.,
                     0., 1.,
                     uv_nodes);

  int child = 0;
  for (int j = 0; j < order; j++) {
    int k1 = 0;
    for (int i = ibeg; i < iend - 1; i++) {

      int idx1 = i;
      int idx2 = i+1;
      int idx3 = iend + 1 + k1;
      int idx4 = iend + 2 + k1;

      double x1 = nodes_coords[3*idx1];
      double y1 = nodes_coords[3*idx1 + 1];
      double z1 = nodes_coords[3*idx1 + 2];

      double x2 = nodes_coords[3*idx2];
      double y2 = nodes_coords[3*idx2 + 1];
      double z2 = nodes_coords[3*idx2 + 2];

      double x3 = nodes_coords[3*idx3];
      double y3 = nodes_coords[3*idx3 + 1];
      double z3 = nodes_coords[3*idx3 + 2];

      double x4 = nodes_coords[3*idx4];
      double y4 = nodes_coords[3*idx4 + 1];
      double z4 = nodes_coords[3*idx4 + 2];

      double __vertex_coords[9] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3};
      double _closest_pointP1[3];
      // double *_uvClosestPointP1;
      double _uvClosestPointP1[2];
      double _uvClosestPointPn[2];
      double _weightsClosestPointP1[3];
      double _dist2;

      PDM_triangle_status_t stat = PDM_triangle_evaluate_position(point_coords,
                                                                  __vertex_coords,
                                                                  _closest_pointP1,
                                                                  &_dist2,
                                                                  _weightsClosestPointP1);

      double _uvPn_sub_tria[6];

      _uvPn_sub_tria[0] = uv_nodes[2*idx1];
      _uvPn_sub_tria[1] = uv_nodes[2*idx1+1];
      _uvPn_sub_tria[2] = uv_nodes[2*idx2];
      _uvPn_sub_tria[3] = uv_nodes[2*idx2+1];
      _uvPn_sub_tria[4] = uv_nodes[2*idx3];
      _uvPn_sub_tria[5] = uv_nodes[2*idx3+1];

      if (stat != PDM_TRIANGLE_DEGENERATED) {

        // _uvClosestPointP1 = _weightsClosestPointP1;
        for (int j1 = 0; j1 < 2; j1++) {
          _uvClosestPointP1[j1] = _weightsClosestPointP1[j1+1];
          _uvClosestPointPn[j1] = 0;
        }
        for (int j1 = 0; j1 < 2; j1++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
          }
        }


        _heap_s_insert (heap,
                        __vertex_coords,
                        _uvPn_sub_tria,
                        _closest_pointP1,
                        _uvClosestPointP1,
                        _uvClosestPointPn,
                        _dist2,
                        child++);
      }

      __vertex_coords[0] = x2;
      __vertex_coords[1] = y2;
      __vertex_coords[2] = z2;
      __vertex_coords[3] = x4;
      __vertex_coords[4] = y4;
      __vertex_coords[5] = z4;
      __vertex_coords[6] = x3;
      __vertex_coords[7] = y3;
      __vertex_coords[8] = z3;

      stat = PDM_triangle_evaluate_position(point_coords,
                                            __vertex_coords,
                                            _closest_pointP1,
                                            &_dist2,
                                            _weightsClosestPointP1);

      _uvPn_sub_tria[0] = uv_nodes[2*idx2];
      _uvPn_sub_tria[1] = uv_nodes[2*idx2+1];
      _uvPn_sub_tria[2] = uv_nodes[2*idx4];
      _uvPn_sub_tria[3] = uv_nodes[2*idx4+1];
      _uvPn_sub_tria[4] = uv_nodes[2*idx3];
      _uvPn_sub_tria[5] = uv_nodes[2*idx3+1];

      if (stat != PDM_TRIANGLE_DEGENERATED) {

        // _uvClosestPointP1 = _weightsClosestPointP1;
        for (int j1 = 0; j1 < 2; j1++) {
          _uvClosestPointP1[j1] = _weightsClosestPointP1[j1+1];
          _uvClosestPointPn[j1] = 0;
        }
        for (int j1 = 0; j1 < 2; j1++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
          }
        }

        _heap_s_insert (heap,
                        __vertex_coords,
                        _uvPn_sub_tria,
                        _closest_pointP1,
                        _uvClosestPointP1,
                        _uvClosestPointPn,
                        _dist2,
                        child++);
      }

      k1++;
    }

    int idx1 = iend - 1;
    int idx2 = iend - 1 + 1;
    int idx3 = iend + 1 + k1;

    double x1 = nodes_coords[3*idx1];
    double y1 = nodes_coords[3*idx1 + 1];
    double z1 = nodes_coords[3*idx1 + 2];

    double x2 = nodes_coords[3*idx2];
    double y2 = nodes_coords[3*idx2 + 1];
    double z2 = nodes_coords[3*idx2 + 2];

    double x3 = nodes_coords[3*idx3];
    double y3 = nodes_coords[3*idx3 + 1];
    double z3 = nodes_coords[3*idx3 + 2];

    double __vertex_coords[9] = {x1, y1, z1,
                                 x2, y2, z2,
                                 x3, y3, z3};

    double _closest_pointP1[3];
    // double *_uvClosestPointP1;
    double _uvClosestPointP1[2];
    double _uvClosestPointPn[2];

    double _weightsClosestPointP1[3];

    double _dist2;

    PDM_triangle_status_t stat = PDM_triangle_evaluate_position(point_coords,
                                                                __vertex_coords,
                                                                _closest_pointP1,
                                                                &_dist2,
                                                                _weightsClosestPointP1);

    double _uvPn_sub_tria[6];

    _uvPn_sub_tria[0] = uv_nodes[2*idx1];
    _uvPn_sub_tria[1] = uv_nodes[2*idx1+1];
    _uvPn_sub_tria[2] = uv_nodes[2*idx2];
    _uvPn_sub_tria[3] = uv_nodes[2*idx2+1];
    _uvPn_sub_tria[4] = uv_nodes[2*idx3];
    _uvPn_sub_tria[5] = uv_nodes[2*idx3+1];

    if (stat != PDM_TRIANGLE_DEGENERATED) {

      // _uvClosestPointP1 = _weightsClosestPointP1;
      for (int j1 = 0; j1 < 2; j1++) {
        _uvClosestPointP1[j1] = _weightsClosestPointP1[j1+1];
        _uvClosestPointPn[j1] = 0;
      }
      for (int j1 = 0; j1 < 2; j1++) {
        for (int k = 0; k < 3; k++) {
          _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
        }
      }

      _heap_s_insert (heap,
                      __vertex_coords,
                      _uvPn_sub_tria,
                      _closest_pointP1,
                      _uvClosestPointP1,
                      _uvClosestPointPn,
                      _dist2,
                      child++);
    }

    ibeg = iend + 1;
    iend += order - j;
  }

  free (uv_nodes);
}


/*----------------------------------------------------------------------------
 *
 * Add sub-triangles of a qn-quadrangle in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_qn_sub_tria
(
 _heap_s_t     *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{
  double *uvNodes = malloc (sizeof(double) * 2 * n_nodes);

  _uv_ho_quad_nodes (order,
                     0., 1.,
                     0., 1.,
                     uvNodes);

  int child = 0;
  int step = order + 1;

  for (int j = 0; j < order; j++) {
    int j1 = j+1;

    for (int i = 0; i < order; i++) {

      int i1 = i+1;

      int idx1 = j*step  + i;
      int idx2 = j*step  + i1;
      int idx3 = j1*step + i;
      int idx4 = j1*step + i1;

      double x1 = nodes_coords[3*idx1];
      double y1 = nodes_coords[3*idx1 + 1];
      double z1 = nodes_coords[3*idx1 + 2];

      double x2 = nodes_coords[3*idx2];
      double y2 = nodes_coords[3*idx2 + 1];
      double z2 = nodes_coords[3*idx2 + 2];

      double x3 = nodes_coords[3*idx3];
      double y3 = nodes_coords[3*idx3 + 1];
      double z3 = nodes_coords[3*idx3 + 2];

      double x4 = nodes_coords[3*idx4];
      double y4 = nodes_coords[3*idx4 + 1];
      double z4 = nodes_coords[3*idx4 + 2];

      double __vertex_coords[9] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3};
      double _closest_pointP1[3];
      // double *_uvClosestPointP1;
      double _uvClosestPointP1[2];
      double _uvClosestPointPn[2];
      double _weightsClosestPointP1[3];
      double _dist2;

      PDM_triangle_status_t stat = PDM_triangle_evaluate_position(point_coords,
                                                                  __vertex_coords,
                                                                  _closest_pointP1,
                                                                  &_dist2,
                                                                  _weightsClosestPointP1);

      double _uvPn_sub_tria[6];

      _uvPn_sub_tria[0] = uvNodes[2*idx1];
      _uvPn_sub_tria[1] = uvNodes[2*idx1+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx2];
      _uvPn_sub_tria[3] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];

      if (stat != PDM_TRIANGLE_DEGENERATED) {

        // _uvClosestPointP1 = _weightsClosestPointP1;
        for (int j2 = 0; j2 < 2; j2++) {
          _uvClosestPointP1[j2] = _weightsClosestPointP1[j2+1];
          _uvClosestPointPn[j2] = 0;
        }
        for (int j2 = 0; j2 < 2; j2++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j2] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j2];
          }
        }

        _heap_s_insert (heap,
                        __vertex_coords,
                        _uvPn_sub_tria,
                        _closest_pointP1,
                        _uvClosestPointP1,
                        _uvClosestPointPn,
                        _dist2,
                        child++);
      }

      __vertex_coords[0] = x2;
      __vertex_coords[1] = y2;
      __vertex_coords[2] = z2;
      __vertex_coords[3] = x4;
      __vertex_coords[4] = y4;
      __vertex_coords[5] = z4;
      __vertex_coords[6] = x3;
      __vertex_coords[7] = y3;
      __vertex_coords[8] = z3;

      stat = PDM_triangle_evaluate_position(point_coords,
                                            __vertex_coords,
                                            _closest_pointP1,
                                            &_dist2,
                                            _weightsClosestPointP1);

      _uvPn_sub_tria[0] = uvNodes[2*idx2];
      _uvPn_sub_tria[1] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx4];
      _uvPn_sub_tria[3] = uvNodes[2*idx4+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];

      if (stat != PDM_TRIANGLE_DEGENERATED) {

        // _uvClosestPointP1 = _weightsClosestPointP1;
        for (int j2 = 0; j2 < 2; j2++) {
          _uvClosestPointP1[j2] = _weightsClosestPointP1[j2+1];
          _uvClosestPointPn[j2] = 0;
        }
        for (int j2 = 0; j2 < 2; j2++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j2] +=
              _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j2];
          }
        }

        _heap_s_insert (heap,
                        __vertex_coords,
                        _uvPn_sub_tria,
                        _closest_pointP1,
                        _uvClosestPointP1,
                        _uvClosestPointPn,
                        _dist2,
                        child++);
      }
    }
  }

  free (uvNodes);
}


/*-----------------------------------------------------------------------------
 *
 * Add children to a heap
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   vtx_tria_current <-- current triangle
 *   uvPn_tria_current<-- uv of current triangle vertices in the on element
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static void
_insert_subtria
(
 _heap_s_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      vtx_tria_current[],
 double                      uvPn_tria_current[]
 )
{
  double _vtx_tria_children[18];
  double _uvPn_tria_children[12];

  const int idx_sub_tria[12] = {0, 3, 5,
                                3, 4, 5,
                                3, 1, 4,
                                5, 4, 2};

  /* Compute middle vertices */

  for (int i = 0; i < 9; i++) {
    _vtx_tria_children[i] = vtx_tria_current[i];
  }

  for (int i = 0; i < 6; i++) {
    _uvPn_tria_children[i] = uvPn_tria_current[i];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      _uvPn_tria_children[6+2*i+j] = 0.5 *
        (uvPn_tria_current[2*i+j] + uvPn_tria_current[2*((i+1)%3)+j]);
    }
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _uvPn_tria_children + 6 + 2*i,
                  weightsPn);

    for (int j = 0; j < 3; j++) {
      _vtx_tria_children[9+3*i+j] = 0;
    }
    for (int k = 0; k < n_nodes; k++) {
      const double *_node_coords = nodes_coords + 3 * k;
      for (int j = 0; j < 3; j++) {
        _vtx_tria_children[9+3*i+j] += weightsPn[k] * _node_coords[j];
      }
    }
  }

  int child = 0;
  for (int i = 0; i < 4; i++) {

    double _vtx_tria_child[9];
    double _uvPn_tria_child[6];

    for (int j = 0; j < 3; j++) {
      int _j = idx_sub_tria[3 * i + j];
      for (int k = 0; k < 3; k++) {
        _vtx_tria_child[3*j + k] = _vtx_tria_children[3*_j + k];
      }
      for (int k = 0; k < 2; k++) {
        _uvPn_tria_child[2*j + k] = _uvPn_tria_children[2*_j + k];
      }
    }

    double _closest_pt_child[3];
    // double *_closest_pt_uvP1_child;
    double _closest_pt_uvP1_child[2];
    double _closest_pt_uvPn_child[2];
    double _dist2_child = 0;
    double _closest_pt_weights_child[3];

    PDM_triangle_status_t stat = PDM_triangle_evaluate_position(point_coords,
                                                                _vtx_tria_child,
                                                                _closest_pt_child,
                                                                &_dist2_child,
                                                                _closest_pt_weights_child);


    if (stat == PDM_TRIANGLE_DEGENERATED) {
      continue;
    }

    // _closest_pt_uvP1_child = _closest_pt_weights_child;
    for (int j = 0; j < 2; j++) {
      _closest_pt_uvP1_child[j] = _closest_pt_weights_child[j+1];
      _closest_pt_uvPn_child[j] = 0;
    }
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 3; k++) {
        _closest_pt_uvPn_child[j] +=
          _closest_pt_weights_child[k] * _uvPn_tria_child[2*k + j];
      }
    }

    _heap_s_insert (heap,
                    _vtx_tria_child,
                    _uvPn_tria_child,
                    _closest_pt_child,
                    _closest_pt_uvP1_child,
                    _closest_pt_uvPn_child,
                    _dist2_child,
                    child++);
  }

}


/*-----------------------------------------------------------------------------
 *
 * compute distance from closest triangle subdivision
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   uncertain_result --> 1 if the result is uncertain
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_closest_tria_subdivision
(
 _heap_s_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj,
 int                        *uncertain_result
 )
{
  *uncertain_result = 0;
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;
  double dist2_pre = HUGE_VAL;
  int distance_extension = 0;

  while (1) {

    double _vtx_tria_current[9];
    double _uvPn_tria_current[6];

    double _closest_pt_current[3];
    double _closest_pt_uvP1_current[2];
    double _closest_pt_uvPn_current[2];
    double _dist2_current;

    /* Get closest triangle stored in the heap */

    int _child;
    int is_empty = _heap_s_top_get (heap,
                                    _vtx_tria_current,
                                    _uvPn_tria_current,
                                    _closest_pt_current,
                                    _closest_pt_uvP1_current,
                                    _closest_pt_uvPn_current,
                                    &_dist2_current,
                                    &_child);


    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty %s\n");
      abort();
    }

    if ((distance_extension == 0) && (_dist2_current > dist2_pre)) {
      distance_extension = 1;
    }

    else if (distance_extension == 1) {
      if (_dist2_current <= dist2_min_min) {
        distance_extension = 0;
      }
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);
    dist2_pre = _dist2_current;

    /* Compute projected from current P1 triangle */
    double weightsP1[3];
    PDM_ho_basis (PDM_MESH_NODAL_TRIAHO,
                  1,
                  3,
                  1,
                  _closest_pt_uvP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_tria_current[3*j+k];
      }
    }

    /* Compute projected from current Pn triangle */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uvPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }

    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */
    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uvPn_current[0];
      uvw[1] = _closest_pt_uvPn_current[1];

      break;
    }

    /*
     * Insert sub-triangles in the heap
     */

    _insert_subtria (heap,
                     order,
                     type,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_tria_current,
                     _uvPn_tria_current);

  }


  if (distance_extension) {
    *uncertain_result = 1;
  }

  return dist2;
}


/*-----------------------------------------------------------------------------
 *
 * compute distance from closest triangle subdivision
 *
 * parameters:
 *   heap1            <-> heap
 *   heap2            <-> work heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_uniform_tria_subdivision
(
 _heap_s_t                  *heap1,
 _heap_s_t                  *heap2,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj
 )
{
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;

  _heap_s_t *heap      = heap1;
  _heap_s_t *next_heap = heap2;

  while (1) {

    double _vtx_tria_current[9];
    double _uvPn_tria_current[6];

    double _closest_pt_current[3];
    double _closest_pt_uvP1_current[2];
    double _closest_pt_uvPn_current[2];
    double _dist2_current;

    /* Get closest triangle stored in the heap */

    int _child;
    int is_empty = _heap_s_top_get (heap,
                                    _vtx_tria_current,
                                    _uvPn_tria_current,
                                    _closest_pt_current,
                                    _closest_pt_uvP1_current,
                                    _closest_pt_uvPn_current,
                                    &_dist2_current,
                                    &_child);


    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty %s\n");
      abort();
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);

    /* Compute projected from current P1 triangle */
    double weightsP1[3];
    PDM_ho_basis (PDM_MESH_NODAL_TRIAHO,
                  1,
                  3,
                  1,
                  _closest_pt_uvP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_tria_current[3*j+k];
      }
    }

    /* Compute projected from current Pn triangle */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uvPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }

    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */
    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uvPn_current[0];
      uvw[1] = _closest_pt_uvPn_current[1];

      break;
    }

    /*
     * Insert sub-triangles in the next heap
     */

    _heap_s_init (next_heap);

    _insert_subtria (next_heap,
                     order,
                     type,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_tria_current,
                     _uvPn_tria_current);

    double _vtx_tria_current2[9];
    double _uvPn_tria_current2[6];

    double _closest_pt_current2[3];
    double _closest_pt_uvP1_current2[2];
    double _dist2_current2;
    int _child_current2;

    while ( !_heap_s_top_get (heap,
                              _vtx_tria_current2,
                              _uvPn_tria_current2,
                              _closest_pt_current2,
                              _closest_pt_uvP1_current2,
                              _closest_pt_uvPn_current,
                              &_dist2_current2,
                              &_child_current2)) {


      _insert_subtria (next_heap,
                       order,
                       type,
                       n_nodes,
                       nodes_coords,
                       point_coords,
                       weightsPn,
                       _vtx_tria_current2,
                       _uvPn_tria_current2);

    }

    _heap_s_t *heap_tmp = heap;
    heap = next_heap;
    next_heap = heap_tmp;

  }

  return dist2;
}


/*-----------------------------------------------------------------------------
 *
 * Default point location on a high order triangle
 *
 * parameters:
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   uvw              --> parametric coordinates in the element
 *
 * return:
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double
_default_location_generic_2d
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const double                char_size,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw,
 _heap_fill_init_sub_tria_t  fill_init_fct
 )
{
  const int n_it_max = 100;
  double err_max = PDM_MAX (char_size * 1e-6, 1e-15);

  double dist2 = HUGE_VAL;

  _heap_s_t heap;
  _heap_s_t heap2;

  double *weightsPn = malloc (sizeof(double) * n_nodes);

  /* Initialize heap */

  _heap_s_init (&heap);

  /* Build initial sub-triangles and store them in the heap */

  (fill_init_fct) (&heap,
                   order,
                   n_nodes,
                   nodes_coords,
                   point_coords);

  /*
   *  While error > error_max
   *    - Get closest triangle in the heap
   *    - Cut it in sub-triangles
   *    - Store them in the heap
   */


  const int method = 0;

  int n_it;
  double err_proj = HUGE_VAL;
  int uncertain_result = 0;;

  if (method == 0) {
    dist2 = _compute_dist2_from_closest_tria_subdivision (&heap,
                                                          order,
                                                          type,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj,
                                                          &uncertain_result);

    if (uncertain_result) {

      /* Initialize heap */

      _heap_s_init (&heap);
      _heap_s_init (&heap2);

      /* Build initial sub-triangles and store them in the heap */

      (fill_init_fct) (&heap,
                       order,
                       n_nodes,
                       nodes_coords,
                       point_coords);


      dist2 = _compute_dist2_from_uniform_tria_subdivision (&heap,
                                                            &heap2,
                                                            order,
                                                            type,
                                                            n_nodes,
                                                            n_it_max,
                                                            err_max,
                                                            nodes_coords,
                                                            point_coords,
                                                            weightsPn,
                                                            projected_coords,
                                                            uvw,
                                                            &n_it,
                                                            &err_proj);

    }
  }

  else {

    _heap_s_init (&heap2);
    dist2 = _compute_dist2_from_uniform_tria_subdivision (&heap,
                                                          &heap2,
                                                          order,
                                                          type,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj);

  }

  free (weightsPn);

  return dist2;

}


/*----------------------------------------------------------------------------
 *
 * Compute the radius of a triangle inscribed circle
 *
 * parameters:
 *   coords           <-- coordinates of vertices
 *
 * return:
 *   radius
 *
 *----------------------------------------------------------------------------*/

static double
_radius_inscribed_circle
(
 const double *coords
 )
{
  double a =
    sqrt ((coords[3*1    ] - coords[3*0    ]) * (coords[3*1    ] - coords[3*0    ]) +
          (coords[3*1 + 1] - coords[3*0 + 1]) * (coords[3*1 + 1] - coords[3*0 + 1]) +
          (coords[3*1 + 2] - coords[3*0 + 2]) * (coords[3*1 + 2] - coords[3*0 + 2]));

  double b =
    sqrt ((coords[3*2    ] - coords[3*1    ]) * (coords[3*2    ] - coords[3*1    ]) +
          (coords[3*2 + 1] - coords[3*1 + 1]) * (coords[3*2 + 1] - coords[3*1 + 1]) +
          (coords[3*2 + 2] - coords[3*1 + 2]) * (coords[3*2 + 2] - coords[3*1 + 2]));

  double c =
    sqrt ((coords[3*0    ] - coords[3*2    ]) * (coords[3*0    ] - coords[3*2    ]) +
          (coords[3*0 + 1] - coords[3*2 + 1]) * (coords[3*0 + 1] - coords[3*2 + 1]) +
          (coords[3*0 + 2] - coords[3*2 + 2]) * (coords[3*0 + 2] - coords[3*2 + 2]));

  double p = a + b + c;
  double S = sqrt (p*(p-a)*(p-b)*(p-c));

  return S/p;
}


/*============================================================================
 * VOLUME
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Init heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *
 *----------------------------------------------------------------------------*/

static void
_heap_v_init
(
 _heap_v_t *heap
 )
{
  heap->idx = -1;

  for (int i = 0; i < S_HEAP; i++) {
    heap->free_idx[i] = i;
  }
}


/*----------------------------------------------------------------------------
 *
 * Get top of the heap
 *
 * parameters:
 *   heap             <-> heap to initialize
 *   order            <-> element order
 *   n_node           <-> number of nodes
 *   ho_vertex_num    <-> high order vertex num (internal ordering)
 *   vertex_coords    <-> vertex coordinates
 *   point_coords     <-> point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static int
_heap_v_top_get
(
 _heap_v_t *heap,
 double    *vtx_tetra_current,
 double    *uvwInPn_tetra_current,
 double    *closest_pt_current,
 double    *closest_pt_uvwP1_current,
 double    *closest_pt_uvwInPn_current,
 double    *dist2_current,
 int       *child
 )
{
  if (heap->idx < 0) {
    return 1;
  }

  int idx = heap->sorted_idx[heap->idx];
  heap->free_idx[heap->idx--]= idx;

  double *_vtx_tetra_current = heap->vtx_tetra + 12 * idx;
  for (int i = 0; i < 12; i++) {
    vtx_tetra_current[i] = _vtx_tetra_current[i];
  }

  double *_uvwInPn_tetra_current = heap->uvwInPn_tetra + 12 *idx;
  for (int i = 0; i < 12; i++) {
    uvwInPn_tetra_current[i] = _uvwInPn_tetra_current[i];
  }

  double *_closest_pt_current = heap->closest_pt + 3 *idx;
  for (int i = 0; i < 3; i++) {
    closest_pt_current[i] = _closest_pt_current[i];
  }

  double *_closest_pt_uvwP1_current = heap->closest_pt_uvwP1 + 3 *idx;
  double *_closest_pt_uvwInPn_current = heap->closest_pt_uvwInPn + 3 *idx;
  for (int i = 0; i < 3; i++) {
    closest_pt_uvwP1_current[i] = _closest_pt_uvwP1_current[i];
    closest_pt_uvwInPn_current[i] = _closest_pt_uvwInPn_current[i];
  }

  *child = heap->child[idx];
  *dist2_current = heap->dist2[idx];

  return 0;
}


/*----------------------------------------------------------------------------
 *
 * Add sub-tetrahedron of a pn-tetrahedron in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_v_insert
(
 _heap_v_t *heap,
 double    *vtx_tetra,
 double    *uvwInPn_tetra,
 double    *closest_pt,
 double    *closest_pt_uvwP1,
 double    *closest_pt_uvwInPn,
 double     dist2,
 int        child
 )
{

  // Look for index (binary search)
  int curr_idx = heap->idx;
  int *sorted_idx = heap->sorted_idx;
  double *sorted_dist2 = heap->dist2;

  int beg = 0;
  int end = curr_idx;

  while (beg <= end) {
    double dist2_beg = sorted_dist2[sorted_idx[beg]];
    double dist2_end = sorted_dist2[sorted_idx[end]];
    if (dist2 >= dist2_beg) {
      end = beg - 1;
    }
    else if (dist2 <= dist2_end) {
      beg = end + 1;
    }
    else {

      const int middle = (end + beg) / 2;
      if (beg == middle) {
        beg = beg + 1;
        end = beg - 1;
      }
      else {
        const double dist2_middle = sorted_dist2[sorted_idx[middle]];
        if (dist2 > dist2_middle) {
          end = middle;
        }
        else {
          beg = middle;
        }
      }
    }
  }


  // If the heap is full remove the most distant
  if (curr_idx >= (S_HEAP - 1)) {
    if (beg == 0) {
      return;
    }
    else {
      const int free_idx = sorted_idx[0];

      for (int i = 1; i < S_HEAP; i++) {
        sorted_idx[i-1] = sorted_idx[i];
      }

      heap->free_idx[heap->idx] = free_idx;
      heap->idx--;

      beg = beg - 1;
      end = end - 1;

    }
  }

  // Add the element to the heap
  heap->idx++;
  assert (heap->free_idx[heap->idx] != -1);

  for (int j = heap->idx; j > beg; j--) {
    sorted_idx[j] = sorted_idx[j-1];
  }

  sorted_idx[beg] = heap->free_idx[heap->idx];

  heap->free_idx[heap->idx] = -1;

  int _idx = sorted_idx[beg];

  for (int j = 0; j < 12; j++) {
    heap->vtx_tetra[12*_idx+j] = vtx_tetra[j];
  }

  for (int j = 0; j < 12; j++) {
    heap->uvwInPn_tetra[12*_idx+j] = uvwInPn_tetra[j];
  }

  for (int j = 0; j < 3; j++) {
    heap->closest_pt[3*_idx+j]  = closest_pt[j];
  }

  for (int j = 0; j < 3; j++) {
    heap->closest_pt_uvwP1[3*_idx+j] = closest_pt_uvwP1[j];
    heap->closest_pt_uvwInPn[3*_idx+j] = closest_pt_uvwInPn[j];
  }

  heap->dist2[_idx] = dist2;
  heap->child[_idx] = child;
}



/*----------------------------------------------------------------------------
 *
 *  Extract 4 tetrahedron from an octahedron
 *
 *----------------------------------------------------------------------------*/

static void _octa_to_tetra
(
 double *vertex_octa,
 double *vertex_tetra,
 double *uvwPn_octa,
 double *uvw_tetra
 )
{

  // 1st tetra
  vertex_tetra[0]  = vertex_octa[0];  uvw_tetra[0]  = uvwPn_octa[0];
  vertex_tetra[1]  = vertex_octa[1];  uvw_tetra[1]  = uvwPn_octa[1];
  vertex_tetra[2]  = vertex_octa[2];  uvw_tetra[2]  = uvwPn_octa[2];

  vertex_tetra[3]  = vertex_octa[3];  uvw_tetra[3]  = uvwPn_octa[3];
  vertex_tetra[4]  = vertex_octa[4];  uvw_tetra[4]  = uvwPn_octa[4];
  vertex_tetra[5]  = vertex_octa[5];  uvw_tetra[5]  = uvwPn_octa[5];

  vertex_tetra[6]  = vertex_octa[6];  uvw_tetra[6]  = uvwPn_octa[6];
  vertex_tetra[7]  = vertex_octa[7];  uvw_tetra[7]  = uvwPn_octa[7];
  vertex_tetra[8]  = vertex_octa[8];  uvw_tetra[8]  = uvwPn_octa[8];

  vertex_tetra[9]  = vertex_octa[15]; uvw_tetra[9]  = uvwPn_octa[15];
  vertex_tetra[10] = vertex_octa[16]; uvw_tetra[10] = uvwPn_octa[16];
  vertex_tetra[11] = vertex_octa[17]; uvw_tetra[11] = uvwPn_octa[17];

  // 2nd tetra
  vertex_tetra[12] = vertex_octa[0];  uvw_tetra[12] = uvwPn_octa[0];
  vertex_tetra[13] = vertex_octa[1];  uvw_tetra[13] = uvwPn_octa[1];
  vertex_tetra[14] = vertex_octa[2];  uvw_tetra[14] = uvwPn_octa[2];

  vertex_tetra[15] = vertex_octa[6];  uvw_tetra[15] = uvwPn_octa[6];
  vertex_tetra[16] = vertex_octa[7];  uvw_tetra[16] = uvwPn_octa[7];
  vertex_tetra[17] = vertex_octa[8];  uvw_tetra[17] = uvwPn_octa[8];

  vertex_tetra[18] = vertex_octa[9];  uvw_tetra[18] = uvwPn_octa[9];
  vertex_tetra[19] = vertex_octa[10]; uvw_tetra[19] = uvwPn_octa[10];
  vertex_tetra[20] = vertex_octa[11]; uvw_tetra[20] = uvwPn_octa[11];

  vertex_tetra[21] = vertex_octa[15]; uvw_tetra[21] = uvwPn_octa[15];
  vertex_tetra[22] = vertex_octa[16]; uvw_tetra[22] = uvwPn_octa[16];
  vertex_tetra[23] = vertex_octa[17]; uvw_tetra[23] = uvwPn_octa[17];

  // 3rd tetra
  vertex_tetra[24] = vertex_octa[0];  uvw_tetra[24] = uvwPn_octa[0];
  vertex_tetra[25] = vertex_octa[1];  uvw_tetra[25] = uvwPn_octa[1];
  vertex_tetra[26] = vertex_octa[2];  uvw_tetra[26] = uvwPn_octa[2];

  vertex_tetra[27] = vertex_octa[9];  uvw_tetra[27] = uvwPn_octa[9];
  vertex_tetra[28] = vertex_octa[10]; uvw_tetra[28] = uvwPn_octa[10];
  vertex_tetra[29] = vertex_octa[11]; uvw_tetra[29] = uvwPn_octa[11];

  vertex_tetra[30] = vertex_octa[12]; uvw_tetra[30] = uvwPn_octa[12];
  vertex_tetra[31] = vertex_octa[13]; uvw_tetra[31] = uvwPn_octa[13];
  vertex_tetra[32] = vertex_octa[14]; uvw_tetra[32] = uvwPn_octa[14];

  vertex_tetra[33] = vertex_octa[15]; uvw_tetra[33] = uvwPn_octa[15];
  vertex_tetra[34] = vertex_octa[16]; uvw_tetra[34] = uvwPn_octa[16];
  vertex_tetra[35] = vertex_octa[17]; uvw_tetra[35] = uvwPn_octa[17];

  // 4th tetra
  vertex_tetra[36] = vertex_octa[0];  uvw_tetra[36] = uvwPn_octa[0];
  vertex_tetra[37] = vertex_octa[1];  uvw_tetra[37] = uvwPn_octa[1];
  vertex_tetra[38] = vertex_octa[2];  uvw_tetra[38] = uvwPn_octa[2];

  vertex_tetra[39] = vertex_octa[3];  uvw_tetra[39] = uvwPn_octa[3];
  vertex_tetra[40] = vertex_octa[4];  uvw_tetra[40] = uvwPn_octa[4];
  vertex_tetra[41] = vertex_octa[5];  uvw_tetra[41] = uvwPn_octa[5];

  vertex_tetra[42] = vertex_octa[12]; uvw_tetra[42] = uvwPn_octa[12];
  vertex_tetra[43] = vertex_octa[13]; uvw_tetra[43] = uvwPn_octa[13];
  vertex_tetra[44] = vertex_octa[14]; uvw_tetra[44] = uvwPn_octa[14];

  vertex_tetra[45] = vertex_octa[15]; uvw_tetra[45] = uvwPn_octa[15];
  vertex_tetra[46] = vertex_octa[16]; uvw_tetra[46] = uvwPn_octa[16];
  vertex_tetra[47] = vertex_octa[17]; uvw_tetra[47] = uvwPn_octa[17];
}


/*----------------------------------------------------------------------------
 *
 * Add sub-tetrahedron of a pn-tetrahedron in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_tetra_sub_tetra
(
 _heap_v_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{

  double *uvw_nodes   = malloc (sizeof(double) * 3 * n_nodes);
  double _vertex_tetra[48];
  double _uvw_vertex_tetra[48];
  _uvw_ho_tetra_nodes (order,
                       0., 1.,
                       0., 1.,
                       0., 1.,
                       uvw_nodes);

  int step = 0;
  int ibeg = 0;
  int child = 0;
  for (int k = 0; k < order; k++) {
    int n_nodes_basis = (order+1-k)*(order+2-k)/2;
    for (int j = 0; j < order - k; j++) {
      step = order - j - k + 1;
      for (int i = ibeg; i < ibeg + step - 1; i++) {
        int i1 = i+1;

        int idx1 = i;
        int idx2 = i1;
        int idx3 = i + step;
        int idx4 = i + n_nodes_basis - j;

        double x1 = nodes_coords[3*idx1];
        double y1 = nodes_coords[3*idx1 + 1];
        double z1 = nodes_coords[3*idx1 + 2];

        double x2 = nodes_coords[3*idx2];
        double y2 = nodes_coords[3*idx2 + 1];
        double z2 = nodes_coords[3*idx2 + 2];

        double x3 = nodes_coords[3*idx3];
        double y3 = nodes_coords[3*idx3 + 1];
        double z3 = nodes_coords[3*idx3 + 2];

        double x4 = nodes_coords[3*idx4];
        double y4 = nodes_coords[3*idx4 + 1];
        double z4 = nodes_coords[3*idx4 + 2];

        double __vertex_coords[12] = {x1, y1, z1,
                                      x2, y2, z2,
                                      x3, y3, z3,
                                      x4, y4, z4};

        double _uvwPn_sub_tetra[12];
        _uvwPn_sub_tetra[0]  = uvw_nodes[3*idx1];
        _uvwPn_sub_tetra[1]  = uvw_nodes[3*idx1+1];
        _uvwPn_sub_tetra[2]  = uvw_nodes[3*idx1+2];
        _uvwPn_sub_tetra[3]  = uvw_nodes[3*idx2];
        _uvwPn_sub_tetra[4]  = uvw_nodes[3*idx2+1];
        _uvwPn_sub_tetra[5]  = uvw_nodes[3*idx2+2];
        _uvwPn_sub_tetra[6]  = uvw_nodes[3*idx3];
        _uvwPn_sub_tetra[7]  = uvw_nodes[3*idx3+1];
        _uvwPn_sub_tetra[8]  = uvw_nodes[3*idx3+2];
        _uvwPn_sub_tetra[9]  = uvw_nodes[3*idx4];
        _uvwPn_sub_tetra[10] = uvw_nodes[3*idx4+1];
        _uvwPn_sub_tetra[11] = uvw_nodes[3*idx4+2];


        double _closest_point_P1[3];
        double _weights_closest_point_P1[4];
        double _dist2;

        PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                               __vertex_coords,
                                                               _closest_point_P1,
                                                               &_dist2,
                                                               _weights_closest_point_P1);


        if (status != PDM_TETRAHEDRON_DEGENERATED) {

          double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
          double _uvw_closest_point_Pn[3] = {0., 0., 0.};
          for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 4; j2++) {
              _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
            }
          }


          _heap_v_insert (heap,
                          __vertex_coords,
                          _uvwPn_sub_tetra,
                          _closest_point_P1,
                          _uvw_closest_point_P1,
                          _uvw_closest_point_Pn,
                          _dist2,
                          child++);

        }

        if ( j != 0) {
          idx1 = i;
          idx2 = i1;
          idx3 = i - step;
          idx4 = i - step + n_nodes_basis - j;
          int idx5 = i + n_nodes_basis - j;
          int idx6 = i - step + n_nodes_basis - j + 1;


          x1 = nodes_coords[3*idx1];
          y1 = nodes_coords[3*idx1 + 1];
          z1 = nodes_coords[3*idx1 + 2];

          x2 = nodes_coords[3*idx2];
          y2 = nodes_coords[3*idx2 + 1];
          z2 = nodes_coords[3*idx2 + 2];

          x3 = nodes_coords[3*idx3];
          y3 = nodes_coords[3*idx3 + 1];
          z3 = nodes_coords[3*idx3 + 2];

          x4 = nodes_coords[3*idx4];
          y4 = nodes_coords[3*idx4 + 1];
          z4 = nodes_coords[3*idx4 + 2];

          double x5 = nodes_coords[3*idx5];
          double y5 = nodes_coords[3*idx5 + 1];
          double z5 = nodes_coords[3*idx5 + 2];

          double x6 = nodes_coords[3*idx6];
          double y6 = nodes_coords[3*idx6 + 1];
          double z6 = nodes_coords[3*idx6 + 2];


          double _vertex_octa[18] = {x1, y1, z1,
                                     x2, y2, z2,
                                     x3, y3, z3,
                                     x4, y4, z4,
                                     x5, y5, z5,
                                     x6, y6, z6};

          double _uvwPn_sub_octa[18];
          _uvwPn_sub_octa[0]  = uvw_nodes[3*idx1];
          _uvwPn_sub_octa[1]  = uvw_nodes[3*idx1+1];
          _uvwPn_sub_octa[2]  = uvw_nodes[3*idx1+2];
          _uvwPn_sub_octa[3]  = uvw_nodes[3*idx2];
          _uvwPn_sub_octa[4]  = uvw_nodes[3*idx2+1];
          _uvwPn_sub_octa[5]  = uvw_nodes[3*idx2+2];
          _uvwPn_sub_octa[6]  = uvw_nodes[3*idx3];
          _uvwPn_sub_octa[7]  = uvw_nodes[3*idx3+1];
          _uvwPn_sub_octa[8]  = uvw_nodes[3*idx3+2];
          _uvwPn_sub_octa[9]  = uvw_nodes[3*idx4];
          _uvwPn_sub_octa[10] = uvw_nodes[3*idx4+1];
          _uvwPn_sub_octa[11] = uvw_nodes[3*idx4+2];
          _uvwPn_sub_octa[12] = uvw_nodes[3*idx5];
          _uvwPn_sub_octa[13] = uvw_nodes[3*idx5+1];
          _uvwPn_sub_octa[14] = uvw_nodes[3*idx5+2];
          _uvwPn_sub_octa[15] = uvw_nodes[3*idx6];
          _uvwPn_sub_octa[16] = uvw_nodes[3*idx6+1];
          _uvwPn_sub_octa[17] = uvw_nodes[3*idx6+2];

          _octa_to_tetra (_vertex_octa,
                          _vertex_tetra,
                          _uvwPn_sub_octa,
                          _uvw_vertex_tetra);

          for (int itetra = 0; itetra < 4; itetra++) {

            for (int n1 = 0; n1 < 12; n1++) {
              const double *__vertex_tetra = _vertex_tetra + 12*itetra + n1;
              const double *__uvw_vertex_tetra = _uvw_vertex_tetra + 12*itetra + n1;
              __vertex_coords[n1] = *__vertex_tetra;
              _uvwPn_sub_tetra[n1] = *__uvw_vertex_tetra;
            }


            status  = PDM_tetrahedron_evaluate_position (point_coords,
                                                         __vertex_coords,
                                                         _closest_point_P1,
                                                         &_dist2,
                                                         _weights_closest_point_P1);

            if (status != PDM_TETRAHEDRON_DEGENERATED) {

              double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
              double _uvw_closest_point_Pn[3] = {0., 0., 0.};
              for (int j1 = 0; j1 < 3; j1++) {
                for (int j2 = 0; j2 < 4; j2++) {
                  _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
                }
              }
              _heap_v_insert (heap,
                              __vertex_coords,
                              _uvwPn_sub_tetra,
                              _closest_point_P1,
                              _uvw_closest_point_P1,
                              _uvw_closest_point_Pn,
                              _dist2,
                              child++);
            }
          }


        }

      }
      ibeg += step;
    }
    ibeg ++;
  }

  free (uvw_nodes);
}



/*----------------------------------------------------------------------------
 *
 *  Extract 2 tetrahedron from a pyramid
 *
 *----------------------------------------------------------------------------*/

static void _pyra_to_tetra
(
 double *vertex_pyra,
 double *vertex_tetra,
 double *uvwPn_pyra,
 double *uvw_tetra
 )
{

  // 1st tetra
  vertex_tetra[0]  = vertex_pyra[0];  uvw_tetra[0]  = uvwPn_pyra[0];
  vertex_tetra[1]  = vertex_pyra[1];  uvw_tetra[1]  = uvwPn_pyra[1];
  vertex_tetra[2]  = vertex_pyra[2];  uvw_tetra[2]  = uvwPn_pyra[2];

  vertex_tetra[3]  = vertex_pyra[3];  uvw_tetra[3]  = uvwPn_pyra[3];
  vertex_tetra[4]  = vertex_pyra[4];  uvw_tetra[4]  = uvwPn_pyra[4];
  vertex_tetra[5]  = vertex_pyra[5];  uvw_tetra[5]  = uvwPn_pyra[5];

  vertex_tetra[6]  = vertex_pyra[9];  uvw_tetra[6]  = uvwPn_pyra[9];
  vertex_tetra[7]  = vertex_pyra[10]; uvw_tetra[7]  = uvwPn_pyra[10];
  vertex_tetra[8]  = vertex_pyra[11]; uvw_tetra[8]  = uvwPn_pyra[11];

  vertex_tetra[9]  = vertex_pyra[12]; uvw_tetra[9]  = uvwPn_pyra[12];
  vertex_tetra[10] = vertex_pyra[13]; uvw_tetra[10] = uvwPn_pyra[13];
  vertex_tetra[11] = vertex_pyra[14]; uvw_tetra[11] = uvwPn_pyra[14];

  // 2nd tetra
  vertex_tetra[12] = vertex_pyra[3];  uvw_tetra[12] = uvwPn_pyra[3];
  vertex_tetra[13] = vertex_pyra[4];  uvw_tetra[13] = uvwPn_pyra[4];
  vertex_tetra[14] = vertex_pyra[5];  uvw_tetra[14] = uvwPn_pyra[5];

  vertex_tetra[15] = vertex_pyra[6];  uvw_tetra[15] = uvwPn_pyra[6];
  vertex_tetra[16] = vertex_pyra[7];  uvw_tetra[16] = uvwPn_pyra[7];
  vertex_tetra[17] = vertex_pyra[8];  uvw_tetra[17] = uvwPn_pyra[8];

  vertex_tetra[18] = vertex_pyra[9];  uvw_tetra[18] = uvwPn_pyra[9];
  vertex_tetra[19] = vertex_pyra[10]; uvw_tetra[19] = uvwPn_pyra[10];
  vertex_tetra[20] = vertex_pyra[11]; uvw_tetra[20] = uvwPn_pyra[11];

  vertex_tetra[21] = vertex_pyra[12]; uvw_tetra[21] = uvwPn_pyra[12];
  vertex_tetra[22] = vertex_pyra[13]; uvw_tetra[22] = uvwPn_pyra[13];
  vertex_tetra[23] = vertex_pyra[14]; uvw_tetra[23] = uvwPn_pyra[14];

}


/*----------------------------------------------------------------------------
 *
 * Add sub-tetrahedron of a pn-pyramid in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_pyra_sub_tetra
(
 _heap_v_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{

  double *uvw_nodes = malloc (sizeof(double) * 3 * n_nodes);
  double _vertex_tetra[24];
  double _uvw_vertex_tetra[24];

  _uvw_ho_pyra_nodes (order,
                      0., 1.,
                      0., 1.,
                      0., 1.,
                      uvw_nodes);


  int step = 0;
  int ibeg = 0;
  int child = 0;
  for (int k = 0; k < order; k++) {
    step = order - k + 1;
    int n_nodes_basis = step*step;
    for (int j = 0; j < order - k; j++) {
      for (int i = ibeg; i < ibeg + step - 1; i++) {

        int idx1 = 3*i;
        int idx2 = 3*(i + 1);
        int idx3 = 3*(i + 1 + step);
        int idx4 = 3*(i + step);
        int idx5 = 3*(i + n_nodes_basis - j);

        double _vertex_pyra[15] = {
          nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
          nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
          nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
          nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2],
          nodes_coords[idx5], nodes_coords[idx5 + 1], nodes_coords[idx5 + 2]
        };

        double _uvwPn_sub_pyra[15] = {
          uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
          uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
          uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
          uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2],
          uvw_nodes[idx5], uvw_nodes[idx5 + 1], uvw_nodes[idx5 + 2]
        };

        _pyra_to_tetra (_vertex_pyra,
                        _vertex_tetra,
                        _uvwPn_sub_pyra,
                        _uvw_vertex_tetra);

        for (int itetra = 0; itetra < 2; itetra++) {

          double *__vertex_coords = _vertex_tetra + 12*itetra;
          double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

          double _closest_point_P1[3];
          double _weights_closest_point_P1[4];
          double _dist2;

          PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                                 __vertex_coords,
                                                                 _closest_point_P1,
                                                                 &_dist2,
                                                                 _weights_closest_point_P1);


          if (status != PDM_TETRAHEDRON_DEGENERATED) {

            double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
            double _uvw_closest_point_Pn[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_coords,
                            _uvwPn_sub_tetra,
                            _closest_point_P1,
                            _uvw_closest_point_P1,
                            _uvw_closest_point_Pn,
                            _dist2,
                            child++);

          }
        }

        if (j != 0 && i != ibeg) { //les points qui ne sont pas sur une face triangulaire de la pyramide
          // 1st tetra
          idx1 = 3*i;
          idx2 = 3*(i - 1);
          idx3 = 3*(i + n_nodes_basis - j - 1);
          idx4 = 3*(i + n_nodes_basis - j - step);

          double __vertex_tetra1[12] = {
            nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
            nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
            nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
            nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2]
          };

          double _uvwPn_sub_tetra1[12] = {
            uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
            uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
            uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
            uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2]
          };


          double _closest_point_P1_tetra1[3];
          double _weights_closest_point_P1_tetra1[4];
          double _dist2_tetra1;

          int is_degenerate_tetra1 = PDM_tetrahedron_evaluate_position (point_coords,
                                                                        __vertex_tetra1,
                                                                        _closest_point_P1_tetra1,
                                                                        &_dist2_tetra1,
                                                                        _weights_closest_point_P1_tetra1);


          if (is_degenerate_tetra1 != -1) {

            double *_uvw_closest_point_P1_tetra1 = _weights_closest_point_P1_tetra1 + 1;
            double _uvw_closest_point_Pn_tetra1[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn_tetra1[j1] += _weights_closest_point_P1_tetra1[j2] * _uvwPn_sub_tetra1[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_tetra1,
                            _uvwPn_sub_tetra1,
                            _closest_point_P1_tetra1,
                            _uvw_closest_point_P1_tetra1,
                            _uvw_closest_point_Pn_tetra1,
                            _dist2_tetra1,
                            child++);
          }


          // 2nd tetra
          idx1 = 3*i;
          idx2 = 3*(i - step);
          idx3 = 3*(i + n_nodes_basis - j - step);
          idx4 = 3*(i + n_nodes_basis - j - step + 1);

          double __vertex_tetra2[12] = {
            nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
            nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
            nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
            nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2]
          };

          double _uvwPn_sub_tetra2[12] = {
            uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
            uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
            uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
            uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2]
          };


          double _closest_point_P1_tetra2[3];
          double _weights_closest_point_P1_tetra2[4];
          double _dist2_tetra2;

          int is_degenerate_tetra2 = PDM_tetrahedron_evaluate_position (point_coords,
                                                                        __vertex_tetra2,
                                                                        _closest_point_P1_tetra2,
                                                                        &_dist2_tetra2,
                                                                        _weights_closest_point_P1_tetra2);

          if (is_degenerate_tetra2 != -1) {

            double *_uvw_closest_point_P1_tetra2 = _weights_closest_point_P1_tetra2 + 1;
            double _uvw_closest_point_Pn_tetra2[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn_tetra2[j1] += _weights_closest_point_P1_tetra2[j2] * _uvwPn_sub_tetra2[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_tetra2,
                            _uvwPn_sub_tetra2,
                            _closest_point_P1_tetra2,
                            _uvw_closest_point_P1_tetra2,
                            _uvw_closest_point_Pn_tetra2,
                            _dist2_tetra2,
                            child++);
          }


          // 3rd tetra
          idx1 = 3*i;
          idx2 = 3*(i + 1);
          idx3 = 3*(i + n_nodes_basis - j - step + 1);
          idx4 = 3*(i + n_nodes_basis - j);

          double __vertex_tetra3[12] = {
            nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
            nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
            nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
            nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2]
          };

          double _uvwPn_sub_tetra3[12] = {
            uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
            uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
            uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
            uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2]
          };


          double _closest_point_P1_tetra3[3];
          double _weights_closest_point_P1_tetra3[4];
          double _dist2_tetra3;

          int is_degenerate_tetra3 = PDM_tetrahedron_evaluate_position (point_coords,
                                                                        __vertex_tetra3,
                                                                        _closest_point_P1_tetra3,
                                                                        &_dist2_tetra3,
                                                                        _weights_closest_point_P1_tetra3);

          if (is_degenerate_tetra3 != -1) {

            double *_uvw_closest_point_P1_tetra3 = _weights_closest_point_P1_tetra3 + 1;
            double _uvw_closest_point_Pn_tetra3[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn_tetra3[j1] += _weights_closest_point_P1_tetra3[j2] * _uvwPn_sub_tetra3[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_tetra3,
                            _uvwPn_sub_tetra3,
                            _closest_point_P1_tetra3,
                            _uvw_closest_point_P1_tetra3,
                            _uvw_closest_point_Pn_tetra3,
                            _dist2_tetra3,
                            child++);
          }


          // 4th tetra
          idx1 = 3*i;
          idx2 = 3*(i + step);
          idx3 = 3*(i + n_nodes_basis - j);
          idx4 = 3*(i + n_nodes_basis - j - 1);

          double __vertex_tetra4[12] = {
            nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
            nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
            nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
            nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2]
          };

          double _uvwPn_sub_tetra4[12] = {
            uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
            uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
            uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
            uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2]
          };


          double _closest_point_P1_tetra4[3];
          double _weights_closest_point_P1_tetra4[4];
          double _dist2_tetra4;

          int is_degenerate_tetra4 = PDM_tetrahedron_evaluate_position (point_coords,
                                                                        __vertex_tetra4,
                                                                        _closest_point_P1_tetra4,
                                                                        &_dist2_tetra4,
                                                                        _weights_closest_point_P1_tetra4);

          if (is_degenerate_tetra4 != -1) {

            double *_uvw_closest_point_P1_tetra4 = _weights_closest_point_P1_tetra4 + 1;
            double _uvw_closest_point_Pn_tetra4[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn_tetra4[j1] += _weights_closest_point_P1_tetra4[j2] * _uvwPn_sub_tetra4[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_tetra4,
                            _uvwPn_sub_tetra4,
                            _closest_point_P1_tetra4,
                            _uvw_closest_point_P1_tetra4,
                            _uvw_closest_point_Pn_tetra4,
                            _dist2_tetra4,
                            child++);
          }


          // inverted pyramid
          idx1 = 3*(i + n_nodes_basis - j - step);
          idx2 = 3*(i + n_nodes_basis - j - 1);
          idx3 = 3*(i + n_nodes_basis - j);
          idx4 = 3*(i + n_nodes_basis - j - step + 1);
          idx5 = 3*i;

          double _vertex_pyrai[15] = {
            nodes_coords[idx1], nodes_coords[idx1 + 1], nodes_coords[idx1 + 2],
            nodes_coords[idx2], nodes_coords[idx2 + 1], nodes_coords[idx2 + 2],
            nodes_coords[idx3], nodes_coords[idx3 + 1], nodes_coords[idx3 + 2],
            nodes_coords[idx4], nodes_coords[idx4 + 1], nodes_coords[idx4 + 2],
            nodes_coords[idx5], nodes_coords[idx5 + 1], nodes_coords[idx5 + 2]
          };

          double _uvwPn_sub_pyrai[15] = {
            uvw_nodes[idx1], uvw_nodes[idx1 + 1], uvw_nodes[idx1 + 2],
            uvw_nodes[idx2], uvw_nodes[idx2 + 1], uvw_nodes[idx2 + 2],
            uvw_nodes[idx3], uvw_nodes[idx3 + 1], uvw_nodes[idx3 + 2],
            uvw_nodes[idx4], uvw_nodes[idx4 + 1], uvw_nodes[idx4 + 2],
            uvw_nodes[idx5], uvw_nodes[idx5 + 1], uvw_nodes[idx5 + 2]
          };

          _pyra_to_tetra(_vertex_pyrai,
                         _vertex_tetra,
                         _uvwPn_sub_pyrai,
                         _uvw_vertex_tetra);

          for (int itetra = 0; itetra < 2; itetra++) {

            double *__vertex_coords = _vertex_tetra + 12*itetra;
            double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

            double _closest_point_P1[3];
            double _weights_closest_point_P1[4];
            double _dist2;

            PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                                   __vertex_coords,
                                                                   _closest_point_P1,
                                                                   &_dist2,
                                                                   _weights_closest_point_P1);

            if (status != PDM_TETRAHEDRON_DEGENERATED) {

              double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
              double _uvw_closest_point_Pn[3] = {0., 0., 0.};
              for (int j1 = 0; j1 < 3; j1++) {
                for (int j2 = 0; j2 < 4; j2++) {
                  _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
                }
              }

              _heap_v_insert (heap,
                              __vertex_coords,
                              _uvwPn_sub_tetra,
                              _closest_point_P1,
                              _uvw_closest_point_P1,
                              _uvw_closest_point_Pn,
                              _dist2,
                              child++);
            }
          }
        }

      }
      ibeg += step;
    }
    ibeg += step;
  }

  free(uvw_nodes);
}



/*----------------------------------------------------------------------------
 *
 *  Extract 3 tetrahedron from a prism
 *
 *----------------------------------------------------------------------------*/

static void _prism_to_tetra
(
 double *vertex_prism,
 double *vertex_tetra,
 double *uvwPn_prism,
 double *uvw_tetra
 )
{
  // 1st tetra
  vertex_tetra[0]  = vertex_prism[0];  uvw_tetra[0]  = uvwPn_prism[0];
  vertex_tetra[1]  = vertex_prism[1];  uvw_tetra[1]  = uvwPn_prism[1];
  vertex_tetra[2]  = vertex_prism[2];  uvw_tetra[2]  = uvwPn_prism[2];

  vertex_tetra[3]  = vertex_prism[3];  uvw_tetra[3]  = uvwPn_prism[3];
  vertex_tetra[4]  = vertex_prism[4];  uvw_tetra[4]  = uvwPn_prism[4];
  vertex_tetra[5]  = vertex_prism[5];  uvw_tetra[5]  = uvwPn_prism[5];

  vertex_tetra[6]  = vertex_prism[6];  uvw_tetra[6]  = uvwPn_prism[6];
  vertex_tetra[7]  = vertex_prism[7];  uvw_tetra[7]  = uvwPn_prism[7];
  vertex_tetra[8]  = vertex_prism[8];  uvw_tetra[8]  = uvwPn_prism[8];

  vertex_tetra[9]  = vertex_prism[9];  uvw_tetra[9]  = uvwPn_prism[9];
  vertex_tetra[10] = vertex_prism[10]; uvw_tetra[10] = uvwPn_prism[10];
  vertex_tetra[11] = vertex_prism[11]; uvw_tetra[11] = uvwPn_prism[11];

  // 2nd tetra
  vertex_tetra[12] = vertex_prism[3];  uvw_tetra[12] = uvwPn_prism[3];
  vertex_tetra[13] = vertex_prism[4];  uvw_tetra[13] = uvwPn_prism[4];
  vertex_tetra[14] = vertex_prism[5];  uvw_tetra[14] = uvwPn_prism[5];

  vertex_tetra[15] = vertex_prism[6];  uvw_tetra[15] = uvwPn_prism[6];
  vertex_tetra[16] = vertex_prism[7];  uvw_tetra[16] = uvwPn_prism[7];
  vertex_tetra[17] = vertex_prism[8];  uvw_tetra[17] = uvwPn_prism[8];

  vertex_tetra[18] = vertex_prism[9];  uvw_tetra[18] = uvwPn_prism[9];
  vertex_tetra[19] = vertex_prism[10]; uvw_tetra[19] = uvwPn_prism[10];
  vertex_tetra[20] = vertex_prism[11]; uvw_tetra[20] = uvwPn_prism[11];

  vertex_tetra[21] = vertex_prism[12]; uvw_tetra[21] = uvwPn_prism[12];
  vertex_tetra[22] = vertex_prism[13]; uvw_tetra[22] = uvwPn_prism[13];
  vertex_tetra[23] = vertex_prism[14]; uvw_tetra[23] = uvwPn_prism[14];

  // 3rd tetra
  vertex_tetra[24] = vertex_prism[6];  uvw_tetra[24] = uvwPn_prism[6];
  vertex_tetra[25] = vertex_prism[7];  uvw_tetra[25] = uvwPn_prism[7];
  vertex_tetra[26] = vertex_prism[8];  uvw_tetra[26] = uvwPn_prism[8];

  vertex_tetra[27] = vertex_prism[9];  uvw_tetra[27] = uvwPn_prism[9];
  vertex_tetra[28] = vertex_prism[10]; uvw_tetra[28] = uvwPn_prism[10];
  vertex_tetra[29] = vertex_prism[11]; uvw_tetra[29] = uvwPn_prism[11];

  vertex_tetra[30] = vertex_prism[12]; uvw_tetra[30] = uvwPn_prism[12];
  vertex_tetra[31] = vertex_prism[13]; uvw_tetra[31] = uvwPn_prism[13];
  vertex_tetra[32] = vertex_prism[14]; uvw_tetra[32] = uvwPn_prism[14];

  vertex_tetra[33] = vertex_prism[15]; uvw_tetra[33] = uvwPn_prism[15];
  vertex_tetra[34] = vertex_prism[16]; uvw_tetra[34] = uvwPn_prism[16];
  vertex_tetra[35] = vertex_prism[17]; uvw_tetra[35] = uvwPn_prism[17];
}


/*----------------------------------------------------------------------------
 *
 * Add sub-tetrahedron of a pn-prism in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_prism_sub_tetra
(
 _heap_v_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{
  int ibeg = 0;
  int iend = order;
  int k1;
  int n_nodes_basis = (order + 1)*(order + 2) / 2;

  double *uvw_nodes = malloc (sizeof(double) * 3 * n_nodes);
  double _vertex_tetra[36];
  double _uvw_vertex_tetra[36];

  _uvw_ho_prism_nodes (order,
                       0., 1.,
                       0., 1.,
                       0., 1.,
                       uvw_nodes);


  int child = 0;
  for (int k = 0; k < order; k++) {
    for (int j = 0; j < order; j++) {
      k1 = 0;
      for (int i = ibeg; i < iend - 1; i++) {
        int idx1 = i;
        int idx2 = i+1;
        int idx3 = iend + 1 + k1;
        int idx4 = iend + 2 + k1;
        int idx5 = n_nodes_basis + i;
        int idx6 = n_nodes_basis + i+1;
        int idx7 = n_nodes_basis + iend + 1 + k1;
        int idx8 = n_nodes_basis + iend + 2 + k1;

        double x1 = nodes_coords[3*idx1];
        double y1 = nodes_coords[3*idx1 + 1];
        double z1 = nodes_coords[3*idx1 + 2];

        double x2 = nodes_coords[3*idx2];
        double y2 = nodes_coords[3*idx2 + 1];
        double z2 = nodes_coords[3*idx2 + 2];

        double x3 = nodes_coords[3*idx3];
        double y3 = nodes_coords[3*idx3 + 1];
        double z3 = nodes_coords[3*idx3 + 2];

        double x4 = nodes_coords[3*idx4];
        double y4 = nodes_coords[3*idx4 + 1];
        double z4 = nodes_coords[3*idx4 + 2];

        double x5 = nodes_coords[3*idx5];
        double y5 = nodes_coords[3*idx5 + 1];
        double z5 = nodes_coords[3*idx5 + 2];

        double x6 = nodes_coords[3*idx6];
        double y6 = nodes_coords[3*idx6 + 1];
        double z6 = nodes_coords[3*idx6 + 2];

        double x7 = nodes_coords[3*idx7];
        double y7 = nodes_coords[3*idx7 + 1];
        double z7 = nodes_coords[3*idx7 + 2];

        double x8 = nodes_coords[3*idx8];
        double y8 = nodes_coords[3*idx8 + 1];
        double z8 = nodes_coords[3*idx8 + 2];

        double _vertex_prism[18] = {x1, y1, z1,
                                    x2, y2, z2,
                                    x3, y3, z3,
                                    x5, y5, z5,
                                    x6, y6, z6,
                                    x7, y7, z7};


        double _uvwPn_sub_prism[18];
        _uvwPn_sub_prism[0]  = uvw_nodes[3*idx1];
        _uvwPn_sub_prism[1]  = uvw_nodes[3*idx1+1];
        _uvwPn_sub_prism[2]  = uvw_nodes[3*idx1+2];
        _uvwPn_sub_prism[3]  = uvw_nodes[3*idx2];
        _uvwPn_sub_prism[4]  = uvw_nodes[3*idx2+1];
        _uvwPn_sub_prism[5]  = uvw_nodes[3*idx2+2];
        _uvwPn_sub_prism[6]  = uvw_nodes[3*idx3];
        _uvwPn_sub_prism[7]  = uvw_nodes[3*idx3+1];
        _uvwPn_sub_prism[8]  = uvw_nodes[3*idx3+2];
        _uvwPn_sub_prism[9]  = uvw_nodes[3*idx5];
        _uvwPn_sub_prism[10] = uvw_nodes[3*idx5+1];
        _uvwPn_sub_prism[11] = uvw_nodes[3*idx5+2];
        _uvwPn_sub_prism[12] = uvw_nodes[3*idx6];
        _uvwPn_sub_prism[13] = uvw_nodes[3*idx6+1];
        _uvwPn_sub_prism[14] = uvw_nodes[3*idx6+2];
        _uvwPn_sub_prism[15] = uvw_nodes[3*idx7];
        _uvwPn_sub_prism[16] = uvw_nodes[3*idx7+1];
        _uvwPn_sub_prism[17] = uvw_nodes[3*idx7+2];

        _prism_to_tetra (_vertex_prism,
                         _vertex_tetra,
                         _uvwPn_sub_prism,
                         _uvw_vertex_tetra);

        for (int itetra = 0; itetra < 3; itetra++) {

          double *__vertex_coords = _vertex_tetra + 12*itetra;
          double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

          double _closest_point_P1[3];
          double _weights_closest_point_P1[4];
          double _dist2;

          PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                                 __vertex_coords,
                                                                 _closest_point_P1,
                                                                 &_dist2,
                                                                 _weights_closest_point_P1);

          if (status != PDM_TETRAHEDRON_DEGENERATED) {

            double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
            double _uvw_closest_point_Pn[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
              }
            }


            _heap_v_insert (heap,
                            __vertex_coords,
                            _uvwPn_sub_tetra,
                            _closest_point_P1,
                            _uvw_closest_point_P1,
                            _uvw_closest_point_Pn,
                            _dist2,
                            child++);

          }

        }

        _vertex_prism[0]  = x2;
        _vertex_prism[1]  = y2;
        _vertex_prism[2]  = z2;
        _vertex_prism[3]  = x3;
        _vertex_prism[4]  = y3;
        _vertex_prism[5]  = z3;
        _vertex_prism[6]  = x4;
        _vertex_prism[7]  = y4;
        _vertex_prism[8]  = z4;
        _vertex_prism[9]  = x6;
        _vertex_prism[10] = y6;
        _vertex_prism[11] = z6;
        _vertex_prism[12] = x7;
        _vertex_prism[13] = y7;
        _vertex_prism[14] = z7;
        _vertex_prism[15] = x8;
        _vertex_prism[16] = y8;
        _vertex_prism[17] = z8;


        _uvwPn_sub_prism[0]  = uvw_nodes[3*idx2];
        _uvwPn_sub_prism[1]  = uvw_nodes[3*idx2+1];
        _uvwPn_sub_prism[2]  = uvw_nodes[3*idx2+2];
        _uvwPn_sub_prism[3]  = uvw_nodes[3*idx3];
        _uvwPn_sub_prism[4]  = uvw_nodes[3*idx3+1];
        _uvwPn_sub_prism[5]  = uvw_nodes[3*idx3+2];
        _uvwPn_sub_prism[6]  = uvw_nodes[3*idx4];
        _uvwPn_sub_prism[7]  = uvw_nodes[3*idx4+1];
        _uvwPn_sub_prism[8]  = uvw_nodes[3*idx4+2];
        _uvwPn_sub_prism[9]  = uvw_nodes[3*idx6];
        _uvwPn_sub_prism[10] = uvw_nodes[3*idx6+1];
        _uvwPn_sub_prism[11] = uvw_nodes[3*idx6+2];
        _uvwPn_sub_prism[12] = uvw_nodes[3*idx7];
        _uvwPn_sub_prism[13] = uvw_nodes[3*idx7+1];
        _uvwPn_sub_prism[14] = uvw_nodes[3*idx7+2];
        _uvwPn_sub_prism[15] = uvw_nodes[3*idx8];
        _uvwPn_sub_prism[16] = uvw_nodes[3*idx8+1];
        _uvwPn_sub_prism[17] = uvw_nodes[3*idx8+2];

        _prism_to_tetra (_vertex_prism,
                         _vertex_tetra,
                         _uvwPn_sub_prism,
                         _uvw_vertex_tetra);

        for (int itetra = 0; itetra < 3; itetra++) {

          double *__vertex_coords = _vertex_tetra + 12*itetra;
          double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

          double _closest_point_P1[3];
          double _weights_closest_point_P1[4];
          double _dist2;

          PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                                 __vertex_coords,
                                                                 _closest_point_P1,
                                                                 &_dist2,
                                                                 _weights_closest_point_P1);

          if (status != PDM_TETRAHEDRON_DEGENERATED) {

            double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
            double _uvw_closest_point_Pn[3] = {0., 0., 0.};
            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 4; j2++) {
                _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_coords,
                            _uvwPn_sub_tetra,
                            _closest_point_P1,
                            _uvw_closest_point_P1,
                            _uvw_closest_point_Pn,
                            _dist2,
                            child++);

          }
        }
        k1++;
      } // End of i-loop


      int idx1 = iend - 1;
      int idx2 = iend - 1 + 1;
      int idx3 = iend + 1 + k1;
      int idx4 = n_nodes_basis + iend - 1;
      int idx5 = n_nodes_basis + iend - 1 + 1;
      int idx6 = n_nodes_basis + iend + 1 + k1;

      double x1 = nodes_coords[3*idx1];
      double y1 = nodes_coords[3*idx1 + 1];
      double z1 = nodes_coords[3*idx1 + 2];

      double x2 = nodes_coords[3*idx2];
      double y2 = nodes_coords[3*idx2 + 1];
      double z2 = nodes_coords[3*idx2 + 2];

      double x3 = nodes_coords[3*idx3];
      double y3 = nodes_coords[3*idx3 + 1];
      double z3 = nodes_coords[3*idx3 + 2];

      double x4 = nodes_coords[3*idx4];
      double y4 = nodes_coords[3*idx4 + 1];
      double z4 = nodes_coords[3*idx4 + 2];

      double x5 = nodes_coords[3*idx5];
      double y5 = nodes_coords[3*idx5 + 1];
      double z5 = nodes_coords[3*idx5 + 2];

      double x6 = nodes_coords[3*idx6];
      double y6 = nodes_coords[3*idx6 + 1];
      double z6 = nodes_coords[3*idx6 + 2];

      double _vertex_prism[18] = {x1, y1, z1,
                                  x2, y2, z2,
                                  x3, y3, z3,
                                  x4, y4, z4,
                                  x5, y5, z5,
                                  x6, y6, z6};

      double _uvwPn_sub_prism[18];
      _uvwPn_sub_prism[0]  = uvw_nodes[3*idx1];
      _uvwPn_sub_prism[1]  = uvw_nodes[3*idx1+1];
      _uvwPn_sub_prism[2]  = uvw_nodes[3*idx1+2];
      _uvwPn_sub_prism[3]  = uvw_nodes[3*idx2];
      _uvwPn_sub_prism[4]  = uvw_nodes[3*idx2+1];
      _uvwPn_sub_prism[5]  = uvw_nodes[3*idx2+2];
      _uvwPn_sub_prism[6]  = uvw_nodes[3*idx3];
      _uvwPn_sub_prism[7]  = uvw_nodes[3*idx3+1];
      _uvwPn_sub_prism[8]  = uvw_nodes[3*idx3+2];
      _uvwPn_sub_prism[9]  = uvw_nodes[3*idx4];
      _uvwPn_sub_prism[10] = uvw_nodes[3*idx4+1];
      _uvwPn_sub_prism[11] = uvw_nodes[3*idx4+2];
      _uvwPn_sub_prism[12] = uvw_nodes[3*idx5];
      _uvwPn_sub_prism[13] = uvw_nodes[3*idx5+1];
      _uvwPn_sub_prism[14] = uvw_nodes[3*idx5+2];
      _uvwPn_sub_prism[15] = uvw_nodes[3*idx6];
      _uvwPn_sub_prism[16] = uvw_nodes[3*idx6+1];
      _uvwPn_sub_prism[17] = uvw_nodes[3*idx6+2];

      _prism_to_tetra (_vertex_prism,
                       _vertex_tetra,
                       _uvwPn_sub_prism,
                       _uvw_vertex_tetra);

      for (int itetra = 0; itetra < 3; itetra++) {

        double *__vertex_coords = _vertex_tetra + 12*itetra;
        double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

        double _closest_point_P1[3];
        double _weights_closest_point_P1[4];
        double _dist2;

        PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                               __vertex_coords,
                                                               _closest_point_P1,
                                                               &_dist2,
                                                               _weights_closest_point_P1);

        if (status != PDM_TETRAHEDRON_DEGENERATED) {

          double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
          double _uvw_closest_point_Pn[3] = {0., 0., 0.};
          for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 4; j2++) {
              _uvw_closest_point_Pn[j1] += _weights_closest_point_P1[j2] * _uvwPn_sub_tetra[3*j2 + j1];
            }
          }

          _heap_v_insert (heap,
                          __vertex_coords,
                          _uvwPn_sub_tetra,
                          _closest_point_P1,
                          _uvw_closest_point_P1,
                          _uvw_closest_point_Pn,
                          _dist2,
                          child++);

        }


      }

      ibeg = iend + 1;
      iend += order - j;
    } // End of j-loop

    ibeg ++;
    iend = ibeg + order;
  } // End of k-loop


  free (uvw_nodes);
}



/*----------------------------------------------------------------------------
 *
 *  Extract 5 tetrahedron from an hexahedron
 *
 *----------------------------------------------------------------------------*/

static void _hexa_to_tetra
(
 double *vertex_hexa,
 double *vertex_tetra,
 double *uvwPn_hexa,
 double *uvw_tetra
 )
{

  // 1st tetra
  vertex_tetra[0]  = vertex_hexa[0];  uvw_tetra[0]  = uvwPn_hexa[0];
  vertex_tetra[1]  = vertex_hexa[1];  uvw_tetra[1]  = uvwPn_hexa[1];
  vertex_tetra[2]  = vertex_hexa[2];  uvw_tetra[2]  = uvwPn_hexa[2];

  vertex_tetra[3]  = vertex_hexa[3];  uvw_tetra[3]  = uvwPn_hexa[3];
  vertex_tetra[4]  = vertex_hexa[4];  uvw_tetra[4]  = uvwPn_hexa[4];
  vertex_tetra[5]  = vertex_hexa[5];  uvw_tetra[5]  = uvwPn_hexa[5];

  vertex_tetra[6]  = vertex_hexa[6];  uvw_tetra[6]  = uvwPn_hexa[6];
  vertex_tetra[7]  = vertex_hexa[7];  uvw_tetra[7]  = uvwPn_hexa[7];
  vertex_tetra[8]  = vertex_hexa[8];  uvw_tetra[8]  = uvwPn_hexa[8];

  vertex_tetra[9]  = vertex_hexa[12]; uvw_tetra[9]  = uvwPn_hexa[12];
  vertex_tetra[10] = vertex_hexa[13]; uvw_tetra[10] = uvwPn_hexa[13];
  vertex_tetra[11] = vertex_hexa[14]; uvw_tetra[11] = uvwPn_hexa[14];

  // 2nd tetra
  vertex_tetra[12] = vertex_hexa[3];  uvw_tetra[12] = uvwPn_hexa[3];
  vertex_tetra[13] = vertex_hexa[4];  uvw_tetra[13] = uvwPn_hexa[4];
  vertex_tetra[14] = vertex_hexa[5];  uvw_tetra[14] = uvwPn_hexa[5];

  vertex_tetra[15] = vertex_hexa[12]; uvw_tetra[15] = uvwPn_hexa[12];
  vertex_tetra[16] = vertex_hexa[13]; uvw_tetra[16] = uvwPn_hexa[13];
  vertex_tetra[17] = vertex_hexa[14]; uvw_tetra[17] = uvwPn_hexa[14];

  vertex_tetra[18] = vertex_hexa[15]; uvw_tetra[18] = uvwPn_hexa[15];
  vertex_tetra[19] = vertex_hexa[16]; uvw_tetra[19] = uvwPn_hexa[16];
  vertex_tetra[20] = vertex_hexa[17]; uvw_tetra[20] = uvwPn_hexa[17];

  vertex_tetra[21] = vertex_hexa[21]; uvw_tetra[21] = uvwPn_hexa[21];
  vertex_tetra[22] = vertex_hexa[22]; uvw_tetra[22] = uvwPn_hexa[22];
  vertex_tetra[23] = vertex_hexa[23]; uvw_tetra[23] = uvwPn_hexa[23];

  // 3rd tetra
  vertex_tetra[24] = vertex_hexa[3];  uvw_tetra[24] = uvwPn_hexa[3];
  vertex_tetra[25] = vertex_hexa[4];  uvw_tetra[25] = uvwPn_hexa[4];
  vertex_tetra[26] = vertex_hexa[5];  uvw_tetra[26] = uvwPn_hexa[5];

  vertex_tetra[27] = vertex_hexa[6];  uvw_tetra[27] = uvwPn_hexa[6];
  vertex_tetra[28] = vertex_hexa[7];  uvw_tetra[28] = uvwPn_hexa[7];
  vertex_tetra[29] = vertex_hexa[8];  uvw_tetra[29] = uvwPn_hexa[8];

  vertex_tetra[30] = vertex_hexa[9];  uvw_tetra[30] = uvwPn_hexa[9];
  vertex_tetra[31] = vertex_hexa[10]; uvw_tetra[31] = uvwPn_hexa[10];
  vertex_tetra[32] = vertex_hexa[11]; uvw_tetra[32] = uvwPn_hexa[11];

  vertex_tetra[33] = vertex_hexa[21]; uvw_tetra[33] = uvwPn_hexa[21];
  vertex_tetra[34] = vertex_hexa[22]; uvw_tetra[34] = uvwPn_hexa[22];
  vertex_tetra[35] = vertex_hexa[23]; uvw_tetra[35] = uvwPn_hexa[23];

  // 4th tetra
  vertex_tetra[36] = vertex_hexa[6];  uvw_tetra[36] = uvwPn_hexa[6];
  vertex_tetra[37] = vertex_hexa[7];  uvw_tetra[37] = uvwPn_hexa[7];
  vertex_tetra[38] = vertex_hexa[8];  uvw_tetra[38] = uvwPn_hexa[8];

  vertex_tetra[39] = vertex_hexa[12]; uvw_tetra[39] = uvwPn_hexa[12];
  vertex_tetra[40] = vertex_hexa[13]; uvw_tetra[40] = uvwPn_hexa[13];
  vertex_tetra[41] = vertex_hexa[14]; uvw_tetra[41] = uvwPn_hexa[14];

  vertex_tetra[42] = vertex_hexa[18]; uvw_tetra[42] = uvwPn_hexa[18];
  vertex_tetra[43] = vertex_hexa[19]; uvw_tetra[43] = uvwPn_hexa[19];
  vertex_tetra[44] = vertex_hexa[20]; uvw_tetra[44] = uvwPn_hexa[20];

  vertex_tetra[45] = vertex_hexa[21]; uvw_tetra[45] = uvwPn_hexa[21];
  vertex_tetra[46] = vertex_hexa[22]; uvw_tetra[46] = uvwPn_hexa[22];
  vertex_tetra[47] = vertex_hexa[23]; uvw_tetra[47] = uvwPn_hexa[23];

  // 5th tetra
  vertex_tetra[48] = vertex_hexa[3];  uvw_tetra[48] = uvwPn_hexa[3];
  vertex_tetra[49] = vertex_hexa[4];  uvw_tetra[49] = uvwPn_hexa[4];
  vertex_tetra[50] = vertex_hexa[5];  uvw_tetra[50] = uvwPn_hexa[5];

  vertex_tetra[51] = vertex_hexa[6];  uvw_tetra[51] = uvwPn_hexa[6];
  vertex_tetra[52] = vertex_hexa[7];  uvw_tetra[52] = uvwPn_hexa[7];
  vertex_tetra[53] = vertex_hexa[8];  uvw_tetra[53] = uvwPn_hexa[8];

  vertex_tetra[54] = vertex_hexa[12]; uvw_tetra[54] = uvwPn_hexa[12];
  vertex_tetra[55] = vertex_hexa[13]; uvw_tetra[55] = uvwPn_hexa[13];
  vertex_tetra[56] = vertex_hexa[14]; uvw_tetra[56] = uvwPn_hexa[14];

  vertex_tetra[57] = vertex_hexa[21]; uvw_tetra[57] = uvwPn_hexa[21];
  vertex_tetra[58] = vertex_hexa[22]; uvw_tetra[58] = uvwPn_hexa[22];
  vertex_tetra[59] = vertex_hexa[23]; uvw_tetra[59] = uvwPn_hexa[23];
}



/*----------------------------------------------------------------------------
 *
 * Add sub-tetrahedron of a pn-hexahedron in the heap
 *
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_hexa_sub_tetra
(
 _heap_v_t    *heap,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords
 )
{
  int step = order + 1;
  int n_nodes_basis = step * step;

  double *uvw_nodes = malloc (sizeof(double) * 3 * n_nodes);
  double _vertex_tetra[60];
  double _uvw_vertex_tetra[60];

  _uvw_ho_hexa_nodes (order,
                      0., 1.,
                      0., 1.,
                      0., 1.,
                      uvw_nodes);

  int child = 0;
  for (int k = 0; k < order; k++) {
    int k1 = k+1;
    for (int j = 0; j < order; j++) {
      int j1 = j+1;
      for (int i = 0; i < order; i++) {
        int i1 = i+1;

        int idx1 = k * n_nodes_basis + j * step + i;
        int idx2 = k * n_nodes_basis + j * step + i1;
        int idx3 = k * n_nodes_basis + j1 * step + i;
        int idx4 = k * n_nodes_basis + j1 * step + i1;
        int idx5 = k1 * n_nodes_basis + j * step + i;
        int idx6 = k1 * n_nodes_basis + j * step + i1;
        int idx7 = k1 * n_nodes_basis + j1 * step + i;
        int idx8 = k1 * n_nodes_basis + j1 * step + i1;

        double x1 = nodes_coords[3*idx1];
        double y1 = nodes_coords[3*idx1 + 1];
        double z1 = nodes_coords[3*idx1 + 2];

        double x2 = nodes_coords[3*idx2];
        double y2 = nodes_coords[3*idx2 + 1];
        double z2 = nodes_coords[3*idx2 + 2];

        double x3 = nodes_coords[3*idx3];
        double y3 = nodes_coords[3*idx3 + 1];
        double z3 = nodes_coords[3*idx3 + 2];

        double x4 = nodes_coords[3*idx4];
        double y4 = nodes_coords[3*idx4 + 1];
        double z4 = nodes_coords[3*idx4 + 2];

        double x5 = nodes_coords[3*idx5];
        double y5 = nodes_coords[3*idx5 + 1];
        double z5 = nodes_coords[3*idx5 + 2];

        double x6 = nodes_coords[3*idx6];
        double y6 = nodes_coords[3*idx6 + 1];
        double z6 = nodes_coords[3*idx6 + 2];

        double x7 = nodes_coords[3*idx7];
        double y7 = nodes_coords[3*idx7 + 1];
        double z7 = nodes_coords[3*idx7 + 2];

        double x8 = nodes_coords[3*idx8];
        double y8 = nodes_coords[3*idx8 + 1];
        double z8 = nodes_coords[3*idx8 + 2];

        double _vertex_hexa[24] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3,
                                   x4, y4, z4,
                                   x5, y5, z5,
                                   x6, y6, z6,
                                   x7, y7, z7,
                                   x8, y8, z8};


        double _uvwPn_sub_hexa[24];
        _uvwPn_sub_hexa[0]  = uvw_nodes[3*idx1];
        _uvwPn_sub_hexa[1]  = uvw_nodes[3*idx1+1];
        _uvwPn_sub_hexa[2]  = uvw_nodes[3*idx1+2];
        _uvwPn_sub_hexa[3]  = uvw_nodes[3*idx2];
        _uvwPn_sub_hexa[4]  = uvw_nodes[3*idx2+1];
        _uvwPn_sub_hexa[5]  = uvw_nodes[3*idx2+2];
        _uvwPn_sub_hexa[6]  = uvw_nodes[3*idx3];
        _uvwPn_sub_hexa[7]  = uvw_nodes[3*idx3+1];
        _uvwPn_sub_hexa[8]  = uvw_nodes[3*idx3+2];
        _uvwPn_sub_hexa[9]  = uvw_nodes[3*idx4];
        _uvwPn_sub_hexa[10] = uvw_nodes[3*idx4+1];
        _uvwPn_sub_hexa[11] = uvw_nodes[3*idx4+2];
        _uvwPn_sub_hexa[12] = uvw_nodes[3*idx5];
        _uvwPn_sub_hexa[13] = uvw_nodes[3*idx5+1];
        _uvwPn_sub_hexa[14] = uvw_nodes[3*idx5+2];
        _uvwPn_sub_hexa[15] = uvw_nodes[3*idx6];
        _uvwPn_sub_hexa[16] = uvw_nodes[3*idx6+1];
        _uvwPn_sub_hexa[17] = uvw_nodes[3*idx6+2];
        _uvwPn_sub_hexa[18] = uvw_nodes[3*idx7];
        _uvwPn_sub_hexa[19] = uvw_nodes[3*idx7+1];
        _uvwPn_sub_hexa[20] = uvw_nodes[3*idx7+2];
        _uvwPn_sub_hexa[21] = uvw_nodes[3*idx8];
        _uvwPn_sub_hexa[22] = uvw_nodes[3*idx8+1];
        _uvwPn_sub_hexa[23] = uvw_nodes[3*idx8+2];


        _hexa_to_tetra (_vertex_hexa,
                        _vertex_tetra,
                        _uvwPn_sub_hexa,
                        _uvw_vertex_tetra);

        for (int itetra = 0; itetra < 5; itetra++) {

          double *__vertex_coords = _vertex_tetra + 12*itetra;
          double *_uvwPn_sub_tetra = _uvw_vertex_tetra + 12*itetra;

          double _closest_point_P1[3];
          double _weights_closest_point_P1[4];
          double _dist2;

          PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                                 __vertex_coords,
                                                                 _closest_point_P1,
                                                                 &_dist2,
                                                                 _weights_closest_point_P1);

          if (status != PDM_TETRAHEDRON_DEGENERATED) {

            double *_uvw_closest_point_P1 = _weights_closest_point_P1 + 1;
            double _uvw_closest_point_Pn[3] = {0., 0., 0.};
            for (int n1 = 0; n1 < 3; n1++) {
              for (int n2 = 0; n2 < 4; n2++) {
                _uvw_closest_point_Pn[n1] += _weights_closest_point_P1[n2] * _uvwPn_sub_tetra[3*n2 + n1];
              }
            }

            _heap_v_insert (heap,
                            __vertex_coords,
                            _uvwPn_sub_tetra,
                            _closest_point_P1,
                            _uvw_closest_point_P1,
                            _uvw_closest_point_Pn,
                            _dist2,
                            child++);
          }
        }
      }
    }
  }

  free (uvw_nodes);
}


/*-----------------------------------------------------------------------------
 *
 * Add children to a heap
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   vtx_tria_current <-- current triangle
 *   uvPn_tria_current<-- uv of current triangle vertices in the on element
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static void
_insert_subtetra
(
 _heap_v_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      vtx_tetra_current[],
 double                      uvwPn_tetra_current[]
 )
{

  double _vtx_tetra_children[30];
  double _uvwPn_tetra_children[30];

  const int idx_sub_tetra[32] = {0, 4, 6, 7,
                                 4, 9, 7, 6,
                                 4, 5, 9, 6,
                                 4, 9, 7, 8,
                                 4, 5, 9, 8,
                                 4, 1, 5, 8,
                                 6, 5, 2, 9,
                                 7, 8, 9, 3};



  // Compute middle vertices
  for (int i = 0; i < 12; i++) {
    _vtx_tetra_children[i] = vtx_tetra_current[i];
  }

  for (int i = 0; i < 12; i++) {
    _uvwPn_tetra_children[i] = uvwPn_tetra_current[i];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {

      _uvwPn_tetra_children[12+3*i+j] = 0.5 *
        (uvwPn_tetra_current[3*i+j] + uvwPn_tetra_current[3*((i+1)%3)+j]);

      _uvwPn_tetra_children[21+3*i+j] = 0.5 *
        (uvwPn_tetra_current[3*i+j] + uvwPn_tetra_current[9+j]);

    }

    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _uvwPn_tetra_children + 12 + 3*i,
                  weightsPn);


    for (int j = 0; j < 3; j++) {
      _vtx_tetra_children[12+3*i+j] = 0;
    }
    for (int k = 0; k < n_nodes; k++) {
      const double *_node_coords = nodes_coords + 3 * k;
      for (int j = 0; j < 3; j++) {
        _vtx_tetra_children[12+3*i+j] += weightsPn[k] * _node_coords[j];
      }
    }








    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _uvwPn_tetra_children + 21 + 3*i,
                  weightsPn);


    for (int j = 0; j < 3; j++) {
      _vtx_tetra_children[21+3*i+j] = 0;
    }
    for (int k = 0; k < n_nodes; k++) {
      const double *_node_coords = nodes_coords + 3 * k;
      for (int j = 0; j < 3; j++) {
        _vtx_tetra_children[21+3*i+j] += weightsPn[k] * _node_coords[j];
      }
    }
  }


  int child = 0;
  for (int i = 0; i < 8; i++) {

    double _vtx_tetra_child[12];
    double _uvwPn_tetra_child[12];

    for (int j = 0; j < 4; j++) {
      int _j = idx_sub_tetra[4 * i + j];
      for (int k = 0; k < 3; k++) {
        _vtx_tetra_child[3*j + k] = _vtx_tetra_children[3*_j + k];
      }
      for (int k = 0; k < 3; k++) {
        _uvwPn_tetra_child[3*j + k] = _uvwPn_tetra_children[3*_j + k];
      }
    }




    double _closest_pt_child[3];
    double _dist2_child = 0;
    double _closest_pt_weights_child[4];

    PDM_tetrahedron_status_t status = PDM_tetrahedron_evaluate_position (point_coords,
                                                           _vtx_tetra_child,
                                                           _closest_pt_child,
                                                           &_dist2_child,
                                                           _closest_pt_weights_child);



    if (status == PDM_TETRAHEDRON_DEGENERATED) {
      continue;
    }

    double *_closest_pt_uvwP1_child = _closest_pt_weights_child + 1;
    double _closest_pt_uvwPn_child[3] = {0., 0., 0.};
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 4; k++) {
        _closest_pt_uvwPn_child[j] +=
          _closest_pt_weights_child[k] * _uvwPn_tetra_child[3*k + j];
      }
    }

    _heap_v_insert (heap,
                    _vtx_tetra_child,
                    _uvwPn_tetra_child,
                    _closest_pt_child,
                    _closest_pt_uvwP1_child,
                    _closest_pt_uvwPn_child,
                    _dist2_child,
                    child++);

  }

}


/*-----------------------------------------------------------------------------
 *
 * compute distance from closest tetrahedron subdivision
 *
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   uncertain_result --> 1 if the result is uncertain
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_closest_tetra_subdivision
(
 _heap_v_t                  *heap,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj,
 int                        *uncertain_result
 )
{
  const int dbg_enabled = 0;
  if (dbg_enabled) {
    log_trace(">> _compute_dist2_from_closest_tetra_subdivision\n");
  }
  *uncertain_result = 0;
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;
  double dist2_pre = HUGE_VAL;
  int distance_extension = 0;

  while (1) {

    double _vtx_tetra_current[12];
    double _uvwPn_tetra_current[12];

    double _closest_pt_current[3];
    double _closest_pt_uvwP1_current[3];
    double _closest_pt_uvwPn_current[3];
    double _dist2_current;

    // Get closest tetrahedron stored in the heap

    int _child;
    int is_empty = _heap_v_top_get (heap,
                                    _vtx_tetra_current,
                                    _uvwPn_tetra_current,
                                    _closest_pt_current,
                                    _closest_pt_uvwP1_current,
                                    _closest_pt_uvwPn_current,
                                    &_dist2_current,
                                    &_child);

    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty %s\n");
      abort();
    }

    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "closest_tetra_it%d.vtk", *n_it);
      int connec[4] = {1, 2, 3, 4};
      PDM_vtk_write_std_elements(filename,
                                 4,
                                 _vtx_tetra_current,
                                 NULL,
                                 PDM_MESH_NODAL_TETRA4,
                                 1,
                                 connec,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }

    if ((distance_extension == 0) && (_dist2_current > dist2_pre)) {
      distance_extension = 1;
    }

    else if (distance_extension == 1) {
      if (_dist2_current <= dist2_min_min) {
        distance_extension = 0;
      }
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);
    dist2_pre = _dist2_current;

    /* Compute projected from current P1 tetrahedron */
    double weightsP1[4];
    PDM_ho_basis (PDM_MESH_NODAL_TETRAHO,
                  1,
                  4,
                  1,
                  _closest_pt_uvwP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 4; k++) {
        _projected_coords_from_p1[j] += weightsP1[k] * _vtx_tetra_current[3*k+j];
      }
    }

    /* Compute projected from current Pn element */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uvwPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }


    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */
    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {
      if (dbg_enabled) {
        log_trace("finished at n_it = %d, err_proj = %e/%e\n", *n_it, sqrt(*err_proj), err_max);
      }

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uvwPn_current[0];
      uvw[1] = _closest_pt_uvwPn_current[1];
      uvw[2] = _closest_pt_uvwPn_current[2];

      break;
    }

    /* Insert sub-tetrahedra in the heap */
    _insert_subtetra (heap,
                      order,
                      type,
                      n_nodes,
                      nodes_coords,
                      point_coords,
                      weightsPn,
                      _vtx_tetra_current,
                      _uvwPn_tetra_current);


  }


  if (distance_extension) {
    *uncertain_result = 1;
  }

  return dist2;
}



/*-----------------------------------------------------------------------------
 *
 * compute distance from uniform tetrahedron subdivision
 *
 * parameters:
 *   heap1            <-> heap
 *   heap2            <-> work heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   _basis_generic   <-- generic basis
 *
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_uniform_tetra_subdivision
(
 _heap_v_t                  *heap1,
 _heap_v_t                  *heap2,
 const int                   order,
 const PDM_Mesh_nodal_elt_t  type,
 const int                   n_nodes,
 const int                   n_it_max,
 const double                err_max,
 const double                nodes_coords[],
 const double                point_coords[],
 double                      weightsPn[],
 double                      projected_coords[],
 double                      uvw[],
 int                        *n_it,
 double                     *err_proj
 )
{
  const int dbg_enabled = 0;
  if (dbg_enabled) {
    log_trace(">> _compute_dist2_from_uniform_tetra_subdivision\n");
  }
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;

  _heap_v_t *heap      = heap1;
  _heap_v_t *next_heap = heap2;

  while (1) {

    double _vtx_tetra_current[12];
    double _uvwPn_tetra_current[12];

    double _closest_pt_current[3];
    double _closest_pt_uvwP1_current[3];
    double _closest_pt_uvwPn_current[3];
    double _dist2_current;

    /* Get closest triangle stored in the heap */
    int _child;
    int is_empty = _heap_v_top_get (heap,
                                    _vtx_tetra_current,
                                    _uvwPn_tetra_current,
                                    _closest_pt_current,
                                    _closest_pt_uvwP1_current,
                                    _closest_pt_uvwPn_current,
                                    &_dist2_current,
                                    &_child);

    if (is_empty) {
      PDM_error(__FILE__, __LINE__, 0,
                "Heap is empty\n");
      abort();
    }

    dist2_min_min = PDM_MIN (dist2_min_min, _dist2_current);

    /* Compute projected from current P1 tetrahedron */
    double weightsP1[4];
    PDM_ho_basis (PDM_MESH_NODAL_TETRAHO,
                  1,
                  4,
                  1,
                  _closest_pt_uvwP1_current,
                  weightsP1);

    double _projected_coords_from_p1[3] = {0., 0., 0.};
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_tetra_current[3*j+k];
      }
    }

    /* Compute projected from current Pn element */
    PDM_ho_basis (type,
                  order,
                  n_nodes,
                  1,
                  _closest_pt_uvwPn_current,
                  weightsPn);

    double _projected_coords_from_pn[3] = {0., 0., 0.};
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k];
      }
    }

    /* Compute distance between two projected */
    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */

    if (sqrt(*err_proj) <= err_max || (*n_it)++ >= n_it_max) {
      if (dbg_enabled) {
        log_trace("finished at n_it = %d, err_proj = %e/%e\n", *n_it, sqrt(*err_proj), err_max);
      }

      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp;
      }

      uvw[0] = _closest_pt_uvwPn_current[0];
      uvw[1] = _closest_pt_uvwPn_current[1];
      uvw[2] = _closest_pt_uvwPn_current[2];

      break;
    }

    /* Insert sub-tetrahedra in the next heap */

    _heap_v_init (next_heap);

    _insert_subtetra (next_heap,
                      order,
                      type,
                      n_nodes,
                      nodes_coords,
                      point_coords,
                      weightsPn,
                      _vtx_tetra_current,
                      _uvwPn_tetra_current);

    double _vtx_tetra_current2[12];
    double _uvwPn_tetra_current2[12];

    double _closest_pt_current2[3];
    double _closest_pt_uvwP1_current2[3];
    double _dist2_current2;
    int _child_current2;

    while ( !_heap_v_top_get (heap,
                              _vtx_tetra_current2,
                              _uvwPn_tetra_current2,
                              _closest_pt_current2,
                              _closest_pt_uvwP1_current2,
                              _closest_pt_uvwPn_current,
                              &_dist2_current2,
                              &_child_current2)) {


      _insert_subtetra (next_heap,
                        order,
                        type,
                        n_nodes,
                        nodes_coords,
                        point_coords,
                        weightsPn,
                        _vtx_tetra_current2,
                        _uvwPn_tetra_current2);

    }

    _heap_v_t *heap_tmp = heap;
    heap = next_heap;
    next_heap = heap_tmp;

  }

  return dist2;
}


/*-----------------------------------------------------------------------------
 *
 * Default point location on a high order tetrahedron
 *
 * parameters:
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   uvw              --> parametric coordinates in the element
 *
 * return:
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double
_default_location_generic_3d
(
 const PDM_Mesh_nodal_elt_t    type,
 const int                     order,
 const double                  char_size,
 const int                     n_nodes,
 const double                 *nodes_coords,
 const double                 *point_coords,
 double                       *projected_coords,
 double                       *uvw,
 _heap_fill_init_sub_tetra_te  fill_init_fct
 )
{
  const int n_it_max = 10000;
  double err_max = PDM_MAX (char_size * 1e-6, 1e-15);

  double dist2 = HUGE_VAL;

  _heap_v_t heap;
  _heap_v_t heap2;

  double *weightsPn = malloc (sizeof(double) * n_nodes);

  // Initialize heap

  _heap_v_init (&heap);

  // Build initial sub-tetrahedron and store them in the heap

  (fill_init_fct) (&heap,
                   order,
                   n_nodes,
                   nodes_coords,
                   point_coords);

  //  While error > error_max    - Get closest triangle in the heap    - Cut it in sub-triangles   - Store them in the heap



  const int method = 0;

  int n_it;
  double err_proj = HUGE_VAL;
  int uncertain_result = 0;;

  if (method == 0) {
    dist2 = _compute_dist2_from_closest_tetra_subdivision (&heap,
                                                           order,
                                                           type,
                                                           n_nodes,
                                                           n_it_max,
                                                           err_max,
                                                           nodes_coords,
                                                           point_coords,
                                                           weightsPn,
                                                           projected_coords,
                                                           uvw,
                                                           &n_it,
                                                           &err_proj,
                                                           &uncertain_result);

    if (uncertain_result) {

      // Initialize heap

      _heap_v_init (&heap);
      _heap_v_init (&heap2);

      // Build initial sub-triangles and store them in the heap

      (fill_init_fct) (&heap,
                       order,
                       n_nodes,
                       nodes_coords,
                       point_coords);


      dist2 = _compute_dist2_from_uniform_tetra_subdivision (&heap,
                                                             &heap2,
                                                             order,
                                                             type,
                                                             n_nodes,
                                                             n_it_max,
                                                             err_max,
                                                             nodes_coords,
                                                             point_coords,
                                                             weightsPn,
                                                             projected_coords,
                                                             uvw,
                                                             &n_it,
                                                             &err_proj);

    }

  }

  else {

    _heap_v_init (&heap2);
    dist2 = _compute_dist2_from_uniform_tetra_subdivision (&heap,
                                                           &heap2,
                                                           order,
                                                           type,
                                                           n_nodes,
                                                           n_it_max,
                                                           err_max,
                                                           nodes_coords,
                                                           point_coords,
                                                           weightsPn,
                                                           projected_coords,
                                                           uvw,
                                                           &n_it,
                                                           &err_proj);

  }

  free (weightsPn);




  return dist2;

}

static double
_reference_length
(
 const double *coords
 )
{
  double a =
    sqrt ((coords[3*1    ] - coords[3*0    ]) * (coords[3*1    ] - coords[3*0    ]) +
          (coords[3*1 + 1] - coords[3*0 + 1]) * (coords[3*1 + 1] - coords[3*0 + 1]) +
          (coords[3*1 + 2] - coords[3*0 + 2]) * (coords[3*1 + 2] - coords[3*0 + 2]));

  double b =
    sqrt ((coords[3*2    ] - coords[3*1    ]) * (coords[3*2    ] - coords[3*1    ]) +
          (coords[3*2 + 1] - coords[3*1 + 1]) * (coords[3*2 + 1] - coords[3*1 + 1]) +
          (coords[3*2 + 2] - coords[3*1 + 2]) * (coords[3*2 + 2] - coords[3*1 + 2]));

  double c =
    sqrt ((coords[3*0    ] - coords[3*2    ]) * (coords[3*0    ] - coords[3*2    ]) +
          (coords[3*0 + 1] - coords[3*2 + 1]) * (coords[3*0 + 1] - coords[3*2 + 1]) +
          (coords[3*0 + 2] - coords[3*2 + 2]) * (coords[3*0 + 2] - coords[3*2 + 2]));

  double d =
    sqrt ((coords[3*3    ] - coords[3*0    ]) * (coords[3*3    ] - coords[3*0    ]) +
          (coords[3*3 + 1] - coords[3*0 + 1]) * (coords[3*3 + 1] - coords[3*0 + 1]) +
          (coords[3*3 + 2] - coords[3*0 + 2]) * (coords[3*3 + 2] - coords[3*0 + 2]));

  double e =
    sqrt ((coords[3*3    ] - coords[3*1    ]) * (coords[3*3    ] - coords[3*1    ]) +
          (coords[3*3 + 1] - coords[3*1 + 1]) * (coords[3*3 + 1] - coords[3*1 + 1]) +
          (coords[3*3 + 2] - coords[3*1 + 2]) * (coords[3*3 + 2] - coords[3*1 + 2]));

  double f =
    sqrt ((coords[3*3    ] - coords[3*2    ]) * (coords[3*3    ] - coords[3*2    ]) +
          (coords[3*3 + 1] - coords[3*2 + 1]) * (coords[3*3 + 1] - coords[3*2 + 1]) +
          (coords[3*3 + 2] - coords[3*2 + 2]) * (coords[3*3 + 2] - coords[3*2 + 2]));



  double max_ab = PDM_MAX (a,b);
  double max_cd = PDM_MAX (c,d);
  double max_ef = PDM_MAX (e,f);

  double max_abcd = PDM_MAX (max_ab,max_cd);
  return PDM_MAX (max_abcd,max_ef);

}


/*----------------------------------------------------------------------------
 *
 * Point location on a high order cell
 *
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates (size = 3 * n_nodes)
 *   point_coords     <-- point to locate coordinates (size = 3)
 *   projected_coords --> projected point coordinates (size = 3)
 *   uvw              --> parametric coordinates of the projected point on the element
 *
 * return:
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double
_default_location
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw
 )
{
  PDM_Mesh_nodal_elt_t _type = type;
  double dist2 = HUGE_VAL;


  switch (_type) {

  case PDM_MESH_NODAL_BARHO: {

    int v1 = 0;
    int v2 = order;

    double p1_coords[6] =
      {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
       nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2]};

    double l2 =
      (p1_coords[0] - p1_coords[3]) * (p1_coords[0] - p1_coords[3]) +
      (p1_coords[1] - p1_coords[4]) * (p1_coords[1] - p1_coords[4]) +
      (p1_coords[2] - p1_coords[5]) * (p1_coords[2] - p1_coords[5]);

    double len_ref = sqrt (l2);

    dist2 = _default_location_generic_1d (type,
                                          order,
                                          len_ref,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_sub_edge);
    break;

  }


  case PDM_MESH_NODAL_TRIAHO: {

    int v1 = 0;
    int v2 = order;
    int v3 = (order + 2) * (order + 1) / 2 - 1;

    double p1_coords[9] =
      {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
       nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
       nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2]};

    double char_size = _radius_inscribed_circle (p1_coords);

    dist2 = _default_location_generic_2d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_sub_tria);
    break;

  }


  case PDM_MESH_NODAL_QUADHO: {

    const int order1 = order + 1;

    const int quadrangle_vertices[4] = {1,
                                        order1,
                                        order1 * order + 1,
                                        order1 * order1};
    int triangle_vertices[6];
    int n_sub_tria = PDM_triangulate_quadrangle (3,
                                                 nodes_coords,
                                                 NULL,
                                                 quadrangle_vertices,
                                                 triangle_vertices);

    double char_size = -1.0;
    for (int i = 0; i < n_sub_tria; i++) {

      int v1 = triangle_vertices[3*i    ] - 1;
      int v2 = triangle_vertices[3*i + 1] - 1;
      int v3 = triangle_vertices[3*i + 2] - 1;

      double p1_coords[9] =
        {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
         nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
         nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2]};

      double _char_size = _radius_inscribed_circle (p1_coords);

      char_size = PDM_MAX (char_size, _char_size);

    }

    dist2 = _default_location_generic_2d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_qn_sub_tria);
    break;
  }


  case PDM_MESH_NODAL_TETRAHO: {

    int v1 = 0;
    int v2 = order;
    int v3 = (order + 2) * (order + 1) / 2 - 1;
    int v4 = (order+1)*(order+2)*(order+3)/6 - 1;

    double p1_coords[12] =
      {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
       nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
       nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2],
       nodes_coords[3*v4], nodes_coords[3*v4+1], nodes_coords[3*v4+2]};

    double char_size = _reference_length (p1_coords);

    dist2 = _default_location_generic_3d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_tetra_sub_tetra);

    break;

  }


  case PDM_MESH_NODAL_PYRAMIDHO: {

    const int pyra_vertices[5] = {1,
                                  order + 1,
                                  order * (order + 1) + 1,
                                  (order + 1) * (order + 1),
                                  (order+1)*(order+2)*(2*order+3)/6};

    int tetra_vertices[8];
    int n_sub_tetra = PDM_triangulate_pyramid (3,
                                               nodes_coords,
                                               NULL,
                                               pyra_vertices,
                                               tetra_vertices);

    double char_size = -1.0;
    for (int i = 0; i < n_sub_tetra; i++) {

      int v1 = tetra_vertices[4*i    ] - 1;
      int v2 = tetra_vertices[4*i + 1] - 1;
      int v3 = tetra_vertices[4*i + 2] - 1;
      int v4 = tetra_vertices[4*i + 3] - 1;

      double p1_coords[12] =
        {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
         nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
         nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2],
         nodes_coords[3*v4], nodes_coords[3*v4+1], nodes_coords[3*v4+2]};

      double _char_size = _reference_length (p1_coords);

      char_size = PDM_MAX (char_size, _char_size);

    }

    dist2 = _default_location_generic_3d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_pyra_sub_tetra);
    break;
  }


  case PDM_MESH_NODAL_PRISMHO: {

    const int prism_vertices[6] = {1,
                                   order + 1,
                                   (order + 2) * (order + 1) / 2,
                                   order * (order + 2) * (order + 1) / 2 + 1,
                                   order * (order + 2) * (order + 1) / 2 + order + 1,
                                   (order+1)*(order+1)*(order+2)/2};

    int tetra_vertices[12];
    int n_sub_tetra = PDM_triangulate_prism (3,
                                             nodes_coords,
                                             NULL,
                                             prism_vertices,
                                             tetra_vertices);

    double char_size = -1.0;
    for (int i = 0; i < n_sub_tetra; i++) {

      int v1 = tetra_vertices[4*i    ] - 1;
      int v2 = tetra_vertices[4*i + 1] - 1;
      int v3 = tetra_vertices[4*i + 2] - 1;
      int v4 = tetra_vertices[4*i + 3] - 1;

      double p1_coords[12] =
        {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
         nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
         nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2],
         nodes_coords[3*v4], nodes_coords[3*v4+1], nodes_coords[3*v4+2]};

      double _char_size = _reference_length (p1_coords);

      char_size = PDM_MAX (char_size, _char_size);

    }

    dist2 = _default_location_generic_3d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_prism_sub_tetra);
    break;
  }


  case PDM_MESH_NODAL_HEXAHO: {

    const int hexa_vertices[8] = {1,
                                  order + 1,
                                  order * (order + 1) + 1,
                                  (order + 1) * (order + 1),
                                  order * (order + 1) * (order + 1) + 1,
                                  order * ((order + 1) * (order + 1) + 1) + 1,
                                  order * (order + 1) * (order + 2) + 1,
                                  (order + 1) * (order + 1) * (order + 1)};

    int tetra_vertices[20];
    int n_sub_tetra = PDM_triangulate_hexahedron (3,
                                                  nodes_coords,
                                                  NULL,
                                                  hexa_vertices,
                                                  tetra_vertices);

    double char_size = -1.0;
    for (int i = 0; i < n_sub_tetra; i++) {

      int v1 = tetra_vertices[4*i    ] - 1;
      int v2 = tetra_vertices[4*i + 1] - 1;
      int v3 = tetra_vertices[4*i + 2] - 1;
      int v4 = tetra_vertices[4*i + 3] - 1;

      double p1_coords[12] =
        {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
         nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
         nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2],
         nodes_coords[3*v4], nodes_coords[3*v4+1], nodes_coords[3*v4+2]};

      double _char_size = _reference_length (p1_coords);


      char_size = PDM_MAX (char_size, _char_size);

    }

    dist2 = _default_location_generic_3d (type,
                                          order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_hexa_sub_tetra);
    break;
  }




  default:
    PDM_error(__FILE__, __LINE__, 0,
              "PDM_ho_location : %d Element type not implemented yet\n", _type);

  }

  return dist2;

}


/*----------------------------------------------------------------------------
 *
 * high order basis
 *
 * parameters:
 *   type            <-- element type
 *
 * return:
 *
 *----------------------------------------------------------------------------*/

static PDM_ho_location_user_elt_t **
_get_user_elt (PDM_Mesh_nodal_elt_t elt_type)
{

  switch(elt_type) {

  case PDM_MESH_NODAL_BARHO:
    return &_user_edge;
    break;

  case PDM_MESH_NODAL_TRIAHO:
    return &_user_tria;
    break;

  case PDM_MESH_NODAL_QUADHO:
    return &_user_quad;
    break;

  case PDM_MESH_NODAL_TETRAHO:
    return &_user_tetra;
    break;

  case PDM_MESH_NODAL_PYRAMIDHO:
    return &_user_pyra;
    break;

  case PDM_MESH_NODAL_PRISMHO:
    return &_user_prism;
    break;

  case PDM_MESH_NODAL_HEXAHO:
    return &_user_hexa;
    break;

  case PDM_MESH_NODAL_POINT:
  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_QUAD4:
  case PDM_MESH_NODAL_TETRA4:
  case PDM_MESH_NODAL_PYRAMID5:
  case PDM_MESH_NODAL_PRISM6:
  case PDM_MESH_NODAL_HEXA8:
    return NULL;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,
              "Unvailable element type %d\n", (int) elt_type);
  }

  return NULL;
}





static inline int
_argmax
(
 const int     n,
 const double *a
 )
{
  double amax = -HUGE_VAL;
  int    imax = -1;
  for (int i = 0; i < n; i++) {
    if (a[i] > amax) {
      amax = a[i];
      imax = i;
    }
  }

  return imax;
}

#define _tria_ij2idx(i, j, order)        ((i) + (j)*((order) + 1) - ((j)-1)*(j)/2)
#define _quad_ij2idx(i, j, order)        ((i) + (j)*((order) + 1))
#define _tetra_ijk2idx(i, j, k, order)   ((i) + (j)*((order) + 1 - (k)) - (j)*((j)-1)/2 + ((k)*((k)*((k) - 3*(order) - 6) + 3*(order)*((order) + 4) + 11)) / 6)
#define _pyramid_ijk2idx(i, j, k, order) ((i) + (j)*((order)+1-(k)) + ((k)*((k)*(2*(k) - 6*(order) - 9) + 6*(order)*((order) + 3) + 13)) / 6)
#define _prism_ijk2idx(i, j, k, order)   ((i) + (j)*((order)+1) - (j)*((j)-1)/2 + (k)*((order)+1)*((order)+2)/2)
#define _hexa_ijk2idx(i, j, k, order)    ((i) + ((order)+1)*((j) + ((order)+1)*(k)))

static int
_boundary_child_element_coord
(
 const PDM_Mesh_nodal_elt_t  parent_elt_type,
 const int                   order,
 const double               *parent_uvw,
 const double               *parent_node_coord,
       PDM_Mesh_nodal_elt_t *child_elt_type,
       int                  *child_n_node,
       double               *child_node_coord
 )
 {
  int i_bnd = -1;

  switch (parent_elt_type) {
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      *child_elt_type = PDM_MESH_NODAL_POINT;
      if (parent_uvw[0] < 0) {
        i_bnd = 0;
        // memcpy(child_node_coord, parent_node_coord, sizeof(double) * 3);
      }
      else if (parent_uvw[0] > 1) {
        i_bnd = 1;
        // memcpy(child_node_coord, parent_node_coord + 3*order, sizeof(double) * 3);
      }
      break;
    }

    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      *child_elt_type = PDM_MESH_NODAL_BARHO;
      if (parent_elt_type == PDM_MESH_NODAL_TRIAHO_BEZIER) {
        *child_elt_type = PDM_MESH_NODAL_BARHO_BEZIER;
      }

      double dist_bnd[3] = {
        -parent_uvw[1],
        parent_uvw[0] + parent_uvw[1] - 1,
        -parent_uvw[0]
      };
      int imax = _argmax(3, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        switch (imax) {
          case 0: {
            memcpy(child_node_coord, parent_node_coord, sizeof(double) * 3 * (order+1));
            break;
          }
          case 1: {
            for (int j = 0; j <= order; j++) {
              int idx = _tria_ij2idx(order-j, j, order);
              memcpy(child_node_coord + 3*j, parent_node_coord + 3*idx, sizeof(double) * 3);
            }
            break;
          }
          case 2: {
            for (int j = 0; j <= order; j++) {
              int idx = _tria_ij2idx(0, order-j, order);
              memcpy(child_node_coord + 3*j, parent_node_coord + 3*idx, sizeof(double) * 3);
            }
            break;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_QUADHO: {
      *child_elt_type = PDM_MESH_NODAL_BARHO;

      double dist_bnd[4] = {
        -parent_uvw[1],
        parent_uvw[0] - 1,
        parent_uvw[1] - 1,
        -parent_uvw[0]
      };
      int imax = _argmax(4, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        switch (imax) {
          case 0: {
            memcpy(child_node_coord, parent_node_coord, sizeof(double) * 3 * (order+1));
            break;
          }
          case 1: {
            for (int j = 0; j <= order; j++) {
              int idx = _quad_ij2idx(order, j, order);
              memcpy(child_node_coord + 3*j, parent_node_coord + 3*idx, sizeof(double) * 3);
            }
            break;
          }
          case 2: {
            for (int j = 0; j <= order; j++) {
              int idx = _quad_ij2idx(order-j, order, order);
              memcpy(child_node_coord + 3*j, parent_node_coord + 3*idx, sizeof(double) * 3);
            }
            break;
          }
          case 3: {
            for (int j = 0; j <= order; j++) {
              int idx = _quad_ij2idx(0, order-j, order);
              memcpy(child_node_coord + 3*j, parent_node_coord + 3*idx, sizeof(double) * 3);
            }
            break;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_TETRAHO: {

      *child_elt_type = PDM_MESH_NODAL_TRIAHO;

      double dist_bnd[4] = {
        -parent_uvw[1],
        parent_uvw[0] + parent_uvw[1] + parent_uvw[2] - 1,
        -parent_uvw[0],
        -parent_uvw[2]
      };
      int imax = _argmax(4, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        switch (imax) {
          case 0: {
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int i = 0; i <= order - k; i++) {
                int idx_read = _tetra_ijk2idx(i, 0, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 1: {
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order - k; j++) {
                int idx_read = _tetra_ijk2idx(order-j-k, j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 2: {
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order - k; j++) {
                int idx_read = _tetra_ijk2idx(0, order-k-j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 3: {
            memcpy(child_node_coord, parent_node_coord, sizeof(double) * 3 * (order+1)*(order+2)/2);
            break;
          }
        }
      }

      break;
    }

    case PDM_MESH_NODAL_PYRAMIDHO: {

      double dist_bnd[5] = {
        -parent_uvw[2],
        -parent_uvw[1],
        parent_uvw[0] + parent_uvw[2] - 1,
        parent_uvw[1] + parent_uvw[2] - 1,
        -parent_uvw[0],
      };
      int imax = _argmax(5, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        switch (imax) {
          case 0: {
            *child_elt_type = PDM_MESH_NODAL_QUADHO;
            memcpy(child_node_coord, parent_node_coord, sizeof(double) * (order+1)*(order+1) * 3);
            break;
          }
          case 1: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int i = 0; i <= order - k; i++) {
                int idx_read = _pyramid_ijk2idx(i, 0, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 2: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order - k; j++) {
                int idx_read = _pyramid_ijk2idx(order-k, j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 3: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int i = 0; i <= order - k; i++) {
                int idx_read = _pyramid_ijk2idx(order-k-i, order-k, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 4: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order - k; j++) {
                int idx_read = _pyramid_ijk2idx(0, order-k-j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
        }
      }

      break;
    }

    case PDM_MESH_NODAL_PRISMHO: {

      double dist_bnd[5] = {
        -parent_uvw[2],
        parent_uvw[2] - 1,
        -parent_uvw[1],
        parent_uvw[0] + parent_uvw[1] - 1,
        -parent_uvw[0],
      };
      int imax = _argmax(5, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        switch (imax) {
          case 0: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            memcpy(child_node_coord, parent_node_coord, sizeof(double) * (order+1)*(order+1) * 3);
            break;
          }
          case 1: {
            *child_elt_type = PDM_MESH_NODAL_TRIAHO;
            int idx_write = 0;
            for (int j = 0; j <= order; j++) {
              for (int i = 0; i <= order - j; i++) {
                int idx_read = _prism_ijk2idx(i, j, order, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 2: {
            *child_elt_type = PDM_MESH_NODAL_QUADHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int i = 0; i <= order; i++) {
                int idx_read = _prism_ijk2idx(i, 0, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 3: {
            *child_elt_type = PDM_MESH_NODAL_QUADHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order; j++) {
                int idx_read = _prism_ijk2idx(order-j, j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
          case 4: {
            *child_elt_type = PDM_MESH_NODAL_QUADHO;
            int idx_write = 0;
            for (int k = 0; k <= order; k++) {
              for (int j = 0; j <= order; j++) {
                int idx_read = _prism_ijk2idx(0, order-j, k, order);
                memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
                idx_write++;
              }
            }
            break;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_HEXAHO: {
      *child_elt_type = PDM_MESH_NODAL_QUADHO;

      double dist_bnd[6];
      for (int i = 0; i < 3; i++) {
        dist_bnd[2*i  ] = -parent_uvw[i];
        dist_bnd[2*i+1] = parent_uvw[i] - 1;
      }

      int imax = _argmax(6, dist_bnd);

      if (dist_bnd[imax] > 0) {
        i_bnd = imax;

        int ivar = imax/2;
        int ival = order * (imax%2);

        int ijk[3];
        int idx_write = 0;
        for (int jj = 0; jj <= order; jj++) {
          for (int ii = 0; ii <= order; ii++) {
            ijk[ivar      ] = ival;
            ijk[(ivar+1)%3] = ii;
            ijk[(ivar+2)%3] = jj;
            int idx_read = _hexa_ijk2idx(ijk[0], ijk[1], ijk[2], order);
            memcpy(child_node_coord + 3*idx_write, parent_node_coord + 3*idx_read, sizeof(double) * 3);
            idx_write++;
          }
        }
      }

      break;
    }

    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid elt_type %d\n", (int) parent_elt_type);

  }

  if (i_bnd >= 0) {
    *child_n_node = PDM_Mesh_nodal_n_vtx_elt_get(*child_elt_type, order);
  }

  return i_bnd;
}

#undef _tria_ij2idx
#undef _quad_ij2idx
#undef _tetra_ijk2idx
#undef _pyramid_ijk2idx
#undef _prism_ijk2idx
#undef _hexa_ijk2idx



static void
_boundary_child_uvw_to_parent_uvw
(
 const PDM_Mesh_nodal_elt_t  parent_elt_type,
 const int                   i_bnd,
 const double               *child_uvw,
       double               *parent_uvw
 )
 {
  switch (parent_elt_type) {
  case PDM_MESH_NODAL_BARHO:
  case PDM_MESH_NODAL_BARHO_BEZIER:
    if (i_bnd == 0) {
      parent_uvw[0] = 0;
    }
    else {
      parent_uvw[0] = 1;
    }
    break;
  case PDM_MESH_NODAL_TRIAHO:
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
    switch (i_bnd) {
    case 0:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = 0;
      break;
    case 1:
      parent_uvw[0] = 1 - child_uvw[0];
      parent_uvw[1] = child_uvw[0];
      break;
    case 2:
      parent_uvw[0] = 0;
      parent_uvw[1] = 1 - child_uvw[0];
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid i_bnd (%d) for elt_type %d\n", i_bnd, (int) parent_elt_type);
    }
    break;
  case PDM_MESH_NODAL_QUADHO:
    switch (i_bnd) {
    case 0:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = 0;
      break;
    case 1:
      parent_uvw[0] = 1;
      parent_uvw[1] = child_uvw[0];
      break;
    case 2:
      parent_uvw[0] = 1 - child_uvw[0];
      parent_uvw[1] = 1;
      break;
    case 3:
      parent_uvw[0] = 0;
      parent_uvw[1] = 1 - child_uvw[0];
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid i_bnd (%d) for elt_type %d\n", i_bnd, (int) parent_elt_type);
    }
    break;
  case PDM_MESH_NODAL_TETRAHO:
    switch (i_bnd) {
    case 0:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = 0;
      parent_uvw[2] = child_uvw[1];
      break;
    case 1:
      parent_uvw[0] = 1 - child_uvw[0] - child_uvw[1];
      parent_uvw[1] = child_uvw[0];
      parent_uvw[2] = child_uvw[1];
      break;
    case 2:
      parent_uvw[0] = 0;
      parent_uvw[1] = 1 - child_uvw[0];
      parent_uvw[2] = child_uvw[1];
      break;
    case 3:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = child_uvw[1];
      parent_uvw[2] = 0;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid i_bnd (%d) for elt_type %d\n", i_bnd, (int) parent_elt_type);
    }
    break;
  case PDM_MESH_NODAL_PYRAMIDHO:
    switch (i_bnd) {
    case 0:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = child_uvw[1];
      parent_uvw[2] = 0;
      break;
    case 1:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = 0;
      parent_uvw[2] = child_uvw[1];
      break;
    case 2:
      parent_uvw[0] = 1 - child_uvw[1];
      parent_uvw[1] = child_uvw[0];
      parent_uvw[2] = child_uvw[1];
      break;
    case 3:
      parent_uvw[0] = 1 - child_uvw[0] - child_uvw[1];
      parent_uvw[1] = 1 - child_uvw[1];
      parent_uvw[2] = child_uvw[1];
      break;
    case 4:
      parent_uvw[0] = 0;
      parent_uvw[1] = 1 - child_uvw[0] - child_uvw[1];
      parent_uvw[2] = child_uvw[1];
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid i_bnd (%d) for elt_type %d\n", i_bnd, (int) parent_elt_type);
    }
    break;
  case PDM_MESH_NODAL_PRISMHO:
    switch (i_bnd) {
    case 0:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = child_uvw[1];
      parent_uvw[2] = 0;
      break;
    case 1:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = child_uvw[1];
      parent_uvw[2] = 1;
      break;
    case 2:
      parent_uvw[0] = child_uvw[0];
      parent_uvw[1] = 0;
      parent_uvw[2] = child_uvw[1];
      break;
    case 3:
      parent_uvw[0] = 1 - child_uvw[0];
      parent_uvw[1] = child_uvw[0];
      parent_uvw[2] = child_uvw[1];
      break;
    case 4:
      parent_uvw[0] = 0;
      parent_uvw[1] = 1 - child_uvw[0];
      parent_uvw[2] = child_uvw[1];
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid i_bnd (%d) for elt_type %d\n", i_bnd, (int) parent_elt_type);
    }
    break;
  case PDM_MESH_NODAL_HEXAHO: {
    int ivar = i_bnd/2;
    parent_uvw[ivar      ] = i_bnd%2;
    parent_uvw[(ivar+1)%3] = child_uvw[0];
    parent_uvw[(ivar+2)%3] = child_uvw[1];
    break;
  }
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid elt_type %d\n", (int) parent_elt_type);
  }
}


static void
_init_uvw
(
 const PDM_Mesh_nodal_elt_t  elt_type,
       double               *uvw
 )
{
  switch (elt_type) {
  case PDM_MESH_NODAL_BARHO:
  case PDM_MESH_NODAL_BARHO_BEZIER:
    uvw[0] = 0.5;
    break;
  case PDM_MESH_NODAL_TRIAHO:
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
    uvw[0] = 1./3.;
    uvw[1] = uvw[0];
    break;
  case PDM_MESH_NODAL_QUADHO:
    uvw[0] = 0.5;
    uvw[1] = 0.5;
    break;
  case PDM_MESH_NODAL_TETRAHO:
    uvw[0] = 0.25;
    uvw[1] = 0.25;
    uvw[2] = 0.25;
    break;
  case PDM_MESH_NODAL_PYRAMIDHO:
    uvw[0] = 0.5;
    uvw[1] = 0.5;
    uvw[2] = 0.25;
    break;
  case PDM_MESH_NODAL_PRISMHO:
    uvw[0] = 1./3.;
    uvw[1] = uvw[0];
    uvw[2] = 0.25;
    break;
  case PDM_MESH_NODAL_HEXAHO:
    uvw[0] = 0.5;
    uvw[1] = 0.5;
    uvw[2] = 0.5;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid elt_type %d\n", (int) elt_type);
  }
}



static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

// static inline double
// _determinant_3x3
// (
//  const double a[3],
//  const double b[3],
//  const double c[3]
//  )
// {
//   return a[0] * (b[1]*c[2] - b[2]*c[1])
//     +    a[1] * (b[2]*c[0] - b[0]*c[2])
//     +    a[2] * (b[0]*c[1] - b[1]*c[0]);
// }

/*---------------------------------------------------------------------------
 * Solve the equation "A.x = b" with Cramer's rule.
 *
 * parameters:
 *   A[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   PDM_FALSE if matrix is singular, PDM_TRUE otherwise
 *----------------------------------------------------------------------------*/

static int
_solve_3x3(double *A,
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0]*(A[4]*A[8] - A[7]*A[5])
    -   A[3]*(A[1]*A[8] - A[7]*A[2])
    +   A[6]*(A[1]*A[5] - A[4]*A[2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return 1;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[4]*A[8] - A[7]*A[5])
          - b[1]*(A[1]*A[8] - A[7]*A[2])
          + b[2]*(A[1]*A[5] - A[4]*A[2])) * det_inv;

  x1 = (  A[0]*(b[1]*A[8] - b[2]*A[5])
          - A[3]*(b[0]*A[8] - b[2]*A[2])
          + A[6]*(b[0]*A[5] - b[1]*A[2])) * det_inv;

  x2 = (  A[0]*(A[4]*b[2] - A[7]*b[1])
          - A[3]*(A[1]*b[2] - A[7]*b[0])
          + A[6]*(A[1]*b[1] - A[4]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return 0;
}


static int
_compute_uvw_ho
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   order,
 const int                   n_node,
 const double                point_coords[3],
 const double                vertex_coords[],
 const double                tolerance,
       double               *uvw,
       double               *work_array
 )
{
  const int dbg_enabled = 0;

  // Take these 2 as input arguments to avoid redundant computations?
  int elt_dim = PDM_Mesh_nodal_elt_dim_get  (elt_type);

  assert(elt_dim > 0);
  // int n_node  = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);

  if (dbg_enabled) {
    log_trace(">> _compute_uvw_ho (type %d, dim %d)\n", (int) elt_type, elt_dim);
  }

  if (1 && dbg_enabled) {
    int *connec = malloc(sizeof(int) * n_node);

    int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                       elt_type,
                                                       order);

    for (int i = 0; i < n_node; i++) {
      connec[ijk_to_user[i]] = i + 1;
    }

    char filename[999];
    sprintf(filename, "_compute_uvw_ho_%d.vtk", (int) elt_type);

    PDM_vtk_write_std_elements_ho(filename,
                                  order,
                                  n_node,
                                  vertex_coords,
                                  NULL,
                                  elt_type,
                                  1,
                                  connec,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);
    free(connec);
  }

  double *weight     = work_array;
  double *diff_weight[3];
  for (int i = 0; i < elt_dim; i++) {
    diff_weight[i] = &work_array[(i+1)*n_node];
  }

  /* TO DO : use smarter finite difference step
   * - adapt to tolerance (or last Newton step size)
   * - beware of cancellation error amplified by division
   * - compute exact Jacobian matrix?
   */
  double finite_difference_dx = 1e-8; // hard-coded -_-'
  double inv_finite_difference_dx = 1./finite_difference_dx;


  const int    max_iter   = 15;
  const double tolerance2 = tolerance * tolerance;

  double _uvw[12];

  double a[9];
  double b[3];

  double duvw[3]; // Newton step in uvw-space

  _init_uvw(elt_type, uvw);

  // PDM_log_trace_array_double(uvw, elt_dim, "init uvw : ");

  PDM_Mesh_nodal_elt_t child_elt_type;
  int child_n_node;
  double child_node_coord[n_node*3];

  for (int iter = 0; iter < max_iter; iter++) {

    /* Evaluate position and approximate gradient */
    memcpy(_uvw, uvw, sizeof(double) * elt_dim);
    for (int i = 0; i < elt_dim; i++) {
      memcpy(_uvw + elt_dim*(i+1), uvw, sizeof(double) * elt_dim);
      _uvw[elt_dim*(i+1)+i] += finite_difference_dx;
    }

    PDM_ho_basis(elt_type,
                 order,
                 n_node,
                 elt_dim+1,
                 _uvw,
                 work_array);

    for (int j = 0; j < elt_dim; j++) {
      for (int i = 0; i < n_node; i++) {
        diff_weight[j][i] = (diff_weight[j][i] - weight[i]) * inv_finite_difference_dx;
      }
    }

    int singular = 0;

    for (int i = 0; i < 3; i++) {
      b[i] = point_coords[i];
    }
    for (int k = 0; k < n_node; k++) {
      for (int i = 0; i < 3; i++) {
        b[i] -= weight[k] * vertex_coords[3*k+i];
      }
    }
    for (int i = 0; i < 3*elt_dim; i++) {
      a[i] = 0;
    }

    for (int k = 0; k < n_node; k++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < elt_dim; j++) {
          a[elt_dim*i+j] += diff_weight[j][k] * vertex_coords[3*k+i];
        }
      }
    }

    switch (elt_dim) {
      case 3: {
        singular = _solve_3x3(a, b, duvw);
        break;
      }

      case 2: {
        double xuxu = a[0]*a[0] + a[2]*a[2] + a[4]*a[4];
        double xuxv = a[0]*a[1] + a[2]*a[3] + a[4]*a[5];
        double xvxv = a[1]*a[1] + a[3]*a[3] + a[5]*a[5];

        double denom = xuxu*xvxv - xuxv*xuxv;

        singular = (PDM_ABS(denom) < _epsilon_denom);

        if (!singular) {
          double rhs[2] = {
            a[0]*b[0] + a[2]*b[1] + a[4]*b[2],
            a[1]*b[0] + a[3]*b[1] + a[5]*b[2]
          };

          double idenom = 1./denom;
          duvw[0] = idenom * (rhs[0]*xvxv - rhs[1]*xuxv);
          duvw[1] = idenom * (rhs[1]*xuxu - rhs[0]*xuxv);
        }

        break;
      }

      case 1: {
        double aa = PDM_DOT_PRODUCT(a, a);

        singular = (aa < _epsilon_denom);

        if (!singular) {
          duvw[0] = PDM_DOT_PRODUCT(a, b) / aa;
        }

        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid elt dimension %d\n", elt_dim);
      }

    if (singular) {
      printf("_compute_uvw_ho : singularity\n");
      return -1;
    }

    double mag2_duvw = 0;
    for (int i = 0; i < elt_dim; i++) {
      mag2_duvw += duvw[i]*duvw[i];
      uvw[i] += duvw[i];
    }

    if (dbg_enabled) {
      log_trace("iter %2d, mag2_duvw %e\n", iter, mag2_duvw);
      PDM_log_trace_array_double(uvw, elt_dim, "uvw : ");
    }


    if (mag2_duvw < tolerance2) {
      int i_bnd = _boundary_child_element_coord(elt_type,
                                                order,
                                                uvw,
                                                vertex_coords,
                                                &child_elt_type,
                                                &child_n_node,
                                                child_node_coord);

      if (dbg_enabled) {
        log_trace("i_bnd = %d\n", i_bnd);
      }

      if (i_bnd < 0) {
        if (dbg_enabled) {
          log_trace("converged\n");
        }
        return 1;
      }
      else {
        if (child_elt_type == PDM_MESH_NODAL_POINT) {
          _boundary_child_uvw_to_parent_uvw(elt_type,
                                            i_bnd,
                                            NULL,
                                            uvw);
          if (dbg_enabled) {
            PDM_log_trace_array_double(uvw, elt_dim, "parent_uvw : ");
          }
          return 1;
        }

        double child_uvw[2];
        _init_uvw(child_elt_type, child_uvw);
        int child_converged = _compute_uvw_ho(child_elt_type,
                                              order,
                                              child_n_node,
                                              point_coords,
                                              child_node_coord,
                                              tolerance,
                                              child_uvw,
                                              work_array);

        if (child_converged) {
          if (dbg_enabled) {
            log_trace("child converged\n");
          }
          // Translate child_uvw to parent_uvw
          _boundary_child_uvw_to_parent_uvw(elt_type,
                                            i_bnd,
                                            child_uvw,
                                            uvw);
          if (dbg_enabled) {
            PDM_log_trace_array_double(child_uvw, elt_dim-1, "child_uvw  : ");
            PDM_log_trace_array_double(uvw,       elt_dim,   "parent_uvw : ");
          }
        }
        else {
          if (dbg_enabled) {
            log_trace("child did not converged\n");
          }
        }
        return child_converged;
      }
    }


  }

  if (dbg_enabled) {
    log_trace("failed to converge\n");
  }
  return 0;
}





/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Point location in a high order element
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   n_nodes           Number of nodes
 * \param [in]   nodes_coords      Coordinates of the nodes (size = 3 * \ref n_nodes)
 * \param [in]   point_coords      Coordinates of the point to locate (size = 3)
 * \param [out]  projected_coords  Coordinates of the projection on the element (size = 3)
 * \param [out]  uvw               Parametric coordinates of the projection on the element
 *
 * \return       Squared distance from the point to the element
 *
 */

double
PDM_ho_location
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw
 )
{
  // PDM_ho_location_user_elt_t *user_elt = *(_get_user_elt (type));
  PDM_ho_location_user_elt_t *user_elt = NULL;
  if (type > PDM_MESH_NODAL_HEXA8) {
    user_elt = *(_get_user_elt (type));
  }

  int entities_dim = PDM_Mesh_nodal_elt_dim_get(type);

  if (user_elt != NULL) {
    if (user_elt->location_in_elt != NULL) {
      return (user_elt->location_in_elt ) (entities_dim,
                                              order,
                                              n_nodes,
                                              nodes_coords,
                                              point_coords,
                                              projected_coords,
                                              uvw);
    }
    else {
      return _default_location (type,
                                order,
                                n_nodes,
                                nodes_coords,
                                point_coords,
                                projected_coords,
                                uvw);
    }
  }

  else {

    return _default_location (type,
                              order,
                              n_nodes,
                              nodes_coords,
                              point_coords,
                              projected_coords,
                              uvw);
  }

  return HUGE_VAL;
}



/**
 * \brief Get the parametric coordinates of the Lagrange nodes
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   umin              Lower bound on u-coordinate
 * \param [in]   umax              Upper bound on u-coordinate
 * \param [in]   vmin              Lower bound on v-coordinate
 * \param [in]   vmax              Upper bound on v-coordinate
 * \param [in]   wmin              Lower bound on w-coordinate
 * \param [in]   wmax              Upper bound on w-coordinate
 * \param [out]  uvw_node          Parametric coordinates of the Lagrange nodes
 *
 */

void
PDM_ho_location_uvw_nodes
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const double                umin,
 const double                umax,
 const double                vmin,
 const double                vmax,
 const double                wmin,
 const double                wmax,
 double                     *uvw_node
 )
{
  switch (type) {

    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER:
    _u_ho_edge_nodes(order, umin, umax, uvw_node);
    break;

    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
    _uv_ho_tria_nodes(order, umin, umax, vmin, vmax, uvw_node);
    break;

    case PDM_MESH_NODAL_QUADHO:
    _uv_ho_quad_nodes(order, umin, umax, vmin, vmax, uvw_node);
    break;

    case PDM_MESH_NODAL_TETRAHO:
    _uvw_ho_tetra_nodes(order, umin, umax, vmin, vmax, wmin, wmax, uvw_node);
    break;

    case PDM_MESH_NODAL_PYRAMIDHO:
    _uvw_ho_pyra_nodes(order, umin, umax, vmin, vmax, wmin, wmax, uvw_node);
    break;

    case PDM_MESH_NODAL_PRISMHO:
    _uvw_ho_prism_nodes(order, umin, umax, vmin, vmax, wmin, wmax, uvw_node);
    break;

    case PDM_MESH_NODAL_HEXAHO:
    _uvw_ho_hexa_nodes(order, umin, umax, vmin, vmax, wmin, wmax, uvw_node);
    break;

    default:
    PDM_error(__FILE__ , __LINE__, 0,
              "Invalid element type %d\n", type);
  }
}


/**
 * \brief Point location in a high order element
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   n_nodes           Number of nodes
 * \param [in]   nodes_coords      Coordinates of the nodes (size = 3 * \ref n_nodes)
 * \param [in]   point_coords      Coordinates of the point to locate (size = 3)
 * \param [in]   tolerance         Tolerance for Newton step in uvw-space
 * \param [out]  projected_coords  Coordinates of the projection on the element (size = 3)
 * \param [out]  uvw               Parametric coordinates of the projection on the element
 * \param [out]  converged         Convergence status (1 if successful, 0 else)
 * \param [in]   work_array        Optional work array (size = (elt_dim+1) * \ref n_nodes or NULL)
 *
 * \return       Distance from the point to the element
 *
 */

double
PDM_ho_location_newton
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 const double                tolerance,
       double               *projected_coords,
       double               *uvw,
       int                  *converged,
       double               *work_array
)
{
  double *_work_array = work_array;
  if (work_array == NULL) {
    int elt_dim = PDM_Mesh_nodal_elt_dim_get(type);
    _work_array = malloc(sizeof(double) * n_nodes * (elt_dim + 1));
  }
  int stat = _compute_uvw_ho(type,
                             order,
                             n_nodes,
                             point_coords,
                             nodes_coords,
                             tolerance,
                             uvw,
                             _work_array);

  *converged = (stat > 0);

  double dist2 = HUGE_VAL;
  if (*converged) {
    PDM_ho_basis(type,
                 order,
                 n_nodes,
                 1,
                 uvw,
                 _work_array);

    for (int j = 0; j < 3; j++) {
      projected_coords[j] = 0;
    }

    for (int i = 0; i < n_nodes; i++) {
      for (int j = 0; j < 3; j++) {
        projected_coords[j] += _work_array[i] * nodes_coords[3*i+j];
      }
    }

    dist2 = 0;
    for (int j = 0; j < 3; j++) {
      dist2 += (projected_coords[j] - point_coords[j])*(projected_coords[j] - point_coords[j]);
    }
  }

  if (work_array == NULL) {
    free(_work_array);
  }

  return dist2;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
