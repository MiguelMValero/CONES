/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_binary_search.h"
#include "pdm_vtk.h"
#include "pdm_geom_elem.h"

#include "pdm_sphere_vol_gen.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define ij2idx(i, j, n) ((i) + ((n)+1)*(j) - ((j)-1)*(j)/2)
#define ijk2idx(i, j, k, n) ((i) + (j)*((n) + 1 - (k)) - (j)*((j)-1)/2 + ((k)*((k)*((k) - 3*(n) - 6) + 3*(n)*((n) + 4) + 11)) / 6)

#define _point_on_edge(ibase, iedge, gnum_vtx, i, j, k) \
  int ibase_edge = base_cell_edge[6*(ibase) + (iedge)]; \
  int sign = PDM_SIGN(ibase_edge); \
  ibase_edge = PDM_ABS(ibase_edge) - 1; \
  int ei; \
  _local_edge_frame((iedge), (i), (j), (k), &ei); \
  if (sign < 0) { \
    (gnum_vtx) = idx_vtx_edge + ibase_edge*n + n - ei; \
  } \
  else { \
    (gnum_vtx) = idx_vtx_edge + ibase_edge*n + ei + 1; \
  }

#define _point_on_face(ibase, iface, gnum_vtx, i, j, k) \
  int ibase_face = base_cell_face[4*(ibase) + (iface)]; \
  int sign = PDM_SIGN(ibase_face); \
  ibase_face = PDM_ABS(ibase_face) - 1; \
  int perm = base_cell_face_perm[4*(ibase) + (iface)]; \
  int fi, fj; \
  _local_face_frame((iface), i, j, k, &fi, &fj); \
  _permute_ij(perm, sign, &fi, &fj, n-2); \
  (gnum_vtx) = 1 + idx_vtx_face + ibase_face*face_int_vtx_n + ij2idx(fi, fj, n-2);



#define _hextet_vtx_0(i, j, k, gnum_vtx) \
  if ((k) == 0) { \
    if ((j) == 0) { \
      if ((i) == 0) { \
        (gnum_vtx) = base_cell_vtx[4*ibase+0]; \
      } \
      else { \
        _point_on_edge(ibase, 0, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 1, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 3, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
    } \
  } \
  else { \
    if ((j) == 0) { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 2, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 2, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_face(ibase, 1, (gnum_vtx), (i)-1, (j)-1, (k)-1); \
      } \
      else { \
        (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i)-1, (j)-1, (k)-1, n-3); \
      } \
    } \
  }

#define _hextet_vtx_1(i, j, k, gnum_vtx) \
  if ((k) == 0) { \
    if ((j) == 0) { \
      if ((i) == n-(k)-(j)) { \
        (gnum_vtx) = base_cell_vtx[4*ibase+1]; \
      } \
      else { \
        _point_on_edge(ibase, 0, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
    } \
    else { \
      if ((i) == n-(k)-(j)) { \
        _point_on_edge(ibase, 3, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 3, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
    } \
  } \
  else { \
    if ((j) == 0) { \
      if ((i) == n-(k)-(j)) { \
        _point_on_edge(ibase, 4, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 2, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
    } \
    else { \
      if ((i) == n-(k)-(j)) { \
        _point_on_face(ibase, 0, (gnum_vtx), (i), (j)-1, (k)-1); \
      } \
      else { \
        (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i), (j)-1, (k)-1, n-3); \
      } \
    } \
  }

#define _hextet_vtx_2(i, j, k, gnum_vtx) \
  if ((k) == 0) { \
    if ((j) == n-(k)-(i)) { \
      if ((i) == 0) { \
        (gnum_vtx) = base_cell_vtx[4*ibase+2]; \
      } \
      else { \
        _point_on_edge(ibase, 3, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 1, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 3, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
    } \
  } \
  else { \
    if ((j) == n-(k)-i) { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 5, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
      else { \
        _point_on_face(ibase, 0, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_face(ibase, 1, (gnum_vtx), (i)-1, (j), (k)-1); \
      } \
      else { \
        (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i)-1, (j), (k)-1, n-3); \
      } \
    } \
  }

#define _hextet_vtx_3(i, j, k, gnum_vtx) \
  if ((k) == 0) {  \
    if ((i) == n-1-(k)-(j)) { \
      _point_on_edge(ibase, 3, (gnum_vtx), (i), (j), (k)-1); \
    } \
    else { \
      _point_on_face(ibase, 3, (gnum_vtx), (i), (j), (k)-1); \
    } \
  } \
  else { \
    if ((i) == n-1-(k)-(j)) { \
      _point_on_face(ibase, 0, (gnum_vtx), (i), (j), (k)-1); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i), (j), (k)-1, n-3); \
    } \
  }

#define _hextet_vtx_4(i, j, k, gnum_vtx) \
  if ((k) == n-(j)-(i)) { \
    if ((j) == 0) { \
      if ((i) == 0) { \
        (gnum_vtx) = base_cell_vtx[4*ibase+3]; \
      } \
      else { \
        _point_on_edge(ibase, 4, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 5, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
      else { \
        _point_on_face(ibase, 0, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
    } \
  } \
  else { \
    if ((j) == 0) { \
      if ((i) == 0) { \
        _point_on_edge(ibase, 2, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
      else { \
        _point_on_face(ibase, 2, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
    } \
    else { \
      if ((i) == 0) { \
        _point_on_face(ibase, 1, (gnum_vtx), (i)-1, (j)-1, (k)); \
      } \
      else { \
        (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i)-1, (j)-1, (k), n-3); \
      } \
    } \
  }

#define _hextet_vtx_5(i, j, k, gnum_vtx) \
  if ((j) == 0) { \
    if ((i) == n-1-(k)-(j)) { \
      _point_on_edge(ibase, 4, (gnum_vtx), (i), (j)-1, (k)); \
    } \
    else { \
      _point_on_face(ibase, 2, (gnum_vtx), (i), (j)-1, (k)); \
    } \
  } \
  else { \
    if ((i) == n-1-(k)-(j)) { \
      _point_on_face(ibase, 0, (gnum_vtx), (i), (j)-1, (k)); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i), (j)-1, (k), n-3); \
    } \
  }

#define _hextet_vtx_6(i, j, k, gnum_vtx) \
  if ((j) == n-1-(k)-i) { \
    if ((i) == 0) { \
      _point_on_edge(ibase, 5, (gnum_vtx), (i)-1, (j), (k)); \
    } \
    else { \
      _point_on_face(ibase, 0, (gnum_vtx), (i)-1, (j), (k)); \
    } \
  } \
  else { \
    if ((i) == 0) { \
      _point_on_face(ibase, 1, (gnum_vtx), (i)-1, (j), (k)); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i)-1, (j), (k), n-3); \
    } \
  }

#define _hextet_vtx_7(i, j, k, gnum_vtx) \
  if ((i) == n-2-(k)-(j)) { \
    _point_on_face(ibase, 0, (gnum_vtx), (i), (j), (k)); \
  } \
  else { \
    (gnum_vtx) = 1 + idx_vtx_cell + ibase*cell_int_vtx_n + ijk2idx((i), (j), (k), n-3); \
  }


#define _point_on_edge2(ibase, iedge, gnum_vtx, i, j) \
  int ibase_edge = base_face_edge[3*(ibase) + (iedge)]; \
  int sign = PDM_SIGN(ibase_edge); \
  ibase_edge = PDM_ABS(ibase_edge) - 1; \
  int ei; \
  _local_edge_frame2((iedge), (i), (j), n, &ei); \
  if (sign < 0) { \
    (gnum_vtx) = idx_vtx_edge + ibase_edge*n + n - ei; \
  } \
  else { \
    (gnum_vtx) = idx_vtx_edge + ibase_edge*n + ei + 1; \
  }

#define _quadtria_vtx_0(i, j, gnum_vtx) \
  if ((j) == 0) { \
    if ((i) == 0) { \
      (gnum_vtx) = base_face_vtx[3*ibase_face]; \
    } \
    else { \
      _point_on_edge2(ibase_face, 0, (gnum_vtx), (i)-1, (j)-1); \
    } \
  } \
  else { \
    if ((i) == 0) { \
      _point_on_edge2(ibase_face, 2, (gnum_vtx), (i)-1, (j)-1); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_face + face_int_vtx_n*ibase_face + ij2idx((i)-1, (j)-1, n-2); \
    } \
  }


#define _quadtria_vtx_1(i, j, gnum_vtx) \
  if ((j) == 0) { \
    if ((i) == n-(j)) { \
      (gnum_vtx) = base_face_vtx[3*ibase_face+1]; \
    } \
    else { \
      _point_on_edge2(ibase_face, 0, (gnum_vtx), (i), (j)-1); \
    } \
  } \
  else { \
    if ((i) == n-(j)) { \
      _point_on_edge2(ibase_face, 1, (gnum_vtx), (i), (j)-1); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_face + face_int_vtx_n*ibase_face + ij2idx((i), (j)-1, n-2); \
    } \
  }


#define _quadtria_vtx_2(i, j, gnum_vtx) \
  if ((j) == n-i) { \
    if ((i) == 0) { \
      (gnum_vtx) = base_face_vtx[3*ibase_face+2]; \
    } \
    else { \
      _point_on_edge2(ibase_face, 1, (gnum_vtx), (i)-1, (j)); \
    } \
  } \
  else { \
    if ((i) == 0) { \
      _point_on_edge2(ibase_face, 2, (gnum_vtx), (i)-1, (j)); \
    } \
    else { \
      (gnum_vtx) = 1 + idx_vtx_face + face_int_vtx_n*ibase_face + ij2idx((i)-1, (j), n-2); \
    } \
  }

#define _quadtria_vtx_3(i, j, gnum_vtx) \
  if ((j) == n-i-1) { \
    _point_on_edge2(ibase_face, 1, (gnum_vtx), (i), (j)); \
  } \
  else { \
    (gnum_vtx) = 1 + idx_vtx_face + face_int_vtx_n*ibase_face + ij2idx((i), (j), n-2); \
  }

/*============================================================================
 * Private function definitions
 *============================================================================*/


static inline double
_smooth_step
(
 const double x
 )
{
  return PDM_MIN(1., PDM_MAX(0., 3*x*x*(2 - x)));
}

/**
 *
 * -1 <= xc, yc, zc <= 1
 *
 */

static void
_cube_to_sphere
(
 const double  xc,
 const double  yc,
 const double  zc,
 double       *xs,
 double       *ys,
 double       *zs
 )
{
  double r2 = xc*xc + yc*yc + zc*zc;
  double r  = sqrt(r2);
  double rinf = PDM_MAX(PDM_MAX(PDM_ABS(xc), PDM_ABS(yc)), PDM_ABS(zc));

  double scale1 = rinf;
  if (r > 1e-15) {
    scale1 /= r;
  }
  double x1 = xc * scale1;
  double y1 = yc * scale1;
  double z1 = zc * scale1;

  double scale2 = sqrt(r2 -
                       xc*xc*yc*yc - yc*yc*zc*zc - zc*zc*xc*xc +
                       xc*xc*yc*yc*zc*zc);
  if (r > 1e-15) {
    scale2 /= r;
  }
  double x2 = xc * scale2;
  double y2 = yc * scale2;
  double z2 = zc * scale2;

  double t = 0;
  if(0) {
    t = 0.7*_smooth_step(r);
  }

  *xs = (1 - t)*x2 + t*x1;
  *ys = (1 - t)*y2 + t*y1;
  *zs = (1 - t)*z2 + t*z1;
}


static inline void
idx2ij
(
 const int  idx,
 const int  n,
 int       *i,
 int       *j
 )
{
  int b = -(2*n + 3);
  int d = b*b - 8*idx;
  int _j = (int) (0.5 * (-b - sqrt(d)));

  *i = idx - ij2idx(0, _j, n);
  *j = _j;
}


static int tetra_face_vtx[4][3] = {
  {1, 2, 3},
  {0, 3, 2},
  {0, 1, 3},
  {0, 2, 1}//{2, 0, 1}
};

static int tetra_edge_vtx[6][2] = {
  {0, 1},
  {0, 2},
  {0, 3},
  {1, 2},
  {1, 3},
  {2, 3}
};

static int tria_edge_vtx[3][2] = {
  {0, 1}, //{1, 2},
  {1, 2}, //{2, 0},
  {2, 0}  //{0, 1}
};

static int _has_same_vtx_triangle
(
 const int tri1[3],
 const int tri2[3]
 )
{
  // We know that sum(tri1) == sum(tri2)
  if (tri1[0] == tri2[0]) {
    if (tri1[1] == tri2[1]) {
      return 1;
    } else if (tri1[1] == tri2[2]) {
      return -1;
    } else {
      return 0;
    }
  }

  else if (tri1[0] == tri2[1]) {
    if (tri1[1] == tri2[2]) {
      return 1;
    } else if (tri1[1] == tri2[0]) {
      return -1;
    } else {
      return 0;
    }
  }

  else if (tri1[0] == tri2[2]) {
    if (tri1[1] == tri2[0]) {
      return 1;
    } else if (tri1[1] == tri2[1]) {
      return -1;
    } else {
      return 0;
    }
  }

  return 0;
}


static int _has_same_vtx_edge
(
 const int edg1[2],
 const int edg2[2]
 )
{
  // We know that sum(edg1) == sum(edg2)
  if (edg1[0] == edg2[0]) {
    if (edg1[1] == edg2[1]) {
      return 1;
    } else {
      return 0;
    }
  }

  else if (edg1[0] == edg2[1]) {
    if (edg1[1] == edg2[0]) {
      return -1;
    } else {
      return 0;
    }
  }

  return 0;
}


static void
_build_edges_from_triangles
(
const int   n_vtx,
const int   n_face,
const int  *face_vtx,
      int  *n_edge,
      int **edge_vtx,
      int **face_edge
 )
{
  *face_edge = malloc(sizeof(int) * n_face * 3);

  int key_max = PDM_MAX(0, 2*(n_vtx-1));
  int n_key = key_max + 1;
  int *key_edge_n = PDM_array_zeros_int(n_key);

  for (int iface = 0; iface < n_face; iface++) {
    const int *_face_vtx = face_vtx + 3*iface;

    for (int idx_edge = 0; idx_edge < 3; idx_edge++) {
      int ev[2];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
        ev[idx_vtx] = _face_vtx[tria_edge_vtx[idx_edge][idx_vtx]];
        key += ev[idx_vtx] - 1;
      }
      key_edge_n[key]++;
    }
  }

  int *key_edge_idx = PDM_array_new_idx_from_sizes_int(key_edge_n, n_key);
  PDM_array_reset_int(key_edge_n, n_key, 0);

  int *key_edge = (int *) malloc(sizeof(int) * key_edge_idx[n_key]);

  int n_edge_max = 3*n_face;
  *n_edge = 0;
  *edge_vtx  = malloc(sizeof(int) * n_edge_max * 2);
  for (int iface = 0; iface < n_face; iface++) {
    const int *_face_vtx  =  face_vtx  + 3*iface;
    int       *_face_edge = *face_edge + 3*iface;

    for (int idx_edge = 0; idx_edge < 3; idx_edge++) {
      int ev[2];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
        ev[idx_vtx] = _face_vtx[tria_edge_vtx[idx_edge][idx_vtx]];
        key += ev[idx_vtx] - 1;
      }


      int edge_id = 0;
      for (int idx = 0; idx < key_edge_n[key]; idx++) {
        int jedge = key_edge[key_edge_idx[key] + idx] - 1;
        int *ev2 = (*edge_vtx) + 2*jedge;
        int sign = _has_same_vtx_edge(ev, ev2);

        if (sign != 0) {
          edge_id = sign*(jedge+1);
          break;
        }
      }

      if (edge_id == 0) {
        // New edge
        (*edge_vtx)[2*(*n_edge)  ] = ev[0];
        (*edge_vtx)[2*(*n_edge)+1] = ev[1];
        edge_id = ++(*n_edge);

        key_edge[key_edge_idx[key] + key_edge_n[key]++] = edge_id;
      }

      _face_edge[idx_edge] = edge_id;
    }
  }
  free(key_edge_n);
  free(key_edge_idx);
  free(key_edge);

  *edge_vtx = realloc(*edge_vtx, sizeof(int) * (*n_edge) * 2);
}

static void
_build_base_edges_and_faces_from_tetra
(
const int   n_vtx,
const int   n_cell,
const int  *cell_vtx,
      int  *n_edge,
      int  *n_face,
      int **edge_vtx,
      int **face_vtx,
      int **face_edge,
      int **cell_face,
      int **cell_edge
 )
{
  *cell_face = malloc(sizeof(int) * n_cell * 4);
  *cell_edge = malloc(sizeof(int) * n_cell * 6);

  /* Build faces */
  int key_max = PDM_MAX(0, 3*(n_vtx-1));
  int n_key = key_max + 1;
  int *key_face_n = PDM_array_zeros_int(n_key);

  for (int icell = 0; icell < n_cell; icell++) {
    const int *_cell_vtx = cell_vtx + 4*icell;

    for (int idx_face = 0; idx_face < 4; idx_face++) {
      int fv[3];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
        fv[idx_vtx] = _cell_vtx[tetra_face_vtx[idx_face][idx_vtx]];
        key += fv[idx_vtx] - 1;
      }
      key_face_n[key]++;
    }
  }

  int *key_face_idx = PDM_array_new_idx_from_sizes_int(key_face_n, n_key);
  PDM_array_reset_int(key_face_n, n_key, 0);

  int *key_face = (int *) malloc(sizeof(int) * key_face_idx[n_key]);

  int n_face_max = 4*n_cell;
  *n_face = 0;
  *face_vtx  = malloc(sizeof(int) * n_face_max * 3);
  for (int icell = 0; icell < n_cell; icell++) {
    const int *_cell_vtx  =  cell_vtx  + 4*icell;
    int       *_cell_face = *cell_face + 4*icell;

    for (int idx_face = 0; idx_face < 4; idx_face++) {
      int fv[3];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
        fv[idx_vtx] = _cell_vtx[tetra_face_vtx[idx_face][idx_vtx]];
        key += fv[idx_vtx] - 1;
      }


      // log_trace("face #%d from cell %d : %d %d %d, key = %d\n",
      //           idx_face, icell, fv[0], fv[1], fv[2], key);

      int face_id = 0;
      for (int idx = 0; idx < key_face_n[key]; idx++) {
        int jface = key_face[key_face_idx[key] + idx] - 1;
        int *fv2 = (*face_vtx) + 3*jface;
        // log_trace("  jface %d: %d %d %d\n", jface, fv2[0], fv2[1], fv2[2]);
        int sign = _has_same_vtx_triangle(fv, fv2);

        if (sign != 0) {
          face_id = sign*(jface+1);
          // log_trace(" >> existing face %d\n", face_id);
          break;
        }
      }

      if (face_id == 0) {
        // New face
        (*face_vtx)[3*(*n_face)  ] = fv[0];
        (*face_vtx)[3*(*n_face)+1] = fv[1];
        (*face_vtx)[3*(*n_face)+2] = fv[2];
        face_id = ++(*n_face);
        // log_trace(" >> new face %d\n", face_id);

        key_face[key_face_idx[key] + key_face_n[key]++] = face_id;
      }

      _cell_face[idx_face] = face_id;
    }
  }
  free(key_face_n);
  free(key_face_idx);
  free(key_face);

  *face_vtx  = realloc(*face_vtx, sizeof(int) * (*n_face) * 3);
  *face_edge = malloc(sizeof(int) * (*n_face) * 3);

  /* Build edges */
  key_max = PDM_MAX(0, 2*(n_vtx-1));
  n_key = key_max + 1;
  int *key_edge_n = PDM_array_zeros_int(n_key);

  for (int iface = 0; iface < *n_face; iface++) {
     int *_face_vtx = *face_vtx + 3*iface;

    for (int idx_edge = 0; idx_edge < 3; idx_edge++) {
      int ev[2];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
        ev[idx_vtx] = _face_vtx[tria_edge_vtx[idx_edge][idx_vtx]];
        key += ev[idx_vtx] - 1;
      }
      key_edge_n[key]++;
    }
  }

  int *key_edge_idx = PDM_array_new_idx_from_sizes_int(key_edge_n, n_key);
  PDM_array_reset_int(key_edge_n, n_key, 0);

  int *key_edge = (int *) malloc(sizeof(int) * key_edge_idx[n_key]);

  int n_edge_max = 3*(*n_face);
  *n_edge = 0;
  *edge_vtx  = malloc(sizeof(int) * n_edge_max * 2);
  for (int iface = 0; iface < *n_face; iface++) {
    int *_face_vtx  = *face_vtx  + 3*iface;
    int *_face_edge = *face_edge + 3*iface;

    for (int idx_edge = 0; idx_edge < 3; idx_edge++) {
      int ev[2];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
        ev[idx_vtx] = _face_vtx[tria_edge_vtx[idx_edge][idx_vtx]];
        key += ev[idx_vtx] - 1;
      }


      int edge_id = 0;
      for (int idx = 0; idx < key_edge_n[key]; idx++) {
        int jedge = key_edge[key_edge_idx[key] + idx] - 1;
        int *ev2 = (*edge_vtx) + 2*jedge;
        int sign = _has_same_vtx_edge(ev, ev2);

        if (sign != 0) {
          edge_id = sign*(jedge+1);
          break;
        }
      }

      if (edge_id == 0) {
        // New edge
        (*edge_vtx)[2*(*n_edge)  ] = ev[0];
        (*edge_vtx)[2*(*n_edge)+1] = ev[1];
        edge_id = ++(*n_edge);

        key_edge[key_edge_idx[key] + key_edge_n[key]++] = edge_id;
      }

      _face_edge[idx_edge] = edge_id;
    }
  }

  *edge_vtx = realloc(*edge_vtx, sizeof(int) * (*n_edge) * 2);



  for (int icell = 0; icell < n_cell; icell++) {
    const int *_cell_vtx  =  cell_vtx  + 4*icell;
    int       *_cell_edge = *cell_edge + 6*icell;

    for (int idx_edge = 0; idx_edge < 6; idx_edge++) {
      int ev[2];
      int key = 0;
      for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
        ev[idx_vtx] = _cell_vtx[tetra_edge_vtx[idx_edge][idx_vtx]];
        key += ev[idx_vtx] - 1;
      }


      int edge_id = 0;
      for (int idx = 0; idx < key_edge_n[key]; idx++) {
        int jedge = key_edge[key_edge_idx[key] + idx] - 1;
        int *ev2 = (*edge_vtx) + 2*jedge;
        int sign = _has_same_vtx_edge(ev, ev2);

        if (sign != 0) {
          edge_id = sign*(jedge+1);
          break;
        }
      }

      if (edge_id == 0) {
        PDM_error(__FILE__, __LINE__, 0, "Error cell_edge edge not in a face\n");
      }

      _cell_edge[idx_edge] = edge_id;
    }
  }
  free(key_edge_n);
  free(key_edge_idx);
  free(key_edge);

}



static inline int
_n_subtetra
(
 const int n
 )
{
  return (n+1)*(n+2)*(n+3)/6;
}


static inline int
_n_subtria
(
 const int n
 )
{
  return (n+1)*(n+2)/2;
}



static inline void
_permute_ij
(
 const int  perm,
 const int  sign,
       int *i,
       int *j,
 const int  n
 )
{
  if (sign > 0) {
    if (perm == 0) {
      // do nothing
    }
    else if (perm == 1) {
      int tmp = *i;
      *i = *j;
      *j = n - tmp;
    }
    else {
      int tmp = *i;
      *i = n - (*j);
      *j = tmp;
    }
  }
  else {
    if (perm == 0) {
      int tmp = *i;
      *i = *j;
      *j = tmp;
    }
    else if (perm == 1) {
      *i = n - (*i);
    }
    else {
      *j = n - (*j);
    }
  }
}


static inline void
_local_edge_frame
(
const int  iedge,
const int  i,
const int  j,
const int  k,
      int *ei
)
{
  if (iedge == 0) {
    *ei = i;
  }
  else if (iedge == 1) {
    *ei = j;
  }
  else if (iedge == 2) {
    *ei = k;
  }
  else if (iedge == 3) {
    *ei = j;
  }
  else if (iedge == 4) {
    *ei = k;
  }
  else {
    *ei = k;
  }
}


static inline void
_local_edge_frame2
(
const int  iedge,
const int  i,
const int  j,
const int  n,
      int *ei
)
{
  if (iedge == 0) {
    *ei = i;
  }
  else if (iedge == 1) {
    *ei = j;
  }
  else {
    *ei = n-j-1;
  }
}


static inline void
_local_face_frame
(
const int  iface,
const int  i,
const int  j,
const int  k,
      int *fi,
      int *fj
)
{
  if (iface == 0) {
    *fi = j;
    *fj = k;
  }
  else if (iface == 1) {
    *fi = k;
    *fj = j;
  }
  else if (iface == 2) {
    *fi = i;
    *fj = k;
  }
  else {
    *fi = j;
    *fj = i;
  }
}




static void
_gen_from_base_mesh
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t    n,
 const int            base_n_vtx,
 const int            base_n_edge,
 const int            base_n_face,
 const int            base_n_cell,
       double        *base_vtx_coord,
       int           *base_edge_vtx,
       int           *base_face_vtx,
       int           *base_face_edge,
       int           *base_cell_face,
       int           *base_cell_edge,
       int           *base_cell_vtx,
       double       **dvtx_coord,
       PDM_g_num_t  **dface_vtx,
       PDM_g_num_t  **dcell_vtx,
       PDM_g_num_t  **distrib_vtx,
       PDM_g_num_t  **distrib_face,
       PDM_g_num_t  **distrib_cell
)
{
  assert(n >= 0);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /* Check orientation */
  for (int i = 0; i < base_n_cell; i++) {

    int *bcv = base_cell_vtx + 4*i;

    double u[3], v[3], w[3];
    for (int j = 0; j < 3; j++) {
      u[j] = base_vtx_coord[3*(bcv[1]-1)+j] - base_vtx_coord[3*(bcv[0]-1)+j];
      v[j] = base_vtx_coord[3*(bcv[2]-1)+j] - base_vtx_coord[3*(bcv[0]-1)+j];
      w[j] = base_vtx_coord[3*(bcv[3]-1)+j] - base_vtx_coord[3*(bcv[0]-1)+j];
    }

    double uv[3];
    PDM_CROSS_PRODUCT(uv, u, v);
    double vol6 = PDM_DOT_PRODUCT(uv, w);

    if (vol6 <= 0) {
      log_trace("!!! base cell %d has volume = %e\n", i, vol6/6.);
    }
  }

  /* Auxiliary arrays */
  int *tria_j_idx = malloc(sizeof(int) * n);
  if (n > 0) {
    tria_j_idx[0] = 0;
    for (int j = 0; j < n-1; j++) {
      tria_j_idx[j+1] = tria_j_idx[j] + n-1-j;
    }
  }

  // PDM_log_trace_array_int(tria_j_idx, n, "tria_j_idx : ");

  int **tetra_j_idx = malloc(sizeof(int *) * PDM_MAX(n-2, 0));
  int  *tetra_k_idx = malloc(sizeof(int  ) * PDM_MAX(n-1, 0));

  if (n > 1) {
    for (int k = 0; k < n-2; k++) {
      tetra_j_idx[k] = malloc(sizeof(int) * (n-1-k));
      tetra_j_idx[k][0] = 0;
      for (int j = 0; j < n-2-k; j++) {
        tetra_j_idx[k][j+1] = tetra_j_idx[k][j] + n-2-k-j;
      }

      // log_trace("tetra_j_idx[%d] : ", k);
      // PDM_log_trace_array_int(tetra_j_idx[k], n-1-k, "");
    }

    tetra_k_idx[0] = 0;
    for (int k = 0; k < n-2; k++) {
      int p = n-3-k;
      tetra_k_idx[k+1] = tetra_k_idx[k] + (p+1)*(p+2)/2;
    }
  }

  // PDM_log_trace_array_int(tetra_k_idx, n-1, "tetra_k_idx : ");

  /*
   *  Vertices
   */
  PDM_g_num_t face_int_vtx_n = n*(n-1) / 2;
  PDM_g_num_t cell_int_vtx_n = n*(n-1)*(n-2) / 6;

  PDM_g_num_t gn_vtx =
  base_n_vtx                 +
  base_n_edge*n              +
  base_n_face*face_int_vtx_n +
  base_n_cell*cell_int_vtx_n;

  *distrib_vtx = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  int dn_vtx = (int) ((*distrib_vtx)[i_rank+1] - (*distrib_vtx)[i_rank]);

  *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);

  PDM_g_num_t idx_vtx_edge = base_n_vtx;
  PDM_g_num_t idx_vtx_face = idx_vtx_edge + base_n_edge*n;
  PDM_g_num_t idx_vtx_cell = idx_vtx_face + base_n_face*face_int_vtx_n;

  // log_trace("idx_vtx_edge/face/cell = "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
  //           idx_vtx_edge, idx_vtx_face, idx_vtx_cell);
  // log_trace("gn_vtx = "PDM_FMT_G_NUM"\n", gn_vtx);

  double step = 1. / (double) (n + 1);

  for (int ivtx = 0; ivtx < dn_vtx; ivtx++) {

    PDM_g_num_t g = (*distrib_vtx)[i_rank] + ivtx;

    if (g < idx_vtx_edge) {
      // Base vertex
      memcpy(*dvtx_coord + 3*ivtx, base_vtx_coord + 3*g, sizeof(double) * 3);
    }

    else if (g < idx_vtx_face) {
      // Base edge
      int ibase = (int) ((g - idx_vtx_edge) / n);
      int i     = (int) (g - idx_vtx_edge - n*ibase);
      // log_trace("base edge %d, i = %d --> gnum "PDM_FMT_G_NUM"\n",
      //           ibase, i, g+1);

      int ivtx1 = base_edge_vtx[2*ibase  ] - 1;
      int ivtx2 = base_edge_vtx[2*ibase+1] - 1;
      double t = (i + 1) * step;

      for (int k = 0; k < 3; k++) {
        (*dvtx_coord)[3*ivtx + k] = (1 - t)*base_vtx_coord[3*ivtx1 + k] + t*base_vtx_coord[3*ivtx2 + k];
      }
    }

    else if (g < idx_vtx_cell) {
      // Base face
      int ibase = (int) ((g - idx_vtx_face) / face_int_vtx_n);
      int idx   = (int) ( g - idx_vtx_face  - face_int_vtx_n*ibase);
      int j = PDM_binary_search_gap_int(idx,
                                        tria_j_idx,
                                        n);
      int i = idx - tria_j_idx[j];
      // log_trace("base face %d, ij = %d %d --> gnum "PDM_FMT_G_NUM"\n",
      //           ibase, i, j, g+1);

      int ivtx1 = base_face_vtx[3*ibase  ] - 1;
      int ivtx2 = base_face_vtx[3*ibase+1] - 1;
      int ivtx3 = base_face_vtx[3*ibase+2] - 1;
      double u = (i + 1) * step;
      double v = (j + 1) * step;

      for (int k = 0; k < 3; k++) {
        (*dvtx_coord)[3*ivtx + k] =
        (1 - u - v) * base_vtx_coord[3*ivtx1 + k] +
        u           * base_vtx_coord[3*ivtx2 + k] +
        v           * base_vtx_coord[3*ivtx3 + k];
      }
    }

    else {
      // Base cell
      int ibase = (int) ((g - idx_vtx_cell) / cell_int_vtx_n);
      int idx   = (int)  (g - idx_vtx_cell  - cell_int_vtx_n*ibase);
      int k = PDM_binary_search_gap_int(idx,
                                        tetra_k_idx,
                                        n-1);
      int j = PDM_binary_search_gap_int(idx - tetra_k_idx[k],
                                        tetra_j_idx[k],
                                        n-1-k);
      int i = idx - tetra_k_idx[k] - tetra_j_idx[k][j];
      // log_trace("base cell %d, idx = %d, ijk = %d %d %d --> gnum "PDM_FMT_G_NUM"\n",
      //           ibase, i, j, k, g+1);

      int ivtx1 = base_cell_vtx[4*ibase  ] - 1;
      int ivtx2 = base_cell_vtx[4*ibase+1] - 1;
      int ivtx3 = base_cell_vtx[4*ibase+2] - 1;
      int ivtx4 = base_cell_vtx[4*ibase+3] - 1;
      double u = (i + 1) * step;
      double v = (j + 1) * step;
      double w = (k + 1) * step;

      for (int l = 0; l < 3; l++) {
        (*dvtx_coord)[3*ivtx + l] =
        (1 - u - v - w) * base_vtx_coord[3*ivtx1 + l] +
        u               * base_vtx_coord[3*ivtx2 + l] +
        v               * base_vtx_coord[3*ivtx3 + l] +
        w               * base_vtx_coord[3*ivtx4 + l];
      }
    }

  } // End of loop on vertices
  free(tria_j_idx);
  for (int k = 0; k < n-2; k++) {
    free(tetra_j_idx[k]);
  }
  free(tetra_j_idx);
  free(tetra_k_idx);




  /*
   *  Cells
   */
  PDM_g_num_t cell_subcell_n = (n+1)*(n+1)*(n+1);

  PDM_g_num_t gn_cell = base_n_cell * cell_subcell_n;

  int hextet_size[6] = {
    n+1, n, n, n, n, n-1
  };

  int subcell_hextet_n[6] = {
    _n_subtetra(n),
    _n_subtetra(n-1),
    _n_subtetra(n-1),
    _n_subtetra(n-1),
    _n_subtetra(n-1),
    _n_subtetra(n-2)
  };

  int subcell_hextet_idx[7];
  PDM_array_idx_from_sizes_int(subcell_hextet_n,
                               6,
                               subcell_hextet_idx);

  // PDM_log_trace_array_int(subcell_hextet_idx,
  //                         7,
  //                         "subcell_hextet_idx : ");

  *distrib_cell = PDM_compute_uniform_entity_distribution(comm, gn_cell);

  int dn_cell = (int) ((*distrib_cell)[i_rank+1] - (*distrib_cell)[i_rank]);

  *dcell_vtx = malloc(sizeof(PDM_g_num_t) * dn_cell * 4);

  // //-->> tmp
  // for (int i = 0; i < 4*dn_cell; i++) {
  //   (*dcell_vtx)[i] = 1;
  // }
  // //<<--


  int *base_cell_face_perm = malloc(sizeof(int) * 4 * base_n_cell);
  for (int icell = 0; icell < base_n_cell; icell++) {
    for (int idx_face = 0; idx_face < 4; idx_face++) {
      int face_id = PDM_ABS(base_cell_face[4*icell+idx_face]) - 1;
      for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
        int vtx_id = base_cell_vtx[4*icell + tetra_face_vtx[idx_face][idx_vtx]];
        if (vtx_id == base_face_vtx[3*face_id]) {
          base_cell_face_perm[4*icell+idx_face] = idx_vtx;
          break;
        }
      }
    }
  }




  int *hextet_k_idx[6];
  hextet_k_idx[0] = malloc(sizeof(int) * (n+2));
  hextet_k_idx[0][0] = 0;
  for (int k = 0; k < n+1; k++) {
    int p = n-k;
    hextet_k_idx[0][k+1] = hextet_k_idx[0][k] + (p+1)*(p+2)/2;
  }
  assert(hextet_k_idx[0][n+1] == subcell_hextet_n[0]);

  hextet_k_idx[1] = malloc(sizeof(int) * (n+1));
  hextet_k_idx[1][0] = 0;
  for (int k = 0; k < n; k++) {
    int p = n-1-k;
    hextet_k_idx[1][k+1] = hextet_k_idx[1][k] + (p+1)*(p+2)/2;
  }
  assert(hextet_k_idx[1][n] == subcell_hextet_n[1]);

  hextet_k_idx[2] = hextet_k_idx[1];
  hextet_k_idx[3] = hextet_k_idx[1];
  hextet_k_idx[4] = hextet_k_idx[1];

  hextet_k_idx[5] = malloc(sizeof(int) * n);
  if (n > 0) {
    hextet_k_idx[5][0] = 0;
    for (int k = 0; k < n-1; k++) {
      int p = n-2-k;
      hextet_k_idx[5][k+1] = hextet_k_idx[5][k] + (p+1)*(p+2)/2;
    }
    assert(hextet_k_idx[5][n-1] == subcell_hextet_n[5]);
  }


  for (int icell = 0; icell < dn_cell; icell++) {

    PDM_g_num_t g = (*distrib_cell)[i_rank] + icell;
    // log_trace("icell = %d/%d, g = "PDM_FMT_G_NUM"\n", icell, dn_cell, g);

    PDM_g_num_t *_dcell_vtx = *dcell_vtx + 4*icell;

    int ibase = (int) (g / cell_subcell_n);
    int idx   = (int) (g - cell_subcell_n*ibase);

    // log_trace("  base cell %d / %d\n", ibase, base_n_cell);

    int hextet = PDM_binary_search_gap_int(idx,
                                         subcell_hextet_idx,
                                         7);

    idx -= subcell_hextet_idx[hextet];

    int k = PDM_binary_search_gap_int(idx,
                                      hextet_k_idx[hextet],
                                      hextet_size[hextet] + 1);
    // log_trace("  hextet = %d, idx = %d, k = %d\n", hextet, idx, k);
    int j, i;
    idx2ij(idx - hextet_k_idx[hextet][k], hextet_size[hextet]-1-k, &i, &j);

    // log_trace("  hextet %d, idx = %5d, ijk = %3d %3d %3d\n",
    //           hextet, idx, i, j, k);


    switch (hextet) {

      case 0: {
        // hextet_vtx : [0, 1, 2, 4]

        // vertex #0
        _hextet_vtx_0(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_1(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_2(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_4(i, j, k, _dcell_vtx[3]);

        break;
      }

      case 1: {
        // hextet_vtx : [1, 5, 3, 4]

        // vertex #0
        _hextet_vtx_1(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_5(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_3(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_4(i, j, k, _dcell_vtx[3]);

        break;
      }

      case 2: {
        // hextet_vtx : [1, 3, 2, 4]

        // vertex #0
        _hextet_vtx_1(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_3(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_2(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_4(i, j, k, _dcell_vtx[3]);

        break;
      }

      case 3: {
        // hextet_vtx : [2, 3, 6, 4]

        // vertex #0
        _hextet_vtx_2(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_3(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_6(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_4(i, j, k, _dcell_vtx[3]);

        break;
      }

      case 4: {
        // hextet_vtx : [3, 4, 5, 6]

        // vertex #0
        _hextet_vtx_3(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_4(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_5(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_6(i, j, k, _dcell_vtx[3]);

        break;
      }

      case 5: {
        // hextet_vtx : [3, 7, 6, 5]

        // vertex #0
        _hextet_vtx_3(i, j, k, _dcell_vtx[0]);

        // vertex #1
        _hextet_vtx_7(i, j, k, _dcell_vtx[1]);

        // vertex #2
        _hextet_vtx_6(i, j, k, _dcell_vtx[2]);

        // vertex #3
        _hextet_vtx_5(i, j, k, _dcell_vtx[3]);

        break;
      }

      default: {
        PDM_error(__FILE__, __LINE__, 0, "Wrong subcell hextet %d\n", hextet);
      }
    }

    if (0) {
      /* Flip */
      PDM_g_num_t tmp = _dcell_vtx[0];
      _dcell_vtx[0] = _dcell_vtx[1];
      _dcell_vtx[1] = tmp;
    }

    // PDM_log_trace_array_long(_dcell_vtx, 4, "_dcell_vtx : ");
    for (int ivtx = 0; ivtx < 4; ivtx++) {
      // _dcell_vtx[ivtx] = PDM_MAX(1, PDM_MIN(gn_vtx, _dcell_vtx[ivtx]));
      assert(_dcell_vtx[ivtx] > 0 && _dcell_vtx[ivtx] <= gn_vtx);
    }

  } // End of loop on cells

  free(base_cell_face_perm);


  free(hextet_k_idx[0]);
  free(hextet_k_idx[1]);
  free(hextet_k_idx[5]);




  /* Boundary faces */
  int *base_face_tag = PDM_array_zeros_int(base_n_face);
  for (int icell = 0; icell < base_n_cell; icell++) {
    for (int idx_face = 4*icell; idx_face < 4*(icell+1); idx_face++) {
      int face_id = PDM_ABS(base_cell_face[idx_face]) - 1;
      base_face_tag[face_id] += PDM_SIGN(base_cell_face[idx_face]);
    }
  }


  int base_n_bdr_face = 0;
  int *base_bdr_face = malloc(sizeof(int) * base_n_face);
  for (int iface = 0; iface < base_n_face; iface++) {
    if (base_face_tag[iface] != 0) {
      base_bdr_face[base_n_bdr_face++] = PDM_SIGN(base_face_tag[iface]) * (iface+1);
    }
  }
  base_bdr_face = realloc(base_bdr_face, sizeof(int) * base_n_bdr_face);


  PDM_g_num_t face_subface_n = (n+1)*(n+1);
  PDM_g_num_t gn_face = base_n_bdr_face*face_subface_n;


  int quadtria_size[2] = {
    n+1, n,
  };

  int subface_quadtria_n[2] = {
    _n_subtria(n),
    _n_subtria(n-1)
  };

  int subface_quadtria_idx[3];
  PDM_array_idx_from_sizes_int(subface_quadtria_n,
                               2,
                               subface_quadtria_idx);


  *distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);

  int dn_face = (int) ((*distrib_face)[i_rank+1] - (*distrib_face)[i_rank]);

  *dface_vtx = malloc(sizeof(PDM_g_num_t) * dn_face * 3);

  // //-->> tmp
  for (int i = 0; i < 3*dn_face; i++) {
    (*dface_vtx)[i] = 1;
  }
  // //<<--


  for (int iface = 0; iface < dn_face; iface++) {

    PDM_g_num_t g = (*distrib_face)[i_rank] + iface;

    PDM_g_num_t *_dface_vtx = *dface_vtx + 3*iface;

    int ibase = (int) (g / face_subface_n);
    int idx   = (int) (g - face_subface_n*ibase);


    int quadtria = PDM_binary_search_gap_int(idx,
                                             subface_quadtria_idx,
                                             3);

    idx -= subface_quadtria_idx[quadtria];

    int j, i;
    idx2ij(idx, quadtria_size[quadtria]-1, &i, &j);

    // log_trace("g = "PDM_FMT_G_NUM", ibase = %d, quadtria = %d, idx = %d, ij = %d %d, \n",
    //           g, ibase, quadtria, idx, i, j);

    int ibase_face = base_bdr_face[ibase];
    int ibase_sign = PDM_SIGN(ibase_face);
    ibase_face = PDM_ABS(ibase_face) - 1;
    // log_trace("  base_id = %d, base_sign = %d\n", ibase_face, ibase_sign);

    int fi;
    int fj = j;
    if (ibase_sign > 0) {
      fi = i;
    }
    else {
      fi = n-i;
    };

    switch (quadtria) {

      case 0: {
        // quadtria_vtx : [0, 1, 2]
        _quadtria_vtx_0(fi, fj, _dface_vtx[0]);
        _quadtria_vtx_1(fi, fj, _dface_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dface_vtx[2]);
        break;
      }

      case 1: {
        // quadtria_vtx : [1, 3, 2]
        _quadtria_vtx_1(fi, fj, _dface_vtx[0]);
        _quadtria_vtx_3(fi, fj, _dface_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dface_vtx[2]);
        break;
      }

      default: {
        PDM_error(__FILE__, __LINE__, 0, "Wrong subface quadtria %d\n", quadtria);
      }
    }

    // PDM_log_trace_array_long(_dface_vtx, 3, "_dface_vtx : ");
    for (int ivtx = 0; ivtx < 3; ivtx++) {
      // _dface_vtx[ivtx] = PDM_MAX(1, PDM_MIN(gn_vtx, _dface_vtx[ivtx]));
      assert(_dface_vtx[ivtx] > 0 && _dface_vtx[ivtx] <= gn_vtx);
    }

  } // End of loop on faces


  free(base_bdr_face);
  free(base_face_tag);

}


static void
_extrude_base_surface_mesh
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t    n,
 const PDM_g_num_t    n_layer,
 const int            base_n_vtx,
 const int            base_n_edge,
 const int            base_n_face,
       double        *base_vtx_coord,
       int           *base_edge_vtx,
       int           *base_face_vtx,
       int           *base_face_edge,
       double       **dvtx_coord,
       PDM_g_num_t  **dface_vtx,
       PDM_g_num_t  **dcell_vtx,
       int           *n_group,
       int          **dgroup_face_idx,
       PDM_g_num_t  **dgroup_face,
       PDM_g_num_t  **distrib_vtx,
       PDM_g_num_t  **distrib_face,
       PDM_g_num_t  **distrib_cell
)
{
  assert(n_layer > 0);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   *  Compute unit normals at base vertices
   */
  double *base_vtx_normal = malloc(sizeof(double) * base_n_vtx * 3);
  for (int i = 0; i < 3*base_n_vtx; i++) {
    base_vtx_normal[i] = 0.;
  }

  for (int iface = 0; iface < base_n_face; iface++) {
    double face_normal[3];

    PDM_geom_elem_tria_surface_vector(1,
                                      base_face_vtx + 3*iface,
                                      base_vtx_coord,
                                      face_normal,
                                      NULL,
                                      NULL);

    for (int ivtx = 0; ivtx < 3; ivtx++) {
      int vtx_id = base_face_vtx[3*iface+ivtx]-1;
      for (int i = 0; i < 3; i++) {
        base_vtx_normal[3*vtx_id+i] += face_normal[i];
      }
    }
  }

  /* Normalize */
  for (int i = 0; i < base_n_vtx; i++) {
    double mag = PDM_MODULE(base_vtx_normal + 3*i);
    assert(mag >= 0.);

    mag = 1./mag;
    for (int j = 0; j < 3; j++) {
      base_vtx_normal[3*i+j] *= mag;
    }
  }

  /* Auxiliary arrays */
  int *tria_j_idx = malloc(sizeof(int) * n);
  if (n > 0) {
    tria_j_idx[0] = 0;
    for (int j = 0; j < n-1; j++) {
      tria_j_idx[j+1] = tria_j_idx[j] + n-1-j;
    }
  }


  /*
   *  Vertices
   */
  PDM_g_num_t face_int_vtx_n = n*(n-1) / 2;

  PDM_g_num_t gn_vtx_layer =
  base_n_vtx                 +
  base_n_edge*n              +
  base_n_face*face_int_vtx_n;

  PDM_g_num_t gn_vtx = (n_layer+1)*gn_vtx_layer;

  *distrib_vtx = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  int dn_vtx = (int) ((*distrib_vtx)[i_rank+1] - (*distrib_vtx)[i_rank]);

  *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);

  PDM_g_num_t idx_vtx_edge = base_n_vtx;
  PDM_g_num_t idx_vtx_face = idx_vtx_edge + base_n_edge*n;

  double step_tangent = 1. / (double) (n + 1);
  double step_normal  = 1. / (double) n_layer;

  for (int ivtx = 0; ivtx < dn_vtx; ivtx++) {

    PDM_g_num_t g = (*distrib_vtx)[i_rank] + ivtx;

    PDM_g_num_t ilayer = g/gn_vtx_layer;
    g = g%gn_vtx_layer;

    if (g < idx_vtx_edge) {
      // Base vertex
      for (int i = 0; i < 3; i++) {
        (*dvtx_coord)[3*ivtx+i] = base_vtx_coord[3*g+i] + step_normal*ilayer*base_vtx_normal[3*g+i];
      }
    }

    else if (g < idx_vtx_face) {
      // Base edge
      int ibase = (int) ((g - idx_vtx_edge) / n);
      int i     = (int) (g - idx_vtx_edge - n*ibase);
      // log_trace("base edge %d, i = %d --> gnum "PDM_FMT_G_NUM"\n",
      //           ibase, i, g+1);

      int ivtx1 = base_edge_vtx[2*ibase  ] - 1;
      int ivtx2 = base_edge_vtx[2*ibase+1] - 1;
      double t = (i + 1) * step_tangent;

      for (int k = 0; k < 3; k++) {
        (*dvtx_coord)[3*ivtx + k] = (1 - t)*base_vtx_coord[3*ivtx1 + k] + t*base_vtx_coord[3*ivtx2 + k] +
        step_normal*ilayer * ((1 - t)*base_vtx_normal[3*ivtx1 + k] + t*base_vtx_normal[3*ivtx2 + k]);
      }
    }

    else {
      // Base face
      int ibase = (int) ((g - idx_vtx_face) / face_int_vtx_n);
      int idx   = (int) ( g - idx_vtx_face  - face_int_vtx_n*ibase);
      int j = PDM_binary_search_gap_int(idx,
                                        tria_j_idx,
                                        n);
      int i = idx - tria_j_idx[j];
      // log_trace("base face %d, ij = %d %d --> gnum "PDM_FMT_G_NUM"\n",
      //           ibase, i, j, g+1);

      int ivtx1 = base_face_vtx[3*ibase  ] - 1;
      int ivtx2 = base_face_vtx[3*ibase+1] - 1;
      int ivtx3 = base_face_vtx[3*ibase+2] - 1;
      double u = (i + 1) * step_tangent;
      double v = (j + 1) * step_tangent;

      for (int k = 0; k < 3; k++) {
        (*dvtx_coord)[3*ivtx + k] =
        (1 - u - v) * base_vtx_coord[3*ivtx1 + k] +
        u           * base_vtx_coord[3*ivtx2 + k] +
        v           * base_vtx_coord[3*ivtx3 + k] +
        step_normal*ilayer * ((1 - u - v) * base_vtx_normal[3*ivtx1 + k] +
                              u           * base_vtx_normal[3*ivtx2 + k] +
                              v           * base_vtx_normal[3*ivtx3 + k]);
      }
    }

  }
  free(tria_j_idx);
  free(base_vtx_normal);


  /*
   *  Cells (prisms)
   */
  int quadtria_size[2] = {
    n+1, n,
  };

  int subface_quadtria_n[2] = {
    _n_subtria(n),
    _n_subtria(n-1)
  };

  int subface_quadtria_idx[3];
  PDM_array_idx_from_sizes_int(subface_quadtria_n,
                               2,
                               subface_quadtria_idx);


  PDM_g_num_t face_subface_n = (n+1)*(n+1);
  PDM_g_num_t gn_cell_layer = face_subface_n * base_n_face;
  PDM_g_num_t gn_cell = n_layer * gn_cell_layer;

  *distrib_cell = PDM_compute_uniform_entity_distribution(comm, gn_cell);

  int dn_cell = (int) ((*distrib_cell)[i_rank+1] - (*distrib_cell)[i_rank]);

  *dcell_vtx = malloc(sizeof(PDM_g_num_t) * dn_cell * 6);

  for (int icell = 0; icell < dn_cell; icell++) {

    PDM_g_num_t *_dcell_vtx = *dcell_vtx + 6*icell;

    PDM_g_num_t g = (*distrib_cell)[i_rank] + icell;

    PDM_g_num_t ilayer = g/gn_cell_layer;
    g = g%gn_cell_layer;

    int ibase_face = (int) (g/face_subface_n);
    int idx        = (int) (g%face_subface_n);

    int quadtria = PDM_binary_search_gap_int(idx,
                                             subface_quadtria_idx,
                                             3);

    idx -= subface_quadtria_idx[quadtria];

    int j, i;
    idx2ij(idx, quadtria_size[quadtria]-1, &i, &j);

    // log_trace("icell %d ("PDM_FMT_G_NUM"), ilayer %d, ibase_face %d, idx %d, quadtria %d, i %d, j %d\n",
    //           icell, g, ilayer, ibase_face, idx, idx, quadtria, i, j);

    int fi = i;
    int fj = j;
    switch (quadtria) {

      case 0: {
        // quadtria_vtx : [0, 1, 2]
        _quadtria_vtx_0(fi, fj, _dcell_vtx[0]);
        _quadtria_vtx_1(fi, fj, _dcell_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dcell_vtx[2]);
        break;
      }

      case 1: {
        // quadtria_vtx : [1, 3, 2]
        _quadtria_vtx_1(fi, fj, _dcell_vtx[0]);
        _quadtria_vtx_3(fi, fj, _dcell_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dcell_vtx[2]);
        break;
      }

      default: {
        PDM_error(__FILE__, __LINE__, 0, "Wrong subface quadtria %d\n", quadtria);
      }
    }

    for (int k = 0; k < 3; k++) {
      _dcell_vtx[k] += gn_vtx_layer*ilayer;
      _dcell_vtx[3+k] = _dcell_vtx[k] + gn_vtx_layer;
    }
    // PDM_log_trace_array_long(_dcell_vtx, 6, "  _dcell_vtx : ");

  }


  /*
   *  Boundary faces
   *
   * /!\ if the base mesh has boundaries, we must generate quad faces as well
   */
  int *base_edge_tag = PDM_array_zeros_int(base_n_edge);
  for (int iface = 0; iface < base_n_face; iface++) {
    for (int idx_edge = 3*iface; idx_edge < 3*(iface+1); idx_edge++) {
      int edge_id = PDM_ABS(base_face_edge[idx_edge]) - 1;
      base_edge_tag[edge_id] += PDM_SIGN(base_face_edge[idx_edge]);
    }
  }


  int base_n_bdr_edge = 0;
  int *base_bdr_edge = malloc(sizeof(int) * base_n_edge);
  for (int iedge = 0; iedge < base_n_edge; iedge++) {
    if (base_edge_tag[iedge] != 0) {
      base_bdr_edge[base_n_bdr_edge++] = PDM_SIGN(base_edge_tag[iedge]) * (iedge+1);
    }
  }
  free(base_edge_tag);
  base_bdr_edge = realloc(base_bdr_edge, sizeof(int) * base_n_bdr_edge);

  assert(base_n_bdr_edge == 0);

  PDM_g_num_t gn_face = 2*gn_cell_layer;
  // + n_layer * base_n_brd_edge * n;

  *distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);

  int dn_face = (int) ((*distrib_face)[i_rank+1] - (*distrib_face)[i_rank]);

  *dface_vtx = malloc(sizeof(PDM_g_num_t) * dn_face * 3); // !!! if quads

  for (int iface = 0; iface < dn_face; iface++) {

    PDM_g_num_t *_dface_vtx = *dface_vtx + 3*iface;

    PDM_g_num_t g = (*distrib_face)[i_rank] + iface;

    int ilayer = (int) (g / gn_cell_layer);
    g = g % gn_cell_layer;

    int ibase_face = (int) (g / face_subface_n);
    int idx        = (int) (g - face_subface_n*ibase_face);


    int quadtria = PDM_binary_search_gap_int(idx,
                                             subface_quadtria_idx,
                                             3);

    idx -= subface_quadtria_idx[quadtria];

    int j, i;
    idx2ij(idx, quadtria_size[quadtria]-1, &i, &j);

    int fi = i;
    int fj = j;
    switch (quadtria) {

      case 0: {
        // quadtria_vtx : [0, 1, 2]
        _quadtria_vtx_0(fi, fj, _dface_vtx[0]);
        _quadtria_vtx_1(fi, fj, _dface_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dface_vtx[2]);
        break;
      }

      case 1: {
        // quadtria_vtx : [1, 3, 2]
        _quadtria_vtx_1(fi, fj, _dface_vtx[0]);
        _quadtria_vtx_3(fi, fj, _dface_vtx[1]);
        _quadtria_vtx_2(fi, fj, _dface_vtx[2]);
        break;
      }

      default: {
        PDM_error(__FILE__, __LINE__, 0, "Wrong subface quadtria %d\n", quadtria);
      }
    }

    if (ilayer == 0) {
      PDM_g_num_t tmp = _dface_vtx[0];
      _dface_vtx[0] = _dface_vtx[1];
      _dface_vtx[1] = tmp;
    }
    else {
      for (int k = 0; k < 3; k++) {
        _dface_vtx[k] += gn_vtx_layer*n_layer;
      }
    }

  }

  /* Groups */
  PDM_g_num_t *distrib_group = PDM_compute_uniform_entity_distribution(comm, gn_cell_layer);
  int dn_group_face = (int) (distrib_group[i_rank+1] - distrib_group[i_rank]);

  *n_group         = 2;
  *dgroup_face_idx = malloc(sizeof(int) * 3);
  (*dgroup_face_idx)[0] = 0;
  (*dgroup_face_idx)[1] = dn_group_face;
  (*dgroup_face_idx)[2] = 2*dn_group_face;
  *dgroup_face     = malloc(sizeof(PDM_g_num_t) * dn_group_face * 2);

  for (int i = 0; i < dn_group_face; i++) {
    (*dgroup_face)[i] = distrib_group[i_rank] + i + 1;
    (*dgroup_face)[dn_group_face+i] = gn_cell_layer + distrib_group[i_rank] + i + 1;
  }

  free(distrib_group);
}


static PDM_dmesh_nodal_t *
_set_dmesh_nodal
(
 const PDM_MPI_Comm          comm,
       double               *dvtx_coord,
       PDM_g_num_t          *dface_vtx,
       PDM_g_num_t          *dcell_vtx,
       PDM_g_num_t          *distrib_vtx,
       PDM_g_num_t          *distrib_face,
       PDM_g_num_t          *distrib_cell,
       int                   n_group,
       int                  *dgroup_face_idx,
       PDM_g_num_t          *dgroup_face,
       PDM_Mesh_nodal_elt_t  cell_type
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  PDM_g_num_t gn_cell = distrib_cell[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];
  int dn_cell = distrib_cell[i_rank+1] - distrib_cell[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  3,
                                                  gn_vtx,
                                                  gn_cell,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  /* Surface */
  dmn->surfacic->n_g_elmts = gn_face;
  int id_section_tria = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                          PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section_tria,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int          _n_group         = n_group;
  int         *_dgroup_face_idx = dgroup_face_idx;
  PDM_g_num_t *_dgroup_face     = dgroup_face;
  if (_dgroup_face == NULL) {
    _n_group = 1;
    _dgroup_face_idx = (int *) malloc(sizeof(int) * (_n_group + 1));
    _dgroup_face_idx[0] = 0;
    _dgroup_face_idx[1] = dn_face;

    _dgroup_face = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _dgroup_face_idx[_n_group]);
    for (int i = 0; i < dn_face; i++) {
      _dgroup_face[i] = distrib_face[i_rank] + i + 1;
    }
  }
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  _n_group,
                                  _dgroup_face_idx,
                                  _dgroup_face,
                                  PDM_OWNERSHIP_KEEP);

  /* Volume */
  dmn->volumic->n_g_elmts = gn_cell;
  int id_section_volume = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                            cell_type);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                        id_section_volume,
                                        dn_cell,
                                        dcell_vtx,
                                        PDM_OWNERSHIP_KEEP);

  // int n_group = 1;
  // int *dgroup_cell_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  // dgroup_cell_idx[0] = 0;
  // dgroup_cell_idx[1] = dn_cell;

  // PDM_g_num_t *dgroup_cell = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_cell_idx[n_group]);
  // for (int i = 0; i < dn_face; i++) {
  //   dgroup_cell[i] = distrib_cell[i_rank] + i + 1;
  // }
  // PDM_DMesh_nodal_elmts_group_set(dmn->volumic,
  //                                 n_group,
  //                                 dgroup_cell_idx,
  //                                 dgroup_cell,
  //                                 PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a volume mesh bounded by a sphere (deformed cube)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n_vtx_x         Number of vertices on segments in x-direction
 * \param[in]  n_vtx_y         Number of vertices on segments in y-direction
 * \param[in]  n_vtx_z         Number of vertices on segments in z-direction
 * \param[in]  radius          Radius of the sphere
 * \param[in]  center_x        x coordinate of the center of the sphere
 * \param[in]  center_y        y coordinate of the center of the sphere
 * \param[in]  center_z        z coordinate of the center of the sphere
 * \param[in]  t_elt           Element type
 * \param[in]  order           Element order
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_gen_nodal
(
 PDM_MPI_Comm           comm,
 const PDM_g_num_t      n_vtx_x,
 const PDM_g_num_t      n_vtx_y,
 const PDM_g_num_t      n_vtx_z,
 const double           radius,
 const double           center_x,
 const double           center_y,
 const double           center_z,
 PDM_Mesh_nodal_elt_t   t_elt,
 const int              order,
 PDM_dmesh_nodal_t    **dmn
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(t_elt);
  if (mesh_dimension != 3) {
    PDM_error(__FILE__, __LINE__, 0,
              "Not implemented yes for dimension %d\n", mesh_dimension);
  }

  /* First: generate a dcube nodal */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_x,
                                                        n_vtx_y,
                                                        n_vtx_z,
                                                        2*radius,
                                                        center_x - radius,
                                                        center_y - radius,
                                                        center_z - radius,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *_dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(_dmn);

  /* Second: "spherify" */
  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(_dmn);
  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(_dmn);

  double center[3] = {center_x, center_y, center_z};

  double *_dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  for (int i = 0; i < dn_vtx; i++) {

    // double r2 = 0., rinf = 0.;
    // for (int j = 0; j < 3; j++) {
    //   double x = dvtx_coord[3*i + j] - center[j];
    //   r2 += x*x;
    //   rinf = PDM_MAX(rinf, PDM_ABS(x));
    // }

    // r2 = sqrt(r2);

    // double scale = rinf/r2;
    // for (int j = 0; j < 3; j++) {
    //   double x = dvtx_coord[3*i + j] - center[j];

    //   _dvtx_coord[3*i + j] = center[j] + scale*x;
    // }
    double xyzc[3];
    for (int j = 0; j < 3; j++) {
      xyzc[j] = (dvtx_coord[3*i + j] - center[j]) / radius;
    }

    _cube_to_sphere(xyzc[0], xyzc[1], xyzc[2],
                    _dvtx_coord + 3*i,
                    _dvtx_coord + 3*i+1,
                    _dvtx_coord + 3*i+2);

    for (int j = 0; j < 3; j++) {
      _dvtx_coord[3*i + j] = center[j] + radius*_dvtx_coord[3*i + j];
    }

  }


  /* Third: get rid of ridges and unify surfaces */

  *dmn = PDM_DMesh_nodal_create(comm,
                                mesh_dimension,
                                distrib_vtx[n_rank],
                                1,
                                0,
                                0);

  PDM_DMesh_nodal_coord_set(*dmn,
                            dn_vtx,
                            _dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  // Volume
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  int n_section    = PDM_DMesh_nodal_n_section_get  (_dmn, geom_kind);
  int *sections_id = PDM_DMesh_nodal_sections_id_get(_dmn, geom_kind);

  PDM_g_num_t gn_elt_vol = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = sections_id[i_section];
    int _order;
    const char *ho_ordering = NULL;
    const PDM_g_num_t    *distrib_elt = PDM_DMesh_nodal_distrib_section_get(_dmn, geom_kind, id_section);
    int                   dn_elt      = PDM_DMesh_nodal_section_n_elt_get  (_dmn, geom_kind, id_section);
    PDM_g_num_t          *delt_vtx    = PDM_DMesh_nodal_section_std_ho_get (_dmn, geom_kind, id_section, &_order, &ho_ordering);
    PDM_Mesh_nodal_elt_t  elt_type    = PDM_DMesh_nodal_section_type_get   (_dmn, geom_kind, id_section);

    int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    PDM_g_num_t *_delt_vtx = malloc(sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);
    memcpy(_delt_vtx, delt_vtx, sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);

    int _id_section = PDM_DMesh_nodal_elmts_section_ho_add((*dmn)->volumic,
                                                           elt_type,
                                                           order,
                                                           ho_ordering);
    PDM_DMesh_nodal_elmts_section_std_set((*dmn)->volumic,
                                          _id_section,
                                          dn_elt,
                                          _delt_vtx,
                                          PDM_OWNERSHIP_KEEP);

    gn_elt_vol += distrib_elt[n_rank];
  }
  (*dmn)->volumic->n_g_elmts = gn_elt_vol;


  // Surface
  geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
  n_section   = PDM_DMesh_nodal_n_section_get  (_dmn, geom_kind);
  sections_id = PDM_DMesh_nodal_sections_id_get(_dmn, geom_kind);

  PDM_g_num_t gn_elt_surf = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = sections_id[i_section];
    int _order;
    const char *ho_ordering = NULL;
    const PDM_g_num_t    *distrib_elt = PDM_DMesh_nodal_distrib_section_get(_dmn, geom_kind, id_section);
    int                   dn_elt      = PDM_DMesh_nodal_section_n_elt_get  (_dmn, geom_kind, id_section);
    PDM_g_num_t          *delt_vtx    = PDM_DMesh_nodal_section_std_ho_get (_dmn, geom_kind, id_section, &_order, &ho_ordering);
    PDM_Mesh_nodal_elt_t  elt_type    = PDM_DMesh_nodal_section_type_get   (_dmn, geom_kind, id_section);

    int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    PDM_g_num_t *_delt_vtx = malloc(sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);
    memcpy(_delt_vtx, delt_vtx, sizeof(PDM_g_num_t) * dn_elt * elt_vtx_n);

    int _id_section = PDM_DMesh_nodal_elmts_section_ho_add((*dmn)->surfacic,
                                                           elt_type,
                                                           order,
                                                           ho_ordering);
    PDM_DMesh_nodal_elmts_section_std_set((*dmn)->surfacic,
                                          _id_section,
                                          dn_elt,
                                          _delt_vtx,
                                          PDM_OWNERSHIP_KEEP);

    gn_elt_surf += distrib_elt[n_rank];
  }
  (*dmn)->surfacic->n_g_elmts = gn_elt_surf;


  // Groups
  PDM_g_num_t *distrib_face = PDM_compute_uniform_entity_distribution(comm,
                                                                      gn_elt_surf);
  int n_group = 1;
  int *dgroup_elt_idx = malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  dgroup_elt_idx[n_group] = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  PDM_g_num_t *dgroup_elt = malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  for (int i = 0; i < dgroup_elt_idx[n_group]; i++) {
    dgroup_elt[i] = distrib_face[i_rank] + i + 1;
  }

  free(distrib_face);

  PDM_DMesh_nodal_elmts_group_set((*dmn)->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(*dmn);


  PDM_dcube_nodal_gen_free(dcube);
}






/**
 * \brief Create a volume mesh bounded by an icosphere
 * (all cells are tetrahedra)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dvtx_coord      Connectivity of distributed vertex to coordinates
 * \param[out] dface_vtx       Connectivity of distributed face to vertex
 * \param[out] dcell_vtx       Connectivity of distributed cell to vertex
 * \param[out] distrib_vtx     Distribution of vertices
 * \param[out] distrib_face    Distribution of faces
 * \param[out] distrib_face    Distribution of cells
 *
 */

void
PDM_sphere_vol_icosphere_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **dcell_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face,
       PDM_g_num_t       **distrib_cell
)
{
  assert(n >= 0);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  // Base vertices
  const int base_n_vtx = 13;
  const double base_vtx_coord[39] = {
    0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    0.0000000000000000,   0.0000000000000000,  -1.0000000000000000,
    0.7236073208526308,  -0.5257253406640403,  -0.4472195337774913,
    -0.2763880040858484,  -0.8506492161065184,  -0.4472198366964583,
    -0.8944261896382816,   0.0000000000000000,  -0.4472155982176212,
    -0.2763880040858484,   0.8506492161065184,  -0.4472198366964583,
    0.7236073208526308,   0.5257253406640403,  -0.4472195337774913,
    0.2763880040858484,  -0.8506492161065184,   0.4472198366964583,
    -0.7236073208526308,  -0.5257253406640403,   0.4472195337774913,
    -0.7236073208526308,   0.5257253406640403,   0.4472195337774913,
    0.2763880040858484,   0.8506492161065184,   0.4472198366964583,
    0.8944261896382816,   0.0000000000000000,   0.4472155982176212,
    0.0000000000000000,   0.0000000000000000,   1.0000000000000000
  };


  // Base cells (tetrahedra)
  const int base_n_cell = 20;
  const int base_cell_vtx[80] = {
    1, 3,  4,  2,
    1, 2,  7,  3,
    1, 4,  5,  2,
    1, 5,  6,  2,
    1, 6,  7,  2,
    1, 7,  12, 3,
    1, 3,  8,  4,
    1, 4,  9,  5,
    1, 5,  10, 6,
    1, 6,  11, 7,
    1, 12, 8,  3,
    1, 8,  9,  4,
    1, 9,  10, 5,
    1, 10, 11, 6,
    1, 11, 12, 7,
    1, 12, 13, 8,
    1, 8,  13, 9,
    1, 9,  13, 10,
    1, 10, 13, 11,
    1, 11, 13, 12
  };


  int  base_n_edge    = 0;
  int  base_n_face    = 0;
  int *base_edge_vtx  = NULL;
  int *base_face_vtx  = NULL;
  int *base_face_edge = NULL;
  int *base_cell_face = NULL;
  int *base_cell_edge = NULL;
  _build_base_edges_and_faces_from_tetra(base_n_vtx,
                                         base_n_cell,
                                         base_cell_vtx,
                                         &base_n_edge,
                                         &base_n_face,
                                         &base_edge_vtx,
                                         &base_face_vtx,
                                         &base_face_edge,
                                         &base_cell_face,
                                         &base_cell_edge);

  if (0 && i_rank == 0) {
    PDM_vtk_write_std_elements("icoball_base_edge.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               base_n_edge,
                               base_edge_vtx,
                               NULL,
                               0, NULL, NULL);

    PDM_vtk_write_std_elements("icoball_base_face.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_TRIA3,
                               base_n_face,
                               base_face_vtx,
                               NULL,
                               0, NULL, NULL);

    PDM_vtk_write_std_elements("icoball_base_cell.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_TETRA4,
                               base_n_cell,
                               base_cell_vtx,
                               NULL,
                               0, NULL, NULL);
  }



  _gen_from_base_mesh(comm,
                      n,
                      base_n_vtx,
                      base_n_edge,
                      base_n_face,
                      base_n_cell,
           (double *) base_vtx_coord,
                      base_edge_vtx,
                      base_face_vtx,
                      base_face_edge,
                      base_cell_face,
                      base_cell_edge,
              (int *) base_cell_vtx,
                      dvtx_coord,
                      dface_vtx,
                      dcell_vtx,
                      distrib_vtx,
                      distrib_face,
                      distrib_cell);


  /* "Spherify" */
  double eps = 1e-10;

  int dn_vtx = (*distrib_vtx)[i_rank+1] - (*distrib_vtx)[i_rank];

  int *ibase_cell = PDM_array_const_int(dn_vtx, -1);
  double *w = malloc(sizeof(double) * dn_vtx);

  for (int icell = 0; icell < base_n_cell; icell++) {

    double tetra_coord[12];
    for (int i = 0; i < 4; i++) {
      int vtx_id = base_cell_vtx[4*icell + i] - 1;
      memcpy(tetra_coord + 3*i,
             base_vtx_coord + 3*vtx_id,
             sizeof(double) * 3);
    }

    double v[3][3];
    for (int ivtx = 0; ivtx < 3; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        v[ivtx][idim] = tetra_coord[3*(ivtx+1) + idim] - tetra_coord[idim];
      }
    }

    double vol6 = v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
    v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
    v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]);

    double ivol6 = 1/vol6;

    double r[3][3];
    for (int i = 0; i < 3; i++) {
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;

      PDM_CROSS_PRODUCT(r[i], v[k], v[j]);

      for (int idim = 0; idim < 3; idim++) {
        r[i][idim] *= ivol6;
      }
    }

    for (int i = 0; i < dn_vtx; i++) {

      if (ibase_cell[i] >= 0) {
        continue;
      }

      double vp0[3] = {
        tetra_coord[0] - (*dvtx_coord)[3*i+0],
        tetra_coord[0] - (*dvtx_coord)[3*i+1],
        tetra_coord[0] - (*dvtx_coord)[3*i+2]
      };

      // Compute barycentrics
      double bc[4];
      bc[0] = 1.;
      for (int j = 0; j < 3; j++) {
        bc[j+1] = PDM_DOT_PRODUCT(vp0, r[j]);
        bc[0] -= bc[j+1];
      }

      if (bc[0] > -eps &&
          bc[1] > -eps &&
          bc[2] > -eps &&
          bc[3] > -eps) {
        ibase_cell[i] = icell;
        w[i] = PDM_MAX(bc[0], 0.);
      }

    } // End of loop on volume vtx

  } // End of loop on base cells

  for (int i = 0; i < dn_vtx; i++) {
    assert(ibase_cell[i] >= 0);

    double r = PDM_MODULE(*dvtx_coord + 3*i);
    double scale = radius * (1 - w[i]);
    if (r > 1.e-15) {
      scale /= r;
    }

    (*dvtx_coord)[3*i  ] = x_center + scale*(*dvtx_coord)[3*i  ];
    (*dvtx_coord)[3*i+1] = y_center + scale*(*dvtx_coord)[3*i+1];
    (*dvtx_coord)[3*i+2] = z_center + scale*(*dvtx_coord)[3*i+2];
  }
  free(ibase_cell);
  free(w);


  free(base_edge_vtx);
  free(base_face_vtx);
  free(base_face_edge);
  free(base_cell_face);
  free(base_cell_edge);
}




/**
 * \brief Create a volume mesh bounded by an icosphere
 * (all cells are tetrahedra)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_icosphere_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **dmn
)
{
  double      *dvtx_coord    = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *dcell_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;
  PDM_g_num_t *distrib_cell  = NULL;
  PDM_sphere_vol_icosphere_gen(comm,
                               n,
                               x_center,
                               y_center,
                               z_center,
                               radius,
                               &dvtx_coord,
                               &dface_vtx,
                               &dcell_vtx,
                               &distrib_vtx,
                               &distrib_face,
                               &distrib_cell);

  /*
   *  Create dmesh nodal
   */
  *dmn = _set_dmesh_nodal(comm,
                          dvtx_coord,
                          dface_vtx,
                          dcell_vtx,
                          distrib_vtx,
                          distrib_face,
                          distrib_cell,
                          0,
                          NULL,
                          NULL,
                          PDM_MESH_NODAL_TETRA4);

  free(distrib_vtx);
  free(distrib_face);
  free(distrib_cell);
}




/**
 * \brief Create a volume mesh bounded by two concentric icospheres
 * (all cells are prisms)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  n_layer         Number of extrusion layers
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius_interior Radius of the interior sphere
 * \param[in]  radius_exterior Radius of the exterior sphere
 * \param[in]  geometric_ratio Geometric ratio for layer thickness
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_hollow_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const PDM_g_num_t         n_layer,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius_interior,
 const double              radius_exterior,
 const double              geometric_ratio,
       PDM_dmesh_nodal_t **dmn
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(n_layer > 0);

  // Base vertices
  const int base_n_vtx = 12;
  const double base_vtx_coord[36] = {
    0.0000000000000000,   0.0000000000000000,  -1.0000000000000000,
    0.7236073208526308,  -0.5257253406640403,  -0.4472195337774913,
    -0.2763880040858484,  -0.8506492161065184,  -0.4472198366964583,
    -0.8944261896382816,   0.0000000000000000,  -0.4472155982176212,
    -0.2763880040858484,   0.8506492161065184,  -0.4472198366964583,
    0.7236073208526308,   0.5257253406640403,  -0.4472195337774913,
    0.2763880040858484,  -0.8506492161065184,   0.4472198366964583,
    -0.7236073208526308,  -0.5257253406640403,   0.4472195337774913,
    -0.7236073208526308,   0.5257253406640403,   0.4472195337774913,
    0.2763880040858484,   0.8506492161065184,   0.4472198366964583,
    0.8944261896382816,   0.0000000000000000,   0.4472155982176212,
    0.0000000000000000,   0.0000000000000000,   1.0000000000000000
  };

  // Base faces (triangles)
  const int base_n_face = 20;
  const int base_face_vtx[60] = {
    2,  3,  1,
    1,  6,  2,
    3,  4,  1,
    4,  5,  1,
    5,  6,  1,
    6,  11, 2,
    2,  7,  3,
    3,  8,  4,
    4,  9,  5,
    5,  10, 6,
    11, 7,  2,
    7,  8,  3,
    8,  9,  4,
    9,  10, 5,
    10, 11, 6,
    11, 12, 7,
    7,  12, 8,
    8,  12, 9,
    9,  12, 10,
    10, 12, 11
  };


  int  base_n_edge    = 0;
  int *base_edge_vtx  = NULL;
  int *base_face_edge = NULL;
  _build_edges_from_triangles(base_n_vtx,
                              base_n_face,
                              base_face_vtx,
                              &base_n_edge,
                              &base_edge_vtx,
                              &base_face_edge);

  if (0 && i_rank == 0) {
    PDM_vtk_write_std_elements("hollow_base_edge.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               base_n_edge,
                               base_edge_vtx,
                               NULL,
                               0, NULL, NULL);

    PDM_vtk_write_std_elements("hollow_base_face.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_TRIA3,
                               base_n_face,
                               base_face_vtx,
                               NULL,
                               0, NULL, NULL);
  }


  double      *dvtx_coord      = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  PDM_g_num_t *dcell_vtx       = NULL;
  int          n_group         = 0;
  int         *dgroup_face_idx = NULL;
  PDM_g_num_t *dgroup_face     = NULL;
  PDM_g_num_t *distrib_vtx     = NULL;
  PDM_g_num_t *distrib_face    = NULL;
  PDM_g_num_t *distrib_cell    = NULL;
  _extrude_base_surface_mesh(comm,
                             n,
                             n_layer,
                             base_n_vtx,
                             base_n_edge,
                             base_n_face,
                  (double *) base_vtx_coord,
                             base_edge_vtx,
                  (int    *) base_face_vtx,
                             base_face_edge,
                             &dvtx_coord,
                             &dface_vtx,
                             &dcell_vtx,
                             &n_group,
                             &dgroup_face_idx,
                             &dgroup_face,
                             &distrib_vtx,
                             &distrib_face,
                             &distrib_cell);
  free(base_edge_vtx);
  free(base_face_edge);



  /* "Spherify", scale and translate */
  double center[3] = {x_center, y_center, z_center};
  double delta_radius = radius_exterior - radius_interior;
  double step0;
  double idenom;

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (geometric_ratio == 1) {
    step0 = delta_radius / (double) n_layer;
  }
  else {
    idenom = 1./(1 - geometric_ratio);
    step0 = delta_radius * (1 - geometric_ratio) / (1 - pow(geometric_ratio, (double) n_layer));
  }

  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  PDM_g_num_t face_int_vtx_n = n*(n-1) / 2;

  PDM_g_num_t gn_vtx_layer =
  base_n_vtx                 +
  base_n_edge*n              +
  base_n_face*face_int_vtx_n;

  for (int i = 0; i < dn_vtx; i++) {
    PDM_g_num_t g = distrib_vtx[i_rank] + i;

    int ilayer = (int) (g / gn_vtx_layer);

    double rlayer = radius_interior;
    if (geometric_ratio == 1) {
      rlayer += step0*ilayer;
    }
    else {
     rlayer += step0*(1 - pow(geometric_ratio, ilayer)) * idenom;
    }
    double r = PDM_MODULE(dvtx_coord + 3*i);
    double scale = rlayer/r;

    for (int j = 0; j < 3; j++) {
      dvtx_coord[3*i+j] = center[j] + scale*dvtx_coord[3*i+j];
    }
  }
PDM_GCC_SUPPRESS_WARNING_POP


  /*
   *  Create dmesh nodal
   */
  *dmn = _set_dmesh_nodal(comm,
                          dvtx_coord,
                          dface_vtx,
                          dcell_vtx,
                          distrib_vtx,
                          distrib_face,
                          distrib_cell,
                          n_group,
                          dgroup_face_idx,
                          dgroup_face,
                          PDM_MESH_NODAL_PRISM6);

  free(distrib_vtx);
  free(distrib_face);
  free(distrib_cell);
}


#undef ij2idx
#undef ijk2idx
#undef _point_on_edge
#undef _point_on_face
#undef _hextet_vtx_0
#undef _hextet_vtx_1
#undef _hextet_vtx_2
#undef _hextet_vtx_3
#undef _hextet_vtx_4
#undef _hextet_vtx_5
#undef _hextet_vtx_6
#undef _hextet_vtx_7

#ifdef __cplusplus
}
#endif /* __cplusplus */
