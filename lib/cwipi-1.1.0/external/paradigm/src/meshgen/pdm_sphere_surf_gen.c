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

#include "pdm_sphere_surf_gen.h"
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
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

static inline void
idx2ij
(
 const int  idx,
 const int  n,
 int       *i,
 int       *j
 )
{
  // int _j = 0;
  // for (_j = 0; _j <= n; _j++) {
  //   if (idx >= ij2idx(0, _j, n) && idx < ij2idx(0, _j+1, n)) {
  //     break;
  //   }
  // }

  int b = -(2*n + 3);
  int d = b*b - 8*idx;
  int _j = (int) (0.5 * (-b - sqrt(d)));

  *i = idx - ij2idx(0, _j, n);
  *j = _j;
}


static PDM_dmesh_nodal_t *
_set_dmesh_nodal
(
 const PDM_MPI_Comm  comm,
       double       *dvtx_coord,
       PDM_g_num_t  *dface_vtx,
       PDM_g_num_t  *distrib_vtx,
       PDM_g_num_t  *distrib_face
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
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                     PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int n_group = 1;
  int *dgroup_elt_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  dgroup_elt_idx[1] = dn_face;

  PDM_g_num_t *dgroup_elt = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  for (int i = 0; i < dn_face; i++) {
    dgroup_elt[i] = distrib_face[i_rank] + i + 1;
  }
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a surface mesh of a sphere
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dvtx_coord      Connectivity of distributed vertex to coordinates
 * \param[out] dface_vtx_idx   Index of distributed face to vertex
 * \param[out] dface_vtx       Connectivity of distributed face to vertex
 * \param[out] distrib_vtx     Distribution of vertices
 * \param[out] distrib_face    Distribution of faces
 *
 */

void
PDM_sphere_surf_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(nu > 1);
  assert(nv > 1);

  /*
   *  Vertices
   */
  PDM_g_num_t gn_vtx = (nu - 1) * nv + 2;

  PDM_g_num_t *_distrib_vtx = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  const double step_u = 2*PDM_PI / (double) (nu - 1);
  const double step_v =   PDM_PI / (double) (nv + 1);

  int dn_vtx = (int) (_distrib_vtx[i_rank+1] - _distrib_vtx[i_rank]);
  double *_dvtx_coord = (double *) malloc(sizeof(double) * dn_vtx * 3);

  for (int i = 0; i < dn_vtx; i++) {

    PDM_g_num_t g = _distrib_vtx[i_rank] + i;

    if (g == gn_vtx-1) {
      // north pole
      _dvtx_coord[3*i  ] = x_center;
      _dvtx_coord[3*i+1] = y_center;
      _dvtx_coord[3*i+2] = z_center + radius;

    } else if (g == gn_vtx-2) {
      // south pole
      _dvtx_coord[3*i  ] = x_center;
      _dvtx_coord[3*i+1] = y_center;
      _dvtx_coord[3*i+2] = z_center - radius;

    } else {

      PDM_g_num_t jj = g / (nu - 1);
      PDM_g_num_t ii = g % (nu - 1);

      double v = -0.5*PDM_PI + (jj+1)*step_v;
      double u = ii*step_u;
      double c = cos(v);

      _dvtx_coord[3*i  ] = x_center + radius * c * cos(u);
      _dvtx_coord[3*i+1] = y_center + radius * c * sin(u);
      _dvtx_coord[3*i+2] = z_center + radius * sin(v);

    }

  }

  /*
   *  Faces
   */
  PDM_g_num_t gn_face = 2*((nu - 1)*(nv - 1) + nu - 1);

  PDM_g_num_t *_distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);

  int           dn_face       = (int) (_distrib_face[i_rank+1] - _distrib_face[i_rank]);
  int         *_dface_vtx_idx = (int *)         malloc(sizeof(int)         * (dn_face + 1));
  PDM_g_num_t *_dface_vtx     = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face * 3);

  _dface_vtx_idx[0] = 0;
  for (int i = 0; i < dn_face; i++) {

    PDM_g_num_t g = _distrib_face[i_rank] + i;

    _dface_vtx_idx[i+1] = _dface_vtx_idx[i] + 3;
    PDM_g_num_t *_fv = _dface_vtx + _dface_vtx_idx[i];


    if (g >= 2*(nu - 1)*(nv - 1) + nu - 1) {

      // north pole cap
      PDM_g_num_t ii = g - (2*(nu - 1)*(nv - 1) + nu - 1);

      _fv[0] = 1 + ii            + (nu-1)*(nv-1);
      _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(nv-1);
      _fv[2] = gn_vtx;
    }

    else if (g >= 2*(nu - 1)*(nv - 1)) {

      // south pole cap
      PDM_g_num_t ii = g - 2*(nu - 1)*(nv - 1);

      _fv[0] = 1 + ii;
      _fv[1] = gn_vtx - 1;
      _fv[2] = 1 + (ii+1)%(nu-1);
    }

    else {

      if (g >= (nu - 1)*(nv - 1)) {
        g -= (nu - 1)*(nv - 1);

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*(jj+1);
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);
      }

      else {

        PDM_g_num_t jj = g / (nu - 1);
        PDM_g_num_t ii = g % (nu - 1);

        _fv[0] = 1 + ii +            (nu-1)*jj;
        _fv[1] = 1 + (ii+1)%(nu-1) + (nu-1)*jj;
        _fv[2] = 1 + ii            + (nu-1)*(jj+1);

      }
    }
  }

  *distrib_vtx  = _distrib_vtx;
  *distrib_face = _distrib_face;

  *dvtx_coord    = _dvtx_coord;
  *dface_vtx_idx = _dface_vtx_idx;
  *dface_vtx     = _dface_vtx;


}


/**
 *
 * \brief Create a surface mesh of a sphere
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] _dmn            Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **_dmn
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;
  PDM_sphere_surf_gen(comm,
                      nu,
                      nv,
                      x_center,
                      y_center,
                      z_center,
                      radius,
                      &dvtx_coord,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &distrib_vtx,
                      &distrib_face);

  /*
   *  Create dmesh nodal
   */
  *_dmn = _set_dmesh_nodal(comm,
                           dvtx_coord,
                           dface_vtx,
                           distrib_vtx,
                           distrib_face);

  free(distrib_vtx);
  free(distrib_face);
  free(dface_vtx_idx);
}


/**
 *
 * \brief Create a surface mesh of a sphere (icosphere)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] dvtx_coord      Connectivity of distributed vertex to coordinates
 * \param[out] dface_vtx_idx   Index of distributed face to vertex
 * \param[out] dface_vtx       Connectivity of distributed face to vertex
 * \param[out] distrib_vtx     Distribution of vertices
 * \param[out] distrib_face    Distribution of faces
 *
 */

void
PDM_sphere_surf_icosphere_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
)
{

  assert(n >= 0);

  const int n_base_vtx  = 12;
  const int n_base_edge = 30;
  const int n_base_face = 20;

  double base_vtx_coord[12*3] = {
    0.00, 0.00, -1.00,
    0.72, -0.53, -0.45,
    -0.28, -0.85, -0.45,
    -0.89, 0.00, -0.45,
    -0.28, 0.85, -0.45,
    0.72, 0.53, -0.45,
    0.28, -0.85, 0.45,
    -0.72, -0.53, 0.45,
    -0.72, 0.53, 0.45,
    0.28, 0.85, 0.45,
    0.89, 0.00, 0.45,
    0.00, 0.00, 1.00
  };

  int base_edge_vtx[30*2] = {
    1, 2,
    2, 3,
    1, 3,
    1, 6,
    2, 6,
    3, 4,
    1, 4,
    4, 5,
    1, 5,
    5, 6,
    6, 11,
    2, 11,
    2, 7,
    3, 7,
    3, 8,
    4, 8,
    4, 9,
    5, 9,
    5, 10,
    6, 10,
    7, 11,
    7, 8,
    8, 9,
    9, 10,
    10, 11,
    11, 12,
    7, 12,
    8, 12,
    9, 12,
    10, 12
  };

  int base_face_vtx[20*3] = {
    1, 2, 3,
    2, 1, 6,
    1, 3, 4,
    1, 4, 5,
    1, 5, 6,
    2, 6, 11,
    3, 2, 7,
    4, 3, 8,
    5, 4, 9,
    6, 5, 10,
    2, 11, 7,
    3, 7, 8,
    4, 8, 9,
    5, 9, 10,
    6, 10, 11,
    7, 11, 12,
    8, 7, 12,
    9, 8, 12,
    10, 9, 12,
    11, 10, 12
  };

  int base_face_edge[20*3] = {
    1, 2, -3,
    -1, 4, -5,
    3, 6, -7,
    7, 8, -9,
    9, 10, -4,
    5, 11, -12,
    -2, 13, -14,
    -6, 15, -16,
    -8, 17, -18,
    -10, 19, -20,
    12, -21, -13,
    14, 22, -15,
    16, 23, -17,
    18, 24, -19,
    20, 25, -11,
    21, 26, -27,
    -22, 27, -28,
    -23, 28, -29,
    -24, 29, -30,
    -25, 30, -26
  };

  double center[3] = {x_center, y_center, z_center};



  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);


  PDM_g_num_t face_int_vtx_n = n*(n-1) / 2;
  PDM_g_num_t face_subface_n = (n+1)*(n+1);
  PDM_g_num_t gn_vtx  = n_base_vtx + n_base_edge*n + n_base_face*face_int_vtx_n;
  PDM_g_num_t gn_face = n_base_face*face_subface_n;

  *distrib_vtx  = PDM_compute_uniform_entity_distribution(comm, gn_vtx);
  *distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);

  int dn_vtx  = (int) ((*distrib_vtx)[i_rank+1]  - (*distrib_vtx)[i_rank]);
  int dn_face = (int) ((*distrib_face)[i_rank+1] - (*distrib_face)[i_rank]);

  *dvtx_coord    = malloc(sizeof(double     ) * dn_vtx  * 3);
  *dface_vtx     = malloc(sizeof(PDM_g_num_t) * dn_face * 3);
  *dface_vtx_idx = PDM_array_new_idx_from_const_stride_int(3, dn_face);

  PDM_g_num_t idx_edge = n_base_vtx;
  PDM_g_num_t idx_face = idx_edge + n_base_edge*n;

  double step = 1. / (double) (n + 1);

  /* Vertices */
  for (int ivtx = 0; ivtx < dn_vtx; ivtx++) {

    PDM_g_num_t g = (*distrib_vtx)[i_rank] + ivtx;

    if (g < idx_edge) {

      // Base vertex
      memcpy(*dvtx_coord + 3*ivtx, base_vtx_coord + 3*g, sizeof(double) * 3);

    }

    else if (g < idx_face) {

      // Base edge
      int ibase = (int) ((g - idx_edge) / n);
      int i     = (int) (g - idx_edge - n*ibase);

      int ivtx1 = base_edge_vtx[2*ibase  ] - 1;
      int ivtx2 = base_edge_vtx[2*ibase+1] - 1;
      double t = (i + 1) * step;

      for (int k = 0; k < 3; k++) {
        (*dvtx_coord)[3*ivtx + k] = (1 - t)*base_vtx_coord[3*ivtx1 + k] + t*base_vtx_coord[3*ivtx2 + k];
      }

    }

    else {

      // Base face
      int ibase = (int) ((g - idx_face) / face_int_vtx_n);
      int idx   = (int) (g - idx_face - face_int_vtx_n*ibase);
      int i, j;
      idx2ij(idx, n-2, &i, &j);

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

    // Translate and scale
    double r = PDM_MODULE(*dvtx_coord + 3*ivtx);

    double scale = radius / r;
    for (int k = 0; k < 3; k++) {
      (*dvtx_coord)[3*ivtx + k] = center[k] + scale*(*dvtx_coord)[3*ivtx + k];
    }
  }




  /* Faces */
  for (int iface = 0; iface < dn_face; iface++) {

    PDM_g_num_t g = (*distrib_face)[i_rank] + iface;

    int ibase = (int) (g / face_subface_n);
    int idx   = (int) (g - face_subface_n*ibase);

    if (idx < 3) {
      // Corner of base face
      if (n == 0) {
        for (int i = 0; i < 3; i++) {
          (*dface_vtx)[3*iface+i] = (PDM_g_num_t) base_face_vtx[3*ibase+i];
        }
      }

      else {
        (*dface_vtx)[3*iface] = base_face_vtx[3*ibase + idx];

        int iedge = base_face_edge[3*ibase + idx];
        if (iedge < 0) {
          iedge = -iedge - 1;
          (*dface_vtx)[3*iface+1] = idx_edge + n*iedge + n;
        } else {
          iedge = iedge - 1;
          (*dface_vtx)[3*iface+1] = idx_edge + n*iedge + 1;
        }

        iedge = base_face_edge[3*ibase + (idx+2)%3];
        if (iedge > 0) {
          iedge = iedge - 1;
          (*dface_vtx)[3*iface+2] = idx_edge + n*iedge + n;
        } else {
          iedge = -iedge - 1;
          (*dface_vtx)[3*iface+2] = idx_edge + n*iedge + 1;
        }
      }

    }

    else if (idx < 6) {
      // Corner of base face
      idx -= 3;

      if (n == 1) {

        for (int i = 0; i < 3; i++) {
          int iedge = PDM_ABS(base_face_edge[3*ibase + i]) - 1;
          (*dface_vtx)[3*iface+i] = 1 + idx_edge + n*iedge;
        }

      }

      else {

        if (idx == 0) {
          (*dface_vtx)[3*iface] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(0,  0,  n-2);
        } else if (idx == 1) {
          (*dface_vtx)[3*iface] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(n-2,0,  n-2);
        } else {
          (*dface_vtx)[3*iface] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(0,  n-2,n-2);
        }

        int iedge = base_face_edge[3*ibase + (idx+2)%3];
        if (iedge > 0) {
          iedge = iedge - 1;
          (*dface_vtx)[3*iface+1] = idx_edge + n*iedge + n;
        } else {
          iedge = -iedge - 1;
          (*dface_vtx)[3*iface+1] = idx_edge + n*iedge + 1;
        }

        iedge = base_face_edge[3*ibase + idx];
        if (iedge < 0) {
          iedge = -iedge - 1;
          (*dface_vtx)[3*iface+2] = idx_edge + n*iedge + n;
        } else {
          iedge = iedge - 1;
          (*dface_vtx)[3*iface+2] = idx_edge + n*iedge + 1;
        }
      }

    }

    else if (idx < 6 + 3*(n-1)) {

      int jedge = (int) ((idx - 6) / (n-1));
      int i     = (int) ((idx - 6) - jedge*(n-1));

      int iedge = base_face_edge[3*ibase + jedge];
      int sign, start;
      if (iedge < 0) {
        iedge = -iedge - 1;
        sign  = -1;
        start = n - 1;
      } else {
        iedge = iedge - 1;
        sign  = 1;
        start = 0;
      }

      (*dface_vtx)[3*iface  ] = 1 + idx_edge + n*iedge + start + sign*i;
      (*dface_vtx)[3*iface+1] = 1 + idx_edge + n*iedge + start + sign*(i+1);
      if (jedge == 0) {
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i, 0, n-2);
      } else if (jedge == 1) {
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(n-2-i, i, n-2);
      } else {
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(0, n-2-i, n-2);
      }

    }

    else if (idx < 6 + 3*(n-1) + 3*(n-2)) {

      int jedge = (int) ((idx - 6 - 3*(n-1)) / (n-2));
      int i     = (int) ((idx - 6 - 3*(n-1)) - jedge*(n-2));

      int iedge = base_face_edge[3*ibase + jedge];
      int sign, start;
      if (iedge < 0) {
        iedge = -iedge - 1;
        sign  = -1;
        start = n - 1;
      } else {
        iedge = iedge - 1;
        sign  = 1;
        start = 0;
      }

      (*dface_vtx)[3*iface  ] = 1 + idx_edge + n*iedge + start + sign*(i+1);
      if (jedge == 0) {
        (*dface_vtx)[3*iface+1] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i+1, 0, n-2);
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i,   0, n-2);
      } else if (jedge == 1) {
        (*dface_vtx)[3*iface+1] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(n-2-i-1, i+1, n-2);
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(n-2-i, i, n-2);
      } else {
        (*dface_vtx)[3*iface+1] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(0, n-2-i-1, n-2);
        (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(0, n-2-i, n-2);
      }

    }

    else if (idx < 6 + 3*(n-1) + 3*(n-2) + (n-2) + (n-3)*(n-2)/2) {

      idx -= (6 + 3*(n-1) + 3*(n-2));

      int i, j;
      idx2ij(idx, n-3, &i, &j);

      (*dface_vtx)[3*iface  ] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i,   j,   n-2);
      (*dface_vtx)[3*iface+1] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i+1, j,   n-2);
      (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i,   j+1, n-2);

    }

    else {

      idx -= (6 + 3*(n-1) + 3*(n-2) + (n-2) + (n-3)*(n-2)/2);

      int i, j;
      idx2ij(idx, n-4, &i, &j);

      (*dface_vtx)[3*iface  ] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i+1, j,   n-2);
      (*dface_vtx)[3*iface+1] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i+1, j+1, n-2);
      (*dface_vtx)[3*iface+2] = 1 + idx_face + face_int_vtx_n*ibase + ij2idx(i,   j+1, n-2);

    }



  }


}


/**
 *
 * \brief Create a surface mesh of a sphere (icosphere)
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Number of icosphere subdivisions
 * \param[in]  x_center        x coordinate of the center of the sphere
 * \param[in]  y_center        y coordinate of the center of the sphere
 * \param[in]  z_center        z coordinate of the center of the sphere
 * \param[in]  radius          Radius of the sphere
 * \param[out] _dmn            Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_icosphere_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **_dmn
)
{
  double      *dvtx_coord    = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;
  PDM_sphere_surf_icosphere_gen(comm,
                                n,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dvtx_coord,
                                &dface_vtx_idx,
                                &dface_vtx,
                                &distrib_vtx,
                                &distrib_face);

  /*
   *  Create dmesh nodal
   */
  *_dmn = _set_dmesh_nodal(comm,
                           dvtx_coord,
                           dface_vtx,
                           distrib_vtx,
                           distrib_face);

  free(distrib_vtx);
  free(distrib_face);
  free(dface_vtx_idx);
}






void
PDM_sphere_surf_icosphere_gen_part
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
       int               **pn_vtx,
       double           ***pvtx_coord,
       PDM_g_num_t      ***pvtx_ln_to_gn,
       int               **pn_face,
       int              ***pface_vtx_idx,
       int              ***pface_vtx,
       PDM_g_num_t      ***pface_ln_to_gn
)
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_surf_icosphere_gen_nodal(comm,
                                      n,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  PDM_multipart_compute(mpart);


  *pn_vtx         = malloc(sizeof(int          ) * n_part);
  *pvtx_coord     = malloc(sizeof(double      *) * n_part);
  *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_face        = malloc(sizeof(int          ) * n_part);
  *pface_vtx_idx  = malloc(sizeof(int         *) * n_part);
  *pface_vtx      = malloc(sizeof(int         *) * n_part);
  *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {

    /* Vertices */
    PDM_g_num_t *_vtx_ln_to_gn;
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
                                                       &_vtx_ln_to_gn,
                                                       PDM_OWNERSHIP_USER);
    (*pvtx_ln_to_gn)[ipart] = _vtx_ln_to_gn;
    // (*pvtx_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*pn_vtx)[ipart]);
    // memcpy((*pvtx_ln_to_gn)[ipart], _vtx_ln_to_gn, sizeof(PDM_g_num_t) * (*pn_vtx)[ipart]);

    double *_vtx_coord;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     ipart,
                                     &_vtx_coord,
                                     PDM_OWNERSHIP_USER);
    (*pvtx_coord)[ipart] = _vtx_coord;
    // (*pvtx_coord)[ipart] = malloc(sizeof(double) * (*pn_vtx)[ipart] * 3);
    // memcpy((*pvtx_coord)[ipart], _vtx_coord, sizeof(double) * (*pn_vtx)[ipart] * 3);

    /* Faces */
    PDM_g_num_t *_face_ln_to_gn;
    (*pn_face)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_FACE,
                                                       &_face_ln_to_gn,
                                                       PDM_OWNERSHIP_USER);
    (*pface_ln_to_gn)[ipart] = _face_ln_to_gn;
    // (*pface_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*pn_face)[ipart]);
    // memcpy((*pface_ln_to_gn)[ipart], _face_ln_to_gn, sizeof(PDM_g_num_t) * (*pn_face)[ipart]);

    int *_face_vtx;
    int *_face_vtx_idx;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &_face_vtx_idx,
                                        &_face_vtx,
                                        PDM_OWNERSHIP_USER);

    if (_face_vtx != NULL) {
      (*pface_vtx_idx)[ipart] = _face_vtx_idx;
      // (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      // memcpy((*pface_vtx_idx)[ipart], _face_vtx_idx, sizeof(int) * ((*pn_face)[ipart]+1));

      (*pface_vtx)[ipart] = _face_vtx;
      // (*pface_vtx)[ipart] = malloc(sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
      // memcpy((*pface_vtx)[ipart], _face_vtx, sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
    }

    else {
      int *_face_edge;
      int *_face_edge_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &_face_edge_idx,
                                          &_face_edge,
                                          PDM_OWNERSHIP_USER);

      int *_edge_vtx;
      int *_edge_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &_edge_vtx_idx,
                                          &_edge_vtx,
                                          PDM_OWNERSHIP_KEEP);

      // (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      // memcpy((*pface_vtx_idx)[ipart], _face_edge_idx, sizeof(int) * ((*pn_face)[ipart]+1));
      (*pface_vtx_idx)[ipart] = _face_edge_idx;

      PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &(*pface_vtx)[ipart]);
      free(_face_edge);
    }
  }

  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
}


#undef ij2idx

#ifdef __cplusplus
}
#endif /* __cplusplus */
