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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_ho_location.h"
#include "pdm_mesh_nodal.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_ho_bezier.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#define ij2idx(i, j, n) ((i) + (j)*((n) + 1) - ((j)-1)*(j)/2)

/*============================================================================
 * Private function definitions
 *============================================================================*/

// TO DO: Newton in arbitrary dimension

static int
_newton_curve
(
 const int     order,
 const double *target,
       double *u,
       double *p,
       double *dp_du,
       double *xyz
)
{
  const int vb = 0;
  if (vb) {
    log_trace(">>> Newton curve\n");
  }

  const int    it_max  = 10;
  // const double tol_res = 1e-6;
  const double tol_u   = 1e-6;

  double _u = *u;

  double dxyz_du[3];

  int converged = 0;

  for (int it = 0; it < it_max; it++) {
    PDM_ho_bezier_de_casteljau_curve(3,
                                     order,
                                     _u,
                                     p,
                                     xyz,
                                     NULL, NULL);

    // log_trace("&&& %f %f %f\n", xyz[0], xyz[1], xyz[2]);

    double vec[3] = {
      target[0] - xyz[0],
      target[1] - xyz[1],
      target[2] - xyz[2]
    };

    // log_trace("dist2 = %f\n", PDM_DOT_PRODUCT(vec, vec));

    // Derivative
    PDM_ho_bezier_de_casteljau_curve(3,
                                     order-1,
                                     _u,
                                     dp_du,
                                     dxyz_du,
                                     NULL, NULL);

    double numer = PDM_DOT_PRODUCT(vec, dxyz_du);
    double denom = PDM_DOT_PRODUCT(dxyz_du, dxyz_du);

    if (vb) {
      log_trace("  it %d, u = %f, v = %f, res = %e\n",
                it, _u, 1-_u, PDM_ABS(numer));
    }

    // if (res < tol_res) {
    //   converged = 1;
    //   if (vb) {
    //     log_trace("  converged (res << 1)\n");
    //   }
    //   break;
    // }

    // TODO: check if Jacobian is singular

    double du = numer / denom;

    _u += du;

    if (_u >= 0 && _u <= 1) {
      if (PDM_ABS(du) < tol_u) {
        converged = 1; // ?
        if (vb) {
          log_trace("  converged (du << 1)\n");
        }
        break;
      }
    }
    else {
      // Outside curve
      converged = 1; // ?
      if (vb) {
        log_trace("  outside domain\n");
      }

      double dist[2] = {0.};
      double cp[3][2];
      for (int i = 0; i < 2; i++) {
        double *_p = p + i*3*order;

        for (int j = 0; j < 3; j++) {
          double dj = _p[j] - xyz[j];
          dist[i] += dj*dj;
        }
      }

      if (dist[0] < dist[1]) {
        // _u = 0.;
        memcpy(xyz, cp[0], sizeof(double)*3);
      } else {
        // _u = 1.;
        memcpy(xyz, cp[1], sizeof(double)*3);
      }

      break;
    }

  } // End of Newton iteration


  if (converged) {
    *u = _u;

    if (_u >= 0 && _u <= 1) {
      PDM_ho_bezier_de_casteljau_curve(3,
                                       order-1,
                                       _u,
                                       dp_du,
                                       dxyz_du,
                                       NULL, NULL);
    }

    // log_trace("u = %f\n", _u);
    // log_trace("&&& %f %f %f\n", xyz[0], xyz[1], xyz[2]);
  }

  return converged;
}


static int
_newton_triangle
(
 const int     order,
       double *target,
       double *u,
       double *v,
       double *p,
       double *dp_du,
       double *dp_dv,
       double *xyz
)
{
  const int vb = 0;
  if (vb) {
    log_trace(">>> Newton triangle\n");
  }

  const int    it_max  = 10;
  const double tol_res = 1e-6;
  const double tol_uv2 = 1e-12;

  double _u = *u;
  double _v = *v;

  double dxyz_du[3], dxyz_dv[3];

  int converged = 0;
  int outside   = 0;

  for (int it = 0; it < it_max; it++) {
    PDM_ho_bezier_de_casteljau_triangle(3,
                                        order,
                                        _u,
                                        _v,
                                        p,
                                        xyz,
                                        NULL, NULL, NULL);
    if (vb) {
      PDM_log_trace_array_double(xyz, 3, "xyz = ");
    }

    double vec[3] = {
      target[0] - xyz[0],
      target[1] - xyz[1],
      target[2] - xyz[2]
    };

    // log_trace("dist2 = %f\n", PDM_DOT_PRODUCT(vec, vec));

    // Jacobian
    PDM_ho_bezier_de_casteljau_triangle(3,
                                        order-1,
                                        _u,
                                        _v,
                                        dp_du,
                                        dxyz_du,
                                        NULL, NULL, NULL);

    PDM_ho_bezier_de_casteljau_triangle(3,
                                        order-1,
                                        _u,
                                        _v,
                                        dp_dv,
                                        dxyz_dv,
                                        NULL, NULL, NULL);

    double a00 = PDM_DOT_PRODUCT(dxyz_du, dxyz_du);
    double a01 = PDM_DOT_PRODUCT(dxyz_du, dxyz_dv);
    double a11 = PDM_DOT_PRODUCT(dxyz_dv, dxyz_dv);

    double b0 = PDM_DOT_PRODUCT(dxyz_du, vec);
    double b1 = PDM_DOT_PRODUCT(dxyz_dv, vec);

    double res = PDM_MAX(PDM_ABS(b0)/a00, PDM_ABS(b1)/a11);

    if (vb) {
      log_trace("  it %d, u = %f, v = %f, w = %f, res = %e\n",
                it, _u, _v, 1-_u-_v, res);
    }

    if (res < tol_res) {
      converged = 1;
      if (vb) {
        log_trace("  converged (res << 1)\n");
      }
      break;
    }

    double det = a00*a11 - a01*a01;

    // TODO: check if Jacobian is singular
    assert(PDM_ABS(det) > 1.e-16);

    double idet = 1. / det;

    double du = (b0*a11 - b1*a01)*idet;
    double dv = (b1*a00 - b0*a01)*idet;

    _u += du;
    _v += dv;
    if (vb) {
      log_trace("  du = %e, dv = %e (--> %f %f %f)\n",
                du, dv, _u, _v, 1-_u-_v);
    }

    if (_u >= 0 && _u <= 1 && _v >= 0 && _v <= 1 && _u + _v <= 1) {
      if (du*du + dv*dv < tol_uv2) {
        converged = 1; // ?
        if (vb) {
          log_trace("  converged (duv << 1)\n");
        }
        break;
      }
    }
    else {
      // Outside triangle
      converged = 1; // ?
      outside   = 1;
      if (vb) {
        log_trace("  outside domain\n");
      }
      break;
    }

  } // End of Newton iteration

  if (converged) {
    *u = _u;
    *v = _v;

    //-->>
    // TO DO: project on the triangle's curved edges
    if (outside) {
      if (1) {
        // Quick and dirty (projection inside uv domain)
        if (_u < 0) {
          _u = 0;
        }

        if (_v < 0) {
          _v = 0;
        }

        if (_u >= 0 && _v >= 0) {
          if (_v > _u + 1) {
            _u = 0;
            _v = 1;
          }
          else if (_v < _u - 1) {
            _u = 1;
            _v = 0;
          }
          else if (_u + _v > 1) {
            double c = _v - _u;
            _v = 0.5*(1 + c);
            _u = 1 - _v;
          }
        }
      }
      else {
        // Projection on the triangle's curved edges (in 3-space)
        double t[3];
        double cp[3][3];
        double dist[3];
        const int n_node_edge = order + 1;
        double edge_node_coord[3*n_node_edge];

        PDM_ho_bezier_curve_location(order,
                                     n_node_edge,
                                     p,
                                     target,
                                     cp[0],
                                     &t[0]);

        for (int j = 0; j <= order; j++) {
          int idx = ij2idx(order-j, j, order);
          memcpy(edge_node_coord + 3*j, p + 3*idx, sizeof(double) * 3);
        }
        PDM_ho_bezier_curve_location(order,
                                     n_node_edge,
                                     edge_node_coord,
                                     target,
                                     cp[1],
                                     &t[1]);

        for (int j = 0; j <= order; j++) {
          int idx = ij2idx(0, j, order);
          memcpy(edge_node_coord + 3*j, p + 3*idx, sizeof(double) * 3);
        }
        PDM_ho_bezier_curve_location(order,
                                     n_node_edge,
                                     edge_node_coord,
                                     target,
                                     cp[2],
                                     &t[2]);

        int jmin = -1;
        double dmin = HUGE_VAL;
        for (int j = 0; j < 3; j++) {
          if (dmin > dist[j]) {
            dmin = dist[j];
            jmin = j;
          }
        }

        // log_trace("closest edge: %d\n", jmin);
        memcpy(xyz, cp[jmin], sizeof(double) * 3);
        if (jmin == 0) {
          _u = t[0];
          _v = 0.;
        }
        else if (jmin == 1) {
          _u = 1 - t[1];
          _v = t[1];
        }
        else {
          _u = 0.;
          _v = t[2];
        }

      }
    }
    else {
      // inside
      PDM_ho_bezier_de_casteljau_triangle(3,
                                          order,
                                          _u,
                                          _v,
                                          p,
                                          xyz,
                                          NULL, NULL, NULL);
    }
    //<<--

    // log_trace("uv = %f %f\n", _u, _v);
    // log_trace("&&& %f %f %f\n", xyz[0], xyz[1], xyz[2]);
  }

  return converged;
}


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
 )
{
  const int n = order + 1;

  double s = 1. - t;

  double p0[n*dim];
  double p1[(n-1)*dim];
  double *p[2] = {p0, p1};

  // initialize
  if (b != NULL) {
    memcpy(p[0], b, sizeof(double) * n * dim);
  }
  else {
    int idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k < dim; k++) {
        if (k == idx) {
          p[0][dim*idx + k] = 1.;
        } else {
          p[0][dim*idx + k] = 0.;
        }
      }

      idx++;
    }
  }

  // subdivision
  if (l != NULL) {
    memcpy(l, b, sizeof(double)*dim);
  }
  if (r != NULL) {
    memcpy(r + (n-1)*dim, b + (n-1)*dim, sizeof(double)*dim);
  }

  for (int j = 1; j <= order; j++) {

    for (int i = 0; i <= order - j; i++) {
      for (int k = 0; k < dim; k++) {
        p[1][dim*i + k] = s*p[0][dim*i + k] + t*p[0][dim*(i+1) + k];
      }
    }

    // subdivision
    if (l != NULL) {
      memcpy(l + dim*j, p[1], sizeof(double)*dim);
    }
    if (r != NULL) {
      memcpy(r + dim*(order-j), p[1] + dim*(order-j), sizeof(double)*dim);
    }

    if (j < order) {
      // swap pointers
      double *tmp = p[1];
      p[1] = p[0];
      p[0] = tmp;
    }

  }

  if (val != NULL) {
    memcpy(val, p[1], sizeof(double) * dim);
  }
}


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
 )
{
  const int n = (order+1)*(order+2)/2;

  double w = 1. - u - v;

  double p0[n*dim];
  double p1[(n-order-1)*dim];
  double *p[2] = {p0, p1};

  // initialize
  if (b != NULL) {
    memcpy(p[0], b, sizeof(double) * n * dim);
  }
  else {
    int idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order - j; i++) {
        for (int k = 0; k < dim; k++) {
          if (k == idx) {
            p[0][dim*idx + k] = 1.;
          } else {
            p[0][dim*idx + k] = 0.;
          }
        }

        idx++;
      }
    }
  }


  // subdivision
  if (atr != NULL || ars != NULL || ast != NULL) {
    int idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        if (atr != NULL && j == 0) {
          int idx_atr = ij2idx(order-i, i, order);
          memcpy(atr + dim*idx_atr, b + dim*idx, sizeof(double)*dim);
        }
        if (ars != NULL) {
          memcpy(ars + dim*idx, b + dim*idx, sizeof(double)*dim);
        }
        if (ast != NULL && i == 0) {
          int idx_ast = ij2idx(j, order-j, order);
          memcpy(ast + dim*idx_ast, b + dim*idx, sizeof(double)*dim);
        }

        idx++;
      }
    }
  }



  for (int l = 1; l <= order; l++) {

    int idx = 0;
    for (int j = 0; j <= order-l; j++) {
      for (int i = 0; i <= order-l-j; i++) {
        int idxu = ij2idx(i+1, j,   order-l+1);
        int idxv = ij2idx(i,   j+1, order-l+1);
        int idxw = ij2idx(i,   j,   order-l+1);

        for (int k = 0; k < dim; k++) {
          p[1][dim*idx + k] =
          u*p[0][dim*idxu + k] +
          v*p[0][dim*idxv + k] +
          w*p[0][dim*idxw + k];
        }


        // subdivision
        if (atr != NULL && j == 0) {
          int idx_atr = ij2idx(order-l-i, i, order);
          memcpy(atr + dim*idx_atr, p[1] + dim*idx, sizeof(double)*dim);
        }
        if (ars != NULL) {
          int idx_ars = ij2idx(i, j, order);
          memcpy(ars + dim*idx_ars, p[1] + dim*idx, sizeof(double)*dim);
        }
        if (ast != NULL && i == 0) {
          int idx_ast = ij2idx(j, order-l-j, order);
          memcpy(ast + dim*idx_ast, p[1] + dim*idx, sizeof(double)*dim);
        }


        idx++;
      }
    }

    if (l < order) {
      // swap pointers
      double *tmp = p[1];
      p[1] = p[0];
      p[0] = tmp;
    }
  }

  if (val != NULL) {
    memcpy(val, p[1], sizeof(double) * dim);
  }
}


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
 )
{
  for (int i = 0; i < order; i++) {
    for (int k = 0; k < dim; k++) {
      db_dt[dim*i + k] = order * (b[dim*(i+1) + k] - b[dim*i + k]);
    }
  }
}


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
 const int     dim,
 const int     order,
       double *b,
       double *bu,
       double *bv
 )
{
  int idx = 0;
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order-j; i++) {

      int idxij  = ij2idx(i, j, order);
      int idxi1j = idxij + 1;            //ij2idx(i+1, j, order);
      int idxij1 = idxij + order - j + 1;//ij2idx(i, j+1, order);

      for (int k = 0; k < dim; k++) {
        bu[dim*idx + k] = order * (b[dim*idxi1j + k] - b[dim*idxij + k]);
        bv[dim*idx + k] = order * (b[dim*idxij1 + k] - b[dim*idxij + k]);
      }

      idx++;
    }
  }
}


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
 )
{
  PDM_UNUSED(n_node);

  double distance, t;
  distance = PDM_line_distance(point_coord,
                               node_coord,
                               node_coord + 3*order,
                               &t,
                               projected_coord);

  *u = t;

  if (order == 1) {
    return distance;
  }

  const int n = order*(order+1)/2;
  double db_du[3*n];
  PDM_ho_bezier_curve_derivative(3,
                                 order,
                                 node_coord,
                                 db_du);

  int converged = _newton_curve(order,
                                point_coord,
                                u,
                                node_coord,
                                db_du,
                                projected_coord);
  // log_trace("projected_coord : %f %f %f\n", projected_coord[0], projected_coord[1], projected_coord[2]);
  PDM_UNUSED(converged);

  distance = 0;
  for (int i = 0; i < 3; i++) {
    double d = point_coord[i] - projected_coord[i];
    distance += d*d;
  }

  return distance;
}


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
 )
{
  int vb = 0;

  double P1_coord[9];
  memcpy(P1_coord,     node_coord,                sizeof(double) * 3);
  memcpy(P1_coord + 3, node_coord + 3*order,      sizeof(double) * 3);
  memcpy(P1_coord + 6, node_coord + 3*(n_node-1), sizeof(double) * 3);

  double distance, weight[3];
  // PDM_triangle_evaluate_position(point_coord,
  //                                P1_coord,
  //                                projected_coord,
  //                                &distance,
  //                                weight);
  PDM_triangle_closest_point(point_coord,
                             P1_coord,
                             projected_coord,
                             &distance,
                             weight);
  uv[0] = weight[1];
  uv[1] = weight[2];
  uv[2] = weight[0];


  if (order == 1) {
    return distance;
  }

  const int n = order*(order+1)/2;
  double db_du[3*n], db_dv[3*n];
  PDM_ho_bezier_triangle_derivatives(3,
                                     order,
                                     node_coord,
                                     db_du,
                                     db_dv);

  int converged = _newton_triangle(order,
                                   point_coord,
                                   &uv[0],
                                   &uv[1],
                                   node_coord,
                                   db_du,
                                   db_dv,
                                   projected_coord);
  if (vb && !converged) {
    log_trace("_newton_triangle failed to convege :\n");
    PDM_log_trace_array_double(point_coord, 3, "point_coord : ");
    for (int i = 0; i < n_node; i++) {
      log_trace("node #%d : ", i);
      PDM_log_trace_array_double(node_coord+3*i, 3, "");
    }
  }

  if (vb) {
    log_trace("projected_coord : %f %f %f\n",
              projected_coord[0], projected_coord[1], projected_coord[2]);
  }
  PDM_UNUSED(converged);

  uv[2] = 1 - uv[0] - uv[1];

  distance = 0;
  for (int i = 0; i < 3; i++) {
    double d = point_coord[i] - projected_coord[i];
    distance += d*d;
  }

  return distance;
}





#undef ij2idx
