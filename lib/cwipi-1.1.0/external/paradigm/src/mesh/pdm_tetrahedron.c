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

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_triangle.h"
#include "pdm_geom_elem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_tetrahedron.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

/*=============================================================================
 * Private function definition
 *============================================================================*/
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
_solve_3x3(double  A[3][3],
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
    -   A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    +   A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
          - b[1]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
          + b[2]*(A[0][1]*A[1][2] - A[1][1]*A[0][2])) * det_inv;

  x1 = (  A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2])
          - A[1][0]*(b[0]*A[2][2] - b[2]*A[0][2])
          + A[2][0]*(b[0]*A[1][2] - b[1]*A[0][2])) * det_inv;

  x2 = (  A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1])
          - A[1][0]*(A[0][1]*b[2] - A[2][1]*b[0])
          + A[2][0]*(A[0][1]*b[1] - A[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return PDM_TRUE;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Evaluates the position in a tetrahedron
 *
 * \param [in]  x                         Point coordinates to evaluate position
 * \param [in]  vtx_coord                 Tetrahedron vertices coordinates
 * \param [out] closest_point             Closest Point in Tetrahedron or NULL
 * \param [out] closest_point_dist2       Square of the distance
 * \param [out] closest_point_weights     Vertices weights or NULL
 *
 * \return      -1 if the tetrahedron is degenerate, 0 else
 *
 */

PDM_tetrahedron_status_t
PDM_tetrahedron_evaluate_position
(
 const double  x[3],
 const double  vtx_coord[12],
 double        closest_point[3],
 double       *closest_point_dist2,
 double        closest_point_weights[4]
 )
{
  double vtx_tria[9];
  double p0p1[3], p0p2[3], p0p3[3], p1p2[3], p1p3[3], p2p3[3];
  double norm2_p0p1, norm2_p0p2, norm2_p0p3, norm2_p1p2, norm2_p1p3, norm2_p2p3;
  double xp0[3], xp1[3], xp2[3], xp3[3];
  double u, v, w;
  double vol6, ivol6;

  const double *pt0 = vtx_coord;
  const double *pt1 = vtx_coord + 3;
  const double *pt2 = vtx_coord + 6;
  const double *pt3 = vtx_coord + 9;


  for (int i = 0; i < 3; i++) {
    p0p1[i] = pt1[i] - pt0[i];
    p0p2[i] = pt2[i] - pt0[i];
    p0p3[i] = pt3[i] - pt0[i];

    p1p2[i] = pt2[i] - pt1[i];
    p1p3[i] = pt3[i] - pt1[i];
    p2p3[i] = pt3[i] - pt2[i];

    xp0[i] = pt0[i] - x[i];
    xp1[i] = pt1[i] - x[i];
    xp2[i] = pt2[i] - x[i];
    xp3[i] = pt3[i] - x[i];
  }

  norm2_p0p1 = PDM_DOT_PRODUCT (p0p1, p0p1);
  norm2_p0p2 = PDM_DOT_PRODUCT (p0p2, p0p2);
  norm2_p0p3 = PDM_DOT_PRODUCT (p0p3, p0p3);
  norm2_p1p2 = PDM_DOT_PRODUCT (p1p2, p1p2);
  norm2_p1p3 = PDM_DOT_PRODUCT (p1p3, p1p3);
  norm2_p2p3 = PDM_DOT_PRODUCT (p2p3, p2p3);

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (norm2_p0p1 == 0.0 ||
      norm2_p0p2 == 0.0 ||
      norm2_p0p3 == 0.0 ||
      norm2_p1p2 == 0.0 ||
      norm2_p1p3 == 0.0 ||
      norm2_p2p3 == 0.0) {
    return PDM_TETRAHEDRON_DEGENERATED;
  }

  vol6 = p0p1[0] * (p0p2[1]*p0p3[2] - p0p2[2]*p0p3[1])
    -    p0p1[1] * (p0p2[0]*p0p3[2] - p0p2[2]*p0p3[0])
    +    p0p1[2] * (p0p2[0]*p0p3[1] - p0p2[1]*p0p3[0]);

  if (vol6 == 0) {
    return PDM_TETRAHEDRON_DEGENERATED;
  }
PDM_GCC_SUPPRESS_WARNING_POP

  ivol6 = 1. / vol6;

  u = xp0[0] * (xp2[1]*xp3[2] - xp2[2]*xp3[1])
    - xp0[1] * (xp2[0]*xp3[2] - xp2[2]*xp3[0])
    + xp0[2] * (xp2[0]*xp3[1] - xp2[1]*xp3[0]);
  u *= -ivol6;

  v = xp0[0] * (xp3[1]*xp1[2] - xp3[2]*xp1[1])
    - xp0[1] * (xp3[0]*xp1[2] - xp3[2]*xp1[0])
    + xp0[2] * (xp3[0]*xp1[1] - xp3[1]*xp1[0]);
  v *= -ivol6;

  w = xp0[0] * (xp1[1]*xp2[2] - xp1[2]*xp2[1])
    - xp0[1] * (xp1[0]*xp2[2] - xp1[2]*xp2[0])
    + xp0[2] * (xp1[0]*xp2[1] - xp1[1]*xp2[0]);
  w *= -ivol6;


  if (u + v + w <= 1 && u >= 0 && v >= 0 && w >= 0) { // point a l'interieur du tetra
    if (closest_point != NULL) {
      for (int i = 0; i < 3; i++) {
        closest_point[i] = x[i];
      }
    }
    *closest_point_dist2 = 0.0;
    if (closest_point_weights != NULL) {
      closest_point_weights[0] = 1 - u - v - w;
      closest_point_weights[1] = u;
      closest_point_weights[2] = v;
      closest_point_weights[3] = w;
    }

    return PDM_TETRAHEDRON_INSIDE;
  }

  else if (u + v + w > 1) {// la face la plus proche est [P1,P2,P3]
    vtx_tria[0] = pt1[0];
    vtx_tria[1] = pt1[1];
    vtx_tria[2] = pt1[2];

    vtx_tria[3] = pt2[0];
    vtx_tria[4] = pt2[1];
    vtx_tria[5] = pt2[2];

    vtx_tria[6] = pt3[0];
    vtx_tria[7] = pt3[1];
    vtx_tria[8] = pt3[2];

    PDM_triangle_evaluate_position(x,
                                   vtx_tria,
                                   closest_point,
                                   closest_point_dist2,
                                   NULL);
    double p0cp[3], p1cp[3], p2cp[3], p3cp[3];

    for (int i = 0; i < 3; i++){
      p0cp[i] = closest_point[i] - pt0[i];
      p1cp[i] = closest_point[i] - pt1[i];
      p2cp[i] = closest_point[i] - pt2[i];
      p3cp[i] = closest_point[i] - pt3[i];
    }

    u = p0cp[0] * (p2cp[1]*p3cp[2] - p2cp[2]*p3cp[1])
      - p0cp[1] * (p2cp[0]*p3cp[2] - p2cp[2]*p3cp[0])
      + p0cp[2] * (p2cp[0]*p3cp[1] - p2cp[1]*p3cp[0]);
    u *= ivol6;

    v = p0cp[0] * (p3cp[1]*p1cp[2] - p3cp[2]*p1cp[1])
      - p0cp[1] * (p3cp[0]*p1cp[2] - p3cp[2]*p1cp[0])
      + p0cp[2] * (p3cp[0]*p1cp[1] - p3cp[1]*p1cp[0]);
    v *= ivol6;

    w = p0cp[0] * (p1cp[1]*p2cp[2] - p1cp[2]*p2cp[1])
      - p0cp[1] * (p1cp[0]*p2cp[2] - p1cp[2]*p2cp[0])
      + p0cp[2] * (p1cp[0]*p2cp[1] - p1cp[1]*p2cp[0]);
    w *= ivol6;

    closest_point_weights[0] = 1 - u - v - w;
    closest_point_weights[1] = u;
    closest_point_weights[2] = v;
    closest_point_weights[3] = w;
  }

  else if (u < 0) {// la face la plus proche est [P0,P3,P2]
    vtx_tria[0] = pt0[0];
    vtx_tria[1] = pt0[1];
    vtx_tria[2] = pt0[2];

    vtx_tria[3] = pt3[0];
    vtx_tria[4] = pt3[1];
    vtx_tria[5] = pt3[2];

    vtx_tria[6] = pt2[0];
    vtx_tria[7] = pt2[1];
    vtx_tria[8] = pt2[2];

    PDM_triangle_evaluate_position(x,
                                   vtx_tria,
                                   closest_point,
                                   closest_point_dist2,
                                   NULL);
    double p0cp[3], p1cp[3], p2cp[3], p3cp[3];
    for (int i = 0; i < 3; i++){
      p0cp[i] = closest_point[i] - pt0[i];
      p1cp[i] = closest_point[i] - pt1[i];
      p2cp[i] = closest_point[i] - pt2[i];
      p3cp[i] = closest_point[i] - pt3[i];
    }

    u = p0cp[0] * (p2cp[1]*p3cp[2] - p2cp[2]*p3cp[1])
      - p0cp[1] * (p2cp[0]*p3cp[2] - p2cp[2]*p3cp[0])
      + p0cp[2] * (p2cp[0]*p3cp[1] - p2cp[1]*p3cp[0]);
    u *= ivol6;

    v = p0cp[0] * (p3cp[1]*p1cp[2] - p3cp[2]*p1cp[1])
      - p0cp[1] * (p3cp[0]*p1cp[2] - p3cp[2]*p1cp[0])
      + p0cp[2] * (p3cp[0]*p1cp[1] - p3cp[1]*p1cp[0]);
    v *= ivol6;

    w = p0cp[0] * (p1cp[1]*p2cp[2] - p1cp[2]*p2cp[1])
      - p0cp[1] * (p1cp[0]*p2cp[2] - p1cp[2]*p2cp[0])
      + p0cp[2] * (p1cp[0]*p2cp[1] - p1cp[1]*p2cp[0]);
    w *= ivol6;

    closest_point_weights[0] = 1 - u - v - w;
    closest_point_weights[1] = u;
    closest_point_weights[2] = v;
    closest_point_weights[3] = w;
  }

  else if (v < 0) {// la face la plus proche est [P0,P3,P1]
    vtx_tria[0] = pt0[0];
    vtx_tria[1] = pt0[1];
    vtx_tria[2] = pt0[2];

    vtx_tria[3] = pt3[0];
    vtx_tria[4] = pt3[1];
    vtx_tria[5] = pt3[2];

    vtx_tria[6] = pt1[0];
    vtx_tria[7] = pt1[1];
    vtx_tria[8] = pt1[2];

    PDM_triangle_evaluate_position(x,
                                   vtx_tria,
                                   closest_point,
                                   closest_point_dist2,
                                   NULL);
    double p0cp[3], p1cp[3], p2cp[3], p3cp[3];
    for (int i = 0; i < 3; i++){
      p0cp[i] = closest_point[i] - pt0[i];
      p1cp[i] = closest_point[i] - pt1[i];
      p2cp[i] = closest_point[i] - pt2[i];
      p3cp[i] = closest_point[i] - pt3[i];
    }

    u = p0cp[0] * (p2cp[1]*p3cp[2] - p2cp[2]*p3cp[1])
      - p0cp[1] * (p2cp[0]*p3cp[2] - p2cp[2]*p3cp[0])
      + p0cp[2] * (p2cp[0]*p3cp[1] - p2cp[1]*p3cp[0]);
    u *= ivol6;

    v = p0cp[0] * (p3cp[1]*p1cp[2] - p3cp[2]*p1cp[1])
      - p0cp[1] * (p3cp[0]*p1cp[2] - p3cp[2]*p1cp[0])
      + p0cp[2] * (p3cp[0]*p1cp[1] - p3cp[1]*p1cp[0]);
    v *= ivol6;

    w = p0cp[0] * (p1cp[1]*p2cp[2] - p1cp[2]*p2cp[1])
      - p0cp[1] * (p1cp[0]*p2cp[2] - p1cp[2]*p2cp[0])
      + p0cp[2] * (p1cp[0]*p2cp[1] - p1cp[1]*p2cp[0]);
    w *= ivol6;

    closest_point_weights[0] = 1 - u - v - w;
    closest_point_weights[1] = u;
    closest_point_weights[2] = v;
    closest_point_weights[3] = w;
  }

  else if (w < 0) {// la face la plus proche est [P0,P1,P2]
    vtx_tria[0] = pt0[0];
    vtx_tria[1] = pt0[1];
    vtx_tria[2] = pt0[2];

    vtx_tria[3] = pt1[0];
    vtx_tria[4] = pt1[1];
    vtx_tria[5] = pt1[2];

    vtx_tria[6] = pt2[0];
    vtx_tria[7] = pt2[1];
    vtx_tria[8] = pt2[2];

    PDM_triangle_evaluate_position(x,
                                   vtx_tria,
                                   closest_point,
                                   closest_point_dist2,
                                   NULL);
    double p0cp[3], p1cp[3], p2cp[3], p3cp[3];
    for (int i = 0; i < 3; i++){
      p0cp[i] = closest_point[i] - pt0[i];
      p1cp[i] = closest_point[i] - pt1[i];
      p2cp[i] = closest_point[i] - pt2[i];
      p3cp[i] = closest_point[i] - pt3[i];
    }

    u = p0cp[0] * (p2cp[1]*p3cp[2] - p2cp[2]*p3cp[1])
      - p0cp[1] * (p2cp[0]*p3cp[2] - p2cp[2]*p3cp[0])
      + p0cp[2] * (p2cp[0]*p3cp[1] - p2cp[1]*p3cp[0]);
    u *= ivol6;

    v = p0cp[0] * (p3cp[1]*p1cp[2] - p3cp[2]*p1cp[1])
      - p0cp[1] * (p3cp[0]*p1cp[2] - p3cp[2]*p1cp[0])
      + p0cp[2] * (p3cp[0]*p1cp[1] - p3cp[1]*p1cp[0]);
    v *= ivol6;

    w = p0cp[0] * (p1cp[1]*p2cp[2] - p1cp[2]*p2cp[1])
      - p0cp[1] * (p1cp[0]*p2cp[2] - p1cp[2]*p2cp[0])
      + p0cp[2] * (p1cp[0]*p2cp[1] - p1cp[1]*p2cp[0]);
    w *= ivol6;

    closest_point_weights[0] = 1 - u - v - w;
    closest_point_weights[1] = u;
    closest_point_weights[2] = v;
    closest_point_weights[3] = w;
  }

  return PDM_TETRAHEDRON_OUTSIDE;
}




/**
 * \brief Computes tetrahedron barycenter
 *
 * \param [in]   pts     Tetrahedron vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_tetrahedron_compute_barycenter
(
 const double pts[12],
       double bary[3]
)
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 4; i++) {
    for (int ipt = 0; ipt < 3; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] *= 0.25;
  }
}



/**
 * \brief Computes the center and radius of a tetrahedron's circumsphere
 *
 * \param [in]   vtx_coord  Tetrahedron vertices coordinates
 * \param [out]  center     Circumsphere center
 * \param [out]  radius     Circumsphere radius
 *
 */

void
PDM_tetrahedron_circumsphere
(
 const double  vtx_coord[12],
 double        center[3],
 double       *radius
 )
{
  /*double d1[3] = {vtx_coord[ 3] - vtx_coord[0],
                  vtx_coord[ 4] - vtx_coord[1],
                  vtx_coord[ 5] - vtx_coord[2]};

  double d2[3] = {vtx_coord[ 6] - vtx_coord[0],
                  vtx_coord[ 7] - vtx_coord[1],
                  vtx_coord[ 8] - vtx_coord[2]};

  double d3[3] = {vtx_coord[ 9] - vtx_coord[0],
                  vtx_coord[10] - vtx_coord[1],
                  vtx_coord[11] - vtx_coord[2]};*/


  double A[3][3];
  double b[3] = {0.};
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      double aij = vtx_coord[3*(i+1) + j] - vtx_coord[j];
      A[i][j] = 2.*aij;
      b[i] += aij*aij;
    }
  }

  _solve_3x3 (A, b, center);

  double radius2 = 0.;
  for (int i = 0; i < 3; i++) {
    radius2 += center[i] * center[i];
    center[i] += vtx_coord[i];
  }

  *radius = sqrt(radius2);
}


/**
 * \brief Build tetrahedron->vertex from tetrahedron->face and face->vertex connectivities.
 *
 * \note In each tetrahedron, face #i is opposite to vertex #i.
 *
 * \param [in]  n_cell     Number of cells
 * \param [in]  cell_face  Cell -> face (signed) connectivity (1-based, size : 4 * \p n_cell)
 * \param [in]  cell_face  Face -> vertex connectivity (1-based, size : 3 * *n_face*)
 * \param [out] cell_vtx   Cell -> vertex (signed) connectivity (size : 4 * \p n_cell)
 *
 */

void
PDM_tetrahedron_ngon_to_nodal
(
 int   n_cell,
 int  *cell_face,
 int  *face_vtx,
 int **cell_vtx
 )
{
  *cell_vtx = malloc(sizeof(int) * n_cell * 4);

  for (int i = 0; i < n_cell; i++) {
    int iface = cell_face[4*i+3];

    if (iface < 0) {
      iface = -iface - 1;
      (*cell_vtx)[4*i  ] = face_vtx[3*iface  ];
      (*cell_vtx)[4*i+1] = face_vtx[3*iface+1];
      (*cell_vtx)[4*i+2] = face_vtx[3*iface+2];
    } else {
      iface = iface - 1;
      (*cell_vtx)[4*i  ] = face_vtx[3*iface+2];
      (*cell_vtx)[4*i+1] = face_vtx[3*iface+1];
      (*cell_vtx)[4*i+2] = face_vtx[3*iface  ];
    }

    iface = PDM_ABS(cell_face[4*i+1]) - 1;
    for (int j = 0; j < 3; j++) {
      int ivtx = face_vtx[3*iface+j];
      if (ivtx != (*cell_vtx)[4*i  ] &&
          ivtx != (*cell_vtx)[4*i+1] &&
          ivtx != (*cell_vtx)[4*i+2] ) {
        (*cell_vtx)[4*i+3] = ivtx;
      }
    }
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
