
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
#include "pdm_vtk.h"
#include "pdm_mesh_nodal.h"
#include "pdm_linear_algebra.h"
#include "pdm_lagrange_to_bezier.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Private function definitions
 *============================================================================*/


static inline double _pow(const double x, const int p) {
  if (p == 0) return 1.;

  double y = x;
  for (int i = 1; i < p; i++) {
    y *= x;
  }
  return y;
}

static void
_bezier_matrix_bar
(
 const int     order,
       double *b
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  // compute binomial coefficients
  int coef[order/2+1];
  coef[0] = 1;
  for (int n = 2; n <= order; n++) {

    if (n%2 == 0) coef[n/2] = coef[n/2-1];

    for (int k = n/2; k > 0; k--) {
      coef[k] += coef[k-1];
    }
  }

  double in = 1. / (double) order;

  b[0] = 1.;
  for (int j = 1; j <= order; j++) {
    b[j] = 0.;
  }

  for (int j = 0; j <= order; j++) {
    int c = coef[PDM_MIN(j,order-j)];
    for (int i = 1; i <= order/2; i++) {
      double u = i * in;
      b[n_nodes*i+j] = c * _pow(u,j) * _pow(1. - u, order-j);
    }
  }

  for (int i = order/2+1; i <= order; i++) {
    for (int j = 0; j <= order; j++) {
      b[n_nodes*i+j] = b[n_nodes*(order-i)+order-j];
    }
  }
}


static void
_bezier_matrix_tria
(
 const int     order,
       double *b
 )
{
#define ij2idx(i, j) ((i) + (j)*(order + 1 - (j)))
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  // compute trinomial coefficients
  const int n_coef = (order/2 + 1)*(order + 1 - order/2);
  int coef[n_coef];
  coef[0] = 1;
  for (int i = 1; i < n_coef; i++) {
    coef[i] = 0;
  }

  for (int n = 1; n <= order; n++) {
    for (int j = n/2; j >=0; j--) {
      int idx = ij2idx(n-j,j);
      for (int i = n-j; i >= j; i--) {
        if (i > 0) {
          if (i > j) {
            coef[idx] += coef[ij2idx(i-1,j)];
          } else {
            coef[idx] += coef[ij2idx(i,j-1)];
          }
        }

        if (j > 0) {
          coef[idx] += coef[ij2idx(i,j-1)];
        }
        idx--;
      }
    }
  }

  double in = 1. / (double) order;

  int icol = 0;
  for (int j = 0; j <= order; j++) {
    for (int i = 0; i <= order-j; i++) {

      int c = coef[ij2idx(PDM_MAX(i,j), PDM_MIN(i,j))];

      int irow = 0;
      for (int l = 0; l <= order; l++) {
        double v = l*in;
        for (int k = 0; k <= order-l; k++) {
          double u = k*in;
          b[n_nodes*irow+icol] = c * _pow(u,i) * _pow(v,j) * _pow(1 - u - v, order - i - j);
          irow++;
        }
      }
      icol++;
    }
  }

#undef ij2idx
}


static void
_bezier_matrix_quad
(
 const int     order,
       double *b
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

  double *a = malloc(sizeof(double) * (order+1)*(order+1));
  _bezier_matrix_bar(order, a);

  for (int i = 0; i <= order; i++) {
    for (int j = 0; j <= order; j++) {
      int k = i + (order+1)*j;
      for (int ii = 0; ii <= order; ii++) {
        for (int jj = 0; jj <= order; jj++) {
          int l = ii + (order+1)*jj;
          b[n_nodes*k+l] = a[(order+1)*i+ii] * a[(order+1)*j+jj];
        }
      }
    }
  }

  free(a);
}


/*----------------------------------------------------------------------------
 *  Get bounding boxes
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order bar
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_bar
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }
  }

  else if (order == 3) {

    double f833 = 5. / 6.;
    double f333 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f833*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f333*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f333*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f833*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }
  }

  else {

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes * n_nodes);
    }
    _bezier_matrix_bar(order, B);

    memcpy(bez, lag, sizeof(double) * n_nodes * 3);
    PDM_linear_algebra_linsolve_gauss(n_nodes, 3, B, bez);

    if (matrix == NULL) {
      free(B);
    }
  }

}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order triangle
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_tria
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -0.5*lag[6+j] + 2*lag[12+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = lag[15+j];
    }
  }

  else if (order == 3) {
    double f5_6 = 5. / 6.;
    double f1_3 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f5_6*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f1_3*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f1_3*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f5_6*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f5_6*lag[j] + 3*lag[12+j] - 1.5*lag[21+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f1_3*lag[j] - 0.75*lag[3+j] - 0.75*lag[6+j] + f1_3*lag[9+j] - 0.75*lag[12+j] + 4.5*lag[15+j] - 0.75*lag[18+j] - 0.75*lag[21+j] - 0.75*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f5_6*lag[9+j] + 3*lag[18+j] - 1.5*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = f1_3*lag[j] - 1.5*lag[12+j] + 3*lag[21+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f1_3*lag[9+j] - 1.5*lag[18+j] + 3*lag[24+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = lag[27+j];
    }
  }

  else {

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes * n_nodes);
    }
    _bezier_matrix_tria(order, B);

    memcpy(bez, lag, sizeof(double) * n_nodes * 3);
    PDM_linear_algebra_linsolve_gauss(n_nodes, 3, B, bez);

    if (matrix == NULL) {
      free(B);
    }
  }

}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order quadrangle
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_quad
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[18+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = 0.25*lag[j] - lag[3+j] + 0.25*lag[6+j] - lag[9+j] + 4*lag[12+j] - lag[15+j] + 0.25*lag[18+j] - lag[21+j] + 0.25*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = -0.5*lag[6+j] + 2*lag[15+j] - 0.5*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = lag[18+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = -0.5*lag[18+j] + 2*lag[21+j] - 0.5*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = lag[24+j];
    }
  }

  else if (order == 3) {

    double f833 = 5. / 6.;
    double f333 = 1. / 3.;
    double f694 = 25./ 36.;
    double f278 = 5. / 18.;
    double f111 = 1. / 9.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f833*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f333*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f333*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f833*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f833*lag[j] + 3*lag[12+j] - 1.5*lag[24+j] + f333*lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f694*lag[j] - 2.5*lag[3+j] + 1.25*lag[6+j] - f278*lag[9+j] - 2.5*lag[12+j] + 9*lag[15+j] - 4.5*lag[18+j] + lag[21+j] + 1.25*lag[24+j] - 4.5*lag[27+j] + 2.25*lag[30+j] - 0.5*lag[33+j] - f278*lag[36+j] + lag[39+j] - 0.5*lag[42+j] + f111*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f278*lag[j] + 1.25*lag[3+j] - 2.5*lag[6+j] + f694*lag[9+j] + lag[12+j] - 4.5*lag[15+j] + 9*lag[18+j] - 2.5*lag[21+j] - 0.5*lag[24+j] + 2.25*lag[27+j] - 4.5*lag[30+j] + 1.25*lag[33+j] + f111*lag[36+j] - 0.5*lag[39+j] + lag[42+j] - f278*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = -f833*lag[9+j] + 3*lag[21+j] - 1.5*lag[33+j] + f333*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f333*lag[j] -1.5*lag[12+j] + 3*lag[24+j] - f833*lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = -f278*lag[j] + lag[3+j] - 0.5*lag[6+j] + f111*lag[9+j] + 1.25*lag[12+j] - 4.5*lag[15+j] + 2.25*lag[18+j] - 0.5*lag[21+j] - 2.5*lag[24+j] + 9*lag[27+j] - 4.5*lag[30+j] + lag[33+j] + f694*lag[36+j] - 2.5*lag[39+j] + 1.25*lag[42+j] - f278*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[30+j] = f111*lag[j] - 0.5*lag[3+j] + lag[6+j] - f278*lag[9+j] - 0.5*lag[12+j] + 2.25*lag[15+j] - 4.5*lag[18+j] + 1.25*lag[21+j] + lag[24+j] - 4.5*lag[27+j] + 9*lag[30+j] - 2.5*lag[33+j] - f278*lag[36+j] + 1.25*lag[39+j] - 2.5*lag[42+j] + f694*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[33+j] = f333*lag[9+j] - 1.5*lag[21+j] + 3*lag[33+j] - f833*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[36+j] = lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[39+j] = -f833*lag[36+j] + 3*lag[39+j] - 1.5*lag[42+j] + f333*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[42+j] = f333*lag[36+j] - 1.5*lag[39+j] + 3*lag[42+j] - f833*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[45+j] = lag[45+j];
    }
  }

  else {

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes * n_nodes);
    }
    _bezier_matrix_quad(order, B);

    memcpy(bez, lag, sizeof(double) * n_nodes * 3);
    PDM_linear_algebra_linsolve_gauss(n_nodes, 3, B, bez);

    if (matrix == NULL) {
      free(B);
    }
  }
}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order tetrahedron
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_tetra
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order + 1 - (k)) - (j)*((j)-1)/2 + ((k)*((k)*((k) - 3*order - 6) + 3*order*(order + 4) + 11)) / 6)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {
    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

    double *work = malloc(sizeof(double) * n_nodes_tria * 3 * 4);

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes_tria * n_nodes_tria);
    }
    _bezier_matrix_tria(order, B);

    // hack for internal nodes
    memcpy(bez, lag, sizeof(double) * n_nodes * 3);


    int idx;

    // face w = 0
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order - j; i++) {
        for (int l = 0; l < 3; l++) {
          work[12*idx+l] = lag[3*idx+l];
        }
        idx++;
      }
    }


    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          work[3+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order - k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          work[6+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face u+v+w = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int i = order - j - k;
        int idx2 = ijk2idx(i,j,k);
        for (int l = 0; l < 3; l++) {
          work[9+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    PDM_linear_algebra_linsolve_gauss(n_nodes_tria, 12, B, work);

    // face w = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        for (int l = 0; l < 3; l++) {
          bez[3*idx+l] = work[12*idx+l];
        }
        idx++;
      }
    }

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[3+12*idx+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order - k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[6+12*idx+l];
        }
        idx++;
      }
    }

    // face u+v+w = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int i = order - j - k;
        int idx2 = ijk2idx(i,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[9+12*idx+l];
        }
        idx++;
      }
    }

    if (matrix == NULL) {
      free(B);
    }
    free(work);
  }
  #undef ijk2idx
}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order pyramid
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_pyramid
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order+1-(k)) + ((k)*((k)*(2*(k) - 6*order - 9) + 6*order*(order + 3) + 13)) / 6)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PYRAMIDHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {

    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *work = malloc(sizeof(double) * n_nodes_tria * 3 * 4);

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes_quad * n_nodes_quad);
    }

    // hack for internal nodes
    memcpy(bez, lag, sizeof(double) * n_nodes * 3);

    int idx;

    /* Quadrangular face */
    // face w = 0
    _bezier_matrix_quad(order, B);
    PDM_linear_algebra_linsolve_gauss(n_nodes_quad, 3, B, bez);


    /* Triangular face */
    _bezier_matrix_tria(order, B);

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          work[12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face u = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(order-k,j,k);
        for (int l = 0; l < 3; l++) {
          work[3+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          work[6+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,order-k,k);
        for (int l = 0; l < 3; l++) {
          work[9+12*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    PDM_linear_algebra_linsolve_gauss(n_nodes_tria, 12, B, work);

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[12*idx+l];
        }
        idx++;
      }
    }

    // face u = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(order-k,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[3+12*idx+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[6+12*idx+l];
        }
        idx++;
      }
    }

    // face v = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,order-k,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[9+12*idx+l];
        }
        idx++;
      }
    }

    if (matrix == NULL) {
      free(B);
    }
    free(work);
  }

  #undef ijk2idx
}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order prism
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_prism
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order+1) - (j)*((j)-1)/2 + (k)*(order+1)*(order+2)/2)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISMHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {

    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *work = malloc(sizeof(double) * n_nodes_quad * 3 * 3);

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes_quad * n_nodes_quad);
    }

    // hack for internal nodes
    memcpy(bez, lag, sizeof(double) * n_nodes * 3);


    int idx;


    /* Triangular faces */
    _bezier_matrix_tria(order, B);

    // face w = 0
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        for (int l = 0; l < 3; l++) {
          work[6*idx+l] = lag[3*idx+l];
        }
        idx++;
      }
    }

    // face w = 1
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        int idx2 = ijk2idx(i,j,order);
        for (int l = 0; l < 3; l++) {
          work[3+6*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    PDM_linear_algebra_linsolve_gauss(n_nodes_tria, 6, B, work);

    // face w = 0
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        for (int l = 0; l < 3; l++) {
          bez[3*idx+l] = work[6*idx+l];
        }
        idx++;
      }
    }

    // face w = 1
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        int idx2 = ijk2idx(i,j,order);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[3+6*idx+l];
        }
        idx++;
      }
    }


    /* Quadrangular faces */
    _bezier_matrix_quad(order, B);

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          work[9*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          work[3+9*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face u+v = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,order-i,k);
        for (int l = 0; l < 3; l++) {
          work[6+9*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    PDM_linear_algebra_linsolve_gauss(n_nodes_quad, 9, B, work);

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[9*idx+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[3+9*idx+l];
        }
        idx++;
      }
    }

    // face u+v = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,order-i,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[6+9*idx+l];
        }
        idx++;
      }
    }


    if (matrix == NULL) {
      free(B);
    }
    free(work);
  }

  #undef ijk2idx
}


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order hexahedron
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_hexa
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
)
{
  #define ijk2idx(i, j, k) ((i) + (order+1)*((j) + (order+1)*(k)))

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *work = malloc(sizeof(double) * n_nodes_quad * 3 * 6);

    double *B = matrix;
    if (matrix == NULL) {
      B = malloc(sizeof(double) * n_nodes_quad * n_nodes_quad);
    }
    _bezier_matrix_quad(order, B);

    // hack for internal nodes
    memcpy(bez, lag, sizeof(double) * n_nodes * 3);


    int idx;

    // face w = 0
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        for (int l = 0; l < 3; l++) {
          work[18*idx+l] = lag[3*idx+l];
        }
        idx++;
      }
    }

    // face w = 1
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,j,order);
        for (int l = 0; l < 3; l++) {
          // work[3+18*idx+l] = lag[order*n_nodes_quad + 3*idx+l];
          work[3+18*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          work[6+18*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face u = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(order,j,k);
        for (int l = 0; l < 3; l++) {
          work[9+18*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          work[12+18*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    // face v = 1
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,order,k);
        for (int l = 0; l < 3; l++) {
          work[15+18*idx+l] = lag[3*idx2+l];
        }
        idx++;
      }
    }

    PDM_linear_algebra_linsolve_gauss(n_nodes_quad, 18, B, work);

    // face w = 0
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        for (int l = 0; l < 3; l++) {
          bez[3*idx+l] = work[18*idx+l];
        }
        idx++;
      }
    }

    // face w = 1
    idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,j,order);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[3+18*idx+l];
          // bez[order*n_nodes_quad + 3*idx+l] = work[3+18*idx+l];
        }
        idx++;
      }
    }

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[6+18*idx+l];
        }
        idx++;
      }
    }

    // face u = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(order,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[9+18*idx+l];
        }
        idx++;
      }
    }

    // face v = 0
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[12+18*idx+l];
        }
        idx++;
      }
    }

    // face v = 1
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,order,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = work[15+18*idx+l];
        }
        idx++;
      }
    }

    if (matrix == NULL) {
      free(B);
    }
    free(work);
  }

}


/* Get bezier bounding box (copied from pdm_t_dcube_nodal_gen.c) */
void
PDM_bezier_bounding_boxes
(
 const PDM_Mesh_nodal_elt_t   t_elt,
 const int                    order,
 const int                    n_nodes,
 int                          n_elt,
 double                      *lagrange_coord,
 double                     **extents
)
{
  double *bezier_coord = malloc (sizeof(double) * n_nodes * 3);

  double *matrix = NULL;
  if (order > 3) {
    matrix = malloc(sizeof(double) * n_nodes * n_nodes);
  }

  for (int i = 0; i < n_elt; i++) {
    double *_min = (*extents) + 6*i;
    double *_max = _min + 3;

    for (int j = 0; j < 3; j++) {
      _min[j] =  1e30;
      _max[j] = -1e30;
    }

    if (t_elt == PDM_MESH_NODAL_BAR2 ||
        t_elt == PDM_MESH_NODAL_BARHO) {
      PDM_lagrange_to_bezier_bar(order, lagrange_coord, bezier_coord, matrix);
    }
    else if (t_elt == PDM_MESH_NODAL_TRIA3 ||
             t_elt == PDM_MESH_NODAL_TRIAHO) {
      PDM_lagrange_to_bezier_tria(order, lagrange_coord, bezier_coord, matrix);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Only implemented for other elements in pdm_t_dcube_nodal_gen.c\n");
    }

    for (int k = 0; k < n_nodes; k++) {
      for (int j = 0; j < 3; j++) {
        _min[j] = PDM_MIN(_min[j], bezier_coord[3*k + j]);
        _max[j] = PDM_MAX(_max[j], bezier_coord[3*k + j]);
      }
    }
  }

  free(bezier_coord);
  if (matrix != NULL) {
    free(matrix);
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
