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
#include "pdm_mesh_nodal.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_line.h"

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

/**
 *
 * \brief Edge Bézier basis
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  u               Parametric coordinates (size = \ref n_pts)
 * \param [out] weights         Weights (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_edge
(
 const int              order,
 const int              n_pts,
 const double *restrict u,
 double       *restrict weights
)
{

  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      double _u =  u[i];

      weights[2*i]   = 1 - _u;
      weights[2*i+1] = _u;

    }
  }

  if (order == 2) {
    for (int i = 0; i < n_pts; i++) {

      double _u =  u[i];

      weights[3*i+0] = (1 - _u) * (1 - _u);
      weights[3*i+1] = 2 * _u * (1 - _u);
      weights[3*i+2] = _u * _u;

    }
  }

  if (order == 3) {
    for (int i = 0; i < n_pts; i++) {

      double _u =  u[i];

      weights[4*i+0] = (1 - _u) * (1 - _u) * (1 - _u);
      weights[4*i+1] = 3 * _u * (1 - _u) * (1 - _u);
      weights[4*i+2] = 3 * _u * _u * (1 - _u);
      weights[4*i+3] = _u * _u * _u;

    }
  }

  if (order > 3) {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }
}

/**
 *
 * \brief Triangle Bézier basis
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  uv              Parametric coordinates (size = \ref n_pts)
 * \param [out] weights         Weights (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_tria
(
 const int              order,
 const int              n_pts,
 const double *restrict uv,
 double       *restrict weights
)
{

  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      double u =  uv[2*i];
      double v =  uv[2*i+1];

      weights[3*i]   = 1 - u - v;
      weights[3*i+1] = u;
      weights[3*i+2] = v;

    }
  }

  if (order == 2) {
    for (int i = 0; i < n_pts; i++) {

      double u =  uv[2*i];
      double v =  uv[2*i+1];

      weights[6*i+0] = (1 - u -v) * (1 - u -v);  // (i,j,k)=(0,0,2)
      weights[6*i+1] = 2 * u * (1 - u -v);       // (i,j,k)=(1,0,1)
      weights[6*i+2] = u * u;                    // (i,j,k)=(2,0,0)
      weights[6*i+3] = 2 * v * (1 - u -v);       // (i,j,k)=(0,1,1)
      weights[6*i+4] = 2 * u * v;                // (i,j,k)=(1,1,0)
      weights[6*i+5] = v * v;                    // (i,j,k)=(0,2,0)

    }
  }

  if (order == 3) {

    for (int i = 0; i < n_pts; i++) {

      double u =  uv[2*i];
      double v =  uv[2*i+1];

      weights[10*i+0] = (1 - u -v) * (1 - u -v) * (1 - u -v);  // (i,j,k)=(003)
      weights[10*i+3] = u * u * u;                             // (i,j,k)=(300)
      weights[10*i+9] = v * v * v;                             // (i,j,k)=(030)
      weights[10*i+1] = 3 * u * (1 - u -v) * (1 - u -v);       // (i,j,k)=(102)
      weights[10*i+2] = 3 * u * u * (1 -u -v);                 // (i,j,k)=(201)
      weights[10*i+6] = 3 * u * u * v;                         // (i,j,k)=(210)
      weights[10*i+8] = 3 * u * v * v;                         // (i,j,k)=(120)
      weights[10*i+7] = 3 * v * v * (1 - u -v);                // (i,j,k)=(021)
      weights[10*i+4] = 3 * v * (1 - u -v) * (1 - u -v);       // (i,j,k)=(012)
      weights[10*i+5] = 6 * u * v * (1 - u -v);                // (i,j,k)=(111)

    }
  }

  if (order > 3) {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }
}

/**
 *
 * \brief Edge Bézier basis derivative with respect to u
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  u               Parametric coordinates (size = \ref n_pts)
 * \param [out] dw              Weights derivative (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_edge_derivative
(
 const int              order,
 const int              n_pts,
 const double *restrict u,
 double       *restrict dw
)
{

  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      dw[2*i]   = -1;
      dw[2*i+1] = 1;

    }
  }

  if (order == 2) {
    for (int i = 0; i < n_pts; i++) {

      double _u =  u[i];

      dw[3*i+0] = 2 * (_u - 1);
      dw[3*i+1] = 2 * (-2 * _u + 1);
      dw[3*i+2] = 2 * _u;
    }
  }

  if (order == 3) {
    for (int i = 0; i < n_pts; i++) {

      double _u =  u[i];

      dw[4*i+0] = 3 * (-1 + 2 * _u - _u * _u);
      dw[4*i+1] = 3 * (1 - 2 * _u + 3 * _u * _u);
      dw[4*i+2] = 3 * (2 * _u - 3 * _u * _u);
      dw[4*i+3] = 3 * _u * _u;
    }
  }

  if (order > 3) {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }
}

/**
 *
 * \brief Triangle Bézier basis derivative
 *
 * \param [in]  order           Order
 * \param [in]  n_pts           Number of points
 * \param [in]  uvw              Parametric coordinates (size = \ref n_pts)
 * \param [out] dw_du           Weights derivative with respect to u (size = n_nodes * \ref n_pts)
 * \param [out] dw_dv           Weights derivative with respect to v (size = n_nodes * \ref n_pts)
 * \param [out] dw_dw           Weights derivative with respect to w = 1 - u - v (size = n_nodes * \ref n_pts)
 *
 */

static void
_basis_bezier_tria_derivative
(
 const int              order,
 const int              n_pts,
 const double *restrict uvw,
 double       *restrict dw_du,
 double       *restrict dw_dv,
 double       *restrict dw_dw
)
{
  PDM_UNUSED(dw_dw);
  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      dw_du[3*i]   = -1;
      dw_du[3*i+1] = 1;
      dw_du[3*i+2] = 0;

      dw_dv[3*i]   = -1;
      dw_dv[3*i+1] = 0;
      dw_dv[3*i+2] = 1;

    }
  }

  if (order == 2) {
    for (int i = 0; i < n_pts; i++) {

      double u =  uvw[2*i];
      double v =  uvw[2*i+1];

      dw_du[6*i+0] = 2 * (-1 + u + v);    // (i,j,k)=(0,0,2)
      dw_du[6*i+1] = 2 * (1 - 2 * u - v); // (i,j,k)=(1,0,1)
      dw_du[6*i+2] = 2 * u;               // (i,j,k)=(2,0,0)
      dw_du[6*i+3] = 2 * (- v);           // (i,j,k)=(0,1,1)
      dw_du[6*i+4] = 2 * v;               // (i,j,k)=(1,1,0)
      dw_du[6*i+5] = 0;                   // (i,j,k)=(0,2,0)

      dw_dv[6*i+0] = 2 * (-1 + v + u);    // (i,j,k)=(0,0,2)
      dw_dv[6*i+1] = 2 * (- u);           // (i,j,k)=(1,0,1)
      dw_dv[6*i+2] = 0;                   // (i,j,k)=(2,0,0)
      dw_dv[6*i+3] = 2 * (1 - u - 2 * v); // (i,j,k)=(0,1,1)
      dw_dv[6*i+4] = 2 * u;               // (i,j,k)=(1,1,0)
      dw_dv[6*i+5] = 2 * v;               // (i,j,k)=(0,2,0)

    }
  }

  if (order == 3) {

    for (int i = 0; i < n_pts; i++) {

      double u =  uvw[2*i];
      double v =  uvw[2*i+1];

      dw_du[10*i+0] = 3 * (- u * u - 1 + 2 * v + 2 * u - v * v - 2 * u * v);   // (i,j,k)=(003)
      dw_du[10*i+3] = 3 * u * u;                                               // (i,j,k)=(300)
      dw_du[10*i+9] = 0;                                                       // (i,j,k)=(030)
      dw_du[10*i+1] = (1 - 4 * u - 2 * v + 4 * u * v + 3 * u * u + v * v);     // (i,j,k)=(102)
      dw_du[10*i+2] = 3 * (2 * u - 3 * u * u - 2 * u * v);                     // (i,j,k)=(201)
      dw_du[10*i+6] = 6 * u * v;                                               // (i,j,k)=(210)
      dw_du[10*i+8] = 3 * v * v;                                               // (i,j,k)=(120)
      dw_du[10*i+7] = -3 * v * v;                                              // (i,j,k)=(021)
      dw_du[10*i+4] = 3 * (1 - 2 * u - 2 * v * v + 2 * u * v);                 // (i,j,k)=(012)
      dw_du[10*i+5] = 6 * (v - 2 * u * v - v * v);                             // (i,j,k)=(111)

      dw_dv[10*i+0] = 3 * (- v * v - 1 + 2 * u + 2 * v - 2 * u * v - u * u);   // (i,j,k)=(003)
      dw_dv[10*i+3] = 0;                                                       // (i,j,k)=(300)
      dw_dv[10*i+9] = 3 * v * v;                                               // (i,j,k)=(030)
      dw_dv[10*i+1] = 3 * (-2 * u + 2 * u * u + 2 * v * u);                    // (i,j,k)=(102)
      dw_dv[10*i+2] = -3 * u * u;                                              // (i,j,k)=(201)
      dw_dv[10*i+6] = 3 * u * u;                                               // (i,j,k)=(210)
      dw_dv[10*i+8] = 6 * u * v;                                               // (i,j,k)=(120)
      dw_dv[10*i+7] = 3 * (2* v - 2 * u * v - 3 * v * v);                      // (i,j,k)=(021)
      dw_dv[10*i+4] = 3 * (1 - 2 * u - 4 * v + 4 * u * v + u * u + 3 * v * v); // (i,j,k)=(012)
      dw_dv[10*i+5] = 6 * (u - u * u - 2 * u * v);                             // (i,j,k)=(111)

    }
  }

  // if (mesh_dimension == 3) {

  //   if (order == 1) {
  //     for (int i = 0; i < n_pts; i++) {

  //       dw_du[3*i]   = 0;
  //       dw_du[3*i+1] = 1;
  //       dw_du[3*i+2] = 0;

  //       dw_dv[3*i]   = 0;
  //       dw_dv[3*i+1] = 0;
  //       dw_dv[3*i+2] = 1;

  //       dw_dw[3*i]   = 1;
  //       dw_dw[3*i+1] = 0;
  //       dw_dw[3*i+2] = 0;

  //     }
  //   }

  //   if (order == 2) {
  //     for (int i = 0; i < n_pts; i++) {

  //       double u =  uv[2*i];
  //       double v =  uv[2*i+1];

  //       dw_du[6*i+0] = 0;              // (i,j,k)=(0,0,2)
  //       dw_du[6*i+1] = 2 * (1 - u -v); // (i,j,k)=(1,0,1)
  //       dw_du[6*i+2] = 2 * u;          // (i,j,k)=(2,0,0)
  //       dw_du[6*i+3] = 0;              // (i,j,k)=(0,1,1)
  //       dw_du[6*i+4] = 2 * v;          // (i,j,k)=(1,1,0)
  //       dw_du[6*i+5] = 0;              // (i,j,k)=(0,2,0)

  //       dw_dv[6*i+0] = 0;              // (i,j,k)=(0,0,2)
  //       dw_dv[6*i+1] = 0;              // (i,j,k)=(1,0,1)
  //       dw_dv[6*i+2] = 0;              // (i,j,k)=(2,0,0)
  //       dw_dv[6*i+3] = 2 * (1 - u -v); // (i,j,k)=(0,1,1)
  //       dw_dv[6*i+4] = 2 * u;          // (i,j,k)=(1,1,0)
  //       dw_dv[6*i+5] = 2 * v;          // (i,j,k)=(0,2,0)

  //       dw_dv[6*i+0] = 2 * (1 - u -v); // (i,j,k)=(0,0,2)
  //       dw_dv[6*i+1] = 2 * u;          // (i,j,k)=(1,0,1)
  //       dw_dv[6*i+2] = 0;              // (i,j,k)=(2,0,0)
  //       dw_dv[6*i+3] = 2 * v;          // (i,j,k)=(0,1,1)
  //       dw_dv[6*i+4] = 0;              // (i,j,k)=(1,1,0)
  //       dw_dv[6*i+5] = 0;              // (i,j,k)=(0,2,0)

  //     }
  //   }

  //   if (order == 3) {

  //     for (int i = 0; i < n_pts; i++) {

  //       double u =  uv[2*i];
  //       double v =  uv[2*i+1];

  //       dw_du[10*i+0] =  0;                           // (i,j,k)=(003)
  //       dw_du[10*i+3] =  3 * u * u;                   // (i,j,k)=(300)
  //       dw_du[10*i+9] =  0;                           // (i,j,k)=(030)
  //       dw_du[10*i+1] =  3 * (1 - u -v) * (1 - u -v); // (i,j,k)=(102)
  //       dw_du[10*i+2] =  6 * u * (1 -u -v);           // (i,j,k)=(201)
  //       dw_du[10*i+6] =  6 * u * v;                   // (i,j,k)=(210)
  //       dw_du[10*i+8] =  3 *  v * v;                  // (i,j,k)=(120)
  //       dw_du[10*i+7] =  0;                           // (i,j,k)=(021)
  //       dw_du[10*i+4] =  0;                           // (i,j,k)=(012)
  //       dw_du[10*i+5] =  6 * v * (1 - u -v);          // (i,j,k)=(111)

  //       dw_dv[10*i+0] =  0;                           // (i,j,k)=(003)
  //       dw_dv[10*i+3] =  0;                           // (i,j,k)=(300)
  //       dw_dv[10*i+9] =  3 * v * v;                   // (i,j,k)=(030)
  //       dw_dv[10*i+1] =  0;                           // (i,j,k)=(102)
  //       dw_dv[10*i+2] =  0;                           // (i,j,k)=(201)
  //       dw_dv[10*i+6] =  3 * u * u;                   // (i,j,k)=(210)
  //       dw_dv[10*i+8] =  6 * u * v;                   // (i,j,k)=(120)
  //       dw_dv[10*i+7] =  6 * v * (1 - u -v);          // (i,j,k)=(021)
  //       dw_dv[10*i+4] =  3 * (1 - u -v) * (1 - u -v); // (i,j,k)=(012)
  //       dw_dv[10*i+5] =  6 * u * (1 - u -v);          // (i,j,k)=(111)

  //       dw_dv[10*i+0] =  3 * (1 - u -v) * (1 - u -v); // (i,j,k)=(003)
  //       dw_dv[10*i+3] =  0;                           // (i,j,k)=(300)
  //       dw_dv[10*i+9] =  0;                           // (i,j,k)=(030)
  //       dw_dv[10*i+1] =  6 * u * (1 - u -v);          // (i,j,k)=(102)
  //       dw_dv[10*i+2] =  3 * u * u;                   // (i,j,k)=(201)
  //       dw_dv[10*i+6] =  0;                           // (i,j,k)=(210)
  //       dw_dv[10*i+8] =  0;                           // (i,j,k)=(120)
  //       dw_dv[10*i+7] =  3 * v * v;                   // (i,j,k)=(021)
  //       dw_dv[10*i+4] =  6 * v * (1 - u -v);          // (i,j,k)=(012)
  //       dw_dv[10*i+5] =  6 * u * v;                   // (i,j,k)=(111)

  //     }
  //   }

  // } // en 3D case



  if (order > 3) {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Evaluate high-order basis Bézier functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */

void
PDM_ho_bezier_basis
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double                     *weights
)
{
 switch (type) {

  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_BARHO:
  case PDM_MESH_NODAL_BARHO_BEZIER:
    _basis_bezier_edge(order, n_pts, uvw, weights);
    break;

  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_TRIAHO:
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
    _basis_bezier_tria(order, n_pts, uvw, weights);
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet type other than BAR and TRIA\n");
    break;
  }
}

/**
 *
 * \brief Evaluate high-order basis Bézier functions derivatives
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] dw_du   Weights derivatives with respect to u
 * \param [out] dw_dv   Weights derivatives with respect to v
 * \param [out] dw_dw   Weights derivatives with respect to w
 *
 */

void
PDM_ho_bezier_basis_derivative
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double            *restrict dw_du,
 double            *restrict dw_dv,
 double            *restrict dw_dw
)
{
 switch (type) {

  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_BARHO:
    _basis_bezier_edge_derivative(order, n_pts, uvw, dw_du);
    break;

  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_TRIAHO:
    _basis_bezier_tria_derivative(order, n_pts, uvw, dw_du, dw_dv, dw_dw);
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet type other than BAR and TRIA\n");
    break;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
