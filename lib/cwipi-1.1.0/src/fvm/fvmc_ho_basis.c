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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bftc_error.h"
#include "bftc_mem.h"
#include "bftc_printf.h"
#include "fvmc_point_location.h"
#include "cwipi_config.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_ho_basis.h"

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

typedef struct _fvmc_ho_basis_user_elt_t {
  fvmc_ho_basis_fct_t elt_basis;
} fvmc_ho_basis_user_elt_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static fvmc_ho_basis_user_elt_t *_user_edge = NULL;

static fvmc_ho_basis_user_elt_t *_user_tria = NULL;
static fvmc_ho_basis_user_elt_t *_user_quad = NULL;

static fvmc_ho_basis_user_elt_t *_user_tetra = NULL;
static fvmc_ho_basis_user_elt_t *_user_hexa = NULL;
static fvmc_ho_basis_user_elt_t *_user_prism = NULL;
static fvmc_ho_basis_user_elt_t *_user_pyra = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Monomonial product
 * Procedure utilisee pour calculer (rapidement) les fonctions de base des simplex
 * dont les noeuds d'interpolation sont equidistants.
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   u               <-- u (size =  n_pts)
 *   fn              --> fn (size =  n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_monomialProduct
(
 const int order,
 const int n_pts,
 const int i_pts,
 const double *restrict u,
       double *restrict fn
)
{

  double constant;

  for (int i = 0; i < n_pts; i++) {
    fn[i]= 1.;
  }

  for (int i = 0; i < i_pts; i++) {

    constant = 1. / (double) (i - i_pts);

    for (int j = 0; j < n_pts; j++) {
      fn[j] *= ((double) i - (double) order * u[j]) * constant;
    }
  }
}



/*----------------------------------------------------------------------------
 *
 * Edge Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   u              <-- u (size = n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_edge_pn
(
 const int order,
 const int n_pts,
 const double *restrict u,
 double *restrict weights
)
{

  const int n_nodes = order + 1;
  double uNodes[n_nodes];


  for (int i = 0; i<n_nodes; i++) {
    uNodes[i] = (double) i / (double) order;
  }


  if (order == 1) {

    for (int i = 0; i < n_pts; i++) {
      double _u =  u[i];
      weights[2*i]   = 1. - _u;
      weights[2*i+1] = _u;
    }
  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {

     double _u = u[i];

     weights[3*i+0] =                    (_u - uNodes[1]) * (_u - uNodes[2]) / ((uNodes[0] - uNodes[1]) * (uNodes[0] - uNodes[2]));
     weights[3*i+1] = (_u - uNodes[0]) *                    (_u - uNodes[2]) / ((uNodes[1] - uNodes[0]) * (uNodes[1] - uNodes[2]));
     weights[3*i+2] = (_u - uNodes[0]) * (_u - uNodes[1])                    / ((uNodes[2] - uNodes[0]) * (uNodes[2] - uNodes[1]));
    }
  }

  else if (order == 3) {

    for (int i = 0;  i < n_pts; i++) {

      double _u = u[i];

      weights[4*i+0] =                     (_u - uNodes[1]) * (_u - uNodes[2]) * (_u - uNodes[3]) / ((uNodes[0] - uNodes[1]) * (uNodes[0] - uNodes[2]) * (uNodes[0] - uNodes[3]));
      weights[4*i+1] =  (_u - uNodes[0]) *                    (_u - uNodes[2]) * (_u - uNodes[3]) / ((uNodes[1] - uNodes[0]) * (uNodes[1] - uNodes[2]) * (uNodes[1] - uNodes[3]));
      weights[4*i+2] =  (_u - uNodes[0]) * (_u - uNodes[1]) *                    (_u - uNodes[3]) / ((uNodes[2] - uNodes[0]) * (uNodes[2] - uNodes[1]) * (uNodes[2] - uNodes[3]));
      weights[4*i+3] =  (_u - uNodes[0]) * (_u - uNodes[1]) * (_u - uNodes[2])                    / ((uNodes[3] - uNodes[0]) * (uNodes[3] - uNodes[1]) * (uNodes[3] - uNodes[2]));

    }
  }

  else if (order == 4) {

    for (int i = 0;  i < n_pts; i++) {

      double _u = u[i];

      weights[5*i+0]= 10.666666666666666*(-1. + _u)*(-0.75 + _u)*(-0.5 + _u)*(-0.25 + _u);
      weights[5*i+1]= -42.666666666666664*(-1. + _u)*(-0.75 + _u)*(-0.5 + _u)*(0. + _u);
      weights[5*i+2]= 64.*(-1. + _u)*(-0.75 + _u)*(-0.25 + _u)*(0. + _u);
      weights[5*i+3]= -42.666666666666664*(-1. + _u)*(-0.5 + _u)*(-0.25 + _u)*(0. + _u);
      weights[5*i+4]= 10.666666666666666*(-0.75 + _u)*(-0.5 + _u)*(-0.25 + _u)*(0. + _u);


    }
  }

  else {

    for (int i = 0; i < n_pts; i++) {

      double _u = u[i];

      for (int i_node = 0; i_node < n_nodes; i_node++) {

        weights[n_nodes * i + i_node] = 1.0;
        for (int n = 0; n < n_nodes; n++) {

          if (n != i_node) {

            double coef = (double) 1 / (double) (uNodes[i_node] - uNodes[n]);
            weights[n_nodes * i + i_node] *= _u - uNodes[n];
            weights[n_nodes * i + i_node] *= coef;
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------
 *
 * Triangle Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_tria_pn
(
 const int order,
 const int n_pts,
 const double *restrict uv,
 double *restrict weights
)
{
  const int n_nodes = (order + 2) * (order + 1) / 2;

  if (order == 1) {

    for (int i = 0; i < n_pts; i++) {
      double u =  uv[2*i];
      double v =  uv[2*i+1];
      weights[3*i]   = 1. - u - v;
      weights[3*i+1] = u;
      weights[3*i+2] = v;
    }
  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {

     double u  = uv[2*i];
     double v  = uv[2*i+1];
     double w  = 1. - u - v;
     double u2 = 2. * u;
     double v2 = 2. * v;
     double w2 = 2. * w;

     weights[6*i+0] = w * (-1. + w2);  /* (i,j,k)=(0,0,2) */
     weights[6*i+1] = u2 * w2;         /* (i,j,k)=(1,0,1) */
     weights[6*i+2] = u * (-1. + u2);  /* (i,j,k)=(2,0,0) */
     weights[6*i+3] = v2 * w2;         /* (i,j,k)=(0,1,1) */
     weights[6*i+4] = u2 * v2;         /* (i,j,k)=(1,1,0) */
     weights[6*i+5] = v * (-1. + v2);  /* (i,j,k)=(0,2,0) */
    }
  }

  else if (order == 3) {

    for (int i = 0;  i < n_pts; i++) {

      double u = uv[2*i];
      double u3 = 3.*u;
      double u3m1 = (u3-1.)*5e-1;

      double v = uv[2*i+1];
      double v3 = 3.*v;
      double v3m1 = (v3-1.)*5e-1;

      double w = 1. - u - v;
      double w3 = 3.*w;
      double w3m1 = (w3-1.)*5e-1;

      weights[10*i+0] = w*w3m1*(w3-2.);   // (i,j,k)=(003)
      weights[10*i+3] = u*u3m1*(u3-2.);   // (i,j,k)=(300)
      weights[10*i+9] = v*v3m1*(v3-2.);   // (i,j,k)=(030)
      weights[10*i+1] = u3*w3*w3m1;       // (i,j,k)=(102)

      double coef = u3*u3m1;
      weights[10*i+2] = coef*w3;           //(i,j,k)=(201)
      weights[10*i+6] = coef*v3;           // (i,j,k)=(210)

      coef=v3*v3m1;
      weights[10*i+8] = coef*u3;           // (i,j,k)=(120)
      weights[10*i+7] = coef*w3;           // (i,j,k)=(021)

      coef=v3*w3;
      weights[10*i+4] = coef*w3m1;         // (i,j,k)=(012)
      weights[10*i+5] = coef*u3;           // (i,j,k)=(111)

    }
  }

  else {

    double *u  = malloc (sizeof(double) *n_pts);
    double *v  = malloc (sizeof(double) *n_pts);
    double *w  = malloc (sizeof(double) *n_pts);
    double *fu = malloc (sizeof(double) *n_pts);
    double *fv = malloc (sizeof(double) *n_pts);
    double *fw = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];
      w[i] = 1. - u[i] - v[i];
    }

    int i_node = 0;
    for (int iv = 0; iv < order + 1; iv++) {
      for (int iu = 0; iu < order + 1 - iv; iu++) {
        int iw = order - iu - iv;

        _monomialProduct(order, n_pts, iu, u, fu);
        _monomialProduct(order, n_pts, iv, v, fv);
        _monomialProduct(order, n_pts, iw, w, fw);

        for (int i = 0; i < n_pts; i++) {
          weights[i * n_nodes + i_node] = fu[i] * fv[i] * fw[i];
        }
        i_node++;
      }
    }

    free (u);
    free (v);
    free (w);
    free (fu);
    free (fv);
    free (fw);

  }
}


/*----------------------------------------------------------------------------
 *
 * Compte uv of edge nodes
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   u               <-- u (size = n_pts)
 *   weights         --> weights (size = (order + 1) * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_uNodesEdges
(
 const int order,
 double *restrict xi
)
{
  const int n_nodes = order + 1;

  for (int i = 0; i < n_nodes; i++) {
    xi[i] = -1. + 2. * (double) i / (double) order;
  }
}

/*----------------------------------------------------------------------------
 *
 * Edge basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   u               <-- u (size = n_pts)
 *   weights         --> weights (size = (order + 1) * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_setL2BasisEqui
(
 const int order,
 const int n_pts,
 const double *restrict u,
       double *restrict weights
 )
{

  const int nMod = order + 1;

  double *xi = malloc (sizeof(double) * nMod);

  _uNodesEdges (order, xi);

  const int sWeights = n_pts * nMod;
  for (int i = 0; i < sWeights; i++) {
    weights[i] = 1.;
  }

  for (int i = 0; i < nMod; i++) {
    for (int j = 0; j < nMod; j++) {

      if (i != j) {
        double var = 1. / (xi[i] - xi[j]);

        for (int i_pts = 0; i_pts < n_pts; i_pts++) {
          weights[nMod * i_pts + i] *=
            (u[i_pts] - xi[j]) * var;
        }
      }
    }
  }

  free (xi);
}


/*----------------------------------------------------------------------------
 *
 * Quadrangle Qn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/


static void
_basis_quad_qn
(
 const int order,
 const int n_pts,
 const double *restrict uv,
 double *restrict weights
)
{

  if (order == 1) {

    for (int i = 0; i < n_pts; i++) {
      double u = uv[2*i];
      double v = uv[2*i+1];
      double u1 = (1 - u);
      double v1 = (1 - v);

      weights[4*i+0] = u1 * v1;
      weights[4*i+1] = u * v1;
      weights[4*i+2] = u1 * v;
      weights[4*i+3] = u * v;
    }

  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {
      double u = uv[2*i];
      double v = uv[2*i+1];

      double uM = 2*(1-u);
      double uP = 2*u;
      double u0 = u-0.5;

      double au1 = -uM * u0;
      double au2 =  uM * uP;
      double au3 =  u0 * uP;

      double vM = 2*(1-v);
      double vP = 2*v;
      double v0 = v-0.5;

      double av1 = -vM * v0;
      double av2 =  vM * vP;
      double av3 =  v0 * vP;

      weights[9*i+0]=au1*av1;
      weights[9*i+1]=au2*av1;
      weights[9*i+2]=au3*av1;
      weights[9*i+3]=au1*av2;
      weights[9*i+4]=au2*av2;
      weights[9*i+5]=au3*av2;
      weights[9*i+6]=au1*av3;
      weights[9*i+7]=au2*av3;
      weights[9*i+8]=au3*av3;
    }

  }

  else {

    const int nMod = order + 1;
    const int n_nodes = nMod * nMod;

    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = 2 * uv[2*i]   - 1;
      v[i] = 2 * uv[2*i+1] - 1;
    }

    double *lagrangeL2_u = malloc (sizeof(double) * nMod * n_pts);
    double *lagrangeL2_v = malloc (sizeof(double) * nMod * n_pts);

    _setL2BasisEqui (order, n_pts, u, lagrangeL2_u);
    _setL2BasisEqui (order, n_pts, v, lagrangeL2_v);

    int i_node = 0;
    for (int iv = 0; iv < nMod; iv++) {
      for (int iu = 0; iu < nMod; iu++) {
        for (int i_pts = 0; i_pts < n_pts; i_pts++) {
          weights[i_pts * n_nodes + i_node] =
            lagrangeL2_u[i_pts * nMod + iu] *
            lagrangeL2_v[i_pts * nMod + iv];
        }
        i_node++;
      }
    }

    free (lagrangeL2_u);
    free (lagrangeL2_v);
    free (u);
    free (v);
  }
}



/*----------------------------------------------------------------------------
 *
 * Tetrahedron Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uvw             <-- u (size = 3 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_tetra_pn
(
  const int order,
  const int n_pts,
  const double *restrict uvw,
  double *restrict weights
)
{

  int n_nodes = (order+1)*(order+2)*(order+3)/6;


  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      double u = uvw[3*i  ];
      double v = uvw[3*i+1];
      double w = uvw[3*i+2];
      double t = 1 - u - v - w;


      weights[4*i  ] = t;  /* (0, 0, 0)*/
      weights[4*i+1] = u;  /* (1, 0, 0)*/
      weights[4*i+2] = v;  /* (0, 1, 0)*/
      weights[4*i+3] = w;  /* (0, 0, 1)*/
    }
  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {


     double u  = uvw[3*i  ];
     double v  = uvw[3*i+1];
     double w  = uvw[3*i+2];
     double t = 1 - u - v - w;


     double u2 = 2. * u;
     double v2 = 2. * v;
     double w2 = 2. * w;
     double t2 = 2. * t;

     weights[10*i+0] = t * (t2 - 1);    // (0  , 0  , 0  )
     weights[10*i+1] = u2 * t2;         // (0.5, 0  , 0  )
     weights[10*i+2] = u * (-1. + u2);  // (1  , 0  , 0  )
     weights[10*i+3] = v2 * t2;         // (0  , 0.5, 0  )
     weights[10*i+4] = u2 * v2;         // (0.5, 0.5, 0  )
     weights[10*i+5] = v * (-1. + v2);  // (0  , 1  , 0  )
     weights[10*i+6] = w2 * t2;         // (0  , 0  , 0.5)
     weights[10*i+7] = u2 * w2;         // (0.5, 0  , 0.5)
     weights[10*i+8] = v2 * w2;         // (0  , 0.5, 0.5)
     weights[10*i+9] = w * (-1. + w2);  // (0  , 0  , 1  )



    }
  }
  else {

    double *u  = malloc (sizeof(double) *n_pts);
    double *v  = malloc (sizeof(double) *n_pts);
    double *w  = malloc (sizeof(double) *n_pts);
    double *t  = malloc (sizeof(double) *n_pts);
    double *fu = malloc (sizeof(double) *n_pts);
    double *fv = malloc (sizeof(double) *n_pts);
    double *fw = malloc (sizeof(double) *n_pts);
    double *ft = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uvw[3*i];
      v[i] = uvw[3*i+1];
      w[i] = uvw[3*i+2];
      t[i] = 1. - u[i] - v[i] - w[i];
    }

    int i_node = 0;
    for (int iw = 0; iw < order + 1; iw++) {
      for (int iv = 0; iv < order + 1 - iw; iv++) {
        for (int iu = 0; iu < order + 1 - iv - iw; iu++) {
          int it = order - iu - iv - iw;

          _monomialProduct(order, n_pts, iu, u, fu);
          _monomialProduct(order, n_pts, iv, v, fv);
          _monomialProduct(order, n_pts, iw, w, fw);
          _monomialProduct(order, n_pts, it, t, ft);

          for (int i = 0; i < n_pts; i++) {
            weights[i * n_nodes + i_node] = fu[i] * fv[i] * fw[i] * ft[i];
          }
          i_node++;
        }
      }
    }
    free (u);
    free (v);
    free (w);
    free (t);
    free (fu);
    free (fv);
    free (fw);
    free (ft);
  }

}



/*----------------------------------------------------------------------------
 *
 * Hexahedron Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uvw             <-- u (size = 3 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_hexa_pn
(
  const int order,
  const int n_pts,
  const double *restrict uvw,
  double *restrict weights
)
{

  int n_nodes = (order+1)*(order+1)*(order+1);

  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      double u = uvw[3*i];
      double v = uvw[3*i+1];
      double w = uvw[3*i+2];
      double u1 = (1 - u);
      double v1 = (1 - v);
      double w1 = (1 - w);

      weights[8*i+0] = u1 * v1 * w1;
      weights[8*i+1] = u  * v1 * w1;
      weights[8*i+2] = u1 * v  * w1;
      weights[8*i+3] = u  * v  * w1;
      weights[8*i+4] = u1 * v1 * w;
      weights[8*i+5] = u  * v1 * w;
      weights[8*i+6] = u1 * v  * w;
      weights[8*i+7] = u  * v  * w;
    }
  }

  if (order == 2) {

    for (int i = 0; i < n_pts; i++) {
      double u = uvw[3*i];
      double v = uvw[3*i+1];
      double w = uvw[3*i+2];

      double uM = 2*(1-u);
      double uP = 2*u;
      double u0 = u-0.5;

      double au1 = -uM * u0;
      double au2 =  uM * uP;
      double au3 =  u0 * uP;

      double vM = 2*(1-v);
      double vP = 2*v;
      double v0 = v-0.5;

      double av1 = -vM * v0;
      double av2 =  vM * vP;
      double av3 =  v0 * vP;

      double wM = 2*(1-w);
      double wP = 2*w;
      double w0 = w-0.5;

      double aw1 = -wM * w0;
      double aw2 =  wM * wP;
      double aw3 =  w0 * wP;

      weights[27*i+ 0]=au1*av1*aw1;
      weights[27*i+ 1]=au2*av1*aw1;
      weights[27*i+ 2]=au3*av1*aw1;
      weights[27*i+ 3]=au1*av2*aw1;
      weights[27*i+ 4]=au2*av2*aw1;
      weights[27*i+ 5]=au3*av2*aw1;
      weights[27*i+ 6]=au1*av3*aw1;
      weights[27*i+ 7]=au2*av3*aw1;
      weights[27*i+ 8]=au3*av3*aw1;
      weights[27*i+ 9]=au1*av1*aw2;
      weights[27*i+10]=au2*av1*aw2;
      weights[27*i+11]=au3*av1*aw2;
      weights[27*i+12]=au1*av2*aw2;
      weights[27*i+13]=au2*av2*aw2;
      weights[27*i+14]=au3*av2*aw2;
      weights[27*i+15]=au1*av3*aw2;
      weights[27*i+16]=au2*av3*aw2;
      weights[27*i+17]=au3*av3*aw2;
      weights[27*i+18]=au1*av1*aw3;
      weights[27*i+19]=au2*av1*aw3;
      weights[27*i+20]=au3*av1*aw3;
      weights[27*i+21]=au1*av2*aw3;
      weights[27*i+22]=au2*av2*aw3;
      weights[27*i+23]=au3*av2*aw3;
      weights[27*i+24]=au1*av3*aw3;
      weights[27*i+25]=au2*av3*aw3;
      weights[27*i+26]=au3*av3*aw3;
    }
  }
  else {

    int nMod = order + 1;

    double *u  = malloc (sizeof(double) *n_pts);
    double *v  = malloc (sizeof(double) *n_pts);
    double *w  = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = 2 * uvw[3*i]   - 1;
      v[i] = 2 * uvw[3*i+1] - 1;
      w[i] = 2 * uvw[3*i+2] - 1;
    }

    double *lagrangeL2_u = malloc (sizeof(double) * nMod * n_pts);
    double *lagrangeL2_v = malloc (sizeof(double) * nMod * n_pts);
    double *lagrangeL2_w = malloc (sizeof(double) * nMod * n_pts);

    _setL2BasisEqui (order, n_pts, u, lagrangeL2_u);
    _setL2BasisEqui (order, n_pts, v, lagrangeL2_v);
    _setL2BasisEqui (order, n_pts, w, lagrangeL2_w);

    int i_node = 0;
    for (int iw = 0; iw < nMod; iw++) {
      for (int iv = 0; iv < nMod; iv++) {
        for (int iu = 0; iu < nMod; iu++) {
          for (int i_pts = 0; i_pts < n_pts; i_pts++) {
            weights[i_pts * n_nodes + i_node] =
              lagrangeL2_u[i_pts * nMod + iu] *
              lagrangeL2_v[i_pts * nMod + iv] *
              lagrangeL2_w[i_pts * nMod + iw];
          }
          i_node++;
        }
      }
    }
    free(u);
    free(v);
    free(w);
    free(lagrangeL2_u);
    free(lagrangeL2_v);
    free(lagrangeL2_w);
  }
}


/*----------------------------------------------------------------------------
 *
 * Prism Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uvw             <-- u (size = 3 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_prism_pn
(
  const int order,
  const int n_pts,
  const double *restrict uvw,
  double *restrict weights
)
{

  int n_nodes = (order+1)*(order+1)*(order+2)/2;

  if (order == 1) {

    for (int i = 0; i < n_pts; i++) {
      double u =  uvw[3*i];
      double v =  uvw[3*i+1];
      double w =  uvw[3*i+2];

      double w1 = 1 - w;

      weights[6*i]   = (1. - u - v) * w1;
      weights[6*i+1] = u * w1;
      weights[6*i+2] = v * w1;
      weights[6*i+3] = (1. - u - v) * w;
      weights[6*i+4] = u * w;
      weights[6*i+5] = v * w;
    }
  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {

     double u  = uvw[3*i];
     double v  = uvw[3*i+1];
     double w  = uvw[3*i+2];


     double t  = 1. - u - v;
     double u2 = 2. * u;
     double v2 = 2. * v;
     double t2 = 2. * t;

     double wM = 2*(1-w);
     double wP = 2*w;
     double w0 = w-0.5;
     double aw1 = -wM * w0;
     double aw2 =  wM * wP;
     double aw3 =  w0 * wP;

     weights[18*i+ 0] = (t * (-1. + t2)) * aw1;
     weights[18*i+ 1] = (u2 * t2)        * aw1;
     weights[18*i+ 2] = (u * (-1. + u2)) * aw1;
     weights[18*i+ 3] = (v2 * t2)        * aw1;
     weights[18*i+ 4] = (u2 * v2)        * aw1;
     weights[18*i+ 5] = (v * (-1. + v2)) * aw1;

     weights[18*i+ 6] = (t * (-1. + t2)) * aw2;
     weights[18*i+ 7] = (u2 * t2)        * aw2;
     weights[18*i+ 8] = (u * (-1. + u2)) * aw2;
     weights[18*i+ 9] = (v2 * t2)        * aw2;
     weights[18*i+10] = (u2 * v2)        * aw2;
     weights[18*i+11] = (v * (-1. + v2)) * aw2;

     weights[18*i+12] = (t * (-1. + t2)) * aw3;
     weights[18*i+13] = (u2 * t2)        * aw3;
     weights[18*i+14] = (u * (-1. + u2)) * aw3;
     weights[18*i+15] = (v2 * t2)        * aw3;
     weights[18*i+16] = (u2 * v2)        * aw3;
     weights[18*i+17] = (v * (-1. + v2)) * aw3;

    }
  }
  else {

    int nMod = order + 1;

    double *u   = malloc (sizeof(double) *n_pts);
    double *v   = malloc (sizeof(double) *n_pts);
    double *w   = malloc (sizeof(double) *n_pts);
    double *t   = malloc (sizeof(double) *n_pts);
    double *fu  = malloc (sizeof(double) *n_pts);
    double *fv  = malloc (sizeof(double) *n_pts);
    double *ft  = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uvw[3*i];
      v[i] = uvw[3*i+1];
      w[i] = 2*uvw[3*i+2]-1;
      t[i] = 1 - u[i] - v[i];
    }

    double *lagrangeL2_w = malloc (sizeof(double) * nMod * n_pts);

    _setL2BasisEqui (order, n_pts, w, lagrangeL2_w);

    int i_node = 0;
    for (int iw = 0; iw < nMod; iw++) {
      for (int iv = 0; iv < nMod; iv++) {
        for (int iu = 0; iu < nMod - iv; iu++) {
          int it = order - iu - iv;

          _monomialProduct(order, n_pts, iu, u, fu);
          _monomialProduct(order, n_pts, iv, v, fv);
          _monomialProduct(order, n_pts, it, t, ft);


          for (int i_pts = 0; i_pts < n_pts; i_pts++) {
            weights[i_pts * n_nodes + i_node] = fu[i_pts] * fv[i_pts] * ft[i_pts] * lagrangeL2_w[i_pts * nMod + iw];
          }
          i_node++;
        }
      }
    }
    free(u);
    free(v);
    free(w);
    free(t);
    free(fu);
    free(fv);
    free(ft);
    free(lagrangeL2_w);
  }
}



/*----------------------------------------------------------------------------
 *
 * Pyramid Pn basis
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uvw             <-- u (size = 3 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_pyra_pn
(
  const int order,
  const int n_pts,
  const double *restrict uvw,
  double *restrict weights
)
{

  if (order == 1) {
    for (int i = 0; i < n_pts; i++) {

      double u = uvw[3*i];
      double v = uvw[3*i+1];
      double w = uvw[3*i+2];

      weights[5*i+0] = (1-u-w)*(1-v-w)/(1-w);
      weights[5*i+1] = u*(1-v-w)/(1-w);
      weights[5*i+2] = (1-u-w)*v/(1-w);
      weights[5*i+3] = u*v/(1-w);
      weights[5*i+4] = w;
    }
  }

  else if (order == 2) {

    for (int i = 0; i < n_pts; i++) {
      double e1 = uvw[3*i];
      double e2 = uvw[3*i+1];
      double e5 = uvw[3*i+2];
      double e3 = 1 - e1 - e5;
      double e4 = 1 - e2 - e5;

      weights[14*i+ 8] = e1*e2*( ((2*e1/(1-e5) - 1) * (2*e2/(1-e5) - 1)) - e5/(1-e5) );
      weights[14*i+ 6] = e2*e3*( ((2*e2/(1-e5) - 1) * (2*e3/(1-e5) - 1)) - e5/(1-e5) );
      weights[14*i+ 0] = e3*e4*( ((2*e3/(1-e5) - 1) * (2*e4/(1-e5) - 1)) - e5/(1-e5) );
      weights[14*i+ 2] = e4*e1*( ((2*e4/(1-e5) - 1) * (2*e1/(1-e5) - 1)) - e5/(1-e5) );

      weights[14*i+ 5] = 4*(e2*e4/(1-e5))*e1*( (2*e1/(1-e5)) - 1 );
      weights[14*i+ 7] = 4*(e3*e1/(1-e5))*e2*( (2*e2/(1-e5)) - 1 );
      weights[14*i+ 3] = 4*(e4*e2/(1-e5))*e3*( (2*e3/(1-e5)) - 1 );
      weights[14*i+ 1] = 4*(e1*e3/(1-e5))*e4*( (2*e4/(1-e5)) - 1 );

      weights[14*i+ 4] = 16*e1*e2*e3*e4/((1-e5)*(1-e5));

      weights[14*i+12] = 4*e1*e2*e5/(1-e5);
      weights[14*i+11] = 4*e2*e3*e5/(1-e5);
      weights[14*i+ 9] = 4*e3*e4*e5/(1-e5);
      weights[14*i+10] = 4*e4*e1*e5/(1-e5);

      weights[14*i+13] = e5*(2*e5-1);



    }
  }
  else {
    bftc_error(__FILE__, __LINE__, 0,
               _("Ordre de la pyramide non implémenté\n"));
  }

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

static fvmc_ho_basis_user_elt_t **
_get_user_elt (fvmc_element_t elt_type)
{

  switch(elt_type) {

  case FVMC_EDGE:
    return &_user_edge;
    break;

  case FVMC_FACE_TRIA:
    return &_user_tria;
    break;

  case FVMC_FACE_QUAD:
    return &_user_quad;
    break;

  case FVMC_CELL_TETRA:
    return &_user_tetra;
    break;

  case FVMC_CELL_PYRAM:
    return &_user_pyra;
    break;

  case FVMC_CELL_PRISM:
    return &_user_prism;
    break;

  case FVMC_CELL_HEXA:
    return &_user_hexa;
    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
               _("fvmc_ho_user_elt_unset : Unvailable element type\n"));
  }

  return NULL;

}

/*----------------------------------------------------------------------------
 *
 * default high order basis
 *
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uvw             <-- uvw (size = elt_dim * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

static void
_default_elt_basis
(
const fvmc_element_t type,
const int order,
const int n_pts,
const double *uvw,
      double *weights
)
{
  switch (type) {

  case FVMC_FACE_TRIA:
    _basis_tria_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_FACE_QUAD:
    _basis_quad_qn (order, n_pts, uvw, weights);
    break;

  case FVMC_EDGE:
    _basis_edge_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_CELL_TETRA:
    _basis_tetra_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_CELL_PRISM:
    _basis_prism_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_CELL_PYRAM:
    _basis_pyra_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_CELL_HEXA:
    _basis_hexa_pn (order, n_pts, uvw, weights);
    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
               _("FVMC_ho_basis : '%d' element type not yet implemented\n"),
               type);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Unset a user element
 *
 *----------------------------------------------------------------------------*/

void
FVMC_ho_basis_user_elt_unset (fvmc_element_t elt_type)
{

  fvmc_ho_basis_user_elt_t **_user_elt = _get_user_elt (elt_type);

  if (*_user_elt != NULL) {
    free (*_user_elt);
    *_user_elt = NULL;
  }

}


void
FVMC_ho_basis_user_elts_unset (void)
{
  FVMC_ho_basis_user_elt_unset (FVMC_EDGE);
  FVMC_ho_basis_user_elt_unset (FVMC_FACE_TRIA);
  FVMC_ho_basis_user_elt_unset (FVMC_FACE_QUAD);
  FVMC_ho_basis_user_elt_unset (FVMC_CELL_TETRA);
  FVMC_ho_basis_user_elt_unset (FVMC_CELL_HEXA);
  FVMC_ho_basis_user_elt_unset (FVMC_CELL_PRISM);
  FVMC_ho_basis_user_elt_unset (FVMC_CELL_PYRAM);
}

/*----------------------------------------------------------------------------
 *
 * Unset a user element
 *
 *----------------------------------------------------------------------------*/


void
FVMC_ho_basis_user_elt_set (fvmc_element_t elt_type,
                            fvmc_ho_basis_fct_t elt_basis)
{
  fvmc_ho_basis_user_elt_t **user_elt = _get_user_elt (elt_type);

  if (*user_elt == NULL) {
    *user_elt = (fvmc_ho_basis_user_elt_t *) malloc (sizeof(fvmc_ho_basis_user_elt_t));
  }

  (*user_elt)->elt_basis = elt_basis;

}


/*----------------------------------------------------------------------------
 *
 * high order basis
 *
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_nodes         <-- number of nodes
 *   n_pts           <-- number of points
 *   uvw             <-- uvw
 *   uvw             <-- uvw (size = n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

void
FVMC_ho_basis
(
const fvmc_element_t type,
const int order,
const int n_nodes,
const int n_pts,
const double *uvw,
      double *weights
)
{
  fvmc_ho_basis_user_elt_t *user_elt = *(_get_user_elt (type));

  int entities_dim = 0;
  switch(type) {
  case FVMC_EDGE:
    entities_dim = 1;
    break;
  case FVMC_FACE_TRIA:          /* Triangle */
  case FVMC_FACE_QUAD:          /* Quadrangle */
    entities_dim = 2;
    break;
  case FVMC_CELL_TETRA:         /* Tetrahedron */
  case FVMC_CELL_PYRAM:         /* Pyramid */
  case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
  case FVMC_CELL_HEXA:          /* Hexahedron (brick) */
    entities_dim = 3;
    break;
  default:
    bftc_error(__FILE__, __LINE__, 0,
               "%d is not high order element type\n", type);

  }

  // printf("FVMC_ho_basis %ld\n",user_elt );

  if (user_elt != NULL) {
    //printf("FVMC_ho_basis 2 %ld\n",user_elt->elt_basis );
    if (user_elt->elt_basis != NULL) {
      (user_elt->elt_basis) (entities_dim,
                             order,
                             n_nodes,
                             n_pts,
                             uvw,
                             weights);
    }
    else {
      _default_elt_basis (type,
                          order,
                          n_pts,
                          uvw,
                          weights);
    }
  }

  else {
    _default_elt_basis (type,
                        order,
                        n_pts,
                        uvw,
                        weights);

  }
}

/*----------------------------------------------------------------------------
 *
 * Free static variables
 *
 *----------------------------------------------------------------------------*/


void
FVMC_ho_basis_free
(
 void
)
{
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
