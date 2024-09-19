/*============================================================================
 * Linear algebra routines
 *============================================================================*/

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

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

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_config.h"
#include "pdm_priv.h"


#include "pdm_linear_algebra.h"


#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) ((a) == 0.0 ? 0.0 : (a) * (a))

/*============================================================================
 * Type definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")


/*=============================================================================
 * Private function definitions
 *============================================================================*/

// calculates sqrt( a^2 + b^2 ) with decent precision
static double pythag(double a, double b) {
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);

  if (absa > absb)
    return (absa * sqrt(1.0 + SQR(absb/absa)));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

/**
 * \brief Perform "back substitution" to solve the linear system of equations
 * after a Singular Value Decomposition.
 *
 * u * diag(w) * v *  x = b
 *
 * u : n_row*n_col matrix (uij = u[n_col*i+j])
 * w : vector of singular values (size = n_col)
 * v : n_col*n_col matrix (vij = v[n_col*i+j])
 * b : n_row*stride matrix (bij = b[stride*i+j])
 * x : n_col*stride matrix (xij = x[stride*i+j])
 *
 * \param [in]    n_row    Number of rows    in matrix u
 * \param [in]    n_col    Number of columns in matrix u
 * \param [in]    stride   Number of columns in matrix b
 * \param [in]    tol      Tolerance for singular value truncation (relative to greatest singular value)
 * \param [in]    u        Matrix U (size = n_row * n_col)
 * \param [in]    w        Array of singular values (size = n_col)
 * \param [in]    v        Matrix Vt (size = n_col * n_col)
 * \param [in]    b        Right-hand side of the system (size = n_row * stride)
 * \param [out]   x        Solution to the system (size = n_col * stride)
 *
 */

static void svbksb
(
 const int     n_row,
 const int     n_col,
 const int     stride,
 const double  tol,
       double *u,
       double *w,
       double *v,
       double *b,
       double *x
)
{
  double wtol = 0;
  for (int i = 0; i < n_col; i++) {
    wtol = PDM_MAX(wtol, tol*PDM_ABS(w[i]));
  }

  double y[n_col*stride];

  for (int i = 0; i < n_col*stride; i++) {
    y[i] = 0;
  }

  for (int j = 0; j < n_col; j++) {
    if (w[j] > wtol) {
      double iw = 1/w[j];
      for (int i = 0; i < n_row; i++) {
        double wuij = u[n_col*i+j] * iw;
        for (int k = 0; k < stride; k++) {
          y[stride*j+k] += wuij * b[stride*i+k];
        }
      }
    }
  }

  for (int i = 0; i < n_col*stride; i++) {
    x[i] = 0;
  }

  for (int j = 0; j < n_col; j++) {
    for (int i = 0; i < n_col; i++) {
      for (int k = 0; k < stride; k++) {
        x[stride*j+k] += v[n_col*j+i] * y[stride*i+k];
      }
    }
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Compute Singular Value Decomposition of a rectangular matrix
 * (Adapted from Numerical Recipes)
 *
 * Given a n_row*n_col matrix A (n_col <= n_row), the SVD A = U * diag(w) * Vt is computed.
 * A is overwritten by U.
 * (matrices are stored in C-order, i.e. Aij = A[n_col*i + j])
 *
 * \param [in]    n_row  Number of rows    in matrix A
 * \param [in]    n_col  Number of columns in matrix A
 * \param [inout] a      Rectangular matrix to decompose, overwritten by U (size = n_row * n_col)
 * \param [out]   w      Array of singular values (size = n_col)
 * \param [out]   v      Matrix Vt (size = n_col * n_col)
 *
 * \return 0 if converged, 1 else
 *
 */

int
PDM_linear_algebra_svd
(
 const int     n_row,
 const int     n_col,
       double *a,
       double *w,
       double *v
 )
{
  double rv1[n_col];
  double g     = 0;
  double scale = 0;
  double anorm = 0;
  double c, f, h, s, x, y, z;

  int l;

  int stat = 0;

  for (int i = 0; i < n_col; i++) {
    l = i+1;

    rv1[i] = scale * g;

    g     = 0;
    s     = 0;
    scale = 0;

    if (i < n_row) {
      for (int k = i; k < n_row; k++) {
        scale += PDM_ABS(a[n_col*k+i]);
      }

      if (scale) {
        double iscale = 1./scale;
        for (int k = i; k < n_row; k++) {
          a[n_col*k+i] *= iscale;
          s += a[n_col*k+i] * a[n_col*k+i];
        }

        f = a[n_col*i+i];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        double ih = 1./h;
        a[n_col*i+i] = f - g;
        for (int j = l; j < n_col; j++) {
          s = 0;
          for (int k = i; k < n_row; k++) {
            s += a[n_col*k+i] * a[n_col*k+j];
          }
          f = s * ih;
          for (int k = i; k < n_row; k++) {
            a[n_col*k+j] += f * a[n_col*k+i];
          }
        }

        for (int k = i; k < n_row; k++) {
          a[n_col*k+i] *= scale;
        }
      }
    }

    w[i] = scale * g;

    g     = 0;
    s     = 0;
    scale = 0;

    if (i < n_row && i != n_col-1) {
      for (int k = l; k < n_col; k++) {
        scale += PDM_ABS(a[n_col*i+k]);
      }

      if (scale) {
        double iscale = 1./scale;
        for (int k = l; k < n_col; k++) {
          a[n_col*i+k] *= iscale;
          s += a[n_col*i+k]*a[n_col*i+k];
        }
        f = a[n_col*i+l];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        a[n_col*i+l] = f - g;

        double ih = 1./h;
        for (int k = l; k < n_col; k++) {
          rv1[k] = a[n_col*i+k] * ih;
        }

        for (int j = l; j < n_row; j++) {
          s = 0;
          for (int k = l; k < n_col; k++ ) {
            s += a[n_col*j+k] * a[n_col*i+k];
          }
          for (int k = l; k < n_col; k++) {
            a[n_col*j+k] += s * rv1[k];
          }
        }

        for (int k = l; k < n_col; k++) {
          a[n_col*i+k] *= scale;
        }
      }

    }
    anorm = PDM_MAX(anorm, (PDM_ABS(w[i]) + PDM_ABS(rv1[i])));
  }

  for (int i = n_col-1; i >= 0; i--) {
    if (i < n_col-1) {
      if (g) {
        for (int j = l; j < n_col; j++) {
          v[n_col*j+i] = (a[n_col*i+j] / a[n_col*i+l]) / g;
        }

        for (int j = l; j < n_col; j++) {
          s = 0;
          for (int k = l; k < n_col; k++) {
            s += a[n_col*i+k] * v[n_col*k+j];
          }
          for (int k = l; k < n_col; k++) {
            v[n_col*k+j] += s * v[n_col*k+i];
          }
        }
      }

      for (int j = l; j < n_col; j++) {
        v[n_col*i+j] = 0;
        v[n_col*j+i] = 0;
      }
    }

    v[n_col*i+i] = 1;
    g = rv1[i];
    l = i;
  }

  for (int i = PDM_MIN(n_row, n_col) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];

    for (int j = l; j < n_col; j++) {
      a[n_col*i+j] = 0;
    }

    if (g) {
      g = 1./g;

      for (int j = l; j < n_col; j++) {
        s = 0;
        for (int k = l; k < n_row; k++) {
          s += a[n_col*k+i] * a[n_col*k+j];
        }

        f = (s / a[n_col*i+i]) * g;

        for (int k = i; k < n_row; k++) {
          a[n_col*k+j] += f * a[n_col*k+i];
        }
      }

      for (int j = i; j < n_row; j++) {
        a[n_col*j+i] *= g;
      }
    }
    else {
      for (int j = i; j < n_row; j++) {
        a[n_col*j+i] = 0;
      }
    }
    a[n_col*i+i] += 1;
  }

  int nm = -1;
  for (int k = n_col-1; k >= 0; k--) {
    for (int its = 0; its < 30; its++) {
      int flag = 1;

      for (l = k; l >= 0; l--) {
        nm = l - 1;

        if ((PDM_ABS(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }

        if ((PDM_ABS(w[nm] + anorm)) == anorm) {
          break;
        }
      }


      if (flag) {
        c = 0;
        s = 1;

        for (int i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];

          if ((PDM_ABS(f) + anorm) == anorm) {
            break;
          }

          g = w[i];
          h = pythag(f, g);
          w[i] = h;
          h = 1./h;
          c =  g * h;
          s = -f * h;

          for (int j = 0; j < n_row; j++) {
            y = a[n_col*j+nm];
            z = a[n_col*j+i];

            a[n_col*j+nm] = y*c + z*s;
            a[n_col*j+i ] = z*c - y*s;
          }
        }
      }

      z = w[k];
      if (l == k) {
        if (z < 0) {
          w[k] = -z;
          for (int j = 0; j < n_col; j++) {
            v[n_col*j+k] = -v[n_col*j+k];
          }
        }
        break;
      }

      if (its == 29) {
        stat = 1;
        printf("no convergence in 30 svdcmp iterations\n");
      }

      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f)))- h)) / x;
      c = 1;
      s = 1;

      for (int j = l; j <= nm; j++) {
        int i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (int jj = 0; jj < n_col; jj++) {
          x = v[n_col*jj+j];
          z = v[n_col*jj+i];
          v[n_col*jj+j] = x*c + z*s;
          v[n_col*jj+i] = z*c - x*s;
        }
        z = pythag(f, h);
        w[j] = z;
        if (z) {
          z = 1./z;
          c = f * z;
          s = h * z;
        }
        f = c*g + s*y;
        x = c*y - s*g;
        for (int jj = 0; jj < n_row; jj++) {
          y = a[n_col*jj+j];
          z = a[n_col*jj+i];
          a[n_col*jj+j] = y * c + z * s;
          a[n_col*jj+i] = z * c - y * s;
        }
      }
      rv1[l] = 0;
      rv1[k] = f;
      w[k]   = x;
    }
  }

  return stat;
}


/**
 * \brief Solve the linear system ax = b using Singular Value Decomposition
 *
 * a : n_row * n_col matrix of the linear system
 * b : n_row * stride matrix (right-hand side term)
 * x : n_col * stride matrix (solution)
 * (a is overwritten in the process)
 *
 * \param [in]    n_row    Number of rows    in matrix A
 * \param [in]    n_col    Number of columns in matrix A
 * \param [in]    stride   Number of columns in matrix b
 * \param [in]    tol      Tolerance for singular value truncation (relative to greatest singular value)
 * \param [inout] a        Rectangular matrix (overwritten) (size = n_row * n_col)
 * \param [in]    b        Right-hand side of the system (size = n_row * stride)
 * \param [out]   x        Solution to the system (if SVD converged) (size = n_col * stride)
 *
 * \return 0 if SVD converged, 1 else
 *
 */

int
PDM_linear_algebra_linsolve_svd
(
 const int     n_row,
 const int     n_col,
 const int     stride,
 const double  tol,
       double *a,
       double *b,
       double *x
 )
{
  double w[n_col];
  double v[n_col*n_col];

  int stat = PDM_linear_algebra_svd(n_row,
                                    n_col,
                                    a,
                                    w,
                                    v);

  if (stat == 0) {
    svbksb(n_row,
           n_col,
           stride,
           tol,
           a,
           w,
           v,
           b,
           x);
  }

  return stat;
}



/**
 * \brief Solve the square linear system Ax = b using Gaussian elimination,
 * where A is a n*n matrix and b, x are n*stride matrices
 * (Aij = A[n*i+j], bij = b[stride*i+j], xij = x[stride*i+j])
 *
 * /!\ Gaussian elimination is performed in place
 * (A and x are used as work arrays)
 *
 * \param [in]    n  Number of rows and columns
 * \param [inout] A  Matrix (overwritten) (size = n * n)
 * \param [inout] x  Right-hand side term at input, solution at output (size = n * stride)
 *
 * \return 1 if A is singular, 0 else
 */

int
PDM_linear_algebra_linsolve_gauss
(
 const int     n,
 const int     stride,
       double *A,
       double *x
 )
{
  const double eps = 1e-15;

  for (int i = 0; i < n; i++) {
    /* Seek best pivot */
    double amax = PDM_ABS(A[n*i+i]);
    int imax = i;
    for (int k = i+1; k < n; k++) {
      double aki = PDM_ABS(A[n*k+i]);
      if (aki > amax) {
        amax = aki;
        imax = k;
      }
    }

    if (amax <= eps) {
      /* matrix A is singular */
      return 1;
    }

    /* Swap rows i and imax */
    if (i != imax) {
      for (int j = 0; j < n; j++) {
        double tmp = A[n*i+j];
        A[n*i   +j] = A[n*imax+j];
        A[n*imax+j] = tmp;
      }

      for (int j = 0; j < stride; j++) {
        double tmp = x[stride*i + j];
        x[stride*i    + j] = x[stride*imax + j];
        x[stride*imax + j] = tmp;
      }
    }

    /* Eliminate subdiagonal terms */
    double inv_amax = 1./A[n*i+i];

    for (int k = i+1; k < n; k++) {
      double r = A[n*k+i] * inv_amax;
      for (int j = i+1; j < n; j++) {
        A[n*k+j] -= r * A[n*i+j];
      }
      A[n*k+i] = 0.;

      for (int j = 0; j < stride; j++) {
        x[stride*k + j] -= r * x[stride*i + j];
      }
    }
  }

  /* Solve triangular system */
  for (int i = n-1; i >= 0; i--) {
    for (int j = i+1; j < n; j++) {
      for (int k = 0; k < stride; k++) {
        x[stride*i + k] -= x[stride*j + k] * A[n*i+j];
      }
    }

    double inv_ai = 1./A[n*i+i];
    for (int k = 0; k < stride; k++) {
      x[stride*i + k] *= inv_ai;
    }
  }

  return 0;
}




PDM_GCC_SUPPRESS_WARNING_POP

#ifdef  __cplusplus
}
#endif
#undef SIGN
#undef SQR
