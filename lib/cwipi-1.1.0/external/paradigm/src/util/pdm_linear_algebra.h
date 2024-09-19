#ifndef __PDM_LINEAR_ALGEBRA_H
#define __PDM_LINEAR_ALGEBRA_H

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/


/*============================================================================
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
 );


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
);


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
 );


#ifdef  __cplusplus
}
#endif

#endif // PDM_LINEAR_ALGEBRA_H


