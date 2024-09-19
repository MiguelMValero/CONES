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
#include "solve_ax_b_4.h"

#define ABS(a)     ((a) <  0  ? -(a) : (a))

/*----------------------------------------------------------------------------
 * Solve Ax = B for a 4x4 system.
 *
 * parameters:
 *   a                  <-- matrix
 *   b                  <-- right hand side
 *   x                  --> solution of Ax = b
 *
 * returns: 0 if solution was found, 1 in case of error (zero pivot)
 *----------------------------------------------------------------------------*/

int solve_ax_b_4(double            a[4][4],
              double  *restrict b,
              double  *restrict x)
{
  int i, j, k, k_pivot;
  
  double abs_pivot, abs_a_ki, factor;
  double _a[4][4], _b[4];
  double _a_swap[4], _b_swap;

  const double _epsilon = 1.0e-15;

  /* Copy array and RHS (avoid modifying a and b) */

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      _a[i][j] = a[i][j];
    }
    _b[i] = b[i];
  }

  /* Forward elimination */

  for (i = 0; i < 4-1; i++) {

    /* Seek best pivot */

    k_pivot = i;
    abs_pivot = ABS(_a[i][i]);

    for (k = i+1; k < 4; k++) {
      abs_a_ki = ABS(_a[k][i]);
      if (abs_a_ki > abs_pivot) {
        abs_pivot = abs_a_ki;
        k_pivot = k;
      }
    }

    /* Swap if necessary */

    if (k_pivot != i) {

      for (j = 0; j < 4; j++) {
        _a_swap[j] = _a[i][j];
        _a[i][j] = _a[k_pivot][j];
        _a[k_pivot][j] = _a_swap[j];
      }

      _b_swap = _b[i];
      _b[i] = _b[k_pivot];
      _b[k_pivot] = _b_swap;

    }

    if (abs_pivot < _epsilon)
      return 1;

    /* Eliminate values */

    for (k = i+1; k < 4; k++) {

      factor = _a[k][i] / _a[i][i];

      _a[k][i] = 0.0;
      for (j = i+1; j < 4; j++) {
        _a[k][j] -= _a[i][j]*factor;
      }
      _b[k] -= _b[i]*factor;

    }

  }

  /* Solve system */

  x[3] =  _b[3]                                                   / _a[3][3];
  x[2] = (_b[2] - _a[2][3]*x[3])                                  / _a[2][2];
  x[1] = (_b[1] - _a[1][3]*x[3]  - _a[1][2]*x[2])                 / _a[1][1];
  x[0] = (_b[0] - _a[0][3]*x[3]  - _a[0][2]*x[2] - _a[0][1]*x[1]) / _a[0][0];

  return 0;
}

