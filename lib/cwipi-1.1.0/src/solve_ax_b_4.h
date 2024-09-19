#ifndef __SOLVE_AX_B_4_H__
#define __SOLVE_AX_B_4_H__
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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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

//int solve_ax_b_4(double a[4][4], double  *restrict b, double  *restrict x);
int solve_ax_b_4(double a[4][4], double  * b, double  * x);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BAR_COORDS_H__ */
