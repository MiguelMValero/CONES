#ifndef __PDM_MORTON_CUH__
#define __PDM_MORTON_CUH__

/*============================================================================
 * Morton encoding for 2D or 3D coordinates.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA
  Copyright (C) 2008-2010  EDF

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

#include <stdio.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   level    <-- level in the grid
 *   extents  <-- coordinate extents for normalization (size: dim*2)
 *   n_coords <-- nomber of coordinates in array
 *   coords   <-- coordinates in the grid (interlaced, not normalized)
 *   m_code   --> array of corresponding Morton codes
 *   d        --> Normalization (dilatation component)
 *   s        --> Normalization (translation component)
 *----------------------------------------------------------------------------*/

void
PDM_morton_encode_coords_CPU(int                dim,
                             PDM_morton_int_t   level,
                             const double       extents[],
                             size_t             n_coords,
                             double             *d_pts,
                             PDM_morton_code_t  pts_code[],
                             double             d[3],
                             double             s[3]);

__device__
void
PDM_morton_encode_coords_GPU(int                dim,
                         PDM_morton_int_t   level,
                         const double       extents[],
                         size_t             n_coords,
                         const double       coords[],
                         PDM_morton_code_t  m_code[],
                         double             d[3],
                         double             s[3]);                            

/*----------------------------------------------------------------------------
 * Get the index associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the array.
 *
 * parameters:
 *   size  <-- size of the array
 *   code  <-- code we are searching for
 *   codes <-- array of Morton codes
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

__device__
int
PDM_morton_binary_search_GPU(int                size,
                             PDM_morton_code_t  code,
                             PDM_morton_code_t *codes);


/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

__device__
_Bool
PDM_morton_a_gt_b_GPU(PDM_morton_code_t  a,
                      PDM_morton_code_t  b);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MORTON_H__ */
