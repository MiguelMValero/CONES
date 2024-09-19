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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_morton.h"
#include "pdm_morton.cuh"
#include "pdm_binary_search.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_cuda_error.cuh"
#include "pdm_cuda.cuh"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
inline static _Bool
_a_gt_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[1] << a_diff;
    code_a.X[2] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = ((code_b.X[0] + 1) << b_diff) - 1;
    code_b.X[1] = ((code_b.X[1] + 1) << b_diff) - 1;
    code_b.X[2] = ((code_b.X[2] + 1) << b_diff) - 1;
  }

  i = l - 1;
  while (i > 0) {

    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

   a =   ((code_a.X[0] >> i) % 2) * 4
       + ((code_a.X[1] >> i) % 2) * 2
       + ((code_a.X[2] >> i) % 2);
   b =   ((code_b.X[0] >> i) % 2) * 4
       + ((code_b.X[1] >> i) % 2) * 2
       + ((code_b.X[2] >> i) % 2);

  return (a > b) ? true : false;
}


/*============================================================================
 * CUDA kernel definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Encode coordinates.
 *
 * parameters:
 *   dim         <-- 1D, 2D or 3D
 *   level       <-- level in the grid
 *   n_coords    <-- nomber of coordinates in array
 *   coords      <-- coordinates in the grid (interlaced, not normalized)
 *   m_code      <-> array of corresponding Morton codes
 *   d           --> Normalization (dilatation component)
 *   s           --> Normalization (translation component)
 *   n           <-> ?
 *   refinement  <-- ?
 *   switch_case <-- dimension case
 *----------------------------------------------------------------------------*/

__global__
static
void
_encode_coords_1
(
  int                dim,
  PDM_morton_int_t   level,
  size_t             n_coords,
  const double       coords[],
  PDM_morton_code_t  m_code[],
  double             d[3],
  double             s[3],
  double             n[3],
  PDM_morton_int_t   refinement
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n_coords)
  {
    return;
  }

  m_code[i].L = level;
  n[0] = (coords[i] - s[0]) / d[0];
  m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
  m_code[i].X[1] = 0;
  m_code[i].X[2] = 0;
}



__global__
static
void
_encode_coords_2
(
  int                dim,
  PDM_morton_int_t   level,
  size_t             n_coords,
  const double       coords[],
  PDM_morton_code_t  m_code[],
  double             d[3],
  double             s[3],
  double             n[3],
  PDM_morton_int_t   refinement
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if ((i >= n_coords) || (j >= 2))
  {
    return;
  }

  m_code[i].L = level;
  n[j] = (coords[i*dim + j] - s[j]) / d[j];
  m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
  m_code[i].X[2] = 0;
}


__global__
static
void
_encode_coords_3
(
  int                dim,
  PDM_morton_int_t   level,
  size_t             n_coords,
  const double       coords[],
  PDM_morton_code_t  m_code[],
  double             d[3],
  double             s[3],
  double             n[3],
  PDM_morton_int_t   refinement
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if ((i >= n_coords) || (j >= 3))
  {
    return;
  }

  m_code[i].L = level;
  
  n[j] = (coords[i*dim + j] - s[j]) / d[j];
  m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
}

// __global__
// static
// void
// _encode_coords
// (
//   int                dim,
//   PDM_morton_int_t   level,
//   size_t             n_coords,
//   const double       coords[],
//   PDM_morton_code_t  m_code[],
//   double             d[3],
//   double             s[3],
//   double             n[3],
//   PDM_morton_int_t   refinement,
//   int                switch_case
// )
// {
//   int i = blockIdx.x * blockDim.x + threadIdx.x;
//   int j = blockIdx.y * blockDim.y + threadIdx.y;

 
//   if (switch_case == 3)
//   {
//     if ((i >= n_coords) || (j >= 3))
//     {
//       return;
//     }

//     m_code[i].L = level;
    
//     n[j] = (coords[i*dim + j] - s[j]) / d[j];
//     m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
//   }
//   else if (switch_case == 2)
//   {
//     if ((i >= n_coords) || (j >= 2))
//     {
//       return;
//     }

//     m_code[i].L = level;
//     n[j] = (coords[i*dim + j] - s[j]) / d[j];
//     m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
//     m_code[i].X[2] = 0;
//   }
//   else if (switch_case == 1)
//   {
//     if (i >= n_coords)
//     {
//       return;
//     }

//     m_code[i].L = level;
//     n[0] = (coords[i] - s[0]) / d[0];
//     m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
//     m_code[i].X[1] = 0;
//     m_code[i].X[2] = 0;
//   }
//   else
//   {
//     printf("Error : invalid switch_case argument\n");
//     __trap();
//   }
  
// }



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
 *   pts_code   --> array of corresponding Morton codes
 *   d        --> Normalization (dilatation component)
 *   s        --> Normalization (translation component)
 *----------------------------------------------------------------------------*/


void
PDM_morton_encode_coords_CPU(int                dim,
                             PDM_morton_int_t   level,
                             const double       extents[],
                             size_t             n_coords,
                             double             *d_pts,
                             PDM_morton_code_t  d_pts_code[],
                             double             d[3],
                             double             s[3])
{
  size_t i;
  double d_max = 0.0;

  PDM_morton_int_t  refinement = 1u << level;

  for (i = 0; i < (size_t)dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
    d_max = PDM_MAX(d_max, d[i]);
  }

  for (i = 0; i < (size_t)dim; i++) { /* Reduce effective dimension */
    if (d[i] < d_max * 1e-10)
      d[i] = d_max * 1e-10;
  }

  double *d_d = NULL;
  gpuErrchk(cudaMalloc(&d_d, sizeof(double) * 3));
  double *d_s = NULL;
  gpuErrchk(cudaMalloc(&d_s, sizeof(double) * 3));
  double *d_n = NULL;
  gpuErrchk(cudaMalloc(&d_n, sizeof(double) * 3));

  gpuErrchk(cudaMemcpy(d_d, d, sizeof(double) * 3, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_s, s, sizeof(double) * 3, cudaMemcpyHostToDevice));

  dim3 n_threads;
  dim3 n_blocks;

  switch(dim) {

    case 3:
      n_threads = set_dim3_value(32, 32, 1);
      n_blocks = set_dim3_value((n_coords + 31)/32, (3 + 31)/32, 1);
      _encode_coords_1<<<n_blocks,n_threads>>>(dim, level, n_coords, d_pts, d_pts_code, d_d, d_s, d_n, refinement);
      break;

    case 2:
      n_threads = set_dim3_value(32, 32, 1);
      n_blocks = set_dim3_value((n_coords + 31)/32, (2 + 31)/32, 1);
      _encode_coords_2<<<n_blocks,n_threads>>>(dim, level, n_coords, d_pts, d_pts_code, d_d, d_s, d_n, refinement);
      break;

    case 1:
      n_threads = set_dim3_value(32, 32, 1);
      n_blocks = set_dim3_value((n_coords + 31)/32, (1 + 31)/32, 1);
      _encode_coords_3<<<n_blocks,n_threads>>>(dim, level, n_coords, d_pts, d_pts_code, d_d, d_s, d_n, refinement);
      break;

    default:
      assert(dim > 0 && dim < 4);
      break;
  }

  gpuErrchk(cudaMemcpy(d, d_d, 3 * sizeof(double), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(s, d_s, 3 * sizeof(double), cudaMemcpyDeviceToHost));

  gpuErrchk(cudaFree(d_d));
  gpuErrchk(cudaFree(d_s));
  gpuErrchk(cudaFree(d_n));
}

__device__
void
PDM_morton_encode_coords_GPU(int                dim,
                             PDM_morton_int_t   level,
                             const double       extents[],
                             size_t             n_coords,
                             const double       coords[],
                             PDM_morton_code_t  m_code[],
                             double             d[3],
                             double             s[3])
{
  size_t i, j;
  double n[3];
  double d_max = 0.0;

  PDM_morton_int_t  refinement = 1u << level;

  for (i = 0; i < (size_t)dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
    d_max = PDM_MAX(d_max, d[i]);
  }

  for (i = 0; i < (size_t)dim; i++) { /* Reduce effective dimension */
    if (d[i] < d_max * 1e-10)
      d[i] = d_max * 1e-10;
  }

  switch(dim) {

  case 3:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 3; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
    }
    break;

  case 2:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 2; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
      m_code[i].X[2] = 0;
    }
    break;

  case 1:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      n[0] = (coords[i] - s[0]) / d[0];
      m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
      m_code[i].X[1] = 0;
      m_code[i].X[2] = 0;
    }
    break;

  default:
    assert(dim > 0 && dim < 4);
    break;
  }
}

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
                             PDM_morton_code_t *codes)
{
  int start = 0;
  int end = size;

  while (end - start > 1) {

    int  middle = (end - start)/2 + start;

    /* PDM_morton_dump( 3, codes[middle]); */
    /* PDM_morton_dump( 3, code); */
    /* printf ("start end middle : %d %d %d\n", start, end, middle); */
    if (_a_gt_b(codes[middle], code))
      end = middle;
    else
      start = middle;

  }

  return start;
}

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
                      PDM_morton_code_t  b)
{
  return  _a_gt_b(a, b);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
