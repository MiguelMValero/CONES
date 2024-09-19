#ifndef __PDM_CUDA_CUH__
#define __PDM_CUDA_CUH__

/*============================================================================
 * CUDA functions
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017 ONERA

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

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

#include "pdm_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/


//Kernel reduction
#define Reduce_kernel(n_threads, n_blocks, shared_size, blockSize, length, data_in, data_out) {reduce_kernel6<<<n_blocks,n_threads,shared_size>>>(data_in, data_out, blockSize, length);}
inline __device__ void warpReduce
(
  volatile int *sdata, 
  unsigned int tid,
  int blockSize
  ) 
{
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
  if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
  if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
  if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

inline __global__ void reduce_kernel6
(
  int *g_idata, 
  int *g_odata,
  int blockSize,
  int length
  ) 
{
  extern __shared__ int sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  if ((i >= length))  
  {
    sdata[tid] = 0;
  }
  else if ((i + blockDim.x) >= length)
  {
    sdata[tid] = g_idata[i];
  }
  else
  {
    sdata[tid] = g_idata[i] + g_idata[i + blockDim.x];
  }
  __syncthreads();


  if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }

  if (tid < 32) warpReduce(sdata, tid, blockSize);

  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0];
  }
}

//CUDA reallocation
#define cudaRealloc(ptr, oldLength, newLength)                \
({                                                            \
  __typeof__(ptr) newptr = NULL;                              \
  if (newLength == 0)                                         \
  {                                                           \
    gpuErrchk(cudaFree(ptr));                                 \
    newptr = NULL;                                            \
  }                                                           \
  else if (!ptr)                                              \
  {                                                           \
    gpuErrchk(cudaMalloc(&ptr, sizeof(newLength)));           \
    newptr = ptr;                                             \
  }                                                           \
  else if (newLength <= oldLength)                            \
  {                                                           \
    newptr = ptr;                                             \
  }                                                           \
  else                                                        \
  {                                                           \
    assert((ptr) && (newLength > oldLength));                 \
    gpuErrchk(cudaMalloc(&newptr, newLength));                \
    if (newptr)                                               \
    {                                                         \
      memcpy(newptr, ptr, oldLength);                         \
      gpuErrchk(cudaFree(ptr));                               \
    }                                                         \
  }                                                           \
  newptr;                                                     \
})

/*============================================================================
 * Public types
 *============================================================================*/


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*!
* \brief set dim3 value, for compatibility with intel compiler
*
*
* \param [in]        x
* \param [in]        y
* \param [in]        z
*
*/

__host__ __device__
dim3
set_dim3_value(int x, int y, int z);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CUDA_CUH_ */