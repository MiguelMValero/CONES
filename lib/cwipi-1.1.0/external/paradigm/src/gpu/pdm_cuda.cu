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

#include "pdm_config.h"

/*
* Standard C library headers
*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
* Optional library and BFT headers
*/

#include "pdm_printf.h"
#include "pdm_cuda_error.cuh"
#include "pdm_cuda.cuh"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*============================================================================
* Public function definitions
*============================================================================*/

/*!
* \brief set dim3 value, for compatibility with intel compiler
*
*
* \param [inout]     value
* \param [in]        x
* \param [in]        y
* \param [in]        z
*
*/

__host__ __device__
dim3
set_dim3_value(int x, int y, int z)
{
  dim3 value;
  value.x = x;
  value.y = y;
  value.z = z;
  return value;
}

/*-----------------------------------------------------------------------------*/

// //To use on a device, need compute capability >= 3.5 (so cudaMalloc and cudaFree can be called from device)
// //Broken function as is, DO NOT use it in the code, or all cuda API calls will return error 999: unknown error
// inline
// __device__
// void*
// cudaRealloc
// (
//   void* ptr, 
//   size_t oldLength, 
//   size_t newLength
//   )
// {

//   if (newLength == 0)
//   {
//     gpuErrchk(cudaFree(ptr));
//     return NULL;
//   }
//   else if (!ptr)
//   {
//     gpuErrchk(cudaMalloc(&ptr, sizeof(newLength)));
//     return ptr;
//   }
//   else if (newLength <= oldLength)
//   {
//     return ptr;
//   }
//   else
//   {
//     assert((ptr) && (newLength > oldLength));
//     void* newptr = NULL;
//     gpuErrchk(cudaMalloc(&newptr, sizeof(newLength)));
//     if (newptr)
//     {
//       memcpy(newptr, ptr, oldLength);
//       gpuErrchk(cudaFree(ptr));
//     }
//     return newptr;
//   }
// }

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

