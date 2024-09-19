#ifndef __FVMC_GATHER_H__
#define __FVMC_GATHER_H__

/*============================================================================
 * Base functions for parallelism
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2007  EDF

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

#include "fvmc_config.h"

#if defined(FVMC_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_io_num.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)

/*
  Pointer to structure keeping track of the status of a series of
  fvmc_gather_...() operations by slices of element global I/O number
  intervals. The structure itself is private, and is defined in fvmc_gather.c
*/

typedef struct _fvmc_gather_slice_t fvmc_gather_slice_t;

#endif /* defined(FVMC_HAVE_MPI) */

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a fvmc_gather_slice_t structure.
 *
 * parameters:
 *   entity_io_num    <-- I/O numbering structure associated with slice entity
 *   slice_size       <-- reference slice size
 *   comm             <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

fvmc_gather_slice_t *
fvmc_gather_slice_create(const fvmc_io_num_t  *entity_io_num,
                        const fvmc_gnum_t     slice_size,
                        MPI_Comm             comm);

/*----------------------------------------------------------------------------
 * Destroy a fvmc_gather_slice_t structure.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_gather_slice_t *
fvmc_gather_slice_destroy(fvmc_gather_slice_t  * this_slice);

/*----------------------------------------------------------------------------
 * Advance a fvmc_gather_slice_t structure to the next start and end values.
 *
 * Elements within this slice will be those for whose global number
 * is >= global_num_start and < global_num_end.
 *
 * parameters:
 *   this_slice        <-- pointer to structure that should be advanced
 *   global_num_start  --> new current global slice start number
 *   global_num_end    --> new current global slice past the end number
 *
 * returns:
 *   0 if the end of the slice has not been reached before this call,
 *   1 if we have already attained the end of the slice.
 *----------------------------------------------------------------------------*/

int
fvmc_gather_slice_advance(fvmc_gather_slice_t  *this_slice,
                         fvmc_gnum_t          *global_num_start,
                         fvmc_gnum_t          *global_num_end);

/*----------------------------------------------------------------------------
 * Reset a fvmc_gather_slice_t structure to its initial state.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be reinitialized
 *----------------------------------------------------------------------------*/

void
fvmc_gather_slice_reinitialize(fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Limit an fvmc_gather_slice_t structure's end value.
 *
 * This allows setting a lower global_num_end value than that previously
 * defined (which may be necessary when buffer size limits require it).
 *
 * parameters:
 *   this_slice        <-- pointer to structure that should be advanced
 *   global_num_end    --> new current global slice past the end number
 *----------------------------------------------------------------------------*/

void
fvmc_gather_slice_limit(fvmc_gather_slice_t  *this_slice,
                       fvmc_gnum_t          *global_num_end);

/*----------------------------------------------------------------------------
 * Build a slice index (0 to n-1 numbering) on rank 0 from local index arrays.
 *
 * This is done by computing the local block lengths from the local
 * index, gathering those lengths to rank 0, and rebuilding a 0 to n-1
 * numbered slice index on rank 0.
 *
 * This function is intended to be used within a loop on subsets of the global
 * lengths array (so as to enable writing to file or sending to an
 * external process without requiring the full array to reside in the process
 * directly handling I/O's memory). As such, it avoids allocating its own
 * working arrays (buffers), so that they may be allocated outside the loop
 * and reused for each call (avoiding the overhead associated with memory
 * allocation).
 *
 * All or most elements in a given portion may belong to a same process rank
 * (depending on mesh numbering and domain splitting). To account for
 * this, for each process rank, the slice_index[] arrays must be large
 * enough to contain (slice_size * stride) values, even though most processes
 * will require less.
 *
 * parameters:
 *   local_index      <-- local index array
 *   slice_index      --> global slice index section for elements
 *                        slice global_num_start to global_num_end
 *                        (output for rank 0, working array only for others)
 *   element_io_num   <-- I/O numbering structure associated with elements
 *   comm             <-- MPI communicator for structures considered
 *   this_slice       <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_slice_index(const fvmc_lnum_t     local_index[],
                       fvmc_gnum_t           slice_index[],
                       const fvmc_io_num_t  *element_io_num,
                       MPI_Comm             comm,
                       fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Recompute maximum value of global_num_end and slice connectivity size for
 * an indexed connectivity slice.
 *
 * Given an initial global connectivity buffer size associated with the
 * slice (global_connect_s_size), this function verifies that the connectivity
 * associated with the slice from global_num_start to global_num_end may fit
 * in this buffer. If this is not the case, global_num_end is reduced to the
 * largest value such that the associated indexed connectivity or values may
 * fit in the indicated buffer size.
 *
 * In any case, slice size will neither be increased above the current
 * slice size, nor be reduced to less than
 * than min(n_g_elements, n_elements_s_min) if initially larger than this.
 * If necessary, global_connect_s_size is increased so that this minimal
 * slice may fit in a buffer of this same size.
 *
 * parameters:
 *   n_elements_s_min      <-- minimum number of elements per slice desired
 *   global_num_end        --> new current global slice past the end number
 *   global_connect_s_size <-> pointer to global connectivity slice size
 *   comm                  <-- associated MPI communicator
 *   slice_index           <-- index of blocks corresponding to a given
 *                             element in the global_connect_s array
 *                             (required for rank 0 only)
 *   this_slice            <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_resize_indexed_slice(const fvmc_gnum_t     n_elements_s_min,
                                fvmc_gnum_t          *global_num_end,
                                fvmc_gnum_t          *global_connect_s_size,
                                MPI_Comm             comm,
                                const fvmc_gnum_t     slice_index[],
                                fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Gather a given portion of an array to rank 0.
 *
 * This function is intended to be used within a loop on subsets of the global
 * array (so as to enable writing to file or sending to an external process
 * without requiring the full array to reside in the process directly
 * handling I/O's memory). As such, it avoids allocating its own working arrays
 * (buffers), so that they may be allocated outside the loop and reused for
 * each call (avoiding the overhead associated with memory allocation).
 *
 * All or most elements in a given portion may belong to a same process rank
 * (depending on mesh numbering and domain splitting). To account for
 * this, for each process rank, the global_array_s[] array must be large
 * enough to contain (slice_size * stride) values, even though most processes
 * will require less.
 *
 * parameters:
 *   local_array      <-- local array (size n_local_elements * stride)
 *   global_array_s   --> global array section for elements
 *                        slice global_num_start to global_num_end
 *                        (output for rank 0, working array only for others)
 *   datatype         <-- MPI datatype of each value
 *   stride           <-- number of (interlaced) values per element
 *   element_io_num   <-- I/O numbering structure associated with elements
 *   comm             <-- MPI communicator for structures considered
 *   this_slice       <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_array(const void          *local_array,
                 void                *global_array_s,
                 MPI_Datatype         datatype,
                 size_t               stride,
                 const fvmc_io_num_t  *element_io_num,
                 MPI_Comm             comm,
                 fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Gather a given portion of an indexed array of to rank 0.
 *
 * A slice_index[] array indicating the index (0 to n-1) of blocks in
 * the slice is required for rank 0. This implies that the block sizes in
 * the slice have already been gathered through the use of
 * fvmc_gather_slice_index() or some similar method, and used to build this
 * slice index.
 *
 * This function is intended to be used within a loop on subsets of the global
 * lengths array (so as to enable writing to file or sending to an
 * external process without requiring the full array to reside in the process
 * directly handling I/O's memory). As such, it avoids allocating its own
 * working arrays (buffers), so that they may be allocated outside the loop
 * and reused for each call (avoiding the overhead associated with memory
 * allocation).
 *
 * All or most elements in a given portion may belong to a same process rank
 * (depending on mesh numbering and domain splitting). To account for
 * this, for each process rank, the global_lengths_s[] arrays must be large
 * enough to contain (slice_index[current_slice_size] - 1) values, even
 * though most processes will require less.
 * Use fvmc_gather_resize_indexed_slice() to adjust current_slice_size.
 *
 * parameters:
 *   local_array      <-- local array
 *                        (size: local_index[n_local_elements] * stride)
 *   global_array_s   --> global array section for elements
 *                        slice global_num_start to global_num_end
 *                        (output for rank 0, working array only for others)
 *   datatype         <-- MPI datatype of each value
 *   local_index      <-- local index array
 *   element_io_num   <-- I/O numbering structure associated with elements
 *   comm             <-- MPI communicator for structures considered
 *   slice_index      <-- index of blocks corresponding to a given
 *                        element in the global_numbers_s array
 *                        (required for rank 0 only)
 *   this_slice       <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_indexed(const void          *local_array,
                   void                *global_array_s,
                   const MPI_Datatype   datatype,
                   const fvmc_lnum_t     local_index[],
                   const fvmc_io_num_t  *element_io_num,
                   MPI_Comm             comm,
                   const fvmc_gnum_t     slice_index[],
                   fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Gather a given portion of a strided (i.e. regular) connectivity array
 * to rank 0. Connectivity values are converted from local to global values
 * (both with 1 to n type numbering).
 *
 * This function is intended to be used within a loop on subsets of the global
 * connectivity array (so as to enable writing to file or sending to an
 * external process without requiring the full array to reside in the process
 * directly handling I/O's memory). As such, it avoids allocating its own
 * working arrays (buffers), so that they may be allocated outside the loop
 * and reused for each call (avoiding the overhead associated with memory
 * allocation).
 *
 * All or most elements in a given portion may belong to a same process rank
 * (depending on mesh numbering and domain splitting). To account for
 * this, for each process rank, the global_connect_s[] array must be large
 * enough to contain (slice_size * stride) values, even though most processes
 * will require less.
 *
 * parameters:
 *   local_connect    <-- local connectivity array (1 to n numbering)
 *   global_connect_s --> global connectivity array section for elements
 *                        slice global_num_start to global_num_end
 *                        (output for rank 0, working array only for others)
 *   stride           <-- number of connected entities (i.e. vertices in
 *                        a nodal connectivity) per element
 *   connected_io_num <-- I/O numbering structure associated with "connected"
 *                        entities (i.e. vertices in a nodal connectivity)
 *   element_io_num   <-- I/O numbering structure associated with elements
 *   comm             <-- MPI communicator for structures considered
 *   this_slice       <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_strided_connect(const fvmc_lnum_t     local_connect[],
                           fvmc_gnum_t           global_connect_s[],
                           const int            stride,
                           const fvmc_io_num_t  *connected_io_num,
                           const fvmc_io_num_t  *element_io_num,
                           MPI_Comm             comm,
                           fvmc_gather_slice_t  *this_slice);

/*----------------------------------------------------------------------------
 * Gather a given portion of an indexed array of numbers to rank 0.
 * If the connected_io_num argument is non-NULL, these numbers
 * are assumed to represent connectivity values, and are converted from
 * local to global values (both with 1 to n type numbering).
 * Otherwise, they are considered to represent any other type of positive
 * integer (such as the number of vertices for each of a polyhedron's faces).
 *
 * A slice_index[] array indicating the index (0 to n-1) of blocks in
 * the slice is required for rank 0. This implies that the block sizes in
 * the slice have already been gathered through the use of
 * fvmc_gather_slice_index() or some similar method, and used to build this
 * slice index.
 *
 * This function is intended to be used within a loop on subsets of the global
 * lengths array (so as to enable writing to file or sending to an
 * external process without requiring the full array to reside in the process
 * directly handling I/O's memory). As such, it avoids allocating its own
 * working arrays (buffers), so that they may be allocated outside the loop
 * and reused for each call (avoiding the overhead associated with memory
 * allocation).
 *
 * All or most elements in a given portion may belong to a same process rank
 * (depending on mesh numbering and domain splitting). To account for
 * this, for each process rank, the global_lengths_s[] arrays must be large
 * enough to contain (slice_index[slice_size] - 1) values, even though most
 * processes will require less.
 * Use fvmc_gather_resize_indexed_slice() to adjust current_slice_size.
 *
 * parameters:
 *   local_index      <-- local index array
 *   local_numbers    <-- local numbers array
 *   global_numbers_s --> global numbers array section for elements
 *                        slice global_num_start to global_num_end
 *                        (output for rank 0, working array only for others)
 *   connected_io_num <-- I/O numbering structure associated with "connected"
 *                        entities (i.e. vertices in a nodal connectivity)
 *   element_io_num   <-- I/O numbering structure associated with elements
 *   comm             <-- MPI communicator for structures considered
 *   slice_index      <-- index of blocks corresponding to a given
 *                        element in the global_numbers_s array
 *                        (required for rank 0 only)
 *   this_slice       <-> structure for management of slice status
 *----------------------------------------------------------------------------*/

void
fvmc_gather_indexed_numbers(const fvmc_lnum_t     local_index[],
                           const fvmc_lnum_t     local_numbers[],
                           fvmc_gnum_t           global_numbers_s[],
                           const fvmc_io_num_t  *connected_io_num,
                           const fvmc_io_num_t  *element_io_num,
                           MPI_Comm             comm,
                           const fvmc_gnum_t     slice_index[],
                           fvmc_gather_slice_t  *this_slice);

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_GATHER_H__ */
