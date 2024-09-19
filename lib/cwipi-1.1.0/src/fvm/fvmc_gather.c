/*============================================================================
 * Functions related to the gathering of local arrays based on global ranks
 * in parallel mode.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2009  EDF

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_io_num.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_gather.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Structure keeping track of the status of a series of fvmc_gather_...()
 * operations by slices of element global I/O number intervals.
 *----------------------------------------------------------------------------*/

struct _fvmc_gather_slice_t {

  int         local_rank;      /* Local rank in communicator */
  int         n_ranks;         /* Number of ranks in communicator */

  /* Initial and final global element numbers for all slices */

  fvmc_gnum_t  global_num_initial;
  fvmc_gnum_t  global_num_final;

  /* First and past the last global element numbers for current slice */

  fvmc_gnum_t  ref_slice_size;
  fvmc_gnum_t  global_num_slice_start;
  fvmc_gnum_t  global_num_slice_end;

  /* Local indexes corresponding to first element of a series with
     global number >= global_num_slice_start initially, and
     global number > global_num_end after a gather operation
     (initialized to 0 before first of a series of calls) */

  fvmc_lnum_t  local_index_start;   /* Current start value */
  fvmc_lnum_t  local_index_last;    /* Last value reached, to advance */

  /* Number of local entities (constant after initialization) */

  fvmc_lnum_t  local_index_end;

  /* Next global element number expected for each rank before and after
     current call (rank 0 only); initialized upon the first call of
     a series of fvmc_gather_...() operations. Keeping track of this
     allows us to avoid having to receive empty messages from
     processors having no element in a given slice, which could
     represent the majority of MPI calls for high processor counts */

  fvmc_gnum_t  *next_global_num;      /* Current values */
  fvmc_gnum_t  *next_global_num_last; /* Next values, to advance */

  /* We may only use the next global element number to avoid
     sending and receiving empty messages when it is up to date;
     this is not yet the case when a slice has just been initialized
     or resized and no operations on it have been done yet */

  _Bool  use_next_global_num;

  /* Buffers associated with MPI Datatype creation */

  size_t       recv_buf_size;  /* Receive buffer size, in bytes */
  void        *recv_buf;       /* Receive buffer */

  int         *blocklengths;   /* Required only for indexed operations */
  fvmc_gnum_t  *displacements;  /* Always required for all ranks; size
                                  slice_size + 1 (last position used to
                                  exchange next_global_num[] values) */

};

#endif /* defined(FVMC_HAVE_MPI) */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Check size of receive buffer for array data, resizing if necessary.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be destroyed
 *   n_entities <-- number of entities associated with receive
 *   stride     <-- number of values per entity
 *   type size  <-- sizeof(datatype)
 *----------------------------------------------------------------------------*/

static void
_slice_recv_buf_size_array(fvmc_gather_slice_t  * this_slice,
                           size_t                n_entities,
                           size_t                stride,
                           size_t                type_size)
{
  size_t size_mult = stride*type_size;

  _Bool reallocate = false;

  /* Base resizing */

  if (this_slice->recv_buf_size < this_slice->ref_slice_size*size_mult) {
    this_slice->recv_buf_size = this_slice->ref_slice_size*size_mult;
    reallocate = true;
  }

  /* Safety resizing (should not be necessary) */

  while (this_slice->recv_buf_size < n_entities*size_mult) {
    this_slice->recv_buf_size *= 2;
    reallocate = true;
  }

  if (reallocate == true)
    BFTC_REALLOC(this_slice->recv_buf, this_slice->recv_buf_size, char);
}

/*----------------------------------------------------------------------------
 * Check size of receive buffer for indexed data, resizing if necessary.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be destroyed
 *   n_values   <-- number of values to receive
 *   type size  <-- sizeof(datatype)
 *----------------------------------------------------------------------------*/

static void
_slice_recv_buf_size_indexed(fvmc_gather_slice_t  * this_slice,
                             size_t                n_values,
                             size_t                type_size)
{
  size_t recv_size = n_values*type_size;
  _Bool reallocate = false;

  /* Base resizing */

  if (  this_slice->recv_buf_size < this_slice->ref_slice_size*type_size) {
    this_slice->recv_buf_size = this_slice->ref_slice_size * type_size;
    reallocate = true;
  }

  /* Safety resizing (may be necessary) */

  while (this_slice->recv_buf_size < recv_size) {
    this_slice->recv_buf_size *= 2;
    reallocate = true;
  }

  if (reallocate == true)
    BFTC_REALLOC(this_slice->recv_buf, this_slice->recv_buf_size, char);
}

#if 0 && defined(DEBUG) && !defined(NDEBUG)

/*----------------------------------------------------------------------------
 * Print information on a fvmc_gather_slice_t structure.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be printed
 *----------------------------------------------------------------------------*/

static void
fvmc_gather_slice_dump(fvmc_gather_slice_t  * this_slice)
{
  if (this_slice != NULL) {

    bftc_printf("slice info:\n"
               "  adress                 = %p\n"
               "  local_rank             = %d\n"
               "  n_ranks                = %d\n"
               "  global_num_initial     = %lu\n"
               "  global_num_final       = %lu\n"
               "  ref_slice_size         = %lu\n"
               "  global_num_slice_start = %lu\n"
               "  global_num_slice_end   = %lu\n"
               "  local_index_start      = %ld\n"
               "  local_index_last       = %ld\n"
               "  local_index_end        = %ld\n"
               "  use_next_global_num    = %d\n",
               this_slice,
               this_slice->local_rank,
               this_slice->n_ranks,
               (long)(this_slice->local_index_end),
               (unsigned long)(this_slice->global_num_initial),
               (unsigned long)(this_slice->global_num_final),
               (unsigned long)(this_slice->ref_slice_size),
               (unsigned long)(this_slice->global_num_slice_start),
               (unsigned long)(this_slice->global_num_slice_end),
               (long)(this_slice->local_index_start),
               (long)(this_slice->local_index_last),
               (long)(this_slice->local_index_end),
               (int)(this_slice->use_next_global_num));

    if (this_slice->next_global_num != NULL) {
      int i;
      bftc_printf("    rank | next_global_num | next_global_num_last\n");
      for (i = 0; i < this_slice->n_ranks; i++) {
        bftc_printf(" %7d |    %12lu |   %12lu\n",
                   i,
                   (unsigned long)(this_slice->next_global_num[i]),
                   (unsigned long)(this_slice->next_global_num_last[i]));
      }
    }
  }
}

#endif /* 0 && defined(DEBUG) && !defined(NDEBUG) */

#endif /* defined(FVMC_HAVE_MPI) */

/*=============================================================================
 * Public function definitions
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
                        MPI_Comm             comm)
{
  int i;
  int  local_rank, n_ranks;
  fvmc_gather_slice_t  *this_slice;

  /* Get local rank and size of the current MPI communicator */
  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &n_ranks);

  /* Allocate and initialize slice structure */

  BFTC_MALLOC(this_slice, 1, fvmc_gather_slice_t);

  this_slice->local_rank = local_rank;
  this_slice->n_ranks = n_ranks;

  this_slice->global_num_initial = 1;
  this_slice->global_num_final = fvmc_io_num_get_global_count(entity_io_num);

  this_slice->ref_slice_size = slice_size;
  this_slice->global_num_slice_start = 1;
  this_slice->global_num_slice_end = 1;

  this_slice->local_index_end = fvmc_io_num_get_local_count(entity_io_num);

  this_slice->local_index_start = 0;
  this_slice->local_index_last = 0;

  /* Allocate and initialize "next expected global number" arrays */

  if (local_rank == 0) {

    BFTC_MALLOC(this_slice->next_global_num, n_ranks, fvmc_gnum_t);
    BFTC_MALLOC(this_slice->next_global_num_last, n_ranks, fvmc_gnum_t);

    for (i = 0; i < n_ranks; i++) {
      this_slice->next_global_num[i] = 0;
      this_slice->next_global_num_last[i] = 0;
    }

  }
  else {
    this_slice->next_global_num = NULL;
    this_slice->next_global_num_last = NULL;
  }

  this_slice->use_next_global_num = false;

  /* Allocated buffers (blocklengths allocated only if needed for ranks > 0) */

  this_slice->recv_buf_size = 0;
  this_slice->recv_buf = NULL;

  this_slice->blocklengths = NULL;

  BFTC_MALLOC(this_slice->displacements, slice_size + 1,  fvmc_gnum_t);

  return this_slice;
}

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
fvmc_gather_slice_destroy(fvmc_gather_slice_t  * this_slice)
{
  if (this_slice != NULL) {

    if (this_slice->next_global_num != NULL) {
      BFTC_FREE(this_slice->next_global_num);
      BFTC_FREE(this_slice->next_global_num_last);
    }

    if (this_slice->recv_buf != NULL)
      BFTC_FREE(this_slice->recv_buf);

    if (this_slice->blocklengths != NULL)
      BFTC_FREE(this_slice->blocklengths);

    BFTC_FREE(this_slice->displacements);

    BFTC_FREE(this_slice);

  }

  return NULL;
}

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
                         fvmc_gnum_t          *global_num_end)
{
  int retval = 0;

  if (this_slice != NULL) {

    if (this_slice->global_num_slice_end > this_slice->global_num_final)
      retval = 1;

    this_slice->global_num_slice_start
      = this_slice->global_num_slice_end;

    this_slice->global_num_slice_end
      = this_slice->global_num_slice_start + this_slice->ref_slice_size;

    if (  this_slice->global_num_slice_end
        > this_slice->global_num_final + 1)
      this_slice->global_num_slice_end = this_slice->global_num_final + 1;

    this_slice->local_index_start = this_slice->local_index_last;

    if (this_slice->next_global_num != NULL) {
      int i;
      for (i = 0; i < this_slice->n_ranks; i++)
        this_slice->next_global_num[i]
          = this_slice->next_global_num_last[i];
    }

    if (   this_slice->global_num_slice_start
        != this_slice->global_num_initial)
      this_slice->use_next_global_num = true;

    *global_num_start = this_slice->global_num_slice_start;
    *global_num_end = this_slice->global_num_slice_end;

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Reset a fvmc_gather_slice_t structure to its initial state.
 *
 * parameters:
 *   this_slice <-- pointer to structure that should be reinitialized
 *----------------------------------------------------------------------------*/

void
fvmc_gather_slice_reinitialize(fvmc_gather_slice_t  *this_slice)
{
  if (this_slice != NULL) {

    this_slice->global_num_slice_start = this_slice->global_num_initial;
    this_slice->global_num_slice_end = this_slice->global_num_initial;

    this_slice->local_index_start = 0;
    this_slice->local_index_last = 0;

    if (this_slice->next_global_num != NULL) {
      int i;
      for (i = 0; i < this_slice->n_ranks; i++) {
        this_slice->next_global_num[i] = 0;
        this_slice->next_global_num_last[i] = 0;
      }
    }

    this_slice->use_next_global_num = false;
  }
}

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
                       fvmc_gnum_t          *global_num_end)
{
  if (*global_num_end != this_slice->global_num_slice_end) {

    if (*global_num_end > this_slice->global_num_final)
      *global_num_end = this_slice->global_num_final;

    this_slice->global_num_slice_end = *global_num_end;

    /* If slice size is changed, the next_global_num[] array of rank 0
       may not be up to date in certain cases, so do not use it */

    this_slice->use_next_global_num = false;

  }
}

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
                       fvmc_gather_slice_t  *this_slice)
{
  int  i, j;
  int  n_local_entities, n_distant_entities;
  fvmc_lnum_t  local_index_start, local_index_stop;

  /* MPI related variables */
  MPI_Status status;
  int  distant_rank;
  const int local_rank = this_slice->local_rank;
  const int n_ranks = this_slice->n_ranks;
  fvmc_gnum_t *const displacements = this_slice->displacements;

  /* Local number of elements */
  const fvmc_lnum_t  local_index_end = this_slice->local_index_end;

  /* Global numbering */
  const fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;
  const fvmc_gnum_t global_num_end = this_slice->global_num_slice_end;
  const fvmc_gnum_t *entity_global_num
                      = fvmc_io_num_get_global_num(element_io_num);

  /* Initialize blocklengths and displacements */

  local_index_start = this_slice->local_index_start;

  assert(   local_index_start >= local_index_end
         || entity_global_num[local_index_start] >= global_num_start);

  for (i = 0, j = local_index_start;
       i < local_index_end && j < local_index_end
                           && entity_global_num[j] < global_num_end;
       i++, j++) {
    displacements[i] = entity_global_num[j] - global_num_start;
  }

  n_local_entities = i;
  local_index_stop = local_index_start + n_local_entities;
  this_slice->local_index_last = local_index_stop;

  if (local_index_stop < local_index_end)
    displacements[n_local_entities] = entity_global_num[j];
  else
    displacements[n_local_entities] = this_slice->global_num_final + 1;

  /*
    Prepare send buffer:
    For rank 0, set "final" values directly;
    values are lengths and not indexes at this stage.
  */

  if (local_rank == 0) {
    for (i = 0, j = local_index_start;
         i < n_local_entities;
         i++, j++) {
      slice_index[displacements[i]] = local_index[j+1] - local_index[j];
    }
  }
  else {
    int n_local_values = n_local_entities;
    for (i = 0, j = local_index_start;
         i < n_local_values;
         i++, j++) {
      slice_index[i] = local_index[j+1] - local_index[j];
    }
  }

  /* Gather lengths information from ranks > 0 */

  if (local_rank == 0) {

    for (distant_rank = 1; distant_rank < n_ranks; distant_rank++) {

      /* Get index from distant rank */

      if (   this_slice->next_global_num[distant_rank] < global_num_end
          || this_slice->use_next_global_num == false) {

        MPI_Send(&distant_rank, 1, MPI_INT, distant_rank, FVMC_MPI_TAG, comm);
        MPI_Recv(&n_distant_entities, 1, MPI_INT,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        MPI_Recv(displacements, n_distant_entities, FVMC_MPI_GNUM,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        n_distant_entities -= 1;
        this_slice->next_global_num_last[distant_rank]
          = displacements[n_distant_entities];

        if (n_distant_entities > 0) {

          fvmc_gnum_t *recv_buf = NULL;

          _slice_recv_buf_size_array(this_slice, n_distant_entities,
                                     1, sizeof(fvmc_gnum_t));

          recv_buf = (fvmc_gnum_t *) this_slice->recv_buf;

          MPI_Recv(this_slice->recv_buf, n_distant_entities,
                   FVMC_MPI_GNUM, distant_rank, FVMC_MPI_TAG, comm, &status);

          for (i = 0 ; i < n_distant_entities ; i++)
            slice_index[displacements[i]] = recv_buf[i];

        }

      }

    }

    /* Transform lengths to global index (0 to n-1)*/

    {
      int l_cur;
      int l_sum = 0;
      int n_entities_max = (int)(global_num_end - global_num_start);
      for (i = 0; i < n_entities_max; i++) {
        l_cur = slice_index[i];
        slice_index[i] = l_sum;
        l_sum += l_cur;
      }
      slice_index[n_entities_max] = l_sum;
    }

  }
  else {

    if (   n_local_entities > 0
        || this_slice->use_next_global_num == false) {

      /* Send local index to rank zero */

      int buf_val;

      MPI_Recv(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm, &status);

      buf_val = n_local_entities + 1;
      MPI_Send(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm);

      MPI_Send(displacements, n_local_entities + 1,
               FVMC_MPI_GNUM, 0, FVMC_MPI_TAG, comm);

      /* Send local portion of lengths array to rank zero */

      if (n_local_entities > 0)
        MPI_Send(slice_index, (int)n_local_entities,
                 FVMC_MPI_GNUM, 0, FVMC_MPI_TAG, comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

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
 * min(n_g_elements, n_elements_s_min) if initially larger than this.
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
                                fvmc_gather_slice_t  *this_slice)
{
  fvmc_gnum_t  i_s;
  fvmc_gnum_t  buf[2];

  fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;

  const int local_rank = this_slice->local_rank;

  *global_num_end = this_slice->global_num_slice_end;

  /* Recompute maximum value of global_num_end for this slice */

  if (local_rank == 0) {

    for (i_s = 1; i_s < *global_num_end - global_num_start + 1; i_s++) {
      if (slice_index[i_s] > *global_connect_s_size) {
        *global_num_end = global_num_start + i_s - 1;
        break;
      }
    }

    /* If possible, size to at least n_elements_s_min (unless
       reference slice size is smaller or already at end of slice) */

    if (*global_num_end - global_num_start < n_elements_s_min) {

      *global_num_end = global_num_start + n_elements_s_min;

      if (*global_num_end - global_num_start > this_slice->ref_slice_size)
        *global_num_end = global_num_start + this_slice->ref_slice_size;

      if (*global_num_end > this_slice->global_num_final + 1)
        *global_num_end = this_slice->global_num_final + 1;

      /* Limit to initial global_num_end */

      if (*global_num_end > this_slice->global_num_slice_end)
        *global_num_end = this_slice->global_num_slice_end;

      *global_connect_s_size =
        FVMC_MAX(*global_connect_s_size,
                slice_index[*global_num_end - global_num_start]);

    }

    buf[0] = *global_num_end, buf[1] = *global_connect_s_size;

  }

  MPI_Bcast(buf, 2, FVMC_MPI_GNUM, 0, comm);

  /* Modify (reduce) slice size limit if necessary */

  fvmc_gather_slice_limit(this_slice, &(buf[0]));

  /* Set output values */

  *global_num_end = buf[0];
  *global_connect_s_size = buf[1];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

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
                 fvmc_gather_slice_t  *this_slice)
{
  int  n_local_entities, n_distant_entities;
  size_t  i, j, k;
  size_t  size_mult;
  fvmc_lnum_t  local_index_start, local_index_stop;

  const char *local_array_val = (const char *) local_array;
  char *global_array_s_val = (char *) global_array_s;

  /* MPI related variables */
  MPI_Status  status;
  int  size, distant_rank;
  const int local_rank = this_slice->local_rank;
  const int n_ranks = this_slice->n_ranks;
  fvmc_gnum_t *const displacements = this_slice->displacements;

  /* Local number of elements */
  const fvmc_lnum_t  local_index_end = this_slice->local_index_end;

  /* Global numbering */
  const fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;
  const fvmc_gnum_t global_num_end = this_slice->global_num_slice_end;
  const fvmc_gnum_t *entity_global_num
                      = fvmc_io_num_get_global_num(element_io_num);

  /* Get info on the current MPI communicator */

  MPI_Type_size(datatype, &size);

  /* Initialize displacements */

  local_index_start = this_slice->local_index_start;

  assert(   local_index_start >= local_index_end
         || entity_global_num[local_index_start] >= global_num_start);

  /* Displacements should be expressed in bytes */

  size_mult = size * stride;

  for (i = 0, j = local_index_start;
       (   i < (size_t)local_index_end
        && j < (size_t)local_index_end
        && entity_global_num[j] < global_num_end);
       i++, j++) {
    displacements[i] = (entity_global_num[j] - global_num_start) * size_mult;
  }

  n_local_entities = i;
  local_index_stop = local_index_start + n_local_entities;
  this_slice->local_index_last = local_index_stop;

  if (local_index_stop < local_index_end)
    displacements[n_local_entities] = entity_global_num[j];
  else
    displacements[n_local_entities] = this_slice->global_num_final + 1;

  /*
    Prepare send buffer (we use a copy to ensure constedness of input)
    For rank 0, set final values directly.
  */

  if (local_rank == 0) {
    for (i = 0, j = (size_t)local_index_start;
         i < (size_t)n_local_entities;
         i++, j++) {
      for (k = 0; k < size_mult; k++) {
        global_array_s_val[displacements[i] + k]
          = local_array_val[j*size_mult + k];
      }
    }
  }
  else
    memcpy(global_array_s_val,
           local_array_val+(local_index_start*size_mult),
           n_local_entities*size_mult);

  /* Gather connectivity information from ranks > 0 */

  if (local_rank == 0) {

    for (distant_rank = 1; distant_rank < n_ranks; distant_rank++) {

      /* Get index from distant rank */

      if (   this_slice->next_global_num[distant_rank] < global_num_end
          || this_slice->use_next_global_num == false) {

        MPI_Send(&distant_rank, 1, MPI_INT, distant_rank, FVMC_MPI_TAG, comm);
        MPI_Recv(&n_distant_entities, 1, MPI_INT,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        MPI_Recv(displacements, n_distant_entities, FVMC_MPI_GNUM,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        n_distant_entities -= 1;
        this_slice->next_global_num_last[distant_rank]
          = displacements[n_distant_entities];

        if (n_distant_entities > 0) {

          char *_recv_buf = NULL;

          _slice_recv_buf_size_array(this_slice, n_distant_entities,
                                     stride, size);

          _recv_buf = (char *) this_slice->recv_buf;

          MPI_Recv(this_slice->recv_buf, n_distant_entities*stride, datatype,
                   distant_rank, FVMC_MPI_TAG, comm, &status);

          for (i = 0 ; i < (size_t)n_distant_entities ; i++) {
            for (j = 0 ; j < size_mult ; j++) {
              global_array_s_val[displacements[i] + j]
                = _recv_buf[i*size_mult + j];
            }
          }

        }

      }

    }

  }
  else {

    if (   n_local_entities > 0
        || this_slice->use_next_global_num ==false) {

      /* Send local index to rank zero */

      int buf_val;

      MPI_Recv(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm, &status);

      buf_val = n_local_entities + 1;
      MPI_Send(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm);

      MPI_Send(displacements, n_local_entities + 1, FVMC_MPI_GNUM, 0,
               FVMC_MPI_TAG, comm);

      /* Send local portion of connectivity array to rank zero */

      if (n_local_entities > 0)
        MPI_Send(global_array_s, (int)(n_local_entities * stride),
                 datatype, 0, FVMC_MPI_TAG, comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

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
                   fvmc_gather_slice_t  *this_slice)
{
  int  i, j, k, l;
  int  n_local_entities, n_distant_entities;
  int  n_values_send = 0;
  fvmc_lnum_t  local_index_start, local_index_stop;

  const char *local_array_val = (const char *) local_array;
  char *global_array_s_val = (char *) global_array_s;

  /* MPI related variables */
  MPI_Status status;
  int  size, distant_rank;
  const int local_rank = this_slice->local_rank;
  const int n_ranks = this_slice->n_ranks;
  int  *blocklengths = this_slice->blocklengths;
  fvmc_gnum_t *const displacements = this_slice->displacements;

  /* Local number of elements */
  const fvmc_lnum_t  local_index_end = this_slice->local_index_end;

  /* Global numbering */
  const fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;
  const fvmc_gnum_t global_num_end = this_slice->global_num_slice_end;
  const fvmc_gnum_t *entity_global_num
                      = fvmc_io_num_get_global_num(element_io_num);

  MPI_Type_size(datatype, &size);

  /* Initialize blocklengths and displacements; this operation
     requires a blocklengths[] array for all ranks (contrary
     to others), so it updates the slice structure accordingly */

  if (blocklengths == NULL) {
    BFTC_MALLOC(this_slice->blocklengths, this_slice->ref_slice_size, int);
    blocklengths = this_slice->blocklengths;
  }

  local_index_start = this_slice->local_index_start;

  assert(   local_index_start >= local_index_end
         || entity_global_num[local_index_start] >= global_num_start);

  /* Displacements are first used to transfer the global slice index position
     of a given entity, and will be set to the index value later on rank 0 */

  for (i = 0, j = local_index_start;
       j < local_index_end && entity_global_num[j] < global_num_end;
       i++, j++)
    displacements[i] = entity_global_num[j] - global_num_start;

  n_local_entities = i;
  local_index_stop = local_index_start + n_local_entities;
  this_slice->local_index_last = local_index_stop;

  if (local_index_stop < local_index_end)
    displacements[n_local_entities] = entity_global_num[j];
  else
    displacements[n_local_entities] = this_slice->global_num_final + 1;

  /*
    Prepare send buffer:
    For rank 0, set final values directly; for others, set blocklengths
    to be sent to rank 0.
  */

  if (local_rank == 0) {
    for (i = 0, j = local_index_start;
         i < n_local_entities;
         i++, j++) {
      const size_t displacement = slice_index[displacements[i]] * size;
      for (k = local_index[j]*(int)size, l = 0;
           k < (local_index[j+1])*(int)size;
           k++, l++)
        global_array_s_val[displacement + l] = local_array_val[k];
    }
  }
  else {
    n_values_send = (  local_index[local_index_start + n_local_entities]
                     - local_index[local_index_start]);
    memcpy(global_array_s_val,
           local_array_val + ((local_index[local_index_start]) * size),
           n_values_send * size);
    for (i = 0, j = local_index_start;
         i < n_local_entities;
         i++, j++) {
      blocklengths[i] = (local_index[j+1] - local_index[j]);
    }
  }

  /* Gather numbers information from ranks > 0 */

  if (local_rank == 0) {

    for (distant_rank = 1; distant_rank < n_ranks; distant_rank++) {

      size_t  recv_size = 0;

      /* Get index from distant rank */

      if (   this_slice->next_global_num[distant_rank] < global_num_end
          || this_slice->use_next_global_num == false) {

        MPI_Send(&distant_rank, 1, MPI_INT, distant_rank, FVMC_MPI_TAG, comm);
        MPI_Recv(&n_distant_entities, 1, MPI_INT,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        /* Get index from distant rank */

        MPI_Recv(displacements, n_distant_entities, FVMC_MPI_GNUM,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        n_distant_entities -= 1;
        this_slice->next_global_num_last[distant_rank]
          = displacements[n_distant_entities];

        for (i = 0; i < n_distant_entities; i++) {
          j = displacements[i];
          blocklengths[i] = (slice_index[j+1] - slice_index[j])*size;
          displacements[i] = slice_index[j] * size;
          recv_size += blocklengths[i];
        }

        if (n_distant_entities > 0) {

          size_t  recv_id;
          char  *_recv_buf = NULL;

          _slice_recv_buf_size_indexed(this_slice,
                                       recv_size,
                                       size);

          MPI_Recv(this_slice->recv_buf, (int)recv_size, datatype,
                   distant_rank, FVMC_MPI_TAG, comm, &status);

          _recv_buf = (char *) this_slice->recv_buf;

          for (i = 0, recv_id = 0 ; i < n_distant_entities ; i++) {
            for (j = 0 ; j < blocklengths[i] ; j++) {
              global_array_s_val[displacements[i] + j]
                = _recv_buf[recv_id++];
            }
          }

        }

      }

    }

  }
  else {

    if (   n_local_entities > 0
        || this_slice->use_next_global_num == false) {

      /* Send local index to rank zero */

      int buf_val;

      MPI_Recv(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm, &status);

      buf_val = n_local_entities + 1;
      MPI_Send(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm);

      MPI_Send(displacements, n_local_entities + 1, FVMC_MPI_GNUM,
               0, FVMC_MPI_TAG, comm);

      /* Send local portion of numbers array to rank zero */

      if (n_local_entities > 0)
        MPI_Send(global_array_s, (int)n_values_send,
                 datatype, 0, FVMC_MPI_TAG, comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

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
fvmc_gather_strided_connect(const fvmc_lnum_t    local_connect[],
                           fvmc_gnum_t          global_connect_s[],
                           const int           stride,
                           const fvmc_io_num_t  *connected_io_num,
                           const fvmc_io_num_t  *element_io_num,
                           MPI_Comm             comm,
                           fvmc_gather_slice_t  *this_slice)
{
  int  i, j, k;
  int  n_local_entities, n_distant_entities;
  fvmc_lnum_t  local_index_start, local_index_stop;

  /* MPI related variables */
  MPI_Status status;
  int  distant_rank;
  const int local_rank = this_slice->local_rank;
  const int n_ranks = this_slice->n_ranks;
  fvmc_gnum_t *const displacements = this_slice->displacements;

  /* Local number of elements */
  const fvmc_lnum_t  local_index_end = this_slice->local_index_end;

  /* Global numbering */
  const fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;
  const fvmc_gnum_t global_num_end = this_slice->global_num_slice_end;
  const fvmc_gnum_t *connected_global_num
                      = fvmc_io_num_get_global_num(connected_io_num);
  const fvmc_gnum_t *entity_global_num
                      = fvmc_io_num_get_global_num(element_io_num);

  /* Initialize displacements */

  local_index_start = this_slice->local_index_start;

  assert(   local_index_start >= local_index_end
         || entity_global_num[local_index_start] >= global_num_start);

  for (i = 0, j = local_index_start;
       i < local_index_end && j < local_index_end
                           && entity_global_num[j] < global_num_end;
       i++, j++) {
    displacements[i] = (entity_global_num[j] - global_num_start) * stride;
  }

  n_local_entities = i;
  local_index_stop = local_index_start + n_local_entities;
  this_slice->local_index_last = local_index_stop;

  if (local_index_stop < local_index_end)
    displacements[n_local_entities] = entity_global_num[j];
  else
    displacements[n_local_entities] = this_slice->global_num_final + 1;

  /*
    Prepare send buffer:
    replace local connected entity numbers with their global counterparts.
    For rank 0, set final values directly.
  */

  if (local_rank == 0) {
    for (i = 0, j = local_index_start;
         i < n_local_entities;
         i++, j++) {
      for (k = 0; k < stride; k++)
        global_connect_s[displacements[i] + k]
          = connected_global_num[local_connect[j*stride + k] - 1];
    }
  }
  else {
    int n_local_values = n_local_entities * stride;
    for (i = 0, j = (fvmc_gnum_t)(local_index_start * stride);
         i < n_local_values;
         i++, j++) {
      global_connect_s[i] = connected_global_num[local_connect[j] - 1];
    }
  }

  /* Gather connectivity information from ranks > 0 */

  if (local_rank == 0) {

    for (distant_rank = 1; distant_rank < n_ranks; distant_rank++) {

      /* Get index from distant rank */

      if (   this_slice->next_global_num[distant_rank] < global_num_end
          || this_slice->use_next_global_num == false) {

        MPI_Send(&distant_rank, 1, MPI_INT, distant_rank, FVMC_MPI_TAG, comm);
        MPI_Recv(&n_distant_entities, 1, MPI_INT,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        MPI_Recv(displacements, n_distant_entities, FVMC_MPI_GNUM,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        n_distant_entities -= 1;
        this_slice->next_global_num_last[distant_rank]
          = displacements[n_distant_entities];

        if (n_distant_entities > 0) {

          fvmc_gnum_t *_recv_buf;

          _slice_recv_buf_size_array(this_slice, n_distant_entities,
                                     stride, sizeof(fvmc_gnum_t));

          _recv_buf = (fvmc_gnum_t *) this_slice->recv_buf;

          MPI_Recv(this_slice->recv_buf, (int)(n_distant_entities*stride),
                   FVMC_MPI_GNUM, distant_rank, FVMC_MPI_TAG, comm, &status);

          for (i = 0 ; i < n_distant_entities ; i++) {
            for (j = 0 ; j < stride ; j++)
              global_connect_s[displacements[i] + j]
                = _recv_buf[i*stride + j];
          }

        }

      }

    }

  }
  else {

    if (   n_local_entities > 0
        || this_slice->use_next_global_num == false) {

      /* Send local index to rank zero */

      int buf_val;

      MPI_Recv(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm, &status);

      buf_val = n_local_entities + 1;
      MPI_Send(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm);

      MPI_Send(displacements, n_local_entities + 1, FVMC_MPI_GNUM, 0,
               FVMC_MPI_TAG, comm);

      /* Send local portion of connectivity array to rank zero */

      if (n_local_entities > 0)
        MPI_Send(global_connect_s, (int)(n_local_entities * stride),
                 FVMC_MPI_GNUM, 0, FVMC_MPI_TAG, comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

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
 * enough to contain (slice_index[current_slice_size] - 1) values, even
 * though most processes will require less.
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
                           fvmc_gather_slice_t  *this_slice)
{
  int  i, j, k, l;
  int  n_local_entities, n_distant_entities;
  int  n_values_send = 0;
  fvmc_lnum_t  local_index_start, local_index_stop;

  /* MPI related variables */
  MPI_Status status;
  int  distant_rank;
  const int local_rank = this_slice->local_rank;
  const int n_ranks = this_slice->n_ranks;
  int  *blocklengths = this_slice->blocklengths;
  fvmc_gnum_t *const displacements = this_slice->displacements;

  /* Local number of elements */
  const fvmc_lnum_t  local_index_end = this_slice->local_index_end;

  /* Global numbering */
  const fvmc_gnum_t global_num_start = this_slice->global_num_slice_start;
  const fvmc_gnum_t global_num_end = this_slice->global_num_slice_end;
  const fvmc_gnum_t *connected_global_num = NULL;
  const fvmc_gnum_t *entity_global_num
                      = fvmc_io_num_get_global_num(element_io_num);

  if (connected_io_num != NULL)
    connected_global_num = fvmc_io_num_get_global_num(connected_io_num);

  /* Initialize blocklengths and displacements; this operation
     requires a blocklengths[] array for all ranks (contrary
     to others), so it updates the slice structure accordingly */

  if (blocklengths == NULL) {
    BFTC_MALLOC(this_slice->blocklengths, this_slice->ref_slice_size, int);
    blocklengths = this_slice->blocklengths;
  }

  local_index_start = this_slice->local_index_start;

  assert(   local_index_start >= local_index_end
         || entity_global_num[local_index_start] >= global_num_start);

  /* Displacements are first used to transfer the global slice index position
     of a given entity, and will be set to the index value later on rank 0 */

  for (i = 0, j = local_index_start;
       i < local_index_end && j < local_index_end
                           && entity_global_num[j] < global_num_end;
       i++, j++)
    displacements[i] = entity_global_num[j] - global_num_start;

  n_local_entities = i;
  local_index_stop = local_index_start + n_local_entities;
  this_slice->local_index_last = local_index_stop;

  if (local_index_stop < local_index_end)
    displacements[n_local_entities] = entity_global_num[j];
  else
    displacements[n_local_entities] = this_slice->global_num_final + 1;

  /*
    Prepare send buffer:
    For rank 0, set final values directly; for others, set blocklengths
    to be sent to rank 0.
  */

  if (connected_io_num == NULL) {

    if (local_rank == 0) {
      for (i = 0, j = local_index_start;
           i < n_local_entities;
           i++, j++) {
        displacements[i] = slice_index[displacements[i]];
        for (k = local_index[j], l = 0; k < local_index[j+1]; k++, l++)
          global_numbers_s[displacements[i] + l] = local_numbers[k];
      }
    }
    else {
      l = 0;
      for (i = 0, j = local_index_start;
           i < n_local_entities;
           i++, j++) {
        blocklengths[i] = local_index[j+1] - local_index[j];
        for (k = local_index[j]; k < local_index[j+1]; k++)
          global_numbers_s[l++] = local_numbers[k];
      }
      n_values_send = l;
    }

  }
  else {

    if (local_rank == 0) {
      for (i = 0, j = local_index_start;
           i < n_local_entities;
           i++, j++) {
        displacements[i] = slice_index[displacements[i]];
        for (k = local_index[j], l = 0; k < local_index[j+1]; k++, l++)
          global_numbers_s[displacements[i] + l]
            = connected_global_num[local_numbers[k] - 1];
      }
    }
    else {
      l = 0;
      for (i = 0, j = local_index_start;
           i < n_local_entities;
           i++, j++) {
        blocklengths[i] = local_index[j+1] - local_index[j];
        for (k = local_index[j]; k < local_index[j+1]; k++)
          global_numbers_s[l++] = connected_global_num[local_numbers[k] - 1];
      }
      n_values_send = l;
    }

  }

  /* Gather numbers information from ranks > 0 */

  if (local_rank == 0) {

    for (distant_rank = 1; distant_rank < n_ranks; distant_rank++) {

      /* Get index from distant rank */

      if (   this_slice->next_global_num[distant_rank] < global_num_end
          || this_slice->use_next_global_num == false) {

        size_t  recv_size = 0;

        MPI_Send(&distant_rank, 1, MPI_INT, distant_rank, FVMC_MPI_TAG, comm);
        MPI_Recv(&n_distant_entities, 1, MPI_INT,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

      /* Get index from distant rank */

        MPI_Recv(displacements, n_distant_entities, FVMC_MPI_GNUM,
                 distant_rank, FVMC_MPI_TAG, comm, &status);

        n_distant_entities -= 1;
        this_slice->next_global_num_last[distant_rank]
          = displacements[n_distant_entities];

        for (i = 0; i < n_distant_entities; i++) {
          j = displacements[i];
          blocklengths[i] = slice_index[j+1] - slice_index[j];
          displacements[i] = slice_index[j];
          recv_size += blocklengths[i];
        }

        if (n_distant_entities > 0) {

          size_t  recv_id;
          fvmc_gnum_t  *_recv_buf = NULL;

          _slice_recv_buf_size_indexed(this_slice,
                                       recv_size,
                                       sizeof(fvmc_gnum_t));

          MPI_Recv(this_slice->recv_buf, recv_size, FVMC_MPI_GNUM,
                   distant_rank, FVMC_MPI_TAG, comm, &status);

          _recv_buf = (fvmc_gnum_t *) this_slice->recv_buf;

          for (i = 0, recv_id = 0 ; i < n_distant_entities ; i++) {
            for (j = 0 ; j < blocklengths[i] ; j++)
              global_numbers_s[displacements[i] + j]
                = _recv_buf[recv_id++];
          }

        }

      }

    }

  }
  else {

    if (   n_local_entities > 0
        || this_slice->use_next_global_num == false) {

      /* Send local index to rank zero */

      int buf_val;

      MPI_Recv(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm, &status);

      buf_val = n_local_entities + 1;
      MPI_Send(&buf_val, 1, MPI_INT, 0, FVMC_MPI_TAG, comm);

      MPI_Send(displacements, n_local_entities + 1, FVMC_MPI_GNUM,
               0, FVMC_MPI_TAG, comm);

      /* Send local portion of numbers array to rank zero */

      if (n_local_entities > 0)
        MPI_Send(global_numbers_s, (int)n_values_send,
                 FVMC_MPI_GNUM, 0, FVMC_MPI_TAG, comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
