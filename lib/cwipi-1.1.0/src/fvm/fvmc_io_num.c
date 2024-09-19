/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
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
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_config_defs.h"
#include "fvmc_morton.h"
#include "fvmc_order.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_io_num.h"

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

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/*
 * Notes:
 *
 * This structure currently only contains a global numbering array containing
 * each entity's global number (a 1 to n index). In the future, it may
 * also contain information relative to ghost zones, for I/O to file formats
 * enabling domain splitting with (with multiple groups of processes writing
 * different subsets, using ghost zones for entities on boundaries between
 * mesh parts assigned to different process groups). In such a case, the
 * main numbering would only be global as regards a process group, at least
 * for use with formats such as that of EnSight Gold.
 *
 * In some cases, a global number may appear on multiple processes (for a face
 * or vertex on a processor boundary). This means that multiple processes may
 * update the corresponding value in a gather-type operation preceding or
 * associated with I/O, in an undefined order. If the associated data is up to
 * date on each process, it should be identical, so this should not be a
 * problem. MPI-IO file writes or MPI-2 one-sided communication PUT operations
 * also authorize this, though for the latter, the MPI-2 standard indicates
 * that this is authorized if processes overlapping PUT operations should use
 * the same predefined datatype, which seems to exclude similar indexed
 * datatypes with different indexes. To avoid problems if we wish to use
 * MPI-2 one-sided communication, one relatively simple solution would be to
 * consider that for processes other than that of lowest rank in which an
 * entity appears, it appears after all entities not occuring in processes
 * of lower rank. In that case, we would have two array sizes:
 * global_num_size_tot defining the array's full size, and global_num_size
 * defining the size of the portion to use in gather type operations.
 */

struct _fvmc_io_num_t {

  fvmc_gnum_t         global_count;    /* Global number of entities */
  fvmc_lnum_t         global_num_size; /* Local size of global numbering array */
  const fvmc_gnum_t  *global_num;      /* Global (possibly shared) entity
                                         numbers (1 to n) */
  fvmc_gnum_t        *_global_num;     /* Global entity numbers if owner,
                                         NULL otherwise */

};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Use bubble sort on an expectedly short sequence of coordinates
 * to ensure lexicographical ordering.
 *
 * parameters:
 *   dim        <-- spatial dimension
 *   start_id   <-- start id in array
 *   end_id     <-- past-the-end id in array
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   order      <-> ordering array base on Morton encoding, or
 *                  lexicographical coordinate ordering for ties
 *----------------------------------------------------------------------------*/

inline static void
_reorder_coords_lexicographic(int                dim,
                              size_t             start_id,
                              size_t             end_id,
                              const fvmc_coord_t  coords[],
                              fvmc_lnum_t         order[])
{
  size_t  i;
  _Bool g_swap;

  do {

    g_swap = false;

    for (i = start_id + 1; i < end_id; i++) {

      size_t j_prev = order[i-1], j = order[i];
      _Bool l_swap = false;

      if (dim == 3) {
        if (coords[j_prev*3] < coords[j*3])
          continue;
        else if (coords[j_prev*3] > coords[j*3])
          l_swap = true;
        else if (coords[j_prev*3 + 1] < coords[j*3 + 1])
          continue;
        else if (   coords[j_prev*3 + 1] > coords[j*3 + 1]
                 || coords[j_prev*3 + 2] > coords[j*3 + 2])
          l_swap = true;
      }
      else if (dim == 2) {
        if (coords[j_prev*2] < coords[j*2 + 1])
          continue;
        else if (   coords[j_prev*2]     > coords[j*2]
                 || coords[j_prev*2 + 1] > coords[j*2 + 1])
          l_swap = true;
      }
      else { /* if (dim == 1) */
        if (coords[j_prev] > coords[j])
          l_swap = true;
      }

      if (l_swap) {
        fvmc_lnum_t o_save = order[i-1];
        order[i-1] = order[i];
        order[i] = o_save;
        g_swap = true;
      }
    }

  } while (g_swap);
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be different, but their order is undetermined.
 *
 * parameters:
 *   dim        <-- spatial dimension
 *   n_entities <-- number of entities considered
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   m_code     <-- Morton code associated with each entity
 *   order      <-> ordering array base on Morton encoding, or
 *                  lexicographical coordinate ordering for ties
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

static void
_check_morton_ordering(int                      dim,
                       size_t                   n_entities,
                       const fvmc_coord_t        coords[],
                       const fvmc_morton_code_t  m_code[],
                       fvmc_lnum_t               order[])
{
  size_t  i_prev = 0, i = 1;

  if (n_entities == 0)
    return;

  /* Check ordering; if two entities have the same Morton codes,
     use lexicographical coordinates ordering to ensure the
     final order is deterministic. */

  for (i = 1; i < n_entities; i++) {

    size_t j_prev = order[i_prev], j = order[i];

    if (   m_code[j_prev].X[0] != m_code[j].X[0]
        || m_code[j_prev].X[1] != m_code[j].X[1]
        || m_code[j_prev].X[2] != m_code[j].X[2]) {

      /* If successive values have the same Morton code,
         order them lexicographically */
      if (i_prev < i - 1)
        _reorder_coords_lexicographic(dim, i_prev, i-1, coords, order);

    }
    i_prev = i;
  }

  if (i_prev < n_entities - 1)
    _reorder_coords_lexicographic(dim, i_prev, n_entities - 1, coords, order);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Copy selected shared global ordering information to private ordering
 * information for an I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to numbering structure
 *----------------------------------------------------------------------------*/

static void
_fvmc_io_num_copy_on_write(fvmc_io_num_t  *const this_io_num)
{
  if (this_io_num->_global_num == NULL) {
    fvmc_lnum_t i;
    BFTC_MALLOC(this_io_num->_global_num,
               this_io_num->global_num_size,
               fvmc_gnum_t);
    for (i = 0; i < this_io_num->global_num_size; i++)
      this_io_num->_global_num[i] = this_io_num->global_num[i];
    this_io_num->global_num = this_io_num->_global_num;
  }
  assert(this_io_num->global_num == this_io_num->_global_num);
}

/*----------------------------------------------------------------------------
 * Copy selected shared global ordering information to private ordering
 * information for an I/O numbering structure.
 *
 * parameters:
 *   this_io_num          <-> pointer to numbering structure
 *   parent_global_number <-- pointer to shared list of global parent
 *                            entity numbers
 *----------------------------------------------------------------------------*/

static void
_fvmc_io_num_try_to_set_shared(fvmc_io_num_t      *const this_io_num,
                              const fvmc_gnum_t         parent_global_number[])
{
  if (this_io_num->_global_num != NULL && parent_global_number != NULL) {
    fvmc_lnum_t i;
    for (i = 0; i < this_io_num->global_num_size; i++)
      if (this_io_num->_global_num[i] != parent_global_number[i])
        break;
    if (i < this_io_num->global_num_size)
      this_io_num->global_num = this_io_num->_global_num;
    else {
      this_io_num->global_num = parent_global_number;
      BFTC_FREE(this_io_num->_global_num);
    }
  }
}

/*----------------------------------------------------------------------------
 * Maximum global number associated with an I/O numbering structure
 *
 * parameters:
 *   this_io_num <-- pointer to partially initialized I/O numbering structure.
 *   comm        <-- associated MPI communicator
 *
 * returns:
 *   maximum global number associated with the I/O numbering
 *----------------------------------------------------------------------------*/

static fvmc_gnum_t
_fvmc_io_num_global_max(const fvmc_io_num_t  *const this_io_num,
                       const MPI_Comm             comm)
{
  fvmc_gnum_t  local_max, global_max;
  size_t      n_ent;

  /* Get maximum global number value */

  n_ent = this_io_num->global_num_size;
  if (n_ent > 0)
    local_max = this_io_num->global_num[n_ent - 1];
  else
    local_max = 0;

  MPI_Allreduce(&local_max, &global_max, 1, FVMC_MPI_GNUM, MPI_MAX, comm);

  return global_max;
}

/*----------------------------------------------------------------------------
 * Global ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * As an option, a number of sub-entities per initial entity may be
 * given, in which case sub-entities of a same entity will have contiguous
 * numbers in the final ordering.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   n_sub_entities <-- optional number of sub-entities per initial entity,
 *                      or NULL if unused
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvmc_io_num_global_order(fvmc_io_num_t       *this_io_num,
                         const fvmc_lnum_t    n_sub_entities[],
                         MPI_Comm            comm)
{

  fvmc_gnum_t  n_ent_recv, num_prev, num_cur;
  size_t      i, j, slice_size;
  fvmc_lnum_t  k;
  int         rank;

  _Bool       may_be_shared = false;

  fvmc_gnum_t  *recv_global_num = NULL;
  fvmc_lnum_t  *recv_n_sub = NULL, *recv_order = NULL;
  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;
  int         have_sub_loc = 0, have_sub_glob = 0;

  int         local_rank, size;

  fvmc_gnum_t  current_global_num = 0, global_num_shift = 0;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  num_prev = 0;    /* true initialization later for slice 0, */

  /* If numbering is shared, we will check later it was changed or if
     it can remain shared (in which case the copy may be discarded) */

  if (this_io_num->global_num != this_io_num->_global_num)
    may_be_shared = true;
  else
    may_be_shared = false;

  /* Get temporary maximum global number value */

  this_io_num->global_count = _fvmc_io_num_global_max(this_io_num, comm);

  /* slice_size = ceil(this_io_num->global_count/size) */

  slice_size = this_io_num->global_count / size;
  if (this_io_num->global_count % size > 0)
    slice_size += 1;

  assert(sizeof(fvmc_gnum_t) >= sizeof(fvmc_lnum_t));

  BFTC_MALLOC(send_count, size, int);
  BFTC_MALLOC(recv_count, size, int);

  BFTC_MALLOC(send_shift, size, int);
  BFTC_MALLOC(recv_shift, size, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
    send_count[(this_io_num->global_num[i] - 1) / slice_size] += 1;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 1; rank < size; rank++) {
    send_shift[rank] = send_shift[rank - 1] + send_count[rank -1];
    recv_shift[rank] = recv_shift[rank - 1] + recv_count[rank -1];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the n_sub_entities[] array to the correct processors */

  n_ent_recv = recv_shift[size - 1] + recv_count[size - 1];

  BFTC_MALLOC(recv_global_num, n_ent_recv, fvmc_gnum_t);
  BFTC_MALLOC(recv_order, n_ent_recv, fvmc_lnum_t);

  MPI_Alltoallv(this_io_num->_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVMC_MPI_GNUM, comm);

  /* Do we have sub-entities ? */

  if (n_sub_entities != NULL)
    have_sub_loc = 1;

  MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX, comm);

  if (have_sub_glob > 0) {

    fvmc_lnum_t  *send_n_sub;

    BFTC_MALLOC(send_n_sub, this_io_num->global_num_size, fvmc_lnum_t);
    BFTC_MALLOC(recv_n_sub, n_ent_recv, fvmc_lnum_t);

    if (n_sub_entities != NULL) {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = n_sub_entities[i];
    }
    else {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = 1;
    }

    MPI_Alltoallv(send_n_sub, send_count, send_shift, FVMC_MPI_LNUM,
                  recv_n_sub, recv_count, recv_shift, FVMC_MPI_LNUM, comm);

    BFTC_FREE(send_n_sub);
  }

  if (n_ent_recv > 0) {

    fvmc_order_local_allocated(NULL,
                              recv_global_num,
                              recv_order,
                              n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    if (have_sub_glob > 0) {

      current_global_num = recv_n_sub[recv_order[0]];
      num_prev = recv_global_num[recv_order[0]];
      recv_global_num[recv_order[0]] = current_global_num;

      for (i = 1; i < n_ent_recv; i++) {
        num_cur = recv_global_num[recv_order[i]];
        if (num_cur > num_prev)
          current_global_num += recv_n_sub[recv_order[i]];
        recv_global_num[recv_order[i]] = current_global_num;
        num_prev = num_cur;
      }

    }
    else { /* if (have_sub_glob == 0) */

      current_global_num = 1;
      num_prev = recv_global_num[recv_order[0]];
      recv_global_num[recv_order[0]] = current_global_num;

      for (i = 1; i < n_ent_recv; i++) {
        num_cur = recv_global_num[recv_order[i]];
        if (num_cur > num_prev)
          current_global_num += 1;
        recv_global_num[recv_order[i]] = current_global_num;
        num_prev = num_cur;
      }

    }

  }

  /* Partial clean-up */

  BFTC_FREE(recv_n_sub);
  BFTC_FREE(recv_order);

  /* At this stage, recv_global_num[] is valid for this process, and
     current_global_num indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_global_num, &global_num_shift, 1, FVMC_MPI_GNUM,
           MPI_SUM, comm);
  global_num_shift -= current_global_num;

  for (i = 0; i < n_ent_recv; i++)
    recv_global_num[i] += global_num_shift;

  /* Return global order to all processors */

  MPI_Alltoallv(recv_global_num, recv_count, recv_shift, FVMC_MPI_GNUM,
                this_io_num->_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                comm);

  /* Free memory */

  BFTC_FREE(recv_global_num);

  BFTC_FREE(send_count);
  BFTC_FREE(recv_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_shift);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvmc_io_num_global_max(this_io_num, comm);

  /* When sub-entities have been added, now switch from a numbering on
     the initial entities (shifted by number of sub-entities) to
     a numbering on the final sub-entities */

  if (n_sub_entities != NULL) {

    fvmc_gnum_t *_global_num;

    for (i = 0, j = 0; i < (size_t)(this_io_num->global_num_size); i++)
      j += n_sub_entities[i];

    BFTC_MALLOC(_global_num, j, fvmc_gnum_t);

    for (i = 0, j = 0; i < (size_t)(this_io_num->global_num_size); i++) {
      for (k = 0; k < n_sub_entities[i]; j++, k++)
        _global_num[j] = this_io_num->_global_num[i] - n_sub_entities[i] + k + 1;
    }

    BFTC_FREE(this_io_num->_global_num);
    this_io_num->_global_num = _global_num;

    if (this_io_num->global_num_size != (fvmc_lnum_t)j) {
      this_io_num->global_num_size = j;
      may_be_shared = false;
    }

    if (may_be_shared == false)
      this_io_num->global_num = this_io_num->_global_num;
  }

  /* If numbering was initially shared, check if it was changed or if it
     may remain shared (in which case the copy may be discarded) */

  if (may_be_shared == true) {
    for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
      if (this_io_num->_global_num[i] != this_io_num->global_num[i])
        break;
    if (i < (size_t)(this_io_num->global_num_size))
      this_io_num->global_num = this_io_num->_global_num;
    else
      BFTC_FREE(this_io_num->_global_num);
  }
}

/*----------------------------------------------------------------------------
 * Global ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   stride         <-- values per entity
 *   global_num     <-- global numbering array
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvmc_io_num_global_order_s(fvmc_io_num_t       *this_io_num,
                           size_t              stride,
                           fvmc_gnum_t          global_num[],
                           MPI_Comm            comm)
{
  fvmc_gnum_t  n_ent_recv;
  size_t      i, j, slice_size;
  int         rank;

  fvmc_gnum_t  *block_global_num = NULL, *recv_global_num = NULL;
  fvmc_lnum_t  *recv_order = NULL;
  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;

  int         local_rank, size;

  fvmc_gnum_t  current_global_num = 0, global_num_shift = 0;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  /* Get maximum global number value for first value of each series
     (does not need to be exact, simply used to define blocks) */

  {
    fvmc_gnum_t  local_max = 0, global_max = 0;
    size_t      n_ent = this_io_num->global_num_size;

    if (n_ent > 0)
      local_max = global_num[(n_ent-1)*stride];
    MPI_Allreduce(&local_max, &global_max, 1, FVMC_MPI_GNUM, MPI_MAX, comm);
    this_io_num->global_count = global_max;
  }

  /* slice_size = ceil(this_io_num->global_count/size) */

  slice_size = this_io_num->global_count / size;
  if (this_io_num->global_count % size > 0)
    slice_size += 1;

  assert(sizeof(fvmc_gnum_t) >= sizeof(fvmc_lnum_t));

  BFTC_MALLOC(send_count, size, int);
  BFTC_MALLOC(recv_count, size, int);

  BFTC_MALLOC(send_shift, size, int);
  BFTC_MALLOC(recv_shift, size, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
    send_count[(global_num[stride*i] - 1) / slice_size] += stride;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 1; rank < size; rank++) {
    send_shift[rank] = send_shift[rank - 1] + send_count[rank -1];
    recv_shift[rank] = recv_shift[rank - 1] + recv_count[rank -1];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the n_sub_entities[] array to the correct processors */

  n_ent_recv = (recv_shift[size - 1] + recv_count[size - 1]) / stride;

  BFTC_MALLOC(recv_global_num, stride*n_ent_recv, fvmc_gnum_t);
  BFTC_MALLOC(recv_order, n_ent_recv, fvmc_lnum_t);

  MPI_Alltoallv(global_num, send_count, send_shift, FVMC_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVMC_MPI_GNUM, comm);

  if (n_ent_recv > 0) {

    size_t prev_id, cur_id;

    fvmc_order_local_allocated_s(NULL,
                                recv_global_num,
                                stride,
                                recv_order,
                                n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    BFTC_MALLOC(block_global_num, n_ent_recv, fvmc_gnum_t);

    current_global_num = 1;
    prev_id = recv_order[0];
    block_global_num[recv_order[0]] = current_global_num;

    for (i = 1; i < n_ent_recv; i++) {
      _Bool greater_than_prev = false;
      cur_id = recv_order[i];
      for (j = 0; j < stride; j++) {
        if (  recv_global_num[cur_id*stride + j]
            > recv_global_num[prev_id*stride + j])
          greater_than_prev = true;
      }
      if (greater_than_prev)
        current_global_num += 1;
      block_global_num[recv_order[i]] = current_global_num;
      prev_id = cur_id;
    }

  }

  /* Partial clean-up */

  BFTC_FREE(recv_order);
  BFTC_FREE(recv_global_num);

  /* At this stage, recv_global_num[] is valid for this process, and
     current_global_num indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_global_num, &global_num_shift, 1, FVMC_MPI_GNUM,
           MPI_SUM, comm);
  global_num_shift -= current_global_num;

  for (i = 0; i < n_ent_recv; i++)
    block_global_num[i] += global_num_shift;

  /* Return global order to all processors */

  for (rank = 0; rank < size; rank++) {
    send_count[rank] /= stride;
    recv_count[rank] /= stride;
  }

  for (rank = 1; rank < size; rank++) {
    send_shift[rank] = send_shift[rank - 1] + send_count[rank -1];
    recv_shift[rank] = recv_shift[rank - 1] + recv_count[rank -1];
  }

  MPI_Alltoallv(block_global_num, recv_count, recv_shift, FVMC_MPI_GNUM,
                this_io_num->_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                comm);

  /* Free memory */

  BFTC_FREE(block_global_num);
  BFTC_FREE(send_count);
  BFTC_FREE(recv_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_shift);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvmc_io_num_global_max(this_io_num, comm);
}

/*----------------------------------------------------------------------------
 * Compare two elements in an indexed list and returns true if element in
 * position i1 is strictly greater than element in position i2.
 *
 * parameters:
 *   i1        <-- position in index for the first element
 *   i2        <-- position in index for the second element
 *   index     <-- number of values to compare for each entity
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_indexed_is_greater(size_t            i1,
                    size_t            i2,
                    const fvmc_lnum_t  index[],
                    const fvmc_gnum_t  number[])
{
  int  i;

  fvmc_lnum_t  i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  fvmc_lnum_t  i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

  if (s1 > s2) {

    for (i = 0; i < s2; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return true;
  }
  else { /* s1 <= s2 */

    for (i = 0; i < s1; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return false;
  }

}

/*----------------------------------------------------------------------------
 * Global indexed ordering associated with an I/O numbering structure.
 *
 * The structure should contain an initial ordering, which should
 * be sorted, but need not be contiguous. On output, the numbering
 * will be contiguous.
 *
 * parameters:
 *   this_io_num    <-> pointer to structure that should be ordered
 *   index          <-- index on entities for global_num[]
 *   global_num     <--
 *   comm           <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_fvmc_io_num_global_order_index(fvmc_io_num_t       *this_io_num,
                               fvmc_lnum_t          index[],
                               fvmc_gnum_t          global_num[],
                               MPI_Comm            comm)
{
  int  rank, local_rank, size;
  size_t  i, shift, slice_size;

  fvmc_gnum_t  n_ent_recv = 0, n_ent_send = 0;
  fvmc_gnum_t  current_global_num = 0, global_num_shift = 0;
  int  *send_count = NULL, *recv_count = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  fvmc_lnum_t  *recv_order = NULL, *recv_sub_index = NULL;
  fvmc_lnum_t  *recv_sub_count = NULL, *send_sub_count = NULL;
  fvmc_gnum_t  *block_global_num = NULL, *recv_global_num = NULL;

  /* Initialization */

  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &size);

  /* Get maximum global number value for first value of each series
     (does not need to be exact, simply used to define blocks) */

  {
    fvmc_gnum_t  local_max = 0, global_max = 0;
    size_t      n_ent = this_io_num->global_num_size;

    if (n_ent > 0)
      local_max = global_num[index[n_ent-1]];
    MPI_Allreduce(&local_max, &global_max, 1, FVMC_MPI_GNUM, MPI_MAX, comm);
    this_io_num->global_count = global_max;
  }

  /* slice_size = ceil(this_io_num->global_count/size) */

  slice_size = this_io_num->global_count / size;
  if (this_io_num->global_count % size > 0)
    slice_size += 1;

  /* Build for each slice, a new ordered indexed list from the received
     elements */

  assert(sizeof(fvmc_gnum_t) >= sizeof(fvmc_lnum_t));

  BFTC_MALLOC(send_count, size, int);
  BFTC_MALLOC(recv_count, size, int);
  BFTC_MALLOC(send_shift, size + 1, int);
  BFTC_MALLOC(recv_shift, size + 1, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / slice_size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  /* Get recv_index */

  n_ent_recv = recv_shift[size];
  n_ent_send = send_shift[size];

  BFTC_MALLOC(recv_sub_count, n_ent_recv, fvmc_lnum_t);
  BFTC_MALLOC(send_sub_count, n_ent_send, fvmc_lnum_t);

  assert(n_ent_send == (fvmc_gnum_t)this_io_num->global_num_size);

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)this_io_num->global_num_size; i++) {
    rank = (global_num[index[i]] - 1) / slice_size;
    shift = send_shift[rank] + send_count[rank];
    send_sub_count[shift] = index[i+1] - index[i];
    send_count[rank] +=  1;
  }

  MPI_Alltoallv(send_sub_count, send_count, send_shift, FVMC_MPI_LNUM,
                recv_sub_count, recv_count, recv_shift, FVMC_MPI_LNUM, comm);

  BFTC_MALLOC(recv_sub_index, n_ent_recv + 1, fvmc_lnum_t);

  recv_sub_index[0] = 0;
  for (i = 0; i < n_ent_recv; i++)
      recv_sub_index[i+1] = recv_sub_index[i] + recv_sub_count[i];

  BFTC_FREE(send_sub_count);
  BFTC_FREE(recv_sub_count);

  /* Get recv_global_num */

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / slice_size;
    send_count[rank] += index[i+1] - index[i];
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  BFTC_MALLOC(recv_global_num, recv_sub_index[n_ent_recv], fvmc_gnum_t);

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the indexed list to the correct processors */

  MPI_Alltoallv(global_num, send_count, send_shift, FVMC_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVMC_MPI_GNUM, comm);

  if (n_ent_recv > 0) { /* Order received elements of the indexed list */

    size_t prev_id, cur_id;

    BFTC_MALLOC(recv_order, n_ent_recv, fvmc_lnum_t);

    fvmc_order_local_allocated_i(NULL,
                                recv_global_num,
                                recv_sub_index,
                                recv_order,
                                n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of elements */

    BFTC_MALLOC(block_global_num, n_ent_recv, fvmc_gnum_t);

    current_global_num = 1;
    prev_id = recv_order[0];
    block_global_num[recv_order[0]] = current_global_num;

    for (i = 1; i < n_ent_recv; i++) {

      cur_id = recv_order[i];

      if (_indexed_is_greater(cur_id, prev_id, recv_sub_index, recv_global_num))
        current_global_num += 1;

      block_global_num[recv_order[i]] = current_global_num;
      prev_id = cur_id;

    }

  } /* End if n_ent_recv > 0 */

  /* Partial clean-up */

  BFTC_FREE(recv_order);
  BFTC_FREE(recv_sub_index);
  BFTC_FREE(recv_global_num);

  /* At this stage, block_global_num[] is valid for this process, and
     current_global_num indicates the total number of entities handled
     by this process; we must now shift global numberings on different
     processes by the cumulative total number of entities handled by
     each process */

  MPI_Scan(&current_global_num, &global_num_shift, 1, FVMC_MPI_GNUM,
           MPI_SUM, comm);
  global_num_shift -= current_global_num;

  for (i = 0; i < n_ent_recv; i++)
    block_global_num[i] += global_num_shift;

  /* Return global order to all processors */

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++) {
    rank = (global_num[index[i]] - 1) / slice_size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < size; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  MPI_Alltoallv(block_global_num, recv_count, recv_shift, FVMC_MPI_GNUM,
                this_io_num->_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                comm);

  /* Free memory */

  BFTC_FREE(block_global_num);
  BFTC_FREE(send_count);
  BFTC_FREE(recv_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_shift);

  /* Get final maximum global number value */

  this_io_num->global_count = _fvmc_io_num_global_max(this_io_num, comm);
}

/*----------------------------------------------------------------------------
 * Return the global number of sub-entities associated with an initial
 * entity whose global numbering is known, given the number of
 * sub-entities per initial entity.
 *
 * parameters:
 *   this_io_num    <-- pointer to base io numbering
 *   n_sub_entities <-- number of sub-entities per initial entity
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   global number of sub-entities
 *----------------------------------------------------------------------------*/

static fvmc_gnum_t
_fvmc_io_num_global_sub_size(const fvmc_io_num_t  *this_io_num,
                            const fvmc_lnum_t     n_sub_entities[],
                            MPI_Comm             comm)
{

  fvmc_gnum_t  global_count, n_ent_recv, num_prev, num_cur;
  size_t      i, slice_size;
  int         rank;

  fvmc_gnum_t  *recv_global_num = NULL;
  fvmc_gnum_t  *send_global_num = NULL;
  fvmc_lnum_t  *recv_n_sub = NULL, *recv_order = NULL;
  int         *send_count = NULL, *recv_count = NULL;
  int         *send_shift = NULL, *recv_shift = NULL;
  int         have_sub_loc = 0, have_sub_glob = 0;

  int         size;

  fvmc_gnum_t  current_global_num = 0;
  fvmc_gnum_t  retval = 0;

  /* Initialization */

  MPI_Comm_size(comm, &size);

  num_prev = 0;    /* true initialization later for slice 0, */

  /* Get temporary maximum global number value */

  global_count = _fvmc_io_num_global_max(this_io_num, comm);

  /* slice_size = ceil(this_io_num->global_count/size) */

  slice_size = global_count / size;
  if (global_count % size > 0)
    slice_size += 1;

  assert(sizeof(fvmc_gnum_t) >= sizeof(fvmc_lnum_t));

  BFTC_MALLOC(send_count, size, int);
  BFTC_MALLOC(recv_count, size, int);

  BFTC_MALLOC(send_shift, size, int);
  BFTC_MALLOC(recv_shift, size, int);

  /* Count number of values to send to each process */

  for (rank = 0; rank < size; rank++)
    send_count[rank] = 0;

  for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
    send_count[(this_io_num->global_num[i] - 1) / slice_size] += 1;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 1; rank < size; rank++) {
    send_shift[rank] = send_shift[rank - 1] + send_count[rank -1];
    recv_shift[rank] = recv_shift[rank - 1] + recv_count[rank -1];
  }

  /* As data is sorted by increasing base global numbering, we do not
     need to build an extra array, but only to send the correct parts
     of the n_sub_entities[] array to the correct processors */

  n_ent_recv = recv_shift[size - 1] + recv_count[size - 1];

  BFTC_MALLOC(recv_global_num, n_ent_recv, fvmc_gnum_t);
  BFTC_MALLOC(recv_order, n_ent_recv, fvmc_lnum_t);

  if (this_io_num->_global_num != NULL)
    send_global_num = this_io_num->_global_num;
  else {
    BFTC_MALLOC(send_global_num,
               this_io_num->global_num_size,
               fvmc_gnum_t);
    memcpy(send_global_num,
           this_io_num->global_num,
           this_io_num->global_num_size * sizeof(fvmc_gnum_t));
  }

  MPI_Alltoallv(send_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                recv_global_num, recv_count, recv_shift, FVMC_MPI_GNUM, comm);

  if (send_global_num != this_io_num->_global_num)
    BFTC_FREE(send_global_num);

  /* Do we have sub-entities ? */

  if (n_sub_entities != NULL)
    have_sub_loc = 1;

  MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX, comm);

  if (have_sub_glob > 0) {

    fvmc_lnum_t  *send_n_sub;

    BFTC_MALLOC(send_n_sub, this_io_num->global_num_size, fvmc_lnum_t);
    BFTC_MALLOC(recv_n_sub, n_ent_recv, fvmc_lnum_t);

    if (n_sub_entities != NULL) {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = n_sub_entities[i];
    }
    else {
      for (i = 0; i < (size_t)(this_io_num->global_num_size); i++)
        send_n_sub[i] = 1;
    }

    MPI_Alltoallv(send_n_sub, send_count, send_shift, FVMC_MPI_LNUM,
                  recv_n_sub, recv_count, recv_shift, FVMC_MPI_LNUM, comm);

    BFTC_FREE(send_n_sub);
  }

  if (n_ent_recv > 0) {

    fvmc_order_local_allocated(NULL,
                              recv_global_num,
                              recv_order,
                              n_ent_recv);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    current_global_num = recv_n_sub[recv_order[0]];
    num_prev = recv_global_num[recv_order[0]];
    recv_global_num[recv_order[0]] = current_global_num;

    for (i = 1; i < n_ent_recv; i++) {
      num_cur = recv_global_num[recv_order[i]];
      if (num_cur > num_prev)
        current_global_num += recv_n_sub[recv_order[i]];
      num_prev = num_cur;
    }

  }

  /* Partial clean-up */

  BFTC_FREE(recv_n_sub);
  BFTC_FREE(recv_order);
  BFTC_FREE(recv_global_num);

  BFTC_FREE(send_count);
  BFTC_FREE(recv_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_shift);

  /* At this stage, current_global_num indicates the total number of
     entities handled by this process; we must now shift global
     numberings on different processes by the cumulative total
     number of entities handled by each process */

  MPI_Allreduce(&current_global_num, &retval, 1, FVMC_MPI_GNUM, MPI_SUM, comm);

  return retval;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   parent_entity_number <-- pointer to list of selected entitie's parent's
 *                            numbers, or NULL if all first n_ent entities
 *                            are used
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   n_entities           <-- number of entities considered
 *   share_parent_global  <-- if non zero, try to share parent_global_number
 *                            instead of using a local copy
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create(const fvmc_lnum_t  parent_entity_number[],
                  const fvmc_gnum_t  parent_global_number[],
                  size_t            n_entities,
                  int               share_parent_global)
{
  fvmc_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (fvmc_parall_get_size() < 2)
    return NULL;

  assert(fvmc_order_local_test(parent_entity_number,
                              parent_global_number,
                              n_entities) == true);

#if defined(FVMC_HAVE_MPI)

  /* Create structure */

  BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFTC_MALLOC(this_io_num->_global_num, n_entities, fvmc_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  if (n_entities > 0) {

    size_t  i;

    /* Assign initial global numbers */

    if (parent_entity_number != NULL) {
      for (i = 0 ; i < n_entities ; i++)
        this_io_num->_global_num[i]
          = parent_global_number[parent_entity_number[i]-1];
    }
    else {
      for (i = 0 ; i < n_entities ; i++)
        this_io_num->_global_num[i] = parent_global_number[i];
    }

  }

  /* Order globally */

  this_io_num->global_count = n_entities;

  _fvmc_io_num_copy_on_write(this_io_num);
  _fvmc_io_num_global_order(this_io_num,
                           NULL,
                           fvmc_parall_get_mpi_comm());

  if (share_parent_global != 0)
    _fvmc_io_num_try_to_set_shared(this_io_num,
                                  parent_global_number);

#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure,
 * sharing a given global numbering array.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   global_number <-- pointer to list of global (i.e. domain splitting
 *                     independent) entity numbers
 *   global_count  <-- global number of entities
 *   n_entities    <-- number of local entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_shared(const fvmc_gnum_t  global_number[],
                         fvmc_gnum_t        global_count,
                         size_t            n_entities)
{
  fvmc_io_num_t  *this_io_num;

  /* Create structure */

  BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

  this_io_num->global_count = global_count;
  this_io_num->global_num_size = n_entities;

  this_io_num->global_num = global_number;
  this_io_num->_global_num = NULL;

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an an initial
 * I/O numbering and a number of new entities per base entity.
 *
 * This is useful for example to create an I/O numbering for
 * triangles based on split polygons, whose I/O numbering is defined.
 *
 * parameters:
 *   base_io_num    <-- pointer to base I/O numbering structure
 *   n_sub_entities <-- number of new entities per base entity
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_from_sub(const fvmc_io_num_t  *base_io_num,
                           const fvmc_lnum_t     n_sub_entities[])
{
  fvmc_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (base_io_num == NULL)
    return NULL;

  assert(fvmc_parall_get_size() > 1); /* Otherwise, base_io_num should be NULL */

#if defined(FVMC_HAVE_MPI)
  {
    fvmc_lnum_t  i, n_ent;

    /* Create structure */

    BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

    n_ent = base_io_num->global_num_size;

    BFTC_MALLOC(this_io_num->_global_num, n_ent, fvmc_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    this_io_num->global_num_size = n_ent;

    /* Assign initial global numbers */

    for (i = 0 ; i < n_ent ; i++)
      this_io_num->_global_num[i] = base_io_num->global_num[i];

    /* Order globally */

    this_io_num->global_count = n_ent;

    _fvmc_io_num_copy_on_write(this_io_num);
    _fvmc_io_num_global_order(this_io_num,
                             n_sub_entities,
                             fvmc_parall_get_mpi_comm());
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a strided adjacency.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   parent_entity_number <-- pointer to list of selected entitie's parent's
 *                            numbers, or NULL if all first n_ent entities
 *                            are used
 *   adjacency            <-- entity adjacency (1 to n global numbering)
 *   n_entities           <-- number of entities considered
 *   stride               <-- values per entity
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_from_adj_s(const fvmc_lnum_t  parent_entity_number[],
                             const fvmc_gnum_t  adjacency[],
                             size_t            n_entities,
                             size_t            stride)
{
  fvmc_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (fvmc_parall_get_size() < 2)
    return NULL;

  assert(fvmc_order_local_test_s(parent_entity_number,
                                adjacency,
                                stride,
                                n_entities) == true);

#if defined(FVMC_HAVE_MPI)
  {
    fvmc_gnum_t *_adjacency = NULL;

    /* Create structure */

    BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

    this_io_num->global_num_size = n_entities;

    BFTC_MALLOC(this_io_num->_global_num, n_entities, fvmc_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    if (n_entities > 0) {

      size_t  i, j;

      /* Assign initial global numbers */

      BFTC_MALLOC(_adjacency, n_entities*stride, fvmc_gnum_t);

      if (parent_entity_number != NULL) {
        for (i = 0 ; i < n_entities ; i++) {
          for (j = 0; j < stride; j++)
            _adjacency[i*stride + j]
              = adjacency[(parent_entity_number[i]-1)*stride + j];
        }
      }
      else
        memcpy(_adjacency, adjacency, n_entities*stride*sizeof(fvmc_gnum_t));

    }

    /* Order globally */

    this_io_num->global_count = n_entities;

    _fvmc_io_num_global_order_s(this_io_num,
                               stride,
                               _adjacency,
                               fvmc_parall_get_mpi_comm());

    BFTC_FREE(_adjacency);
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an indexed adjacency.
 *
 * The corresponding entities do not need to be locally ordered.
 *
 * parameters:
 *  parent_entity_number <-- pointer to list of selected entitie's parent's
 *                           numbers, or NULL if all first n_ent entities
 *                           are used
 *  index                <-- index on entities for adjacency
 *  adjacency            <-- entity adjacency (1 to n global numbering)
 *  n_entities           <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_from_adj_i(const fvmc_lnum_t  parent_entity_number[],
                             const fvmc_lnum_t  index[],
                             const fvmc_gnum_t  adjacency[],
                             fvmc_lnum_t        n_entities)
{
  fvmc_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (fvmc_parall_get_size() < 2)
    return NULL;

#if defined(FVMC_HAVE_MPI)
  {
    fvmc_lnum_t  *_index = NULL;
    fvmc_gnum_t *_adjacency = NULL;

#if defined(DEBUG) && !defined(NDEBUG)
    const char no_adjacent_elt_msg[]
      = " Error during the creation of a fvmc_io_num_t structure\n"
        " for an indexed adjacency.\n\n"
        " At least one entity has no adjacent element, and\n"
        " this case is not currently handled.\n\n"
        " Check if this definition is correct.";
#endif

    /* Create structure */

    BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

    this_io_num->global_num_size = n_entities;

    BFTC_MALLOC(this_io_num->_global_num, n_entities, fvmc_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    if (n_entities > 0) {

      fvmc_lnum_t  i, j, k, ent_id, _shift;

      /* Assign initial global numbers */

      BFTC_MALLOC(_index, n_entities + 1, fvmc_lnum_t);
      _index[0] = 0;

      if (parent_entity_number != NULL) {

#if defined(DEBUG) && !defined(NDEBUG)
        for (i = 0 ; i < n_entities ; i++) {
          ent_id = parent_entity_number[i]-1;
          if ((index[ent_id+1] - index[ent_id]) == 0)
            bftc_error(__FILE__, __LINE__, 0, no_adjacent_elt_msg);
        }
#endif

        /* Count reduced size */

        for (i = 0 ; i < n_entities ; i++) {
          ent_id = parent_entity_number[i]-1;
          _index[i+1] = index[ent_id+1] - index[ent_id];
        }

        for (i = 0 ; i < n_entities ; i++)
          _index[i+1] += _index[i];

        BFTC_MALLOC(_adjacency, _index[n_entities], fvmc_gnum_t);

        /* Define reduced index and adjacency */

        for (i = 0 ; i < n_entities ; i++) {

          ent_id = parent_entity_number[i]-1;
          _shift = _index[i];

          for (j = index[ent_id], k = 0; j < index[ent_id+1]; j++, k++)
            _adjacency[_shift + k] = adjacency[j];

        }

      }
      else {

#if defined(DEBUG) && !defined(NDEBUG)
        for (i = 0 ; i < n_entities ; i++) {
          if ((index[i+1] - index[i]) == 0)
            bftc_error(__FILE__, __LINE__, 0, no_adjacent_elt_msg);
        }
#endif

        BFTC_MALLOC(_adjacency, index[n_entities], fvmc_gnum_t);

        memcpy(_index, index, (n_entities+1)*sizeof(fvmc_lnum_t));
        memcpy(_adjacency, adjacency, index[n_entities]*sizeof(fvmc_gnum_t));

      }

    }

    /* Order globally */

    this_io_num->global_count = n_entities;

    _fvmc_io_num_global_order_index(this_io_num,
                                   _index,
                                   _adjacency,
                                   fvmc_parall_get_mpi_comm());

    if (_adjacency != NULL)
      BFTC_FREE(_adjacency);
    if (_index != NULL)
      BFTC_FREE(_index);
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   dim        <-- spatial dimension
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_from_coords(const fvmc_coord_t  coords[],
                              int                dim,
                              size_t             n_entities)
{
  size_t  i;
  fvmc_coord_t extents[6];
  fvmc_lnum_t  *order = NULL;
  fvmc_morton_code_t *m_code = NULL;

#if defined(FVMC_HAVE_MPI)
  MPI_Comm comm = fvmc_parall_get_mpi_comm();
#endif

  const int level = sizeof(fvmc_morton_int_t)*8 - 1;
  const int n_ranks = fvmc_parall_get_size();

  fvmc_io_num_t  *this_io_num = NULL;

  /* Create structure */

  BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

  this_io_num->global_num_size = n_entities;

  BFTC_MALLOC(this_io_num->_global_num, n_entities, fvmc_gnum_t);
  this_io_num->global_num = this_io_num->_global_num;

  /* Build Morton encoding and order it */

#if defined(FVMC_HAVE_MPI)
  fvmc_morton_get_coord_extents(dim, n_entities, coords, extents, comm);
#else
  fvmc_morton_get_coord_extents(dim, n_entities, coords, extents);
#endif

  BFTC_MALLOC(m_code, n_entities, fvmc_morton_code_t);
  BFTC_MALLOC(order, n_entities, fvmc_lnum_t);

  fvmc_morton_encode_coords(dim, level, extents, n_entities, coords, m_code);

  fvmc_morton_local_order(n_entities, m_code, order);

#if defined(FVMC_HAVE_MPI)

  if (n_ranks > 1) {

    int rank_id;
    fvmc_lnum_t j, shift;

    size_t n_block_ents = 0;
    fvmc_gnum_t current_global_num = 0, global_num_shift = 0;

    int *c_rank = NULL;
    int *send_count = NULL, *send_shift = NULL;
    int *recv_count = NULL, *recv_shift = NULL;
    fvmc_coord_t *send_coords = NULL, *recv_coords = NULL;
    fvmc_lnum_t *weight = NULL;
    fvmc_gnum_t *block_global_num = NULL, *part_global_num = NULL;
    fvmc_morton_code_t *morton_index = NULL;

    BFTC_MALLOC(weight, n_entities, fvmc_lnum_t);
    BFTC_MALLOC(morton_index, n_ranks + 1, fvmc_morton_code_t);

    for (i = 0; i < n_entities; i++)
      weight[i] = 1;

    fvmc_morton_build_rank_index(dim,
                                 level,
                                 n_entities,
                                 m_code,
                                 weight,
                                 order,
                                 morton_index,
                                 comm);

    BFTC_FREE(order);
    BFTC_FREE(weight);
    BFTC_MALLOC(c_rank, n_entities, int);

    for (i = 0; i < n_entities; i++)
      c_rank[i] = fvmc_morton_quantile_search(n_ranks,
                                             m_code[i],
                                             morton_index);

    BFTC_FREE(morton_index);
    BFTC_FREE(m_code);

    /* Build send_buf, send_count and send_shift
       to build a rank to coords indexed list */

    BFTC_MALLOC(send_count, n_ranks, int);
    BFTC_MALLOC(recv_count, n_ranks, int);
    BFTC_MALLOC(send_shift, n_ranks + 1, int);
    BFTC_MALLOC(recv_shift, n_ranks + 1, int);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (i = 0; i < n_entities; i++)
      send_count[c_rank[i]] += dim;

    /* Exchange number of coords to send to each process */

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
      recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
    }

    /* Build send and receive buffers */

    BFTC_MALLOC(send_coords, send_shift[n_ranks], fvmc_coord_t);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (i = 0; i < n_entities; i++) {
      rank_id = c_rank[i];
      shift = send_shift[rank_id] + send_count[rank_id];
      for (j = 0; j < dim; j++)
        send_coords[shift + j] = coords[i*dim + j];
      send_count[rank_id] += dim;
    }

    BFTC_MALLOC(recv_coords, recv_shift[n_ranks], fvmc_coord_t);

    /* Exchange coords between processes */

    MPI_Alltoallv(send_coords, send_count, send_shift, FVMC_MPI_COORD,
                  recv_coords, recv_count, recv_shift, FVMC_MPI_COORD,
                  comm);

    BFTC_FREE(send_coords);

    /* Now re-build Morton codes on block distribution */

    n_block_ents = recv_shift[n_ranks] / dim;

    BFTC_MALLOC(m_code, n_block_ents, fvmc_morton_code_t);
    BFTC_MALLOC(order, n_block_ents, fvmc_lnum_t);

    fvmc_morton_encode_coords(dim,
                             level,
                             extents,
                             n_block_ents,
                             recv_coords,
                             m_code);

    fvmc_morton_local_order(n_block_ents, m_code, order);

    /* Check ordering; if two entities have the same Morton codes,
       use lexicographical coordinates ordering to ensure the
       final order is deterministic. */

    _check_morton_ordering(dim, n_block_ents, recv_coords, m_code, order);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    BFTC_FREE(m_code);
    BFTC_FREE(recv_coords);
    BFTC_MALLOC(block_global_num, n_block_ents, fvmc_gnum_t);

    for (i = 0; i < n_block_ents; i++)
      block_global_num[order[i]] = i+1;

    BFTC_FREE(order);

    current_global_num = n_block_ents;

    /* At this stage, block_global_num[] is valid for this process, and
       current_global_num indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    MPI_Scan(&current_global_num, &global_num_shift, 1, FVMC_MPI_GNUM,
             MPI_SUM, comm);
    global_num_shift -= current_global_num;

    for (i = 0; i < n_block_ents; i++)
      block_global_num[i] += global_num_shift;

    /* Return global order to all processors */

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] /= dim;
      recv_count[rank_id] /= dim;
      send_shift[rank_id] /= dim;
      recv_shift[rank_id] /= dim;
    }

    send_shift[n_ranks] /= dim;

    BFTC_MALLOC(part_global_num, send_shift[n_ranks], fvmc_gnum_t);

    MPI_Alltoallv(block_global_num, recv_count, recv_shift, FVMC_MPI_GNUM,
                  part_global_num, send_count, send_shift, FVMC_MPI_GNUM,
                  comm);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (i = 0; i < n_entities; i++) {
      rank_id = c_rank[i];
      shift = send_shift[rank_id] + send_count[rank_id];
      this_io_num->_global_num[i] = part_global_num[shift];
      send_count[rank_id] += 1;
    }

    /* Free memory */

    BFTC_FREE(c_rank);

    BFTC_FREE(block_global_num);
    BFTC_FREE(part_global_num);

    BFTC_FREE(send_count);
    BFTC_FREE(recv_count);
    BFTC_FREE(send_shift);
    BFTC_FREE(recv_shift);

    /* Get final maximum global number value */

    this_io_num->global_count = _fvmc_io_num_global_max(this_io_num, comm);

  }

#endif /* FVMC_HAVE_MPI */

  if (n_ranks == 1) {

    _check_morton_ordering(dim, n_entities, coords, m_code, order);

    BFTC_FREE(m_code);

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[order[i]] = i+1;

    BFTC_FREE(order);

    this_io_num->global_count = n_entities;

  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a simple accumulation
 * (i.e. scan) of counts on successive ranks.
 *
 * parameters:
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_create_from_scan(size_t  n_entities)
{
  fvmc_io_num_t  *this_io_num = NULL;

  /* Initial checks */

  if (fvmc_parall_get_size() < 2)
    return NULL;

#if defined(FVMC_HAVE_MPI)
  {
    size_t  i;
    fvmc_gnum_t gnum_base = n_entities;
    fvmc_gnum_t gnum_sum = n_entities;
    fvmc_gnum_t gnum_shift = 0;

    MPI_Comm comm = fvmc_parall_get_mpi_comm();

    /* Create structure */

    BFTC_MALLOC(this_io_num, 1, fvmc_io_num_t);

    BFTC_MALLOC(this_io_num->_global_num, n_entities, fvmc_gnum_t);
    this_io_num->global_num = this_io_num->_global_num;

    this_io_num->global_num_size = n_entities;

    MPI_Scan(&gnum_base, &gnum_shift, 1, FVMC_MPI_GNUM, MPI_SUM, comm);

    gnum_base = gnum_shift - gnum_base + 1;

    for (i = 0; i < n_entities; i++)
      this_io_num->_global_num[i] = gnum_base + i;

    gnum_base = n_entities;

    MPI_Allreduce(&gnum_base, &gnum_sum, 1, FVMC_MPI_GNUM, MPI_SUM, comm);

    this_io_num->global_count = gnum_sum;
  }
#endif

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Destruction of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_io_num_t *
fvmc_io_num_destroy(fvmc_io_num_t  * this_io_num)
{
  if (this_io_num != NULL) {
    BFTC_FREE(this_io_num->_global_num);
    BFTC_FREE(this_io_num);
  }

  return this_io_num;
}

/*----------------------------------------------------------------------------
 * Return local number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  local number of associated entities
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_io_num_get_local_count(const fvmc_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_num_size;
}

/*----------------------------------------------------------------------------
 * Return global number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  global number of associated entities
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_io_num_get_global_count(const fvmc_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_count;
}

/*----------------------------------------------------------------------------
 * Return global numbering associated with an I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  pointer to array of global numbers associated with local entities
 *  (1 to n numbering)
 *----------------------------------------------------------------------------*/

const fvmc_gnum_t *
fvmc_io_num_get_global_num(const fvmc_io_num_t  *const this_io_num)
{
  assert(this_io_num != NULL);

  return this_io_num->global_num;
}

/*----------------------------------------------------------------------------
 * Return the global number of sub-entities associated with an initial
 * entity whose global numbering is known, given the number of
 * sub-entities per initial entity.
 *
 * parameters:
 *   this_io_num    <-> pointer to base io numbering
 *   n_sub_entities <-- number of sub-entities per initial entity
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   global number of sub-entities
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_io_num_global_sub_size(const fvmc_io_num_t  *this_io_num,
                           const fvmc_lnum_t     n_sub_entities[])
{
  fvmc_gnum_t  retval = 0;

  /* Initial checks */

  if (this_io_num == NULL)
    return retval;

  assert(fvmc_parall_get_size() > 1); /* Otherwise, base_io_num should be NULL */

#if defined(FVMC_HAVE_MPI)
 {
   int  have_sub_loc = 0, have_sub_glob = 0;

   /* Caution: we may have sub-entities on some ranks and not on others */

   if (n_sub_entities != NULL)
     have_sub_loc = 1;

   MPI_Allreduce(&have_sub_loc, &have_sub_glob, 1, MPI_INT, MPI_MAX,
                 fvmc_parall_get_mpi_comm());

   if (have_sub_glob > 0)
     retval = _fvmc_io_num_global_sub_size(this_io_num,
                                          n_sub_entities,
                                          fvmc_parall_get_mpi_comm());
 }
#endif

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump printout of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvmc_io_num_dump(const fvmc_io_num_t  *const this_io_num)
{
  fvmc_lnum_t i;

  if (this_io_num == NULL) {
    bftc_printf("  global numbering: nil\n");
    return;
  }

  bftc_printf("  global numbering size:            %u\n",
             (unsigned)this_io_num->global_num_size);

  bftc_printf("\n"
             "  pointer to shareable array:\n"
             "    global_num:                     %p\n",
             this_io_num->global_num);

  bftc_printf("\n"
             "  pointer to local array:\n"
             "    _global_num:                    %p\n",
             this_io_num->global_num);

  if (this_io_num->global_num_size > 0) {

    bftc_printf("\n  global number:\n\n");
    for (i = 0 ; i < this_io_num->global_num_size ; i++)
      bftc_printf("  %10u : %10lu\n",
                 (unsigned)i + 1,
                 (unsigned long)this_io_num->global_num[i]);
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
