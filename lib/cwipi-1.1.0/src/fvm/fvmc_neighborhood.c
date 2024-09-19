/*============================================================================
 * Determine and update geometrical neighborhood information.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2008-2009  EDF

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
#include <limits.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>
#include <bftc_timer.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_order.h"
#include "fvmc_parall.h"
#include "fvmc_part_to_block.h"

#include "fvmc_box.h"
#include "fvmc_box_tree.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_neighborhood.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro and Type definitions
 *============================================================================*/

/* Neighborhood statistics (using box tree) */
/*------------------------------------------*/

typedef struct {

  int         dim;                     /* Layout dimension */

  /* The following fields have 3 global values:
     mean on ranks, minimum on ranks, and maximum on ranks */

  int         depth[3];                /* Tree depth */
  fvmc_lnum_t  n_leaves[3];             /* Number of leaves */
  fvmc_lnum_t  n_boxes[3];              /* Number of associated boxes */
  fvmc_lnum_t  n_threshold_leaves[3];   /* Number of leaves over threshold */
  fvmc_lnum_t  n_leaf_boxes[3];         /* Number of boxes per leaf */
  size_t      mem_used[3];             /* Memory used */
  size_t      mem_required[3];         /* Memory temporarily required */

} _box_tree_stats_t;

/* Main neighborhood structure */
/*-----------------------------*/

struct _fvmc_neighborhood_t {

  fvmc_lnum_t        n_elts;          /* Number of elements */

  fvmc_gnum_t       *elt_num;         /* Global numbers associated with
                                        elements in local block
                                        (size: n_elts) */
  fvmc_lnum_t       *neighbor_index;  /* Start index of neighbors
                                        (size: n_elts + 1) */
  fvmc_gnum_t       *neighbor_num;    /* Global element neighbor numbers
                                        (size: neighbor_index[n_elts]) */

#if defined(FVMC_HAVE_MPI)
  MPI_Comm  comm;                    /* Associated MPI communicator */
#endif

  /* Algorithm-related options */

  int  max_tree_depth;               /* Maximum search tree depth */
  int  leaf_threshold;               /* Maximum number of boxes which can
                                        be related to a leaf of the tree if
                                        level < max_tree_depth */
  int  max_box_ratio;                /* Stop adding levels to tree when
                                        (  n_linked_boxes
                                         > max_box_ratio*n_initial_boxes) */

  _box_tree_stats_t  bt_stats;       /* Statistics associated with the
                                        box-trees used for search */

  /* Timings */

  double  cpu_time[2];   /* CPU time for tree construction and query */
  double  wtime[2];      /* Wall clock time for tree construction and query */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize box_tree statistics
 *
 * parameters:
 *   bts --> pointer to box tree statistics structure
 *---------------------------------------------------------------------------*/

static void
_init_bt_statistics(_box_tree_stats_t  *bts)
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}
/*----------------------------------------------------------------------------
 * Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * parameters:
 *   bts   <-> pointer to box tree statistics structure
 *   boxes <-- pointer to box tree structure
 *---------------------------------------------------------------------------*/

static void
_update_bt_statistics(_box_tree_stats_t     *bts,
                      const fvmc_box_tree_t  *bt)
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = fvmc_box_tree_get_stats(bt,
                               bts->depth,
                               bts->n_leaves,
                               bts->n_boxes,
                               bts->n_threshold_leaves,
                               bts->n_leaf_boxes,
                               bts->mem_used,
                               mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = FVMC_MAX(bts->mem_required[i], mem_required[i]);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Distribute bounding boxes over the ranks according to a Morton encoding
 * index. Try to get a well-balanced distribution and spatially coherent.
 *
 * parameters:
 *   n     <-> pointer to neighborhood management structure
 *   boxes <-> box set to redistribute
 *---------------------------------------------------------------------------*/

static void
_redistribute_boxes(fvmc_neighborhood_t  *n,
                    fvmc_box_set_t       *boxes)
{
  fvmc_box_tree_t  *coarse_tree = NULL;
  fvmc_box_distrib_t  *distrib = NULL;

  const int  max_box_ratio = FVMC_MIN(6, n->max_box_ratio);

  /* Sanity checks */

  assert(boxes != NULL);

  coarse_tree = fvmc_box_tree_create(n->max_tree_depth,
                                    n->leaf_threshold,
                                    max_box_ratio);

  /* Build a tree and associate boxes */

  fvmc_box_tree_set_boxes(coarse_tree,
                         boxes,
                         FVMC_BOX_TREE_SYNC_LEVEL);

  _update_bt_statistics(&(n->bt_stats), coarse_tree);

  /*
    Compute an index based on Morton encoding to ensure a good distribution
    of bounding boxes among the ranks.
  */

  distrib = fvmc_box_tree_get_distrib(coarse_tree, boxes);

  fvmc_box_tree_destroy(&coarse_tree);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  fvmc_box_distrib_dump_statistics(distrib, n->comm);
#endif

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  fvmc_box_set_redistribute(distrib, boxes);

  /* Delete intermediate structures */

  fvmc_box_distrib_destroy(&distrib);
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * using a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

static inline void
_gnum_shellsort(fvmc_lnum_t  l,
                fvmc_lnum_t  r,
                fvmc_gnum_t  a[])
{
  fvmc_lnum_t i, j, h;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1);

  /* Sort array */
  for (; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      fvmc_gnum_t  v = a[i];

      j = i;
      while ((j >= l+h) && (v < a[j-h])) {
        a[j] = a[j-h];
        j -= h;
      }
      a[j] = v;

    } /* Loop on array elements */

  } /* End of loop on stride */
}

/*----------------------------------------------------------------------------
 * Remove numbers with no neighbors
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

static void
_remove_nums_without_neighbors(fvmc_neighborhood_t  *n)
{
  fvmc_lnum_t  i, start_id, end_id, saved_id, n_elts;

  fvmc_lnum_t  e_count = 0;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;

  /* Remove elements with no neighbors */

  saved_id = n->neighbor_index[0];

  for (i = 0; i < n_elts; i++) {

    start_id = saved_id;
    end_id = n->neighbor_index[i+1];

    if (end_id - start_id > 0) {

      n->elt_num[e_count] = n->elt_num[i];

      saved_id = end_id;
      n->neighbor_index[e_count+1] = end_id;

      e_count += 1;

    }

  }

  if (e_count < n_elts) {
    n->n_elts = e_count;
    BFTC_REALLOC(n->elt_num, e_count, fvmc_gnum_t);
    BFTC_REALLOC(n->neighbor_index, e_count + 1, fvmc_lnum_t);
  }
}

/*----------------------------------------------------------------------------
 * Remove multiple element neighbor occurences
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

static void
_clean_neighbor_nums(fvmc_neighborhood_t  *n)
{
  fvmc_lnum_t  i, j, start_id, end_id, saved_id, n_elts, n_neighbors;

  fvmc_lnum_t  n_count = 0;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;
  n_neighbors = n->neighbor_index[n_elts];

  /* Remove redundant elements */

  saved_id = n->neighbor_index[0];

  for (i = 0; i < n_elts; i++) {

    start_id = saved_id;
    end_id = n->neighbor_index[i+1];

    if (end_id - start_id > 0) {

      _gnum_shellsort(start_id, end_id, n->neighbor_num);

      n->neighbor_num[n_count++] = n->neighbor_num[start_id];

      for (j = start_id + 1; j < end_id; j++) {
        if (n->neighbor_num[j] != n->neighbor_num[j-1])
          n->neighbor_num[n_count++] = n->neighbor_num[j];
      }

    }

    saved_id = end_id;
    n->neighbor_index[i+1] = n_count;

  } /* End of loop on elements */

  if (n_count < n_neighbors)
    BFTC_REALLOC(n->neighbor_num, n_count, fvmc_gnum_t);
}

/*----------------------------------------------------------------------------
 * Order a neighborhood based on the global numbering of elements.
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

static void
_order_neighborhood(fvmc_neighborhood_t  *n)
{
  fvmc_lnum_t  i, j, k, order_id, shift, e_count;
  fvmc_lnum_t  n_elts, n_neighbors, n_elt_neighbors;
  fvmc_gnum_t  prev_num, cur_num;

  fvmc_lnum_t  *order = NULL, *old_index = NULL;
  fvmc_gnum_t  *old_e_num = NULL, *old_n_num = NULL;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;
  n_neighbors = n->neighbor_index[n_elts];

  BFTC_MALLOC(order, n_elts, fvmc_lnum_t);
  BFTC_MALLOC(old_e_num, n_elts, fvmc_gnum_t);
  BFTC_MALLOC(old_index, n_elts + 1, fvmc_lnum_t);
  BFTC_MALLOC(old_n_num, n_neighbors, fvmc_gnum_t);

  memcpy(old_e_num, n->elt_num, n_elts*sizeof(fvmc_gnum_t));
  memcpy(old_index, n->neighbor_index, (n_elts + 1)*sizeof(fvmc_gnum_t));
  memcpy(old_n_num, n->neighbor_num, n_neighbors*sizeof(fvmc_gnum_t));

  /* Order elt_num */

  fvmc_order_local_allocated(NULL, old_e_num, order, n_elts);

  /* Reshape according to the new ordering */

  /* Add first element */

  order_id = order[0];
  shift = 0;

  n->elt_num[0] = old_e_num[order_id];
  prev_num = n->elt_num[0];

  n->neighbor_index[0] = 0;
  n->neighbor_index[1] = old_index[order_id+1] - old_index[order_id];

  /* Loop on second-to last elements, merging data if an element has
     already appeared */

  for (i = 1, e_count = 1; i < n_elts; i++) {

    order_id = order[i];

    n_elt_neighbors = old_index[order_id+1] - old_index[order_id];

    shift = n->neighbor_index[i];

    cur_num = old_e_num[order_id];

    if (cur_num != prev_num) {
      n->elt_num[e_count] = cur_num;
      n->neighbor_index[e_count+1] = (  n->neighbor_index[e_count]
                                      + n_elt_neighbors);
      e_count += 1;
      prev_num = cur_num;
    }
    else
      n->neighbor_index[e_count] += n_elt_neighbors;

    for (j = old_index[order_id], k = 0; k < n_elt_neighbors; j++, k++)
      n->neighbor_num[shift + k] = old_n_num[j];

  } /* End of loop on elements */

  assert(n->neighbor_index[e_count] == n_neighbors);

  /* Free temporary memory */

  BFTC_FREE(order);
  BFTC_FREE(old_e_num);
  BFTC_FREE(old_index);
  BFTC_FREE(old_n_num);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Synchronize a neighborhood management structure and distribute the
 * resulting data over ranks by block.
 *
 * parameters:
 *   n        <-- pointer to neighborhood management structure
 *   n_g_elts <--  global number of elements
 *----------------------------------------------------------------------------*/

static void
_sync_by_block(fvmc_neighborhood_t  *n,
               fvmc_gnum_t           n_g_elts)
{
  fvmc_lnum_t  i, j;
  int  rank_id, n_ranks, n_recv_elts, n_sub_elts, shift;

  int  *send_count = NULL, *recv_count = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  fvmc_lnum_t  *counter = NULL;
  fvmc_gnum_t  *send_buf = NULL, *recv_buf = NULL;

  assert(n != NULL);

  if (n_g_elts == 0 || n->comm == MPI_COMM_NULL)
    return;

  /* Initialization */

  MPI_Comm_rank(n->comm, &rank_id);
  MPI_Comm_size(n->comm, &n_ranks);

  fvmc_part_to_block_info_t bi = fvmc_part_to_block_compute_sizes(rank_id,
                                                                n_ranks,
                                                                0,
                                                                0,
                                                                n_g_elts);

  /* Allocate buffers */

  BFTC_MALLOC(send_count, n_ranks, int);
  BFTC_MALLOC(recv_count, n_ranks, int);
  BFTC_MALLOC(send_shift, n_ranks + 1, int);
  BFTC_MALLOC(recv_shift, n_ranks + 1, int);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  /* Synchronize definition for each global element */

  for (i = 0; i < n->n_elts; i++) {
    long int g_ent_id = n->elt_num[i] -1;
    int send_rank = (g_ent_id/bi.block_size)*bi.rank_step;
    n_sub_elts = n->neighbor_index[i+1] - n->neighbor_index[i];
    send_count[send_rank] += 2 + n_sub_elts;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, n->comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (i = 0; i < n_ranks; i++) {
    send_shift[i + 1] = send_shift[i] + send_count[i];
    recv_shift[i + 1] = recv_shift[i] + recv_count[i];
  }

  /* Fill send_buf: global number and number of elements in index */

  BFTC_MALLOC(send_buf, send_shift[n_ranks], fvmc_gnum_t);
  BFTC_MALLOC(recv_buf, recv_shift[n_ranks], fvmc_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n->n_elts; i++) {

    long int g_ent_num = n->elt_num[i];
    int send_rank = ((g_ent_num - 1)/bi.block_size)*bi.rank_step;

    shift = send_shift[send_rank] + send_count[send_rank];
    n_sub_elts = n->neighbor_index[i+1] - n->neighbor_index[i];

    send_buf[shift++] = g_ent_num;
    send_buf[shift++] = n_sub_elts;

    for (j = 0; j < n_sub_elts; j++)
      send_buf[shift + j] = n->neighbor_num[n->neighbor_index[i] + j];

    send_count[send_rank] += 2 + n_sub_elts;
  }

  MPI_Alltoallv(send_buf, send_count, send_shift, FVMC_MPI_GNUM,
                recv_buf, recv_count, recv_shift, FVMC_MPI_GNUM,
                n->comm);

  n_recv_elts = recv_shift[n_ranks];

  /* Free what we may */

  BFTC_FREE(send_buf);
  BFTC_FREE(send_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_count);
  BFTC_FREE(recv_shift);

  /* Build arrays corresponding to block distribution of neighborhood */

  n->n_elts = bi.gnum_range[1] - bi.gnum_range[0];

  BFTC_FREE(n->elt_num);
  BFTC_FREE(n->neighbor_index);
  BFTC_FREE(n->neighbor_num);

  BFTC_MALLOC(n->elt_num, n->n_elts, fvmc_gnum_t);
  BFTC_MALLOC(n->neighbor_index, n->n_elts + 1, fvmc_lnum_t);

  for (i = 0; i < n->n_elts; i++) {
    n->elt_num[i] = bi.gnum_range[0] + i;
    n->neighbor_index[i] = 0;
  }
  n->neighbor_index[n->n_elts] = 0;

  /* Count element neighbors in block distribution */

  i = 0; /* start position in recv_buf */

  while (i < n_recv_elts) {

    fvmc_lnum_t elt_id = recv_buf[i++] - bi.gnum_range[0];

    assert(n->elt_num[elt_id] == recv_buf[i-1]);

    n_sub_elts = recv_buf[i++];
    n->neighbor_index[elt_id + 1] += n_sub_elts;
    i += n_sub_elts;
  }

  /* Transform element neighbor count to index */

  n->neighbor_index[0] = 0;
  for (i = 0; i < n->n_elts; i++)
    n->neighbor_index[i+1] += n->neighbor_index[i];

  BFTC_MALLOC(n->neighbor_num, n->neighbor_index[n->n_elts], fvmc_gnum_t);

  /* Fill element neighbors in block distribution */

  BFTC_MALLOC(counter, n->n_elts, fvmc_lnum_t);

  for (i = 0; i < n->n_elts; i++)
    counter[i] = 0;

  i = 0; /* start position in recv_buf */

  while (i < n_recv_elts) {

    fvmc_lnum_t elt_id = recv_buf[i++] - bi.gnum_range[0];

    n_sub_elts = recv_buf[i++];

    shift = n->neighbor_index[elt_id] + counter[elt_id];

    for (j = 0; j < n_sub_elts; j++)
      n->neighbor_num[j + shift] = recv_buf[i++];

    counter[elt_id] += n_sub_elts;

  } /* End of loop on ranks */

  BFTC_FREE(recv_buf);
  BFTC_FREE(counter);

  /* Remove redundant data */

  _clean_neighbor_nums(n);
}

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a neighborhood_t structure and initialize it.
 *
 * parameters:
 *   comm  <-- associated MPI communicator
 *
 * returns:
 *   pointer to an empty fvmc_box_tree_t structure.
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)
fvmc_neighborhood_t *
fvmc_neighborhood_create(MPI_Comm  comm)
#else
fvmc_neighborhood_t *
fvmc_neighborhood_create(void)
#endif
{
  double  w_start, w_end, cpu_start, cpu_end;

  fvmc_neighborhood_t *n = NULL;

  /* Timer start */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  /* Allocate and initialize */

  BFTC_MALLOC(n, 1, fvmc_neighborhood_t);

  n->n_elts = 0;
  n->elt_num = NULL;
  n->neighbor_index = NULL;
  n->neighbor_num = NULL;

#if defined(FVMC_HAVE_MPI)
  n->comm = comm;
#endif

  /* Algorithm options */

  n->max_tree_depth = 30; /* defaults */
  n->leaf_threshold = 30;
  n->max_box_ratio = 10;

  _init_bt_statistics(&(n->bt_stats));

  /* Timer end */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  n->cpu_time[0] = cpu_end - cpu_start;  /* build time */
  n->wtime[0] = w_end - w_start;
  n->cpu_time[1] = 0.0;                  /* query */
  n->wtime[1] = 0.0;

  /* Return structure */

  return n;
}

/*----------------------------------------------------------------------------
 * Destroy a neighborhood_t structure.
 *
 * parameters:
 *   n <-> pointer to pointer to fvmc_neighborhood_t structure to destroy.
 *----------------------------------------------------------------------------*/

void
fvmc_neighborhood_destroy(fvmc_neighborhood_t  **n)
{
  if (n != NULL) {
    fvmc_neighborhood_t *_n = *n;
    if (_n != NULL) {
      if (_n->elt_num != NULL)
        BFTC_FREE(_n->elt_num);
      if (_n->neighbor_index != NULL)
        BFTC_FREE(_n->neighbor_index);
      if (_n->neighbor_num != NULL)
        BFTC_FREE(_n->neighbor_num);
    }
    BFTC_FREE(*n);
  }
}

/*----------------------------------------------------------------------------
 * Set non-default algorithm parameters for neighborhood management structure.
 *
 * parameters:
 *   n              <-> pointer to neighborhood management structure
 *   max_tree_depth <-- maximum search tree depth
 *   leaf_threshold <-- maximum number of boxes which can be related to
 *                      a leaf of the tree if level < max_tree_depth
 *   max_box_ratio  <-- stop adding levels to tree when
 *                      (n_linked_boxes > max_box_ratio*n_initial_boxes)
 *---------------------------------------------------------------------------*/

void
fvmc_neighborhood_set_options(fvmc_neighborhood_t  *n,
                             int                  max_tree_depth,
                             int                  leaf_threshold,
                             int                  max_box_ratio)
{
  if (n == NULL)
    return;

  n->max_tree_depth = max_tree_depth;
  n->leaf_threshold = leaf_threshold;
  n->max_box_ratio = max_box_ratio;
}

/*----------------------------------------------------------------------------
 * Retrieve pointers to of arrays from a neighborhood_t structure.
 *
 * Arrays remain the property of the neighborhood_t structure, and must not
 * be modified by the caller.
 *
 * parameters:
 *   n              <-> pointer to fvmc_neighborhood_t structure.
 *   n_elts         --> number of elements with neighbors in block
 *                      associated with local rank
 *   elt_num        --> global element numbers in local block (size: n_elts)
 *   neighbor_index --> start index of neighbors (size: n_elts + 1)
 *   neighbor_num   --> global element neighbor numbers
 *                      (size: neighbor_index[n_elts])
 *----------------------------------------------------------------------------*/

void
fvmc_neighborhood_get_data(const fvmc_neighborhood_t         *n,
                          fvmc_lnum_t                       *n_elts,
                          fvmc_gnum_t                **const elt_num,
                          fvmc_lnum_t                **const neighbor_index,
                          fvmc_gnum_t                **const neighbor_num)
{
  if (n != NULL) {

    if (n_elts != NULL)
      *n_elts = n->n_elts;

    if (elt_num != NULL)
      *elt_num = n->elt_num;

    if (neighbor_index != NULL)
      *neighbor_index = n->neighbor_index;

    if (neighbor_num != NULL)
      *neighbor_num = n->neighbor_num;
  }
}

/*----------------------------------------------------------------------------
 * Transfer ownership of arrays from a neighborhood_t structure to the
 * calling program.
 *
 * Arrays that are transferred are removed from the structure, so its
 * use after calling this function is limited to querying timing information
 * until a new neighborhood is computed.
 *
 * parameters:
 *   n              <-> pointer to fvmc_neighborhood_t structure.
 *   n_elts         --> number of elements with neighbors in block
 *                      associated with local rank
 *   elt_num        --> global element numbers in local block (size: n_elts)
 *   neighbor_index --> start index of neighbors (size: n_elts + 1)
 *   neighbor_num   --> global element neighbor numbers
 *                      (size: neighbor_index[n_elts])
 *----------------------------------------------------------------------------*/

void
fvmc_neighborhood_transfer_data(fvmc_neighborhood_t   *n,
                               fvmc_lnum_t           *n_elts,
                               fvmc_gnum_t          **elt_num,
                               fvmc_lnum_t          **neighbor_index,
                               fvmc_gnum_t          **neighbor_num)
{
  if (n != NULL) {

    if (n_elts != NULL)
      *n_elts = n->n_elts;

    if (elt_num != NULL) {
      *elt_num = n->elt_num;
      n->elt_num = NULL;
    }
    if (neighbor_index != NULL) {
      *neighbor_index = n->neighbor_index;
      n->neighbor_index = NULL;
    }
    if (neighbor_num != NULL) {
      *neighbor_num = n->neighbor_num;
      n->neighbor_num = NULL;
    }
  }
}

/*----------------------------------------------------------------------------
 * Determine intersecting boxes.
 *
 * Box global numbers and extents may be either copied for the structure's
 * internal use from the caller, or tranferred to the neighborhood management
 * structure: both the box_gnum and extents arguments have an "assigned"
 * variant, in which cas a pointer to a pointer is provided, and the
 * argument's property is transferred to the neighborhod management structure.
 * The unused variant of an argument should be set to NULL.
 *
 * Boxes may be distributed among processors, so their intersections are
 * determined using a block distribution, and defined using their
 * global numbers.
 *
 * parameters:
 *   n                 <-> pointer to neighborhood management structure
 *   dim               <-- spatial dimension
 *   n_boxes           <-- local number of boxes
 *   box_gnum          <-- global numbering of boxes
 *   extents           <-- coordinate extents (size: n_boxes*dim*2, as
 *                         xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   box_gnum_transfer <-> as box_gnum, ownership transferred (NULL on return)
 *   extents_transfer  <-> as extents, ownership transferred (NULL on return)
 *   comm       <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
fvmc_neighborhood_by_boxes(fvmc_neighborhood_t  *n,
                          int                  dim,
                          fvmc_lnum_t           n_boxes,
                          const fvmc_gnum_t    *box_gnum,
                          const fvmc_coord_t   *extents,
                          fvmc_gnum_t         **box_gnum_assigned,
                          fvmc_coord_t        **extents_assigned)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvmc_box_tree_t  *bt = NULL;
  fvmc_box_set_t  *boxes = NULL;

  const fvmc_gnum_t  *_box_gnum = box_gnum;
  const fvmc_coord_t  *_extents = extents;

  int  n_ranks = 1;

  clock_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  /* Transfer data if necessary */

  if (box_gnum_assigned != NULL)
    _box_gnum = *box_gnum_assigned;
  if (extents_assigned != NULL)
    _extents = *extents_assigned;

  /* Reset structure if necessary */

  n->n_elts = 0;
  if (n->elt_num != NULL)
    BFTC_FREE(n->elt_num);
  if (n->neighbor_index != NULL)
    BFTC_FREE(n->neighbor_index);
  if (n->neighbor_num != NULL)
    BFTC_FREE(n->neighbor_num);

  /* Allocate fvmc_box_set_t structure and initialize it */

#if defined(FVMC_HAVE_MPI)

  if (n->comm != MPI_COMM_NULL)
    MPI_Comm_size(n->comm, &n_ranks);

  boxes = fvmc_box_set_create(dim,
                             1,  /* normalize */
                             1,  /*allow_projection */
                             n_boxes,
                             _box_gnum,
                             _extents,
                             n->comm);

  if (n_ranks > 1)
    _redistribute_boxes(n, boxes);

#else

  boxes = fvmc_box_set_create(dim,
                             1,  /* normalize */
                             1,  /*allow_projection */
                             n_boxes,
                             _box_gnum,
                             _extents);

#endif

  /* Free transferred data if applicable */

  if (box_gnum_assigned != NULL) {
    _box_gnum = NULL;
    BFTC_FREE(*box_gnum_assigned);
  }

  if (extents_assigned != NULL) {
    _extents = NULL;
    BFTC_FREE(*extents_assigned);
  }

  /* Build a tree structure and use it to order bounding boxes */

  /* Create and initialize a box tree structure */

  bt = fvmc_box_tree_create(n->max_tree_depth,
                           n->leaf_threshold,
                           n->max_box_ratio);

  /* Build a tree and put bounding boxes */

  fvmc_box_tree_set_boxes(bt,
                         boxes,
                         FVMC_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics((&n->bt_stats), bt);

  /* Update construction times. */

  clock_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  n->cpu_time[0] = cpu_end - cpu_start;
  n->wtime[0] = clock_end - clock_start;

  clock_start = clock_end;
  cpu_start = cpu_end;

  /* Allocate structure to store intersections between boxes */

  n->n_elts = fvmc_box_set_get_size(boxes);

  BFTC_MALLOC(n->elt_num, n->n_elts, fvmc_gnum_t);
  memcpy(n->elt_num,
         fvmc_box_set_get_g_num(boxes),
         n->n_elts*sizeof(fvmc_gnum_t));

  fvmc_box_tree_get_intersects(bt,
                              boxes,
                              &(n->neighbor_index),
                              &(n->neighbor_num));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  fvmc_box_tree_dump(bt);
  fvmc_box_set_dump(boxes, 1);
#endif

  /* Destroy the associated box tree */

  fvmc_box_tree_destroy(&bt);

  /* Compact intersections list, delete redundancies and order intersections */

  _remove_nums_without_neighbors(n);
  _order_neighborhood(n);

#if defined(FVMC_HAVE_MPI)

  /* Synchronize list of intersections for each element of the list
     and distribute it by block over the ranks */

  if (n_ranks > 1)
    _sync_by_block(n, fvmc_box_set_get_global_size(boxes));

#endif /* HAVE_MPI */

  /* Destroy the box set structures */

  fvmc_box_set_destroy(&boxes);

  _clean_neighbor_nums(n);

  /* Update query times. */

  clock_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  n->cpu_time[1] = cpu_end - cpu_start;
  n->wtime[1] = clock_end - clock_start;
}

/*----------------------------------------------------------------------------
 * Get global statistics relative to the search structures used
 * by fvmc_neighborhood_by_boxes().
 *
 * All fields returned are optional: if their argument is set to NULL,
 * the corresponding information will not be returned.
 *
 * For each field not set to NULL, 3 values are always returned:
 * the mean on all ranks (rounded to the closest integer), the minimum,
 * and the maximum value respectively.
 *
 * In serial mode, the mean, minimum, and maximum will be identical for most
 * fields, but all 3 values are returned nonetheless.
 *
 * Note that the final memory use is only relative to the final search
 * structure, and the comparison of the total (or average) with the minima
 * and maxima may give an indication on load balancing.

 * The mem_required field only takes into account the theoretical maximum
 * memory footprint of the main search structure during its construction phase,
 * and of that of temporary structures before load balancing in parallel mode,
 * but does not include minor work arrays or buffers used during the algorithm.
 *
 * Neither of the 2 memory fields include the footprint of the arrays
 * containing the query results.
 *
 * parameters:
 *   n                  <-- pointer to neighborhood management structure
 *   dim                --> layout dimension (3, 2, or 1)
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_final          --> theoretical memory for final search structure
 *   mem_required       --> theoretical maximum memory for main structures
 *                          used during the algorithm
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
fvmc_neighborhood_get_box_stats(const fvmc_neighborhood_t  *n,
                               int                        depth[3],
                               fvmc_lnum_t                 n_leaves[3],
                               fvmc_lnum_t                 n_boxes[3],
                               fvmc_lnum_t                 n_threshold_leaves[3],
                               fvmc_lnum_t                 n_leaf_boxes[3],
                               size_t                     mem_final[3],
                               size_t                     mem_required[3])
{
  size_t i;

  if (n == NULL)
    return 0;

  for (i = 0; i < 3; i++) {

    if (depth != NULL)
      depth[i] = n->bt_stats.depth[i];

    if (n_leaves != NULL)
      n_leaves[i] = n->bt_stats.n_leaves[i];

    if (n_boxes != NULL)
      n_boxes[i] = n->bt_stats.n_boxes[i];

    if (n_threshold_leaves != NULL)
      n_threshold_leaves[i] = n->bt_stats.n_threshold_leaves[i];

    if (n_leaf_boxes != NULL)
      n_leaf_boxes[i] = n->bt_stats.n_leaf_boxes[i];

    if (mem_final != NULL)
      mem_final[i] = n->bt_stats.mem_used[i];

    if (mem_required != NULL)
      mem_required[i] = n->bt_stats.mem_required[i];
  }
  return n->bt_stats.dim;
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * parameters:
 *   n              <-- pointer to neighborhood management structure
 *   build_wtime    --> initialization Wall-clock time (or NULL)
 *   build_cpu_time --> initialization CPU time (or NULL)
 *   query_wtime    --> query Wall-clock time (or NULL)
 *   query_cpu_time --> query CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvmc_neighborhood_get_times(const fvmc_neighborhood_t  *n,
                           double                    *build_wtime,
                           double                    *build_cpu_time,
                           double                    *query_wtime,
                           double                    *query_cpu_time)
{
  if (n == NULL)
    return;

  if (build_wtime != NULL)
    *build_wtime = n->wtime[0];
  if (build_cpu_time != NULL)
    *build_cpu_time = n->cpu_time[0];

  if (query_wtime != NULL)
    *query_wtime = n->wtime[1];
  if (query_cpu_time != NULL)
    *query_cpu_time = n->cpu_time[1];
}

/*----------------------------------------------------------------------------
 * Dump a neighborhood management structure.
 *
 * parameters:
 *   n <-- pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

void
fvmc_neighborhood_dump(const fvmc_neighborhood_t  *n)
{
  fvmc_lnum_t  i, j;

  bftc_printf("\n"
             "Neighborhood information: %p\n\n", n);

  if (n == NULL)
    return;

  bftc_printf("number of elements: %10d\n"
             "list size:          %10d\n\n",
             (int)(n->n_elts), (int)(n->neighbor_index[n->n_elts]));

  bftc_printf("max tree depth:     %d\n"
             "leaf threshold:     %d\n"
             "max box ratio       %d\n\n",
             n->max_tree_depth, n->leaf_threshold, n->max_box_ratio);

#if defined(FVMC_HAVE_MPI)
  if (n->comm != MPI_COMM_NULL)
    bftc_printf("\n"
               "Associated MPI communicator: %ld\n",
               (long)(n->comm));
#endif

  bftc_printf("CPU time:           %f\n"
             "Wall-clock time:    %f\n\n",
             n->cpu_time, n->wtime);

  for (i = 0; i < n->n_elts; i++) {

    int n_neighbors = (n->neighbor_index[i+1] - n->neighbor_index[i]);

    bftc_printf("global num.: %10u | n_neighbors : %3d |",
               n->elt_num[i], n_neighbors);

    for (j = n->neighbor_index[i]; j < n->neighbor_index[i+1]; j++)
      bftc_printf("  %10u ", n->neighbor_num[j]);
    bftc_printf("\n");

  }

  bftc_printf_flush();
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
