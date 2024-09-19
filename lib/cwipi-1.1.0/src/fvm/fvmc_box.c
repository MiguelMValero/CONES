/*============================================================================
 * Handle boxes aligned with Cartesian axes.
 *===========================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

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
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_timer.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "fvmc_box.h"
#include "fvmc_box_priv.h"

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type and structure definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Display a histogram on leaves associated to the boxes and several
 * other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib          <-- pointer to the fvmc_box_distrib_t structure
 *   n_quantiles      <-> number of quantiles requested (or NULL);
 *                        may be lower upon return
 *   quantile_start   --> start of quantile (size: n_quantiles + 1)
 *   n_quantile_boxes --> number of boxes per quantile (size: n_quantiles)
 *   imbalance        --> distribution imbalance measure (or NULL)
 *   n_ranks          --> number of ranks with boxes (or NULL)
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static void
_get_distrib_statistics(const fvmc_box_distrib_t  *distrib,
                        fvmc_lnum_t               *n_quantiles,
                        fvmc_lnum_t                quantile_start[],
                        fvmc_lnum_t                n_quantile_boxes[],
                        double                   *imbalance,
                        int                      *n_ranks,
                        MPI_Comm                  comm)
{
  fvmc_lnum_t  i, j, k, step, delta, _n_rank_boxes;

  int  _n_ranks = 0;
  fvmc_lnum_t  _min = INT_MAX, _max = 0, gmin = 0, gmax = 0;

  /* Sanity checks */

  assert(distrib != NULL);
  assert(distrib->index != NULL);

  if (n_quantiles != NULL) {

    fvmc_lnum_t _n_quantiles = 0;

    /* Get global min and max number of boxes */

    for (i = 0; i < distrib->n_ranks; i++) {

      _n_rank_boxes = distrib->index[i+1] - distrib->index[i];
      _min = FVMC_MIN(_min, _n_rank_boxes);
      _max = FVMC_MAX(_max, _n_rank_boxes);

      if (_n_rank_boxes > 0)
        _n_ranks += 1;

    }

    gmin = _min;
    gmax = _max;

    MPI_Allreduce(&_min, &gmin, 1, FVMC_MPI_LNUM, MPI_MIN, comm);
    MPI_Allreduce(&_max, &gmax, 1, FVMC_MPI_LNUM, MPI_MAX, comm);

    /* Build a histogram for the distribution of boxes */

    delta = gmax - gmin;
    if (delta < _n_quantiles)
      _n_quantiles = delta;

    /* Define quantiles */

    step = delta / _n_quantiles;
    if (delta % _n_quantiles > 0)
      step++;

    for (i = 0; i < _n_quantiles; i++)
      quantile_start[i] = gmin + i*step;

    quantile_start[_n_quantiles] = gmax + 1;

    /* Count for quantiles */

    for (j = 0; j < _n_quantiles; j++)
      n_quantile_boxes[j] = 0;

    if (delta > 0) {  /* Loop on boxes */

      for (i = 0; i < distrib->n_ranks; i++) {

        _n_rank_boxes = distrib->index[i+1] - distrib->index[i];

        for (k = 1; k < _n_quantiles; k++) {
          if (_n_rank_boxes < gmin + k*step)
            break;
        }
        n_quantile_boxes[k-1] += 1;

      } /* End of loop on boxes */

    }

    *n_quantiles = _n_quantiles;
  }

  /* Set other return values */

  if (imbalance != NULL)
    *imbalance = distrib->fit;

  if (n_ranks != NULL)
    *n_ranks = _n_ranks;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create a set of boxes and initialize it.
 *
 * parameters:
 *   dim              <-- spatial dimension
 *   normalize        <-- 1 if boxes are to be normalized, 0 otherwize
 *   allow_projection <-- if 1, project to lower dimension if all boxes
 *                        are cut by the median plane of the set.
 *   n_boxes          <-- number of elements to create
 *   box_gnum         <-- global numbering of boxes
 *   extents          <-- coordinate extents (size: n_boxes*dim*2, as
 *                        xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   a new allocated pointer to a fvmc_box_set_t structure.
 *---------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)
fvmc_box_set_t *
fvmc_box_set_create(int                 dim,
                   int                 normalize,
                   int                 allow_projection,
                   fvmc_lnum_t          n_boxes,
                   const fvmc_gnum_t   *box_gnum,
                   const fvmc_coord_t  *box_extents,
                   MPI_Comm            comm)
#else
fvmc_box_set_t *
fvmc_box_set_create(int                 dim,
                   int                 normalize,
                   int                 allow_projection,
                   fvmc_lnum_t          n_boxes,
                   const fvmc_gnum_t   *box_gnum,
                   const fvmc_coord_t  *box_extents)
#endif
{
  int j, k;
  fvmc_lnum_t  i;
  fvmc_gnum_t  n_g_boxes = n_boxes;
  fvmc_coord_t  g_min[3], g_max[3], g_extents[6];

  fvmc_box_set_t  *boxes = NULL;

  /* Get global min/max coordinates */

#if defined(FVMC_HAVE_MPI)
  fvmc_morton_get_global_extents(dim, n_boxes, box_extents, g_extents, comm);
#else
  fvmc_morton_get_global_extents(dim, n_boxes, box_extents, g_extents);
#endif

  for (j = 0; j < 3; j++) {
    g_min[j] = g_extents[j];
    g_max[j] = g_extents[j+dim];
  }

#if defined(FVMC_HAVE_MPI)

  if (comm != MPI_COMM_NULL) {

    fvmc_gnum_t  box_max = 0;

    for (i = 0; i < n_boxes; i++)
      box_max = FVMC_MAX(box_max, box_gnum[i]);

    MPI_Allreduce(&box_max, &n_g_boxes, 1, FVMC_MPI_GNUM, MPI_MAX, comm);

  }

#endif

  /* Allocate box set structure and initialize it */

  BFTC_MALLOC(boxes, 1, fvmc_box_set_t);

  boxes->dim = dim;
  boxes->n_boxes = n_boxes;
  boxes->n_g_boxes = n_g_boxes;

  for (j = 0; j < 3; j++) {
    boxes->dimensions[j] = j;
    boxes->gmin[j] = g_min[j];
    boxes->gmax[j] = g_max[j];
  }

  boxes->g_num = NULL;
  boxes->extents = NULL;

#if defined(FVMC_HAVE_MPI)
  boxes->comm = comm;
#endif

  /* Optionally allow and detect a layout of lower
     dimension than the spatial dimension */

  if (allow_projection) {

    double g_mid[3];
    int proj[] = {1, 1, 1};

    for (j = 0; j < dim; j++)
      g_mid[j] = (g_min[j] + g_max[j]) * 0.5;

    for (i = 0; i < n_boxes; i++) {
      for (j = 0; j < dim; j++) {
        if (   box_extents[i*dim*2 + j]     > g_mid[j]
            || box_extents[i*dim*2 + j+dim] < g_mid[j])
          proj[j] = 0;
      }
    }

#if defined(FVMC_HAVE_MPI)
    if (comm != MPI_COMM_NULL) {
      int l_proj[3];
      for (j = 0; j < dim; j++)
        l_proj[j] = proj[j];
      MPI_Allreduce(l_proj, proj, dim, MPI_INT, MPI_MIN, comm);
    }
#endif

    boxes->dim = 0;
    for (j = 0; j < dim; j++) {
      if (proj[j] == 0) {
        boxes->dimensions[boxes->dim] = j;
        boxes->dim += 1;
      }
    }


  }

  for (j = boxes->dim; j < 3; j++) /* ensure all is initialized */
    boxes->dimensions[j] = -1;

  /* Now assign values */

  BFTC_MALLOC(boxes->g_num, n_boxes, fvmc_gnum_t);
  BFTC_MALLOC(boxes->extents, n_boxes*boxes->dim*2, fvmc_coord_t);

  for (i = 0; i < n_boxes; i++) {

    fvmc_coord_t *_min = boxes->extents + (boxes->dim*2*i);
    fvmc_coord_t *_max = _min + boxes->dim;

    boxes->g_num[i] = box_gnum[i];

    for (j = 0; j < boxes->dim; j++) {
      k = boxes->dimensions[j];
      _min[j] = box_extents[i*dim*2 + k];
      _max[j] = box_extents[i*dim*2 + k+dim];
      assert(_min[j] <= _max[j]);
    }
  }

  /* Define the normalized min/max coordinates of the box */

  if (normalize) {

    fvmc_coord_t  d[3], s[3];

    for (j = 0; j < boxes->dim; j++) {
      k = boxes->dimensions[j];
      s[j] = g_min[k];
      d[j] = g_max[k] - g_min[k];
    }

    for (i = 0; i < n_boxes; i++) {

      fvmc_coord_t *_min = boxes->extents + (boxes->dim*2*i);
      fvmc_coord_t *_max = _min + boxes->dim;

      for (j = 0; j < boxes->dim; j++) {
        _min[j] = (_min[j] - s[j]) / d[j];
        _max[j] = (_max[j] - s[j]) / d[j];
      }
    }

  }

  /* Return pointer to structure */

  return boxes;
}

/*----------------------------------------------------------------------------
 * Delete a fvmc_box_set_t structure.
 *
 * parameters:
 *   boxes <-> pointer to the fvmc_box_set_t structure to delete
 *---------------------------------------------------------------------------*/

void
fvmc_box_set_destroy(fvmc_box_set_t  **boxes)
{
  if (boxes != NULL) {

    fvmc_box_set_t  *_boxes = *boxes;

    if (_boxes == NULL)
      return;

    BFTC_FREE(_boxes->g_num);
    BFTC_FREE(_boxes->extents);
    BFTC_FREE(_boxes);
  }
}

/*----------------------------------------------------------------------------
 * Return the dimension associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   associated spatial dimension
 *---------------------------------------------------------------------------*/

int
fvmc_box_set_get_dim(const fvmc_box_set_t  *boxes)
{
  int retval = 0;

  if (boxes != NULL)
    retval = boxes->dim;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the local number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_box_set_get_size(const fvmc_box_set_t  *boxes)
{
  fvmc_lnum_t retval = 0;

  if (boxes != NULL)
    retval = boxes->n_boxes;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the global number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   global number of boxes
 *---------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_box_set_get_global_size(const fvmc_box_set_t  *boxes)
{
  fvmc_gnum_t retval = 0;

  if (boxes != NULL)
    retval = boxes->n_g_boxes;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return extents associated with a set of boxes.
 *
 * The extents array is organized in the following fashion:
 * {x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
 *  x_min_n, y_min_n, ..., x_max_n, y_max_n, ...}
 *
 * Its size is thus: n_boxes * dim * 2.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to extents array
 *---------------------------------------------------------------------------*/

const fvmc_coord_t *
fvmc_box_set_get_extents(fvmc_box_set_t  *boxes)
{
  assert(boxes != NULL);

  return boxes->extents;
}

/*----------------------------------------------------------------------------
 * Return global numbers associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to global box numbers array
 *---------------------------------------------------------------------------*/

const fvmc_gnum_t *
fvmc_box_set_get_g_num(fvmc_box_set_t  *boxes)
{
  assert(boxes != NULL);

  return boxes->g_num;
}

/*----------------------------------------------------------------------------
 * Build a Morton_index to get a well-balanced distribution of the boxes.
 *
 * parameters:
 *  boxes      <-- pointer to associated fvmc_box_set_t structure
 *  distrib    <-> pointer to a fvmc_box_distrib_t structure
 *  n_leaves   <-- number of leaves with weight > 0
 *  leaf_codes <-- Morton code for each leaf
 *  weight     <-- number of boxes related to each leaf
 *---------------------------------------------------------------------------*/

void
fvmc_box_set_build_morton_index(const fvmc_box_set_t  *boxes,
                               fvmc_box_distrib_t    *distrib,
                               fvmc_lnum_t            n_leaves,
                               fvmc_morton_code_t    *leaf_codes,
                               fvmc_lnum_t           *weight)
{
#if defined(FVMC_HAVE_MPI)

  fvmc_lnum_t  *order = NULL;

  assert(distrib != NULL);
  assert(distrib->morton_index != NULL);

  BFTC_MALLOC(order, n_leaves, fvmc_lnum_t);

  /* Locally order Morton encoding */

  fvmc_morton_local_order(n_leaves,
                         leaf_codes,
                         order);

  /* Compute a Morton index on ranks and return the associated fit */

  if (boxes->comm != MPI_COMM_NULL)
    distrib->fit = fvmc_morton_build_rank_index(boxes->dim,
                                               distrib->max_level,
                                               n_leaves,
                                               leaf_codes,
                                               weight,
                                               order,
                                               distrib->morton_index,
                                               boxes->comm);
  /* Free memory */

  BFTC_FREE(order);

#endif
}

/*----------------------------------------------------------------------------
 * Redistribute boxes over the ranks according to the Morton index to
 * assume a better balanced distribution of the boxes.
 *
 * parameters:
 *  distrib <--  data structure on box distribution
 *  boxes   <->  pointer to the structure to redistribute
 *---------------------------------------------------------------------------*/

void
fvmc_box_set_redistribute(const fvmc_box_distrib_t  *distrib,
                         fvmc_box_set_t            *boxes)
{
#if defined(FVMC_HAVE_MPI)

  int  rank_id;

  fvmc_lnum_t  i, j;
  int  *send_count = NULL, *send_shift = NULL;
  int  *recv_count = NULL, *recv_shift = NULL;
  fvmc_gnum_t  *send_g_num = NULL;
  fvmc_coord_t  *send_extents = NULL;

  const int stride = boxes->dim * 2;

  /* Sanity checks */

  assert(distrib != NULL);
  assert(boxes != NULL);
  assert(distrib->n_ranks > 1);

  /* Build send_buf, send_count and send_shift
     to build a rank to boxes indexed list */

  BFTC_MALLOC(send_count, distrib->n_ranks, int);
  BFTC_MALLOC(recv_count, distrib->n_ranks, int);
  BFTC_MALLOC(send_shift, distrib->n_ranks + 1, int);
  BFTC_MALLOC(recv_shift, distrib->n_ranks + 1, int);

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    send_count[rank_id]
      = distrib->index[rank_id+1] - distrib->index[rank_id];

  /* Exchange number of boxes to send to each process */

  MPI_Alltoall(send_count, 1, MPI_INT,
               recv_count, 1, MPI_INT, boxes->comm);

  for (i = 0; i < distrib->n_ranks; i++)
    send_shift[i] = distrib->index[i];

  recv_shift[0] = 0;
  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];

  /* Build send_buffers */

  BFTC_MALLOC(send_g_num, distrib->index[distrib->n_ranks], fvmc_gnum_t);
  BFTC_MALLOC(send_extents,
             distrib->index[distrib->n_ranks] * boxes->dim * 2,
             fvmc_coord_t);

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++)
    send_count[rank_id] = 0;

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++) {

    for (i = distrib->index[rank_id];
         i < distrib->index[rank_id+1];
         i++) {

      fvmc_lnum_t  box_id = distrib->list[i];
      fvmc_lnum_t  shift = distrib->index[rank_id] + send_count[rank_id];

      send_g_num[shift] = boxes->g_num[box_id];

      for (j = 0; j < stride; j++)
        send_extents[shift*stride + j] = boxes->extents[box_id*stride + j];

      send_count[rank_id] += 1;

    }

  } /* End of loop on ranks */

  /* Prepare to replace the local arrays */

  boxes->n_boxes = recv_shift[distrib->n_ranks];
  BFTC_FREE(boxes->g_num);
  BFTC_FREE(boxes->extents);

  BFTC_MALLOC(boxes->g_num, boxes->n_boxes, fvmc_gnum_t);
  BFTC_MALLOC(boxes->extents, boxes->n_boxes*stride, fvmc_coord_t);

  /* Exchange boxes between processes */

  MPI_Alltoallv(send_g_num, send_count, send_shift, FVMC_MPI_GNUM,
                boxes->g_num, recv_count, recv_shift, FVMC_MPI_GNUM,
                boxes->comm);

  for (rank_id = 0; rank_id < distrib->n_ranks; rank_id++) {
    send_count[rank_id] *= stride;
    send_shift[rank_id] *= stride;
    recv_count[rank_id] *= stride;
    recv_shift[rank_id] *= stride;
  }

  MPI_Alltoallv(send_extents, send_count, send_shift, FVMC_MPI_COORD,
                boxes->extents, recv_count, recv_shift, FVMC_MPI_COORD,
                boxes->comm);

  /* Free buffers */

  BFTC_FREE(send_g_num);
  BFTC_FREE(send_extents);
  BFTC_FREE(send_count);
  BFTC_FREE(send_shift);
  BFTC_FREE(recv_count);
  BFTC_FREE(recv_shift);

#endif /* FVMC_HAVE_MPI */
}

/*----------------------------------------------------------------------------
 * Dump a fvmc_box_set_t structure.
 *
 * parameters:
 *   boxes     <-- pointer to the fvmc_box_t structure
 *   verbosity <-- verbosity level (0 or 1)
 *----------------------------------------------------------------------------*/

void
fvmc_box_set_dump(const fvmc_box_set_t  *boxes,
                 int                   verbosity)
{
  fvmc_lnum_t  i;

  const char  XYZ[4] = "XYZ";

  if (boxes == NULL)
    return;

  /* Print basic information */

  if (boxes->dim == 3)
    bftc_printf("\nBox set (3D layout):\n\n"
               "global min/max on selected faces:\n"
               "  [%7.5e %7.5e %7.5e] --> [%7.5e %7.5e %7.5e]\n",
               boxes->gmin[0], boxes->gmin[1], boxes->gmin[2],
               boxes->gmax[0], boxes->gmax[1], boxes->gmax[2]);

  else if (boxes->dim == 2) {
    bftc_printf("\nBox set (2D layout, selected axes [%c, %c]\n\n",
               XYZ[boxes->dimensions[0]],
               XYZ[boxes->dimensions[1]]);
    bftc_printf("global min/max on selected faces:\n"
               "  [%7.5e %7.5e] --> [%7.5e %7.5e]\n",
               boxes->gmin[boxes->dimensions[0]],
               boxes->gmin[boxes->dimensions[1]],
               boxes->gmax[boxes->dimensions[0]],
               boxes->gmax[boxes->dimensions[1]]);
  }

  else if (boxes->dim == 1) {
    bftc_printf("\nBox set (1D layout, selected axis [%c]\n\n",
               XYZ[boxes->dimensions[0]]);
    bftc_printf("global min/max on selected faces:\n"
               "  [%7.5e %7.5e] --> [%7.5e %7.5e]\n",
               boxes->gmin[boxes->dimensions[0]],
               boxes->gmin[boxes->dimensions[1]],
               boxes->gmax[boxes->dimensions[0]],
               boxes->gmax[boxes->dimensions[1]]);
  }
  bftc_printf_flush();

  /* Print detailed box information */

  if (verbosity < 1)
    return;

  if (boxes->dim == 3) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const fvmc_coord_t *bmin = boxes->extents + i*6;
      const fvmc_coord_t *bmax = boxes->extents + i*6 + 3;
      bftc_printf("  id %8d, num %9lu: "
                 "[%7.5e %7.5e %7.5e] --> [%7.5e %7.5e %7.5e]\n",
                 i, (unsigned long)(boxes->g_num[i]),
                 bmin[0], bmin[1], bmin[2],
                 bmax[0], bmax[1], bmax[2]);
    }
  }

  else if (boxes->dim == 2) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const fvmc_coord_t *bmin = boxes->extents + i*4;
      const fvmc_coord_t *bmax = boxes->extents + i*4 + 2;
      bftc_printf("  id %8d, num %9lu: "
                 "[%7.5e %7.5e] --> [%7.5e %7.5e]\n",
                 i, (unsigned long)(boxes->g_num),
                 bmin[0], bmin[1], bmax[0], bmax[1]);
    }
  }

  else if (boxes->dim == 1) {
    for (i = 0; i < boxes->n_boxes; i++) {
      const fvmc_coord_t *bmin = boxes->extents + i*2;
      const fvmc_coord_t *bmax = boxes->extents + i*2 + 1;
      bftc_printf("  id %8d, num %9lu: "
                 "[%7.5e] --> [%7.5e]\n",
                 i, (unsigned long)(boxes->g_num),
                 bmin[0], bmax[0]);
    }
  }

  /* Sanity check */

  for (i = 0; i < boxes->n_boxes; i++) {
    int j;
    const fvmc_coord_t *bmin = boxes->extents + boxes->dim*2*i;
    const fvmc_coord_t *bmax = boxes->extents + boxes->dim*(2*i + 1);
    for (j = 0; j < boxes->dim; j++) {
      if (bmin[j] > bmax[j])
        bftc_error(__FILE__, __LINE__, 0,
                  _("Inconsistent box found (min > max):\n"
                    "  global number:  %u\n"
                    "  min       :  %10.4g\n"
                    "  max       :  %10.4g\n"),
                  boxes->g_num[i], bmin[j], bmax[j]);
    }
  }

}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a fvmc_box_distrib_t structure.
 *
 * parameters:
 *   n_boxes   <-- number of boxes
 *   n_g_boxes <-- global number of boxes
 *   max_level <-- max level reached locally in the related tree
 *   comm      <-- MPI communicator. on which the distribution takes place
 *
 * returns:
 *   a pointer to a new allocated fvmc_box_distrib_t structure.
 *---------------------------------------------------------------------------*/

fvmc_box_distrib_t *
fvmc_box_distrib_create(fvmc_lnum_t  n_boxes,
                       fvmc_gnum_t  n_g_boxes,
                       int         max_level,
                       MPI_Comm    comm)
{
  int  i, n_ranks, gmax_level;

  fvmc_box_distrib_t  *new_distrib = NULL;

  if (n_g_boxes == 0)
    return NULL;

  BFTC_MALLOC(new_distrib, 1, fvmc_box_distrib_t);

  /* Parallel parameters */

  MPI_Comm_size(comm, &n_ranks);

  new_distrib->n_ranks = n_ranks;

  new_distrib->n_boxes = n_boxes;

  assert(n_ranks > 1);

  BFTC_MALLOC(new_distrib->morton_index, n_ranks + 1, fvmc_morton_code_t);

  MPI_Allreduce(&max_level, &gmax_level, 1, MPI_INT, MPI_MAX, comm);

  new_distrib->max_level = gmax_level;
  new_distrib->fit = 999.0;

  BFTC_MALLOC(new_distrib->index, n_ranks + 1, fvmc_lnum_t);

  for (i = 0; i < n_ranks + 1; i++)
    new_distrib->index[i] = 0;

  new_distrib->list = NULL;

  return  new_distrib;
}

/*----------------------------------------------------------------------------
 * Destroy a fvmc_box_distrib_t structure.
 *
 * parameters:
 *   distrib <-> pointer to pointer to the structure to destroy
 *---------------------------------------------------------------------------*/

void
fvmc_box_distrib_destroy(fvmc_box_distrib_t  **distrib)
{
  if (distrib != NULL) {

    fvmc_box_distrib_t  *d = *distrib;

    if (d == NULL)
      return;

    BFTC_FREE(d->index);
    BFTC_FREE(d->list);
    BFTC_FREE(d->morton_index);

    BFTC_FREE(d);
  }
}

/*----------------------------------------------------------------------------
 * Delete redundancies in box distribution
 *
 * parameters:
 *   distrib <->  pointer to the fvmc_box_distrib_t structure
 *---------------------------------------------------------------------------*/

void
fvmc_box_distrib_clean(fvmc_box_distrib_t  *distrib)
{
  int  i, rank;

  fvmc_lnum_t  *counter = NULL, *new_index = NULL;

  BFTC_MALLOC(counter, distrib->n_boxes, fvmc_lnum_t);
  BFTC_MALLOC(new_index, distrib->n_ranks + 1, fvmc_lnum_t);

  for (i = 0; i < distrib->n_ranks + 1; i++)
    new_index[i] = 0;

  for (rank = 0; rank < distrib->n_ranks; rank++) {

    fvmc_lnum_t  shift = new_index[rank];
    fvmc_lnum_t  start = distrib->index[rank];
    fvmc_lnum_t  end = distrib->index[rank + 1];

    if (end - start > 0) {

      for (i = 0; i < distrib->n_boxes; i++)
        counter[i] = 0;

      for (i = start; i < end; i++)
        counter[distrib->list[i]] += 1;

      for (i = 0; i < distrib->n_boxes; i++) {

        if (counter[i] > 0)
          distrib->list[shift++] = i;

      }

    } /* end if end - start > 0 */

    new_index[rank+1] = shift;

  } /* End of loop on ranks */

  /* Memory management */

  BFTC_FREE(distrib->index);
  BFTC_REALLOC(distrib->list, new_index[distrib->n_ranks], fvmc_lnum_t);

  distrib->index = new_index;

  BFTC_FREE(counter);
}

/*----------------------------------------------------------------------------
 * Display a histogram on leaves associated to the boxes and several
 * other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib <-- pointer to the fvmc_box_distrib_t structure
 *   comm    <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
fvmc_box_distrib_dump_statistics(const fvmc_box_distrib_t  *distrib,
                                MPI_Comm                  comm)
{
  fvmc_lnum_t  i;

  int          n_ranks = 0;
  fvmc_lnum_t   n_quantiles = 5;
  fvmc_lnum_t   quantile_start[6];
  fvmc_lnum_t   n_boxes[5];

  /* Sanity checks */

  assert(distrib != NULL);
  assert(distrib->index != NULL);

  _get_distrib_statistics(distrib,
                          &n_quantiles,
                          quantile_start,
                          n_boxes,
                          NULL,
                          &n_ranks,
                          comm);

  bftc_printf("\n"
             "- Box distribution statistics -\n\n");

  bftc_printf("   Distribution imbalance:              %10.4g\n",
             distrib->fit);
  bftc_printf("   Number of ranks in distribution:     %8d\n\n",
             n_ranks);

  /* Print histogram to show the distribution of boxes */

  if (n_quantiles > 0) {

    for (i = 0; i < n_quantiles - 1; i++)
      bftc_printf("    %3d : [ %10d ; %10d [ = %10d\n",
                 i+1, quantile_start[i], quantile_start[i+1], n_boxes[i]);

    i = n_quantiles -1;
    bftc_printf("    %3d : [ %10d ; %10d ] = %10d\n",
               i+1, quantile_start[i], quantile_start[i+1] - 1, n_boxes[i]);

  }
  bftc_printf_flush();
}

#endif /* defined(FVMC_HAVE_MPI) */

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
