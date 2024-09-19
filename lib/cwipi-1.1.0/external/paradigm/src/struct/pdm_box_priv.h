#ifndef __PDM_BOX_PRIV_H__
#define __PDM_BOX_PRIV_H__

/*============================================================================
 * Handle bounding boxes.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_morton.h"

#include "pdm_box.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/* Structure use to manage box distribution on a tree structure */

struct _PDM_box_distrib_t {

  int                 n_ranks;      /* Number of associated ranks */

  int                 n_boxes;      /* Number of bounding boxes */

  int                 max_level;    /* Global max level used to compute the
                                       distribution */
  double              fit;          /* Evaluation of the distribution
                                       (lower is better) */

  /* Morton code array defining an index on ranks = resulting distribution */

  PDM_morton_code_t  *morton_index; /* size = n_ranks + 1 */

  /* Indexed list on ranks to list related bounding boxes */

  int   *index;   /* Index on ranks (size = n_ranks + 1) */
  int   *list;    /* List of bounding boxes associated to each rank */

};



struct _PDM_boxes_t {
  int          n_boxes;        /* Number of bounding boxes */
  PDM_g_num_t *g_num;          /* Array of associated global numbers */
  double      *extents;        /* Extents associated with each box:
                                  * x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
                                  * x_min_n, y_min_n, ..., x_max_n, y_max_n,
                                  * (size: n_boxes * dim * 2) */

  int          n_part_orig;
  int         *n_boxes_orig;   /* Number of origin bounding boxes */
  int         *origin;         /* Initial location :
                                  * iproc, i_part, local_num */
};

typedef struct {

  PDM_mpi_win_shared_t *w_g_num;
  PDM_mpi_win_shared_t *w_extents;
  PDM_mpi_win_shared_t *w_n_boxes_orig;
  PDM_mpi_win_shared_t *w_origin;

} _w_boxes_data_t;


/* Set of bounding boxes */

struct _PDM_box_set_t {

  PDM_MPI_Comm    comm;                 /* Associated MPI communicator */
  int             dim;                  /* Spatial dimension (1, 2 or 3) */
  int             dimensions[3];        /* Only used in 1 or 2D: X = 0, Y = 1, Z = 2 */

  PDM_g_num_t  n_g_boxes;               /* Global number of bounding boxes */

  PDM_boxes_t *local_boxes;             /* Local boxes */

  int n_copied_ranks;                   /* Number of copies from other ranks */
  int *copied_ranks;                    /* Copied ranks */
  PDM_boxes_t *rank_boxes;              /* Boxes copied from other ranks */

  double      gmin[3];                  /* Global minima of the coordinates */
  double      gmax[3];                  /* Global maxima of the coordinates */

  int         normalized;               /* 1 if normalized, 0 otherwise */

  double      s[3];                     /* Translation for the normalization */
  double      d[3];                     /* Dilatation for the normalization */

  int                  n_rank_in_shm;
  PDM_MPI_Comm         comm_shared;
  PDM_boxes_t         *shm_boxes;       /* Boxes shared from other ranks */
  _w_boxes_data_t     *wboxes_data;

};

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_BOX_PRIV_H__ */
