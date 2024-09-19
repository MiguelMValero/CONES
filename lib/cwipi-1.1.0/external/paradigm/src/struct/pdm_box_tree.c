/*============================================================================
 * Search octrees and quadtrees of boxes.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_box_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_morton.h"
#include "pdm_vtk.h"
#include "pdm_plane.h"
#include "pdm_mesh_nodal.h"
#include "pdm_linear_programming.h"
#include "pdm_binary_search.h"
#include "pdm_unique.h"

#include "pdm_point_tree_seq.h"
#include "pdm_point_tree_seq_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_box_tree.h"
#include "pdm_box_tree_priv.h"
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro and Type definitions
 *============================================================================*/

#define PDM_BOX_TREE_MAX_BUILD_LOOPS 50

/*=============================================================================
 * Static global variables
 *============================================================================*/

static int iappel = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Build minimal octree between two octants
 *
 * \param [in]     octree    Current octree
 * \param [in]     code      Morton code
 * \param [inout]  extents   Extents associated to the Morton code
 *
 */

static void
_extents
(
 const int dim,
 PDM_morton_code_t code,
 double    extents[]
 )
{
  double side = 1. / (double) (1 << code.L);
  for (int i = 0; i < dim; i++) {
    extents[i]       = (double) code.X[i] * side;
    extents[dim + i] = extents[i] + side;
  }
}


static void
_extents_real
(
 const int               dim,
       PDM_morton_code_t code,
       double            s[],
       double            d[],
       double            extents[]
 )
{
  double side = 1. / (double) (1 << code.L);
  for (int i = 0; i < dim; i++) {
    extents[i]       = s[i] +  (double) code.X[i] * side * d[i];
    extents[dim + i] = extents[i] + side * d[i];
  }
}


/**
 *
 * \brief Compute distance to a box
 *
 * \param [in]   dim        Dimension
 * \param [in]   extents    Box extents
 * \param [in]   coords     Point coords
 * \param [out]  min_dist2  Square of minimum distance
 * \param [out]  max_dist2  Sqaure of maximum distance
 *
 * \return 1 if point is in the box, 0 otherwise
 *
 */

inline static int
_box_dist2_max
(
 const int              dim,
 const int              normalized,
 const double          *restrict d,
 const double          *restrict extents,
 const double           *restrict coords,
 double                *restrict max_dist2
 )
{

  int inbox = 0;
  *max_dist2 = 0.;

  PDM_UNUSED(normalized);

  //  if (normalized) {

  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim]) {
      double _max_dist2 = d[i] * (coords[i] - extents[i]);
      *max_dist2 += _max_dist2 * _max_dist2;
    }

    else if (coords[i] < extents[i]) {
      double _max_dist2 = d[i] * (coords[i] - extents[dim+i]);
      *max_dist2 += _max_dist2 * _max_dist2;
    }

    else {
      inbox += 1;
      double val1 = d[i] * (coords[i] - extents[i]);
      double val2 = d[i] * (coords[i] - extents[dim+i]);
      *max_dist2 += PDM_MAX (val1 * val1, val2 * val2);
    }
  }

  /* } */

  /* else { */

  /*   for (int i = 0; i < dim; i++) { */
  /*     if (coords[i] > extents[i+dim]) { */
  /*       double _max_dist2 = coords[i] - extents[i]; */
  /*       *max_dist2 += _max_dist2 * _max_dist2; */
  /*     } */

  /*     else if (coords[i] < extents[i]) { */
  /*       double _max_dist2 = coords[i] - extents[dim+i]; */
  /*       *max_dist2 += _max_dist2 * _max_dist2; */
  /*     } */

  /*     else { */
  /*       inbox += 1; */
  /*       double val1 = coords[i] - extents[i]; */
  /*       double val2 = coords[i] - extents[dim+i]; */
  /*       *max_dist2 += PDM_MAX (val1 * val1, val2 * val2); */
  /*     } */
  /*   } */
  /* } */

  return inbox == dim;

}


/**
 *
 * \brief Compute distance to a box
 *
 * \param [in]   dim        Dimension
 * \param [in]   extents    Box extents
 * \param [in]   coords     Point coords
 * \param [out]  min_dist2  Square of minimum distance
 * \param [out]  max_dist2  Sqaure of maximum distance
 *
 * \return 1 if point is in the box, 0 otherwise
 *
 */

inline static int
_box_dist2_min
(
 const int              dim,
 const int              normalized,
 const double          *restrict d,
 const double          *restrict extents,
 const double          *restrict coords,
 double                *restrict min_dist2
 )
{

  int inbox = 0;
  *min_dist2 = 0.;

  PDM_UNUSED(normalized);

  /* if (normalized) { */

    for (int i = 0; i < dim; i++) {
      if (coords[i] > extents[i+dim]) {
        double _min_dist2 = d[i] * (coords[i] - extents[dim+i]);
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else if (coords[i] < extents[i]) {
        double _min_dist2 = d[i] * (coords[i] - extents[i]);
        *min_dist2 += _min_dist2 * _min_dist2;
      }

      else {
        inbox += 1;
      }
    }

  /* } */

  /* else { */

  /*   for (int i = 0; i < dim; i++) { */
  /*     if (coords[i] > extents[i+dim]) { */
  /*       double _min_dist2 = coords[i] - extents[dim+i]; */
  /*       *min_dist2 += _min_dist2 * _min_dist2; */
  /*     } */

  /*     else if (coords[i] < extents[i]) { */
  /*       double _min_dist2 = coords[i] - extents[i]; */
  /*       *min_dist2 += _min_dist2 * _min_dist2; */
  /*     } */

  /*     else { */
  /*       inbox += 1; */
  /*     } */
  /*   } */
  /* } */

  return inbox == dim;

}





inline static int
_point_inside_box
(
 const int              dim,
 const double *restrict extents,
 const double *restrict coords
 )
{
  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim] || coords[i] < extents[i]) {
      return 0;
    }
  }

  return 1;
}


inline static int
_point_inside_ellipsoid
(
 const int              dim,
 const double *restrict extents,
 const double *restrict coords
 )
{
  double d2 = 0.;
  for (int i = 0; i < dim; i++) {
    double d = (2*coords[i] - extents[i+dim] - extents[i]) / (extents[i+dim] - extents[i]);
    d2 += d * d;
  }

  return d2 <= 1.;
}



inline static int
_intersect_line_box
(
 const int              dim,
 const double *restrict box_extents,
 const double *restrict line_origin,
 const double *restrict line_invdir
 )
{
  double tmin = 0.;
  double tmax = 1.;

  for (int i = 0; i < dim; i++) {
    double t1 = (box_extents[i]       - line_origin[i]) * line_invdir[i];
    double t2 = (box_extents[i + dim] - line_origin[i]) * line_invdir[i];

    if (line_invdir[i] < 0) {
      double tmp = t1;
      t1 = t2;
      t2 = tmp;
    }

    if (tmin > t2 || tmax < t1) {
      return 0;
    } else {
      tmin = PDM_MAX (tmin, t1);
      tmax = PDM_MIN (tmax, t2);
    }
  }

  return 1;
}



inline static int
_intersect_box_box
(
 const int              dim,
 const double *restrict box_extents_a,
 const double *restrict box_extents_b
 )
{
  // printf("box_extents_a = %12.5e/%12.5e/%12.5e  %12.5e/%12.5e/%12.5e \n",
  //        box_extents_a[0],
  //        box_extents_a[1],
  //        box_extents_a[2],
  //        box_extents_a[3],
  //        box_extents_a[4],
  //        box_extents_a[5]);
  // printf("box_extents_b = %12.5e/%12.5e/%12.5e  %12.5e/%12.5e/%12.5e \n",
  //        box_extents_b[0],
  //        box_extents_b[1],
  //        box_extents_b[2],
  //        box_extents_b[3],
  //        box_extents_b[4],
  //        box_extents_b[5]);

  // int intersect = 1;
  for (int i = 0; i < dim; i++) {
    if (box_extents_a[i] > box_extents_b[i+dim] || box_extents_b[i] > box_extents_a[i+dim]) {
      return 0;
    }
  }

  // printf("intersect --> %i \n", intersect);

  return 1;
}




/**
 *
 * \brief Add children nodes into stack
 *
 * \param [in]    bt            Box tree
 * \param [in]    dim           Dimension
 * \param [in]    id_curr_node  Identifier of current node
 * \param [in]    upper_bound   Upper_bound criteria to store a child in stack
 * \param [in]    pt            Distance to this point must be lesser van upper_bound
 * \param [inout] pos_stack     Position in the stack
 * \param [inout] stack         Stack
 *
 */

inline static void
_push_child_in_stack_v0
(
 PDM_box_tree_t *bt,
 const int       dim,
 const int        normalized,
 const double    *restrict d,
 const int       id_curr_node,
 const double    upper_bound,
 const double    *restrict pt,
 int             *restrict pos_stack,
 int             *restrict stack,
 int             *restrict inbox_stack,
 double          *restrict min_dist2_stack,
 int             flag,
 int             sorted
 )
{
  PDM_UNUSED(flag);

  PDM_box_tree_data_t *_local_data = bt->local_data;
  int sort_child[bt->n_children];
  double dist_child[bt->n_children];
  int inbox_child[bt->n_children];

  for (int i = 0; i < bt->n_children; i++) {
    dist_child[i] = HUGE_VAL;
  }

  /* Sort children and store them into the stack */

  const int *_child_ids = _local_data->child_ids + id_curr_node*bt->n_children;

  int _n_push = 0;

  double child_extents2[2*dim];

  for (int j = 0; j < bt->n_children; j++) {

    double child_min_dist2;

    int child_id = _child_ids[j];

    //const double *child_extents = bt->extents + dim * 2 * child_id;

    _node_t *curr_node = &(_local_data->nodes[child_id]);

    if (curr_node->n_boxes == 0) {
      continue;
    }

    _extents (dim, curr_node->morton_code, child_extents2);

    int inbox = _box_dist2_min (dim,
                                normalized,
                                d,
                                child_extents2,
                                pt,
                                &child_min_dist2);

    if (sorted) {
      int i1 = 0;
      for (i1 = _n_push; (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
        dist_child[i1] = dist_child[i1-1];
        sort_child[i1] = sort_child[i1-1];
        inbox_child[i1] = inbox_child[i1-1];
      }

      sort_child[i1] = _child_ids[j];
      dist_child[i1] =  child_min_dist2;
      inbox_child[i1] =  inbox;

      _n_push += 1;
    }
    else {
      if (((child_min_dist2 < upper_bound) || (inbox == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox;
        min_dist2_stack[*pos_stack] = child_min_dist2;

        (*pos_stack)++;
      }
    }

  }

  if (sorted) {
    for (int j =  _n_push - 1; j >= 0; j--) {
      int child_id = sort_child[j];

      if (((dist_child[j] < upper_bound) || (inbox_child[j] == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox_child[j];
        min_dist2_stack[*pos_stack] = dist_child[j];

        (*pos_stack)++;
      }
    }
  }

}

/**
 *
 * \brief Add children nodes into stack (extended version: local/rank tree data)
 *
 * \param [in]    bt            Box tree
 * \param [in]    i_rank
 * \param [in]    dim           Dimension
 * \param [in]    id_curr_node  Identifier of current node
 * \param [in]    upper_bound   Upper_bound criteria to store a child in stack
 * \param [in]    pt            Distance to this point must be lesser van upper_bound
 * \param [inout] pos_stack     Position in the stack
 * \param [inout] stack         Stack
 *
 */

inline static void
_push_child_in_stack_v2
(
 PDM_box_tree_t       *bt,
 PDM_box_tree_data_t  *box_tree_data,
 const int             dim,
 const int             normalized,
 const double         *restrict d,
 const int             id_curr_node,
 const double          upper_bound,
 const double          *restrict pt,
 int                   *restrict pos_stack,
 int                   *restrict stack,
 int                   *restrict inbox_stack,
 double                *restrict min_dist2_stack,
 int                   flag,
 int                   sorted
 )
{
  PDM_UNUSED(flag);

  PDM_box_tree_data_t *_tree_data = box_tree_data;
  int sort_child[bt->n_children];
  double dist_child[bt->n_children];
  int inbox_child[bt->n_children];

  for (int i = 0; i < bt->n_children; i++) {
    dist_child[i] = HUGE_VAL;
  }

  /* Sort children and store them into the stack */

  const int *_child_ids = _tree_data->child_ids + id_curr_node*bt->n_children;

  int _n_push = 0;

  double child_extents2[2*dim];

  for (int j = 0; j < bt->n_children; j++) {

    double child_min_dist2;

    int child_id = _child_ids[j];

    _node_t *curr_node = &(_tree_data->nodes[child_id]);

    if (curr_node->n_boxes == 0) {
      continue;
    }

    _extents (dim, curr_node->morton_code, child_extents2);

    int inbox = _box_dist2_min (dim,
                                normalized,
                                d,
                                child_extents2,
                                pt,
                                &child_min_dist2);

    if (sorted) {
      int i1 = 0;
      for (i1 = _n_push; (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
        dist_child[i1] = dist_child[i1-1];
        sort_child[i1] = sort_child[i1-1];
        inbox_child[i1] = inbox_child[i1-1];
      }

      sort_child[i1] = _child_ids[j];
      dist_child[i1] =  child_min_dist2;
      inbox_child[i1] =  inbox;

      _n_push += 1;
    }
    else {
      if (((child_min_dist2 < upper_bound) || (inbox == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox;
        min_dist2_stack[*pos_stack] = child_min_dist2;

	(*pos_stack)++;
      }
    }

  }

  if (sorted) {
    for (int j =  _n_push - 1; j >= 0; j--) {
      int child_id = sort_child[j];

      if (((dist_child[j] < upper_bound) || (inbox_child[j] == 1))) {

        stack[*pos_stack]           = child_id; /* push child in th stack */
        inbox_stack[*pos_stack]     = inbox_child[j];
        min_dist2_stack[*pos_stack] = dist_child[j];

        (*pos_stack)++;
      }
    }
  }
}




/*----------------------------------------------------------------------------
 * Get minimum coordinates for a given box.
 *
 * parameters:
 *   box_set <-- pointer to box set structure
 *   box_id  <-- id of box
 *
 * returns:
 *   pointer to minimum box coordinates
 *---------------------------------------------------------------------------*/

inline static const double *
_box_min(const PDM_box_set_t   *boxes,
         int              box_id)
{
  return boxes->local_boxes->extents + box_id*boxes->dim*2;
}

/*----------------------------------------------------------------------------
 * Get maxmum coordinates for a given box.
 *
 * parameters:
 *   box_set <-- pointer to box set structure
 *   box_id  <-- id of box
 *
 * returns:
 *   pointer to maximum box coordinates
 *---------------------------------------------------------------------------*/

inline static const double *
_box_max(const PDM_box_set_t   *boxes,
         int              box_id)
{
  return boxes->local_boxes->extents + box_id*boxes->dim*2 + boxes->dim;
}


/*----------------------------------------------------------------------------
 * Test for intersection between two bounding boxes.
 *
 * parameters:
 *   extents <-- array of box extents
 *   id_0    <-- id of first box in extents
 *   extentsB<-- box extents to intersect
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

inline static _Bool
_boxes_intersect_3d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*6;

  if (   e0[0] > extentsB[3] || extentsB[0] > e0[3]
	 || e0[1] > extentsB[4] || extentsB[1] > e0[4]
	 || e0[2] > extentsB[5] || extentsB[2] > e0[5])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_2d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*4;

  if (   e0[0] > extentsB[2] || extentsB[0] > e0[2]
	 || e0[1] > extentsB[3] || extentsB[1] > e0[3])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_1d(const double  *extents,
                    int            id_0,
                    const double  *extentsB)
{
  const double *e0 = extents + id_0*2;

  if (   e0[0] > extentsB[1] || extentsB[0] > e0[1])
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------------
 * Test for intersection between two bounding boxes.
 *
 * parameters:
 *   extents <-- array of box extents
 *   id_0    <-- id of first box
 *   id_1    <-- id of second box
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

inline static _Bool
_boxes_intern_intersect_3d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*6;
  const double *e1 = extents + id_1*6;

  if (   e0[0] > e1[3] || e1[0] > e0[3]
	 || e0[1] > e1[4] || e1[1] > e0[4]
	 || e0[2] > e1[5] || e1[2] > e0[5])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intern_intersect_2d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*4;
  const double *e1 = extents + id_1*4;

  if (   e0[0] > e1[2] || e1[0] > e0[2]
	 || e0[1] > e1[3] || e1[1] > e0[3])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intern_intersect_1d(const double  *extents,
                           int          id_0,
                           int          id_1)
{
  const double *e0 = extents + id_0*2;
  const double *e1 = extents + id_1*2;

  if (   e0[0] > e1[1] || e1[0] > e0[1])
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------
 * Update octree stat structure (min, max, mean, box ratio, ...)
 *
 * parameters:
 *   bt      <-> pointer on the PDM_box_tree_t structure to deal with
 *   node_id <-- node on which we collect data
 *----------------------------------------------------------------------------*/

static void
_update_tree_stats(PDM_box_tree_t  *bt,
                   int        node_id)
{
  int  i;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  *node = _local_data->nodes + node_id;

  if (node->is_leaf == false) {
    int n_children = bt->n_children;
    const int *_child_ids = _local_data->child_ids + node_id*bt->n_children;
    for (i = 0; i < n_children; i++)
      _update_tree_stats(bt, _child_ids[i]);
  }

  else { /* leaf node */

    PDM_box_tree_stats_t  s = bt->stats;

    s.n_leaves += 1;
    s.n_linked_boxes += node->n_boxes;

    if (node->n_boxes > bt->threshold)
      s.n_spill_leaves += 1;

    s.min_linked_boxes = PDM_MIN(s.min_linked_boxes, node->n_boxes);
    s.max_linked_boxes = PDM_MAX(s.max_linked_boxes, node->n_boxes);
    s.max_level_reached = PDM_MAX(s.max_level_reached, node->morton_code.L);

    bt->stats = s;
  }
}

/*----------------------------------------------------------------------------
 * Define box_tree->stat structure (min, max, mean, box ratio, ...)
 *
 * parameters:
 *   bt <-> pointer to the box-tree structure
 *----------------------------------------------------------------------------*/

static void
_get_box_tree_stats(PDM_box_tree_t  *bt)
{
  if (bt == NULL)
    return;

  /* Initialize statistics */

  bt->stats.max_level_reached = 0;

  bt->stats.n_leaves = 0;
  bt->stats.n_spill_leaves = 0;
  bt->stats.n_linked_boxes = 0;

  bt->stats.min_linked_boxes = INT_MAX;
  bt->stats.max_linked_boxes = 0;

  /* Recursively update stats, starting from root */
  PDM_box_tree_data_t *_local_data = bt->local_data;
  if ( _local_data == NULL )
    return;

  if (_local_data->nodes != NULL)
    _update_tree_stats(bt, 0);
}

/*----------------------------------------------------------------------------
 * Get the coordinates in the grid for the current point at this level.
 *
 * parameters:
 *   level      <--   level on which we want the coordinates
 *   coords     <--   coords of the point to translate in the octree grid.
 *   XYZ        <->   pointer to the X, Y, Z coordinates in the grid
 *----------------------------------------------------------------------------*/

inline static void
_get_grid_coords_3d(PDM_morton_int_t  level,
                    const double      coords[3],
                    double            XYZ[])
{
  PDM_morton_int_t  refinement = 1 << level;

  XYZ[0] = coords[0] * refinement;
  XYZ[1] = coords[1] * refinement;
  XYZ[2] = coords[2] * refinement;
}

inline static void
_get_grid_coords_2d(PDM_morton_int_t  level,
                    const double      coords[2],
                    double            XY[])
{
  PDM_morton_int_t  refinement = 1 << level;

  XY[0] = coords[0] * refinement;
  XY[1] = coords[1] * refinement;
}

inline static void
_get_grid_coords_1d(PDM_morton_int_t  level,
                    const double      coords[1],
                    double            x_tmp[])
{
  PDM_morton_int_t  refinement = 1 << level;

  x_tmp[0] = coords[0] * refinement;
}

/*----------------------------------------------------------------------------
 * Return true if a leaf intersects a box, false otherwise.
 *
 * parameters:
 *   morton_code <-- Morton code of the leaf
 *   min_box     <-- coordinates of min. point of the bounding box
 *   max_box     <-- coordinates of max. point of the bounding box
 *
 * returns:
 *   true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_node_intersect_box_3d(PDM_morton_code_t  morton_code,
                       const double   min_box[3],
                       const double   max_box[3])
{
  int  i;
  double  min_oct[3], max_oct[3];

  for (i = 0; i < 3; i++) {
    min_oct[i] = (double)morton_code.X[i];
    max_oct[i] = (double)(morton_code.X[i] + 1);
  }

  /* printf ("min_oct : %12.5e %12.5e %12.5e\n", min_oct[0], min_oct[1], min_oct[2]); */
  /* printf ("max_oct : %12.5e %12.5e %12.5e\n", max_oct[0], max_oct[1], max_oct[2]); */

  /* printf ("min_box : %12.5e %12.5e %12.5e\n", min_box[0], min_box[1], min_box[2]); */
  /* printf ("max_box : %12.5e %12.5e %12.5e\n\n", max_box[0], max_box[1], max_box[2]); */


  if (   min_box[0] > max_oct[0] || min_oct[0] > max_box[0]
	 || min_box[1] > max_oct[1] || min_oct[1] > max_box[1]
         || min_box[2] > max_oct[2] || min_oct[2] > max_box[2]){
    return false;
  }
  else{
    return true;
  }
}

inline static _Bool
_node_intersect_box_2d(PDM_morton_code_t  morton_code,
                       const double   min_box[2],
                       const double   max_box[2])
{
  int  i;
  double  min_oct[2], max_oct[2];

  for (i = 0; i < 2; i++) {
    min_oct[i] = (double)morton_code.X[i];
    max_oct[i] = (double)(morton_code.X[i] + 1);
  }

  if (   min_box[0] > max_oct[0] || min_oct[0] > max_box[0]
	 || min_box[1] > max_oct[1] || min_oct[1] > max_box[1])
    return false;
  else
    return true;
}

inline static _Bool
_node_intersect_box_1d(PDM_morton_code_t  morton_code,
                       const double   min_box[1],
                       const double   max_box[1])
{
  double  min_oct, max_oct;

  min_oct = (double)morton_code.X[0];
  max_oct = (double)(morton_code.X[0] + 1);

  if (min_box[0] > max_oct || min_oct > max_box[0])
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------------
 * Split a node into its children and evaluate the box distribution.
 *
 * parameters:
 *   bt       <-> pointer to the box tree being built
 *   boxes    <-- pointer to the associated box set structure
 *   node_id  <-- id of the node to split
 *
 * returns:
 *   the number of associations between nodes and their children
 *----------------------------------------------------------------------------*/

static int
_evaluate_splitting_cartesian(PDM_box_tree_t       *bt,
                             const PDM_box_set_t  *boxes,
                                   int             node_id)
{
  PDM_morton_code_t  children[bt->n_children];
  double             child_extents[2*boxes->dim*bt->n_children];

  int  n_linked_boxes = 0;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  node = _local_data->nodes[node_id];

  /* Define a Morton code for each child */
  PDM_morton_get_children(boxes->dim, node.morton_code, children);

  /* Compute extents of each child */
  for(int i = 0; i < bt->n_children; ++i) {
    _extents(boxes->dim, children[i], child_extents + 6*i);
  }

  /* Loop on boxes associated to the node_id */
  for (int j = 0; j < node.n_boxes; j++) {
    int   box_id = _local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    /* Brute force */
    for (int i = 0; i < bt->n_children; i++) {
      if(_intersect_box_box(boxes->dim, box_min, child_extents + 2*boxes->dim*i)){
        n_linked_boxes += 1;
      }
    }
  } /* End of loop on boxes */
  return n_linked_boxes;
}


/*----------------------------------------------------------------------------
 * Split a node into its children and evaluate the box distribution.
 *
 * parameters:
 *   bt       <-> pointer to the box tree being built
 *   boxes    <-- pointer to the associated box set structure
 *   node_id  <-- id of the node to split
 *
 * returns:
 *   the number of associations between nodes and their children
 *----------------------------------------------------------------------------*/

static int
_evaluate_splitting_3d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[8];

  int  n_linked_boxes = 0;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  node = _local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 3);

  /* Define a Morton code for each child */

  PDM_morton_get_children(3, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = _local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);


    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   PDM_morton_compare(3, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_2d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[4];

  int  n_linked_boxes = 0;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  node = _local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 2);

  /* Define a Morton code for each child */

  PDM_morton_get_children(2, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = _local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   PDM_morton_compare(2, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_1d(PDM_box_tree_t       *bt,
                       const PDM_box_set_t  *boxes,
                       int             node_id)
{
  int   i, j;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[2];

  int  n_linked_boxes = 0;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  node = _local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  /* Define a Morton code for each child */

  PDM_morton_get_children(1, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[1], max_grid_coord[1];

    int   box_id = _local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(1, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   PDM_morton_compare(1, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}


/*----------------------------------------------------------------------------
 * Evaluate the intersection between a box and a node tree
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxesB          <--  pointer boxes that intersect associated tree boxes
 *   boxB_id         <--  id of the current boxB in boxesB
 *   node_id         <--  id of the starting node
 *----------------------------------------------------------------------------*/

static _Bool
_boxes_intersect_node_3d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 3);

  double  min_grid_coord[3], max_grid_coord[3];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(3, level, box_min);
  max_code = PDM_morton_encode(3, level, box_max);

  if (   PDM_morton_compare(3, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_3d(level, box_min, min_grid_coord);
    _get_grid_coords_3d(level, box_max, max_grid_coord);

    return _node_intersect_box_3d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(3, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(3, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


static _Bool
_boxes_intersect_node_2d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 2);

  double  min_grid_coord[2], max_grid_coord[2];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(2, level, box_min);
  max_code = PDM_morton_encode(2, level, box_max);

  if (   PDM_morton_compare(2, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_2d(level, box_min, min_grid_coord);
    _get_grid_coords_2d(level, box_max, max_grid_coord);

    return _node_intersect_box_2d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(2, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(2, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


static _Bool
_boxes_intersect_node_1d (const PDM_box_tree_t     *bt,
                          const PDM_box_set_t      *boxesB,
                          int                       boxB_id,
                          int                       node_id)
{
  PDM_morton_code_t  min_code, max_code;

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  level = node.morton_code.L;

  assert(boxesB->dim == 1);

  double  min_grid_coord[1], max_grid_coord[1];

  const double  *box_min = _box_min(boxesB, boxB_id);
  const double  *box_max = _box_max(boxesB, boxB_id);

  min_code = PDM_morton_encode(1, level, box_min);
  max_code = PDM_morton_encode(1, level, box_max);

  if (   PDM_morton_compare(1, min_code, max_code)
         == PDM_MORTON_DIFFERENT_ID) {

    _get_grid_coords_1d(level, box_min, min_grid_coord);
    _get_grid_coords_1d(level, box_max, max_grid_coord);

    return _node_intersect_box_1d (node.morton_code, min_grid_coord, max_grid_coord);

  }

  else { /* Box is included in the same octant */

    assert(   PDM_morton_compare(1, max_code, min_code)
              == PDM_MORTON_EQUAL_ID);

    return (PDM_morton_compare(1, min_code, node.morton_code) == PDM_MORTON_EQUAL_ID);

  }
}


/*----------------------------------------------------------------------------
 * Evaluate the intersection between a box and a node tree
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxesB          <--  pointer boxes that intersect associated tree boxes
 *   boxB_id         <--  id of the current boxB in boxesB
 *   node_id         <--  id of the starting node
 *----------------------------------------------------------------------------*/

static _Bool
_boxes_intersect_node (const PDM_box_tree_t     *bt,
                       const PDM_box_set_t      *boxesB,
                       int                       boxB_id,
                       int                       node_id)
{
  _Bool is_intersect = false;
  if (boxesB->dim == 3) {
    is_intersect = _boxes_intersect_node_3d (bt, boxesB, boxB_id, node_id);
  }
  else if (boxesB->dim == 2) {
    is_intersect = _boxes_intersect_node_2d (bt, boxesB, boxB_id, node_id);
  }
  else if (boxesB->dim == 1) {
    is_intersect = _boxes_intersect_node_1d (bt, boxesB, boxB_id, node_id);
  }
  return is_intersect;
}



/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree to help
 * determine if we should add a level to the tree structure.
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxes           <--  pointer to the associated box set structure
 *   node_id         <--  id of the starting node
 *   build_type      <--  layout variant for building the tree structure
 *   next_level_size -->  size of box_ids for the next level
 *----------------------------------------------------------------------------*/

static void
_count_next_level(PDM_box_tree_t           *bt,
                  const PDM_box_set_t      *boxes,
                  int                       node_id,
                  PDM_box_tree_sync_t       build_type,
                  int                      *next_level_size)
{
  int   i;

  PDM_box_tree_data_t *_local_data = bt->local_data;
  _node_t  *node = _local_data->nodes + node_id;

  if (node->is_leaf == false) {

    assert(_local_data->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _count_next_level(bt,
                        boxes,
                        _local_data->child_ids[bt->n_children*node_id + i],
                        build_type,
                        next_level_size);

  }

  else { /* if (node->is_leaf == true) */
    if (   node->n_boxes < bt->threshold
	   && node_id != 0                    /* Root node is always divided */
	   && build_type == PDM_BOX_TREE_ASYNC_LEVEL)
      *next_level_size += node->n_boxes;

    else { /* Split node and evaluate box distribution between its children */
      if (boxes->dim == 3)
        *next_level_size += _evaluate_splitting_3d(bt, boxes, node_id);
      else if (boxes->dim == 2)
        *next_level_size += _evaluate_splitting_2d(bt, boxes, node_id);
      else if (boxes->dim == 1)
        *next_level_size += _evaluate_splitting_1d(bt, boxes, node_id);
    }
  }
}


/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree to help
 * determine if we should add a level to the tree structure.
 *
 * parameters:
 *   bt              <->  pointer to the box tree being built
 *   boxes           <--  pointer to the associated box set structure
 *   node_id         <--  id of the starting node
 *   build_type      <--  layout variant for building the tree structure
 *   next_level_size -->  size of box_ids for the next level
 *----------------------------------------------------------------------------*/

static void
_count_next_level_cartesian(PDM_box_tree_t           *bt,
                           const PDM_box_set_t      *boxes,
                           int                       node_id,
                           PDM_box_tree_sync_t       build_type,
                           int                      *next_level_size)
{
  PDM_box_tree_data_t *_local_data = bt->local_data;
  _node_t  *node = _local_data->nodes + node_id;

  if (node->is_leaf == false) {

    assert(_local_data->child_ids[bt->n_children*node_id] > 0);

    for (int i = 0; i < bt->n_children; i++)
      _count_next_level_cartesian(bt,
                                 boxes,
                                 _local_data->child_ids[bt->n_children*node_id + i],
                                 build_type,
                                 next_level_size);

  }

  else { /* if (node->is_leaf == true) */
    if (   node->n_boxes < bt->threshold
     && node_id != 0                    /* Root node is always divided */
     && build_type == PDM_BOX_TREE_ASYNC_LEVEL)
      *next_level_size += node->n_boxes;

    else { /* Split node and evaluate box distribution between its children */
      *next_level_size += _evaluate_splitting_cartesian(bt, boxes, node_id);
    }
  }
}


/*----------------------------------------------------------------------------
 * Test if we have to continue the building of the box tree.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   build_type <--  layout variant for building the tree structure
 *   next_size  -->  size of box_ids for the next tree if required
 *
 * returns:
 *   true if we should continue, false otherwise.
 *----------------------------------------------------------------------------*/

static _Bool
_recurse_tree_build(      PDM_box_tree_t      *bt,
                    const PDM_box_set_t       *boxes,
                          PDM_box_tree_sync_t  build_type,
                          int                  build_cartesian,
                          int                 *next_size)
{
  int  state = 0;
  int  _next_size = 0;

  _Bool retval = false;

  int  n_ranks = 1;
  PDM_MPI_Comm comm = boxes->comm;

  if (comm != PDM_MPI_COMM_NULL)
    PDM_MPI_Comm_size(comm, &n_ranks);

  PDM_box_tree_data_t *_local_data = bt->local_data; //*
  _local_data->n_build_loops += 1; //*

  if (bt == NULL) // if (_local_data == NULL)
    state = 1;

  /* To avoid infinite loop on tree building */

  if (_local_data->n_build_loops > PDM_BOX_TREE_MAX_BUILD_LOOPS) //*
    state = 1;

  /* A sufficient accuracy has been reached */

  if ((int)(bt->stats.max_level_reached) == bt->max_level)
    state = 1;

  /* Algorithm is converged. No need to go further */

  if (   bt->stats.n_spill_leaves == 0
	 && bt->stats.max_level_reached > 0)
    state = 1;

  if (n_ranks > 1 && build_type == PDM_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    PDM_MPI_Allreduce(&state, &global_state, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);
    state = global_state; /* Stop if all ranks require it */
  }

  if (state == 0) {

    float box_ratio;

    /* Limit, to avoid excessive memory usage */
    if(build_cartesian == 0) {
      _count_next_level(bt,
                        boxes,
                        0,  /* Starts from root */
                        build_type,
                        &_next_size);
    } else {
      _count_next_level_cartesian(bt,
                                  boxes,
                                  0,  /* Starts from root */
                                  build_type,
                                  &_next_size);
    }

    if (bt->stats.n_boxes > 0)
      box_ratio = (float) ((_next_size*1.0)/bt->stats.n_boxes);
    else
      box_ratio = 0;

    if (bt->stats.max_level_reached > 0 && box_ratio > bt->max_box_ratio)
      state = 1;

  }

  if (n_ranks > 1 && build_type == PDM_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    PDM_MPI_Allreduce(&state, &global_state, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
    state = global_state; /* Stop as as soon as any rank requires it */
  }

  /* If no condition is encoutered, we have to continue */

  *next_size = _next_size;

  if (state == 0)
    retval = true;

  return retval;
}


/*----------------------------------------------------------------------------
 * Create a box tree data by copying from another.
 *
 * parameters:
 *   dest <-> pointer to destination box tree data
 *   src  <-- pointer to source box tree data
 *----------------------------------------------------------------------------*/
static void
_copy_tree_data(PDM_box_tree_data_t        *dest,
                const PDM_box_tree_data_t  *src,
                const int                   n_children,
                const int                   n_linked_boxes)
{
  assert(dest != NULL && src != NULL);

  memcpy(dest, src, sizeof(PDM_box_tree_data_t));

  dest->nodes = (_node_t *) malloc(dest->n_max_nodes * sizeof(_node_t));

  dest->child_ids = (int *) malloc(dest->n_max_nodes*n_children * sizeof(int));

  dest->box_ids = (int *) malloc(n_linked_boxes * sizeof(int));

  memcpy(dest->nodes,
	 src->nodes,
	 dest->n_nodes * sizeof(_node_t));

  memcpy(dest->child_ids,
         src->child_ids,
         dest->n_nodes * n_children * sizeof(int));

  memcpy(dest->box_ids,
         src->box_ids,
         n_linked_boxes * sizeof(int));
}


/*----------------------------------------------------------------------------
 * Create a box tree by copying from another.
 *
 * parameters:
 *   dest <-> pointer to destination box tree
 *   src  <-- pointer to source box tree
 *----------------------------------------------------------------------------*/

static void
_copy_tree(PDM_box_tree_t        *dest,
           const PDM_box_tree_t  *src)
{
  assert(dest != NULL && src != NULL);

  memcpy(dest, src, sizeof(PDM_box_tree_t));

  dest->local_data = (PDM_box_tree_data_t *) malloc(sizeof(PDM_box_tree_data_t));


  /* -->> _copy_tree_data
     dest->nodes = (_node_t *) malloc(dest->n_max_nodes * sizeof(_node_t));
     dest->child_ids = (int *) malloc(dest->n_max_nodes*dest->n_children * sizeof(int));

     dest->box_ids = (int *) malloc((dest->stats).n_linked_boxes * sizeof(int));

     memcpy(dest->nodes, src->nodes, dest->n_nodes * sizeof(_node_t));
     memcpy(dest->child_ids,
     src->child_ids,
     dest->n_nodes * src->n_children * sizeof(int));

     memcpy(dest->box_ids,
     src->box_ids,
     (dest->stats).n_linked_boxes * sizeof(int));
     <<-- */

  _copy_tree_data (dest->local_data,
                   src->local_data,
                   src->n_children,
                   (dest->stats).n_linked_boxes);

}


/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_data_t structure.
 *
 * parameters:
 *   btd <-- pointer to pointer to PDM_box_tree_data_t structure to destroy
 *----------------------------------------------------------------------------*/

static void
_free_tree_data_arrays(PDM_box_tree_data_t  *btd)
{
  assert(btd != NULL);

  free(btd->nodes);
  free(btd->child_ids);
  free(btd->box_ids);
}




/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to PDM_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

static void
_free_tree_arrays(PDM_box_tree_t  *bt)
{
  assert(bt != NULL);

  /*PDM_box_tree_data_t *_local_data = bt->local_data;
    assert(_local_data != NULL);
    free(_local_data->nodes);
    free(_local_data->child_ids);
    free(_local_data->box_ids);*/
  _free_tree_data_arrays(bt->local_data);
}

/*----------------------------------------------------------------------------
 * Create a new node from a Morton code and associate it to a tree.
 *
 * parameters:
 *   bt          <->  pointer to the box tree being built
 *   morton_code <--  Morton identification number related to the node
 *   node_id     <--  id of the starting node
 *----------------------------------------------------------------------------*/

static inline void
_new_node(PDM_box_tree_t     *bt,
          PDM_morton_code_t   morton_code,
          int                 node_id)
{
  int  i;
  _node_t *node;

  assert(bt != NULL);

  PDM_box_tree_data_t *_local_data = bt->local_data;
  assert(_local_data != NULL);

  node = _local_data->nodes + node_id;

  if ((int)(morton_code.L) > bt->max_level) {
    PDM_error(__FILE__, __LINE__, 0,
	      "Error adding a new node in box tree (%p).\n"
	      "Max level reached. Current level: %u and Max level: %d\n",
	      (void *)bt, morton_code.L, bt->max_level);
    abort();
  }

  node->is_leaf = true;
  node->morton_code = morton_code;

  node->n_boxes = 0;
  node->start_id = -1; /* invalid value by default */
  node->extra_weight = 0;

  for (i = 0; i < bt->n_children; i++) {
    _local_data->child_ids[node_id*bt->n_children + i] = -1;
  }

}

/*----------------------------------------------------------------------------
 * Split a node into its children and define the new box distribution.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_split_node_3d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[8];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->local_data->n_nodes;
  _node_t  split_node = next_bt->local_data->nodes[node_id];

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 8);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 8 > next_bt->local_data->n_max_nodes) {
    assert(next_bt->local_data->n_max_nodes > 0);
    next_bt->local_data->n_max_nodes += PDM_MAX (9, next_bt->local_data->n_max_nodes/3);
    next_bt->local_data->nodes = (_node_t *) realloc((void *) next_bt->local_data->nodes, next_bt->local_data->n_max_nodes * sizeof(_node_t));
    next_bt->local_data->child_ids = (int *) realloc((void *) next_bt->local_data->child_ids, next_bt->local_data->n_max_nodes*8 * sizeof(int));

  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(3, node.morton_code, children);

  for (i = 0; i < 8; i++) {

    const int   new_id = n_init_nodes + i;
    next_bt->local_data->child_ids[node_id*8 + i] = new_id;

    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = node.n_boxes;
  split_node.is_leaf = false;
  split_node.extra_weight = 0;

  next_bt->local_data->nodes[node_id] = split_node;
  next_bt->local_data->n_nodes = n_init_nodes + 8;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {
    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord)) {
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;

        }
      }

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   PDM_morton_compare(3, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 8; i++) {

    (next_bt->local_data->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 8; i++)
    (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[3], max_grid_coord[3];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(3, next_level, box_min);
    max_code = PDM_morton_encode(3, next_level, box_max);

    if (   PDM_morton_compare(3, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {

        if (_node_intersect_box_3d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
	    + (next_bt->local_data->nodes[sub_id]).start_id;
          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same octant */

      assert(   PDM_morton_compare(3, max_code, min_code)
                == PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {

        if (   PDM_morton_compare(3, min_code, children[i])
               == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
            + (next_bt->local_data->nodes[sub_id]).start_id;

          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 8; i++) {
    _node_t  n1 = next_bt->local_data->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->local_data->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
	    == (  next_bt->local_data->nodes[n_init_nodes + 8 - 1].start_id
		  + next_bt->local_data->nodes[n_init_nodes + 8 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}


/*----------------------------------------------------------------------------
 * Split a node into its children and define the new box distribution.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_split_node_cartesian_3d(PDM_box_tree_t       *bt,
                        PDM_box_tree_t       *next_bt,
                        const PDM_box_set_t  *boxes,
                        int                   node_id,
                        int                  *shift_ids)
{
  int j, i;
  PDM_morton_code_t  children[8];
  double             child_extents[6*8];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->local_data->n_nodes;
  _node_t  split_node = next_bt->local_data->nodes[node_id];

  const _node_t  node = bt->local_data->nodes[node_id];

  assert(bt->n_children == 8);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 8 > next_bt->local_data->n_max_nodes) {
    assert(next_bt->local_data->n_max_nodes > 0);
    next_bt->local_data->n_max_nodes += PDM_MAX (9, next_bt->local_data->n_max_nodes/3);
    next_bt->local_data->nodes = (_node_t *) realloc((void *) next_bt->local_data->nodes, next_bt->local_data->n_max_nodes * sizeof(_node_t));
    next_bt->local_data->child_ids = (int *) realloc((void *) next_bt->local_data->child_ids, next_bt->local_data->n_max_nodes*8 * sizeof(int));

  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(3, node.morton_code, children);
  /* Compute extents of each child */
  for(i = 0; i < 8; ++i) {
    _extents(3, children[i], child_extents + 6*i);
  }

  for (i = 0; i < 8; i++) {

    const int   new_id = n_init_nodes + i;
    next_bt->local_data->child_ids[node_id*8 + i] = new_id;

    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = node.n_boxes;
  split_node.is_leaf = false;
  split_node.extra_weight = 0;

  next_bt->local_data->nodes[node_id] = split_node;
  next_bt->local_data->n_nodes = n_init_nodes + 8;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);

    /* Brute force */
    for (i = 0; i < 8; i++) {
      if(_intersect_box_box(3, box_min, child_extents + 6*i)){
        (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
      }
    }
  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 8; i++) {
    (next_bt->local_data->nodes[n_init_nodes + i]).start_id = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes;
  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 8; i++) {
    (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes = 0;
  }

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);

    /* Brute force */
    for (i = 0; i < 8; i++) {
      if(_intersect_box_box(3, box_min, child_extents + 6*i)){
        const int sub_id = n_init_nodes + i;
        const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes + (next_bt->local_data->nodes[sub_id]).start_id;
        next_bt->local_data->box_ids[shift] = box_id;
        (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
      }
    }

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 8; i++) {
    _node_t  n1 = next_bt->local_data->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->local_data->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
      == (  next_bt->local_data->nodes[n_init_nodes + 8 - 1].start_id
      + next_bt->local_data->nodes[n_init_nodes + 8 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}






static void
_split_node_2d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[4];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->local_data->n_nodes;
  _node_t  split_node = next_bt->local_data->nodes[node_id];

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 4);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 4 > next_bt->local_data->n_max_nodes) {
    assert(next_bt->local_data->n_max_nodes > 0);
    next_bt->local_data->n_max_nodes += PDM_MAX(5, next_bt->local_data->n_max_nodes/3);
    next_bt->local_data->nodes = (_node_t *) realloc((void *) next_bt->local_data->nodes, next_bt->local_data->n_max_nodes * sizeof(_node_t));
    next_bt->local_data->child_ids = (int *) realloc((void *) next_bt->local_data->child_ids, next_bt->local_data->n_max_nodes*4 * sizeof(int));
  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(2, node.morton_code, children);

  for (i = 0; i < 4; i++) {
    const int   new_id = n_init_nodes + i;
    next_bt->local_data->child_ids[node_id*4 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;
  split_node.extra_weight = 0;

  next_bt->local_data->nodes[node_id] = split_node;
  next_bt->local_data->n_nodes = n_init_nodes + 4;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   PDM_morton_compare(2, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 4; i++) {

    (next_bt->local_data->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 4; i++)
    (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(2, next_level, box_min);
    max_code = PDM_morton_encode(2, next_level, box_max);

    if (   PDM_morton_compare(2, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {

        if (_node_intersect_box_2d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
	    + (next_bt->local_data->nodes[sub_id]).start_id;
          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same quadrant */

      assert(   PDM_morton_compare(2, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {

        if (   PDM_morton_compare(2, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
	    + (next_bt->local_data->nodes[sub_id]).start_id;

          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 4; i++) {
    _node_t  n1 = next_bt->local_data->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->local_data->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
	    == (  next_bt->local_data->nodes[n_init_nodes + 4 - 1].start_id
		  + next_bt->local_data->nodes[n_init_nodes + 4 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

static void
_split_node_1d(PDM_box_tree_t       *bt,
               PDM_box_tree_t       *next_bt,
               const PDM_box_set_t  *boxes,
               int             node_id,
               int            *shift_ids)
{
  int j, i;
  PDM_morton_code_t  min_code, max_code;
  PDM_morton_code_t  children[2];

  int   n_linked_boxes = 0;
  int   _shift_ids = *shift_ids;
  int   n_init_nodes = next_bt->local_data->n_nodes;
  _node_t  split_node = next_bt->local_data->nodes[node_id];

  const _node_t  node = bt->local_data->nodes[node_id];
  const PDM_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 2);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 2 > next_bt->local_data->n_max_nodes) {
    assert(next_bt->local_data->n_max_nodes > 0);
    next_bt->local_data->n_max_nodes += PDM_MAX (3, next_bt->local_data->n_max_nodes/3);
    next_bt->local_data->nodes = (_node_t *) realloc((void *) next_bt->local_data->nodes, next_bt->local_data->n_max_nodes * sizeof(_node_t));
    next_bt->local_data->child_ids = (int *) realloc((void *) next_bt->local_data->child_ids, next_bt->local_data->n_max_nodes*2 * sizeof(int));
  }

  /* Define a Morton code for each child and create the children nodes */

  PDM_morton_get_children(1, node.morton_code, children);

  for (i = 0; i < 2; i++) {

    const int   new_id = n_init_nodes + i;
    next_bt->local_data->child_ids[node_id*2 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;
  split_node.extra_weight = 0;

  next_bt->local_data->nodes[node_id] = split_node;
  next_bt->local_data->n_nodes = n_init_nodes + 2;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same segment */

      assert(   PDM_morton_compare(1, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   PDM_morton_compare(1, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {
          (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 2; i++) {

    (next_bt->local_data->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 2; i++)
    (next_bt->local_data->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    double  min_grid_coord[2], max_grid_coord[2];

    int   box_id = bt->local_data->box_ids[node.start_id + j];
    const double  *box_min = _box_min(boxes, box_id);
    const double  *box_max = _box_max(boxes, box_id);

    min_code = PDM_morton_encode(1, next_level, box_min);
    max_code = PDM_morton_encode(1, next_level, box_max);

    if (   PDM_morton_compare(1, min_code, max_code)
	   == PDM_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {

        if (_node_intersect_box_1d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
	    + (next_bt->local_data->nodes[sub_id]).start_id;
          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same segment */

      assert(   PDM_morton_compare(1, max_code, min_code)
		== PDM_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {

        if (   PDM_morton_compare(1, min_code, children[i])
	       == PDM_MORTON_EQUAL_ID) {

          const int sub_id = n_init_nodes + i;
          const int shift =   (next_bt->local_data->nodes[sub_id]).n_boxes
	    + (next_bt->local_data->nodes[sub_id]).start_id;

          next_bt->local_data->box_ids[shift] = box_id;
          (next_bt->local_data->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 2; i++) {
    _node_t  n1 = next_bt->local_data->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->local_data->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
	    == (  next_bt->local_data->nodes[n_init_nodes + 2 - 1].start_id
		  + next_bt->local_data->nodes[n_init_nodes + 2 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree when adding
 * a level to the tree structure.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   build_type <--  layout variant for building the tree structure
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_build_next_level(PDM_box_tree_t       *bt,
                  PDM_box_tree_t       *next_bt,
                  const PDM_box_set_t  *boxes,
                  int             node_id,
                  PDM_box_tree_sync_t   build_type,
                  int            *shift_ids)
{
  int   i;

  int   _shift_ids = *shift_ids;
  const _node_t  *cur_node = bt->local_data->nodes + node_id;

  if (cur_node->is_leaf == false) {

    assert(bt->local_data->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _build_next_level(bt,
                        next_bt,
                        boxes,
                        bt->local_data->child_ids[bt->n_children*node_id + i],
                        build_type,
                        &_shift_ids);
  }

  else { /* if (node->is_leaf == true) */

    if (   cur_node->n_boxes < bt->threshold
           && node_id != 0                    /* Root node is always divided */
           && build_type == PDM_BOX_TREE_ASYNC_LEVEL) {

      /* Copy related box_ids in the new next_ids */

      _node_t *next_node = next_bt->local_data->nodes + node_id;

      next_node->n_boxes = cur_node->n_boxes;
      next_node->start_id = _shift_ids;
      next_node->extra_weight = 0;

      for (i = 0; i < cur_node->n_boxes; i++)
        next_bt->local_data->box_ids[_shift_ids++]
          = bt->local_data->box_ids[cur_node->start_id + i];
    }
    else {  /* Split node and evaluate box distribution between its children */

      if (boxes->dim == 3)
        _split_node_3d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
      else if (boxes->dim == 2)
        _split_node_2d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
      else if (boxes->dim == 1)
        _split_node_1d(bt,
                       next_bt,
                       boxes,
                       node_id,
                       &_shift_ids);
    }

  }

  /* Prepare return values */

  *shift_ids = _shift_ids;
}


/*----------------------------------------------------------------------------
 * Evaluate the box distribution over the leaves of the box tree when adding
 * a level to the tree structure.
 *
 * parameters:
 *   bt         <->  pointer to the box tree being built
 *   next_bt    <->  pointer to the next box tree being built
 *   boxes      <--  pointer to the associated box set structure
 *   node_id    <--  id of the starting node
 *   build_type <--  layout variant for building the tree structure
 *   shift_ids  <->  first free position free in new box_ids
 *----------------------------------------------------------------------------*/

static void
_build_next_level_cartesian(PDM_box_tree_t       *bt,
                           PDM_box_tree_t       *next_bt,
                           const PDM_box_set_t  *boxes,
                           int                   node_id,
                           PDM_box_tree_sync_t   build_type,
                           int                  *shift_ids)
{
  int   i;

  int   _shift_ids = *shift_ids;
  const _node_t  *cur_node = bt->local_data->nodes + node_id;

  if (cur_node->is_leaf == false) {

    assert(bt->local_data->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _build_next_level_cartesian(bt,
                                 next_bt,
                                 boxes,
                                 bt->local_data->child_ids[bt->n_children*node_id + i],
                                 build_type,
                                 &_shift_ids);
  }

  else { /* if (node->is_leaf == true) */

    if (   cur_node->n_boxes < bt->threshold
           && node_id != 0                    /* Root node is always divided */
           && build_type == PDM_BOX_TREE_ASYNC_LEVEL) {

      /* Copy related box_ids in the new next_ids */

      _node_t *next_node = next_bt->local_data->nodes + node_id;

      next_node->n_boxes = cur_node->n_boxes;
      next_node->start_id = _shift_ids;
      next_node->extra_weight = 0;

      for (i = 0; i < cur_node->n_boxes; i++)
        next_bt->local_data->box_ids[_shift_ids++]
          = bt->local_data->box_ids[cur_node->start_id + i];
    }
    else {  /* Split node and evaluate box distribution between its children */
      _split_node_cartesian_3d(bt,
                              next_bt,
                              boxes,
                              node_id,
                              &_shift_ids);
    }

  }

  /* Prepare return values */

  *shift_ids = _shift_ids;
}


/*----------------------------------------------------------------------------
 * Loop on all nodes of the box tree to define an array with Morton codes
 * and weights (= number of linked boxes) associated to each leaf
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   node_id    <-- id of the current node (to traverse)
 *   n_leaves   <-> current number of leaves in the tree with n_boxes > 0
 *   leaf_codes <-> Morton code associated to each leaf
 *   weight     <-> number of boxes attached to each leaf
 *----------------------------------------------------------------------------*/

static void
_build_leaf_weight(const PDM_box_tree_t    *bt,
                         int                node_id,
                         int                parent_weight,
                         int               *n_leaves,
                         PDM_morton_code_t *leaf_codes,
                         double            *weight)
{

  PDM_UNUSED (parent_weight);

  int  i;

  int _n_leaves = *n_leaves;

  const _node_t  *node = bt->local_data->nodes + node_id;

  if (node->is_leaf == false) {
    int repart_weight = (int) ceil((double) node->extra_weight/bt->n_children);
    // int repart_weight = bt->n_children;
    for (i = 0; i < bt->n_children; i++) {
      _build_leaf_weight(bt,
                         bt->local_data->child_ids[bt->n_children*node_id + i],
                         repart_weight,
                         &_n_leaves,
                         leaf_codes,
                         weight);
    }
  } else { /* node is a leaf */

    if (node->n_boxes > 0 ) {
      leaf_codes[_n_leaves]  = node->morton_code;
      weight    [_n_leaves]  = node->n_boxes;
      _n_leaves += 1;
    }
  }

  *n_leaves = _n_leaves;
}

/*----------------------------------------------------------------------------
 * Loop on all nodes of the tree to define an index of ranks related
 * to each box.
 *
 * parameters:
 *   bt           <-> pointer to PDM_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_index(const PDM_box_tree_t  *bt,
                         PDM_box_distrib_t     *distrib,
                         int                    dim,
                         int              node_id,
                         size_t                 size,
                         PDM_morton_code_t      search_index[],
                         int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->local_data->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_index(bt,
                               distrib,
                               dim,
                               bt->local_data->child_ids[bt->n_children*node_id + i],
                               size,
                               search_index,
                               id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = PDM_morton_binary_search((int) size,
                                         node->morton_code,
                                         search_index);
      int  rank = id_rank[id];

      distrib->index[rank + 1] += node->n_boxes;

    }
  }

}

/*----------------------------------------------------------------------------
 * Loop on all nodes of the tree to define a list of ranks related
 * to each box.
 *
 * parameters:
 *   bt           <-> pointer to PDM_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   counter      <-> counter array used to build the list
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_list(const PDM_box_tree_t  *bt,
                        PDM_box_distrib_t     *distrib,
                        int                    dim,
                        int              node_id,
                        int              counter[],
                        size_t                 size,
                        PDM_morton_code_t      search_index[],
                        int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->local_data->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_list(bt,
                              distrib,
                              dim,
                              bt->local_data->child_ids[bt->n_children*node_id + i],
                              counter,
                              size,
                              search_index,
                              id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = PDM_morton_binary_search((int) size,
                                         node->morton_code,
                                         search_index);
      int  rank = id_rank[id];

      for (i = 0; i < node->n_boxes; i++) {

        int   box_id = bt->local_data->box_ids[node->start_id + i];
        int   shift = distrib->index[rank] + counter[rank];

        distrib->list[shift] = box_id;
        counter[rank] += 1;

      }
    }
  }

}


static void
_count_box_ids(const PDM_box_tree_t  *bt,
                     int              dim,
                     int              node_id,
                     int             *n_box_ids)
{
  int  i;

  const _node_t  *node = bt->local_data->nodes + node_id;

  int _n_box_ids = *n_box_ids;
  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _count_box_ids(bt,
                     dim,
                     bt->local_data->child_ids[bt->n_children*node_id + i],
                     &_n_box_ids);
  }
  else {
    if (node->n_boxes > 0) {
      _n_box_ids += node->n_boxes;
    }
  }

  *n_box_ids = _n_box_ids;
}

static void
_collect_box_ids(const PDM_box_tree_t  *bt,
                       int              dim,
                       int              node_id,
                       int             *is_visited,
                       int             *box_ids,
                       int             *n_box_ids)
{
  int  i;

  const _node_t  *node = bt->local_data->nodes + node_id;

  int _n_box_ids = *n_box_ids;
  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _collect_box_ids(bt,
                       dim,
                       bt->local_data->child_ids[bt->n_children*node_id + i],
                       is_visited,
                       box_ids,
                       &_n_box_ids);
  }
  else {
    if (node->n_boxes > 0) {

      for (i = 0; i < node->n_boxes; i++) {

        int   box_id = bt->local_data->box_ids[node->start_id + i];

        if(is_visited[box_id] == 0) {
          is_visited[box_id] = 1;
          box_ids[_n_box_ids++] = box_id;
        }
      }
    }
  }

  *n_box_ids = _n_box_ids;
}

/*----------------------------------------------------------------------------
 * Recursively build an index on boxes which intersect.
 *
 * parameters:
 *   bt      <-> pointer to PDM_box_tree_t structure.
 *   boxesB  <-- pointer boxes that intersect associated tree boxes
 *   node_id <-- id of the current node (to traverse)
 *   count   <-> intersection count
 *----------------------------------------------------------------------------*/

static void
_count_boxes_intersections(const PDM_box_tree_t  *bt,
                           int *ncall,
                           int *n_boxes,
                           const PDM_box_set_t  *boxesB,
                           int                   boxB_id,
                           int                   node_id,
                           int                   count[])
{
  const PDM_box_tree_data_t *_local_data = bt->local_data;
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->local_boxes->extents;
  const _node_t  *node = _local_data->nodes + node_id;

  const double *boxB_extents = boxesB->local_boxes->extents + 2 * boxes->dim * boxB_id;

  *ncall += 1;

  /*assert (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 node_id));*/
  if (!_boxes_intersect_node (bt,
                              boxesB,
                              boxB_id,
                              node_id)) {
    printf("error _count_boxes_intersections :\n");
    const double  *box_min = _box_min(boxesB, boxB_id);
    const double  *box_max = _box_max(boxesB, boxB_id);
    printf("box : %f %f %f  / %f %f %f\n",
           box_min[0], box_min[1], box_min[2],
           box_max[0], box_max[1], box_max[2]);

    const PDM_morton_code_t code = bt->local_data->nodes[node_id].morton_code;
    printf("node : L = %d, X = %d %d %d\n",
           code.L, code.X[0], code.X[1], code.X[2]);

    const PDM_morton_int_t level = code.L;
    double min_grid_coord[3], max_grid_coord[3];
    _get_grid_coords_3d(level, box_min, min_grid_coord);
    _get_grid_coords_3d(level, box_max, max_grid_coord);
    printf("box grid coords: %f %f %f  / %f %f %f\n",
           min_grid_coord[0], min_grid_coord[1], min_grid_coord[2],
           max_grid_coord[0], max_grid_coord[1], max_grid_coord[2]);
    abort();
  }

  if (node->is_leaf == false) {

    for (int i = 0; i < bt->n_children; i++) { /* traverse downwards */

      if (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 _local_data->child_ids[bt->n_children*node_id + i])) {

        _count_boxes_intersections(bt,
                                   ncall,
                                   n_boxes,
                                   boxesB,
                                   boxB_id,
                                   _local_data->child_ids[bt->n_children*node_id + i],
                                   count);
      }
    }
  }

  else { /* node is a leaf */

    if (boxes->dim == 3) {

      *n_boxes += node->n_boxes;

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_3d (box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 2) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_2d (box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 1) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_1d(box_extents, id0, boxB_extents)) {
          count[id0] += 1;
        }
      }

    }
  }
}


/*----------------------------------------------------------------------------
 * Recursively build an index on boxes which intersect.
 *
 * parameters:
 *   bt      <-> pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *   count   <-> intersection count
 *----------------------------------------------------------------------------*/

static void
_count_intern_intersections(const PDM_box_tree_t  *bt,
                            int                   node_id,
                            int                   count[]) //***
{
  int   i, j;

  const PDM_box_tree_data_t *_local_data = bt->local_data;
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->local_boxes->extents;
  const _node_t  *node = _local_data->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _count_intern_intersections(bt,
                                  _local_data->child_ids[bt->n_children*node_id + i],
                                  count);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_3d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_2d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];
          if (_boxes_intern_intersect_1d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

  }

}


/*----------------------------------------------------------------------------
 * Recursively build a list on bounding boxes which intersect together.
 *
 * parameters:
 *   bt        <-> pointer to PDM_box_tree_t structure.
 *   boxesB    <-- pointer boxes that intersect associated tree boxes
 *   node_id   <-- id of the current node (to traverse)
 *   count     <-> intersection count (workin array)
 *   index     <-- index on intersections
 *   box_g_num <-> global number of intersection boxes
 *----------------------------------------------------------------------------*/

static void
_get_boxes_intersections(const PDM_box_tree_t  *bt,
                         const PDM_box_set_t  *boxesB,
                         int                   boxB_id,
                         int                   node_id,
                         int                   count[],
                         int                   box_index[],
                         int                   box_l_num[])
{
  const PDM_box_tree_data_t *_local_data = bt->local_data;
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->local_boxes->extents;
  const _node_t  *node = _local_data->nodes + node_id;

  const double *boxB_extents = boxesB->local_boxes->extents + 2 * boxes->dim * boxB_id;

  assert (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 node_id));

  if (node->is_leaf == false) {

    for (int i = 0; i < bt->n_children; i++) { /* traverse downwards */

      if (_boxes_intersect_node (bt,
                                 boxesB,
                                 boxB_id,
                                 _local_data->child_ids[bt->n_children*node_id + i])) {

        _get_boxes_intersections (bt,
                                  boxesB,
                                  boxB_id,
                                  _local_data->child_ids[bt->n_children*node_id + i],
                                  count,
                                  box_index,
                                  box_l_num);
      }
    }
  }

  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_3d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 2) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_2d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }

    else if (boxes->dim == 1) {

      for (int i = 0; i < node->n_boxes; i++) {
        int  id0 = _local_data->box_ids[node->start_id + i];
        if (_boxes_intersect_1d (box_extents, id0, boxB_extents)) {
          int   shift0 = box_index[id0] + count[id0];
          box_l_num[shift0] = boxB_id;
          count[id0] += 1;
        }
      }

    }
  } /* End if node is a leaf */
}


/*----------------------------------------------------------------------------
 * Recursively build a list on bounding boxes which intersect together.
 *
 * parameters:
 *   bt        <-> pointer to PDM_box_tree_t structure.
 *   node_id   <-- id of the current node (to traverse)
 *   count     <-> intersection count (workin array)
 *   index     <-- index on intersections
 *   box_g_num <-> global number of intersection boxes
 *----------------------------------------------------------------------------*/

static void
_get_intern_intersections(const PDM_box_tree_t  *bt,
                          int                    node_id,
                          int                    count[],
                          int                    box_index[],
                          PDM_g_num_t             box_g_num[])
{
  int  i, j;

  const PDM_box_tree_data_t *_local_data = bt->local_data;
  const PDM_box_set_t   *boxes = bt->boxes;
  const double  *box_extents = boxes->local_boxes->extents;
  const _node_t  *node = _local_data->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _get_intern_intersections(bt,
                                _local_data->child_ids[bt->n_children*node_id + i],
                                count,
                                box_index,
                                box_g_num);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_3d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->local_boxes->g_num[id1];
            box_g_num[shift1] = boxes->local_boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }
    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_2d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->local_boxes->g_num[id1];
            box_g_num[shift1] = boxes->local_boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          int   id0 = _local_data->box_ids[node->start_id + i];
          int   id1 = _local_data->box_ids[node->start_id + j];

          if (_boxes_intern_intersect_1d(box_extents, id0, id1)) {
            int   shift0 = box_index[id0] + count[id0];
            int   shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->local_boxes->g_num[id1];
            box_g_num[shift1] = boxes->local_boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }

  } /* End if node is a leaf */
}

/*----------------------------------------------------------------------------
 * Recursively define a counter array on the number of bounding boxes
 * associated to a leaf.
 *
 * This will be used for displaying a histogram.
 *
 * parameters:
 *   bt      <-- pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *   n_steps <-- number of steps in histogram
 *   step    <-- steps of the histogram
 *   h_min   <-- min. value of the histogram
 *   counter <-> counter (working array)
 *----------------------------------------------------------------------------*/

static void
_build_histogram(const PDM_box_tree_t  *bt,
                 int              node_id,
                 int              n_steps,
                 int              step,
                 int              h_min,
                 PDM_g_num_t              count[])
{
  int  i, j;

  const _node_t  *node = bt->local_data->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _build_histogram(bt,
                       bt->local_data->child_ids[bt->n_children*node_id + i],
                       n_steps,
                       step,
                       h_min,
                       count);
  }
  else {
    for (i = 0, j = 1; j < n_steps; i++, j++)
      if (node->n_boxes < h_min + j*step)
        break;
    count[i] += 1;
  }
}

/*----------------------------------------------------------------------------
 * Dump a box tree node.
 *
 * parameters:
 *   bt      <-- pointer to PDM_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *----------------------------------------------------------------------------*/

static void
_dump_node(const PDM_box_tree_t  *bt,
           int              node_id)
{
  int  i;

  const char *node_type[] = {"node", "leaf"};

  const PDM_box_tree_data_t *_local_data = bt->local_data;
  const _node_t  *node = _local_data->nodes + node_id;
  const PDM_morton_code_t  m_code = node->morton_code;

  PDM_printf("\n"
             "  node %10d (%s)\n"
             "    level:   %3u - anchor: [ %10u %10u %10u ]\n"
             "    n_boxes: %3d - start_id: %u\n"
             "    boxes:\n",
             node_id, node_type[(int)(node->is_leaf)],
	     m_code.L, m_code.X[0], m_code.X[1], m_code.X[2],
             node->n_boxes, node->start_id);

  for (i = 0; i < node->n_boxes; i++)
    PDM_printf("        %d\n", (int)(_local_data->box_ids[node->start_id + i]));

  if (node->is_leaf == false) {

    const int *c_id = _local_data->child_ids + bt->n_children*node_id;

    if (bt->n_children == 8) {
      PDM_printf("  children_id:  %d %d %d %d %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3],
                 (int)c_id[4], (int)c_id[5], (int)c_id[6], (int)c_id[7]);
    }
    else if (bt->n_children == 4) {
      PDM_printf("  children_id:  %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3]);
    }
    else if (bt->n_children == 2) {
      PDM_printf("  children_id:  %d %d\n",
                 (int)c_id[0], (int)c_id[1]);
    }

    for (i = 0; i < bt->n_children; i++)
      _dump_node(bt, c_id[i]);
  }
}

/*
 * \brief Determine if a box intersects a given volume region
 *
 * \param [in]  n_planes      Number of planes difining the volume
 * \param [in]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [in]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 * \param [in]  box_extents   Points to determine
 *
 * \return 1 if the plane_pt is on the side of the plane where the normal points to, 0 otherwise
 */


static int
_box_intersect_volume
(
int      n_planes,
double  *n,
double  *plane_pt,
double  *box_extents
)
{
  double box_pt[3];
  double vect[3];
  double n_iplane[3];
  double plane_pt_iplane[3];

  int count_intersected_planes, count_points_not_intersect_plane;


  // All planes for one point
  for (int x = 0; x < 4; x += 3) {
    box_pt[0] = box_extents[x];
    for (int y = 1; y < 5; y += 3) {
      box_pt[1] = box_extents[y];
      for (int z = 2; z < 6; z += 3) {
        box_pt[2] = box_extents[z];

        count_intersected_planes = 0;

        for (int iplane = 0; iplane < n_planes; iplane++) {

          n_iplane[0] = n[3*iplane];
          n_iplane[1] = n[3*iplane+1];
          n_iplane[2] = n[3*iplane+2];

          plane_pt_iplane[0] = plane_pt[3*iplane];
          plane_pt_iplane[1] = plane_pt[3*iplane+1];
          plane_pt_iplane[2] = plane_pt[3*iplane+2];

          vect[0] = box_pt[0] - plane_pt_iplane[0]; vect[1] = box_pt[1] - plane_pt_iplane[1]; vect[2] = box_pt[2] - plane_pt_iplane[2];

          if (PDM_DOT_PRODUCT(vect, n_iplane) >= 0) {
            count_intersected_planes++;
          }

        } // end loop on planes

        if (count_intersected_planes == n_planes) {
          return 1;
        }

      }
    }
  }

  // All points for one plane
  for (int iplane = 0; iplane < n_planes; iplane++) {

    count_points_not_intersect_plane = 0;

    for (int x = 0; x < 4; x += 3) {
      box_pt[0] = box_extents[x];
      for (int y = 1; y < 5; y += 3) {
        box_pt[1] = box_extents[y];
        for (int z = 2; z < 6; z += 3) {
          box_pt[2] = box_extents[z];

          n_iplane[0] = n[3*iplane];
          n_iplane[1] = n[3*iplane+1];
          n_iplane[2] = n[3*iplane+2];

          plane_pt_iplane[0] = plane_pt[3*iplane];
          plane_pt_iplane[1] = plane_pt[3*iplane+1];
          plane_pt_iplane[2] = plane_pt[3*iplane+2];

          vect[0] = box_pt[0] - plane_pt_iplane[0]; vect[1] = box_pt[1] - plane_pt_iplane[1]; vect[2] = box_pt[2] - plane_pt_iplane[2];

          if (PDM_DOT_PRODUCT(vect, n_iplane) < 0) {
            count_points_not_intersect_plane++;
          }

        }
      }
    }

    if (count_points_not_intersect_plane == 8) {
      return 0;
    }

  }

  // Undefined case
  return PDM_lp_intersect_volume_box(n_planes, plane_pt, n, box_extents);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a \ref PDM_box_tree_data_t structure and initialize it
 *
 * \return Pointer to an empty \ref PDM_box_tree_data_t structure
 *
 */

PDM_box_tree_data_t *
PDM_box_tree_data_create(void)
{
  PDM_box_tree_data_t  *btd = NULL;

  btd = (PDM_box_tree_data_t *) malloc(sizeof(PDM_box_tree_data_t));

  btd->n_max_nodes = 0;
  btd->n_nodes = 0;

  btd->nodes = NULL;

  btd->box_ids = NULL;

  btd->n_build_loops = 0;

  btd->stack     = NULL;
  btd->pos_stack = NULL;

  return btd;
}



/**
 *
 * \brief Create a \ref PDM_box_tree_t structure and initialize it
 *
 * \param [in]   max_level      Max possible level
 * \param [in]   threshold      Max number of boxes linked to an octant if \ref max_level is not reached
 * \param [in]   max_box_ratio  Max n_linked_boxes / n_boxes ratio
 *
 * \return Pointer to an empty \ref PDM_box_tree_t structure
 *
 */

PDM_box_tree_t *
PDM_box_tree_create(int    max_level,
                    int    threshold,
                    float  max_box_ratio)
{
  PDM_box_tree_t  *bt = NULL;

  bt = (PDM_box_tree_t *) malloc(sizeof(PDM_box_tree_t));

  /* Sanity checks */

  if (max_level < 0) {
    PDM_error(__FILE__, __LINE__, 0,
	      "  Forbidden max_level value (%d) in the tree structure\n",
	      max_level);
    abort();
  }

  if (threshold < 1) {
    PDM_error(__FILE__, __LINE__, 0,
	      "  Forbidden threshold value (%d) in the tree structure\n",
	      threshold);
    abort();
  }

  if (max_box_ratio < 1.0) {
    PDM_error(__FILE__, __LINE__, 0,
	      "  Forbidden max_box_ratio value (%f) in the tree structure\n",
	      (double)max_box_ratio);
    abort();
  }

  /* Create and initialize tree structure according to its type */

  bt->max_level = max_level;
  bt->threshold = threshold;
  bt->max_box_ratio = max_box_ratio;

  bt->comm = PDM_MPI_COMM_NULL;

  /* Set stats */

  bt->stats.max_level_reached = 0;

  bt->stats.n_leaves = 0;
  bt->stats.n_spill_leaves = 0;
  bt->stats.n_linked_boxes = 0;

  bt->stats.min_linked_boxes = INT_MAX;
  bt->stats.max_linked_boxes = 0;

  /* Initialize nodes */

  /*//-->> fonction PDM_box_tree_data_create
    bt->n_max_nodes = 0;
    bt->n_nodes = 0;

    bt->nodes = NULL;

    bt->box_ids = NULL;

    bt->n_build_loops = 0;

    bt->stack = NULL;
    bt->pos_stack = NULL;
    //<<---*/
  bt->local_data = PDM_box_tree_data_create();
  //bt->local_data = NULL;

  bt->n_copied_ranks = 0;
  bt->copied_ranks   = NULL;
  bt->rank_data      = NULL;

  /* Shared memory */
  bt->n_rank_in_shm  = 0;
  bt->shm_data       = NULL;
  bt->wbox_tree_data = NULL;

  return bt;
}



/**
 *
 * \brief Destroy a \ref PDM_box_tree_data_t structure
 *
 * \param [in]   btd  Pointer to pointer to \ref PDM_box_tree_data_t structure to destroy
 *
 */

void
PDM_box_tree_data_destroy(PDM_box_tree_data_t  *btd)
{
  _free_tree_data_arrays(btd);

  if (btd->stack != NULL) {
    free(btd->stack);
  }

  if (btd->pos_stack != NULL) {
    free(btd->pos_stack);
  }
}



/**
 *
 * \brief Destroy a \ref PDM_box_tree_t structure
 *
 * \param [in]   bt   Pointer to pointer to \ref PDM_box_tree_t structure to destroy
 *
 */

void
PDM_box_tree_destroy(PDM_box_tree_t  **bt)
{
  PDM_box_tree_t  *_bt = *bt;

  if (_bt == NULL) {
    return;
  }

  if(_bt->shm_data != NULL) {
    for(int i = 0; i < _bt->n_rank_in_shm; ++i) {
      PDM_mpi_win_shared_unlock_all (_bt->wbox_tree_data[i].w_nodes    );
      PDM_mpi_win_shared_unlock_all (_bt->wbox_tree_data[i].w_child_ids);
      PDM_mpi_win_shared_unlock_all (_bt->wbox_tree_data[i].w_box_ids  );

      PDM_mpi_win_shared_free(_bt->wbox_tree_data[i].w_nodes    );
      PDM_mpi_win_shared_free(_bt->wbox_tree_data[i].w_child_ids);
      PDM_mpi_win_shared_free(_bt->wbox_tree_data[i].w_box_ids  );
    }
    free(_bt->wbox_tree_data);
    free(_bt->shm_data);
  }

  // Free local tree data
  PDM_box_tree_data_destroy(_bt->local_data);
  free(_bt->local_data);

  // Free copied ranks tree data
  // PDM_box_tree_free_copies(_bt);

  free(_bt);
  *bt = NULL;
}



/**
 *
 * \brief Get the deepest level allowed by the tree structure
 *
 * \param [in]   bt   Pointer to pointer to PDM_box_tree_t structure to destroy
 *
 * \return  Deepest allowed level of the tree
 *
 */

int
PDM_box_tree_get_max_level(const PDM_box_tree_t  *bt)
{
  return bt->max_level;
}



/**
 *
 * \brief Assign a set of boxes to an empty \ref PDM_box_tree_t structure
 *
 * The box tree structure must have been created using to PDM_tree_create().
 *
 * The depth of the tree is adjusted so that a maximum of max_n_elts boxes
 * will be assigned to each leaf, unless this would require going beyond
 * the tree's maximum level.
 *
 * If max_level = -1, the highest level reachable is PDM_TREE_MAX_LEVEL but
 * there is no defined target level.
 *
 * \param [in]   bt           Pointer to \ref PDM_box_tree_t structure
 * \param [in]   boxes        Pointer to the associated box set structure
 * \param [in]   build_type   Layout variant for building the tree structure
 *
 */

void
PDM_box_tree_set_boxes(PDM_box_tree_t      *bt,
                       PDM_box_set_t       *boxes,
                       PDM_box_tree_sync_t  build_type)
{
  //printf("  -->> PDM_box_tree_set_boxes\n");
  PDM_boxes_t *_local_boxes = boxes->local_boxes;

  double t1 = PDM_MPI_Wtime();

  int   box_id;

  PDM_box_tree_t  tmp_bt;

  int   next_box_ids_size = 0, shift = 0;
  double anchor[3] = {0., 0., 0.};

  /* Initialization */

  assert(bt != NULL);

  PDM_box_tree_data_t *_local_data = bt->local_data;
  assert( _local_data != NULL );

  _local_data->n_build_loops = 0;//*

  bt->comm = boxes->comm;
  bt->boxes = boxes;

  /* Preallocate for the two first levels of a tree */

  if (boxes->dim == 3) {
    bt->n_children = 8;
    _local_data->n_max_nodes = 73;//*
  }
  else if (boxes->dim == 2) {
    bt->n_children = 4;
    _local_data->n_max_nodes = 21;//*
  }
  else if (boxes->dim == 1) {
    bt->n_children = 2;
    _local_data->n_max_nodes = 7;//*
  }

  _local_data->n_nodes = 1;//*

  _local_data->nodes = (_node_t *) malloc(_local_data->n_max_nodes * sizeof(_node_t));//*
  _local_data->child_ids = (int *) malloc(_local_data->n_max_nodes*bt->n_children * sizeof(int));//*

  /* Define root node */

  _new_node(bt, PDM_morton_encode(boxes->dim, 0, anchor), 0);

  /* Initialize bt by assigning all boxes to the root leaf */

  _local_data->box_ids = (int *) malloc (_local_boxes->n_boxes * sizeof(int));

  for (box_id = 0; box_id < _local_boxes->n_boxes; box_id++)
    _local_data->box_ids[box_id] = box_id;

  (_local_data->nodes[0]).is_leaf = true;
  (_local_data->nodes[0]).n_boxes = _local_boxes->n_boxes;
  (_local_data->nodes[0]).start_id = 0;
  (_local_data->nodes[0]).extra_weight = 0;

  bt->stats.n_boxes = _local_boxes->n_boxes;

  _get_box_tree_stats(bt);

  /* printf ("extents set_boxes :"); */
  /* for (int i = 0; i < _local_boxes->n_boxes; i++) { */
  /*   printf("%12.5e < %12.5e / %12.5e < %12.5e / %12.5e < %12.5e\n", */
  /*          _local_boxes->extents[6*i  ], */
  /*          _local_boxes->extents[6*i+3], */
  /*          _local_boxes->extents[6*i+1], */
  /*          _local_boxes->extents[6*i+4], */
  /*          _local_boxes->extents[6*i+2], */
  /*          _local_boxes->extents[6*i+5] */
  /*          ); */
  /* } */

  char *env_var_oct = getenv ("BOX_TREE_CARTESIAN_BUILD");
  int build_cartesian = 0;
  if (env_var_oct != NULL) {
    build_cartesian = atoi(env_var_oct);
  }

  /* Build local tree structure by adding boxes from the root */
  while (_recurse_tree_build(bt,
                             boxes,
                             build_type,
                             build_cartesian,
                             &next_box_ids_size)) {
    /* Initialize next_bt: copy of bt */
    _copy_tree(&tmp_bt, bt);
    /* Optimize memory usage */

    bt->local_data->n_max_nodes = bt->local_data->n_nodes;
    bt->local_data->nodes = (_node_t *) realloc((void *) bt->local_data->nodes, bt->local_data->n_nodes * sizeof(_node_t));
    bt->local_data->child_ids = (int *) realloc((void *) bt->local_data->child_ids,
						bt->local_data->n_max_nodes*bt->n_children * sizeof(int));

    /* Define a box ids list for the next level of the boxtree */
    tmp_bt.local_data->box_ids = (int*) realloc((void *) tmp_bt.local_data->box_ids, next_box_ids_size * sizeof(int));
    shift = 0;

    if(build_cartesian == 0) {
      _build_next_level(bt,
                        &tmp_bt,
                        boxes,
                        0, /* Starts from root */
                        build_type,
                        &shift);
     } else {
      _build_next_level_cartesian(bt,
                                 &tmp_bt,
                                 boxes,
                                 0, /* Starts from root */
                                 build_type,
                                 &shift);
    }
    assert(shift == next_box_ids_size);

    /* replace current tree by the tree computed at a higher level */
    _free_tree_arrays(bt);
    free(bt->local_data);
    *bt = tmp_bt; /* Overwrite bt members with those of next_bt */
    _get_box_tree_stats(bt);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    PDM_printf("  - New box tree level -\n");
    PDM_box_tree_dump_statistics(bt);
#endif

  } /* While building should continue */
  //printf("  <<-- PDM_box_tree_set_boxes\n");


  double dt = PDM_MPI_Wtime() - t1;
  if(0) {
    log_trace("PDM_box_tree_set_boxes : %12.5e \n", dt);
  }


}



/**
 *
 * \brief Compute an index based on Morton encoding to ensure a good distribution
 * of boxes among the participating ranks
 *
 * \param [in]   bt           Pointer to \ref PDM_box_tree_t structure
 * \param [in]   boxes        Pointer to the associated box set structure
 *
 * \return   Pointer to newly created \ref PDM_box_distrib_t structure
 *
 */

PDM_box_distrib_t *
PDM_box_tree_get_distrib(PDM_box_tree_t        *bt,
                         const PDM_box_set_t   *boxes)
{
  int  i;

  int  reduce_size = 0;
  int   n_leaves = 0;
  int  *reduce_ids = NULL;
  PDM_morton_code_t  *leaf_codes = NULL, *reduce_index = NULL;
  double   *weight = NULL;
  int      *counter = NULL;

  PDM_box_distrib_t  *distrib = NULL;

  assert(bt != NULL);
  assert(boxes != NULL);

  /* Compute basic box distribution */

  distrib = PDM_box_distrib_create(boxes->local_boxes->n_boxes,
                                   boxes->n_g_boxes,
                                   (bt->stats).max_level_reached,
                                   boxes->comm);

  if (distrib == NULL)
    return NULL;

  leaf_codes = (PDM_morton_code_t *) malloc(bt->stats.n_leaves * sizeof(PDM_morton_code_t));
  weight = (double *) malloc(bt->stats.n_leaves * sizeof(double));

  /* Build index for boxes */
  int repart_weight = 0;
  _build_leaf_weight(bt,
                     0,
                     repart_weight,
                     &n_leaves,
                     leaf_codes,
                     weight);

  assert(n_leaves <= bt->stats.n_leaves);

  // PDM_log_trace_array_int(weight, n_leaves, "weight : ");

  leaf_codes = (PDM_morton_code_t *) realloc((void *) leaf_codes, n_leaves * sizeof(PDM_morton_code_t));
  weight = (double *) realloc((void *) weight, n_leaves * sizeof(double));

  /* Compute the resulting Morton index */

  PDM_box_set_build_morton_index(boxes,
                                 distrib,
                                 n_leaves,
                                 leaf_codes,
                                 weight);

  free(leaf_codes);
  free(weight);

  /* Compact Morton_index to get an array without "0 element" */

  for (i = 0; i < distrib->n_ranks; i++)
    if (PDM_morton_a_gt_b(distrib->morton_index[i+1],
                          distrib->morton_index[i]))
      reduce_size++;

  reduce_index = (PDM_morton_code_t *) malloc((reduce_size + 1) * sizeof(PDM_morton_code_t));
  reduce_ids  = (int *) malloc(reduce_size * sizeof(int));

  reduce_size = 0;
  reduce_index[0] = distrib->morton_index[0];

  for (i = 0; i < distrib->n_ranks; i++) {

    if (PDM_morton_a_gt_b(distrib->morton_index[i+1],
                          distrib->morton_index[i])) {

      reduce_index[reduce_size + 1] = distrib->morton_index[i+1];
      reduce_ids[reduce_size++] = i;

    }

  }

  /* Define a rank -> box indexed list */

  _build_rank_to_box_index(bt,
                           distrib,
                           boxes->dim,
                           0,  /* starts from root */
                           reduce_size,
                           reduce_index,
                           reduce_ids);

  for (i = 0; i < distrib->n_ranks; i++)
    distrib->index[i+1] += distrib->index[i];

 distrib->list = (int *) malloc(distrib->index[distrib->n_ranks] * sizeof(int));

  counter = PDM_array_zeros_int(distrib->n_ranks);

  _build_rank_to_box_list(bt,
                          distrib,
                          boxes->dim,
                          0,  /* starts from root */
                          counter,
                          reduce_size,
                          reduce_index,
                          reduce_ids);

  /* Free memory */

  free(counter);
  free(reduce_ids);
  free(reduce_index);

  /* Define the final index (without redundancies) and realloc list */

  PDM_box_distrib_clean(distrib);

  return distrib;
}



/**
 *
 * \brief Build an indexed list on associated bt boxes to list boxes B intersections
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * \param [in]   bt          Pointer to box tree structure to query
 * \param [in]   boxes       Pointer to boxes that intersect associated tree boxes
 * \param [out]  box_index   Pointer to the index array on bounding boxes
 * \param [out]  box_g_num   Pointer to the list of intersecting bounding boxes
 *
 */

void
PDM_box_tree_get_boxes_intersects(PDM_box_tree_t       *bt,
                                  const PDM_box_set_t  *boxesB,
                                  int                  *box_index[],
                                  int                  *box_l_num[])
{
  int  i, list_size;

  int  *counter = NULL;
  int  *_index = NULL;
  int         *_l_num = NULL;
  const PDM_box_set_t  *boxes = bt->boxes;

  int myRank;
  PDM_MPI_Comm_rank(bt->comm, &myRank);

  /* Build index */

  _index = PDM_array_zeros_int(boxes->local_boxes->n_boxes + 1);

  int ncall=0;
  int n_boxes =0;

  for (int k = 0; k < boxesB->local_boxes->n_boxes; k++) {
    ncall = 0;
    n_boxes = 0;
    _count_boxes_intersections(bt,
                               &ncall,
                               &n_boxes,
                               boxesB,
                               k,
                               0, /* start from root */
                               _index + 1);

  // log_trace("ncall   =  %i \n", ncall);
  // log_trace("n_boxes =  %i \n", n_boxes);

  }

  /* Build index from counts */

  for (i = 0; i < boxes->local_boxes->n_boxes; i++)
    _index[i+1] += _index[i];

  list_size = _index[boxes->local_boxes->n_boxes];

  _l_num = (int *) malloc(list_size * sizeof(int));

  counter = PDM_array_zeros_int(boxes->local_boxes->n_boxes);

  /* Build list */

  for (int k = 0; k < boxesB->local_boxes->n_boxes; k++) {
    _get_boxes_intersections(bt,
                             boxesB,
                             k,
                             0, /* start from root */
                             counter,
                             _index,
                             _l_num);
  }

  /* Remove duplicate boxes */

  int idx=0;
  for (i = 0; i < boxes->local_boxes->n_boxes; i++) {
    counter[i] = 0;
    int *_l_num_box = _l_num + _index[i];
    int n_elt = _index[i+1] -_index[i];
    counter[i] = 0;
    _index[i] = 0;
    PDM_sort_int(_l_num_box, NULL, n_elt);

    int pre = -1;
    for (int j = 0; j < n_elt; j++) {
      if (pre < _l_num_box[j]) {
        pre = _l_num_box[j];
        _l_num[idx++] = pre;
        counter[i]++;
      }
    }
  }

  for (i = 0; i < boxes->local_boxes->n_boxes; i++) {
    _index[i+1] = _index[i] + counter[i];
  }

  _l_num = realloc (_l_num, sizeof(int) * _index[boxes->local_boxes->n_boxes]);

  free(counter);

  /* Return pointers */

  *box_index = _index;
  *box_l_num = _l_num;
}



/**
 *
 * \brief Build an indexed list on boxes to list intersections
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * \param [in]   bt          Pointer to box tree structure to query
 * \param [out]  box_index   Pointer to the index array on bounding boxes
 * \param [out]  box_g_num   Pointer to the list of intersecting bounding boxes
 *
 */

void
PDM_box_tree_get_intern_intersects(PDM_box_tree_t       *bt,
                                   int                  *box_index[],
                                   PDM_g_num_t           *box_g_num[])
{
  int  i, list_size;

  int  *counter = NULL;
  int  *_index = NULL;
  PDM_g_num_t  *_g_num = NULL;
  const PDM_box_set_t  *boxes = bt->boxes;
  const PDM_boxes_t *_local_boxes = boxes->local_boxes;

  /* Build index */

  _index = PDM_array_zeros_int(_local_boxes->n_boxes + 1);

  _count_intern_intersections(bt,
                              0, /* start from root */
                              _index + 1);

  /* Build index from counts */

  for (i = 0; i < _local_boxes->n_boxes; i++)
    _index[i+1] += _index[i];

  list_size = _index[_local_boxes->n_boxes];

  _g_num = (PDM_g_num_t *) malloc(list_size * sizeof(PDM_g_num_t));

  counter = PDM_array_zeros_int(_local_boxes->n_boxes);

  /* Build list */

  _get_intern_intersections(bt,
                            0, /* start from root */
                            counter,
                            _index,
                            _g_num);

  free(counter);

  /* Return pointers */

  *box_index = _index;
  *box_g_num = _g_num;
}



/**
 *
 * \brief Get global box tree statistics
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
 * Note that the theoretical memory use includes that of the associated
 * box set.
 *
 * \param [in]  bt                  Pointer to box tree structure
 * \param [in]  depth               Tree depth (max level used)
 * \param [out] n_leaves            Number of leaves in the tree
 * \param [out] n_boxes             Number of boxes in the tree
 * \param [out] n_threshold_leaves  Number of leaves where n_boxes > threshold
 * \param [out] n_leaf_boxes        Number of boxes for a leaf
 * \param [out] mem_used            Theoretical used memory
 * \param [out] mem_allocated       Theoretical allocated memory
 *
 * \return The spatial dimension associated with the box tree layout (3, 2, or 1)
 *
 */

int
PDM_box_tree_get_stats(const PDM_box_tree_t *bt,
                       int                   depth[3],
                       int                   n_leaves[3],
                       int                   n_boxes[3],
                       int                   n_threshold_leaves[3],
                       int                   n_leaf_boxes[3],
                       size_t                mem_used[3],
                       size_t                mem_allocated[3])
{
  int i;
  uint64_t mem_per_node;
  uint64_t s_mean[7], s_min[7], s_max[7];
  PDM_box_tree_stats_t s;

  int dim = 3;

  if (bt == NULL)
    return 0;

  s = bt->stats;

  if (bt->n_children == 4)
    dim = 2;
  else if (bt->n_children == 2)
    dim = 1;

  /* Prepare array of local values; prior to or in the absence of
     MPI communication, mean values are set to local values. */

  s_mean[0] = s.n_linked_boxes / s.n_leaves;
  /* Round to nearest integer, and not floor */
  if (s.n_linked_boxes % s.n_leaves >= s.n_leaves/2)
    s_mean[0] += 1;

  s_min[0] = s.min_linked_boxes;
  s_max[0] = s.max_linked_boxes;

  s_mean[1] = s.max_level_reached;
  s_mean[2] = s.n_leaves;
  s_mean[3] = s.n_boxes;
  s_mean[4] = s.n_spill_leaves;

  /* Estimate theoretical memory usage */

  mem_per_node = sizeof(_node_t) + bt->n_children*sizeof(int);

  s_mean[5] = sizeof(PDM_box_tree_t);
  s_mean[5] += bt->local_data->n_nodes * mem_per_node;
  s_mean[5] += s.n_linked_boxes * sizeof(int);

  s_mean[5] += sizeof(PDM_box_set_t);
  s_mean[5] += s.n_boxes * (  sizeof(PDM_g_num_t)
			      + (dim * 2 * sizeof(double)));

  s_mean[6] = s_mean[5] + (bt->local_data->n_max_nodes - bt->local_data->n_nodes)*mem_per_node;

  /* Pre-synchronize for serial cases (i = 0 already handled) */

  for (i = 1; i < 7; i++) {
    s_min[i] = s_mean[i];
    s_max[i] = s_mean[i];
  }

  /* In parallel mode, synchronize values */

  if (bt->comm != PDM_MPI_COMM_NULL) {

    int n_ranks;
    uint64_t s_l_sum[14], s_g_sum[14];

    PDM_MPI_Comm_size(bt->comm, &n_ranks);

    if (n_ranks > 1) { /* Should always be the case, bat play it safe) */

      /* Split value to avoid exceeding PDM_g_num_t limits
         (especially when it is a 32-bit value) */

      s_l_sum[0] = s.n_linked_boxes/n_ranks;
      s_l_sum[7] = s.n_linked_boxes%n_ranks;
      for (i = 1; i < 7; i++) {
        s_l_sum[i] = ((uint64_t) s_mean[i])/n_ranks;
        s_l_sum[i+7] = ((uint64_t) s_mean[i])%n_ranks;
      }

      PDM_MPI_Allreduce(s_l_sum, s_g_sum, 14, PDM_MPI_UINT64_T, PDM_MPI_SUM, bt->comm);

      s_mean[0] = s.min_linked_boxes;
      PDM_MPI_Allreduce(s_mean, s_min, 7, PDM_MPI_UINT64_T, PDM_MPI_MIN, bt->comm);
      s_mean[0] = s.max_linked_boxes;
      PDM_MPI_Allreduce(s_mean, s_max, 7, PDM_MPI_UINT64_T, PDM_MPI_MAX, bt->comm);

      /* Specific handling for linked boxes, so as to ensure correct
         total using large integers even if we do not know the
         corresponding MPI type */
      {
        uint64_t s_n = s_g_sum[0]*n_ranks + s_g_sum[7]; /* linked boxes */
        uint64_t s_d = s_g_sum[2]*n_ranks + s_g_sum[9]; /* leaves */
        uint64_t s_m = s_n / s_d;
        /* Round to nearest integer, and not floor */
        if (s_n % s_d >= s_d/2)
          s_m += 1;
        s_mean[0] = s_m;
      }
      for (i = 1; i < 7; i++) {
        s_mean[i] = s_g_sum[i] + s_g_sum[i+7]/n_ranks;
        /* Round to nearest integer, and not floor */
        if ((PDM_g_num_t)(s_g_sum[i+7]%n_ranks) >= (PDM_g_num_t) (n_ranks/2))
          s_mean[i] += 1;
      }
    }

  }

  /* Set values already in stats */

  if (depth != NULL) {
    depth[0] = (int) s_mean[1];
    depth[1] = (int) s_min[1];
    depth[2] = (int) s_max[1];
  }

  if (n_leaves != NULL) {
    n_leaves[0] = (int) s_mean[2];
    n_leaves[1] = (int) s_min[2];
    n_leaves[2] = (int) s_max[2];
  }

  if (n_boxes != NULL) {
    n_boxes[0] = (int) s_mean[3];
    n_boxes[1] = (int) s_min[3];
    n_boxes[2] = (int) s_max[3];
  }

  if (n_threshold_leaves != NULL) {
    n_threshold_leaves[0] = (int) s_mean[4];
    n_threshold_leaves[1] = (int) s_min[4];
    n_threshold_leaves[2] = (int) s_max[4];
  }

  if (n_leaf_boxes != NULL) {
    n_leaf_boxes[0] = (int) s_mean[0];
    n_leaf_boxes[1] = (int) s_min[0];
    n_leaf_boxes[2] = (int) s_max[0];
  }

  if (mem_used != NULL) {
    mem_used[0] = (int) s_mean[5];
    mem_used[1] = (int) s_min[5];
    mem_used[2] = (int) s_max[5];
  }

  if (mem_allocated != NULL) {
    mem_allocated[0] = (int) s_mean[6];
    mem_allocated[1] = (int) s_min[6];
    mem_allocated[2] = (int) s_max[6];
  }

  return dim;
}



/**
 *
 * \brief Display local statistics about a \ref PDM_box_tree_t structure
 *
 * \param [in]  bt   Pointer to box tree structure
 *
 */

void
PDM_box_tree_dump_statistics(const PDM_box_tree_t  *bt)
{
  int i, j;
  PDM_box_tree_stats_t s;
  unsigned g_max_level_reached;
  PDM_g_num_t n_g_leaves, n_g_boxes, n_g_linked_boxes, n_g_spill_leaves;
  int g_min_linked_boxes, g_max_linked_boxes;
  double mean_linked_boxes, box_ratio;

  PDM_g_num_t count[5];

  int step = 0, delta = 0;
  const int n_steps = 5;

  if (bt == NULL)
    return;

  s = bt->stats;

  g_max_level_reached = s.max_level_reached;
  n_g_leaves = s.n_leaves;
  n_g_boxes = s.n_boxes;
  n_g_linked_boxes = s.n_linked_boxes;
  n_g_spill_leaves = s.n_spill_leaves;
  g_min_linked_boxes = s.min_linked_boxes;
  g_max_linked_boxes = s.max_linked_boxes;

  if (bt->comm != PDM_MPI_COMM_NULL) {

    PDM_g_num_t l_min[1], g_min[1];
    PDM_g_num_t l_max[2], g_max[2];
    PDM_g_num_t l_sum[3], g_sum[3];

    l_sum[0] = n_g_leaves;
    l_sum[1] = n_g_spill_leaves;
    l_sum[2] = n_g_linked_boxes;

    l_min[0] = g_min_linked_boxes;
    l_max[0] = s.max_level_reached;
    l_max[1] = g_max_linked_boxes;

    PDM_MPI_Allreduce(l_sum, g_sum, 3, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, bt->comm);
    PDM_MPI_Allreduce(l_min, g_min, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MIN, bt->comm);
    PDM_MPI_Allreduce(l_max, g_max, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, bt->comm);

    n_g_leaves = l_sum[0];
    n_g_spill_leaves = l_sum[1];
    n_g_linked_boxes = l_sum[2];

    g_min_linked_boxes = (int) g_min[0];
    g_max_level_reached = (int) g_max[0];
    g_max_linked_boxes = (int) g_max[1];
  }

  /* Redefine final statistics */

  double _n_g_linked_boxes = (double) n_g_linked_boxes;
  double _n_g_leaves = (double) n_g_leaves;
  double _n_g_boxes = (double) n_g_boxes;

  mean_linked_boxes = _n_g_linked_boxes / _n_g_leaves;
  box_ratio = _n_g_linked_boxes / _n_g_boxes;

  /* Define axis subdivisions */
  PDM_array_reset_gnum(count, n_steps, 0);

  delta = g_max_linked_boxes - g_min_linked_boxes;

  if (delta > 0) {

    step = delta/n_steps;

    _build_histogram(bt,
                     0, /* start from root */
                     n_steps,
                     step,
                     g_min_linked_boxes,
                     count);

  } /* max - min > 0 */

  /* Print statistics and bounding boxes histogram */

  PDM_printf("\n"
             "Box tree statistics:\n\n");
  PDM_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (final/init):     %f\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (double)(bt->max_box_ratio), (int)(bt->max_level));

  PDM_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10llu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10llu\n"
             "  Initial number of boxes:            %10llu\n"
             "  Number of linked boxes:             %10llu\n"
             "  Mean number of leaves per box:      %10.4g\n\n",
             g_max_level_reached, (unsigned long long)n_g_leaves,
             (unsigned long long)n_g_spill_leaves, (unsigned long long)n_g_boxes,
             (unsigned long long)n_g_linked_boxes,  box_ratio);

  PDM_printf("Number of linked boxes per box tree leaf:\n"
             "  Mean value:         %10.4g\n"
             "  min. value:         %10llu\n"
             "  max. value:         %10llu\n\n",
             mean_linked_boxes,
             (unsigned long long)(s.min_linked_boxes),
             (unsigned long long)(s.max_linked_boxes));

  if (delta > 0) { /* Number of elements in each subdivision */

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      PDM_printf("    %3d : [ %10llu; %10llu [ = %10llu\n",
                 i+1,
                 (unsigned long long)(g_min_linked_boxes + i*step),
                 (unsigned long long)(g_min_linked_boxes + j*step),
                 (unsigned long long)(count[i]));

    PDM_printf("    %3d : [ %10llu; %10llu ] = %10llu\n",
               n_steps,
               (unsigned long long)(g_min_linked_boxes + (n_steps - 1)*step),
               (unsigned long long)(g_max_linked_boxes),
               (unsigned long long)(count[n_steps - 1]));

  }
}



/**
 *
 * \brief Dump a \ref PDM_box_tree_t structure
 *
 * \param [in]  bt   Pointer to box tree structure
 *
 */

void
PDM_box_tree_dump(PDM_box_tree_t  *bt)
{
  PDM_box_tree_stats_t s;

  if (bt == NULL) {
    PDM_printf("\nBox tree: nil\n");
    return;
  }

  PDM_printf("\nBox tree: %p\n\n", (void *)bt);

  PDM_printf("  n_max_nodes:  %d\n\n"
             "  n_nodes:      %d\n",
             (int)(bt->local_data->n_max_nodes), (int)(bt->local_data->n_nodes));

  s = bt->stats;

  /* Print statistics and bounding boxes histogram */

  PDM_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (linked/init):    %f\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (double)(bt->max_box_ratio), (int)(bt->max_level));

  PDM_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10llu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10llu\n"
             "  Initial number of boxes:            %10llu\n"
             "  Number of linked boxes:             %10llu\n",
             s.max_level_reached,
             (unsigned long long)(s.n_leaves),
             (unsigned long long)(s.n_spill_leaves),
             (unsigned long long)(s.n_boxes),
             (unsigned long long)(s.n_linked_boxes));

  PDM_printf("Bounding boxes related to each leaf of the box tree.\n"
             "  min. value:         %10llu\n"
             "  max. value:         %10llu\n\n",
             (unsigned long long)(s.min_linked_boxes),
             (unsigned long long)(s.max_linked_boxes));

  _dump_node(bt, 0);
}



/**
 *
 * \brief Get minimum of maximum distance of boxes
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   n_pts          Number of points
 * \param [in]   pts            Point coordinates
 * \param [out]  box_id         Leaf box with the minimum of maximum distance
 * \param [out]  box_max_dist   Maximum distance to box_id
 *
 */

void
PDM_box_tree_min_dist_max_box
(
 PDM_box_tree_t  *bt,
 const int        n_pts,
 double          *pts,
 int             *box_id,
 double          *box_max_dist
 )
{

  int s_pt_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);
  int pos_stack = 0;

  int dim = bt->boxes->dim;

  int normalized = bt->boxes->normalized;
  const double *d = bt->boxes->d;

  double *_pts = pts;
  /* if (normalized) { */
    _pts = malloc (sizeof(double) * 3 * n_pts);
    /*for (int i = 0; i < n_pts; i++) {
      const double *_pt_origin =  pts + 3 * i;
      double *_pt        = _pts + 3 * i;
      PDM_box_set_normalize ((PDM_box_set_t *) bt->boxes, _pt_origin, _pt);
      log_trace("pt = %f %f %f >> %f %f %f\n",
                _pt_origin[0],
                _pt_origin[1],
                _pt_origin[2],
                _pt[0],
                _pt[1],
                _pt[2]);
                }*/
    PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                  n_pts,
                                  pts,
                                  _pts);
    // for (int i = 0; i < n_pts; i++) {
    //   const double *_pt_origin =  pts + 3 * i;
    //   double *_pt              = _pts + 3 * i;
    //   log_trace("pt = %f %f %f >> %f %f %f\n",
    //             _pt_origin[0],
    //             _pt_origin[1],
    //             _pt_origin[2],
    //             _pt[0],
    //             _pt[1],
    //             _pt[2]);
    // }
  /* } */

  double extents2[2*dim];

  for (int i = 0; i < n_pts; i++) {

    const double *_pt = _pts + 3 * i;

    /* Init stack */

    box_id[i] = -1;
    box_max_dist[i] = HUGE_VAL;

    pos_stack = 0;
    stack[pos_stack] = 0; /* push root in th stack */

    _extents (dim, bt->local_data->nodes[0].morton_code, extents2);

    inbox_stack[pos_stack] = _box_dist2_min (dim,
                                             normalized,
                                             d,
                                             extents2,
                                             _pt,
                                             min_dist2_stack);

    pos_stack++;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];

      _node_t *curr_node = &(bt->local_data->nodes[id_curr_node]);

      if (curr_node->n_boxes == 0)
        continue;

      _extents (dim, curr_node->morton_code, extents2);

      double max_dist2;

      _box_dist2_max (dim,
                      normalized,
                      d,
                      extents2,
                      _pt,
                      &max_dist2);

      if (max_dist2 <= box_max_dist[i]) {

        if (!curr_node->is_leaf) {

          _push_child_in_stack_v0 (bt,
                                   dim,
                                   normalized,
                                   d,
                                   id_curr_node,
                                   box_max_dist[i],
                                   _pt,
                                   &pos_stack,
                                   stack,
                                   inbox_stack,
                                   min_dist2_stack,
                                   0,
                                   1);

        }

        else {

          for (int j = 0; j < curr_node->n_boxes; j++) {

            double box_max_dist2;

            int   _box_id = bt->local_data->box_ids[curr_node->start_id + j];
            const double *_box_extents =  bt->boxes->local_boxes->extents + _box_id*dim*2;

            _box_dist2_max (dim,
                            normalized,
                            d,
                            _box_extents,
                            _pt,
                            &box_max_dist2);

            if (box_max_dist2 < box_max_dist[i]) {
              box_id[i] = _box_id;
              box_max_dist[i] = box_max_dist2;
            }
          }

        }
      }
    }
  }

  free (stack);
  free (inbox_stack);
  free (min_dist2_stack);

  if (_pts != pts) {
    free (_pts);
  }

}



/**
 *
 * \brief Get an indexed list of all boxes within an upper bound distance
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_pts              Number of points
 * \param [in]   pts                Point coordinates
 * \param [in]   upper_bound_dist2  Squared upper bound distance (size = \ref n_pts)
 * \param [out]  pts_box_idx            Pointer to the index array on bounding boxes
 * \param [out]  boxes              Pointer to the list of bounding boxes
 *
 */

void
PDM_box_tree_closest_upper_bound_dist_boxes_get
(
 PDM_box_tree_t  *bt,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *pts_box_idx[],
 int             *boxes[]
 )
{
  //printf ("PDM_box_tree_closest_upper_bound_dist_boxes_get\n");

  int normalized = bt->boxes->normalized;
  const double *d = bt->boxes->d;

  double *_pts = pts;
  /* if (normalized) { */
  _pts = malloc (sizeof(double) * 3 * n_pts);
  /*for (int i = 0; i < n_pts; i++) {
    const double *_pt_origin = pts + 3 * i;
    double *_pt        = _pts + 3 * i;
    PDM_box_set_normalize ((PDM_box_set_t *)bt->boxes, _pt_origin, _pt);
    }*/
  PDM_box_set_normalize_robust ((PDM_box_set_t *)bt->boxes,
                                n_pts,
                                pts,
                                _pts);
  /* } */

  int s_pt_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);

  *pts_box_idx = malloc (sizeof(int) * (n_pts + 1));
  int *_pts_box_idx = *pts_box_idx;

  PDM_array_reset_int(_pts_box_idx, n_pts + 1, 0);

  int tmp_s_boxes = 4 * n_pts;
  *boxes = malloc (sizeof(int) * tmp_s_boxes);
  int *_boxes = *boxes;

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);

  int pos_stack = 0;

  int dim = bt->boxes->dim;

  int idx_box = 0;

  int n_boxes = bt->boxes->local_boxes->n_boxes;

  int *tag = PDM_array_zeros_int(n_boxes);

  int n_visited_boxes = 0;

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  size_t n_node = 0;
  size_t n_node_vid = 0;

  double extents2[2*dim];

  for (int i = 0; i < n_pts; i++) {
    const double *_pt = _pts + 3 * i;
    int flag = 0;

    /* Init stack :  push root */

    pos_stack = 0;
    stack[pos_stack] = 0; /* push root in th stack */
    _extents (dim, bt->local_data->nodes[0].morton_code, extents2);
    inbox_stack[pos_stack] = _box_dist2_min (dim,
                                             normalized,
                                             d,
                                             extents2,
                                             _pt,
                                             min_dist2_stack);


    pos_stack++;
    n_visited_boxes = 0;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];

      _node_t *curr_node = &(bt->local_data->nodes[id_curr_node]);

      if (curr_node->n_boxes == 0)
        continue;

      n_node++;
      if (curr_node->n_boxes == 0)
        n_node_vid++;

      double min_dist2 = min_dist2_stack[pos_stack];

      int inbox = inbox_stack[pos_stack];

      if ((min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
        if (!curr_node->is_leaf) {

          _push_child_in_stack_v0 (bt,
				   dim,
				   normalized,
				   d,
				   id_curr_node,
				   upper_bound_dist2[i],
				   _pt,
				   &pos_stack,
				   stack,
				   inbox_stack,
				   min_dist2_stack,
				   flag,
				   0);

        }

        else {

          for (int j = 0; j < curr_node->n_boxes; j++) {

            double box_min_dist2;

            int   _box_id = bt->local_data->box_ids[curr_node->start_id + j];

            if (tag[_box_id] == 0) {

              const double *_box_extents =  bt->boxes->local_boxes->extents + _box_id*dim*2;

              inbox = _box_dist2_min (dim,
                                      normalized,
                                      d,
                                      _box_extents,
                                      _pt,
                                      &box_min_dist2);

              if ((box_min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
                if (idx_box >= tmp_s_boxes) {
                  tmp_s_boxes += PDM_MAX(1, tmp_s_boxes/3);
                  *boxes = realloc (*boxes, sizeof(int) * tmp_s_boxes);
                  _boxes = *boxes;
                }
                _boxes[idx_box++] = _box_id;
                _pts_box_idx[i+1]++;
              }
              visited_boxes[n_visited_boxes++] = _box_id;
              tag[_box_id] = 1;
            }
          }
        }
      }
    }

    for (int j = 0; j < n_visited_boxes; j++) {
      tag[visited_boxes[j]] = 0;
    }

  }


  //printf ("[%d] Parours arbre : %ld \n", iappel, n_node);
  iappel+=1;

  for (int i = 0; i < n_pts; i++) {
    _pts_box_idx[i+1] += _pts_box_idx[i];
  }

  *boxes = realloc (*boxes, sizeof(int) * _pts_box_idx[n_pts]);
  if (pts != _pts) {
    free (_pts);
  }

  free (tag);
  free (stack);
  free (inbox_stack);
  free (min_dist2_stack);
  free (visited_boxes);
}



/**
 *
 * \brief Get an indexed list of all boxes within an upper bound distance
 *
 * The search can be performed either in the local box tree (\ref i_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_rank]
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   i_rank             Copied rank
 * \param [in]   n_pts              Number of points
 * \param [in]   pts                Point coordinates
 * \param [in]   upper_bound_dist2  Squared upper bound distance (size = \ref n_pts)
 * \param [out]  pts_box_idx            Pointer to the index array on bounding boxes
 * \param [out]  boxes              Pointer to the list of bounding boxes
 * \param [in]   d_opt              Normalization vector (or NULL)
 *
 */
static
void
_box_tree_closest_upper_bound_dist_boxes_impl
(
 PDM_box_tree_t       *bt,
 PDM_boxes_t          *boxes,
 PDM_box_tree_data_t  *box_tree_data,
 const int             n_pts,
 double                pts[],
 double                upper_bound_dist2[],
 int                  *pts_box_idx[],
 int                  *pts_box[],
 const double         *d_opt
)
{

  int normalized = bt->boxes->normalized;

  /*
     Fix inconsistency between dbbt->d and bt->d
     in the case where  non-NULL global extents
     were passed to PDM_dbbtree_create.
  */
  const double *d;
  if (d_opt == NULL) {
    d = bt->boxes->d;
  } else {
    d = d_opt;
  }

  double *_pts = pts;
  /* if (normalized) { */
  /*  _pts = malloc (sizeof(double) * 3 * n_pts);
    for (int i = 0; i < n_pts; i++) {
      const double *_pt_origin = pts + 3 * i;
      double *_pt              = _pts + 3 * i;
      PDM_box_set_normalize ((PDM_box_set_t *) bt->boxes, _pt_origin, _pt);
      }*/
  PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                n_pts,
                                pts,
                                _pts);
  /* } */

  int s_pt_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);

  *pts_box_idx = malloc (sizeof(int) * (n_pts + 1));
  int *_pts_box_idx = *pts_box_idx;

  PDM_array_reset_int(_pts_box_idx, n_pts + 1, 0);

  int tmp_s_boxes = 4 * n_pts;
  *pts_box = malloc (sizeof(int) * tmp_s_boxes);
  int *_pts_box = *pts_box;

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);

  int pos_stack = 0;

  int dim = bt->boxes->dim;

  int idx_box = 0;

  int n_boxes = boxes->n_boxes;

  int *tag = PDM_array_zeros_int(n_boxes);

  int n_visited_boxes = 0;
  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  double extents2[2*dim];
  _node_t *nodes         = box_tree_data->nodes;
  int     *box_ids       = box_tree_data->box_ids;
  double  *boxes_extents = boxes->extents;
  _extents (dim, nodes[0].morton_code, extents2);

  _node_t *curr_node = NULL;

  for (int i = 0; i < n_pts; i++) {
    const double *_pt = _pts + 3 * i;
    int flag = 0;

    pos_stack = 0;
    stack[pos_stack] = 0; /* push root in stack */
    inbox_stack[pos_stack] = _box_dist2_min (dim,
                                             normalized,
                                             d,
                                             extents2,
                                             _pt,
                                             min_dist2_stack);


    pos_stack++;
    n_visited_boxes = 0;

    while (pos_stack > 0) {

      int id_curr_node = stack[--pos_stack];
      curr_node = nodes + id_curr_node;

      if (curr_node->n_boxes == 0)
        continue;

      double min_dist2 = min_dist2_stack[pos_stack];

      int inbox = inbox_stack[pos_stack];

      if ((min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
        if (!curr_node->is_leaf) {
          _push_child_in_stack_v2 (bt,
                                   box_tree_data,
                                   dim,
                                   normalized,
                                   d,
                                   id_curr_node,
                                   upper_bound_dist2[i],
                                   _pt,
                                   &pos_stack,
                                   stack,
                                   inbox_stack,
                                   min_dist2_stack,
                                   flag,
                                   0);
        }

        else {

          for (int j = 0; j < curr_node->n_boxes; j++) {

            double box_min_dist2;

            int _box_id = box_ids[curr_node->start_id + j];

            if (tag[_box_id] == 0) {

              double *_box_extents = boxes_extents + _box_id*dim*2;

              inbox = _box_dist2_min (dim,
                                      normalized,
                                      d,
                                      _box_extents,
                                      _pt,
                                      &box_min_dist2);

              if ((box_min_dist2 <= upper_bound_dist2[i]) || (inbox == 1)) {
                if (idx_box >= tmp_s_boxes) {
                  tmp_s_boxes += PDM_MAX (1, tmp_s_boxes/3);
                  *pts_box = realloc (*pts_box, sizeof(int) * tmp_s_boxes);
                  _pts_box = *pts_box;
                }
                _pts_box[idx_box++] = _box_id;
                _pts_box_idx[i+1]++;
              }
              visited_boxes[n_visited_boxes++] = _box_id;
              tag[_box_id] = 1;
            }
          }
        }
      }
    }

    for (int j = 0; j < n_visited_boxes; j++) {
      tag[visited_boxes[j]] = 0;
    }

  }

  //printf ("[%d] Parours arbre : %ld \n", iappel, n_node);
  iappel+=1;

  for (int i = 0; i < n_pts; i++) {
    _pts_box_idx[i+1] += _pts_box_idx[i];
  }

  *pts_box = realloc (*pts_box, sizeof(int) * _pts_box_idx[n_pts]);
  if (pts != _pts) {
    free (_pts);
  }

  free (tag);
  free (stack);
  free (inbox_stack);
  free (min_dist2_stack);
  free (visited_boxes);
}



void
PDM_box_tree_closest_upper_bound_dist_boxes_get_v2
(
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *pts_box_idx[],
 int             *pts_box[],
 const double    *d_opt
)
{
  assert(i_copied_rank < bt->n_copied_ranks);
  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  _box_tree_closest_upper_bound_dist_boxes_impl(bt,
                                                boxes,
                                                box_tree_data,
                                                n_pts,
                                                pts,
                                                upper_bound_dist2,
                                                pts_box_idx,
                                                pts_box,
                                                d_opt);
}


void
PDM_box_tree_closest_upper_bound_dist_boxes_get_shared
(
 PDM_box_tree_t  *bt,
 const int        i_shm,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *pts_box_idx[],
 int             *pts_box[],
 const double    *d_opt
)
{
  PDM_boxes_t         *boxes         = &bt->boxes->shm_boxes[i_shm];
  PDM_box_tree_data_t *box_tree_data = &bt->shm_data        [i_shm];
  _box_tree_closest_upper_bound_dist_boxes_impl(bt,
                                                boxes,
                                                box_tree_data,
                                                n_pts,
                                                pts,
                                                upper_bound_dist2,
                                                pts_box_idx,
                                                pts_box,
                                                d_opt);
}


void
PDM_box_tree_closest_upper_bound_dist_boxes_get_shared_box_pov
(
 PDM_box_tree_t  *bt,
 const int        i_shm,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *box_pts_idx[],
 int             *box_pts[],
 const double    *d_opt
)
{

  PDM_boxes_t         *boxes         = &bt->boxes->shm_boxes[i_shm];
  PDM_box_tree_data_t *box_tree_data = &bt->shm_data        [i_shm];
  int *pts_box_idx = NULL;
  int *pts_box     = NULL;
  _box_tree_closest_upper_bound_dist_boxes_impl(bt,
                                                boxes,
                                                box_tree_data,
                                                n_pts,
                                                pts,
                                                upper_bound_dist2,
                                                &pts_box_idx,
                                                &pts_box,
                                                d_opt);

  /* Transpose from pts->box to box->line */
  int n_boxes = boxes->n_boxes;

  int* box_pts_n = PDM_array_zeros_int(n_boxes);
  for(int i = 0; i < n_pts; ++i) {
    for(int idx_box = pts_box_idx[i]; idx_box <  pts_box_idx[i+1]; ++idx_box ) {
      box_pts_n[pts_box[idx_box]]++;
    }
  }

  *box_pts_idx = PDM_array_new_idx_from_sizes_int(box_pts_n, n_boxes);
  *box_pts = (int *) malloc(sizeof(int) * (*box_pts_idx)[n_boxes]);
  PDM_array_reset_int(box_pts_n, n_boxes, 0);
  for (int ipts = 0; ipts < n_pts; ipts++) {
    for (int idx_box = pts_box_idx[ipts]; idx_box < pts_box_idx[ipts+1]; idx_box++) {
      int ibox = pts_box[idx_box];
      (*box_pts)[(*box_pts_idx)[ibox] + box_pts_n[ibox]++] = ipts;
    }
  }
  free(box_pts_n);
  free(pts_box_idx);
  free(pts_box);
}

void
PDM_box_tree_closest_upper_bound_dist_boxes_get_v2_box_pov
(
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *box_pts_idx[],
 int             *box_pts[],
 const double    *d_opt
)
{
  assert(i_copied_rank < bt->n_copied_ranks);
  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int *pts_box_idx = NULL;
  int *pts_box     = NULL;
  _box_tree_closest_upper_bound_dist_boxes_impl(bt,
                                                boxes,
                                                box_tree_data,
                                                n_pts,
                                                pts,
                                                upper_bound_dist2,
                                                &pts_box_idx,
                                                &pts_box,
                                                d_opt);
  /* Transpose from pts->box to box->line */
  int n_boxes = boxes->n_boxes;

  int* box_pts_n = PDM_array_zeros_int(n_boxes);
  for(int i = 0; i < n_pts; ++i) {
    for(int idx_box = pts_box_idx[i]; idx_box <  pts_box_idx[i+1]; ++idx_box ) {
      box_pts_n[pts_box[idx_box]]++;
    }
  }

  *box_pts_idx = PDM_array_new_idx_from_sizes_int(box_pts_n, n_boxes);
  *box_pts = (int *) malloc(sizeof(int) * (*box_pts_idx)[n_boxes]);
  PDM_array_reset_int(box_pts_n, n_boxes, 0);
  for (int ipts = 0; ipts < n_pts; ipts++) {
    for (int idx_box = pts_box_idx[ipts]; idx_box < pts_box_idx[ipts+1]; idx_box++) {
      int ibox = pts_box[idx_box];
      (*box_pts)[(*box_pts_idx)[ibox] + box_pts_n[ibox]++] = ipts;
    }
  }
  free(box_pts_n);
  free(pts_box_idx);
  free(pts_box);
}


/**
 *
 * \brief Copy the local box tree of some ranks on all other ranks
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_copied_ranks     Number of copied ranks
 * \param [in]   copied_ranks       List of copied ranks
 * \param [out]  rank_copy_num      Transpose of list of copied ranks (-1 for non-copied ranks)
 *
 */

void
PDM_box_tree_copy_to_ranks
(
 PDM_box_tree_t *bt,
 int            *n_copied_ranks,
 int            *copied_ranks,
 int            *rank_copy_num
)
{
  // copy boxes
  PDM_box_copy_boxes_to_ranks ((PDM_box_set_t *) bt->boxes, *n_copied_ranks, copied_ranks);

  // set (copy) tree structure
  int my_rank;
  PDM_MPI_Comm_rank(bt->comm, &my_rank);
  int n_ranks;
  PDM_MPI_Comm_size(bt->comm, &n_ranks);

  bt->copied_ranks = (int *) malloc (sizeof(int) * (*n_copied_ranks));
  if ( rank_copy_num == NULL ) {
    rank_copy_num = (int *) malloc (sizeof(int) * n_ranks);
  }
  int i = 0;
  PDM_array_reset_int(rank_copy_num, n_ranks, -1);


  bt->n_copied_ranks = 0;
  int i_rank = 0;

  for (i = 0; i < *n_copied_ranks; i++) {
    i_rank = copied_ranks[i];
    if ( my_rank != i_rank ) {
      rank_copy_num[copied_ranks[i]]         = bt->n_copied_ranks;
      bt->copied_ranks[bt->n_copied_ranks++] = copied_ranks[i];
    }
  }
  bt->copied_ranks = (int *) realloc (bt->copied_ranks, sizeof(int) * bt->n_copied_ranks);

  bt->rank_data = (PDM_box_tree_data_t *) malloc (sizeof(PDM_box_tree_data_t) * bt->n_copied_ranks);
  for (i = 0; i < bt->n_copied_ranks; i++) {
    bt->rank_data[i].stack     = NULL;
    bt->rank_data[i].pos_stack = NULL;
  }

  int n_max_nodes    = 0;
  int n_nodes        = 0;
  int n_build_loops  = 0;
  int l_box_ids      = 0;
  //int int_buffer[4] = {0};

  int *nodes_is_leaf  = NULL;
  int *nodes_morton   = NULL;
  int *nodes_n_boxes  = NULL;
  int *nodes_start_id = NULL;

  int *child_ids      = NULL;
  int *box_ids        = NULL;


  int icopied = 0;
  int j = 0, k = 0;
  for (i = 0; i < *n_copied_ranks; i++) {
    i_rank = copied_ranks[i];
    if ( my_rank == i_rank ) {
      n_max_nodes   = bt->local_data->n_max_nodes;
      n_nodes       = bt->local_data->n_nodes;
      n_build_loops = bt->local_data->n_build_loops;
      l_box_ids     = bt->stats.n_linked_boxes;

      /*
	int_buffer[0] = n_max_nodes;
	int_buffer[1] = n_nodes;
	int_buffer[2] = n_build_loops;
	int_buffer[3] = l_box_ids;
      */
    }

    PDM_MPI_Bcast(&n_max_nodes,   1, PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(&n_nodes,       1, PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(&n_build_loops, 1, PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(&l_box_ids,     1, PDM_MPI_INT, i_rank, bt->comm);
    /*
      PDM_MPI_Bcast(int_buffer, 4, PDM_MPI_INT, i_rank, bt->comm);
      n_max_nodes   = int_buffer[0];
      n_nodes       = int_buffer[1];
      n_build_loops = int_buffer[2];
      l_box_ids     = int_buffer[3];
    */

    // prepare buffers
    nodes_is_leaf  = (int *) malloc (sizeof(int) * n_max_nodes);
    nodes_n_boxes  = (int *) malloc (sizeof(int) * n_max_nodes);
    nodes_start_id = (int *) malloc (sizeof(int) * n_max_nodes);
    nodes_morton   = (int *) malloc (sizeof(int) * n_max_nodes*4);

    child_ids      = (int *) malloc (sizeof(int) * n_max_nodes*bt->n_children);
    box_ids        = (int *) malloc (sizeof(int) * l_box_ids);

    if ( my_rank == i_rank ) {
      // set buffers
      for (j = 0; j < n_max_nodes; j++) {
        nodes_is_leaf[j]  = (int) bt->local_data->nodes[j].is_leaf;
        nodes_n_boxes[j]  =       bt->local_data->nodes[j].n_boxes;
        nodes_start_id[j] =       bt->local_data->nodes[j].start_id;
        nodes_morton[4*j] = (int) bt->local_data->nodes[j].morton_code.L;
        for (k = 0; k < 3; k++) {
          nodes_morton[4*j+k+1] = (int) bt->local_data->nodes[j].morton_code.X[k];
        }
      }

      memcpy(child_ids, bt->local_data->child_ids, sizeof(int) * n_max_nodes*bt->n_children);
      memcpy(box_ids,   bt->local_data->box_ids,   sizeof(int) * l_box_ids);
    }

    // broadcast buffers
    PDM_MPI_Bcast(nodes_is_leaf,  n_max_nodes,   PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(nodes_n_boxes,  n_max_nodes,   PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(nodes_start_id, n_max_nodes,   PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(nodes_morton,   n_max_nodes*4, PDM_MPI_INT, i_rank, bt->comm);

    PDM_MPI_Bcast(child_ids, n_max_nodes*bt->n_children, PDM_MPI_INT, i_rank, bt->comm);
    PDM_MPI_Bcast(box_ids,   l_box_ids,                  PDM_MPI_INT, i_rank, bt->comm);

    if  ( my_rank != i_rank ) {
      bt->rank_data[icopied].n_max_nodes   = n_max_nodes;
      bt->rank_data[icopied].n_nodes       = n_nodes;
      bt->rank_data[icopied].n_build_loops = n_build_loops;

      bt->rank_data[icopied].nodes     = (_node_t *) malloc (sizeof(_node_t) * n_max_nodes);
      bt->rank_data[icopied].child_ids = (int *)     malloc (sizeof(int)     * n_max_nodes*bt->n_children);
      bt->rank_data[icopied].box_ids   = (int *)     malloc (sizeof(int)     * l_box_ids);

      for (j = 0; j < n_max_nodes; j++) {
        bt->rank_data[icopied].nodes[j].is_leaf       = (_Bool)            nodes_is_leaf[j];
        bt->rank_data[icopied].nodes[j].n_boxes       =                    nodes_n_boxes[j];
        bt->rank_data[icopied].nodes[j].start_id      =                    nodes_start_id[j];
        bt->rank_data[icopied].nodes[j].morton_code.L = (PDM_morton_int_t) nodes_morton[4*j];
        for (k = 0; k < 3; k++) {
          bt->rank_data[icopied].nodes[j].morton_code.X[k] = (PDM_morton_int_t) nodes_morton[4*j+k+1];
        }
      }

      memcpy(bt->rank_data[icopied].child_ids, child_ids, sizeof(int) * n_max_nodes*bt->n_children);
      memcpy(bt->rank_data[icopied].box_ids,   box_ids,   sizeof(int) * l_box_ids);

      icopied++;
    }

    free(nodes_is_leaf);
    free(nodes_n_boxes);
    free(nodes_start_id);
    free(nodes_morton);
    free(child_ids);
    free(box_ids);
  }

  *n_copied_ranks = bt->n_copied_ranks;
}


/**
 *
 * \brief Copy the local box tree of some ranks on all other ranks
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_copied_ranks     Number of copied ranks
 * \param [in]   copied_ranks       List of copied ranks
 * \param [out]  rank_copy_num      Transpose of list of copied ranks (-1 for non-copied ranks)
 *
 */

void
PDM_box_tree_copy_to_shm
(
 PDM_box_tree_t *bt
)
{
  // copy boxes
  PDM_box_copy_boxes_to_shm ((PDM_box_set_t *) bt->boxes);

  int n_rank, i_rank;
  PDM_MPI_Comm_rank(bt->comm, &i_rank);
  PDM_MPI_Comm_size(bt->comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(bt->comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  bt->n_rank_in_shm = n_rank_in_shm;
  // log_trace("PDM_box_tree_copy_to_shm ; n_rank_in_shm = %i \n", n_rank_in_shm);


  int s_shm_data_in_rank[4] = {0};
  s_shm_data_in_rank[0] = bt->local_data->n_max_nodes;
  s_shm_data_in_rank[1] = bt->local_data->n_nodes;
  s_shm_data_in_rank[2] = bt->local_data->n_build_loops;
  s_shm_data_in_rank[3] = bt->stats.n_linked_boxes;
  int *s_shm_data_in_all_nodes = malloc(4 * n_rank_in_shm * sizeof(int));

  PDM_MPI_Allgather(s_shm_data_in_rank     , 4, PDM_MPI_INT,
                    s_shm_data_in_all_nodes, 4, PDM_MPI_INT, comm_shared);

  bt->shm_data       = (PDM_box_tree_data_t *) malloc(n_rank_in_shm * sizeof(PDM_box_tree_data_t));
  bt->wbox_tree_data = (_w_box_tree_data_t  *) malloc(n_rank_in_shm * sizeof(_w_box_tree_data_t ));

  /* Creation mémoire des windows */
  for(int i = 0; i < n_rank_in_shm; ++i) {
    int n_max_nodes    = s_shm_data_in_all_nodes[4*i  ];
    int n_nodes        = s_shm_data_in_all_nodes[4*i+1];
    int n_build_loops  = s_shm_data_in_all_nodes[4*i+2];
    int n_linked_boxes = s_shm_data_in_all_nodes[4*i+3];
    bt->wbox_tree_data[i].n_max_nodes    = n_max_nodes;
    bt->wbox_tree_data[i].n_nodes        = n_nodes;
    bt->wbox_tree_data[i].n_linked_boxes = n_linked_boxes;
    bt->wbox_tree_data[i].w_nodes        = PDM_mpi_win_shared_create(n_max_nodes               , sizeof(_node_t), comm_shared);
    bt->wbox_tree_data[i].w_child_ids    = PDM_mpi_win_shared_create(n_max_nodes*bt->n_children, sizeof(int    ), comm_shared);
    bt->wbox_tree_data[i].w_box_ids      = PDM_mpi_win_shared_create(n_linked_boxes            , sizeof(int    ), comm_shared);

    PDM_mpi_win_shared_lock_all (0, bt->wbox_tree_data[i].w_nodes    );
    PDM_mpi_win_shared_lock_all (0, bt->wbox_tree_data[i].w_child_ids);
    PDM_mpi_win_shared_lock_all (0, bt->wbox_tree_data[i].w_box_ids  );

    bt->shm_data[i].n_max_nodes   = n_max_nodes;
    bt->shm_data[i].n_nodes       = n_nodes;
    bt->shm_data[i].n_build_loops = n_build_loops;

    // Setup alias
    bt->shm_data[i].nodes     = PDM_mpi_win_shared_get (bt->wbox_tree_data[i].w_nodes    );
    bt->shm_data[i].child_ids = PDM_mpi_win_shared_get (bt->wbox_tree_data[i].w_child_ids);
    bt->shm_data[i].box_ids   = PDM_mpi_win_shared_get (bt->wbox_tree_data[i].w_box_ids  );

  }
  PDM_MPI_Barrier (comm_shared);

  /* Copy from local to shared (After windows creation bcause window call is collective ) */
  for (int j = 0; j < bt->local_data->n_max_nodes; j++) {
    bt->shm_data[i_rank_in_shm].nodes[j].is_leaf          = bt->local_data->nodes[j].is_leaf;
    bt->shm_data[i_rank_in_shm].nodes[j].n_boxes          = bt->local_data->nodes[j].n_boxes;
    bt->shm_data[i_rank_in_shm].nodes[j].start_id         = bt->local_data->nodes[j].start_id;
    bt->shm_data[i_rank_in_shm].nodes[j].morton_code.L    = bt->local_data->nodes[j].morton_code.L;
    bt->shm_data[i_rank_in_shm].nodes[j].morton_code.X[0] = bt->local_data->nodes[j].morton_code.X[0];
    bt->shm_data[i_rank_in_shm].nodes[j].morton_code.X[1] = bt->local_data->nodes[j].morton_code.X[1];
    bt->shm_data[i_rank_in_shm].nodes[j].morton_code.X[2] = bt->local_data->nodes[j].morton_code.X[2];
  }

  memcpy(bt->shm_data[i_rank_in_shm].child_ids, bt->local_data->child_ids, sizeof(int) * bt->local_data->n_max_nodes*bt->n_children);
  memcpy(bt->shm_data[i_rank_in_shm].box_ids,   bt->local_data->box_ids  , sizeof(int) * bt->stats.n_linked_boxes);

  free(s_shm_data_in_all_nodes);

  PDM_MPI_Barrier (comm_shared);
  // PDM_MPI_Comm_free(&comm_shared);
}

void
PDM_box_tree_free_copies
(
 PDM_box_tree_t *bt
 )
 {
  PDM_box_set_free_copies(&bt->boxes);

  if (bt->copied_ranks != NULL) {
    free (bt->copied_ranks);
    bt->copied_ranks = NULL;
  }

  if (bt->rank_data != NULL) {
    for (int i = 0; i < bt->n_copied_ranks; i++) {
      PDM_box_tree_data_destroy(&(bt->rank_data[i]));
    }
    free(bt->rank_data);
    bt->rank_data = NULL;
  }

  bt->n_copied_ranks = 0;
 }


/**
 *
 * \brief Get an indexed list of all points inside the boxes of a box tree
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_pts              Number of points
 * \param [in]   pts_g_num          Point global ids
 * \param [in]   pts_coord          Point coordinates
 * \param [out]  pts_in_box_idx     Pointer to the index array on boxes
 * \param [out]  pts_in_box_g_num   Pointer to the list of global ids of points inside boxes
 * \param [out]  pts_in_box_coord   Pointer to the list of coordinates of points inside boxes
 *
 */

void
PDM_box_tree_points_inside_boxes
(
 PDM_box_tree_t     *bt,
 const int           n_pts,
 const PDM_g_num_t   pts_g_num[],
 const double        pts_coord[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 )
{
  const int dim = bt->boxes->dim;
  int n_boxes = bt->boxes->local_boxes->n_boxes;

  double *_pts_coord = NULL;
  if (bt->boxes->normalized) {
    _pts_coord = malloc (sizeof(double) * 3 * n_pts);
    /*for (int i = 0; i < n_pts; i++) {
      const double *_pt_origin = pts_coord + 3 * i;
      double *_pt = _pts_coord + 3 * i;
      PDM_box_set_normalize ((PDM_box_set_t *)bt->boxes, _pt_origin, _pt);
      }*/
    PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                  n_pts,
                                  (double *) pts_coord,
                                  _pts_coord);
  } else {
    _pts_coord = (double *) pts_coord;
  }


  int *pts_in_box_count = PDM_array_zeros_int(n_boxes);

  int *boxes_idx = malloc (sizeof(int) * (n_pts + 1));
  boxes_idx[0] = 0;

  int tmp_s_boxes = 4 * n_pts;
  int *boxes_l_num = malloc (sizeof(int) * tmp_s_boxes);


  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  double node_extents[2*dim];


  /* Loop over points */
  for (int ipt = 0; ipt < n_pts; ipt++) {

    boxes_idx[ipt+1] = boxes_idx[ipt];
    const double *_pt = _pts_coord + dim*ipt;

    /* Init stack */
    pos_stack = 0;
    _extents (dim,
              bt->local_data->nodes[0].morton_code,
              node_extents);

    if (_point_inside_box (dim, node_extents, _pt)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {
      int node_id = stack[--pos_stack];

      _node_t *node = &(bt->local_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = bt->local_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            const double *box_extents = bt->boxes->local_boxes->extents + box_id*2*dim;

            if (_point_inside_box (dim, box_extents, _pt)) {
              if (boxes_idx[ipt+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                boxes_l_num = realloc (boxes_l_num, sizeof(int) * tmp_s_boxes);
              }
              boxes_l_num[boxes_idx[ipt+1]++] = box_id;
              pts_in_box_count[box_id]++;
            }
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          }

        }
      }

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = bt->local_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          _extents (dim, bt->local_data->nodes[child_id].morton_code, node_extents);

          if (_point_inside_box (dim, node_extents, _pt)) {
            stack[pos_stack++] = child_id;
          }
        }
      }

    } // While stack not empty

    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // Loop over points

  if (pts_coord != _pts_coord) {
    free (_pts_coord);
  }

  free (is_visited_box);
  free (stack);
  free (visited_boxes);

  /* From {point -> boxes} to {box -> point} */

  *pts_in_box_idx = malloc (sizeof(int) * (n_boxes + 1));
  (*pts_in_box_idx)[0] = 0;
  for (int ibox = 0; ibox < n_boxes; ibox++) {
    (*pts_in_box_idx)[ibox+1] = (*pts_in_box_idx)[ibox] + pts_in_box_count[ibox];
    pts_in_box_count[ibox] = 0;
  }

  int size = (*pts_in_box_idx)[n_boxes];
  *pts_in_box_g_num = malloc (sizeof(PDM_g_num_t) * size);
  *pts_in_box_coord = malloc (sizeof(double)      * size * 3);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    for (int j = boxes_idx[ipt]; j < boxes_idx[ipt+1]; j++) {
      int ibox = boxes_l_num[j];
      int idx = (*pts_in_box_idx)[ibox] + pts_in_box_count[ibox];

      (*pts_in_box_g_num)[idx] = pts_g_num[ipt];
      for (int idim = 0; idim < 3; idim++) {
        (*pts_in_box_coord)[dim*idx + idim] = pts_coord[dim*ipt + idim];
      }

      pts_in_box_count[ibox]++;
    }
  }
  free (pts_in_box_count);
  free (boxes_idx);
  free (boxes_l_num);
}


/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_pts          Number of points
 * \param [in]   pts_coord      Point coordinates
 * \param [out]  box_idx        Pointer to the index array on points (size = \ref n_pts + 1)
 * \param [out]  box_l_num      Pointer to the list boxes containing points (size = \ref box_idx[\ref n_pts])
 *
 */
static
void
_box_tree_boxes_containing_points_impl
(
 PDM_box_tree_t       *bt,
 PDM_boxes_t          *boxes,
 PDM_box_tree_data_t  *box_tree_data,
 const int             n_pts,
 const double         *pts_coord,
 int                 **box_idx,
 int                 **box_l_num
)
{
  const int dim = bt->boxes->dim;
  //int normalized = bt->boxes->normalized;

  /* if (normalized) { */
  double *_pts_coord = malloc (sizeof(double) * dim * n_pts);
  PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                n_pts,
                                (double *) pts_coord,
                                _pts_coord);

  int n_boxes = boxes->n_boxes;

  *box_idx = malloc (sizeof(int) * (n_pts + 1));
  int *_box_idx = *box_idx;
  _box_idx[0] = 0;

  int tmp_s_boxes = 4 * n_pts;
  *box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  int *_box_l_num = *box_l_num;

  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  double node_extents[2*dim];

  for (int ipt = 0; ipt < n_pts; ipt++) {

    _box_idx[ipt+1] = _box_idx[ipt];
    const double *_pt = _pts_coord + dim * ipt;

    /* Init stack */
    pos_stack = 0;
    _extents (dim,
              box_tree_data->nodes[0].morton_code,
              node_extents);

    if (_point_inside_box (dim, node_extents, _pt)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      _node_t *node = &(box_tree_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = box_tree_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            const double *box_extents = boxes->extents + box_id*2*dim;

            if (_point_inside_box (dim, box_extents, _pt)) {
              if (_box_idx[ipt+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                *box_l_num = realloc (*box_l_num, sizeof(int) * tmp_s_boxes);
                _box_l_num = *box_l_num;
              }
              _box_l_num[_box_idx[ipt+1]++] = box_id;
            }
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          }

        }
      }

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          _extents (dim, box_tree_data->nodes[child_id].morton_code, node_extents);

          if (_point_inside_box (dim, node_extents, _pt)) {
            stack[pos_stack++] = child_id;
          }
        }
      }

    } // End while stack not empty

    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // End of loop on points

  if (pts_coord != _pts_coord) {
    free (_pts_coord);
  }

  free (is_visited_box);
  free (stack);
  free (visited_boxes);

  *box_l_num = realloc (*box_l_num, sizeof(int) * _box_idx[n_pts]);
}

/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_pts          Number of points
 * \param [in]   pts_coord      Point coordinates
 * \param [out]  box_idx        Pointer to the index array on points (size = \ref n_pts + 1)
 * \param [out]  box_l_num      Pointer to the list boxes containing points (size = \ref box_idx[\ref n_pts])
 *
 */

void
PDM_box_tree_boxes_containing_points
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_pts,
 const double   *pts_coord,
 int           **box_idx,
 int           **box_l_num
)
{
  assert(i_copied_rank < bt->n_copied_ranks);

  // const int dim = bt->boxes->dim;
  //int normalized = bt->boxes->normalized;

  // Done inside _box_tree_boxes_containing_points_impl
  // /* if (normalized) { */
  // double *_pts_coord = malloc (sizeof(double) * dim * n_pts);
  // /*for (int i = 0; i < n_pts; i++) {
  //   const double *_pt_origin = pts_coord + dim * i;
  //   double *_pt              = _pts_coord + dim * i;
  //   PDM_box_set_normalize ((PDM_box_set_t *) bt->boxes, _pt_origin, _pt);
  //   }*/
  // PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
  //                               n_pts,
  //                               (double *) pts_coord,
  //                               _pts_coord);
  // /* } */

  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  _box_tree_boxes_containing_points_impl(bt,
                                         boxes,
                                         box_tree_data,
                                         n_pts,
                                         pts_coord,
                                         box_idx,
                                         box_l_num);

  // int n_boxes = boxes->n_boxes;

  // *box_idx = malloc (sizeof(int) * (n_pts + 1));
  // int *_box_idx = *box_idx;
  // _box_idx[0] = 0;

  // int tmp_s_boxes = 4 * n_pts;
  // *box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  // int *_box_l_num = *box_l_num;

  // int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  // int *stack = malloc ((sizeof(int)) * s_stack);
  // int pos_stack = 0;

  // int n_visited_boxes = 0;
  // PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  // for (int i = 0; i < n_boxes; i++) {
  //   is_visited_box[i] = PDM_FALSE;
  // }

  // int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  // double node_extents[2*dim];

  // for (int ipt = 0; ipt < n_pts; ipt++) {

  //   _box_idx[ipt+1] = _box_idx[ipt];
  //   const double *_pt = _pts_coord + dim * ipt;

  //   /* Init stack */
  //   pos_stack = 0;
  //   _extents (dim,
  //             box_tree_data->nodes[0].morton_code,
  //             node_extents);

  //   if (_point_inside_box (dim, node_extents, _pt)) {
  //     stack[pos_stack++] = 0;
  //   }

  //   n_visited_boxes = 0;

  //   /* Traverse box tree */
  //   while (pos_stack > 0) {

  //     int node_id = stack[--pos_stack];

  //     _node_t *node = &(box_tree_data->nodes[node_id]);

  //     if (node->n_boxes == 0) {
  //       continue;
  //     }

  //     /* Leaf node */
  //     if (node->is_leaf) {
  //       /* inspect boxes contained in current leaf node */
  //       for (int ibox = 0; ibox < node->n_boxes; ibox++) {
  //         int box_id = box_tree_data->box_ids[node->start_id + ibox];

  //         if (is_visited_box[box_id] == PDM_FALSE) {
  //           const double *box_extents = boxes->extents + box_id*2*dim;

  //           if (_point_inside_box (dim, box_extents, _pt)) {
  //             if (_box_idx[ipt+1] >= tmp_s_boxes) {
  //               tmp_s_boxes *= 2;
  //               *box_l_num = realloc (*box_l_num, sizeof(int) * tmp_s_boxes);
  //               _box_l_num = *box_l_num;
  //             }
  //             _box_l_num[_box_idx[ipt+1]++] = box_id;
  //           }
  //           visited_boxes[n_visited_boxes++] = box_id;
  //           is_visited_box[box_id] = PDM_TRUE;
  //         }

  //       }
  //     }

  //     /* Internal node */
  //     else {
  //       /* inspect children of current node */
  //       int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
  //       for (int ichild = 0; ichild < bt->n_children; ichild++) {
  //         int child_id = child_ids[ichild];
  //         _extents (dim, box_tree_data->nodes[child_id].morton_code, node_extents);

  //         if (_point_inside_box (dim, node_extents, _pt)) {
  //           stack[pos_stack++] = child_id;
  //         }
  //       }
  //     }

  //   } // End while stack not empty

  //   for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
  //     is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
  //   }

  // } // End of loop on points

  // if (pts_coord != _pts_coord) {
  //   free (_pts_coord);
  // }

  // free (is_visited_box);
  // free (stack);
  // free (visited_boxes);

  // *box_l_num = realloc (*box_l_num, sizeof(int) * _box_idx[n_pts]);
}


/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_pts          Number of points
 * \param [in]   pts_coord      Point coordinates
 * \param [out]  box_idx        Pointer to the index array on points (size = \ref n_pts + 1)
 * \param [out]  box_l_num      Pointer to the list boxes containing points (size = \ref box_idx[\ref n_pts])
 *
 */

void
PDM_box_tree_boxes_containing_points_shared
(
 PDM_box_tree_t *bt,
 const int       i_shm,
 const int       n_pts,
 const double   *pts_coord,
 int           **box_idx,
 int           **box_l_num
)
{
  PDM_boxes_t         *boxes         = &bt->boxes->shm_boxes[i_shm];
  PDM_box_tree_data_t *box_tree_data = &bt->shm_data        [i_shm];
  _box_tree_boxes_containing_points_impl(bt,
                                         boxes,
                                         box_tree_data,
                                         n_pts,
                                         pts_coord,
                                         box_idx,
                                         box_l_num);
}

/**
 *
 * \brief Get an indexed list of all ellipsoids containing points
 *
 * The boxes are treated as axis-aligned ellipsoids with same extents.
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_pts          Number of points
 * \param [in]   pts_coord      Point coordinates
 * \param [out]  box_idx        Pointer to the index array on points (size = \ref n_pts + 1)
 * \param [out]  box_l_num      Pointer to the list of ellipsoids containing points (size = \ref box_idx[\ref n_pts])
 *
 */

void
PDM_box_tree_ellipsoids_containing_points
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_pts,
 const double   *pts_coord,
 int           **box_idx,
 int           **box_l_num
 )
{
  assert(i_copied_rank < bt->n_copied_ranks);

  const int dim = bt->boxes->dim;
  //int normalized = bt->boxes->normalized;

  /* if (normalized) { */
  double *_pts_coord = malloc (sizeof(double) * dim * n_pts);
  /*for (int i = 0; i < n_pts; i++) {
    const double *_pt_origin = pts_coord + dim * i;
    double *_pt              = _pts_coord + dim * i;
    PDM_box_set_normalize ((PDM_box_set_t *) bt->boxes, _pt_origin, _pt);
    }*/
  PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                n_pts,
                                (double *) pts_coord,
                                _pts_coord);
  /* } */

  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int n_boxes = boxes->n_boxes;


  *box_idx = malloc (sizeof(int) * (n_pts + 1));
  int *_box_idx = *box_idx;
  _box_idx[0] = 0;

  int tmp_s_boxes = 4 * n_pts;
  *box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  int *_box_l_num = *box_l_num;

  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  double node_extents[2*dim];

  for (int ipt = 0; ipt < n_pts; ipt++) {

    _box_idx[ipt+1] = _box_idx[ipt];
    const double *_pt = _pts_coord + dim * ipt;

    /* Init stack */
    pos_stack = 0;
    _extents (dim,
              box_tree_data->nodes[0].morton_code,
              node_extents);

    if (_point_inside_box (dim, node_extents, _pt)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      _node_t *node = &(box_tree_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = box_tree_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            const double *box_extents = boxes->extents + box_id*2*dim;

            if (_point_inside_ellipsoid (dim, box_extents, _pt)) {
              if (_box_idx[ipt+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                *box_l_num = realloc (*box_l_num, sizeof(int) * tmp_s_boxes);
                _box_l_num = *box_l_num;
              }
              _box_l_num[_box_idx[ipt+1]++] = box_id;
            }
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          }

        }
      }

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          _extents (dim, box_tree_data->nodes[child_id].morton_code, node_extents);

          if (_point_inside_box (dim, node_extents, _pt)) {
            stack[pos_stack++] = child_id;
          }
        }
      }

    } // End while stack not empty

    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // End of loop on points

  if (pts_coord != _pts_coord) {
    free (_pts_coord);
  }

  free (is_visited_box);
  free (stack);
  free (visited_boxes);

  *box_l_num = realloc (*box_l_num, sizeof(int) * _box_idx[n_pts]);
}

static
void
_box_tree_intersect_lines_boxes_impl
(
 PDM_box_tree_t       *bt,
 PDM_boxes_t          *boxes,
 PDM_box_tree_data_t  *box_tree_data,
 const int             n_line,
 const double         *line_coord,
 int                 **box_idx,
 int                 **box_l_num
)
{
  // double t1 = PDM_MPI_Wtime();
  const int dim = bt->boxes->dim;
  const int two_dim = 2*dim;
  //int normalized = bt->boxes->normalized;

  double *_line_coord = malloc(two_dim * n_line * sizeof(double));
  PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                n_line*2,
                     (double *) line_coord,
                                _line_coord);

  if(0 == 1){
    int i_rank;
    PDM_MPI_Comm_rank(bt->comm, &i_rank);
    char filename[999];
    sprintf(filename, "box_tree_intersect_lines_boxes_impl_line_coord_%i.vtk", i_rank);
    PDM_vtk_write_lines(filename,
                        n_line,
                        _line_coord,
                        NULL,
                        NULL);
  }


  int n_boxes = boxes->n_boxes;

  *box_idx = malloc (sizeof(int) * (n_line + 1));
  int *_box_idx = *box_idx;
  _box_idx[0] = 0;

  int tmp_s_boxes = 4 * n_line;
  *box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  int *_box_l_num = *box_l_num;

  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  // double node_extents[2*dim];

  /*
   * Compute once extents
   */
  double *node_extents = malloc(2 * dim * box_tree_data->n_nodes * sizeof(double));
  for(int i = 0; i < box_tree_data->n_nodes; ++i) {
    _extents (dim,
              box_tree_data->nodes[i].morton_code,
              &node_extents[i*2*dim]);
  }

  if(1 == 0){
    char filename[999];
    int i_rank;
    PDM_MPI_Comm_rank(bt->comm, &i_rank);
    sprintf(filename, "box_tree_intersect_lines_boxes_impl_node_extents_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        box_tree_data->n_nodes,
                        node_extents,
                        NULL);
  }

  double invdir[3];
  for (int iline = 0; iline < n_line; iline++) {

    _box_idx[iline+1] = _box_idx[iline];
    const double *origin = _line_coord + two_dim * iline;
    const double *destination = origin + dim;
    for (int i = 0; i < 3; i++) {
      double d = destination[i] - origin[i];
      if (PDM_ABS(d) < 1e-15) {
         invdir[i] = copysign(1, d) * HUGE_VAL;
      } else {
        invdir[i] = 1. / d;
      }
    }

    /* Init stack */
    pos_stack = 0;
    // _extents (dim, box_tree_data->nodes[0].morton_code, node_extents);
    // if (_intersect_line_box (dim, node_extents, origin, invdir)) {
    //   stack[pos_stack++] = 0;
    // }

    double *_node_extents = node_extents + 2 * dim * 0;
    if (_intersect_line_box (dim, _node_extents, origin, invdir)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      _node_t *node = &(box_tree_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = box_tree_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            const double *box_extents = boxes->extents + box_id*2*dim;

            if (_intersect_line_box (dim, box_extents, origin, invdir)) {
              if (_box_idx[iline+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                *box_l_num = realloc (*box_l_num, sizeof(int) * tmp_s_boxes);
                _box_l_num = *box_l_num;
              }
              _box_l_num[_box_idx[iline+1]++] = box_id;
            }
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          }

        }
      }

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          // _extents (dim, box_tree_data->nodes[child_id].morton_code, node_extents);
          // if (_intersect_line_box (dim, node_extents, origin, invdir)) {
          //   stack[pos_stack++] = child_id;
          // }

          double *_node_extents_child = node_extents + 2 * dim * child_id;
          if (_intersect_line_box (dim, _node_extents_child, origin, invdir)) {
            stack[pos_stack++] = child_id;
          }

        }
      }

    } // End while stack not empty

    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // End of loop on points

  if (line_coord != _line_coord) {
    free (_line_coord);
  }

  free (is_visited_box);
  free (stack);
  free (visited_boxes);

  *box_l_num = realloc (*box_l_num, sizeof(int) * _box_idx[n_line]);
  free(node_extents);

  // double t2 = PDM_MPI_Wtime();
  // log_trace("_box_tree_intersect_lines_boxes_impl = %12.5e \n", t2-t1);
}

/**
 *
 * \brief Get an indexed list of all boxes intersecting lines
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_line         Number of lines
 * \param [in]   line_coord     Lines coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  box_idx        Pointer to the index array on lines (size = \ref n_line + 1)
 * \param [out]  box_l_num      Pointer to the list of boxes intersecting lines (size = \ref box_idx[\ref n_line])
 *
 */

void
PDM_box_tree_intersect_lines_boxes
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_line,
 const double   *line_coord,
 int           **box_idx,
 int           **box_l_num
)
{
  assert(i_copied_rank < bt->n_copied_ranks);
  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  _box_tree_intersect_lines_boxes_impl(bt,
                                       boxes,
                                       box_tree_data,
                                       n_line,
                                       line_coord,
                                       box_idx,
                                       box_l_num);

}


void
PDM_box_tree_intersect_lines_boxes_shared
(
 PDM_box_tree_t *bt,
 const int       i_shm,
 const int       n_line,
 const double   *line_coord,
 int           **box_idx,
 int           **box_l_num
)
{
  PDM_boxes_t         *boxes         = &bt->boxes->shm_boxes[i_shm];
  PDM_box_tree_data_t *box_tree_data = &bt->shm_data        [i_shm];
  _box_tree_intersect_lines_boxes_impl(bt,
                                       boxes,
                                       box_tree_data,
                                       n_line,
                                       line_coord,
                                       box_idx,
                                       box_l_num);
}


/**
 *
 * \brief Get an indexed list of all lines intersecting boxes
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_line         Number of lines
 * \param [in]   line_coord     Lines coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  box_line_idx   Pointer to the index array on boxes (size = \ref n_box + 1)
 * \param [out]  box_line       Pointer to the list of lines intersecting boxes (size = \ref box_line_idx[\ref n_box])
 *
 */

void
PDM_box_tree_intersect_boxes_lines
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_line,
 const double   *line_coord,
 int           **box_line_idx,
 int           **box_line
)
{
  int *line_box_idx = NULL;
  int *line_box     = NULL;
  PDM_box_tree_intersect_lines_boxes(bt,
                                     i_copied_rank,
                                     n_line,
                                     line_coord,
                                     &line_box_idx,
                                     &line_box);

  /* Transpose from line->box to box->line */
  PDM_boxes_t *boxes;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
  }
  int n_boxes = boxes->n_boxes;

  int* box_line_n = PDM_array_zeros_int(n_boxes);
  for(int i = 0; i < n_line; ++i) {
    for(int idx_box = line_box_idx[i]; idx_box <  line_box_idx[i+1]; ++idx_box ) {
      box_line_n[line_box[idx_box]]++;
    }
  }

  *box_line_idx = PDM_array_new_idx_from_sizes_int(box_line_n, n_boxes);
  *box_line = (int *) malloc(sizeof(int) * (*box_line_idx)[n_boxes]);
  PDM_array_reset_int(box_line_n, n_boxes, 0);
  for (int iline = 0; iline < n_line; iline++) {
    for (int idx_box = line_box_idx[iline]; idx_box < line_box_idx[iline+1]; idx_box++) {
      int ibox = line_box[idx_box];
      (*box_line)[(*box_line_idx)[ibox] + box_line_n[ibox]++] = iline;
    }
  }
  free(box_line_n);
  free(line_box_idx);
  free(line_box);
}



void
PDM_box_tree_intersect_boxes_lines_shared
(
 PDM_box_tree_t *bt,
 const int       i_shm,
 const int       n_line,
 const double   *line_coord,
 int           **box_line_idx,
 int           **box_line
)
{
  int *line_box_idx = NULL;
  int *line_box     = NULL;
  PDM_box_tree_intersect_lines_boxes_shared(bt,
                                            i_shm,
                                            n_line,
                                            line_coord,
                                            &line_box_idx,
                                            &line_box);

  /* Transpose from line->box to box->line */
  PDM_boxes_t *boxes = &bt->boxes->shm_boxes[i_shm];
  int n_boxes = boxes->n_boxes;

  int* box_line_n = PDM_array_zeros_int(n_boxes);
  for(int i = 0; i < n_line; ++i) {
    for(int idx_box = line_box_idx[i]; idx_box <  line_box_idx[i+1]; ++idx_box ) {
      box_line_n[line_box[idx_box]]++;
    }
  }

  *box_line_idx = PDM_array_new_idx_from_sizes_int(box_line_n, n_boxes);
  *box_line = (int *) malloc(sizeof(int) * (*box_line_idx)[n_boxes]);
  PDM_array_reset_int(box_line_n, n_boxes, 0);
  for (int iline = 0; iline < n_line; iline++) {
    for (int idx_box = line_box_idx[iline]; idx_box < line_box_idx[iline+1]; idx_box++) {
      int ibox = line_box[idx_box];
      (*box_line)[(*box_line_idx)[ibox] + box_line_n[ibox]++] = iline;
    }
  }
  free(box_line_n);
  free(line_box_idx);
  free(line_box);
}
/**
 *
 * \brief Get an indexed list of all boxes intersecting lines
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_line         Number of lines
 * \param [in]   line_coord     Lines coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  box_idx        Pointer to the index array on lines (size = \ref n_line + 1)
 * \param [out]  box_l_num      Pointer to the list of boxes intersecting lines (size = \ref box_idx[\ref n_line])
 *
 */

void
PDM_box_tree_intersect_boxes_boxes
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_tgt_box,
 const double   *tgt_box_extents,
 int           **tgt_box_idx,
 int           **tgt_box_l_num
)
{

  assert(i_copied_rank < bt->n_copied_ranks);

  const int dim = bt->boxes->dim;
  const int two_dim = 2*dim;
  // int normalized = bt->boxes->normalized;

  // if (normalized) {
  double *_tgt_box_extents = malloc (sizeof(double) * two_dim * n_tgt_box);
  PDM_box_set_normalize_robust ((PDM_box_set_t *) bt->boxes,
                                n_tgt_box*2,
                                (double *) tgt_box_extents,
                                _tgt_box_extents);

  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int n_boxes = boxes->n_boxes;

  *tgt_box_idx = malloc (sizeof(int) * (n_tgt_box + 1));
  int *_tgt_box_idx = *tgt_box_idx;
  _tgt_box_idx[0] = 0;

  int tmp_s_boxes = 4 * n_tgt_box;
  *tgt_box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  int *_tgt_box_l_num = *tgt_box_l_num;

  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes); // A optimiser

  double node_extents[2*dim];

  for (int itgt_box = 0; itgt_box < n_tgt_box; itgt_box++) {

    _tgt_box_idx[itgt_box+1] = _tgt_box_idx[itgt_box];
    const double *lextents   = _tgt_box_extents + two_dim * itgt_box;

    /* Init stack */
    pos_stack = 0;
    _extents (dim,
              box_tree_data->nodes[0].morton_code,
              node_extents);

    // printf(" intersect with roots \n");
    if (_intersect_box_box (dim, node_extents, lextents)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      _node_t *node = &(box_tree_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = box_tree_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            const double *box_extents = boxes->extents + box_id*2*dim;

            if (_intersect_box_box (dim, box_extents, lextents)) {
              if (_tgt_box_idx[itgt_box+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                *tgt_box_l_num = realloc (*tgt_box_l_num, sizeof(int) * tmp_s_boxes);
                _tgt_box_l_num = *tgt_box_l_num;
              }
              _tgt_box_l_num[_tgt_box_idx[itgt_box+1]++] = box_id;
            }
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          }

        }
      }

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          _extents (dim, box_tree_data->nodes[child_id].morton_code, node_extents);

          // printf(" intersect with internal child_id = %i \n", child_id);
          if (_intersect_box_box (dim, node_extents, lextents)) {
            stack[pos_stack++] = child_id;
          }
        }
      }

    } // End while stack not empty

    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // End of loop on points

  if (tgt_box_extents != _tgt_box_extents) {
    free (_tgt_box_extents);
  }

  free (is_visited_box);
  free (stack);
  free (visited_boxes);

  *tgt_box_l_num = realloc (*tgt_box_l_num, sizeof(int) * _tgt_box_idx[n_tgt_box]);

}


void
PDM_box_tree_intersect_boxes_boxes2
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_tgt_box,
 const double   *tgt_box_extents,
 int           **tbox_box_idx,
 int           **tbox_box
)
{
  // tbox = tree_box
  int *box_tbox_idx = NULL;
  int *box_tbox     = NULL;
  PDM_box_tree_intersect_boxes_boxes(bt,
                                     i_copied_rank,
                                     n_tgt_box,
                                     tgt_box_extents,
                                     &box_tbox_idx,
                                     &box_tbox);

  /* Transpose from line->box to box->line */
  PDM_boxes_t *boxes;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
  }
  int n_boxes = boxes->n_boxes;

  // PDM_log_trace_connectivity_int(box_tbox_idx, box_tbox, n_tgt_box, "box_tbox ::");

  int* tbox_box_n = PDM_array_zeros_int(n_boxes);
  for(int i = 0; i < n_tgt_box; ++i) {
    for(int idx_box = box_tbox_idx[i]; idx_box <  box_tbox_idx[i+1]; ++idx_box ) {
      tbox_box_n[box_tbox[idx_box]]++;
    }
  }

  *tbox_box_idx = PDM_array_new_idx_from_sizes_int(tbox_box_n, n_boxes);
  *tbox_box = (int *) malloc(sizeof(int) * (*tbox_box_idx)[n_boxes]);
  PDM_array_reset_int(tbox_box_n, n_boxes, 0);
  for (int iline = 0; iline < n_tgt_box; iline++) {
    for (int idx_box = box_tbox_idx[iline]; idx_box < box_tbox_idx[iline+1]; idx_box++) {
      int ibox = box_tbox[idx_box];
      (*tbox_box)[(*tbox_box_idx)[ibox] + tbox_box_n[ibox]++] = iline;
    }
  }
  free(tbox_box_n);
  free(box_tbox_idx);
  free(box_tbox);

}


/**
 *
 * \brief Get an indexed list of all boxes inside given volumes
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt                    Pointer to box tree structure
 * \param [in]   i_copied_rank         Copied rank
 * \param [in]   n_volumes             Number of volumes
 * \param [in]   n_planes_per_volume   Index of the number of planes per volume
 * \param [in]   plane_normal          Oriented normal vector for a given plane (oriented toward the interior of the volume)
 * \param [in]   plane_pt_coord        Point on plane coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  volume_box_idx        Pointer to the index array on lines (size = \ref n_line + 1)
 * \param [out]  volume_box_l_num      Pointer to the list of boxes intersecting lines (size = \ref box_idx[\ref n_line])
 *
 */

void
PDM_box_tree_intersect_volume_boxes
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_volumes,
 const int      *volume_plane_idx,
 double         *plane_normal,
 double         *plane_pt_coord,
 int           **volume_box_idx,
 int           **volume_box_l_num
 )
 {
  if (i_copied_rank >= bt->n_copied_ranks) {
    PDM_error(__FILE__, __LINE__, 0, "Copied rank %d >= Number of copied ranks\n", (int) i_copied_rank);
  }

  const int dim = bt->boxes->dim;

  // Total number of planes
  int n_planes = volume_plane_idx[n_volumes];

  double *plane_pt_coord_normalized = malloc(sizeof(double) * 3 * n_planes);
  double *plane_normal_normalized   = malloc(sizeof(double) * 3 * n_planes);

  // Normalize coordinates
  PDM_box_set_normalize_robust((PDM_box_set_t *) bt->boxes,
                                                 n_planes,
                                                 plane_pt_coord,
                                                 plane_pt_coord_normalized);

  PDM_box_set_normalize_normal_vector((PDM_box_set_t *) bt->boxes,
                                                        n_planes,
                                                        plane_normal,
                                                        plane_normal_normalized);

  // Different pointer depending if it is a copied rank or not
  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
    box_tree_data = bt->local_data;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int n_boxes = boxes->n_boxes;

  // Set up output
  *volume_box_idx      = PDM_array_const_int((n_volumes + 1), 0);
  int *_volume_box_idx = *volume_box_idx;

  int tmp_s_boxes = 4 * n_volumes;
  *volume_box_l_num = malloc (sizeof(int) * tmp_s_boxes);
  int *_volume_box_l_num = *volume_box_l_num;

  // Set up stack
  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  // Set up tracking if boxes have been dealt with yet
  int n_visited_boxes = 0;
  PDM_bool_t *is_visited_box = malloc(sizeof(PDM_bool_t) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    is_visited_box[i] = PDM_FALSE;
  }

  int *visited_boxes = malloc(sizeof(int) * n_boxes);

  double *current_plane_normals = NULL;
  double *current_pt_planes     = NULL;
  int current_n_plane;

  double current_node_box_extents[2*dim];

  for (int ivolume = 0; ivolume < n_volumes; ivolume++) {

    _volume_box_idx[ivolume+1] = _volume_box_idx[ivolume];

    /* Init stack */
    pos_stack = 0;

    // Get current volume information
    current_plane_normals = plane_normal_normalized   + 3*volume_plane_idx[ivolume];
    current_pt_planes     = plane_pt_coord_normalized + 3*volume_plane_idx[ivolume];
    current_n_plane       = volume_plane_idx[ivolume+1] - volume_plane_idx[ivolume];

    // Get current box information
    _extents (dim,
              box_tree_data->nodes[0].morton_code,
              current_node_box_extents);

    if (_box_intersect_volume(current_n_plane, current_plane_normals, current_pt_planes, current_node_box_extents)) {
      stack[pos_stack++] = 0;
    }

    n_visited_boxes = 0;

    /* Traverse box tree */
    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      _node_t *node = &(box_tree_data->nodes[node_id]);

      if (node->n_boxes == 0) {
        continue;
      }

      /* Leaf node */
      if (node->is_leaf) {
        /* inspect boxes contained in current leaf node */
        for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          int box_id = box_tree_data->box_ids[node->start_id + ibox];

          if (is_visited_box[box_id] == PDM_FALSE) {
            double *box_extents = boxes->extents + box_id*2*dim;

            if (_box_intersect_volume(current_n_plane, current_plane_normals, current_pt_planes, box_extents)) {
              // Avoid table overflow
              if (_volume_box_idx[ivolume+1] >= tmp_s_boxes) {
                tmp_s_boxes *= 2;
                *volume_box_l_num = realloc (*volume_box_l_num, sizeof(int) * tmp_s_boxes);
                _volume_box_l_num = *volume_box_l_num;
              }
              _volume_box_l_num[_volume_box_idx[ivolume+1]++] = box_id;

            } // end if box is in volume
            visited_boxes[n_visited_boxes++] = box_id;
            is_visited_box[box_id] = PDM_TRUE;
          } // end if not yet dealt with this box

        } // end loop on leaf boxes
      } // end if is a leaf node

      /* Internal node */
      else {
        /* inspect children of current node */
        int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
        for (int ichild = 0; ichild < bt->n_children; ichild++) {
          int child_id = child_ids[ichild];
          _extents (dim, box_tree_data->nodes[child_id].morton_code, current_node_box_extents);

          if (_box_intersect_volume(current_n_plane, current_plane_normals, current_pt_planes, current_node_box_extents)) {
            stack[pos_stack++] = child_id;
          }
        } // end loop on children of curent node
      } // end if is an internal node

    } // end while dealing with stack

    // Clean up visited information of current volume
    for (int ibox = 0; ibox < n_visited_boxes; ibox++) {
      is_visited_box[visited_boxes[ibox]] = PDM_FALSE;
    }

  } // end loop on volumes

  // Free's
  if (plane_pt_coord != plane_pt_coord_normalized) {
    free (plane_pt_coord_normalized);
  }
  if (plane_normal != plane_normal_normalized) {
    free (plane_normal_normalized);
  }
  free (is_visited_box);
  free (stack);
  free (visited_boxes);


  *volume_box_l_num = realloc (*volume_box_l_num, sizeof(int) * _volume_box_idx[n_volumes]);
}


void
PDM_box_tree_write_vtk
(
 const char     *filename,
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       normalized
 )
{
  assert(bt != NULL);

  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    box_tree_data = bt->local_data;
  } else {
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int n_nodes = box_tree_data->n_nodes;
  int dim     = dim = bt->boxes->dim;

  double *node_extents = (double *) malloc(sizeof(double) * n_nodes * 6);
  int    *node_depth   = PDM_array_zeros_int(n_nodes);


  /* Extents */
  for (int i = 0; i < n_nodes; i++) {
    double *e = node_extents + 6*i;
    _extents(dim,
             box_tree_data->nodes[i].morton_code,
             e);

    if (!normalized) {
      double en[6] = {e[0], e[1], e[2], e[3], e[4], e[5]};
      PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en,   e);
      PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en+3, e+3);
    }
  }

  /* Depth */
  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  stack[pos_stack++] = 0;
  while (pos_stack > 0) {

    int node_id = stack[--pos_stack];

    int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
    _node_t *node = &(box_tree_data->nodes[node_id]);

    if (node->is_leaf) {
      continue;
    }

    for (int ichild = 0; ichild < bt->n_children; ichild++) {
      int child_id = child_ids[ichild];
      node_depth[child_id] = node_depth[node_id] + 1;
      stack[pos_stack++] = child_id;
    }

  }
  free(stack);


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nbox_tree\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    double *e = node_extents + 6*i;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_nodes, 9*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "12\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_nodes);
  fprintf(f, "FIELD node_field 2\n");

  fprintf(f, "depth 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "%d\n", node_depth[i]);
  }

  fprintf(f, "is_leaf 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    _node_t *node = &(box_tree_data->nodes[i]);
    fprintf(f, "%d\n", node->is_leaf);
  }

  fclose(f);

  free(node_extents);
  free(node_depth);

}


void
PDM_box_tree_write_vtk2
(
 const char     *filename,
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 double *s,
 double *d
 )
{
  assert(bt != NULL);

  PDM_box_tree_data_t *box_tree_data;
  if (i_copied_rank < 0) {
    box_tree_data = bt->local_data;
  } else {
    box_tree_data = bt->rank_data + i_copied_rank;
  }

  int n_nodes = box_tree_data->n_nodes;
  int dim     = dim = bt->boxes->dim;

  double *node_extents = (double *) malloc(sizeof(double) * n_nodes * 6);
  int    *node_depth   = PDM_array_zeros_int(n_nodes);


  /* Extents */
  for (int i = 0; i < n_nodes; i++) {
    double *e = node_extents + 6*i;
    _extents(dim,
             box_tree_data->nodes[i].morton_code,
             e);

    if (1) {//!normalized) {
      double en[6] = {e[0], e[1], e[2], e[3], e[4], e[5]};
      // PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en,   e);
      // PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en+3, e+3);
      for (int j = 0; j < 3; j++) {
        e[j] = en[j] * d[j] + s[j];
        e[3+j] = en[3+j] * d[j] + s[j];
      }
    }
  }

  /* Depth */
  int s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  stack[pos_stack++] = 0;
  while (pos_stack > 0) {

    int node_id = stack[--pos_stack];

    int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
    _node_t *node = &(box_tree_data->nodes[node_id]);

    if (node->is_leaf) {
      continue;
    }

    for (int ichild = 0; ichild < bt->n_children; ichild++) {
      int child_id = child_ids[ichild];
      node_depth[child_id] = node_depth[node_id] + 1;
      stack[pos_stack++] = child_id;
    }

  }
  free(stack);


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nbox_tree\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    double *e = node_extents + 6*i;
    fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
    fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
    fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
    fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
  }

  fprintf(f, "CELLS %d %d\n", n_nodes, 9*n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    int j = 8*i;
    fprintf(f, "8 %d %d %d %d %d %d %d %d\n", j, j+1, j+2, j+3, j+4, j+5, j+6, j+7);
  }

  fprintf(f, "CELL_TYPES %d\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "12\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_nodes);
  fprintf(f, "FIELD node_field 2\n");

  fprintf(f, "depth 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "%d\n", node_depth[i]);
  }

  fprintf(f, "is_leaf 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    _node_t *node = &(box_tree_data->nodes[i]);
    fprintf(f, "%d\n", node->is_leaf);
  }

  fclose(f);

  free(node_extents);
  free(node_depth);

}


void
PDM_box_tree_extract_extents
(
 PDM_box_tree_t  *bt,
 const int        normalized,
 const int        depth_max,
       int       *n_extract_boxes,
       double   **extract_extents,
       int       *n_extract_child,
       int      **extract_child_id
)
{
  /*
   * Il faudrait aussi sortir les child_id des noeuds extrait
   *   Si on a le node_id -> On peut rajouter un poids fonction de la solicitation
   *   box_tree_data->nodes[child_id].extra_weight = 0
   */
  assert(bt != NULL);

  PDM_box_tree_data_t *box_tree_data = bt->local_data;

  int n_nodes = box_tree_data->n_nodes;
  int dim     = bt->boxes->dim;

  /* Depth */
  int  s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int    *node_depth        = PDM_array_zeros_int(n_nodes);
  double *_extract_extents  = malloc(n_nodes * 6 * sizeof(double));
  int    *_extract_child_id = malloc(n_nodes     * sizeof(int   ));
  int     _n_extract_boxes  = 0;
  int     _n_extract_child  = 0;

  stack[pos_stack++] = 0;
  while (pos_stack > 0) {

    int node_id = stack[--pos_stack];

    int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
    _node_t *node = &(box_tree_data->nodes[node_id]);

    if (node_depth[node_id] < depth_max && node->is_leaf) {
      continue;
    }

    for (int ichild = 0; ichild < bt->n_children; ichild++) {
      int child_id = child_ids[ichild];
      node_depth[child_id] = node_depth[node_id] + 1;


      if(node_depth[child_id] < depth_max) {
        stack[pos_stack++] = child_id;
      }

      // log_trace("node->n_boxes = %i\n", box_tree_data->nodes[child_id].n_boxes);
      if(box_tree_data->nodes[child_id].n_boxes == 0) {
        continue;
      }

      double *e = _extract_extents + 6*_n_extract_boxes;
      _extents(dim, box_tree_data->nodes[child_id].morton_code, e);

      _extract_child_id[_n_extract_child++] = child_id;

      // log_trace("\n");

      if (!normalized) {
        double en[6] = {e[0], e[1], e[2], e[3], e[4], e[5]};
        PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en,   e  );
        PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en+3, e+3);
      }
      _n_extract_boxes++;
    }
  }
  free(stack);
  free(node_depth);

  _extract_extents  = realloc(_extract_extents , _n_extract_boxes * 6 * sizeof(double));
  _extract_child_id = realloc(_extract_child_id, _n_extract_child     * sizeof(int   ));

  *n_extract_child  = _n_extract_child;
  *n_extract_boxes  = _n_extract_boxes;
  *extract_extents  = _extract_extents;
  *extract_child_id = _extract_child_id;
}



void
PDM_box_tree_extract_leaves
(
 PDM_box_tree_t  *bt,
       int       *n_leaf,
       int      **leaf_id
)
{
  /*
   * Il faudrait aussi sortir les child_id des noeuds extrait
   *   Si on a le node_id -> On peut rajouter un poids fonction de la solicitation
   *   box_tree_data->nodes[child_id].extra_weight = 0
   */
  assert(bt != NULL);

  PDM_box_tree_data_t *box_tree_data = bt->local_data;

  int n_nodes = box_tree_data->n_nodes;

  /* Depth */
  int  s_stack = ((bt->n_children - 1) * (bt->max_level - 1) + bt->n_children);
  int *stack = malloc ((sizeof(int)) * s_stack);
  int pos_stack = 0;

  int    *_leaf_id = malloc(n_nodes     * sizeof(int   ));
  int     _n_leaf  = 0;

  stack[pos_stack++] = 0;
  while (pos_stack > 0) {

    int node_id = stack[--pos_stack];

    int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
    _node_t *node = &(box_tree_data->nodes[node_id]);

    if (node->is_leaf) {
      _leaf_id[_n_leaf++] = node_id;
    }
    else {
      for (int ichild = 0; ichild < bt->n_children; ichild++) {
        int child_id = child_ids[ichild];
        if (child_id < 0 || box_tree_data->nodes[child_id].n_boxes == 0) {
          continue;
        }
        stack[pos_stack++] = child_id;
      }
    }
  }
  free(stack);

  // _extract_extents  = realloc(_extract_extents , _n_extract_boxes * 6 * sizeof(double));
  // _extract_child_id = realloc(_extract_child_id, _n_extract_child     * sizeof(int   ));

  // *n_extract_child  = _n_extract_child;
  // *n_extract_boxes  = _n_extract_boxes;
  // *extract_extents  = _extract_extents;
  // *extract_child_id = _extract_child_id;

  _leaf_id = realloc(_leaf_id, _n_leaf * sizeof(int));
  *n_leaf  = _n_leaf;
  *leaf_id = _leaf_id;
}


void
PDM_box_tree_extract_node_extents
(
 PDM_box_tree_t  *bt,
 int              n_node,
 int             *node_id,
 double          *node_extents,
 const int        normalized
 )
{
  int dim = 3;
  PDM_box_tree_data_t *box_tree_data = bt->local_data;
  double *btree_s, *btree_d;
  PDM_box_set_normalization_get((PDM_box_set_t *) bt->boxes,
                                &btree_s,
                                &btree_d);

  for (int inode = 0; inode < n_node; inode++) {
    if (normalized) {
      _extents(dim,
               box_tree_data->nodes[node_id[inode]].morton_code,
               node_extents + 6*inode);
    }
    else {
      _extents_real(dim,
                    box_tree_data->nodes[node_id[inode]].morton_code,
                    btree_s,
                    btree_d,
                    node_extents + 6*inode);
    }
  }
}


void
PDM_box_tree_extract_extents_by_child_ids
(
 PDM_box_tree_t  *bt,
 const int        normalized,
 const int        n_child_to_extract,
 const int       *child_ids_to_extract,
       int       *n_extract_boxes,
       double   **extract_extents,
       int       *n_extract_child,
       int      **extract_child_id,
       int      **extract_is_leaf
)
{
  /*
   * Il faudrait aussi sortir les child_id des noeuds extrait
   *   Si on a le node_id -> On peut rajouter un poids fonction de la solicitation
   *   box_tree_data->nodes[child_id].extra_weight = 0
   */
  assert(bt != NULL);

  PDM_box_tree_data_t *box_tree_data = bt->local_data;

  int dim     = bt->boxes->dim;

  /* Depth */
  double *_extract_extents  = malloc(n_child_to_extract * bt->n_children * 6 * sizeof(double));
  int    *_extract_child_id = malloc(n_child_to_extract * bt->n_children     * sizeof(int   ));
  int    *_extract_is_leaf  = malloc(n_child_to_extract * bt->n_children     * sizeof(int   ));
  int     _n_extract_boxes  = 0;
  int     _n_extract_child  = 0;

  for(int i = 0; i < n_child_to_extract; ++i) {
    int node_id = child_ids_to_extract[i];

    int *child_ids = box_tree_data->child_ids + node_id*bt->n_children;
    _node_t *node = &(box_tree_data->nodes[node_id]);

    // log_trace("node->n_boxes = %i \n", node->n_boxes);

    if (node->is_leaf) {
      // log_trace("is leaf --> %i \n", node_id);
      continue;
    }

    for (int ichild = 0; ichild < bt->n_children; ichild++) {
      int child_id = child_ids[ichild];

      // log_trace("\t box_tree_data->nodes[child_id].n_boxes = %i \n", box_tree_data->nodes[child_id].n_boxes);
      // log_trace("node->n_boxes = %i\n", box_tree_data->nodes[child_id].n_boxes);
      if(box_tree_data->nodes[child_id].n_boxes == 0) {
        continue;
      }

      double *e = _extract_extents + 6*_n_extract_boxes;
      _extents(dim, box_tree_data->nodes[child_id].morton_code, e);

      _extract_child_id[_n_extract_child  ] = child_id;
      if(box_tree_data->nodes[child_id].is_leaf)  {
        _extract_is_leaf [_n_extract_child++] = 1;
      } else {
        _extract_is_leaf [_n_extract_child++] = 0;
      }


      if (!normalized) {
        double en[6] = {e[0], e[1], e[2], e[3], e[4], e[5]};
        PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en,   e  );
        PDM_box_set_normalize_inv((PDM_box_set_t *) bt->boxes, en+3, e+3);
      }
      _n_extract_boxes++;
    }
  }

  _extract_extents  = realloc(_extract_extents , _n_extract_boxes * 6 * sizeof(double));
  _extract_child_id = realloc(_extract_child_id, _n_extract_child     * sizeof(int   ));
  _extract_is_leaf  = realloc(_extract_is_leaf , _n_extract_child     * sizeof(int   ));

  *n_extract_child  = _n_extract_child;
  *n_extract_boxes  = _n_extract_boxes;
  *extract_extents  = _extract_extents;
  *extract_child_id = _extract_child_id;
  *extract_is_leaf  = _extract_is_leaf;
}


void
PDM_box_tree_assign_weight
(
 PDM_box_tree_t  *bt,
 const int        n_node,
 const int       *nodes_id,
       int       *weight
)
{

  PDM_box_tree_data_t *box_tree_data = bt->local_data;


  for(int i = 0; i < n_node; ++i) {
    int node_id = nodes_id[i];
    _node_t *node = &(box_tree_data->nodes[node_id]);
    node->extra_weight += weight[i];
  }
}


int
PDM_box_tree_get_box_ids
(
 PDM_box_tree_t  *bt,
 int              node_id,
 int            **box_ids
)
{

  PDM_boxes_t *boxes = bt->boxes->local_boxes;

  int n_boxes = boxes->n_boxes;

  int *is_visited = PDM_array_zeros_int(n_boxes);


  // PDM_box_tree_data_t *box_tree_data = bt->local_data;

  int dim              = bt->boxes->dim;
  int n_box_ids_approx = 0;
  _count_box_ids(bt, dim, node_id, &n_box_ids_approx);

  // log_trace("n_box_ids_approx = %i \n", n_box_ids_approx);
  /* Allocate */
  int* _box_ids = malloc( n_box_ids_approx * sizeof(int));

  /* Collect */
  int n_box_ids = 0;
  _collect_box_ids(bt, dim, node_id, is_visited, _box_ids, &n_box_ids);

  // log_trace("n_box_ids = %i \n", n_box_ids);

  free(is_visited);


  _box_ids = realloc(_box_ids, n_box_ids * sizeof(int));
  *box_ids = _box_ids;

  return n_box_ids;
}



int
PDM_box_tree_box_extents_get
(
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 double         **extents
 )
{
  assert(bt != NULL);

  PDM_boxes_t *boxes;
  if (i_copied_rank < 0) {
    boxes = bt->boxes->local_boxes;
  } else {
    boxes = bt->boxes->rank_boxes + i_copied_rank;
  }

  *extents = boxes->extents;

  return boxes->n_boxes;
}


static void
_visu_pair
(
 char                  *filename,
 PDM_box_tree_t        *btree,
 PDM_point_tree_seq_t  *ptree,
 int                    btree_node_id,
 int                    ptree_node_id
 )
{

  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;

  boxes         = btree->boxes->local_boxes;
  box_tree_data = btree->local_data;


  _node_t *node = &(box_tree_data->nodes[btree_node_id]);

  int n_box_in_leaf = 0;
  if (node->is_leaf) {
    n_box_in_leaf = node->n_boxes;
  }

  int n_pts_in_leaf = 0;
  if (ptree->nodes->is_leaf[ptree_node_id]) {
    n_pts_in_leaf = ptree->nodes->n_points[ptree_node_id];
  }

  double *btree_s, *btree_d;
  PDM_box_set_normalization_get((PDM_box_set_t *) btree->boxes,
                                &btree_s,
                                &btree_d);

  double btree_extents[6];
  _extents_real(3,
                box_tree_data->nodes[btree_node_id].morton_code,
                btree_s,
                btree_d,
                btree_extents);


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nintersection_btree_ptree\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 16 + 8*n_box_in_leaf + n_pts_in_leaf);

  double *e;

  e = btree_extents;
  fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
  fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
  fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
  fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
  fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
  fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
  fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
  fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);

  e = ptree->nodes->extents + 6*ptree_node_id;
  fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
  fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
  fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
  fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
  fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
  fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
  fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
  fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);

  if (node->is_leaf) {
    for (int ibox = 0; ibox < n_box_in_leaf; ibox++) {
      int box_id = box_tree_data->box_ids[node->start_id + ibox];

      const double *box_extents = boxes->extents + box_id*6;
      for (int i = 0; i < 3; i++) {
        btree_extents[i  ] = btree_s[i] + btree_d[i]*box_extents[i  ];
        btree_extents[i+3] = btree_s[i] + btree_d[i]*box_extents[i+3];
      }
      e = btree_extents;
      fprintf(f, "%f %f %f\n", e[0], e[1], e[2]);
      fprintf(f, "%f %f %f\n", e[3], e[1], e[2]);
      fprintf(f, "%f %f %f\n", e[3], e[4], e[2]);
      fprintf(f, "%f %f %f\n", e[0], e[4], e[2]);
      fprintf(f, "%f %f %f\n", e[0], e[1], e[5]);
      fprintf(f, "%f %f %f\n", e[3], e[1], e[5]);
      fprintf(f, "%f %f %f\n", e[3], e[4], e[5]);
      fprintf(f, "%f %f %f\n", e[0], e[4], e[5]);
    }
  }

  if (ptree->nodes->is_leaf[ptree_node_id]) {
    for (int i = ptree->nodes->range[2*ptree_node_id]; i < ptree->nodes->range[2*ptree_node_id+1]; i++) {
      double *p = ptree->_pts_coord + 3*i;
      fprintf(f, "%f %f %f\n", p[0], p[1], p[2]);
    }
  }

  fprintf(f, "CELLS %d %d\n",
          2 + n_box_in_leaf + n_pts_in_leaf,
          9*(2+n_box_in_leaf) + 2*n_pts_in_leaf);
  for (int i = 0; i < 2; i++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      int k = 8*i + j;
      fprintf(f, "%d ", k);
    }
    fprintf(f, "\n");
  }

  for (int i = 0; i < n_box_in_leaf; i++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      int k = 16 + 8*i + j;
      fprintf(f, "%d ", k);
    }
    fprintf(f, "\n");
  }

  for (int i = 0; i < n_pts_in_leaf; i++) {
    fprintf(f, "1 %d\n", 16 + 8*n_box_in_leaf + i);
  }

  fprintf(f, "CELL_TYPES %d\n", 2 + n_box_in_leaf + n_pts_in_leaf);
  for (int i = 0; i < 2 + n_box_in_leaf; i++) {
    fprintf(f, "12\n");
  }
  for (int i = 0; i < n_pts_in_leaf; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", 2 + n_box_in_leaf + n_pts_in_leaf);
  fprintf(f, "FIELD field 2\n");

  fprintf(f, "itree 1 %d int\n", 2 + n_box_in_leaf + n_pts_in_leaf);
  fprintf(f, "0\n1\n");
  for (int i = 0; i < n_box_in_leaf; i++) {
    fprintf(f, "0\n");
  }
  for (int i = 0; i < n_pts_in_leaf; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "id 1 %d int\n", 2 + n_box_in_leaf + n_pts_in_leaf);
  fprintf(f, "-1\n-1\n");

  if (node->is_leaf) {
    for (int i = 0; i < n_box_in_leaf; i++) {
      int box_id = box_tree_data->box_ids[node->start_id + i];
      fprintf(f, "%d\n", box_id);
    }
  }

  if (ptree->nodes->is_leaf[ptree_node_id]) {
    for (int i = ptree->nodes->range[2*ptree_node_id]; i < ptree->nodes->range[2*ptree_node_id+1]; i++) {
      fprintf(f, "%d\n", ptree->new_to_old[i]);
    }
  }

  fclose(f);
}


static int
_binary_search
(
 const int  elem,
 const int *array,
 const int  n
 )
{
  int l = 0;
  int r = n;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] < elem)
    return l + 1;
  else
    return l;
}

static inline void
_insertion_sort
(
 const int   point_id,
 int        *box_pts_n,
 int        *box_pts_s,
 int       **box_pts
 )
{
  int dbg = 0;

  if ((*box_pts_n) == 0) {
    (*box_pts)[(*box_pts_n)++] = point_id;
  }
  else {
    int i = _binary_search(point_id,
                           *box_pts,
                           *box_pts_n);
    if (dbg) {
      log_trace("%d at pos %d in array ", point_id, i);
      PDM_log_trace_array_int((*box_pts), (*box_pts_n), "");
    }

    if ((*box_pts)[i] == point_id) {
      return;
    }

    if ((*box_pts_s) <= (*box_pts_n)) {
      (*box_pts_s) *= 2;
      (*box_pts)    = realloc((*box_pts),
                              sizeof(int) * (*box_pts_s));
    }

    for (int j = (*box_pts_n); j > i; j--) {
      (*box_pts)[j] = (*box_pts)[j-1];
    }

    (*box_pts)[i] = point_id;
    (*box_pts_n)++;

    if (dbg) {
      PDM_log_trace_array_int(*box_pts, *box_pts_n, "after insertion : ");
    }
  }
}


typedef enum {
  SUBDIVISION_CRITERION_VOLUME,
  SUBDIVISION_CRITERION_LENGTH,
  SUBDIVISION_CRITERION_DEPTH,
} _subdivision_criterion_t;

void
PDM_tree_intersection_point_box
(
 PDM_box_tree_t        *btree,
 PDM_point_tree_seq_t  *ptree,
 int                  **box_pts_idx,
 int                  **box_pts
 )
{
  _subdivision_criterion_t subdiv_crit = SUBDIVISION_CRITERION_VOLUME;

  int dbg  = 0;
  int visu = 0;

  PDM_boxes_t *boxes;
  PDM_box_tree_data_t *box_tree_data;

  boxes         = btree->boxes->local_boxes;
  box_tree_data = btree->local_data;

  int n_boxes = boxes->n_boxes;

  double *btree_s, *btree_d;
  PDM_box_set_normalization_get((PDM_box_set_t *) btree->boxes,
                                &btree_s,
                                &btree_d);


  /* Get point_tree data (use gets!!!) */
  int ptree_n_children = PDM_point_tree_n_children_get(ptree);
  int    *ptree_depth       = ptree->nodes->depth;
  int    *ptree_is_leaf     = ptree->nodes->is_leaf;
  int    *ptree_range       = ptree->nodes->range;
  int    *ptree_children_id = ptree->nodes->children_id;
  double *ptree_extents     = ptree->nodes->extents;

  // int n_pts = ptree->n_pts;
  double *ptree_pts_coord;
  PDM_point_tree_seq_sorted_points_get(ptree,
                                       &ptree_pts_coord);

  double *_pts_coord = ptree_pts_coord;
  // double *_pts_coord = malloc (sizeof(double) * 3 * n_pts);
  // PDM_box_set_normalize_robust ((PDM_box_set_t *) btree->boxes,
  //                               n_pts,
  //                    (double *) ptree_pts_coord,
  //                               _pts_coord);
  // log_trace("normalization %f %f %f  --->  %f %f %f\n",
  //           ptree_pts_coord[0], ptree_pts_coord[1], ptree_pts_coord[2],
  //           _pts_coord[0], _pts_coord[1], _pts_coord[2]);

  int *ptree_new_to_old = NULL;
  PDM_point_tree_seq_point_new_to_old_get(ptree,
                                          &ptree_new_to_old);



  double btree_extents[6];

  /* Start from both roots */
  int btree_node_id = 0;
  int ptree_node_id = 0;

  _extents_real(3,
                box_tree_data->nodes[btree_node_id].morton_code,
                btree_s,
                btree_d,
                btree_extents);
  // log_trace("btree_extents : %f %f %f  %f %f %f",
  //           btree_extents[0], btree_extents[1], btree_extents[2],
  //           btree_extents[3], btree_extents[4], btree_extents[5]);
  int intersect = _intersect_box_box(3,
                                     btree_extents,
                                     &ptree_extents[6*ptree_node_id]);

  if (!intersect) {
    *box_pts_idx = PDM_array_zeros_int(n_boxes + 1);
    if (dbg) {
      log_trace("roots do not intersect\n");
    }
    return;
  }


  int s_queue = 1000; // ?
  int *queue0 = malloc(sizeof(int) * s_queue * 2);
  int *queue1 = malloc(sizeof(int) * s_queue * 2);
  int *queues[2] = {queue0, queue1};

  int *queue0_depth = NULL;
  int *queue1_depth = NULL;
  if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
    queue0_depth = malloc(sizeof(int) * s_queue);
    queue1_depth = malloc(sizeof(int) * s_queue);
  }

  int *queues_depth[2] = {queue0_depth, queue1_depth};


  int n_queue = 0;
  queues[0][2*n_queue  ] = btree_node_id;
  queues[0][2*n_queue+1] = ptree_node_id;
  if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
    queues_depth[0][n_queue] = 0;
  }
  n_queue++;



  int  *__box_pts_n = PDM_array_zeros_int(n_boxes);
  int  *__box_pts_s = PDM_array_const_int(n_boxes, 4);
  int **__box_pts   = malloc(sizeof(int *) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    __box_pts[i] = malloc(sizeof(int) * __box_pts_s[i]);
  }


  int istep = -1;
  while (n_queue > 0) {

    int new_n_queue = 0;

    if (dbg) {
      PDM_log_trace_array_int(queues[0], 2*n_queue, "queue : ");
    }

    for (int ipair = 0; ipair < n_queue; ipair++) {

      istep++;

      btree_node_id = queues[0][2*ipair  ];
      ptree_node_id = queues[0][2*ipair+1];
      int btree_depth;
      if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
        btree_depth = queues_depth[0][ipair];
      }

      if (dbg) {
        log_trace("  Step %d\n", istep);
        log_trace("  node ids %d %d\n", btree_node_id, ptree_node_id);
      }

      if (visu) {
        char filename[999];
        sprintf(filename, "intersection_btree_ptree_step_%4.4d.vtk", istep);
        _visu_pair(filename,
                   btree,
                   ptree,
                   btree_node_id,
                   ptree_node_id);
      }

      _node_t *node = &(box_tree_data->nodes[btree_node_id]);

      int btree_is_leaf = node->is_leaf;

      int isubdiv;

      /* Which tree do we subdivide? */
      if (btree_is_leaf) {

        if (ptree_is_leaf[ptree_node_id]) {
          /* Both leaves */
          isubdiv = -1;

          if (dbg) {
            log_trace("  both leaves\n");
          }

          /* inspect boxes contained in current leaf node */
          for (int ibox = 0; ibox < node->n_boxes; ibox++) {

            int box_id = box_tree_data->box_ids[node->start_id + ibox];

            if (dbg) {
              log_trace("    box_id = %d\n", box_id);
            }

            const double *box_extents = boxes->extents + box_id*6;
            for (int i = 0; i < 3; i++) {
              btree_extents[i  ] = btree_s[i] + btree_d[i]*box_extents[i  ];
              btree_extents[i+3] = btree_s[i] + btree_d[i]*box_extents[i+3];
            }
            // log_trace(" box_extents : %f %f %f  %f %f %f\n",
            //           box_extents[0], box_extents[1], box_extents[2],
            //           box_extents[3], box_extents[4], box_extents[5]);
            // log_trace("_box_extents : %f %f %f  %f %f %f\n",
            //           btree_extents[0], btree_extents[1], btree_extents[2],
            //           btree_extents[3], btree_extents[4], btree_extents[5]);

            for (int ipt = ptree_range[2*ptree_node_id]; ipt < ptree_range[2*ptree_node_id+1]; ipt++) {

              double *pt = _pts_coord + 3*ipt;
              int point_id = ptree_new_to_old[ipt];

              if (dbg) {
                log_trace("      point_id = %d (%f %f %f)\n",
                          point_id, pt[0], pt[1], pt[2]);
              }

              if (_point_inside_box(3, btree_extents, pt)) {
                if (dbg) {
                  log_trace("        inside box\n");
                }

                _insertion_sort(point_id,
                                &__box_pts_n[box_id],
                                &__box_pts_s[box_id],
                                &__box_pts  [box_id]);
              }

            } // End of loop on current ptree leaf's points
          } // End of loop on current btree leaf's boxes


        }
        else {
          /* Subdivide point tree */
          isubdiv = 1;
        }

      }

      else if (ptree_is_leaf[ptree_node_id]) {
        /* Subdivide box tree */
        isubdiv = 0;
      }

      else {

        // Decide which tree is subdivided
        _extents_real(3,
                      box_tree_data->nodes[btree_node_id].morton_code,
                      btree_s,
                      btree_d,
                      btree_extents);

        double btree_crit = 1.;
        double ptree_crit = 1.;
        switch (subdiv_crit) {

          case SUBDIVISION_CRITERION_VOLUME: {
            btree_crit = 1.;
            ptree_crit = 1.;
            for (int i = 0; i < 3; i++) {
              btree_crit *= (btree_extents[i+3] - btree_extents[i]);
              ptree_crit *= (ptree_extents[6*ptree_node_id + i+3] - ptree_extents[6*ptree_node_id + i]);
            }
            break;
          }
          case SUBDIVISION_CRITERION_LENGTH: {
            btree_crit = 0.;
            ptree_crit = 0.;
            for (int i = 0; i < 3; i++) {
              btree_crit = PDM_MAX(btree_crit, (btree_extents[i+3] - btree_extents[i]));
              ptree_crit = PDM_MAX(ptree_crit, (ptree_extents[6*ptree_node_id + i+3] - ptree_extents[6*ptree_node_id + i]));
            }
            break;
          }
          case SUBDIVISION_CRITERION_DEPTH: {
            btree_crit = btree_depth;
            ptree_crit = ptree_depth[ptree_node_id];
            break;
          }
          default: {
            PDM_error(__FILE__, __LINE__, 0,
                      "Subdivision criterion %d not implemented\n", (int) subdiv_crit);
            break;
          }
        }

        if (btree_crit > ptree_crit) {
          /* Subdivide box tree */
          isubdiv = 0;
        }
        else {
          /* Subdivide point tree */
          isubdiv = 1;
        }

      }


      /* Add children to new queue */
      if (isubdiv == 0) {
        /* Subdivide box tree */
        if (dbg) {
          log_trace("  subdivide box tree\n");
        }

        int *children_id = box_tree_data->child_ids + btree_node_id*btree->n_children;
        for (int ichild = 0; ichild < btree->n_children; ichild++) {
          int child_id = children_id[ichild];

          if (child_id < 0) continue;

          _extents_real(3,
                        box_tree_data->nodes[child_id].morton_code,
                        btree_s,
                        btree_d,
                        btree_extents);
          // log_trace("btree_extents : %f %f %f  %f %f %f",
          //           btree_extents[0], btree_extents[1], btree_extents[2],
          //           btree_extents[3], btree_extents[4], btree_extents[5]);

          intersect = _intersect_box_box(3,
                                         btree_extents,
                                         &ptree_extents[6*ptree_node_id]);

          if (dbg) {
            log_trace("    child %d, intersect? %d\n", child_id, intersect);
          }

          if (intersect) {
            // Check size!!!
            if (new_n_queue >= s_queue) {
              s_queue *= 2;
              queues[0] = realloc(queues[0], sizeof(int) * s_queue * 2);
              queues[1] = realloc(queues[1], sizeof(int) * s_queue * 2);
              if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
                queues_depth[0] = realloc(queues_depth[0], sizeof(int) * s_queue);
                queues_depth[1] = realloc(queues_depth[1], sizeof(int) * s_queue);
              }
            }

            queues[1][2*new_n_queue  ] = child_id;
            queues[1][2*new_n_queue+1] = ptree_node_id;
            if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
              queues_depth[1][new_n_queue] = btree_depth + 1;
            }

            new_n_queue++;
          }
        }

      }
      else if (isubdiv == 1) {
        /* Subdivide point tree */
        if (dbg) {
          log_trace("  subdivide point tree\n");
        }

        _extents_real(3,
                      box_tree_data->nodes[btree_node_id].morton_code,
                      btree_s,
                      btree_d,
                      btree_extents);
        // log_trace("btree_extents : %f %f %f  %f %f %f",
        //           btree_extents[0], btree_extents[1], btree_extents[2],
        //           btree_extents[3], btree_extents[4], btree_extents[5]);

        int *children_id = ptree_children_id + ptree_node_id*ptree_n_children;
        for (int ichild = 0; ichild < ptree_n_children; ichild++) {
          int child_id = children_id[ichild];

          if (child_id < 0) continue;

          intersect = _intersect_box_box(3,
                                         btree_extents,
                                         &ptree_extents[6*child_id]);

          if (dbg) {
            log_trace("    child %d, intersect? %d\n", child_id, intersect);
          }

          if (intersect) {
            // Check size!!!
            if (new_n_queue >= s_queue) {
              s_queue *= 2;
              queues[0] = realloc(queues[0], sizeof(int) * s_queue * 2);
              queues[1] = realloc(queues[1], sizeof(int) * s_queue * 2);
              if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
                queues_depth[0] = realloc(queues_depth[0], sizeof(int) * s_queue);
                queues_depth[1] = realloc(queues_depth[1], sizeof(int) * s_queue);
              }
            }

            queues[1][2*new_n_queue  ] = btree_node_id;
            queues[1][2*new_n_queue+1] = child_id;
            if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
              queues_depth[1][new_n_queue] = btree_depth;
            }

            new_n_queue++;
          }
        }

      }


    } // End of loop in current queue

    if (dbg) {
      PDM_log_trace_array_int(queues[1], 2*new_n_queue, "new_queue : ");
    }

    /* Swap queues */
    int *tmp = queues[1];
    queues[1] = queues[0];
    queues[0] = tmp;
    if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
      tmp = queues_depth[1];
      queues_depth[1] = queues_depth[0];
      queues_depth[0] = tmp;
    }

    n_queue = new_n_queue;

  } // End of while loop
  free(__box_pts_s);
  // free(_pts_coord);
  free(queues[0]);
  free(queues[1]);
  if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
    free(queues_depth[0]);
    free(queues_depth[1]);
  }

  /* Re-arrange result */
  *box_pts_idx = PDM_array_new_idx_from_sizes_int(__box_pts_n, n_boxes);

  *box_pts = malloc(sizeof(int) * (*box_pts_idx)[n_boxes]);
  for (int i = 0; i < n_boxes; i++) {
    int *bp = *box_pts + (*box_pts_idx)[i];

    for (int j = 0; j < __box_pts_n[i]; j++) {
      bp[j] = __box_pts[i][j];
    }

    free(__box_pts[i]);
  }
  free(__box_pts_n);
  free(__box_pts);

  if(queue0_depth != NULL) {
    free(queue0_depth);
  }
  if(queue1_depth != NULL) {
    free(queue1_depth);
  }


}


#ifdef __cplusplus
}
#endif /* __cplusplus */
