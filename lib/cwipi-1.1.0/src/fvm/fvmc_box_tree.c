/*============================================================================
 * Search octrees and quadtrees of boxes.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_parall.h"

#include "fvmc_box_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_box_tree.h"

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

#define FVMC_BOX_TREE_MAX_BUILD_LOOPS 50

/* Structures for each octant or quadrant */
/*----------------------------------------*/

/* If the type is BOX_TREE_NODE, the ordering of children is defined as follows,
   using notation B: bottom, U: up, E: east, W: west, S: south,  N: north.

   octant:   0: BSW, 1: BSE, 2: BNW, 3: BNE, 4: USW, 5: USE, 6: UNW, 7: UNE
   quadrant: 0:  SW, 1:  SE, 2:  NW, 3:  NE
   segment:  0:   W, 1:   E
 */

typedef struct {

  _Bool              is_leaf;      /* True for leaf nodes */

  fvmc_morton_code_t  morton_code;  /* Level and coordinates in the grid
                                      according to Morton encoding */

  fvmc_lnum_t  n_boxes;             /* Number of associated bounding boxes */
  fvmc_lnum_t  start_id;            /* Position of the first box_id */

} _node_t;

/* Structure used to manage statistics */

typedef struct {

  unsigned    max_level_reached;  /* Max level number reached */

  fvmc_lnum_t  n_leaves;           /* Number of leaves in the tree */
  fvmc_lnum_t  n_boxes;            /* Number of boxes to locate in the tree */
  fvmc_lnum_t  n_linked_boxes;     /* Number of linked boxes in the tree */
  fvmc_lnum_t  n_spill_leaves;     /* Number of leaves where n_boxes > threshold */

  fvmc_lnum_t  min_linked_boxes;   /* Minimum number of boxes for a leaf */
  fvmc_lnum_t  max_linked_boxes;   /* Maximum number of boxes for a leaf */

} fvmc_box_tree_stats_t;

/* Main box tree structure */
/*-------------------------*/

struct _fvmc_box_tree_t {

  int               n_children;      /* 8, 4, or 2 (2^dim) */

  fvmc_morton_int_t  max_level;       /* Max. possible level */
  fvmc_lnum_t        threshold;       /* Max number of boxes linked to a
                                        node if max_level is not reached */
  fvmc_lnum_t        max_box_ratio;   /* Max n_linked_boxes / n_boxes value */

  fvmc_box_tree_stats_t stats;        /* Statistics related to the structure */

  fvmc_lnum_t        n_max_nodes;     /* Current max. allocated nodes */
  fvmc_lnum_t        n_nodes;         /* Number of nodes (including leaves) */

  _node_t          *nodes;           /* Array of nodes (root at index 0) */

  fvmc_lnum_t       *child_ids;       /* Ids of associated children
                                        (size: 2^dim * n_max_nodes) */
  fvmc_lnum_t       *box_ids;         /* List of associated box ids.
                                        size = stat.n_linked_boxes */

  int     n_build_loops;             /* Number of loops required to build */

#if defined(FVMC_HAVE_MPI)
  MPI_Comm          comm;            /* Associated MPI communicator */
#endif
};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

inline static const fvmc_coord_t *
_box_min(const fvmc_box_set_t   *boxes,
         fvmc_lnum_t             box_id)
{
  return boxes->extents + box_id*boxes->dim*2;
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

inline static const fvmc_coord_t *
_box_max(const fvmc_box_set_t   *boxes,
         fvmc_lnum_t             box_id)
{
  return boxes->extents + box_id*boxes->dim*2 + boxes->dim;
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
_boxes_intersect_3d(const fvmc_coord_t  *extents,
                    fvmc_lnum_t          id_0,
                    fvmc_lnum_t          id_1)
{
  const fvmc_coord_t *e0 = extents + id_0*6;
  const fvmc_coord_t *e1 = extents + id_1*6;

  if (   e0[0] > e1[3] || e1[0] > e0[3]
      || e0[1] > e1[4] || e1[1] > e0[4]
      || e0[2] > e1[5] || e1[2] > e0[5])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_2d(const fvmc_coord_t  *extents,
                    fvmc_lnum_t          id_0,
                    fvmc_lnum_t          id_1)
{
  const fvmc_coord_t *e0 = extents + id_0*4;
  const fvmc_coord_t *e1 = extents + id_1*4;

  if (   e0[0] > e1[2] || e1[0] > e0[2]
      || e0[1] > e1[3] || e1[1] > e0[3])
    return false;
  else
    return true;
}

inline static _Bool
_boxes_intersect_1d(const fvmc_coord_t  *extents,
                    fvmc_lnum_t          id_0,
                    fvmc_lnum_t          id_1)
{
  const fvmc_coord_t *e0 = extents + id_0*2;
  const fvmc_coord_t *e1 = extents + id_1*2;

  if (   e0[0] > e1[1] || e1[0] > e0[1])
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------
 * Update octree stat structure (min, max, mean, box ratio, ...)
 *
 * parameters:
 *   bt      <-> pointer on the fvmc_box_tree_t structure to deal with
 *   node_id <-- node on which we collect data
 *----------------------------------------------------------------------------*/

static void
_update_tree_stats(fvmc_box_tree_t  *bt,
                   fvmc_lnum_t       node_id)
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {
    int n_children = bt->n_children;
    const fvmc_lnum_t *_child_ids = bt->child_ids + node_id*bt->n_children;
    for (i = 0; i < n_children; i++)
      _update_tree_stats(bt, _child_ids[i]);
  }

  else { /* leaf node */

    fvmc_box_tree_stats_t  s = bt->stats;

    s.n_leaves += 1;
    s.n_linked_boxes += node->n_boxes;

    if (node->n_boxes > bt->threshold)
      s.n_spill_leaves += 1;

    s.min_linked_boxes = FVMC_MIN(s.min_linked_boxes, node->n_boxes);
    s.max_linked_boxes = FVMC_MAX(s.max_linked_boxes, node->n_boxes);
    s.max_level_reached = FVMC_MAX(s.max_level_reached, node->morton_code.L);

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
_get_box_tree_stats(fvmc_box_tree_t  *bt)
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

  if (bt->nodes != NULL)
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
_get_grid_coords_3d(fvmc_morton_int_t  level,
                    const double      coords[3],
                    double            XYZ[])
{
  fvmc_morton_int_t  refinement = 1 << level;

  XYZ[0] = coords[0] * refinement;
  XYZ[1] = coords[1] * refinement;
  XYZ[2] = coords[2] * refinement;
}

inline static void
_get_grid_coords_2d(fvmc_morton_int_t  level,
                    const double      coords[2],
                    double            XY[])
{
  fvmc_morton_int_t  refinement = 1 << level;

  XY[0] = coords[0] * refinement;
  XY[1] = coords[1] * refinement;
}

inline static void
_get_grid_coords_1d(fvmc_morton_int_t  level,
                    const double      coords[1],
                    double            X[])
{
  fvmc_morton_int_t  refinement = 1 << level;

  X[0] = coords[0] * refinement;
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
_node_intersect_box_3d(fvmc_morton_code_t   morton_code,
                       const fvmc_coord_t   min_box[3],
                       const fvmc_coord_t   max_box[3])
{
  int  i;
  double  min_oct[3], max_oct[3];

  for (i = 0; i < 3; i++) {
    min_oct[i] = (double)morton_code.X[i];
    max_oct[i] = (double)(morton_code.X[i] + 1);
  }

  if (   min_box[0] > max_oct[0] || min_oct[0] > max_box[0]
      || min_box[1] > max_oct[1] || min_oct[1] > max_box[1]
      || min_box[2] > max_oct[2] || min_oct[2] > max_box[2])
    return false;
  else
    return true;
}

inline static _Bool
_node_intersect_box_2d(fvmc_morton_code_t   morton_code,
                       const fvmc_coord_t   min_box[2],
                       const fvmc_coord_t   max_box[2])
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
_node_intersect_box_1d(fvmc_morton_code_t   morton_code,
                       const fvmc_coord_t   min_box[1],
                       const fvmc_coord_t   max_box[1])
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
_evaluate_splitting_3d(fvmc_box_tree_t       *bt,
                       const fvmc_box_set_t  *boxes,
                       fvmc_lnum_t            node_id)
{
  fvmc_lnum_t  i, j;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[8];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 3);

  /* Define a Morton code for each child */

  fvmc_morton_get_children(3, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[3], max_grid_coord[3];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(3, next_level, box_min);
    max_code = fvmc_morton_encode(3, next_level, box_max);

    if (   fvmc_morton_compare(3, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same octant */

      assert(   fvmc_morton_compare(3, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   fvmc_morton_compare(3, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_2d(fvmc_box_tree_t       *bt,
                       const fvmc_box_set_t  *boxes,
                       fvmc_lnum_t            node_id)
{
  fvmc_lnum_t  i, j;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[4];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  assert(boxes->dim == 2);

  /* Define a Morton code for each child */

  fvmc_morton_get_children(2, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[2], max_grid_coord[2];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(2, next_level, box_min);
    max_code = fvmc_morton_encode(2, next_level, box_max);

    if (   fvmc_morton_compare(2, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   fvmc_morton_compare(2, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   fvmc_morton_compare(2, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
}

static int
_evaluate_splitting_1d(fvmc_box_tree_t       *bt,
                       const fvmc_box_set_t  *boxes,
                       fvmc_lnum_t            node_id)
{
  fvmc_lnum_t  i, j;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[2];

  int  n_linked_boxes = 0;

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  /* Define a Morton code for each child */

  fvmc_morton_get_children(1, node.morton_code, children);

  /* Loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[1], max_grid_coord[1];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(1, next_level, box_min);
    max_code = fvmc_morton_encode(1, next_level, box_max);

    if (   fvmc_morton_compare(1, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          n_linked_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   fvmc_morton_compare(1, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   fvmc_morton_compare(1, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          n_linked_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  return n_linked_boxes;
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
_count_next_level(fvmc_box_tree_t           *bt,
                  const fvmc_box_set_t      *boxes,
                  fvmc_lnum_t                node_id,
                  fvmc_box_tree_sync_t       build_type,
                  fvmc_lnum_t               *next_level_size)
{
  fvmc_lnum_t  i;

  _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    assert(bt->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _count_next_level(bt,
                        boxes,
                        bt->child_ids[bt->n_children*node_id + i],
                        build_type,
                        next_level_size);

  }

  else { /* if (node->is_leaf == true) */
    if (   node->n_boxes < bt->threshold
        && node_id != 0                    /* Root node is always divided */
        && build_type == FVMC_BOX_TREE_ASYNC_LEVEL)
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
_recurse_tree_build(fvmc_box_tree_t       *bt,
                    const fvmc_box_set_t  *boxes,
                    fvmc_box_tree_sync_t   build_type,
                    fvmc_lnum_t           *next_size)
{
  int  state = 0;
  fvmc_lnum_t  _next_size = 0;

  _Bool retval = false;

#if defined(FVMC_HAVE_MPI)

  int  n_ranks = 1;
  MPI_Comm comm = boxes->comm;

  if (comm != MPI_COMM_NULL)
    MPI_Comm_size(comm, &n_ranks);

#endif

  bt->n_build_loops += 1;

  if (bt == NULL)
    state = 1;

  /* To avoid infinite loop on tree building */

  if (bt->n_build_loops > FVMC_BOX_TREE_MAX_BUILD_LOOPS)
    state = 1;

  /* A sufficient accuracy has been reached */

  if (bt->stats.max_level_reached == bt->max_level)
    state = 1;

  /* Algorithm is converged. No need to go further */

  if (   bt->stats.n_spill_leaves == 0
      && bt->stats.max_level_reached > 0)
    state = 1;

#if defined(FVMC_HAVE_MPI)
  if (n_ranks > 1 && build_type == FVMC_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    MPI_Allreduce(&state, &global_state, 1, MPI_INT, MPI_MIN, comm);
    state = global_state; /* Stop if all ranks require it */
  }
#endif

  if (state == 0) {

    int  box_ratio;

    /* Limit, to avoid excessive memory usage */

    _count_next_level(bt,
                      boxes,
                      0,  /* Starts from root */
                      build_type,
                      &_next_size);

    if (bt->stats.n_boxes > 0)
      box_ratio = _next_size/bt->stats.n_boxes;
    else
      box_ratio = 0;

    if (box_ratio > bt->max_box_ratio)
      state = 1;

  }

#if defined(FVMC_HAVE_MPI)
  if (n_ranks > 1 && build_type == FVMC_BOX_TREE_SYNC_LEVEL) {
    int global_state;
    MPI_Allreduce(&state, &global_state, 1, MPI_INT, MPI_MAX, comm);
    state = global_state; /* Stop as as soon as any rank requires it */
  }
#endif

  /* If no condition is encoutered, we have to continue */

  *next_size = _next_size;

  if (state == 0)
    retval = true;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a box tree by copying from another.
 *
 * parameters:
 *   dest <-> pointer to destination box tree
 *   src  <-- pointer to source box tree
 *----------------------------------------------------------------------------*/

static void
_copy_tree(fvmc_box_tree_t        *dest,
           const fvmc_box_tree_t  *src)
{
  assert(dest != NULL && src != NULL);

  memcpy(dest, src, sizeof(fvmc_box_tree_t));

  BFTC_MALLOC(dest->nodes, dest->n_max_nodes, _node_t);
  BFTC_MALLOC(dest->child_ids, dest->n_max_nodes*dest->n_children, fvmc_lnum_t);
  BFTC_MALLOC(dest->box_ids, (dest->stats).n_linked_boxes, fvmc_lnum_t);

  memcpy(dest->nodes, src->nodes, dest->n_nodes * sizeof(_node_t));
  memcpy(dest->child_ids,
         src->child_ids,
         dest->n_nodes * src->n_children * sizeof(fvmc_lnum_t));

  memcpy(dest->box_ids,
         src->box_ids,
         (dest->stats).n_linked_boxes * sizeof(fvmc_lnum_t));
}

/*----------------------------------------------------------------------------
 * Destroy a fvmc_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to fvmc_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

static void
_free_tree_arrays(fvmc_box_tree_t  *bt)
{
  assert(bt != NULL);

  BFTC_FREE(bt->nodes);
  BFTC_FREE(bt->child_ids);
  BFTC_FREE(bt->box_ids);
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
_new_node(fvmc_box_tree_t     *bt,
          fvmc_morton_code_t   morton_code,
          fvmc_lnum_t          node_id)
{
  int  i;
  _node_t *node;

  assert(bt != NULL);

  node = bt->nodes + node_id;

  if (morton_code.L > bt->max_level)
    bftc_error(__FILE__, __LINE__, 0,
              _("Error adding a new node in box tree (%p).\n"
                "Max level reached. Current level: %u and Max level: %d\n"),
              bt, morton_code.L, bt->max_level);

  node->is_leaf = true;
  node->morton_code = morton_code;

  node->n_boxes = 0;
  node->start_id = -1; /* invalid value by default */

  for (i = 0; i < bt->n_children; i++)
    bt->child_ids[node_id*bt->n_children + i] = -1;
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
_split_node_3d(fvmc_box_tree_t       *bt,
               fvmc_box_tree_t       *next_bt,
               const fvmc_box_set_t  *boxes,
               fvmc_lnum_t            node_id,
               fvmc_lnum_t           *shift_ids)
{
  int j, i;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[8];

  fvmc_lnum_t  n_linked_boxes = 0;
  fvmc_lnum_t  _shift_ids = *shift_ids;
  fvmc_lnum_t  n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 8);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 8 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    BFTC_REALLOC(next_bt->nodes, next_bt->n_max_nodes, _node_t);
    BFTC_REALLOC(next_bt->child_ids, next_bt->n_max_nodes*8, fvmc_lnum_t);
  }

  /* Define a Morton code for each child and create the children nodes */

  fvmc_morton_get_children(3, node.morton_code, children);

  for (i = 0; i < 8; i++) {
    const fvmc_lnum_t  new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*8 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 8;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[3], max_grid_coord[3];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(3, next_level, box_min);
    max_code = fvmc_morton_encode(3, next_level, box_max);

    if (   fvmc_morton_compare(3, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {
        if (_node_intersect_box_3d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same octant */

      assert(   fvmc_morton_compare(3, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {
        if (   fvmc_morton_compare(3, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 8; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 8; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[3], max_grid_coord[3];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(3, next_level, box_min);
    max_code = fvmc_morton_encode(3, next_level, box_max);

    if (   fvmc_morton_compare(3, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_3d(next_level, box_min, min_grid_coord);
      _get_grid_coords_3d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 8; i++) {

        if (_node_intersect_box_3d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same octant */

      assert(   fvmc_morton_compare(3, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 8; i++) {

        if (   fvmc_morton_compare(3, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 8; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 8 - 1].start_id
             + next_bt->nodes[n_init_nodes + 8 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

static void
_split_node_2d(fvmc_box_tree_t       *bt,
               fvmc_box_tree_t       *next_bt,
               const fvmc_box_set_t  *boxes,
               fvmc_lnum_t            node_id,
               fvmc_lnum_t           *shift_ids)
{
  int j, i;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[4];

  fvmc_lnum_t  n_linked_boxes = 0;
  fvmc_lnum_t  _shift_ids = *shift_ids;
  fvmc_lnum_t  n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 4);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 4 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    BFTC_REALLOC(next_bt->nodes, next_bt->n_max_nodes, _node_t);
    BFTC_REALLOC(next_bt->child_ids, next_bt->n_max_nodes*4, fvmc_lnum_t);
  }

  /* Define a Morton code for each child and create the children nodes */

  fvmc_morton_get_children(2, node.morton_code, children);

  for (i = 0; i < 4; i++) {
    const fvmc_lnum_t  new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*4 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 4;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[2], max_grid_coord[2];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(2, next_level, box_min);
    max_code = fvmc_morton_encode(2, next_level, box_max);

    if (   fvmc_morton_compare(2, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {
        if (_node_intersect_box_2d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same quadrant */

      assert(   fvmc_morton_compare(2, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {
        if (   fvmc_morton_compare(2, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 4; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 4; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[2], max_grid_coord[2];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(2, next_level, box_min);
    max_code = fvmc_morton_encode(2, next_level, box_max);

    if (   fvmc_morton_compare(2, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_2d(next_level, box_min, min_grid_coord);
      _get_grid_coords_2d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 4; i++) {

        if (_node_intersect_box_2d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same quadrant */

      assert(   fvmc_morton_compare(2, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 4; i++) {

        if (   fvmc_morton_compare(2, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 4; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 4 - 1].start_id
             + next_bt->nodes[n_init_nodes + 4 - 1].n_boxes));
#endif

  /* Return pointers */

  *shift_ids = _shift_ids;
}

static void
_split_node_1d(fvmc_box_tree_t       *bt,
               fvmc_box_tree_t       *next_bt,
               const fvmc_box_set_t  *boxes,
               fvmc_lnum_t            node_id,
               fvmc_lnum_t           *shift_ids)
{
  int j, i;
  fvmc_morton_code_t  min_code, max_code;
  fvmc_morton_code_t  children[2];

  fvmc_lnum_t  n_linked_boxes = 0;
  fvmc_lnum_t  _shift_ids = *shift_ids;
  fvmc_lnum_t  n_init_nodes = next_bt->n_nodes;
  _node_t  split_node = next_bt->nodes[node_id];

  const _node_t  node = bt->nodes[node_id];
  const fvmc_morton_int_t  next_level = node.morton_code.L + 1;

  assert(bt->n_children == 2);

  /* Add the leaves to the next_bt structure */

  if (n_init_nodes + 2 > next_bt->n_max_nodes) {
    assert(next_bt->n_max_nodes > 0);
    next_bt->n_max_nodes *= 2;
    BFTC_REALLOC(next_bt->nodes, next_bt->n_max_nodes, _node_t);
    BFTC_REALLOC(next_bt->child_ids, next_bt->n_max_nodes*2, fvmc_lnum_t);
  }

  /* Define a Morton code for each child and create the children nodes */

  fvmc_morton_get_children(1, node.morton_code, children);

  for (i = 0; i < 2; i++) {
    const fvmc_lnum_t  new_id = n_init_nodes + i;
    next_bt->child_ids[node_id*2 + i] = new_id;
    _new_node(next_bt, children[i], new_id);
  }

  split_node.start_id = 0;
  split_node.n_boxes = 0;
  split_node.is_leaf = false;

  next_bt->nodes[node_id] = split_node;
  next_bt->n_nodes = n_init_nodes + 2;

  /* Counting loop on boxes associated to the node_id */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[2], max_grid_coord[2];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(1, next_level, box_min);
    max_code = fvmc_morton_encode(1, next_level, box_max);

    if (   fvmc_morton_compare(1, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {
        if (_node_intersect_box_1d(children[i], min_grid_coord, max_grid_coord))
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
      }

    }
    else { /* Box is included in the same segment */

      assert(   fvmc_morton_compare(1, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {
        if (   fvmc_morton_compare(1, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {
          (next_bt->nodes[n_init_nodes + i]).n_boxes += 1;
          break;
        }
      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

  /* Build index */

  for (i = 0; i < 2; i++) {

    (next_bt->nodes[n_init_nodes + i]).start_id
      = _shift_ids + n_linked_boxes;
    n_linked_boxes += (next_bt->nodes[n_init_nodes + i]).n_boxes;

  }

  _shift_ids += n_linked_boxes;

  for (i = 0; i < 2; i++)
    (next_bt->nodes[n_init_nodes + i]).n_boxes = 0;

  /* Second loop on boxes associated to the node_id: fill */

  for (j = 0; j < node.n_boxes; j++) {

    fvmc_coord_t  min_grid_coord[2], max_grid_coord[2];

    fvmc_lnum_t  box_id = bt->box_ids[node.start_id + j];
    const fvmc_coord_t  *box_min = _box_min(boxes, box_id);
    const fvmc_coord_t  *box_max = _box_max(boxes, box_id);

    min_code = fvmc_morton_encode(1, next_level, box_min);
    max_code = fvmc_morton_encode(1, next_level, box_max);

    if (   fvmc_morton_compare(1, min_code, max_code)
        == FVMC_MORTON_DIFFERENT_ID) {

      _get_grid_coords_1d(next_level, box_min, min_grid_coord);
      _get_grid_coords_1d(next_level, box_max, max_grid_coord);

      for (i = 0; i < 2; i++) {

        if (_node_intersect_box_1d(children[i],
                                   min_grid_coord,
                                   max_grid_coord)) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;
          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
        }

      } /* End of loop on children*/

    }
    else { /* Box is included in the same segment */

      assert(   fvmc_morton_compare(1, max_code, min_code)
             == FVMC_MORTON_EQUAL_ID);

      for (i = 0; i < 2; i++) {

        if (   fvmc_morton_compare(1, min_code, children[i])
            == FVMC_MORTON_EQUAL_ID) {

          const fvmc_lnum_t sub_id = n_init_nodes + i;
          const fvmc_lnum_t shift =   (next_bt->nodes[sub_id]).n_boxes
                                   + (next_bt->nodes[sub_id]).start_id;

          next_bt->box_ids[shift] = box_id;
          (next_bt->nodes[sub_id]).n_boxes += 1;
          break;
        }

      } /* End of loop on children */

    } /* If min_code and max_code in the same leaf */

  } /* End of loop on boxes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < 2; i++) {
    _node_t  n1 = next_bt->nodes[n_init_nodes + i - 1];
    _node_t  n2 = next_bt->nodes[n_init_nodes + i];
    assert(n1.n_boxes == (n2.start_id - n1.start_id));
  }
  assert(   _shift_ids
         == (  next_bt->nodes[n_init_nodes + 2 - 1].start_id
             + next_bt->nodes[n_init_nodes + 2 - 1].n_boxes));
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
_build_next_level(fvmc_box_tree_t       *bt,
                  fvmc_box_tree_t       *next_bt,
                  const fvmc_box_set_t  *boxes,
                  fvmc_lnum_t            node_id,
                  fvmc_box_tree_sync_t   build_type,
                  fvmc_lnum_t           *shift_ids)
{
  fvmc_lnum_t  i;

  fvmc_lnum_t  _shift_ids = *shift_ids;
  const _node_t  *cur_node = bt->nodes + node_id;

  if (cur_node->is_leaf == false) {

    assert(bt->child_ids[bt->n_children*node_id] > 0);

    for (i = 0; i < bt->n_children; i++)
      _build_next_level(bt,
                        next_bt,
                        boxes,
                        bt->child_ids[bt->n_children*node_id + i],
                        build_type,
                        &_shift_ids);
  }

  else { /* if (node->is_leaf == true) */

    if (   cur_node->n_boxes < bt->threshold
        && node_id != 0                    /* Root node is always divided */
        && build_type == FVMC_BOX_TREE_ASYNC_LEVEL) {

      /* Copy related box_ids in the new next_ids */

      _node_t *next_node = next_bt->nodes + node_id;

      next_node->n_boxes = cur_node->n_boxes;
      next_node->start_id = _shift_ids;

      for (i = 0; i < cur_node->n_boxes; i++)
        next_bt->box_ids[_shift_ids++]
          = bt->box_ids[cur_node->start_id + i];
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

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Loop on all nodes of the box tree to define an array with Morton codes
 * and weights (= number of linked boxes) associated to each leaf
 *
 * parameters:
 *   bt         <-> pointer to fvmc_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   node_id    <-- id of the current node (to traverse)
 *   n_leaves   <-> current number of leaves in the tree with n_boxes > 0
 *   leaf_codes <-> Morton code associated to each leaf
 *   weight     <-> number of boxes attached to each leaf
 *----------------------------------------------------------------------------*/

static void
_build_leaf_weight(const fvmc_box_tree_t  *bt,
                   fvmc_lnum_t             node_id,
                   fvmc_lnum_t            *n_leaves,
                   fvmc_morton_code_t     *leaf_codes,
                   fvmc_lnum_t            *weight)
{
  int  i;

  fvmc_lnum_t _n_leaves = *n_leaves;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false)
    for (i = 0; i < bt->n_children; i++)
      _build_leaf_weight(bt,
                         bt->child_ids[bt->n_children*node_id + i],
                         &_n_leaves,
                         leaf_codes,
                         weight);

  else { /* node is a leaf */

    if (node->n_boxes > 0) {
      leaf_codes[_n_leaves] = node->morton_code;
      weight[_n_leaves] = node->n_boxes;
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
 *   bt           <-> pointer to fvmc_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_index(const fvmc_box_tree_t  *bt,
                         fvmc_box_distrib_t     *distrib,
                         int                    dim,
                         fvmc_lnum_t             node_id,
                         size_t                 size,
                         fvmc_morton_code_t      search_index[],
                         int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_index(bt,
                               distrib,
                               dim,
                               bt->child_ids[bt->n_children*node_id + i],
                               size,
                               search_index,
                               id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = fvmc_morton_binary_search(size,
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
 *   bt           <-> pointer to fvmc_box_tree_t structure.
 *   distrib      <-- structure holding box distribution data
 *   dim          <-- box tree layout dimension (1, 2, or 3)
 *   node_id      <-- id of the current node (to traverse)
 *   counter      <-> counter array used to build the list
 *   size         <-- size of index in which we search
 *   search_index <-- index on which box distribution is made
 *   id_rank      <-- relation between id and rank
 *----------------------------------------------------------------------------*/

static void
_build_rank_to_box_list(const fvmc_box_tree_t  *bt,
                        fvmc_box_distrib_t     *distrib,
                        int                    dim,
                        fvmc_lnum_t             node_id,
                        fvmc_lnum_t             counter[],
                        size_t                 size,
                        fvmc_morton_code_t      search_index[],
                        int                    id_rank[])
{
  int  i;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++)
      _build_rank_to_box_list(bt,
                              distrib,
                              dim,
                              bt->child_ids[bt->n_children*node_id + i],
                              counter,
                              size,
                              search_index,
                              id_rank);
  }
  else {

    if (node->n_boxes > 0) {

      int  id = fvmc_morton_binary_search(size,
                                         node->morton_code,
                                         search_index);
      int  rank = id_rank[id];

      for (i = 0; i < node->n_boxes; i++) {

        fvmc_lnum_t  box_id = bt->box_ids[node->start_id + i];
        fvmc_lnum_t  shift = distrib->index[rank] + counter[rank];

        distrib->list[shift] = box_id;
        counter[rank] += 1;

      }
    }
  }

}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Recursively build an index on boxes which intersect.
 *
 * parameters:
 *   bt      <-> pointer to fvmc_box_tree_t structure.
 *   boxes   <-- pointer to associated boxes structure
 *   node_id <-- id of the current node (to traverse)
 *   count   <-> intersection count
 *----------------------------------------------------------------------------*/

static void
_count_intersections(const fvmc_box_tree_t  *bt,
                     const fvmc_box_set_t   *boxes,
                     fvmc_lnum_t             node_id,
                     fvmc_lnum_t             count[])
{
  fvmc_lnum_t  i, j;

  const fvmc_coord_t  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _count_intersections(bt,
                           boxes,
                           bt->child_ids[bt->n_children*node_id + i],
                           count);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intersect_3d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intersect_2d(box_extents, id0, id1)) {
            count[id0] += 1;
            count[id1] += 1;
          }

        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];
          if (_boxes_intersect_1d(box_extents, id0, id1)) {
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
 *   bt        <-> pointer to fvmc_box_tree_t structure.
 *   boxes     <-- pointer to associated boxes structure
 *   node_id   <-- id of the current node (to traverse)
 *   count     <-> intersection count (workin array)
 *   index     <-- index on intersections
 *   box_g_num <-> global number of intersection boxes
 *----------------------------------------------------------------------------*/

static void
_get_intersections(const fvmc_box_tree_t  *bt,
                   const fvmc_box_set_t   *boxes,
                   fvmc_lnum_t             node_id,
                   fvmc_lnum_t             count[],
                   fvmc_lnum_t             box_index[],
                   fvmc_gnum_t             box_g_num[])
{
  int  i, j;

  const fvmc_coord_t  *box_extents = boxes->extents;
  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _get_intersections(bt,
                         boxes,
                         bt->child_ids[bt->n_children*node_id + i],
                         count,
                         box_index,
                         box_g_num);
  }
  else { /* node is a leaf */

    if (boxes->dim == 3) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intersect_3d(box_extents, id0, id1)) {
            fvmc_lnum_t  shift0 = box_index[id0] + count[id0];
            fvmc_lnum_t  shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }
    else if (boxes->dim == 2) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intersect_2d(box_extents, id0, id1)) {
            fvmc_lnum_t  shift0 = box_index[id0] + count[id0];
            fvmc_lnum_t  shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
            count[id0] += 1;
            count[id1] += 1;
          }
        }
      }
    }

    else if (boxes->dim == 1) {

      for (i = 0; i < node->n_boxes - 1; i++) {
        for (j = i+1; j < node->n_boxes; j++) {
          fvmc_lnum_t  id0 = bt->box_ids[node->start_id + i];
          fvmc_lnum_t  id1 = bt->box_ids[node->start_id + j];

          if (_boxes_intersect_1d(box_extents, id0, id1)) {
            fvmc_lnum_t  shift0 = box_index[id0] + count[id0];
            fvmc_lnum_t  shift1 = box_index[id1] + count[id1];
            box_g_num[shift0] = boxes->g_num[id1];
            box_g_num[shift1] = boxes->g_num[id0];
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
 *   bt      <-- pointer to fvmc_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *   n_steps <-- number of steps in histogram
 *   step    <-- steps of the histogram
 *   h_min   <-- min. value of the histogram
 *   counter <-> counter (working array)
 *----------------------------------------------------------------------------*/

static void
_build_histogram(const fvmc_box_tree_t  *bt,
                 fvmc_lnum_t             node_id,
                 fvmc_lnum_t             n_steps,
                 fvmc_lnum_t             step,
                 fvmc_lnum_t             h_min,
                 fvmc_gnum_t             count[])
{
  int  i, j;

  const _node_t  *node = bt->nodes + node_id;

  if (node->is_leaf == false) {

    for (i = 0; i < bt->n_children; i++) /* traverse downwards */
      _build_histogram(bt,
                       bt->child_ids[bt->n_children*node_id + i],
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
 *   bt      <-- pointer to fvmc_box_tree_t structure.
 *   node_id <-- id of the current node (to traverse)
 *----------------------------------------------------------------------------*/

static void
_dump_node(const fvmc_box_tree_t  *bt,
           fvmc_lnum_t             node_id)
{
  int  i;

  const char *node_type[] = {"node", "leaf"};

  const _node_t  *node = bt->nodes + node_id;
  const fvmc_morton_code_t  m_code = node->morton_code;

  bftc_printf("\n"
             "  node %10d (%s)\n"
             "    level:   %3hu - anchor: [ %10u %10u %10u ]\n"
             "    n_boxes: %3hd - start_id: %u\n"
             "    boxes:\n",
             node_id, node_type[(int)(node->is_leaf)],
             m_code.L, m_code.X[0], m_code.X[1], m_code.X[2],
             node->n_boxes, node->start_id);

  for (i = 0; i < node->n_boxes; i++)
    bftc_printf("        %d\n", (int)(bt->box_ids[node->start_id + i]));

  if (node->is_leaf == false) {

    const fvmc_lnum_t *c_id = bt->child_ids + bt->n_children*node_id;

    if (bt->n_children == 8)
      bftc_printf("  children_id:  %d %d %d %d %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3],
                 (int)c_id[4], (int)c_id[5], (int)c_id[6], (int)c_id[7]);
    else if (bt->n_children == 4)
      bftc_printf("  children_id:  %d %d %d %d\n",
                 (int)c_id[0], (int)c_id[1], (int)c_id[2], (int)c_id[3]);
    else if (bt->n_children == 2)
      bftc_printf("  children_id:  %d %d\n",
                 (int)c_id[0], (int)c_id[1]);

    for (i = 0; i < bt->n_children; i++)
      _dump_node(bt, c_id[i]);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a fvmc_box_tree_t structure and initialize it.
 *
 * parameters:
 *  max_level     <-- max possible level
 *  threshold     <-- max number of  boxes linked to an octant if
 *                    max_level is not reached
 *  max_box_ratio <-- max n_linked_boxes / n_boxes value
 *
 * returns:
 *   pointer to an empty fvmc_box_tree_t structure.
 *----------------------------------------------------------------------------*/

fvmc_box_tree_t *
fvmc_box_tree_create(int  max_level,
                    int  threshold,
                    int  max_box_ratio)
{
  fvmc_box_tree_t  *bt = NULL;

  BFTC_MALLOC(bt, 1, fvmc_box_tree_t);

  /* Sanity checks */

  if (max_level < 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("  Forbidden max_level value (%d) in the tree structure\n"),
              max_level);

  if (threshold < 1)
    bftc_error(__FILE__, __LINE__, 0,
              _("  Forbidden threshold value (%d) in the tree structure\n"),
              threshold);

  if (max_box_ratio < 1)
    bftc_error(__FILE__, __LINE__, 0,
              _("  Forbidden max_box_ratio value (%d) in the tree structure\n"),
              max_box_ratio);

  /* Create and initialize tree structure according to its type */

  bt->max_level = max_level;
  bt->threshold = threshold;
  bt->max_box_ratio = max_box_ratio;

#if defined(FVMC_HAVE_MPI)
  bt->comm = MPI_COMM_NULL;
#endif

  /* Set stats */

  bt->stats.max_level_reached = 0;

  bt->stats.n_leaves = 0;
  bt->stats.n_spill_leaves = 0;
  bt->stats.n_linked_boxes = 0;

  bt->stats.min_linked_boxes = INT_MAX;
  bt->stats.max_linked_boxes = 0;

  /* Initialize nodes */

  bt->n_max_nodes = 0;
  bt->n_nodes = 0;

  bt->nodes = NULL;

  bt->box_ids = NULL;

  bt->n_build_loops = 0;

  return bt;
}

/*----------------------------------------------------------------------------
 * Destroy a fvmc_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to fvmc_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

void
fvmc_box_tree_destroy(fvmc_box_tree_t  **bt)
{
  fvmc_box_tree_t  *_bt = *bt;

  if (_bt != NULL) {

    BFTC_FREE(_bt->nodes);
    BFTC_FREE(_bt->child_ids);
    BFTC_FREE(_bt->box_ids);

    BFTC_FREE(_bt);
    *bt = _bt;
  }
}

/*----------------------------------------------------------------------------
 * Get the deepest level allowed by the tree structure.
 *
 * parameters:
 *   bt <-- pointer to fvmc_box_tree_t structure.
 *
 * returns:
 *   deepest allowed level of the tree
 *----------------------------------------------------------------------------*/

int
fvmc_box_tree_get_max_level(const fvmc_box_tree_t  *bt)
{
  return bt->max_level;
}

/*----------------------------------------------------------------------------
 * Assign a set of boxes to an empty fvmc_box_tree_t structure.
 *
 * The box tree structure must have been created using to fvmc_tree_create().
 *
 * The depth of the tree is adjusted so that a maximum of max_n_elts boxes
 * will be assigned to each leaf, unless this would require going beyond
 * the tree's maximum level.
 *
 * If max_level = -1, the highest level reachable is FVMC_TREE_MAX_LEVEL but
 * there is no defined target level.
 *
 * parameters:
 *   bt         <-> pointer to fvmc_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   build_type <-- layout variant for building the tree structure
 *----------------------------------------------------------------------------*/

void
fvmc_box_tree_set_boxes(fvmc_box_tree_t       *bt,
                       const fvmc_box_set_t  *boxes,
                       fvmc_box_tree_sync_t   build_type)
{
  fvmc_lnum_t  box_id;

  fvmc_box_tree_t  tmp_bt;

  fvmc_lnum_t  next_box_ids_size = 0, shift = 0;
  fvmc_coord_t anchor[3] = {0., 0., 0.};

  /* Initialization */

  assert(bt != NULL);

  bt->n_build_loops = 0;

#if defined(FVMC_HAVE_MPI)
  bt->comm = boxes->comm;
#endif

  /* Preallocate for the two first levels of a tree */

  if (boxes->dim == 3) {
    bt->n_children = 8;
    bt->n_max_nodes = 73;
  }
  else if (boxes->dim == 2) {
    bt->n_children = 4;
    bt->n_max_nodes = 21;
  }
  else if (boxes->dim == 1) {
    bt->n_children = 2;
    bt->n_max_nodes = 7;
  }

  bt->n_nodes = 1;

  BFTC_MALLOC(bt->nodes, bt->n_max_nodes, _node_t);
  BFTC_MALLOC(bt->child_ids,
             bt->n_max_nodes*bt->n_children,
             fvmc_lnum_t);

  /* Define root node */

  _new_node(bt, fvmc_morton_encode(boxes->dim, 0, anchor), 0);

  /* Initialize bt by assigning all boxes to the root leaf */

  BFTC_MALLOC(bt->box_ids, boxes->n_boxes, fvmc_lnum_t);

  for (box_id = 0; box_id < boxes->n_boxes; box_id++)
    bt->box_ids[box_id] = box_id;

  (bt->nodes[0]).is_leaf = true;
  (bt->nodes[0]).n_boxes = boxes->n_boxes;
  (bt->nodes[0]).start_id = 0;

  bt->stats.n_boxes = boxes->n_boxes;

  _get_box_tree_stats(bt);

  /* Build local tree structure by adding boxes from the root */

  while (_recurse_tree_build(bt,
                             boxes,
                             build_type,
                             &next_box_ids_size)) {

    /* Initialize next_bt: copy of bt */

    _copy_tree(&tmp_bt, bt);

    /* Optimize memory usage */

    bt->n_max_nodes = bt->n_nodes;
    BFTC_REALLOC(bt->nodes, bt->n_nodes, _node_t);
    BFTC_REALLOC(bt->child_ids,
                bt->n_max_nodes*bt->n_children,
                fvmc_lnum_t);

    /* Define a box ids list for the next level of the boxtree */

    BFTC_REALLOC(tmp_bt.box_ids, next_box_ids_size, fvmc_lnum_t);
    shift = 0;

    _build_next_level(bt,
                      &tmp_bt,
                      boxes,
                      0, /* Starts from root */
                      build_type,
                      &shift);

    assert(shift == next_box_ids_size);

    /* replace current tree by the tree computed at a higher level */

    _free_tree_arrays(bt);
    *bt = tmp_bt; /* Overwrite bt members with those of next_bt */

    _get_box_tree_stats(bt);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bftc_printf("  - New box tree level -\n");
    fvmc_box_tree_dump_statistics(bt);
#endif

  } /* While building should continue */
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute an index based on Morton encoding to ensure a good distribution
 * of boxes among the participating ranks.
 *
 * parameters:
 *   bt         <-> pointer to fvmc_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *
 * returns:
 *   pointer to newly created fvmc_box_distrib_t structure.
 *----------------------------------------------------------------------------*/

fvmc_box_distrib_t *
fvmc_box_tree_get_distrib(fvmc_box_tree_t        *bt,
                         const fvmc_box_set_t   *boxes)
{
  int  i;

  int  reduce_size = 0;
  fvmc_lnum_t  n_leaves = 0;
  int  *reduce_ids = NULL;
  fvmc_morton_code_t  *leaf_codes = NULL, *reduce_index = NULL;
  fvmc_lnum_t  *weight = NULL, *counter = NULL;

  fvmc_box_distrib_t  *distrib = NULL;

  assert(bt != NULL);
  assert(boxes != NULL);

  /* Compute basic box distribution */

  distrib = fvmc_box_distrib_create(boxes->n_boxes,
                                   boxes->n_g_boxes,
                                   (bt->stats).max_level_reached,
                                   boxes->comm);

  if (distrib == NULL)
    return NULL;

  BFTC_MALLOC(leaf_codes, bt->stats.n_leaves, fvmc_morton_code_t);
  BFTC_MALLOC(weight, bt->stats.n_leaves, fvmc_lnum_t);

  /* Build index for boxes */

  _build_leaf_weight(bt,
                     0,
                     &n_leaves,
                     leaf_codes,
                     weight);

  assert(n_leaves <= bt->stats.n_leaves);

  BFTC_REALLOC(leaf_codes, n_leaves, fvmc_morton_code_t);
  BFTC_REALLOC(weight, n_leaves, fvmc_lnum_t);

  /* Compute the resulting Morton index */

  fvmc_box_set_build_morton_index(boxes,
                                 distrib,
                                 n_leaves,
                                 leaf_codes,
                                 weight);

  BFTC_FREE(leaf_codes);
  BFTC_FREE(weight);

  /* Compact Morton_index to get an array without "0 element" */

  for (i = 0; i < distrib->n_ranks; i++)
    if (fvmc_morton_a_gt_b(distrib->morton_index[i+1],
                          distrib->morton_index[i]))
      reduce_size++;

  BFTC_MALLOC(reduce_index, reduce_size + 1, fvmc_morton_code_t);
  BFTC_MALLOC(reduce_ids, reduce_size, int);

  reduce_size = 0;
  reduce_index[0] = distrib->morton_index[0];

  for (i = 0; i < distrib->n_ranks; i++) {

    if (fvmc_morton_a_gt_b(distrib->morton_index[i+1],
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

  BFTC_MALLOC(distrib->list,
             distrib->index[distrib->n_ranks], int);

  BFTC_MALLOC(counter, distrib->n_ranks, fvmc_lnum_t);

  for (i = 0; i < distrib->n_ranks; i++)
    counter[i] = 0;

  _build_rank_to_box_list(bt,
                          distrib,
                          boxes->dim,
                          0,  /* starts from root */
                          counter,
                          reduce_size,
                          reduce_index,
                          reduce_ids);

  /* Free memory */

  BFTC_FREE(counter);
  BFTC_FREE(reduce_ids);
  BFTC_FREE(reduce_index);

  /* Define the final index (without redundancies) and realloc list */

  fvmc_box_distrib_clean(distrib);

  return distrib;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Build an indexed list on boxes to list intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   boxes     <-- pointer to a associated box set
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxes
 *----------------------------------------------------------------------------*/

void
fvmc_box_tree_get_intersects(fvmc_box_tree_t        *bt,
                            const fvmc_box_set_t   *boxes,
                            fvmc_lnum_t            *box_index[],
                            fvmc_gnum_t            *box_g_num[])
{
  fvmc_lnum_t  i, list_size;

  fvmc_lnum_t  *counter = NULL;
  fvmc_lnum_t  *_index = NULL;
  fvmc_gnum_t  *_g_num = NULL;

  /* Build index */

  BFTC_MALLOC(_index, boxes->n_boxes + 1, fvmc_lnum_t);

  for (i = 0; i < boxes->n_boxes + 1; i++)
    _index[i] = 0;

  _count_intersections(bt,
                       boxes,
                       0, /* start from root */
                       _index + 1);

  /* Build index from counts */

  for (i = 0; i < boxes->n_boxes; i++)
    _index[i+1] += _index[i];

  list_size = _index[boxes->n_boxes];

  BFTC_MALLOC(_g_num, list_size, fvmc_gnum_t);

  BFTC_MALLOC(counter, boxes->n_boxes, fvmc_lnum_t);

  for (i = 0; i < boxes->n_boxes; i++)
    counter[i] = 0;

  /* Build list */

  _get_intersections(bt,
                     boxes,
                     0, /* start from root */
                     counter,
                     _index,
                     _g_num);

  BFTC_FREE(counter);

  /* Return pointers */

  *box_index = _index;
  *box_g_num = _g_num;
}

/*----------------------------------------------------------------------------
 * Get global box tree statistics.
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
 * parameters:
 *   bt                 <-- pointer to box tree structure
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_used           --> theoretical used memory
 *   mem_allocated      --> theoretical allocated memory
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
fvmc_box_tree_get_stats(const fvmc_box_tree_t  *bt,
                       int                    depth[3],
                       fvmc_lnum_t             n_leaves[3],
                       fvmc_lnum_t             n_boxes[3],
                       fvmc_lnum_t             n_threshold_leaves[3],
                       fvmc_lnum_t             n_leaf_boxes[3],
                       size_t                 mem_used[3],
                       size_t                 mem_allocated[3])
{
  int i;
  size_t mem_per_node;
  size_t s_mean[7], s_min[7], s_max[7];
  fvmc_box_tree_stats_t s;

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

  mem_per_node = sizeof(_node_t) + bt->n_children*sizeof(fvmc_lnum_t);

  s_mean[5] = sizeof(fvmc_box_tree_t);
  s_mean[5] += bt->n_nodes * mem_per_node;
  s_mean[5] += s.n_linked_boxes * sizeof(fvmc_lnum_t);

  s_mean[5] += sizeof(fvmc_box_set_t);
  s_mean[5] += s.n_boxes * (  sizeof(fvmc_gnum_t)
                           + (dim * 2 * sizeof(fvmc_coord_t)));

  s_mean[6] = s_mean[5] + (bt->n_max_nodes - bt->n_nodes)*mem_per_node;

 /* Pre-synchronize for serial cases (i = 0 already handled) */

  for (i = 1; i < 7; i++) {
    s_min[i] = s_mean[i];
    s_max[i] = s_mean[i];
  }

  /* In parallel mode, synchronize values */

#if defined(FVMC_HAVE_MPI)

  if (bt->comm != MPI_COMM_NULL) {

    int n_ranks;
    fvmc_gnum_t s_l_sum[14], s_g_sum[14];

    MPI_Comm_size(bt->comm, &n_ranks);

    if (n_ranks > 1) { /* Should always be the case, bat play it safe) */

      /* Split value to avoid exceeding fvmc_gnum_t limits
         (especially when it is a 32-bit value) */

      s_l_sum[0] = s.n_linked_boxes/n_ranks;
      s_l_sum[7] = s.n_linked_boxes%n_ranks;
      for (i = 1; i < 7; i++) {
        s_l_sum[i] = s_mean[i]/n_ranks;
        s_l_sum[i+7] = s_mean[i]%n_ranks;
      }

      MPI_Allreduce(s_l_sum, s_g_sum, 14, FVMC_MPI_GNUM, MPI_SUM, bt->comm);

      s_mean[0] = s.min_linked_boxes;
      MPI_Allreduce(s_mean, s_min, 7, FVMC_MPI_GNUM, MPI_MIN, bt->comm);
      s_mean[0] = s.max_linked_boxes;
      MPI_Allreduce(s_mean, s_max, 7, FVMC_MPI_GNUM, MPI_MAX, bt->comm);

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
        if (s_g_sum[i+7]%n_ranks >= (fvmc_gnum_t)n_ranks/2)
          s_mean[i] += 1;
      }
    }

  }
#endif

  /* Set values already in stats */

  if (depth != NULL) {
    depth[0] = s_mean[1];
    depth[1] = s_min[1];
    depth[2] = s_max[1];
  }

  if (n_leaves != NULL) {
    n_leaves[0] = s_mean[2];
    n_leaves[1] = s_min[2];
    n_leaves[2] = s_max[2];
  }

  if (n_boxes != NULL) {
    n_boxes[0] = s_mean[3];
    n_boxes[1] = s_min[3];
    n_boxes[2] = s_max[3];
  }

  if (n_threshold_leaves != NULL) {
    n_threshold_leaves[0] = s_mean[4];
    n_threshold_leaves[1] = s_min[4];
    n_threshold_leaves[2] = s_max[4];
  }

  if (n_leaf_boxes != NULL) {
    n_leaf_boxes[0] = s_mean[0];
    n_leaf_boxes[1] = s_min[0];
    n_leaf_boxes[2] = s_max[0];
  }

  if (mem_used != NULL) {
    mem_used[0] = s_mean[5];
    mem_used[1] = s_min[5];
    mem_used[2] = s_max[5];
  }

  if (mem_allocated != NULL) {
    mem_allocated[0] = s_mean[6];
    mem_allocated[1] = s_min[6];
    mem_allocated[2] = s_max[6];
  }

  return dim;
}

/*----------------------------------------------------------------------------
 * Display local statistics about a fvmc_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
fvmc_box_tree_dump_statistics(const fvmc_box_tree_t  *bt)
{
  int i, j;
  fvmc_box_tree_stats_t s;
  unsigned g_max_level_reached;
  fvmc_gnum_t n_g_leaves, n_g_boxes, n_g_linked_boxes, n_g_spill_leaves;
  fvmc_lnum_t g_min_linked_boxes, g_max_linked_boxes;
  double mean_linked_boxes, box_ratio;

  fvmc_gnum_t count[5];

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

#if defined(FVMC_HAVE_MPI)

  if (bt->comm != MPI_COMM_NULL) {

    fvmc_gnum_t l_min[1], g_min[1];
    fvmc_gnum_t l_max[2], g_max[2];
    fvmc_gnum_t l_sum[3], g_sum[3];

    l_sum[0] = n_g_leaves;
    l_sum[1] = n_g_spill_leaves;
    l_sum[2] = n_g_linked_boxes;

    l_min[0] = g_min_linked_boxes;
    l_max[0] = s.max_level_reached;
    l_max[1] = g_max_linked_boxes;

    MPI_Allreduce(l_sum, g_sum, 3, FVMC_MPI_GNUM, MPI_SUM, bt->comm);
    MPI_Allreduce(l_min, g_min, 1, FVMC_MPI_GNUM, MPI_MIN, bt->comm);
    MPI_Allreduce(l_max, g_max, 2, FVMC_MPI_GNUM, MPI_MAX, bt->comm);

    n_g_leaves = l_sum[0];
    n_g_spill_leaves = l_sum[1];
    n_g_linked_boxes = l_sum[2];

    g_min_linked_boxes = g_min[0];
    g_max_level_reached = g_max[0];
    g_max_linked_boxes = g_max[1];
  }

#endif

  /* Redefine final statistics */

  mean_linked_boxes = (double)n_g_linked_boxes / (double)n_g_leaves;
  box_ratio = (double)n_g_linked_boxes / (double)n_g_boxes;

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

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

  bftc_printf("\n"
             "Box tree statistics:\n\n");
  bftc_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (final/init):     %d\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (int)(bt->max_box_ratio), (int)(bt->max_level));

  bftc_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10lu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10lu\n"
             "  Initial number of boxes:            %10lu\n"
             "  Number of linked boxes:             %10lu\n"
             "  Mean number of leaves per box:      %10.4g\n\n",
             g_max_level_reached, (unsigned long)n_g_leaves,
             (unsigned long)n_g_spill_leaves, (unsigned long)n_g_boxes,
             (unsigned long)n_g_linked_boxes,  box_ratio);

  bftc_printf("Number of linked boxes per box tree leaf:\n"
             "  Mean value:         %10.4g\n"
             "  min. value:         %10lu\n"
             "  max. value:         %10lu\n\n",
             mean_linked_boxes,
             (unsigned long)(s.min_linked_boxes),
             (unsigned long)(s.max_linked_boxes));

  if (delta > 0) { /* Number of elements in each subdivision */

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bftc_printf("    %3d : [ %10u; %10u [ = %10u\n",
                 i+1,
                 (unsigned long)(g_min_linked_boxes + i*step),
                 (unsigned long)(g_min_linked_boxes + j*step),
                 (unsigned long)(count[i]));

    bftc_printf("    %3d : [ %10u; %10u ] = %10u\n",
               n_steps,
               (unsigned long)(g_min_linked_boxes + (n_steps - 1)*step),
               (unsigned long)(g_max_linked_boxes),
               (unsigned long)(count[n_steps - 1]));

  }
}

/*----------------------------------------------------------------------------
 * Dump an fvmc_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
fvmc_box_tree_dump(fvmc_box_tree_t  *bt)
{
  fvmc_box_tree_stats_t s;

  if (bt == NULL) {
    bftc_printf("\nBox tree: nil\n");
    return;
  }

  bftc_printf("\nBox tree: %p\n\n", bt);

  bftc_printf("  n_max_nodes:  %d\n\n"
             "  n_nodes:      %d\n",
             (int)(bt->n_max_nodes), (int)(bt->n_nodes));

  s = bt->stats;

  /* Print statistics and bounding boxes histogram */

  bftc_printf("  Number of children per leaf:              %d\n"
             "  Max number of bounding boxes for a leaf:  %d\n"
             "  Max value for box ratio (linked/init):    %d\n"
             "  Max level allowed:                        %d\n\n",
             bt->n_children, (int)(bt->threshold),
             (int)(bt->max_box_ratio), (int)(bt->max_level));

  bftc_printf("  Max level reached:                  %5u\n"
             "  Number of leaves:                   %10lu\n"
             "  Leaves with n_boxes > max_n_boxes:  %10lu\n"
             "  Initial number of boxes:            %10lu\n"
             "  Number of linked boxes:             %10lu\n",
             s.max_level_reached,
             (unsigned long)(s.n_leaves),
             (unsigned long)(s.n_spill_leaves),
             (unsigned long)(s.n_boxes),
             (unsigned long)(s.n_linked_boxes));

  bftc_printf("Bounding boxes related to each leaf of the box tree.\n"
             "  min. value:         %10lu\n"
             "  max. value:         %10lu\n\n",
             (unsigned long)(s.min_linked_boxes),
             (unsigned long)(s.max_linked_boxes));

  _dump_node(bt, 0);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
