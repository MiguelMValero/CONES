/*
 * \file
 */

#ifndef __PDM_BOX_TREE_H__
#define __PDM_BOX_TREE_H__

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

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_box.h"

#include "pdm_point_tree_seq.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _PDM_box_tree_data_t PDM_box_tree_data_t;
typedef struct _PDM_box_tree_t PDM_box_tree_t;

typedef enum {

  PDM_BOX_TREE_ASYNC_LEVEL,  /* Boxes are placed according to tree parameters,
                                and potentially at different levels */
  PDM_BOX_TREE_SYNC_LEVEL    /* All boxes are placed for all ranks at the
                                same level */

} PDM_box_tree_sync_t;

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
PDM_box_tree_data_create(void);


/**
 *
 * \brief Create a \ref PDM_box_tree_t structure and initialize it
 *
 * \param [in]   max_level      Max possible level
 * \param [in]   threshold      Max number of boxes linked to an octant if \ref max_level is not reached
 * \param [in]   max_box_ratio  Max n_linked_boxes / n_boxes ratio
 *
 * \return Pointer to an empty PDM_box_tree_t structure
 *
 */

PDM_box_tree_t *
PDM_box_tree_create(int    max_level,
                    int    threshold,
                    float  max_box_ratio);


/**
 *
 * \brief Destroy a \ref PDM_box_tree_data_t structure
 *
 * \param [in]   btd  Pointer to pointer to PDM_box_tree_data_t structure to destroy
 *
 */

void
PDM_box_tree_data_destroy(PDM_box_tree_data_t  *btd);


/**
 *
 * \brief Destroy a \ref PDM_box_tree_t structure
 *
 * \param [in]   bt   Pointer to pointer to PDM_box_tree_t structure to destroy
 *
 */

void
PDM_box_tree_destroy(PDM_box_tree_t  **bt);


/**
 *
 * \brief Get the deepest level allowed by the tree structure
 *
 * \param [in]   bt   Pointer to pointer to \ref PDM_box_tree_t structure to destroy
 *
 * \return  Deepest allowed level of the tree
 *
 */

int
PDM_box_tree_get_max_level(const PDM_box_tree_t  *bt);


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
PDM_box_tree_set_boxes(PDM_box_tree_t       *bt,
                       PDM_box_set_t        *boxes,
                       PDM_box_tree_sync_t   build_type);


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
                         const PDM_box_set_t   *boxes);


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
                                  int                  *box_l_num[]);


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
                                   PDM_g_num_t          *box_g_num[]);


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
                       size_t                mem_allocated[3]);


/**
 *
 * \brief Display local statistics about a \ref PDM_box_tree_t structure
 *
 * \param [in]  bt   Pointer to box tree structure
 *
 */

void
PDM_box_tree_dump_statistics(const PDM_box_tree_t  *bt);


/**
 *
 * \brief Dump a \ref PDM_box_tree_t structure
 *
 * \param [in]  bt   Pointer to box tree structure
 *
 */

void
PDM_box_tree_dump(PDM_box_tree_t  *bt);


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
 );


/**
 *
 * \brief Get an indexed list of all boxes within an upper bound distance
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_pts              Number of points
 * \param [in]   pts                Point coordinates
 * \param [in]   upper_bound_dist2  Squared upper bound distance (size = \ref n_pts)
 * \param [out]  i_boxes            Pointer to the index array on bounding boxes
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
 int             *i_boxes[],
 int             *boxes[]
 );


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
 * \param [out]  i_boxes            Pointer to the index array on bounding boxes
 * \param [out]  boxes              Pointer to the list of bounding boxes
 * \param [in]   d_opt              Normalization vector (or NULL)
 *
 */

void
PDM_box_tree_closest_upper_bound_dist_boxes_get_v2
(
 PDM_box_tree_t  *bt,
 const int        i_rank,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *i_boxes[],
 int             *boxes[],
 const double    *d_opt
 );

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
);


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
);

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
);

/**
 *
 * \brief Copy the local box tree of some ranks on all other ranks
 *
 * \param [in]   bt                 Pointer to box tree structure
 * \param [in]   n_copied_ranks     Number of copied ranks
 * \param [in]   copied_ranks       List of copied ranks
 * \param [out]  rank_copy_num      Transpose of list of copied ranks (-1 for non-copied ranks) (size = n_rank)
 *
 */

void
PDM_box_tree_copy_to_ranks
(
 PDM_box_tree_t *bt,
 int            *n_copied_ranks,
 int            *copied_ranks,
 int            *rank_copy_num
);

/**
 *
 * \brief Make box_tree shared among nodes
 *
 * \param [in]   bt                 Pointer to box tree structure
 *
 */
void
PDM_box_tree_copy_to_shm
(
 PDM_box_tree_t *bt
);

void
PDM_box_tree_free_copies
(
  PDM_box_tree_t *bt
);

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
 );


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
 );


void
PDM_box_tree_boxes_containing_points_shared
(
 PDM_box_tree_t *bt,
 const int       i_shm,
 const int       n_pts,
 const double   *pts_coord,
 int           **box_idx,
 int           **box_l_num
);

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
 );


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
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 const int        n_line,
 const double    *line_coord,
 int            **box_idx,
 int            **box_l_num
 );

void
PDM_box_tree_intersect_lines_boxes_shared
(
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 const int        n_line,
 const double    *line_coord,
 int            **box_idx,
 int            **box_l_num
 );

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
 * \param [out]  box_line_l_num Pointer to the list of lines intersecting boxes (size = \ref box_line_idx[\ref n_box])
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
 int           **box_line_l_num
 );

void
PDM_box_tree_intersect_boxes_lines_shared
(
 PDM_box_tree_t *bt,
 const int       i_shm,
 const int       n_line,
 const double   *line_coord,
 int           **box_line_idx,
 int           **box_line_l_num
 );


/**
 *
 * \brief Get an indexed list of all boxes intersecting boxes
 *
 * The search can be performed either in the local box tree (\ref i_copied_rank < 0) or in
 * any distant box tree copied locally from rank bt->copied_rank[\ref i_copied_rank]
 *
 * \param [in]   bt             Pointer to box tree structure
 * \param [in]   i_copied_rank  Copied rank
 * \param [in]   n_line         Number of boxes
 * \param [in]   line_coord     Boxes coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  box_idx        Pointer to the index array on boxes (size = \ref n_line + 1)
 * \param [out]  box_l_num      Pointer to the list of boxes intersecting boxes (size = \ref box_idx[\ref n_line])
 *
 */

void
PDM_box_tree_intersect_boxes_boxes
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_box,
 const double   *box_coord,
 int           **box_idx,
 int           **box_l_num
);

void
PDM_box_tree_intersect_boxes_boxes2
(
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       n_tgt_box,
 const double   *tgt_box_extents,
 int           **tbox_box_idx,
 int           **tbox_box
);

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
 );

void
PDM_box_tree_write_vtk
(
 const char     *filename,
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 const int       normalized
 );

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
);

void
PDM_box_tree_extract_leaves
(
 PDM_box_tree_t  *bt,
 int             *n_leaf,
 int            **leaf_id
 );

void
PDM_box_tree_extract_node_extents
(
 PDM_box_tree_t  *bt,
 int              n_node,
 int             *node_id,
 double          *node_extents,
 const int        normalized
 );

void
PDM_box_tree_extract_extents_by_child_ids
(
       PDM_box_tree_t  *bt,
 const int              normalized,
 const int              n_child_to_extract,
 const int             *child_ids_to_extract,
       int             *n_extract_boxes,
       double         **extract_extents,
       int             *n_extract_child,
       int            **extract_child_id,
       int            **extract_is_leaf
);

void
PDM_box_tree_assign_weight
(
 PDM_box_tree_t  *bt,
 const int        n_node,
 const int       *nodes_id,
       int       *weight
);


int
PDM_box_tree_get_box_ids
(
 PDM_box_tree_t  *bt,
 int              node_id,
 int            **box_ids
);

int
PDM_box_tree_box_extents_get
(
 PDM_box_tree_t  *bt,
 const int        i_copied_rank,
 double         **extents
 );


void
PDM_tree_intersection_point_box
(
 PDM_box_tree_t        *btree,
 PDM_point_tree_seq_t  *ptree,
 int                  **box_pts_idx,
 int                  **box_pts
 );

void
PDM_box_tree_write_vtk2
(
 const char     *filename,
 PDM_box_tree_t *bt,
 const int       i_copied_rank,
 double *s,
 double *d
 );
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_BOX_TREE_H__ */
