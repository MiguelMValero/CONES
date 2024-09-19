/*
 * \file
 */

#ifndef PDM_KDTREE_SEQ_H
#define PDM_KDTREE_SEQ_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_kdtree_seq_t     PDM_kdtree_seq_t;
typedef struct _pdm_kdtree_seq_shm_t PDM_kdtree_seq_shm_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a kdtree structure
 *
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_kdtree_seq object
 */

PDM_kdtree_seq_t *
PDM_kdtree_seq_create
(
 const int    depth_max,
 const int    points_in_leaf_max,
 const double tolerance
);


/**
 *
 * \brief Free a kdtree structure
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_free
(
 PDM_kdtree_seq_t *kdtree
);



/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 *
 */


void
PDM_kdtree_seq_point_cloud_set
(
 PDM_kdtree_seq_t *kdtree,
 const int         n_points,
 const double     *coords
);



/**
 *
 * \brief Build kdtree
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_build
(
 PDM_kdtree_seq_t *kdtree
);



/**
 *
 * \brief Write kdtree nodes in a VTK file
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [in]   filename               Output file name
 *
 */

void PDM_kdtree_seq_write_nodes
(
 PDM_kdtree_seq_t *kdtree,
 const char       *filename
 );


/**
 *
 * \brief Look for points inside at set of balls
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  n_ball               Number of balls
 * \param [in]  ball_center          Center of balls (size = \ref n_ball * 3)
 * \param [in]  ball_radius2         Squared radius of balls (size = \ref n_ball)
 * \param [out] ball_pts_idx         Index for ball->points graph (size \ref n_ball + 1)
 * \param [out] ball_pts_l_num       Ball->points graph (cloud_id, point_id)
 * \param [out] ball_pts_dist2       Distance from points to ball centers
 *
 */

void
PDM_kdtree_seq_points_inside_balls
(
 const PDM_kdtree_seq_t  *kdtree,
 const int                n_ball,
 double                  *ball_center,
 double                  *ball_radius2,
 int                    **ball_pts_idx,
 int                    **ball_pts_l_num,
 double                 **ball_pts_dist2
 );


/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_ids             IDs of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_kdtree_seq_extract_extent
(
  PDM_kdtree_seq_t  *kdtree,
  int                root_id,
  int                n_depth,
  int               *n_node,
  int              **node_ids,
  double           **node_extents,
  int              **node_weight
);


/**
 *
 * \brief Get points located inside a set of boxes
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [in]   n_box                  Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [out]  pts_idx                Index of points located in boxes
 * \param [out]  pts_l_num              Local ids of points located in boxes
 *
 */

void
PDM_kdtree_seq_points_inside_boxes
(
       PDM_kdtree_seq_t   *kdtree,
 const int                 n_box,
 const double              box_extents[],
       int               **pts_idx,
       int               **pts_l_num
);


/**
 *
 * \brief Look for closest points stored inside a kdtree
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_kdtree_pt_id   Closest point in kdtree index
 * \param [out]  closest_kdtree_pt_dist Closest point in kdtree distance
 *
 */

void
PDM_kdtree_seq_closest_point
(
PDM_kdtree_seq_t *kdtree,
const int         n_pts,
double           *pts,
int              *closest_kdtree_pt_id,
double           *closest_kdtree_pt_dist2
);


/**
 *
 * \brief Get point order in kdtree
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [out]  new_to_old             New to old order of points in kdtree
 *
 */

void
PDM_kdtree_seq_point_new_to_old_get
(
 PDM_kdtree_seq_t  *kdtree,
 int              **new_to_old
);


/**
 *
 * \brief Get point order in kdtree
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [out]  old_to_new             Old to new order of points in kdtree
 *
 */

void
PDM_kdtree_seq_point_old_to_new_get
(
 PDM_kdtree_seq_t  *kdtree,
 int              **old_to_new
);


void
PDM_kdtree_seq_sorted_points_get
(
 PDM_kdtree_seq_t  *kdtree,
 double           **pts_coord
);


int
PDM_kdtree_seq_point_range_get
(
       PDM_kdtree_seq_t  *kdtree,
 const int                node_id,
       int               *point_range
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_KDTREE_SEQ_H */

