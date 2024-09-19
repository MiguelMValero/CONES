/*
 * \file
 */

#ifndef PDM_POINT_TREE_SEQ_H
#define PDM_POINT_TREE_SEQ_H

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

typedef struct _pdm_point_tree_seq_t     PDM_point_tree_seq_t;
typedef struct _pdm_point_tree_seq_shm_t PDM_point_tree_seq_shm_t;

/**
 * \enum PDM_point_tree_seq_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_POINT_TREE_SEQ_NORTH_WEST_NADIR,
  PDM_POINT_TREE_SEQ_NORTH_WEST_ZENITH,
  PDM_POINT_TREE_SEQ_NORTH_EAST_NADIR,
  PDM_POINT_TREE_SEQ_NORTH_EAST_ZENITH,
  PDM_POINT_TREE_SEQ_SOUTH_WEST_NADIR,
  PDM_POINT_TREE_SEQ_SOUTH_WEST_ZENITH,
  PDM_POINT_TREE_SEQ_SOUTH_EAST_NADIR,
  PDM_POINT_TREE_SEQ_SOUTH_EAST_ZENITH,
} PDM_point_tree_seq_child_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a point_tree structure
 *
 * \param [in]   tree_type          Tree type
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_point_tree_seq object
 */

PDM_point_tree_seq_t *
PDM_point_tree_seq_create
(
 const PDM_doctree_local_tree_t tree_type,
 const int                      depth_max,
 const int                      points_in_leaf_max,
 const double                   tolerance
);


/**
 *
 * \brief Free a point_tree structure
 *
 * \param [in]   point_tree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_free
(
 PDM_point_tree_seq_t *ptree
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_pts             Number of points in cloud
 * \param [in]   pts_coord         Point coordinates
 *
 */

void
PDM_point_tree_seq_point_cloud_set
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_pts,
 const double               *pts_coord
);


/**
 *
 * \brief Build a point_tree
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_build
(
 PDM_point_tree_seq_t *ptree
);


void
PDM_point_tree_seq_build_from_boxes
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_box,
       double               *box_extents
);


/**
 *
 * \brief Write point_tree nodes in a VTK file
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   filename              Output file name
 *
 */

void PDM_point_tree_seq_write_nodes
(
       PDM_point_tree_seq_t *ptree,
 const char                 *filename
 );


/**
 *
 * \brief Get number of children per node in  point_tree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 *
 * \return   Number of children per node in point_tree
 */

int
PDM_point_tree_n_children_get
(
 PDM_point_tree_seq_t *ptree
 );


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_new_to_old_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **new_to_old
);


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  old_to_new             Old to new order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_old_to_new_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **old_to_new
);


/**
 *
 * \brief Get point coords in point_tree's order
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  pts_coord             Point coordinates
 *
 */

void
PDM_point_tree_seq_sorted_points_get
(
 PDM_point_tree_seq_t  *ptree,
 double               **pts_coord
);


/**
 *
 * \brief Get point range of a point_tree node
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   node_id               Current node ID (zero-based)
 * \param [out]  point_range           Point range of current node
 *
 * \return   Number of points inside current node
 *
 */

int
PDM_point_tree_seq_point_range_get
(
       PDM_point_tree_seq_t *ptree,
 const int                   node_id,
       int                  *point_range
);


/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  ptree                Pointer to \ref PDM_point_tree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_ids             IDs of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_point_tree_seq_extract_nodes
(
  PDM_point_tree_seq_t  *ptree,
  int                    root_id,
  int                    n_depth,
  int                   *n_node,
  int                  **node_ids,
  double               **node_extents,
  int                  **node_weight
);

/**
 *
 * \brief Get node extents of subtree of given a list of node_id
 *
 * \param [in]  ptree                Pointer to \ref PDM_point_tree_seq object
 * \param [in]  n_child_to_extract   Number of node to extract
 * \param [in]  child_ids_to_extract Id of node to extract
 * \param [out] n_extract_child      Number of child extracted
 * \param [out] extract_extents      Extents of all child
 * \param [out] extract_child_id     Id of all child
 * \param [out] extract_is_leaf      Tag to know if child is a leaf or not
 *
 */
void
PDM_point_tree_seq_extract_extents_by_child_ids
(
  PDM_point_tree_seq_t  *ptree,
  const int              n_node_to_extract,
  const int             *node_ids_to_extract,
        int             *n_extract_child,
        int            **node_to_child_idx,
        int            **extract_child_id,
        int            **extract_is_leaf,
        double         **extract_extents
);

/**
 *
 * \brief Look for closest points stored inside a point_tree
 *
 * \param [in]   ptree                  Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_kdtree_pt_id   Closest point in kdtree ID (zero-based)
 * \param [out]  closest_kdtree_pt_dist Closest point in kdtree distance
 *
 */

void
PDM_point_tree_seq_closest_point
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_pts,
       double               *pts,
       int                  *closest_kdtree_pt_id,
       double               *closest_kdtree_pt_dist2
 );


/**
 *
 * \brief Get points located inside a set of boxes
 *
 * \param [in]   ptree                  Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_box                  Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [out]  box_pts_idx            Index of points located in boxes
 * \param [out]  box_pts                Local ids of points located in boxes (zero-based)
 *
 */

void
PDM_point_tree_seq_points_inside_boxes
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_box,
 const double                 box_extents[],
       int                  **box_pts_idx,
       int                  **box_pts
);


/**
 *
 * \brief Look for points inside at set of balls
 *
 * \param [in]  ptree                Pointer to \ref PDM_point_tree_seq object
 * \param [in]  n_ball               Number of balls
 * \param [in]  ball_center          Center of balls (size = \ref n_ball * 3)
 * \param [in]  ball_radius2         Squared radius of balls (size = \ref n_ball)
 * \param [out] ball_pts_idx         Index for ball->points graph (size \ref n_ball + 1)
 * \param [out] ball_pts             Ball->points graph (zero-based IDs)
 * \param [out] ball_pts_dist2       Distance from points to ball centers
 *
 */

void
PDM_point_tree_seq_points_inside_balls
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_ball,
       double                *ball_center,
       double                *ball_radius2,
       int                  **ball_pts_idx,
       int                  **ball_pts,
       double               **ball_pts_dist2
);



PDM_point_tree_seq_shm_t *
PDM_point_tree_make_shared
(
  PDM_point_tree_seq_t *local_ptree,
  PDM_MPI_Comm          comm_shared
);

void
PDM_point_tree_seq_shm_free
(
 PDM_point_tree_seq_shm_t *shm_ptree
);


/**
 *
 * \brief Get point coords in point_tree's order
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  pts_coord             Point coordinates
 *
 */
void
PDM_point_tree_seq_shm_sorted_points_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 double                   **pts_coord
);

/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_shm_point_new_to_old_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 int                      **new_to_old
);

void
PDM_point_tree_seq_shm_point_old_to_new_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 int                      **old_to_new
);

/**
 *
 * \brief Get points located inside a set of boxes (search in shared-memory tree)
 *
 * \param [in]   ptree                  Pointer to \ref PDM_point_tree_seq object
 * \param [in]   i_shm_rank             Shared-memory rank to explore
 * \param [in]   n_box                  Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [out]  box_pts_idx            Index of points located in boxes
 * \param [out]  box_pts                Local ids of points located in boxes (zero-based)
 *
 */

void
PDM_point_tree_seq_points_inside_boxes_shared
(
       PDM_point_tree_seq_shm_t  *shm_ptree,
 const int                        i_shm_rank,
 const int                        n_box,
 const double                     box_extents[],
       int                      **box_pts_idx,
       int                      **box_pts
 );


void PDM_point_tree_seq_write_nodes_shared
(
       PDM_point_tree_seq_shm_t *shm_ptree,
 const int                       i_shm_rank,
 const char                     *filename
 );


void
PDM_point_tree_seq_closest_point_shared
(
       PDM_point_tree_seq_shm_t *shm_ptree,
 const int                       i_shm_rank,
 const int                       n_pts,
       double                   *pts_coord,
       int                      *closest_point_id,
       double                   *closest_point_dist2
 );

void
PDM_tree_intersection_point_box2
(
 PDM_point_tree_seq_t  *btree,
 PDM_point_tree_seq_t  *ptree,
 double                *box_extents,
 int                  **box_pts_idx,
 int                  **box_pts
 );


void
PDM_point_tree_seq_intersect_box_leaf
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_box,
 const double                 box_extents[],
       int                  **box_leaf_idx,
       int                  **box_leaf
 );

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_POINT_TREE_SEQ_H */

