#ifndef __PDM_POINT_TREE_SEQ_PRIV_H__
#define __PDM_POINT_TREE_SEQ_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */


/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  int                        *ancestor_id;          /*!< Ids of ancestor in point tree array */
  int                        *is_leaf;              /*!< IS a leaf >*/
  PDM_point_tree_seq_child_t *location_in_ancestor; /*!< Location in ancestor */
  int                        *depth;                /*!< Depth in the tree */
  int                        *children_id;          /*!< Ids of children in point tree array */
  int                        *range;                /*!< Ids of children in point tree array */
  int                        *idx;                  /*!< Start index of point list for each node */
  int                        *n_points;             /*!< Number of points in node*/
  double                     *extents;              /*!< Extents of the node */

} _l_nodes_t;


/**
 * \struct _point_tree_seq_t
 * \brief  Define a point tree
 *
 */

struct _pdm_point_tree_seq_t{

  PDM_doctree_local_tree_t tree_type;

  double         extents[6];         /*!< Extents of current process */
  int            depth_max;          /*!< Maximum depth of the three */
  int            points_in_leaf_max; /*!< Maximum number of points in a leaf */
  double         tolerance;          /*!< Relative geometric tolerance */
  int            n_nodes;            /*!< Current number of nodes in point_tree */
  int            n_nodes_max;        /*!< Maximum number of nodes in point_tree */

  int            n_pts;
  const double  *pts_coord;
  double        *_pts_coord;
  int           *new_to_old;  /*< Zero-based numbering */
  int           *old_to_new;  /*< Zero-based numbering */

  _l_nodes_t    *nodes;

  int n_leaf;
  int n_leaf_max;
  int n_leaf_box_max;
  int *leaf_box_ids;
  int *leaf_box_idx;
  int *leaf_ids;

};


struct _pdm_point_tree_seq_shm_t {

  PDM_MPI_Comm          comm_shared;

  PDM_doctree_local_tree_t tree_type;

  PDM_point_tree_seq_t *ptrees;

  int                  *shared_pts_idx;
  int                  *shared_nodes_idx;

  PDM_mpi_win_shared_t *w_is_leaf;
  PDM_mpi_win_shared_t *w_children_id;
  PDM_mpi_win_shared_t *w_range;
  PDM_mpi_win_shared_t *w_n_points;
  PDM_mpi_win_shared_t *w_extents;

  PDM_mpi_win_shared_t *w_pts_coord;
  PDM_mpi_win_shared_t *w_new_to_old;
  PDM_mpi_win_shared_t *w_old_to_new;

  int                  *shm_n_nodes;     /*!< Number of nodes shm from other ranks */
  int                 **shm_is_leaf;
  int                 **shm_children_id;
  int                 **shm_range;
  int                 **shm_n_points;
  double              **shm_extents;
  int                  *shm_n_pts;     /*!< Number of points shm from other ranks */
  double              **shm_pts_coord;
  int                 **shm_new_to_old;
  int                 **shm_old_to_new;


  // PDM_mpi_win_shared_t *w_idx;
  // PDM_mpi_win_shared_t *w_ancestor_id;
  // PDM_mpi_win_shared_t *w_depth;
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_TREE_SEQ_PRIV_H__ */
