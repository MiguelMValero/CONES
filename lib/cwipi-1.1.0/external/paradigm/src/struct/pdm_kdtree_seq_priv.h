#ifndef __PDM_KDTREE_SEQ_PRIV_H__
#define __PDM_KDTREE_SEQ_PRIV_H__

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

  int                    *ancestor_id;          /*!< Ids of ancestor in kdtree array */
  int                    *is_leaf;              /*!< IS a leaf >*/
  int                    *depth;                /*!< Depth in the tree */
  int                    *children_id;          /*!< Ids of children in kdtree array */
  int                    *range;                /*!< Ids of children in kdtree array */
  int                    *idx;                  /*!< Start index of point list for each node */
  int                    *n_points;             /*!< Number of points in node*/
  double                 *extents;              /*!< Extents of the node */

} _l_nodes_t;


/**
 * \struct _kdtree_seq_t
 * \brief  Define a kdtree
 *
 */

struct _pdm_kdtree_seq_t{

  double         extents[6];            /*!< Extents of current process */
  int            depth_max;             /*!< Maximum depth of the three */
  int            points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double         tolerance;             /*!< Relative geometric tolerance */
  int            n_nodes;               /*!< Current number of nodes in kdtree */
  int            n_nodes_max;           /*!< Maximum number of nodes in kdtree */

  int            n_pts;
  const double  *pts_coord;
  double        *_pts_coord;
  int           *new_to_old;  /*< Zero-based numbering */
  int           *old_to_new;  /*< Zero-based numbering */

  _l_nodes_t    *nodes;

};


struct _pdm_kdtree_seq_shm_t {
  PDM_MPI_Comm       comm_shared;
  PDM_kdtree_seq_t  *kdtrees;

  PDM_mpi_win_shared_t *w_is_leaf;
  PDM_mpi_win_shared_t *w_children_id;
  PDM_mpi_win_shared_t *w_range;
  PDM_mpi_win_shared_t *w_n_points;
  PDM_mpi_win_shared_t *w_extents;

  PDM_mpi_win_shared_t *w_point_ids;
  PDM_mpi_win_shared_t *w_point_clouds;

  // PDM_mpi_win_shared_t *w_idx;
  // PDM_mpi_win_shared_t *w_ancestor_id;
  // PDM_mpi_win_shared_t *w_depth;
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_KDTREE_SEQ_PRIV_H__ */
