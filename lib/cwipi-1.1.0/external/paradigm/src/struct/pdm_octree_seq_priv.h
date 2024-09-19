#ifndef __PDM_OCTREE_SEQ_PRIV_H__
#define __PDM_OCTREE_SEQ_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_timer.h"
#include "pdm_overlay.h"
#include "pdm_error.h"
#include "pdm_printf.h"

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


/**
 * \struct _octant_t
 * \brief  Define an octant
 *
 */

//typedef struct  {
//
//  int                     ancestor_id;          /*!< Ids of ancestor in octree array */
//  int                     is_leaf;              /*!< IS a leaf >*/
//  PDM_octree_seq_child_t  location_in_ancestor; /*!< Location in ancestor */
//  int                     depth;                /*!< Depth in the tree */
//  int                     children_id[8];       /*!< Ids of children in octree array */
//  int                     range[2];             /*!< Ids of children in octree array */
//  int                     idx[9];               /*!< Start index of point list for each octant */
//  int                     n_points;             /*!< Number of points in octant*/
//  double                  extents[6];           /*!< Extents of the node */
//
//} _octant_t;



typedef struct {

  int                    *ancestor_id;          /*!< Ids of ancestor in octree array */
  int                    *is_leaf;              /*!< IS a leaf >*/
  PDM_octree_seq_child_t *location_in_ancestor; /*!< Location in ancestor */
  int                    *depth;                /*!< Depth in the tree */
  int                    *children_id;          /*!< Ids of children in octree array */
  int                    *range;                /*!< Ids of children in octree array */
  int                    *idx;                  /*!< Start index of point list for each octant */
  int                    *n_points;             /*!< Number of points in octant*/
  double                 *extents;              /*!< Extents of the node */

} _l_nodes_t;


/**
 * \struct _octree_seq_t
 * \brief  Define an octree
 *
 */

struct _pdm_octree_seq_t {

  double         extents[6];            /*!< Extents of current process */
  int            depth_max;             /*!< Maximum depth of the three */
  int            points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double         tolerance;             /*!< Relative geometric tolerance */
  int            n_nodes;               /*!< Current number of nodes in octree */
  int            n_nodes_max;           /*!< Maximum number of nodes in octree */
  int           *n_points;              /*!< Number of points in each cloud */
  int            t_n_points;            /*!< total number of points */
  int            n_point_clouds;        /*!< Number of point cloud */
  const double **point_clouds;          /*!< points cloud */
  int           *point_ids;             /*!< Id's of points in it cloud sorted by octree (size: n_points + 1) */
  int           *point_icloud;          /*!< Cloud's of points sorted by octree (size: n_points + 1) */

  _l_nodes_t    *nodes;

} ;


struct _pdm_octree_seq_shm_t {
  PDM_MPI_Comm       comm_shared;
  PDM_octree_seq_t  *octrees;

  PDM_mpi_win_shared_t *w_is_leaf;
  PDM_mpi_win_shared_t *w_children_id;
  PDM_mpi_win_shared_t *w_range;
  PDM_mpi_win_shared_t *w_n_points;
  PDM_mpi_win_shared_t *w_extents;

  PDM_mpi_win_shared_t *w_point_ids;
  PDM_mpi_win_shared_t *w_point_clouds;

  // PDM_mpi_win_shared_t *w_idx;
  // PDM_mpi_win_shared_t *w_ancestor_id;
  // PDM_mpi_win_shared_t *w_location_in_ancestor;
  // PDM_mpi_win_shared_t *w_depth;

};



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_OCTREE_SEQ_PRIV_H__ */
