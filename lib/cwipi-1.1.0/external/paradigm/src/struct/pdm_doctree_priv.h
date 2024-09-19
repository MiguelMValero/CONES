#ifndef __PDM_DOCTREE_PRIV_H__
#define __PDM_DOCTREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_point_tree_seq.h"
#include "pdm_part_to_block.h"

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

#define NTIMER_DOCTREE 14

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                                  = 0,
  REDISTRIBUTE_PTS_HILBERT               = 1,
  BUILD_COARSE_TREE_AND_EXTRACT          = 2,
  BUILD_BBOX_COARSE                      = 3,
  BBOX_COARSE_SOLICITATE                 = 4,
  EQUILIBRATE_WITH_SOLICITATON           = 5,
  EQUILIBRATE_WITH_SOLICITATON_TRANSFERT = 6,
  UPDATE_SOLICITATION_SEND               = 7,
  BUILD_LOCAL_TREE                       = 8,
  BUILD_SHARED_LOCAL_TREE                = 9,
  UPDATE_SOLICITATION_WAIT               = 10,
  LOCAL_SOLICITATE                       = 11,
  EQUILIBRATE_PB                         = 12,
  END                                    = 13

} _doctree_timer_step_t;

/*============================================================================
 * Type definitions
 *============================================================================*/

struct _pdm_doctree_t {

  PDM_MPI_Comm               comm;                       /*!< MPI communicator                          */
  int                        dim;                        /*!< Dimension                                 */
  int                        n_part_cloud;               /*!< Dimension                                 */
  PDM_doctree_local_tree_t   local_tree_kind;


  int                        coarse_depth_max;           /*!< global_octree depth_max                   */
  int                        coarse_points_in_leaf_max;  /*!< global_octree max pts in leaf             */

  int                        local_depth_max;
  int                        local_points_in_leaf_max;
  double                     local_tolerance;

  PDM_point_tree_seq_t      *coarse_tree;              /*! coarse tree to orient among procs        */
  PDM_point_tree_seq_t      *local_tree;               /*! Local tree                               */
  PDM_point_tree_seq_shm_t  *shmem_tree;               /*! Shared tree among cores in current nodes */

  PDM_MPI_Comm               comm_dist_graph;
  int                        n_degree_in;
  int*                       neighbor_in;

  PDM_MPI_Comm               comm_shared;

  /* Cloud - Just reference */
  int          *n_point_cloud;
  PDM_g_num_t **pts_g_num;
  double      **pts_coords;
  int         **pts_init_location;

  /* Solicitation */
  PDM_tree_solicitation_t    solicitation_kind;
  int                        n_part;
  int                       *n_entity;
  int                      **init_location_entity;
  PDM_g_num_t              **entity_gnum;
  double                   **entity_coords;

  /* Equilibrate results */
  PDM_ownership_t            ownership;
  PDM_part_to_block_t       *ptb_unit_op_equi;
  int                       *block_pts_in_box_n;
  PDM_g_num_t               *block_pts_in_box_g_num;
  double                    *block_pts_in_box_coord;

  /* Equilibrate results */
  int                        dn_pts_equi;
  int                       *equi_pts_init_location;
  int                       *equi_pts_init_location_idx;


  /* Misc */
  PDM_timer_t *timer; /*!< Timer */
  double times_elapsed[NTIMER_DOCTREE]; /*!< Elapsed time */
  double times_cpu    [NTIMER_DOCTREE]; /*!< CPU time */
  double times_cpu_u  [NTIMER_DOCTREE]; /*!< User CPU time */
  double times_cpu_s  [NTIMER_DOCTREE]; /*!< System CPU time */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOCTREE_PRIV_H__ */
