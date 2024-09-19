#ifndef __PDM_POINTS_MERGE_PRIV_H__
#define __PDM_POINTS_MERGE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _point_merge_t
 * \brief  Define a point merge structures
 *
 */


struct _pdm_points_merge_t {

  PDM_MPI_Comm     comm;                    /*!< MPI communicator                           */
  PDM_ownership_t  owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t       results_is_getted;       /*!< Flags to indicate if result is getted      */
  double           tolerance;               /*!< Relative geometric tolerance               */
  int              n_point_clouds;          /*!< Number of point cloud                      */
  int              depth_max;               /*!< Maximum depth of internal octrees          */
  int              points_in_leaf_max;      /*!< Maximum number of point in a leaf
                                             *   of internal octrees                        */
  int              *n_points;                /*!< Number of points in each cloud             */
  int               max_n_points;            /*!< Maximum number of points in each cloud     */
  const double    **point_clouds;            /*!< points cloud                               */
  const double    **char_length;             /*!< Characteristic length of points (optionnal)*/
  PDM_octree_t     *octree;                  /*!< Octree pointer                          */
  int             **candidates_idx;          /*!< Candidates indexes for each cloud          */
  int             **candidates_desc;         /*!< Candidates description for each cloud      */

} ;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINTS_MERGE_PRIV_H__ */
