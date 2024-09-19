#ifndef __PDM_DIST_CLOUD_SURF_PRIV_H__
#define __PDM_DIST_CLOUD_SURF_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_surf_mesh.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_timer.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

#define NTIMER 8

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  UPPER_BOUND_DIST              = 1,
  CANDIDATE_SELECTION           = 2,
  LOAD_BALANCING_ELEM_DIST      = 3,
  COMPUTE_ELEM_DIST             = 4,
  RESULT_TRANSMISSION           = 5,
  END                           = 6,
  BBTREE_CREATE                 = 7,

} _ol_timer_step_t;


/**
 * \struct PDM_dist_cloud_surf_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  double      **dist;
  double      **proj;
  PDM_g_num_t **closest_elt_gnum;

} _points_cloud_t;

/**
 * \struct PDM_dist_cloud_surf_t
 * \brief  Distance to a mesh surface structure
 *
 */

struct _pdm_dist_cloud_surf_t {

  PDM_MPI_Comm      comm;                    /*!< MPI communicator */
  PDM_ownership_t   owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t        results_is_getted;       /*!< Flags to indicate if result is getted      */
  int               n_point_cloud;           /*!< Number of point clouds */

  PDM_mesh_nature_t mesh_nature;             /*!< Nature of the mesh */

  PDM_surf_mesh_t  *surf_mesh;               /*!< Surface mesh pointer */
  PDM_surf_mesh_t  *_surf_mesh;              /*!< Surface mesh pointer */

  // PDM_Mesh_nodal_t*  mesh_nodal;             /*!< Surface mesh identifier */
  PDM_part_mesh_nodal_t*  mesh_nodal;       /*!< Nodal mesh identifier */

  _points_cloud_t *points_cloud;             /*!< Point clouds */

  PDM_timer_t *timer;                        /*!< Timer */

  double times_elapsed[NTIMER];              /*!< Elapsed time */

  double times_cpu[NTIMER];                  /*!< CPU time */

  double times_cpu_u[NTIMER];                /*!< User CPU time */

  double times_cpu_s[NTIMER];                /*!< System CPU time */


} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/
#undef NTIMER
#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_DIST_CLOUD_SURF_PRIV_H__ */
