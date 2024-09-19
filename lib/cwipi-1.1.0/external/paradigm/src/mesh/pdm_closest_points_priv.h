#ifndef __PDM_CLOSEST_POINTS_PRIV_H__
#define __PDM_CLOSEST_POINTS_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_part_to_part.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

#define NTIMER 2

/**
 * \struct _tgt_point_cloud_t
 * \brief  Target point cloud structure
 *
 */

typedef struct {

  int           n_part;            /*!< Number of partition */
  int          *n_points;          /*!< Number of points of each partition */
  double      **coords;            /*!< Point coordinates points of each partition */
  PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */
  PDM_g_num_t **closest_src_gnum;  /*!< Global numbering of the n_closest source points
                                     for each point of each partition  */
  double      **closest_src_dist; /*!< Distance to the n_closest source points
                                    for each point of each partition  */

} _tgt_point_cloud_t;


/**
 * \struct _src_point_cloud_t
 * \brief  Src point cloud structure
 *
 */

typedef struct {

  int           n_part;            /*!< Number of partition */
  int          *n_points;          /*!< Number of points of each partition */
  double      **coords;            /*!< Point coordinates points of each partition */
  PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */

  int         **tgt_in_src_idx;    /*!< Reverse results */
  PDM_g_num_t **tgt_in_src;        /*!< Reverse results */
  double      **tgt_in_src_dist;   /*!< Reverse results */

} _src_point_cloud_t ;


/**
 * \struct _PDM_closest_t
 * \brief  Closest points structure
 *
 */

struct _pdm_closest_point_t {

  PDM_MPI_Comm    comm;                         /*!< MPI communicator */
  PDM_ownership_t owner;                        /*!< Which have the responsabilities of results */
  PDM_bool_t      results_is_getted;            /*!< Flags to indicate if result is getted      */
  PDM_bool_t      tgt_in_src_results_is_getted; /*!< Flags to indicate if result is getted      */
  PDM_bool_t      tgt_in_src_results_is_getted_d; /*!< Flags to indicate if result is getted      */

  int n_closest;                                /*!< Number of closest source points to find for each
                                                  target point  */

  _src_point_cloud_t *src_cloud;                /*!< Source point cloud */

  _tgt_point_cloud_t *tgt_cloud;                /*!< Target point cloud */

  PDM_timer_t *timer;                           /*!< Timer */

  double times_elapsed[NTIMER];                 /*!< Elapsed time */

  double times_cpu[NTIMER];                     /*!< CPU time */

  double times_cpu_u[NTIMER];                   /*!< User CPU time */

  double times_cpu_s[NTIMER];                   /*!< System CPU time */


  PDM_part_to_part_t *ptp; /*!< To exchange data between src and tgt point clouds (both in user frame) */
  PDM_ownership_t     ptp_ownership;


} ;


#undef NTIMER
#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_CLOSEST_POINTS_PRIV_H__ */
