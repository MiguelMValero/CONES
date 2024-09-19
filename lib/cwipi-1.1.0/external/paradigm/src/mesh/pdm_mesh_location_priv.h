#ifndef __PDM_MESH_LOCATION_PRIV_H__
#define __PDM_MESH_LOCATION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh_nodal.h"
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

#define NTIMER_MESH_LOCATION 11

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;

  int          *n_located;
  int          *n_un_located;

  int         **located;
  int         **un_located;

  PDM_g_num_t **location;
  double      **dist2;
  double      **uvw;
  int         **weights_idx;
  double      **weights; /*!< Barycentric coordinates */
  double      **projected_coords;

} _point_cloud_t;


/**
 * \struct _points_in_element_t
 * \brief
 *
 */

typedef struct {

  int           n_part;
  int          *n_elts; // Aredescendre dans la structure parente
  int         **pts_inside_idx;
  PDM_g_num_t **gnum;
  double      **coords;
  double      **uvw;
  double      **projected_coords;
  int         **weights_idx;
  double      **weights;
  double      **dist2;

} _points_in_element_t;



/**
 * \struct _pdm_mesh_location_t
 * \brief  Structure to locate point clouds inside a mesh
 *
 */
struct _pdm_mesh_location_t {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */

  int mesh_dimension;

  int  shared_nodal;   /*!< 1 if mesh nodal is shared, 0 otherwise */
  PDM_part_mesh_nodal_t*  mesh_nodal;  /*!< Mesh identifier */
  PDM_part_mesh_nodal_t* _mesh_nodal;
  PDM_l_num_t **face_vtx_n; /* Mandatory to build mesh nodal */
  PDM_l_num_t **cell_face_n; /* Mandatory to build mesh nodal */
  PDM_l_num_t **cell_vtx_idx;
  PDM_l_num_t **cell_vtx;

  _point_cloud_t *point_clouds; /*!< Point clouds */

  int        use_user_extract;
  int      **is_elmt_select_by_user;

  double tolerance;

  PDM_mesh_location_method_t method;

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER_MESH_LOCATION]; /*!< Elapsed time */

  double times_cpu[NTIMER_MESH_LOCATION];     /*!< CPU time */

  double times_cpu_u[NTIMER_MESH_LOCATION];  /*!< User CPU time */

  double times_cpu_s[NTIMER_MESH_LOCATION];  /*!< System CPU time */

  _points_in_element_t *points_in_elements;

  int  reverse_result; /*!< Enable reverse results */

  PDM_ownership_t owner;       /*!< Ownership */
  int  tag_unlocated_get;      /*!< Tag call to unlocated_get function */ 
  int  tag_located_get;        /*!< Tag call to located_get function */ 
  int  tag_point_location_get; /*!< Tag call to point_location_get function */ 
  int  tag_points_in_elt_get;  /*!< Tag call to points_in_elt_get function */ 
  int  tag_cell_vtx_get;       /*!< Tag call to cell_vtx_get function */ 
  // int *tag_pts_in_elt_get;     /*!< Tag call to points_in_elt_get function */


  PDM_part_to_part_t **ptp; /*!< To exchange data between elt and points (both in user frame) */
  PDM_ownership_t     *ptp_ownership;
} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/
#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_MESH_LOCATION_PRIV_H__ */
