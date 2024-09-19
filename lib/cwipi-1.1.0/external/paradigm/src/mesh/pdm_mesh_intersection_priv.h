#ifndef __PDM_MESH_INTERSECTION_PRIV_H__
#define __PDM_MESH_INTERSECTION_PRIV_H__

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
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_to_part.h"
#include "pdm_extract_part.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 23

/*=============================================================================
 * Static global variables
 *============================================================================*/

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

} _mesh_intersection_timer_step_t;


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_ol_t
 * \brief  Overlay type
 *
 * _PDM_ol_t defines a overlaying structure
 *
 */
struct _pdm_mesh_intersection_t {

  double   project_coef;         /*!< Projection coefficient to define the overlay
                                      surface projection :
                                      If value == 0, the surface projection is MeshA
                                      If value == 1, the surface projection is MeshB
                                      If 0 < value < 1, the projection surface is an
                                      intermediate surface */
  double   vtx_car_length_tol;   /*!< Absolute tolerance used to define local geometric
                                      tolerance for vertex caracteristic lenght
                                      (tolerance > 0) */

  double   extents_tol;          /*!< Absolute tolerance used to define local geometric
                                      tolerance for vertex caracteristic lenght
                                      (tolerxance > 0) */

  double   same_plane_tol;       /*!< Absolute tolerance used to check if 2 surfaces
                                      are the same plane surface */

  double bbox_tolerance;


  PDM_MPI_Comm comm;             /*!< MPI communicator */

  PDM_mesh_intersection_kind_t intersect_kind;

  int              n_part_mesh[2];
  int              dim_mesh[2];
  PDM_part_mesh_t *mesh[2];

  PDM_part_mesh_nodal_t *mesh_nodal[2]; /*!< Mesh nodal A/B */


  PDM_extract_part_t *extrp_mesh[2];

  // _ol_mesh_t  *olMeshA;       /*!< Overlay Mesh A */
  // _ol_mesh_t  *olMeshB;       /*!< Overlay Mesh B */

  PDM_timer_t *timer;


  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


  /* Results */
  PDM_ownership_t owner;
  int tag_extrp_mesh;
  int tag_elt_a_elt_b_get;
  int tag_box_a_box_b_get;
  int tag_elt_b_elt_a_get;
  int tag_elt_volume_get[2];
  int         **elt_a_elt_b_idx;
  PDM_g_num_t **elt_a_elt_b;
  int          *box_a_box_b_idx;
  int          *box_a_box_b;
  double      **elt_a_elt_b_volume;

  double      **elt_b_elt_a_volume;
  double      **elt_volume[2];

  PDM_ownership_t     ptp_ownership;
  PDM_part_to_part_t *ptp;




  /* vol_vol */
  int     tetraisation_pt_type;
  double *tetraisation_pt_coord;

  /* debug */
  double local_vol_A_B;
  double global_vol_A_B;
  double global_vol_A;
};

/*=============================================================================
 * Static global variables
 *============================================================================*/

// static PDM_Handles_t *olArray = NULL; /*!< Array to storage overlay identifiers */

/*=============================================================================
 * Static function definitions
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_INTERSECTION_PRIV_H__ */
