#ifndef __PDM_PART_DOMAIN_INTERFACE_PRIV_H__
#define __PDM_PART_DOMAIN_INTERFACE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct pdm_part_domain_interface_t
 * \brief  Connnectivity between domains
 *
 * \ref pdm_part_domain_interface_t defines an interface structure
 *
 */

struct _pdm_part_domain_interface_t {

  int                            n_interface;
  int                            n_domain;
  int                           *n_part;

  int                            interface_describe_by_face;
  int                            interface_describe_by_edge;
  int                            interface_describe_by_vtx;

  PDM_domain_interface_mult_t     multidomain_intrf;

  int                          ***interface_pn_face;
  PDM_g_num_t                 ****interface_face_ln_to_gn;
  int                         ****interface_sgn_face;
  int                         ****interface_sens_face;     // Comment la face opposé est rangé par rapport à moi 1 = même sens  | -1 :  sens opposé
  int                         ****interface_ids_face;      // (i_proc, i_part, i_face)
  int                         ****interface_ids_face_idx;  // Size = interface_pn_face
  int                         ****interface_dom_face;      // (i_dom_cur, i_dom_opp)

  int                          ***interface_pn_edge;
  PDM_g_num_t                 ****interface_edge_ln_to_gn;
  int                         ****interface_sgn_edge;
  int                         ****interface_sens_edge;
  int                         ****interface_ids_edge;      // (i_proc, i_part, i_edge)
  int                         ****interface_ids_edge_idx;  // Size = interface_pn_edge
  int                         ****interface_dom_edge;      // (i_dom_cur, i_dom_opp)
  // ptr is = (i_domain, i_part, i_interface)

  int                          ***interface_pn_vtx;
  PDM_g_num_t                 ****interface_vtx_ln_to_gn;
  int                         ****interface_sgn_vtx;
  int                         ****interface_sens_vtx;    // Pas de sens pour les sommets always = 1
  int                         ****interface_ids_vtx;     // (i_proc, i_part, i_vtx)
  int                         ****interface_ids_vtx_idx; // Size = interface_pn_vtx
  int                         ****interface_dom_vtx;     // (i_dom_cur, i_dom_opp)

  double                      **translation_vect;
  double                      **rotation_direction;
  double                      **rotation_center;
  double                       *rotation_angle;

  PDM_ownership_t ownership;
  int is_result[PDM_BOUND_TYPE_MAX];

  PDM_MPI_Comm    comm;
};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_DOMAIN_INTERFACE_PRIV_H__ */
