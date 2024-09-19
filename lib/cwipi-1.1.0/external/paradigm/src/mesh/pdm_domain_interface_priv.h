#ifndef __PDM_DOMAIN_INTERFACE_PRIV_H__
#define __PDM_DOMAIN_INTERFACE_PRIV_H__

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
 * \struct pdm_domain_interface_t
 * \brief  Connnectivity between domains
 *
 * \ref pdm_domain_interface_t defines an interface structure
 *
 */

struct _pdm_domain_interface_t{

  int                           n_interface;
  int                           n_domain;
  PDM_domain_interface_mult_t   multidomain_intrf;
  int                          *interface_dn_face;
  PDM_g_num_t                 **interface_ids_face;
  int                         **interface_dom_face;

  int                          *interface_dn_edge;
  PDM_g_num_t                 **interface_ids_edge;
  int                         **interface_dom_edge;

  int                          *interface_dn_vtx;
  PDM_g_num_t                 **interface_ids_vtx;
  int                         **interface_dom_vtx;

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

#endif /* __PDM_DOMAIN_INTERFACE_PRIV_H__ */
