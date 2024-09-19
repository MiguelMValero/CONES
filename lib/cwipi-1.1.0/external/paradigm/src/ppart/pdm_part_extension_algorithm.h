/*
 * \file
 */

#ifndef __PDM_PART_EXTENSION_ALGORITHM_H__
#define __PDM_PART_EXTENSION_ALGORITHM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_to_part.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


void
PDM_part_extension_interface_by_entity1_to_interface_by_entity2
(
  PDM_part_domain_interface_t  *pdi,
  PDM_bound_type_t              entity1_bound,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_entity2,
  int                          *n_part,
  int                         **pn_entity1_in,
  PDM_g_num_t                ***pentity1_ln_to_gn_in,
  int                        ***pentity1_hint_in,
  int                         **pn_entity2_in,
  PDM_g_num_t                ***pentity2_ln_to_gn_in,
  int                        ***pentity2_entity1_idx_in,
  int                        ***pentity2_entity1_in,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                  comm
);



void
PDM_part_extension_pconnectivity_to_extented_pconnectivity
(
  PDM_part_domain_interface_t    *pdi,
  PDM_bound_type_t                entity2_bound,
  int                             n_domain,
  PDM_g_num_t                    *shift_by_domain_entity2,
  int                            *n_part,
  int                           **pn_entity1_in,
  PDM_g_num_t                  ***pentity1_ln_to_gn_in,
  int                           **pn_entity2_in,
  PDM_g_num_t                  ***pentity2_ln_to_gn_in,
  int                          ***pentity1_entity2_idx_in,
  int                          ***pentity1_entity2_in,
  int                            *pn_entity1_extented,
  PDM_g_num_t                   **pentity1_extented_ln_to_gn,
  int                           **pentity1_extented_to_pentity1_idx,
  int                           **pentity1_extented_to_pentity1_triplet,
  int                           **pentity1_extented_to_pentity1_interface,
  int                           **pn_entity2_extented_out,
  PDM_g_num_t                  ***pentity2_extented_ln_to_gn_out,
  int                          ***pextented_entity1_entity2_idx_out,
  int                          ***pextented_entity1_entity2_out,
  int                          ***pentity2_extented_to_pentity2_idx_out,
  int                          ***pentity2_extented_to_pentity2_triplet_out,
  int                          ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                    comm
);



// Prevoir un concatenate




#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_ALGORITHM_H__ */
