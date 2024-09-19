/*
 * \file
 */

#ifndef __PDM_DOMAIN_UTILS_H__
#define __PDM_DOMAIN_UTILS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Shift all ln_to_gn of all current partition by the last max last domain
 *
 * \param [in]     n_domain             Number of domain
 * \param [in]     n_part               Number of partition for each domain       (size : n_domain )
 * \param [in]     pn_entity            Number of entity by domain and partition  (size : n_domaine, n_part)
 * \param [inout]  pentity_ln_to_gn     Number of entity by domain and partition  (size : n_domaine, n_part)
 * \param [in]     comm                 MPI Communicator
 *
 * \returnd        Array of shift (size : n_domain+1)
 */
PDM_g_num_t*
PDM_compute_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_MPI_Comm     comm
);


/**
 * \brief Shift all ln_to_gn of all current partition by the last max last domain
 *
 * \param [in]     n_domain             Number of domain
 * \param [in]     n_part               Number of partition for each domain       (size : n_domain )
 * \param [in]     pn_entity            Number of entity by domain and partition  (size : n_domaine, n_part)
 * \param [inout]  pentity_ln_to_gn     Number of entity by domain and partition  (size : n_domaine, n_part)
 * \param [in]     shift_by_domain      For each domain shift to apply (can be compute by PDM_compute_offset_ln_to_gn_by_domain)
 * \param [in]     sens                 Manage sens of shift (1 or - 1)
 *
 */
void
PDM_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_g_num_t     *shift_by_domain,
  int              sens
);





#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOMAIN_UTILS_H__ */
