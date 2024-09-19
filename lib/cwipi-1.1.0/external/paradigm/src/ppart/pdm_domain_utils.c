/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_domain_utils.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
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
)
{

  PDM_g_num_t *shift_by_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *shift_by_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        shift_by_domain_loc[i_domain] = PDM_MAX(shift_by_domain_loc[i_domain], _pentity_ln_to_gn[i]);
      }
    }
  }

  shift_by_domain[0] = 0;
  PDM_MPI_Allreduce(shift_by_domain_loc, &shift_by_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(shift_by_domain, n_domain+1);

  free(shift_by_domain_loc);

  return shift_by_domain;
}

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
)
{
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        _pentity_ln_to_gn[i] = _pentity_ln_to_gn[i] + sens * shift_by_domain[i_domain];
      }
    }
  }
}




#ifdef __cplusplus
}
#endif /* __cplusplus */

