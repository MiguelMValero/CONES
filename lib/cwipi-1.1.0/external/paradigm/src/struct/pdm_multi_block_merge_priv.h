#ifndef __PDM_MULTI_BLOCK_MERGE_PRIV_H__
#define __PDM_MULTI_BLOCK_MERGE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_to_block.h"
#include "pdm_multi_block_to_part.h"

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

struct _pdm_multi_block_merge_t {

  PDM_MPI_Comm            comm;             /*!< MPI Communicator */

  int                     n_rank;
  int                     i_rank;

  int                     n_block;
  int                    *blocks_ini_dn;
  PDM_g_num_t            *multi_block_distrib;
  PDM_g_num_t            *distrib_merge;

  /*
   * Results
   */
  PDM_g_num_t            *old_distrib;
  int                    *dold_to_new_idx;
  PDM_g_num_t            *dold_to_new;

  int                     dn_new_block;

  PDM_multi_block_to_part_t* mbtp;

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTI_BLOCK_MERGE_PRIV_H__ */
