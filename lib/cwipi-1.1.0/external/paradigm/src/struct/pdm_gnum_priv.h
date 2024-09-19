#ifndef __PDM_GNUM_PRIV_H__
#define __PDM_GNUM_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
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

/**
 * \struct _pdm_gen_gnum_t
 * \brief  Define a global numberring
 *
 */

struct _pdm_gen_gnum_t {

  PDM_MPI_Comm     comm;              /*!< MPI communicator */
  PDM_ownership_t  owner;             /*!< Which have the responsabilities of results */
  PDM_bool_t       results_is_getted; /*!< Flags to indicate if result is getted      */
  int              n_part;            /*!< Number of partitions */
  int              dim;               /*!< Spatial dimension */
  int              nuplet;            /*!< Size  of tuple in parent */
  PDM_bool_t       merge;             /*!< Merge double point status */
  double           tolerance;         /*!< Geometric tolerance */
  PDM_g_num_t      n_g_elt;           /*!< Global number of elements */
  int              *n_elts;           /*!< Number of elements in partitions */
  PDM_g_num_t     **g_nums;           /*!< Global numbering of elements */
  double          **coords;           /*!< Coordinates of elements */
  double          **char_length;      /*!< Characteristic length */
  int             **index;            /*!< Index : used if merge is activated */
  PDM_g_num_t     **parent;           /*!< Global n */

} ;

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_GNUM_PRIV_H__ */
