#ifndef __PDM_GLOBAL_REDUCE_PRIV_H__
#define __PDM_GLOBAL_REDUCE_PRIV_H__

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_global_reduce_t
 * \brief  Define a global reduce
 *
 */

struct _pdm_global_reduce_t {

  PDM_MPI_Comm   comm;            /*!< MPI communicator */
  int            n_part;          /*!< Number of partitions */
  int           *n_elts;          /*!< Number of elements in partitions */
  PDM_g_num_t  **g_nums;          /*!< Global numbering of elements */

  PDM_reduce_op_t operation;      /*!< Type of reduction operation */
  double  **local_field;          /*!< Local field */
  double  **global_reduced_field; /*!< Global reduced field */

  PDM_part_to_block_t *ptb;       /*!< Part to block structure */
  PDM_block_to_part_t *btp;       /*!< Block to part structure */
  int                  stride;    /*!< Current Field stride */
  int                **strides;   /*!< Strides array storage
                                   *   (In ptb strides are variable) */
} ;

/*=============================================================================
 * Static global variables
 *============================================================================*/
#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_GLOBAL_REDUCE_PRIV_H__ */
