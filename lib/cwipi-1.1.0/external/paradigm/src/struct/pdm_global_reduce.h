/*
 * \file
 */

#ifndef __PDM_GLOBAL_REDUCE_H__
#define __PDM_GLOBAL_REDUCE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_global_reduce_t PDM_global_reduce_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that computes a global reduction
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_global_reduce object
 */

PDM_global_reduce_t *
PDM_global_reduce_create
(
 const int          n_part,
 const PDM_MPI_Comm comm
);


/**
 *
 * \brief Free a global point reduce structure
 *
 * \param [in]   gre          Pointer to \ref PDM_global_reduce object
 *
 */

void
PDM_global_reduce_free
(
 PDM_global_reduce_t *gre
);


/**
 *
 * \brief Set global ids
 *
 * \param [in]   gre           Pointer to \ref PDM_global_reduce object
 * \param [in]   i_part        Current partition
 * \param [in]   n_pts         Number of points in the partition
 * \param [in]   pts_ln_to_gn  Global ids of points in the partition
 *
 */

void
PDM_global_reduce_g_num_set
(
 PDM_global_reduce_t *gre,
 const int            i_part,
 const int            n_pts,
 const PDM_g_num_t   *pts_ln_to_gn
 );


/**
 *
 * \brief Set reduction operation
 *
 * \param [in]   gre                       Pointer to \ref PDM_global_reduce object
 * \param [in]   operation                 Type of reduction operation
 */

void
PDM_global_reduce_operation_set
(
 PDM_global_reduce_t   *gre,
 const PDM_reduce_op_t  operation
);


/**
 *
 * \brief Set local field
 *
 * \param [in]   gre                       Pointer to \ref PDM_global_reduce object
 * \param [in]   i_part                    Current partition
 * \param [in]   stride                    Stride of the field
 * \param [in]   local_field               Local value of field
 *                                         (can be NULL for any other reduction operation)
 * \param [in]   global_reduced_field_ptr  Pointer where global reduced field
 *                                         will be stored after computing
 */

void
PDM_global_reduce_field_set
(
 PDM_global_reduce_t *gre,
 const int            i_part,
 const int            stride,
 const double        *local_field,
 double              *global_reduced_field_ptr
);


/**
 *
 * \brief Compute the global reduced field
 *
 * \param [in]   gre          Pointer to \ref PDM_global_reduce object
 *
 */

void
PDM_global_reduce_field_compute
(
 PDM_global_reduce_t *gre
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GLOBAL_REDUCE_H__ */
