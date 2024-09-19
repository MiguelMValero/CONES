/*
 * \file
 * \author equemera
 *
 * \date January 18, 2018, 9:48 AM
 */

#ifndef PDM_BLOCK_TO_BLOCK_H
#define	PDM_BLOCK_TO_BLOCK_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_block_to_block_t
 * \brief  Block-to-Block redistribution
 *
 */

typedef struct _pdm_block_to_block_t PDM_block_to_block_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a block-to-block redistribution
 *
 * \param [in]   block_distrib_ini_idx Origin block distribution (size : *n_rank* + 1)
 * \param [in]   block_distrib_end_idx Destination block distribution (size : *n_rank* + 1)
 * \param [in]   comm                  MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_block_t instance
 *
 */

PDM_block_to_block_t *
PDM_block_to_block_create
(
 PDM_g_num_t   *block_distrib_ini_idx,
 PDM_g_num_t   *block_distrib_end_idx,
 PDM_MPI_Comm   comm
);


/**
 *
 * \brief Exchange data from origin block to destination block
 *
 * \warning Variable stride is not yet available
 *
 * \param [in]   btb              Block-to-Block structure
 * \param [in]   s_data           Data size
 * \param [in]   t_stride         Stride type
 * \param [in]   cst_stride       Constant stride
 * \param [in]   block_stride_ini Origin block stride for each block element for \ref PDM_STRIDE_VAR
 *                                Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data_ini   Origin block data
 * \param [out]  block_stride_end Destination block stride
 * \param [out]  block_data_end   Destination block data
 *
 */
int
PDM_block_to_block_exch
(
 PDM_block_to_block_t  *btb,
 size_t                 s_data,
 PDM_stride_t           t_stride,
 int                    cst_stride,
 int                   *block_stride_ini,
 void                  *block_data_ini,
 int                   *block_stride_end,
 void                 **block_data_end
);



/**
 *
 * \brief Exchange data from origin block to destination block
 *
 * \warning Variable stride is not yet available
 *
 * \param [in]   btb              Block-to-Block structure
 * \param [in]   mpi_type         Mpi kind setup by all possible MPI_type_create
 * \param [in]   block_stride_ini Origin block stride for each block element for \ref PDM_STRIDE_VAR
 *                                Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data_ini   Origin block data
 * \param [out]  block_stride_end Destination block stride
 * \param [out]  block_data_end   Destination block data
 *
 */
int
PDM_block_to_block_exch_with_mpi_type
(
 PDM_block_to_block_t  *btb,
 PDM_stride_t           t_stride,
 PDM_MPI_Datatype       mpi_type,
 int                   *block_stride_ini,
 void                  *block_data_ini,
 int                   *block_stride_end,
 void                 **block_data_end
);


/**
 *
 * \brief Free a Block-to-Block structure
 *
 * \param [inout] btb  Block-to-Block structure
 *
 * \return       NULL pointer
 */

PDM_block_to_block_t *
PDM_block_to_block_free
(
 PDM_block_to_block_t *btb
);


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_BLOCK_H */

