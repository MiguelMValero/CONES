/*
 * \file
 */

#ifndef __PDM_MULTI_BLOCK_MERGE_H__
#define __PDM_MULTI_BLOCK_MERGE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_multi_block_merge_t      PDM_multi_block_merge_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a redistribute structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_multi_block_merge_t*
PDM_multi_block_merge_create
(
       PDM_g_num_t  **block_distrib_idx,
 const int            n_block,
       int           *n_selected,
       PDM_g_num_t  **selected_g_num,
       int            graph_size,
       int           *dmerge_idx,
       int           *dmerge_block_id,
       PDM_g_num_t   *dmerge_g_num,
       PDM_MPI_Comm   comm
);


int
PDM_multi_block_merge_get_n_block
(
 PDM_multi_block_merge_t* mbm
);


PDM_g_num_t*
PDM_multi_block_merge_get_distrib
(
 PDM_multi_block_merge_t* mbm
);

PDM_g_num_t*
PDM_multi_block_merge_get_multi_distrib
(
 PDM_multi_block_merge_t* mbm
);


void
PDM_multi_block_merge_exch
(
 PDM_multi_block_merge_t     *mbm,
 size_t                       s_data,
 PDM_stride_t                 t_stride,
 int                        **block_stride,
 void                       **block_data,
 int                        **merge_block_stride,
 void                       **merge_block_data
);


void
PDM_multi_block_merge_get_old_to_new
(
 PDM_multi_block_merge_t  *mbm,
 PDM_g_num_t             **old_distrib,
 int                     **dold_to_new_idx,
 PDM_g_num_t             **dold_to_new
);

void
PDM_multi_block_merge_exch_and_update
(
 PDM_multi_block_merge_t  *mbm_cur,
 PDM_multi_block_merge_t  *mbm_for_update,
 PDM_stride_t              t_stride,
 int                     **block_stride,
 PDM_g_num_t             **block_data,
 int                     **update_domain,
 int                     **merge_block_stride,
 PDM_g_num_t             **merge_block_data
);

// void
// PDM_multi_block_merge_reorder_block_hilbert
// (
//  PDM_multi_block_merge_t  *mbm,
//  double                   *block_coord,
//  double                   *padd_coord
// );

// PDM_g_num_t*
// PDM_multi_block_merge_compute_old_to_new
// (
//  PDM_multi_block_merge_t  *mbm
// );

// void
// PDM_multi_block_merge_get_old_to_new
// (
//  PDM_multi_block_merge_t  *mbm,
//  int                **dold_to_new_idx,
//  PDM_g_num_t        **dold_to_new
// );

// int
// PDM_multi_block_merge_get_n_block
// (
//  PDM_multi_block_merge_t  *mbm
// );

// PDM_g_num_t*
// PDM_multi_block_merge_get_parent_blk_g_num
// (
//  PDM_multi_block_merge_t  *mbm
// );

/**
 * \brief Create a redistribute structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
void
PDM_multi_block_merge_free
(
 PDM_multi_block_merge_t* mbm
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTI_BLOCK_MERGE_H__ */
