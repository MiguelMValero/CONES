/*
 * \file
 */

#ifndef __PDM_PART_DISTANT_NEIGHBOR_H__
#define __PDM_PART_DISTANT_NEIGHBOR_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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
 * Type definitions
 *============================================================================*/

typedef struct _pdm_distant_neighbor_t PDM_distant_neighbor_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an initialized \ref _distant_neighbor_t structure
 *
 * This function returns an initialized \ref _distant_neighbor_t structure
 *
 * \param [in]   comm          MPI communicator
 * \param [in]   n_part        Number of partitions
 * \param [in]   n_entity      Number of entities for each partition
 * \param [out]  neighbor_idx  Indexes of candidate for each current part point
 *                              (size = number of entities in the current part + 1)
 * \param [out]  neighbor_desc Candidates description (process,
 *                                                     part in the process,
 *                                                     entitiy in the part)
 *
 * \return      A new initialized \ref PDM_distant_neighbor structure
 *
 */

PDM_distant_neighbor_t*
PDM_distant_neighbor_create
(
const PDM_MPI_Comm   comm,
const int            n_part,
const int           *n_entity,
      int          **neighbor_idx,
      int          **neighbor_desc
);


/**
 *
 * \brief Free an distant negihtbor structure
 *
 * \param [in]   id                 Identifier
 *
 */
void
PDM_distant_neighbor_free
(
 PDM_distant_neighbor_t* dn
);

/**
 * \brief Exchange data between \ref _distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *
 *
 */
void
PDM_distant_neighbor_exch
(
 PDM_distant_neighbor_t   *dn,
 size_t                    s_data,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                     **send_entity_stride,
 void                    **send_entity_data,
 int                    ***recv_entity_stride,
 void                   ***recv_entity_data
);

/**
 * \brief Exchange data between \ref _distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *
 *
 */
void
PDM_distant_neighbor_exch_int
(
 PDM_distant_neighbor_t    *dn,
 size_t                    s_data,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                     **send_entity_stride,
 int                     **send_entity_data,
 int                    ***recv_entity_stride,
 int                    ***recv_entity_data
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_DISTANT_NEIGHBOR_H__ */
