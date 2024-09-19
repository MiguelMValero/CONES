/*
 * \file
 */

#ifndef PDM_DCONNECTIVITY_TRANSFORM_H_
#define PDM_DCONNECTIVITY_TRANSFORM_H_

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_block_to_part.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
);


/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       Distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       Distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx  Index of dentitiy1->dentity2 connectivity
 * \param [in]   dentity1_entity2      Connectivity of dentitiy1->dentity2
 * \param [in]   is_signed             If connectivity is signed
 * \param [in]   dentity2_entity1_idx  Index of dentitiy2->dentity1 connectivity
 * \param [in]   dentity2_entity1      Connectivity of dentitiy2->dentity1
 */
void
PDM_dconnectivity_transpose
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
       PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
       int              is_signed,
       int            **dentity2_entity1_idx,
       PDM_g_num_t    **dentity2_entity1
);

void PDM_dfacecell_to_dcellface
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const PDM_g_num_t* dface_cell,
  int**              dcell_face_idx,
  PDM_g_num_t**      dcell_face,
  PDM_MPI_Comm       comm
);
void PDM_dcellface_to_dfacecell
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const int*         dcell_face_idx,
  const PDM_g_num_t* dcell_face,
  PDM_g_num_t**      dface_cell,
  PDM_MPI_Comm       comm
);

void
PDM_deduce_combine_connectivity_dual
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       PDM_g_num_t    **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
);


void
PDM_dorder_reverse
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity_distrib,
 const PDM_g_num_t     *dentity1_entity2,
       PDM_g_num_t    **dentity2_entity1
);


void
PDM_dgroup_entity_transpose
(
 int            n_group,
 int           *dgroup_entity_idx,
 PDM_g_num_t   *dgroup_entity,
 PDM_g_num_t   *distrib_entity,
 int          **dentity_group_idx,
 int          **dentity_group,
 PDM_MPI_Comm   comm
);

/* Version of PDM_dentity_group_transpose valid for signed groups (is_signed = 1 for signed groups) */

void
PDM_dentity_group_signed_transpose
(
 int            n_group,
 int           *dentity_group_idx,
 int           *dentity_group,
 PDM_g_num_t   *distrib_entity,
 int          **dgroup_entity_idx,
 PDM_g_num_t  **dgroup_entity,
 PDM_MPI_Comm   comm,
 const int      is_signed
);

/**
 *
 * \brief Create group->entity connectivity
 *
 * \param [in]   n_group               Number of groups
 * \param [in]   dentity_group_idx     Group->Entity connectivity index
 * \param [in]   dentity_group         Group->Entity connectivity
 * \param [in]   distrib_entity        Entity distribution
 * \param [out]  dgroup_entity_idx     Entity->Group connectivity index
 * \param [out]  dgroup_entity         Entity->Group connectivity
 * \param [in]   comm                  PDM_MPI communicator
 *
 */

void
PDM_dentity_group_transpose
(
 int            n_group,
 int           *dentity_group_idx,
 int           *dentity_group,
 PDM_g_num_t   *distrib_entity,
 int          **dgroup_entity_idx,
 PDM_g_num_t  **dgroup_entity,
 PDM_MPI_Comm   comm
);

void
PDM_dconnectivity_to_extract_dconnectivity_bis
(
 const PDM_MPI_Comm    comm,
       int             n_selected_entity1,
       PDM_g_num_t    *select_entity1,
       PDM_g_num_t    *entity1_distribution,
       int            *dentity1_entity2_idx,
       PDM_g_num_t    *dentity1_entity2,
       PDM_g_num_t   **extract_entity1_distribution,
       PDM_g_num_t   **extract_entity2_distribution,
       int           **dextract_entity1_entity2_idx,
       PDM_g_num_t   **dextract_entity1_entity2,
       PDM_g_num_t   **dparent_entity1_g_num,
       PDM_g_num_t   **dparent_entity2_g_num,
       PDM_g_num_t   **entity1_old_to_new
);


void
PDM_dconnectivity_to_extract_dconnectivity_block
(
 const PDM_MPI_Comm          comm,
       int                   dn_extract_entity1,
       PDM_g_num_t          *dextract_gnum_entity1,
       PDM_g_num_t          *entity1_distribution,
       int                  *dentity1_entity2_idx,
       PDM_g_num_t          *dentity1_entity2,
       int                 **dextract_entity1_entity2_idx,
       PDM_g_num_t         **dextract_entity1_entity2,
       PDM_block_to_part_t **btp_entity1_to_extract_entity1,
       PDM_g_num_t         **extract_entity2_distribution,
       PDM_g_num_t         **dparent_entity2_g_num
);

void
PDM_dconnectivity_to_extract_dconnectivity
(
 const PDM_MPI_Comm          comm,
       int                   n_selected_entity1,
       PDM_g_num_t          *select_entity1,
       PDM_g_num_t          *entity1_distribution,
       int                  *dentity1_entity2_idx,
       PDM_g_num_t          *dentity1_entity2,
       PDM_g_num_t         **extract_entity1_distribution,
       PDM_g_num_t         **dparent_entity1_g_num,
       int                 **dextract_entity1_entity2_idx,
       PDM_g_num_t         **dextract_entity1_entity2,
       PDM_block_to_part_t **btp_entity1_to_extract_entity1,
       PDM_g_num_t         **extract_entity2_distribution,
       PDM_g_num_t         **dparent_entity2_g_num
);

void
PDM_dconnectivity_dface_vtx_from_face_and_edge
(
 const PDM_MPI_Comm    comm,
       PDM_g_num_t    *distrib_face,
       PDM_g_num_t    *distrib_edge,
       int            *dface_edge_idx,
       PDM_g_num_t    *dface_edge,
       PDM_g_num_t    *dedge_vtx,
       PDM_g_num_t   **dface_vtx
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DCONNECTIVITY_TRANSFORM_H_ */
