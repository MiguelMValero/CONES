/*
 * \file
 */

#ifndef __PDM_DMESH_NODAL_TO_DMESH_H__
#define __PDM_DMESH_NODAL_TO_DMESH_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"

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

typedef struct _pdm_dmesh_nodal_to_dmesh_t PDM_dmesh_nodal_to_dmesh_t;


typedef enum {
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE = 0,
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE = 1,
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_VTX  = 2,
} PDM_dmesh_nodal_to_dmesh_transform_t;


typedef enum {
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE    = 0,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE = 1,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE = 2,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX  = 3,
} PDM_dmesh_nodal_to_dmesh_translate_group_t;


/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Generate the deduce connectivity in distribute manner
 *
 * \param [in]   comm     MPI communicator
 *
 */
void
PDM_generate_entitiy_connectivity_raw
(
PDM_MPI_Comm   comm,
PDM_g_num_t    n_vtx_abs,
int            n_entity_elt_tot,
PDM_g_num_t   *delmt_entity,
int           *delmt_entity_vtx_idx,
PDM_g_num_t   *delmt_entity_vtx,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt
);

/**
 * \brief Generate the deduce connectivity in distribute manner
 *
 * \param [in]   comm     MPI communicator
 *
 */
void
PDM_generate_entitiy_connectivity
(
PDM_MPI_Comm   comm,
PDM_g_num_t    n_vtx_abs,
PDM_g_num_t    n_g_child,
int            n_part,
int           *n_entity_elt_tot,
PDM_g_num_t  **delmt_entity,
int          **delmt_entity_vtx_idx,
PDM_g_num_t  **delmt_entity_vtx,
int          **dparent_elmt_position,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt,
int          **dentity_parent_element_position,
int          **dparent_idx,
PDM_g_num_t  **dparent_gnum,
int          **dparent_sign,
PDM_g_num_t  **delmt_child_distrib,
PDM_g_num_t  **distrib_missing_child,
PDM_g_num_t  **dmissing_child_parent_g_num
);

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
);

void
PDM_dmesh_nodal_to_dmesh_set_post_treat_result
(
 PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm,
 int                         post_treat_result
);

void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         i_mesh,
        PDM_dmesh_nodal_t          *dmn
);

void
PDM_dmesh_nodal_to_dmesh_compute
(
        PDM_dmesh_nodal_to_dmesh_t                 *dmn_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
);

void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmn_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
);

void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
);


void
PDM_dmesh_nodal_to_dmesh_get_missing
(
 PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
 const int                   i_mesh,
 PDM_geometry_kind_t         geom_kind,
 PDM_g_num_t               **distrib_missing,
 PDM_g_num_t               **dmissing_parent_g_num
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
