/*
 * \file
 */

#ifndef __PDM_FIELD_CELL_TO_VTX_H__
#define __PDM_FIELD_CELL_TO_VTX_H__


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_mesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

typedef struct _pdm_field_cell_to_vtx_t PDM_field_cell_to_vtx_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// Function set pour le inverse distance weighting
// Regarder  dans cwipi pour des arguments variadic

// JC: interterp_kind dans le exch ? Vu que de toute manière les infos les mêmes pour la construction


PDM_field_cell_to_vtx_t*
PDM_field_cell_to_vtx_create
(
 const int                            n_domain,
 const int                           *n_part,
 const int                           *n_group,
 const PDM_cell_to_vtx_interp_kind_t  interp_kind, // IDW(p), RBF, LSQ, USER
 const int                            n_depth,
 const PDM_MPI_Comm                   comm
);


/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   fctv
 *
 */
void
PDM_field_cell_to_vtx_compute
(
  PDM_field_cell_to_vtx_t *fctv
);



void
PDM_field_cell_to_vtx_part_set
(
  PDM_field_cell_to_vtx_t   *mi,
  int                       i_domain,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);


void
PDM_field_cell_to_vtx_part_mesh_nodal_set
(
  PDM_field_cell_to_vtx_t   *mi,
  int                        i_domain,
  PDM_part_mesh_nodal_t     *pmn
);


void
PDM_field_cell_to_vtx_graph_comm_set
(
  PDM_field_cell_to_vtx_t   *mi,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                      *entity_part_bound_proc_idx,
  int                      *entity_part_bound_part_idx,
  int                      *entity_part_bound
);

void
PDM_field_cell_to_vtx_part_group_set
(
  PDM_field_cell_to_vtx_t   *mi,
  int                       i_domain,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity
);

void
PDM_field_cell_to_vtx_part_domain_interface_shared_set
(
  PDM_field_cell_to_vtx_t      *mi,
  PDM_part_domain_interface_t *pdi
);

void
PDM_field_cell_to_vtx_inverse_distance_weighting_p_set
(
  PDM_field_cell_to_vtx_t *fctv,
  int                      p
);

void
PDM_field_cell_to_vtx_exch
(
        PDM_field_cell_to_vtx_t      *mi,
        PDM_field_kind_t            field_kind,
        double                   ***local_field,
        double                  ****bound_field,
        double                  ****result_field
);


void
PDM_field_cell_to_vtx_set_cell_center
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  double                   *cell_center
);


void
PDM_field_cell_to_vtx_set_vtx_center
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  double                   *vtx_center
);



void
PDM_field_cell_to_vtx_free
(
 PDM_field_cell_to_vtx_t *mi
);



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_field_cell_to_vtx_H__ */

