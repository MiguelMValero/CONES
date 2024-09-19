/*
 * \file
 */

#ifndef __PDM_PART_DOMAIN_INTERFACE_H__
#define __PDM_PART_DOMAIN_INTERFACE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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


typedef struct _pdm_part_domain_interface_t PDM_part_domain_interface_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_part_domain_interface_t*
PDM_part_domain_interface_create
(
const int                    n_interface,
const int                    n_domain,
const int                   *n_part,
PDM_domain_interface_mult_t  multidomain_interface,
PDM_ownership_t              ownership,
PDM_MPI_Comm                 comm
);


void
PDM_part_domain_interface_set
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind,
 int                           i_domain,
 int                           i_part,
 int                           i_interface,
 int                           interface_pn,
 PDM_g_num_t                  *interface_ln_to_gn,
 int                          *interface_sgn,
 int                          *interface_sens,
 int                          *interface_ids,
 int                          *interface_ids_idx,
 int                          *interface_dom
);


void
PDM_part_domain_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind,
 int                            i_domain,
 int                            i_part,
 int                            i_interface,
 int                           *interface_pn,
 PDM_g_num_t                  **interface_ln_to_gn,
 int                          **interface_sgn,
 int                          **interface_sens,
 int                          **interface_ids,
 int                          **interface_ids_idx,
 int                          **interface_dom
);

int
PDM_part_domain_interface_n_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf
);

int
PDM_part_domain_interface_exist_get
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind
);

void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
);

void
PDM_part_domain_interface_translation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
  const double                       *vect
);

void
PDM_part_domain_interface_rotation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
  const double                       *direction,
  const double                       *center,
  const double                        angle
);

void
PDM_part_domain_interface_translation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
        double                      **vect
);

void
PDM_part_domain_interface_rotation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
);


void
PDM_part_domain_interface_as_graph
(
  PDM_part_domain_interface_t    *dom_intrf,
  PDM_bound_type_t                interface_kind,
  int                           **n_entity,
  PDM_g_num_t                  ***entity_ln_to_gn,
  int                          ***neighbor_entity_idx,
  int                          ***neighbor_entity_desc,
  int                            *n_g_interface,
  int                           **composed_id_idx,
  int                           **composed_id,
  PDM_g_num_t                   **composed_ln_to_gn_sorted
);

void
PDM_part_domain_interface_translate
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind1,
 PDM_bound_type_t               interface_kind2,
 int                           *n_part,
 int                          **pn_entity1,
 int                          **pn_entity2,
 PDM_g_num_t                 ***entity2_ln_to_gn,
 int                         ***entity1_entity2_idx,
 int                         ***entity1_entity2
);


/**
 *
 * \brief Create from all interface information a view by partition
 *
 * \param [in]    pdi                       part_domain_interface structure
 * \param [in]    interface_kind            Kind of interface
 * \param [out]   pn_entity_num             Number of entity by interface (size = n_part along all domain )
 * \param [out]   pentity_num               Index of entity by partition
 * \param [out]   pentity_opp_interface_idx For each index, idx of connected triplet
 * \param [out]   pentity_opp_location      For each index, all opposite triplet
 * \param [out]   pentity_opp_interface     For each index, by which interface the entity is connected
 *
 */
void
PDM_part_domain_interface_view_by_part
(
  PDM_part_domain_interface_t   *pdi,
  PDM_bound_type_t               interface_kind,
  int                           *pn_entity,
  PDM_g_num_t                  **pentity_ln_to_gn,
  int                          **pn_entity_num,
  int                         ***pentity_num,
  int                         ***pentity_opp_location_idx,
  int                         ***pentity_opp_location,
  int                         ***pentity_opp_interface_idx,
  int                         ***pentity_opp_interface,
  int                         ***pentity_opp_sens,
  PDM_g_num_t                 ***pentity_opp_gnum_out
);

void
PDM_part_domain_interface_add
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind1,
 PDM_bound_type_t               interface_kind2,
 int                           *n_part,
 int                          **pn_entity1,
 PDM_g_num_t                 ***entity1_ln_to_gn,
 int                          **pn_entity2,
 PDM_g_num_t                 ***entity2_ln_to_gn,
 int                         ***entity2_entity1_idx,
 int                         ***entity2_entity1,
 int                            connectivity_is_signed
);

void
PDM_part_domain_interface_face2vtx
(
 PDM_part_domain_interface_t   *dom_intrf,
 int                           *n_part,
 int                          **pn_face,
 PDM_g_num_t                 ***pface_ln_to_gn,
 int                          **pn_vtx,
 PDM_g_num_t                 ***pvtx_ln_to_gn,
 int                         ***pface_vtx_idx,
 int                         ***pface_vtx
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_DOMAIN_INTERFACE_H__ */
