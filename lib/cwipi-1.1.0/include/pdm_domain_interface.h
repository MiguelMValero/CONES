/*
 * \file
 */

#ifndef __PDM_DOMAIN_INTERFACE_H__
#define __PDM_DOMAIN_INTERFACE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_to_block.h"

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

typedef struct _pdm_domain_interface_t PDM_domain_interface_t;
// typedef struct _pdm_part_domain_interface_t PDM_part_domain_interface_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_domain_interface_t *
PDM_domain_interface_create
(
 const int                   n_interface,
 const int                   n_domain,
 PDM_domain_interface_mult_t multidomain_interface,
 PDM_ownership_t             ownership,
 PDM_MPI_Comm                comm
);

void
PDM_domain_interface_set
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                    *interface_dn,
 PDM_g_num_t           **interface_ids,
 int                   **interface_dom
);

void
PDM_domain_interface_get
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_dn,
 PDM_g_num_t          ***interface_ids,
 int                  ***interface_dom
);

int PDM_domain_interface_get_as_graph
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_graph_idx,
 PDM_g_num_t           **interface_graph_ids,
 int                   **interface_graph_dom
);

void
PDM_domain_interface_translate_face2vtx
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
);


void
PDM_domain_interface_translate_vtx2edge
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_edge,
 int                    **dedge_vtx_idx,
 PDM_g_num_t            **dedge_vtx
);

void
PDM_domain_interface_translate_vtx2face
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
);

void
PDM_domain_interface_free
(
 PDM_domain_interface_t *dom_intrf
);


void
PDM_domain_interface_translate_entity1_entity2
(

 int                      n_domain,
 int                      n_interface,
 int                     *dn_entity1,
 int                     *dn_entity2,
 int                     *dn_interface,
 int                    **interface_dom,
 PDM_g_num_t            **interface_ids,
 int                    **dentity2_entity1_idx,
 PDM_g_num_t            **dentity2_entity1,
 int                      connectivity_is_signed,
 PDM_MPI_Comm             comm,
 int                    **interface_dn_entity2,
 PDM_g_num_t           ***interface_ids_entity2,
 int                   ***interface_dom_entity2
);

void
PDM_ddomain_interface_to_pdomain_interface
(
 PDM_MPI_Comm                   comm,
 int                            n_interface,
 int                            n_domain,
 PDM_domain_interface_mult_t    multidomain_intrf,
 PDM_bound_type_t               interface_kind,
 int                           *dn_interface,
 PDM_g_num_t                  **interface_ids,
 int                          **interface_dom,
 int                           *n_part,
 int                          **pn_entity,
 PDM_g_num_t                 ***entity_ln_to_gn,
 PDM_part_domain_interface_t   *pditrf
);


PDM_part_domain_interface_t*
PDM_domain_interface_to_part_domain_interface
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *n_part,
 int                    **pn_face,
 int                    **pn_edge,
 int                    **pn_vtx,
 PDM_g_num_t           ***face_ln_to_gn,
 PDM_g_num_t           ***edge_ln_to_gn,
 PDM_g_num_t           ***vtx_ln_to_gn
);


void
PDM_domain_interface_translation_set
(
        PDM_domain_interface_t  *dom_intrf,
        int                      i_interface,
  const double                  *vect
);

void
PDM_domain_interface_rotation_set
(
        PDM_domain_interface_t  *dom_intrf,
  const int                      i_interface,
  const double                  *direction,
  const double                  *center,
  const double                   angle
);

void
PDM_domain_interface_translation_get
(
        PDM_domain_interface_t       *dom_intrf,
        int                           i_interface,
        double                      **vect
);

void
PDM_domain_interface_rotation_get
(
        PDM_domain_interface_t       *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
);

void
PDM_part_domain_interface_to_domain_interface
(
  PDM_part_domain_interface_t    *dom_intrf,
  PDM_bound_type_t                interface_kind1,
  int                            *n_part,
  int                           **pn_entity1,
  PDM_g_num_t                  ***entity1_ln_to_gn,
  PDM_domain_interface_t        **ditrf_out,
  int                          ***is_entity1_on_itrf_out
);


void
PDM_domain_interface_make_flat_view
(
  PDM_domain_interface_t  *dom_intrf,
  PDM_bound_type_t         interface_kind,
  PDM_g_num_t             *shift_by_domain,
  PDM_part_to_block_t   ***ptb_interface_out,
  PDM_g_num_t           ***entity_opp_gnum_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOMAIN_INTERFACE_H__ */
