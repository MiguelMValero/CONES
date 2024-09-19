/*
 * \file
 */

#ifndef __PDM_DGEOM_ELEM_H__
#define __PDM_DGEOM_ELEM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
PDM_compute_center_from_descending_connectivity
(
  const int         *dentity1_entity2_idx,
  const PDM_g_num_t *dentity1_entity2,
  const int          dn_entity1,
  const PDM_g_num_t *dentity2_distrib,
  double            *dentity1_coord,
  double            *dentity2_coord,
  PDM_MPI_Comm       comm
);


void
PDM_compute_dface_normal
(
  const int         *dface_vtx_idx,
  const PDM_g_num_t *dface_vtx,
  const int          dn_face,
  const PDM_g_num_t *dvtx_distrib,
  double            *dvtx_coord,
  double            *dface_normal,
  PDM_MPI_Comm       comm
);


void
PDM_compute_vtx_characteristic_length
(
 PDM_MPI_Comm    comm,
 int             dn_face,
 int             dn_edge,
 int             dn_vtx,
 int            *dface_vtx_idx,
 PDM_g_num_t    *dface_vtx,
 PDM_g_num_t    *dedge_vtx,
 double         *dvtx_coord,
 double        **dchar_length_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DGEOM_ELEM_H__ */
