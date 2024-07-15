/*
 * \file
 */

#ifndef __PDM_DCUBE_GEN_H__
#define __PDM_DCUBE_GEN_H__

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
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dcube_t PDM_dcube_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a distributed cube
 *
 * \param [in]   comm           Communicator
 * \param [in]   n_vtx_seg      Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin
 * \param [in]   zero_y         Coordinates of the origin
 * \param [in]   zero_z         Coordinates of the origin
 * \param [in]   ownership      Choice of ownership of the results \ref PDM_ownership_t
 *
 * \return      Pointer to \ref PDM_dcube_t object
 *
 */

PDM_dcube_t*
PDM_dcube_gen_init
(
      PDM_MPI_Comm     comm,
const PDM_g_num_t      n_vtx_seg,
const double           length,
const double           zero_x,
const double           zero_y,
const double           zero_z,
      PDM_ownership_t  owner
);


/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   pdm_dcube     Pointer to \ref PDM_dcube_t object
 * \param [out]  n_face_group  Number of faces groups
 * \param [out]  dn_cell       Number of cells stored in this process
 * \param [out]  dn_face       Number of faces stored in this process
 * \param [out]  dn_vtx        Number of vertices stored in this process
 * \param [out]  sface_vtx     Length of dface_vtx array
 * \param [out]  sface_group   Length of dface_group array
 *
 */

void
PDM_dcube_gen_dim_get
(
 PDM_dcube_t        *pdm_dcube,
 int                *n_face_group,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
);

/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  pdm_dcube       Pointer to \ref PDM_dcube_t object
 * \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
 * \param [out] dface_vtx_idx   Faces from vertices connectivity index (size = dn_face + 1)
 * \param [out] dface_vtx       Faces from vertices connectivity (size = dface_vtxL)
 * \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
 * \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
 * \param [out] dface_group     Faces groups (size = dFacegroupL)
 *
 */

void
PDM_dcube_gen_data_get
(
 PDM_dcube_t       *pdm_dcube,
 PDM_g_num_t      **dface_cell,
 int              **dface_vtx_idx,
 PDM_g_num_t      **dface_vtx,
 double           **dvtx_coord,
 int              **dface_group_idx,
 PDM_g_num_t      **dface_group
);

/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  pdm_dcube     Pointer to \ref PDM_dcube_t object
 *
 */

void
PDM_dcube_gen_free
(
 PDM_dcube_t       *pdm_dcube
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PART_DCUBE_H__ */
