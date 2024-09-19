/*
 * \file
 */

#ifndef __PDM_SURF_PART_H__
#define __PDM_SURF_PART_H__

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


/**
 * \struct PDM_surf_part_t
 * \brief  Surface partition
 *
 *  PDM_surf_part_t defines a surface partition
 *
 */

typedef struct _pdm_surf_part_t PDM_surf_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Return an intialized \ref PDM_surf_part_t structure
 *
 * This function returns an initialized \ref PDM_surf_part_t structure
 *
 * \param [in]  n_face       Number of faces
 * \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]  face_vtx_idx  face -> vertex connectivity
 * \param [in]  face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]  n_vtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 *
 * \return      A new initialized \ref _part_t structure
 *
 */

PDM_surf_part_t *
PDM_surf_part_create
(
const int         n_face,
const int        *face_vtx_idx,
const int        *face_vtx,
const PDM_g_num_t *face_ln_to_gn,
const int         n_vtx,
const double     *coords,
const PDM_g_num_t *vtx_ln_to_gn
 );


/**
 * \brief Delete a \ref PDM_surf_part_t structure
 *
 * This function deletes a  PDM_surf_part_t structure
 *
 * \param [in]  part      part to delete
 *
 * \return     Null pointer
 */

PDM_surf_part_t *
PDM_surf_part_free
(
 PDM_surf_part_t * part
);


/**
 * \brief Compute partition edge entities
 *
 * This function defines edges of an initial partition and
 * computes edge connectivities
 *
 * \param [in]  part      Partition to compute
 *
 */

void
PDM_surf_part_build_edges
(
PDM_surf_part_t *part
);


/**
 * \brief Return face_ln_to_gn
 *
 *
 * \param [in]  part      Partition to compute
 *
 */

const PDM_g_num_t *
PDM_surf_part_faceLnToGn_get
(
PDM_surf_part_t *part
);

/**
 * \brief Dump a surf_part_t object
 *
 * This function dumps a PDM_surf_part_t structure
 *
 * \param [in]  part       to dump
 *
 */

void
PDM_surf_part_dump
(
PDM_surf_part_t *part
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_PART_H__ */
