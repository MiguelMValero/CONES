/*
 * \file
 * \author equemera
 *
 * \date August 31, 2017, 1:24 PM
 */

#ifndef PDM_FACECELL_ORIENT_H
#define	PDM_FACECELL_ORIENT_H


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

#ifdef	__cplusplus
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


/**
 * \brief Orient cell->face connectivity
 *
 * At the output of th function, a face number in \ref cell_face is positive
 * if surface normal is inside the cell, negative otherwise. \ref face_cell is
 * oriented in the same way
 *
 * \param [in]      n_cell         Number of cells
 * \param [in]      n_face         Number of faces
 * \param [in]      n_vtx          Number of vertices
 * \param [in]      coords         Vertices coordinates
 * \param [in]      cell_face_idx  Cell to face connectivity index (size = \ref n_cell + 1)
 * \param [in, out] cell_face      Cell to face connectivity (size = cell_face_idx[n_cell])
 * \param [in, out] face_cell      Face to cell connectivity (size = 2 * \ref n_face) or NULL
 * \param [in]      face_vtx_idx   Face to vertex connectivity index (size = \ref n_face + 1)
 * \param [in]      face_vtx       Face to vertex connectivity (size = face_vtx_idx[n_face])
 *
 */

void
PDM_cellface_orient
(
const int      n_cell,
const int      n_face,
const int      n_vtx,
const double  *coords,
const int     *cell_face_idx,
int           *cell_face,
int           *face_cell,
const int     *face_vtx_idx,
const int     *face_vtx
);


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_FACECELL_ORIENT_H */
