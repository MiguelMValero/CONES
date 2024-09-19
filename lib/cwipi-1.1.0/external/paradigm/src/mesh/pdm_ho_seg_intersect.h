/*
 * \file
 */

#ifndef __PDM_HO_SEG_INTERSECT_H__
#define __PDM_HO_SEG_INTERSECT_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"

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

#define _MIN(a,b) ((a) < (b) ? (a) : (b))
#define _MAX(a,b) ((a) > (b) ? (a) : (b))

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Determine intersection between P1 element and segment from line intersecting ho face bounding box
 *
 * \param [in]   t_elt                           Element type
 * \param [in]   n_back_face_to_intersect        Number of background mesh faces to intersect
 * \param [in]   back_face_to_intersect_ln_to_gn Set of background mesh faces to intersect
 * \param [in]   back_face_vtx                   Background face -> vertex connectivity
 * \param [in]   back_vtx_coord                  Coordinates of vertices of the background mesh
 * \param [in]   back_face_line_idx              Index of background face -> line connectivity
 * \param [in]   back_face_line                  Background face -> line connectivity
 * \param [in]   line_boxes_idx                  Index of line -> box connectivity
 * \param [in]   box_line_intersect_points       Line -> box intersection points connectivity
 * \param [out]  newton_initial_point            Initial point for Newthon Method
 *
 */

void
PDM_ho_seg_intersect_P1_line
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 int                         n_back_face_to_intersect,
 PDM_g_num_t                *back_face_to_intersect_ln_to_gn,
 PDM_g_num_t                *back_face_vtx,
 double                     *back_vtx_coord,
 int                        *back_face_line_idx,
 PDM_g_num_t                *back_face_line,
 int                        *line_boxes_idx,
 double                     *box_line_intersect_points,
 double                    **newton_initial_point
 );

/**
 * \brief Compute intersection between line and ho element using Newton Method
 *
 * \param [in]   type                            Element type
 * \param [in]   n_line                          Number of lines
 * \param [in]   line_coords                     Coordinates of lines
 * \param [in]   n_back_face_to_intersect        Number of background mesh faces to intersect
 * \param [in]   back_face_to_intersect_ln_to_gn Set of background mesh faces to intersect
 * \param [in]   back_face_vtx                   Background face -> vertex connectivity
 * \param [in]   back_vtx_coord                  Coordinates of vertices of the background mesh
 * \param [in]   back_face_line_idx              Index of background face -> line connectivity
 * \param [in]   back_face_line                  Background face -> line connectivity
 * \param [out]  line_box_intersection_point     Line->Box->Point connectivity
 *
 */

void
PDM_ho_seg_intersect_compute
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   n_line,
 double                     *line_coords,
 int                         n_back_face_to_intersect,
 PDM_g_num_t                *back_face_to_intersect_ln_to_gn,
 PDM_g_num_t                *back_face_vtx,
 double                     *back_vtx_coord,
 int                        *back_face_line_idx,
 PDM_g_num_t                *back_face_line,
 double                    **line_box_intersection_point
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_SEG_INTERSECT_H__ */
