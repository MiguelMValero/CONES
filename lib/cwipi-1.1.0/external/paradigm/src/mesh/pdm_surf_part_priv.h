#ifndef __PDM_SURF_PART_PRIV_H__
#define __PDM_SURF_PART_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_bound.h"

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
 * \struct _surf_part_t
 * \brief  Mesh partition type
 *
 * _surf_part_t defines a mesh partition structure
 *
 */

struct _pdm_surf_part_t {

  int               n_face;       /*!< Number of faces  */
  double           *faceNormal;  /*!< Face normal  */

  int               nGhostFace;  /*!< Number of ghost faces  */
  int               nTotalFace;  /*!< Total number of faces  */
  int               sface_vtx;    /*!< Size of \ref face_vtx */
  const int        *face_vtx_idx;  /*!< Index in \ref face_vtx */
  const int        *face_vtx;     /*!< face -> vertex connectivity */
  int              *faceEdgeIdx; /*!< face -> edge connectivity (NULL if not computed) */
  int              *faceEdge;    /*!< face -> edge connectivity (NULL if not computed) */
  const PDM_g_num_t *face_ln_to_gn;  /*!< Local face numbering to global face numbering */
  int               n_vtx;        /*!< Number of vertices */
  const double     *coords;      /*!< Vertex coordinates */
  int              *vtxEdgeIdx;  /*!< Index in \ref vtxEdge (NULL if not computed) */
  int              *vtxEdge;     /*!< Vertex -> edge connectivity
                                      (NULL if not computed) */
  const PDM_g_num_t *vtx_ln_to_gn;   /*!< Local vertex numbering to global vertex numbering */
  int               nEdge;       /*!< Number of edges */
  int               nGhostEdge;  /*!< Number of ghost edges */
  int               nTotalEdge;  /*!< Total number of edges */
  int              *edgeFace;    /*!< Edge -> face connectivity
                                      (NULL if not computed) */
  int              *edgeVtx;     /*!< Edge -> vertex connectivity
                                      (NULL if not computed) */
  PDM_g_num_t       *edgeLnToGn;  /*!< Local edge numbering to global edge numbering */

  PDM_part_bound_t *edgePartBound; /*!< Inter partition edges */
  PDM_part_bound_t *vtxPartBound;  /*!< Inter partition vertices */

  double           *carLgthVtx;   /*!< Caracteristic length of vertices */
  double           *extents;      /*!< xmin, ymin, zmin, xmax, ymax, zmax for each face
                                       (size = dim * n_face * 2) */

} ;


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_PART_PRIV_H__ */
