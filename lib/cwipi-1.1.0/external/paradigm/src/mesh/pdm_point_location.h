/*
 * \file
 */

#ifndef __PDM_POINT_LOCATION_H__
#define __PDM_POINT_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"

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
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Locate a set points inside a set of elements
 *
 * Elements are ordered by type (points, lines, triangles, quadrangles,
 * polygons, tetrahedra, pyramids, prisms, hexahedra, polyhedra).
 *
 * \param [in]   type_idx           Index for the element types (size = 11)
 * \param [in]   elt_vtx_idx        Index of the element-vertex connectivity
 * \param [in]   elt_vtx_coord      Coordinates of the elements' vertices
 * \param [in]   poly3d_face_idx    Index of the element-face connectivity (only for polyhedra)
 * \param [in]   face_vtx_idx       Index for the face-vertex connectivity
 * \param [in]   face_vtx           Face-vertex connectivity
 * \param [in]   face_orientation   Orientation of the faces
 * \param [in]   pts_idx            Index of points (size = n_elt + 1)
 * \param [in]   pts_coord          Coordinates of the points to locate
 * \param [in]   tolerance          Geometric tolerance
 * \param [out]  distance           Distance from points to elements (< 0 if inside, > 0 if outside)
 * \param [out]  projected_coord    Coordinates of the projection of the points on the elements
 * \param [out]  bar_coord_idx      Index for the mean-value coordinates of the projections
 * \param [out]  bar_coord          Mean-value coordinates of the projections
 *
 */

void
PDM_point_location_nodal
(
 PDM_part_mesh_nodal_elmts_t   *pmne,
 const int                      n_part,
 const double                 **pvtx_coord,
 const int                    **pts_idx,
 const double                 **pts_coord,
 const double                   tolerance,
 double                      ***distance,
 double                      ***projected_coord,
 int                         ***bar_coord_idx,
 double                      ***bar_coord,
 double                      ***uvw
 );


/**
 * \brief Compute quadrangle, hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * \param [in]   elt_type       Type of element
 * \param [in]   point_coords   Point coordinates
 * \param [in]   vertex_coords  Pointer to element vertex coordinates
 * \param [in]   tolerance      Location tolerance factor
 * \param [out]  uvw            Parametric coordinates of point in element
 * \param [in]   init_uvw       Initial uvw guess for Newton method (or NULL)
 *
 *  \return Convergence status of Newton method
 */

PDM_bool_t
PDM_point_location_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
       double               uvw[3],
       double               init_uvw[3]
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
