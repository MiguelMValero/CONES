/*
 * \file
 */

#ifndef __PDM_TRIANGULATE_H__
#define __PDM_TRIANGULATE_H__

/*============================================================================
 * Triangulation of a polygon
 *============================================================================*/
/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2005-2007  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*
 * Pointer to structure maintaining the state of the current triangulation;
 * the structure itself is private.
 */

typedef struct _PDM_triangulate_state_t PDM_triangulate_state_t;

/*
 * Describe how the resulting triangle connectivity is defined.
 */

typedef enum {

  PDM_TRIANGULATE_MESH_DEF,      /* Definition by mesh vertex numbers */
  PDM_TRIANGULATE_ELT_DEF        /* Definition by local (element) vertex
                                    position (1 to n) */

} PDM_triangulate_def_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Create a structure necessary to the polygon triangulation algorithm
 *
 * \param [in]   n_vertices_max  Maximum expected number of vertices per polygon
 *
 * \return  Pointer to polygon triangulation state structure
 *
 */

PDM_triangulate_state_t *
PDM_triangulate_state_create(const int  n_vertices_max);


/**
 * \brief Destroy a structure necessary to the polygon triangulation algorithm
 *
 * \param [in]   this_state  Pointer to structure that should be destroyed
 *
 * \return  NULL pointer
 *
 */

PDM_triangulate_state_t *
PDM_triangulate_state_destroy(PDM_triangulate_state_t  *this_state);


/**
 * \brief Triangulate a polygonal face
 *
 * For a polygon with n vertices, we should obtain a triangluation with
 * (n-2) triangles and (2n-3) edges. If the polygon_vertices argument
 * is NULL, 1, 2, ...,n local numbering is implied.
 *
 * \param [in]   dim                Spatial dimension (2 or 3)
 * \param [in]   n_vertices         Number of vertices defining the polygon
 * \param [in]   coords             Coordinates of the triangulation's vertices
 * \param [in]   parent_vertex_num  Optional indirection to vertex coordinates (1 to n)
 * \param [in]   polygon_vertices   Polygon connectivity; size: n_vertices or empty
 * \param [in]   mode               Triangles connectivity by vertex number or
 *                                  polygon vertex index (1 to n)
 * \param [out]   triangle_vertices Triangles connectivity
 *                                  (size = (\ref n_vertices - 2) * 3)
 * \param [in]   state              Associated triangulation state structure
 *
 * \return  Number of resulting triangles
 *
 */

int
PDM_triangulate_polygon(int                             dim,
                        int                             n_vertices,
                        const double                    coords[],
                        const PDM_l_num_t               parent_vertex_num[],
                        const PDM_l_num_t               polygon_vertices[],
                        PDM_triangulate_def_t           mode,
                        PDM_l_num_t                     triangle_vertices[],
                        PDM_triangulate_state_t  *const state);


/**
 * \brief Triangulate a quadrangle
 *
 * A convex quadrangle is divided into two triangles along its shortest
 * diagonal. A non-convex quadrangle may only be divided along the diagonal
 * which lies inside the quadrangle.
 *
 * \param [in]   dim                  Spatial dimension (2 or 3)
 * \param [in]   coords               Coordinates of the triangulation's vertices
 * \param [in]   parent_vertex_num    Optional indirection to vertex coordinates (1 to 4)
 * \param [in]   quadrangle_vertices  Quadrangle connectivity (size = 4 or empty)
 * \param [out]  triangle_vertices    Triangles connectivity (size = 2 * 3)
 * \param [in]   state                Associated triangulation state structure
 *
 * \return  Number of resulting triangles
 *
 */

int
PDM_triangulate_quadrangle(int                dim,
                           const double       coords[],
                           const PDM_l_num_t  parent_vertex_num[],
                           const PDM_l_num_t  quadrangle_vertices[],
                           PDM_l_num_t        triangle_vertices[]);


/**
 * \brief Tetrahedralize a pyramid
 *
 * \param [in]   dim               Spatial dimension (3)
 * \param [in]   coords            Coordinates of the vertices
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 5)
 * \param [in]   pyramid_vertices  Pyramid connectivity (size = 5)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 2 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_pyramid (int               dim,
                         const double      coords[],
                         const PDM_l_num_t parent_vertex_num[],
                         const PDM_l_num_t pyramid_vertices[],
                         PDM_l_num_t       tetra_vertices[]);


/**
 * \brief Tetrahedralize a prism
 *
 * A simple look-up table is currently used,
 * thus no validity check is performed on the tetrahedra.
 *
 * \param [in]   dim               Spatial dimension (3) (unused)
 * \param [in]   coords            Coordinates of the vertices (unused)
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 6) (unused)
 * \param [in]   prism_vertices    Prism connectivity (size = 6)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 3 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_prism (int               dim,
                       const double      coords[],
                       const PDM_l_num_t parent_vertex_num[],
                       const PDM_l_num_t prism_vertices[],
                       PDM_l_num_t       tetrahedron_vertices[]);


/**
 * \brief Tetrahedralize a hexahedron
 *
 * A simple look-up table is currently used,
 * thus no validity check is performed on the tetrahedra.
 *
 * \param [in]   dim               Spatial dimension (3) (unused)
 * \param [in]   coords            Coordinates of the vertices (unused)
 * \param [in]   parent_vertex_num Optional indirection to vertex coordinates (1 to 8) (unused)
 * \param [in]   hexa_vertices     Hexahedron connectivity (size = 8)
 * \param [out]  tetra_vertices    Tetrahedra connectivity (size = 5 * 4)
 *
 * \return  Number of resulting tetrahedra
 *
 */

int
PDM_triangulate_hexahedron (int               dim,
                            const double      coords[],
                            const PDM_l_num_t parent_vertex_num[],
                            const PDM_l_num_t hexa_vertices[],
                            PDM_l_num_t       tetrahedron_vertices[]);
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_TRIANGULATE_H__ */
