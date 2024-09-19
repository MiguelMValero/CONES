#ifndef __FVMC_NODAL_EXTRACT_H__
#define __FVMC_NODAL_EXTRACT_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2006-2008  EDF

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

#include "fvmc_defs.h"
#include "fvmc_nodal.h"

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

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Copy global vertex numbers to an array.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   g_vtx_num  --> global vertex numbers (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_global_vertex_num(const fvmc_nodal_t  *this_nodal,
                                fvmc_gnum_t         *g_vtx_num);

/*----------------------------------------------------------------------------
 * Copy global element numbers of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the global element numbers are continued from one section to the next,
 * so to the user, all is as if the sections were concatenated.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   element_type  <-- type of elements to deal with
 *   g_elt_num     <-> pointer to global_element_num array (pre-allocated)
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_global_element_num(const fvmc_nodal_t  *this_nodal,
                                 fvmc_element_t       element_type,
                                 fvmc_gnum_t         *g_elt_num);

/*----------------------------------------------------------------------------
 * Copy vertex coordinates to an array.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   vertex_coords  --> vertices coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_vertex_coords(const fvmc_nodal_t  *this_nodal,
                            fvmc_interlace_t     interlace,
                            fvmc_coord_t        *vertex_coords);

/*----------------------------------------------------------------------------
 * Copy element centers to an array.
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   entity_dim     <-- dimension of entities we want to count (0 to 3)
 *   cell_centers   --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_element_centers(const fvmc_nodal_t  *this_nodal,
                              fvmc_interlace_t     interlace,
                              int                 entity_dim,
                              fvmc_coord_t        *cell_centers);

/*----------------------------------------------------------------------------
 * Copy element -> vertex connectivity of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the connectivity spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * Return local connectivity for sections of a given element_type.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   element_type  <-- type of elements of the section to deal with
 *   connectivity  <-> pointer to connectvity (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_strided_connect(const fvmc_nodal_t  *this_nodal,
                              fvmc_element_t       element_type,
                              fvmc_lnum_t         *connectivity);

/*----------------------------------------------------------------------------
 * Build inverse vertex -> element connectivity.
 *
 * The user is responsible for freeing the returned arrays.
 * The size of element_index[] should be n_vertices + 1, where n_vertices
 * is the value returned by fvmc_nodal_get_n_entities(this_nodal, entity_dim).
 * The size of element_id[] should be element_index[n_vertices].
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * Note also that both the vertex -> element index and connectivity arrays
 * returned use 0 to n numbering.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   entity_dim    <-- dimension of entities we want to count (1 to 3)
 *   element_index --> vertex -> element connectivity index (O to n-1)
 *   element_id    --> vertex -> element connectivity (0 to n-1)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_vertex_elements(const fvmc_nodal_t   *this_nodal,
                              int                  entity_dim,
                              fvmc_lnum_t         **element_index,
                              fvmc_lnum_t         **element_id);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_NODAL_EXTRACT_H__ */
