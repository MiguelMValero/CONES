#ifndef __FVMC_HO_LOCATION_H__
#define __FVMC_HO_LOCATION_H__

/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

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

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"

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
 * 
 * Callback to define location in a high order element
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point to locate coordinates
 *   projected_coords  --> projected point coordinates (if point is outside) 
 *   projected_uvw     --> parametric coordinates of the projected point
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

typedef double (*fvmc_ho_location_fct_t)
(const int entities_dim,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *uvw);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 
 * Set elementary functions
 * 
 * parameters:
 *   location_tetra    <-- Location in a tetrahedron
 *   location_prism    <-- Location in a prism
 *   location_pyramid  <-- Location in a pyramid
 *   location_hexa     <-- Location in a hexaedron
 *   location_tria     <-- Location on a triangle
 *   location_quad     <-- Location on a quandragle
 *   location_edge     <-- Location on a edge
 *   interp_tetra       <-- Interpolation in a tetrahedron
 *   interp_prism       <-- Interpolation in a prism
 *   interp_pyramid     <-- Interpolation in a pyramid
 *   interp_hexa        <-- Interpolation in a hexaedron
 *   interp_tria        <-- Interpolation on a triangle
 *   interp_quad        <-- Interpolation on a quandragle
 *   interp_edge        <-- Interpolation on a edge
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_location_user_elt_set (fvmc_element_t elt_type,
                               fvmc_ho_location_fct_t location_in_elt);


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_location_user_elt_unset (fvmc_element_t elt_type);

void
fvmc_ho_location_user_elts_unset (void);

/*----------------------------------------------------------------------------
 * 
 * Point location in a high order cell
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates (size = 3 * n_nodes)
 *   point_coords     <-- point to locate coordinates (size = 3)
 *   projected_coords --> projected point coordinates (size = 3)
 *   uvw              --> parametric coordinates of the projected point on the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location
(
 const fvmc_element_t type,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *uvw
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_HO_H__ */
