#ifndef __FVMC_HO_BASIS_H__
#define __FVMC_HO_BASIS_H__

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
 * Callback to define the basis functions of an high order 
 * element
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   n_pts             <-- number of points 
 *   uvw               <-- Parametric coordinates of points
 *   projected_uvw     --> Interpolation weights associated to uvw coordinates
 * 
 *----------------------------------------------------------------------------*/

typedef void (*fvmc_ho_basis_fct_t)
(const int entities_dim,
 const int order,
 const int n_nodes,
 const int n_pts,
 const double *uvw,
 double *weights);

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
FVMC_ho_basis_user_elt_set (fvmc_element_t elt_type,
                            fvmc_ho_basis_fct_t location_in_elt);


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
FVMC_ho_basis_user_elt_unset (fvmc_element_t elt_type);

void
FVMC_ho_basis_user_elts_unset (void);


/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_nodes         <-- number of nodes
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw (size = elt_dim * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

void
FVMC_ho_basis
(
const fvmc_element_t type,
const int order,
const int n_nodes,
const int n_pts,
const double *uvw,
      double *weights 
);

/*----------------------------------------------------------------------------
 * 
 * Free static variables
 * 
 *----------------------------------------------------------------------------*/


void
FVMC_ho_basis_free
(
 void
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_HO_BASIS_H__ */
