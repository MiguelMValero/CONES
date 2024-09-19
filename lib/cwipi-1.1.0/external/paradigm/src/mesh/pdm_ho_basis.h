/*
 * \file
 */

#ifndef __PDM_HO_BASIS_H__
#define __PDM_HO_BASIS_H__

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

#include "pdm_mesh_nodal.h"

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

typedef void (*PDM_ho_basis_fct_t)
(const int     entities_dim,
 const int     order,
 const int     n_nodes,
 const int     n_pts,
 const double *uvw,
 double       *weights);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Evaluate high-order basis functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_nodes   Number of nodes
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */

void
PDM_ho_basis
(
 const PDM_Mesh_nodal_elt_t type,
 const int                  order,
 const int                  n_nodes,
 const int                  n_pts,
 const double              *uvw,
 double                    *weights
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_BASIS_H__ */
