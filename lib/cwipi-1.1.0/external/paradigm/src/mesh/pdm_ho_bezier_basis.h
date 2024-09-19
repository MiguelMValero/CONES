/*
 * \file
 */

#ifndef __PDM_HO_BEZIER_BASIS_H__
#define __PDM_HO_BEZIER_BASIS_H__

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

/**
 *
 * \brief Evaluate high-order basis Bézier functions
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] weights   Weights (size = \ref n_pts * \ref n_nodes)
 *
 */

void
PDM_ho_bezier_basis
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double                     *weights
);

/**
 *
 * \brief Evaluate high-order basis Bézier functions derivatives
 *
 *
 * \param [in]  type      Element type structure
 * \param [in]  order     Element order
 * \param [in]  n_pts     Number of points
 * \param [in]  uvw       Parametric coordinates of the points (size = elt_dim * \ref n_pts)
 * \param [out] dw_du   Weights derivatives with respect to u
 * \param [out] dw_dv   Weights derivatives with respect to v
 * \param [out] dw_dw   Weights derivatives with respect to w
 *
 */

void
PDM_ho_bezier_basis_derivative
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_pts,
 const double               *uvw,
 double            *restrict dw_du,
 double            *restrict dw_dv,
 double            *restrict dw_dw
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_BEZIER_BASIS_H__ */
