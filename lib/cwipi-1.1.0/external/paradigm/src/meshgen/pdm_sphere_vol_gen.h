/*
 * \file
 */

#ifndef __PDM_SPHERE_VOL_GEN_H__
#define __PDM_SPHERE_VOL_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/**
 *
 * \brief Create a volume, nodal mesh bounded by a sphere (deformed cube)
 *
 * \note The output mesh is *block-distributed*,
 * in the form of a \ref PDM_dmesh_nodal_t object
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n_vtx_x         Number of vertices on segments in x-direction
 * \param[in]  n_vtx_y         Number of vertices on segments in y-direction
 * \param[in]  n_vtx_z         Number of vertices on segments in z-direction
 * \param[in]  radius          Sphere radius
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  t_elt           Element type
 * \param[in]  order           Element order
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal_t object
 *
 */

void
PDM_sphere_vol_gen_nodal
(
 PDM_MPI_Comm           comm,
 const PDM_g_num_t      n_vtx_x,
 const PDM_g_num_t      n_vtx_y,
 const PDM_g_num_t      n_vtx_z,
 const double           radius,
 const double           x_center,
 const double           y_center,
 const double           z_center,
 PDM_Mesh_nodal_elt_t   t_elt,
 const int              order,
 PDM_dmesh_nodal_t    **dmn
 );


/**
 * \brief Create a volume mesh bounded by an icosphere
 * (all cells are tetrahedra, all faces are triangles)
 *
 * \note The output mesh is *block-distributed*.
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dvtx_coord      Vertex coordinates
 * \param[out] dface_vtx       Face -> vertex connectivity (global ids)
 * \param[out] dcell_vtx       Cell -> face connectivity (global ids)
 * \param[out] distrib_vtx     Vertex distribution (size : *n_rank* + 1)
 * \param[out] distrib_face    Face distribution (size : *n_rank* + 1)
 * \param[out] distrib_face    Cell distribution (size : *n_rank* + 1)
 *
 */

void
PDM_sphere_vol_icosphere_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **dcell_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face,
       PDM_g_num_t       **distrib_cell
);


/**
 * \brief Create a volume, nodal mesh bounded by an icosphere
 * (all cells are tetrahedra)
 *
 * \note The output mesh is *block-distributed*,
 * in the form of a \ref PDM_dmesh_nodal_t object
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_icosphere_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **dmn
);


/**
 * \brief Create a volume, nodal mesh bounded by two concentric icospheres
 * (all cells are prisms)
 *
 * \note The output mesh is *block-distributed*,
 * in the form of a \ref PDM_dmesh_nodal_t object
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  n_layer         Number of extrusion layers
 * \param[in]  x_center        x-coordinate of the spheres center
 * \param[in]  y_center        y-coordinate of the spheres center
 * \param[in]  z_center        z-coordinate of the spheres center
 * \param[in]  radius_interior Interior sphere radius
 * \param[in]  radius_exterior Exterior sphere radius
 * \param[in]  geometric_ratio Geometric ratio for layer thickness
 * \param[out] dmn             Pointer to a \ref PDM_dmesh_nodal object
 *
 */

void
PDM_sphere_vol_hollow_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const PDM_g_num_t         n_layer,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius_interior,
 const double              radius_exterior,
 const double              geometric_ratio,
       PDM_dmesh_nodal_t **dmn
);

#ifdef __cplusplus
}
#endif
#endif
