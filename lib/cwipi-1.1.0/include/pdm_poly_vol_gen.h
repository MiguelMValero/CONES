/*
 * \file
 */

#ifndef __PDM_POLY_VOL_GEN_H__
#define __PDM_POLY_VOL_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2020-2023  ONERA

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
 * \brief Generate a distributed polyhedral mesh
 *
 * \param [in]   pdm_comm         MPI communicator
 * \param [in]   xmin             Minimum x-coordinate
 * \param [in]   ymin             Minimum y-coordinate
 * \param [in]   zmin             Minimum z-coordinate
 * \param [in]   lengthx          Length in the x-direction
 * \param [in]   lengthy          Length in the y-direction
 * \param [in]   lengthz          Length in the z-direction
 * \param [in]   nx               Number of vertices in the x-direction
 * \param [in]   ny               Number of vertices in the y-direction
 * \param [in]   nz               Number of vertices in the z-direction
 * \param [in]   randomize        Enable/disable randomization
 * \param [in]   random_seed      Random seed
 * \param [out]  ng_cell          Global number of cells
 * \param [out]  ng_face          Global number of faces
 * \param [out]  ng_vtx           Global number of vertices
 * \param [out]  n_face_group     Number of face groups
 * \param [out]  dn_cell          Local number of cells
 * \param [out]  dn_face          Local number of faces
 * \param [out]  dn_vtx           Local number of vertices
 * \param [out]  dcell_face_idx   Index of cell-face connectivity (size = \p dn_cell + 1)
 * \param [out]  dcell_face       Distributed cell-face connectivity (size = \p dcell_face_idx[\p dn_cell])
 * \param [out]  dface_cell       Distributed face-cell connectivity (size = 2 * \p dn_face)
 * \param [out]  dface_vtx_idx    Index of face-vertex connectivity (size = \p dn_face + 1)
 * \param [out]  dface_vtx        Distributed face-vertex connectivity (size = \p dface_vtx_idx[\p dn_face])
 * \param [out]  dvtx_coord       Coordinates of local vertices (size = 3 * \p dn_vtx)
 * \param [out]  dface_group_idx  Index of dface_group (size = \p n_face_group + 1)
 * \param [out]  dface_group      Distributed lists of faces in each group (size = \p dface_group_idx[\p n_face_group])
 *
 */

void
PDM_poly_vol_gen
(
 PDM_MPI_Comm  comm,
 double        xmin,
 double        ymin,
 double        zmin,
 double        lengthx,
 double        lengthy,
 double        lengthz,
 PDM_g_num_t   nx,
 PDM_g_num_t   ny,
 PDM_g_num_t   nz,
 int           randomize,
 int           random_seed,
 PDM_g_num_t  *ng_cell,
 PDM_g_num_t  *ng_face,
 PDM_g_num_t  *ng_vtx,
 int          *n_face_group,
 int          *dn_cell,
 int          *dn_face,
 int          *dn_vtx,
 int         **dcell_face_idx,
 PDM_g_num_t **dcell_face,
 PDM_g_num_t **dface_cell,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 double      **dvtx_coord,
 int         **dface_group_idx,
 PDM_g_num_t **dface_group
 );

#ifdef __cplusplus
}
#endif
#endif
