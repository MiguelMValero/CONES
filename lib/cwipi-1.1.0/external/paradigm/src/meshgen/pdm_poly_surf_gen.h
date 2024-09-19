/*
 * \file
 */

#ifndef __PDM_POLY_SURF_GEN_H__
#define __PDM_POLY_SURF_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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
 * \brief Generate a distributed polygonal surface mesh
 *
 * \param [in]   pdm_comm         MPI communicator
 * \param [in]   xmin             Minimum x-coordinate
 * \param [in]   xmax             Maximum x-coordinate
 * \param [in]   ymin             Minimum y-coordinate
 * \param [in]   ymax             Maximum y-coordinate
 * \param [in]   have_random      Enable/disable randomization
 * \param [in]   init_random      Random seed
 * \param [in]   nx               Number of vertices in the x-direction
 * \param [in]   ny               Number of vertices in the y-direction
 * \param [out]  ng_face          Global number of faces
 * \param [out]  ng_vtx           Global number of vertices
 * \param [out]  ng_edge          Global number of edges
 * \param [out]  dn_vtx           Local number of vertices
 * \param [out]  dvtx_coord       Coordinates of local vertices (size = 3 * \ref dn_vtx)
 * \param [out]  dn_face          Local number of faces
 * \param [out]  dface_vtx_idx    Index of face-vertex connectivity (size = \ref dn_face + 1)
 * \param [out]  dface_vtx        Distributed face-vertex connectivity (size = \ref dface_vtx_idx[\ref dn_face])
 * \param [out]  dn_edge          Local number of edges
 * \param [out]  dedge_vtx        Distributed edge-vertex connectivity (size = 2 * \ref dn_edge)
 * \param [out]  dedge_face       Distributed edge-face connectivity (size = 2 * \ref dn_edge)
 * \param [out]  n_edge_group     Number of edge groups
 * \param [out]  dedge_group_idx  Index of dedge_group (size = \ref n_edge_group + 1)
 * \param [out]  dedge_group      Distributed lists of edges in each group (size = \ref dedge_group_idx[\ref n_edge_group])
 *
 */

void PDM_poly_surf_gen
(
PDM_MPI_Comm     pdm_comm,
double           xmin,
double           xmax,
double           ymin,
double           ymax,
int              have_random,
int              init_random,
PDM_g_num_t      nx,
PDM_g_num_t      ny,
PDM_g_num_t     *ng_face,
PDM_g_num_t     *ng_vtx,
PDM_g_num_t     *ng_edge,
int             *dn_vtx,
double         **dvtx_coord,
int             *dn_face,
int            **dface_vtx_idx,
PDM_g_num_t    **dface_vtx,
PDM_g_num_t    **dface_edge,
int             *dn_edge,
PDM_g_num_t    **dedge_vtx,
PDM_g_num_t    **dedge_face,
int             *n_edge_group,
int            **dedge_group_idx,
PDM_g_num_t    **dedge_group
);

#ifdef __cplusplus
}
#endif
#endif
