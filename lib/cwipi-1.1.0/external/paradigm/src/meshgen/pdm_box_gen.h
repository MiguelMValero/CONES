/*
 * \file
 */

#ifndef __PDM_BOX_GEN_H__
#define __PDM_BOX_GEN_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2022-2023       ONERA

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Generate a random set of boxes
 *
 * \param [in]   comm                   MPI Communicator
 * \param [in]   seed                   Random seed
 * \param [in]   geometric_g_num        Compute global ids from coordinates
 * \param [in]   gn_box                 Global number of boxes
 * \param [in]   min_size               Minimal box size
 * \param [in]   max_size               Maximal box size
 * \param [in]   x_min                  Minimal X-coordinate for box centers
 * \param [in]   y_min                  Minimal Y-coordinate for box centers
 * \param [in]   z_min                  Minimal Z-coordinate for box centers
 * \param [in]   x_max                  Maximal X-coordinate for box centers
 * \param [in]   y_max                  Maximal Y-coordinate for box centers
 * \param [in]   z_max                  Maximal Z-coordinate for box centers
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes (size : 6 * \p n_box)
 * \param [out]  box_ln_to_gn           Global ids of the local boxes (size : \p n_box)
 *
 */

void
PDM_box_gen_random
(
 PDM_MPI_Comm   comm,
 int            seed,
 int            geometric_g_num,
 PDM_g_num_t    gn_box,
 double         min_size,
 double         max_size,
 double         x_min,
 double         y_min,
 double         z_min,
 double         x_max,
 double         y_max,
 double         z_max,
 int           *n_box,
 double       **box_extents,
 PDM_g_num_t  **box_ln_to_gn
 );



/**
 *
 * \brief Generate a cartesian set of boxes
 *
 * \param [in]   comm                   MPI Communicator
 * \param [in]   nx                     Number of points in X-direction
 * \param [in]   ny                     Number of points in Y-direction
 * \param [in]   nz                     Number of points in Z-direction
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes (size : 6 * \p n_box)
 * \param [out]  box_ln_to_gn           Global ids of the local boxes (size : \p n_box)
 *
 */

void
PDM_box_gen_cartesian
(
 PDM_MPI_Comm        comm,
 const int           nx,
 const int           ny,
 const int           nz,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 int                *n_box,
 double            **box_extents,
 PDM_g_num_t       **box_ln_to_gn
);

#ifdef __cplusplus
}
#endif

#endif // PDM_BOX_GEN_H
