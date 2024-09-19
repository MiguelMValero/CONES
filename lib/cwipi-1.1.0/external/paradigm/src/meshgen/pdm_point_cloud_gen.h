/*
 * \file
 */

#ifndef __PDM_POINT_CLOUD_GEN_H__
#define __PDM_POINT_CLOUD_GEN_H__

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
 * \brief Generate a uniformly random point cloud inside a cuboid.
 *
 * \param [in]   comm                   MPI Communicator
 * \param [in]   seed                   Random seed
 * \param [in]   gn_pts                 Global number of points in the cloud
 * \param [in]   geometric_g_num        Compute global ids from coordinates (Morton ordering)
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  ln_pts                 Local number of points in the cloud
 * \param [out]  coord                  XYZ-coordinates of the local points
 * \param [out]  g_num                  Global ids of the local points
 *
 */

void
PDM_point_cloud_gen_random
(
 PDM_MPI_Comm        comm,
 const int           seed,
 const int           geometric_g_num,
 const PDM_g_num_t   gn_pts,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 int                *ln_pts,
 double            **coord,
 PDM_g_num_t       **g_num
 );


/**
 *
 * \brief Generate a cartesian point cloud inside a cuboid.
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
 * \param [out]  n_pts                  Local number of points in the cloud
 * \param [out]  pts_coord              XYZ-coordinates of the local points
 * \param [out]  pts_ln_to_gn           Global ids of the local points
 *
 */

void
PDM_point_cloud_gen_cartesian
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
 int                *n_pts,
 double            **pts_coord,
 PDM_g_num_t       **pts_ln_to_gn
);



/**
 *
 * \brief Generate a uniformly random point cloud inside a cuboid.
 *
 * \param [in]   comm                   MPI Communicator
 * \param [in]   seed                   Random seed
 * \param [in]   gn_pts                 Global number of points in the cloud
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  dpts_coord             XYZ-coordinates of the local points
 * \param [out]  distrib_pts            Distribution of points
 *
 */

void
PDM_dpoint_cloud_gen_random
(
 PDM_MPI_Comm        comm,
 const int           seed,
 const PDM_g_num_t   gn_pts,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 double            **dpts_coord,
 PDM_g_num_t       **distrib_pts
);

/**
 *
 * \brief Generate a cartesian point cloud inside a cuboid.
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
 * \param [out]  dpts_coord             XYZ-coordinates of the local points
 * \param [out]  distrib_pts            Distribution of points
 *
 */

void
PDM_dpoint_cloud_gen_cartesian
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
 double            **dpts_coord,
 PDM_g_num_t       **distrib_pts
);

#ifdef __cplusplus
}
#endif

#endif // PDM_POINT_CLOUD_GEN_H
