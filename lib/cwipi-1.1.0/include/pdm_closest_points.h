/*
 * \file
 */

#ifndef PDM_CLOSEST_POINTS_H
#define PDM_CLOSEST_POINTS_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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
#include "pdm_part_to_part.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_closest_point_t PDM_closest_point_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to look for the closest points of a point cloud
 * (target cloud) in an other point cloud (source cloud)
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_closest      Number of closest source points to find for each
 *                              target point
 * \param [in ] owner           Ownership for \ref PDM_closest_point_t
 *
 * \return     Pointer to \ref PDM_closest_points_t object
 *
 */

PDM_closest_point_t*
PDM_closest_points_create
(
 const PDM_MPI_Comm    comm,
 const int             n_closest,
 const PDM_ownership_t owner
);


/**
 *
 * \brief Set the number of partitions of both point clouds
 *
 * \param [in]   cls               Pointer to \ref PDM_closest_points_t object
 * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
 * \param [in]   n_part_cloud_tgt  Number of partitions of the target cloud
 *
 */

void
PDM_closest_points_n_part_cloud_set
(
       PDM_closest_point_t* cls,
 const int                  n_part_cloud_src,
 const int                  n_part_cloud_tgt
);


/**
 *
 * \brief Set the target point cloud
 *
 * \param [in]   cls             Pointer to \ref PDM_closest_points_t object
 * \param [in]   i_part          Partition identifier
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates (size : 3 * \p n_points)
 * \param [in]   gnum            Point global ids (size : \p n_points)
 *
 */

void
PDM_closest_points_tgt_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);


/**
 *
 * \brief Set the source point cloud
 *
 * \param [in]   cls             Pointer to \ref PDM_closest_points_t object
 * \param [in]   i_part          Partition identifier
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates (size : 3 * \p n_points)
 * \param [in]   gnum            Point global ids (size : \p n_points)
 *
 */

void
PDM_closest_points_src_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);

/**
 *
 * \brief Look for closest points
 *
 * \param [in]   cls Pointer to \ref PDM_closest_points_t object
 *
 */

void
PDM_closest_points_compute
(
 PDM_closest_point_t *cls
);


/**
 *
 * \brief Get closest source points global ids and (squared) distance
 *
 * \param [in]   cls                   Pointer to \ref PDM_closest_points_t object
 * \param [in]   i_part_tgt            Partition identifier
 * \param [out]  closest_src_gnum      Global ids of the closest source points (size : *n_closest* * *n_tgt_points*)
 * \param [out]  closest_src_distance  Squared distance from closest source points (size : *n_closest* * *n_tgt_points*)
 *
 */

void
PDM_closest_points_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_tgt,
       PDM_g_num_t         **closest_src_gnum,
       double              **closest_src_distance
);


/**
 *
 * \brief Free a \ref PDM_closest_points_t instance
 *
 * \param [in]  cls      Pointer to \ref PDM_closest_points_t object
 *
 */

void
PDM_closest_points_free
(
 PDM_closest_point_t  *cls
);


/**
 *
 * \brief  Dump elapsed and CPU times
 *
 * \param [in]  cls      Pointer to \ref PDM_closest_points_t object
 *
 */

void
PDM_closest_points_dump_times
(
 PDM_closest_point_t  *cls
);

/**
 *
 * \brief Get source-to-target mapping
 *
 * \param [in]   cls                Pointer to \ref PDM_closest_points_t object
 * \param [in]   i_part_src         Partition identifier
 * \param [out]  src_to_tgt_idx     Index for source-to-target mapping (size : *n_src_points* + 1)
 * \param [out]  src_to_tgt         Source-to-target mapping (global ids) (size : \p src_to_tgt_idx[ *n_src_points* ] )
 *
 */

void
PDM_closest_points_tgt_in_src_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_src,
       int                 **src_to_tgt_idx,
       PDM_g_num_t         **src_to_tgt
);
/**
 *
 * \brief Get source-to-target distance
 *
 * \param [in]   cls                Pointer to \ref PDM_closest_points_t object
 * \param [in]   i_part_src         Partition identifier
 * \param [out]  src_to_tgt_idx     Index for source-to-target mapping (size : *n_src_points* + 1)
 * \param [out]  src_to_tgt_dist    Squared source-to-target distance (size : \p src_to_tgt_idx[ *n_src_points* ] )
 *
 */
void
PDM_closest_points_tgt_in_src_dist_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_src,
       int                 **src_to_tgt_idx,
       double              **src_to_tgt_dist
);



/**
 *
 * \brief WTF?
 *
 */
void
PDM_transform_to_parent_gnum
(
 const int           n_part_initial,
 const int          *n_elmt_initial,
 const PDM_g_num_t **child_ln_to_gn,
 const PDM_g_num_t **parent_ln_to_gn,
 const int           n_part_to_transform,
 const int          *n_elmt_to_transform,
       PDM_g_num_t **gnum_to_transform,
       PDM_MPI_Comm  comm
);


/**
 *
 * \brief  transfert _closest_pts var as it seems this static var is not readable
 *          when we switch to the nvcc compiler
 *
 * **??**
 *
 */

PDM_closest_point_t*
PDM_closest_points_closest_transfert
(
  PDM_closest_point_t  *cls
);


/**
 *
 * \brief  Get the number of target points in a partition
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points_t object
 * \param [in]  i_part  Partition identifier
 *
 * \return   Number of target points in partition \p i_part
 *
 */

int
PDM_closest_points_n_tgt_get
(
  PDM_closest_point_t  *cls,
  const int             i_part
);


/**
 *
 * \brief  Get the number of source points in a partition
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points_t object
 * \param [in]  i_part  Partition identifier
 *
 * \return   Number of source points in partition \p i_part
 *
 */
int
PDM_closest_points_n_src_get
(
  PDM_closest_point_t  *cls,
  const int             i_part
);


/**
 *
 * \brief  Get the number of closest points
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points_t object
 *
 * \return   Number of closest points
 *
 */

int
PDM_closest_points_n_closest_get
(
  PDM_closest_point_t  *cls
);

/**
 * \brief Get \ref PDM_part_to_part_t object to exchange data between
 * the source and target point clouds (both in user frame)
 *
 * \param [in ] cls        Pointer to \ref PDM_closest_point_t object
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for \p ptp
 *
 */

void
PDM_closest_points_part_to_part_get
(
 PDM_closest_point_t  *cls,
 PDM_part_to_part_t  **ptp,
 PDM_ownership_t       ownership
 );

#ifdef	__cplusplus
}
#endif

#endif // PDM_CLOSEST_POINTS_H
