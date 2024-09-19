#ifndef __PDM_DMESH_EXTRACT_PRIV_H__
#define __PDM_DMESH_EXTRACT_PRIV_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_block_to_part.h"

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

/**
 * \struct _pdm_extract_part_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_dmesh_extract_t
{

  int                     dim;
  PDM_MPI_Comm            comm;

  /* Input */
  int                     n_selected;
  PDM_g_num_t            *selected_gnum;

  PDM_dmesh_t            *dmesh;
  PDM_dmesh_t            *dmesh_shared;
  PDM_dmesh_nodal_t      *dmesh_nodal;

  PDM_dmesh_t            *dmesh_extract;
  PDM_dmesh_nodal_t      *dmesh_nodal_extract;
  PDM_ownership_t         dmesh_extract_ownership;

  PDM_g_num_t            *distrib_extract              [PDM_MESH_ENTITY_MAX];
  PDM_g_num_t            *parent_extract_gnum          [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t         distrib_extract_ownership    [PDM_MESH_ENTITY_MAX];
  PDM_ownership_t         parent_extract_gnum_ownership[PDM_MESH_ENTITY_MAX];

  PDM_block_to_part_t    *btp_entity_to_extract_entity[PDM_MESH_ENTITY_MAX];
  PDM_ownership_t         btp_ownership               [PDM_MESH_ENTITY_MAX];

  int                     n_bound                           [PDM_BOUND_TYPE_MAX];
  PDM_block_to_part_t   **btp_bound_entity_to_extract_entity[PDM_BOUND_TYPE_MAX];
  PDM_ownership_t        *btp_bound_ownership               [PDM_BOUND_TYPE_MAX];


};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_EXTRACT_PRIV_H__ */
