/*
 * \file
 */

#ifndef __PDM_PART_MESH_H__
#define __PDM_PART_MESH_H__

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

typedef struct _pdm_part_mesh_t PDM_part_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   dn_bnd              Number of boundaries
 * \param [in]   n_join              Number of interfaces with other domains
 *
 * \return     Identifier
 */

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
);

void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *dmesh
);


void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       pn_entity
);


int
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type
);

void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                      *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_vtx_coord_set
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double            *vtx_coord,
 PDM_ownership_t    ownership
);


void
PDM_part_mesh_vtx_coord_get
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
);

void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 int                        i_part,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_entity_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t              *pentity_ln_to_gn,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_entity_ln_to_gn_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **pentity_ln_to_gn,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_entity_color_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                      *pentity_color,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_entity_color_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                     **pentity_color,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_n_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound
);

int
PDM_part_mesh_n_bound_get
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type
);

int
PDM_part_mesh_tn_part_get
(
 PDM_part_mesh_t          *pmesh
);

int
PDM_part_mesh_n_part_get
(
 PDM_part_mesh_t          *pmesh
);

void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                       pn_bound,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                      *pn_bound,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_bound_concat_compute
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type
);

void
PDM_part_mesh_bound_concat_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                      *pbound_idx,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_bound_concat_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **pbound_idx,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_part_graph_comm_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **ppart_bound_proc_idx,
 int                     **ppart_bound_part_idx,
 int                     **ppart_bound,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_part_graph_comm_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                      *ppart_bound_proc_idx,
 int                      *ppart_bound_part_idx,
 int                      *ppart_bound,
 PDM_ownership_t           ownership
);

/**
 * \brief Export a partitioned mesh in Ensight format
 *
 * \param [in] pmesh          Pointer to \ref PDM_part_mesh_t object
 * \param [in] directory      Output directory
 * \param [in] name           Output name
 * \param [in] export_bounds  Option to export bounds
 *
 */

void
PDM_part_mesh_dump_ensight
(
 PDM_part_mesh_t *pmesh,
 const char      *directory,
 const char      *name,
 PDM_bool_t       export_bounds
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_H__ */
