#ifndef __PDM_DMESH_EXTRACT_H__
#define __PDM_DMESH_EXTRACT_H__

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
#include "pdm_block_to_part.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"

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

typedef struct _pdm_dmesh_extract_t PDM_dmesh_extract_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build an extract_part struct
 *
 * \param [in]   dim                 Extraction dimension
 * \param [in]   comm                MPI communicator
 *
 * \return   Initialized \ref PDM_extract_part_t instance
 */
PDM_dmesh_extract_t*
PDM_dmesh_extract_create
(
 const int                     dim,
       PDM_MPI_Comm            comm
);


/**
 *
 * \brief Compute extraction
 *
 * \param [in]   extrp      PDM_extract_part_t
 *
 */
void
PDM_dmesh_extract_compute
(
 PDM_dmesh_extract_t *dme
);



/**
 *
 * \brief Set the extract number
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   entity_type          Entity kind to be extracted (\ref PDM_mesh_entities_t)
 * \param [in]   n_selected           Number of entity to select
 * \param [in]   selected_gnum        List of global id to extract
 *
 */
void
PDM_dmesh_extract_selected_gnum_set
(
 PDM_dmesh_extract_t *dme,
 PDM_mesh_entities_t  entity_type,
 int                  n_selected,
 PDM_g_num_t         *selected_gnum
);

/**
 *
 * \brief Set the dn_entity of entity_type
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   entity_type          Entity kind (\ref PDM_mesh_entities_t)
 * \param [in]   dn_entity            Number of entity in current process
 *
 */
void
PDM_dmesh_extract_dn_entity_set
(
 PDM_dmesh_extract_t *dme,
 PDM_mesh_entities_t  entity_type,
 int                  dn_entity
);


/**
 *
 * \brief Set vertices coordinates
 *
 * \param [in]   dme           PDM_dmesh_extract_t
 * \param [in]   dvtx_coord    Distributed vertex coordinates (size = 3 * dn_vtx )
 *
 */
void
PDM_dmesh_extract_vtx_coord_set
(
 PDM_dmesh_extract_t *dme,
 double              *dvtx_coord
);

/**
 *
 * \brief Set mesh bound for one bound_type
 *
 * \param [in]   dme           PDM_dmesh_extract_t
 * \param [in]   bound_type    Bound kind (\ref PDM_bound_type_t)
 * \param [in]   n_bound       Number of bound in for bound_type
 * \param [in]   connect       Connectivity between group and entity (size = connect_idx[n_bound])
 * \param [in]   connect_idx   Connectivity index between group and entity (size = n_bound+1)
 *
 */
void
PDM_dmesh_extract_dmesh_bound_set
(
 PDM_dmesh_extract_t *dme,
 PDM_bound_type_t     bound_type,
 int                  n_bound,
 PDM_g_num_t         *connect,
 int                 *connect_idx
);


/**
 *
 * \brief Set connectivity by kind (\ref PDM_connectivity_type_t )
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   connectivity_type    Connectivity kind (\ref PDM_connectivity_type_t)
 * \param [in]   dconnect             Connectivity (size = dconnect_idx[dn_entity])
 * \param [in]   dconnect_idx         Connectivity (size = dn_entity+1)
 *
 */
void
PDM_dmesh_extract_dconnectivity_set
(
       PDM_dmesh_extract_t     *dme,
       PDM_connectivity_type_t  connectivity_type,
       PDM_g_num_t             *dconnect,
       int                     *dconnect_idx
);


/**
 *
 * \brief Set a dmesh object correspond to extraction
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [out]  dmesh                PDM_dmesh_t who need to be extracted
 *
 */
void
PDM_dmesh_extract_dmesh_set
(
 PDM_dmesh_extract_t     *dme,
 PDM_dmesh_t             *dmesh
);

/**
 *
 * \brief Get a dmesh object correspond to extraction
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [out]  dmesh_nodal          PDM_dmesh_nodal_t who need to be extracted
 *
 */
void
PDM_dmesh_extract_dmesh_nodal_set
(
 PDM_dmesh_extract_t     *dme,
 PDM_dmesh_nodal_t       *dmesh_nodal
);

/**
 *
 * \brief Get a dmesh object correspond to extraction
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [out]  dmesh_extract        Current extraction direclty inside a PDM_dmesh_t
 * \param [in]   ownership            KEEP or USER
 *
 */
void
PDM_dmesh_extract_dmesh_get
(
 PDM_dmesh_extract_t     *dme,
 PDM_dmesh_t            **dmesh_extract,
 PDM_ownership_t          ownership
);

/**
 *
 * \brief Get a dmesh object correspond to extraction
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [out]  dmesh_nodal_extract  Current extraction direclty inside a PDM_dmesh_nodal_t
 * \param [in]   ownership            KEEP or USER
 *
 */
void
PDM_dmesh_extract_dmesh_nodal_get
(
 PDM_dmesh_extract_t     *dme,
 PDM_dmesh_nodal_t      **dmesh_nodal_extract,
 PDM_ownership_t          ownership
);


/**
 *
 * \brief Get the redistributed parent_gnum (in block frame)
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   entity_type          Entity type (cell, face, edge, vtx)
 * \param [out]  dn_entity            Size of block of current entity
 * \param [out]  parent_gnum          Parent gnum redistributed
 * \param [in]   ownership            KEEP or USER
 *
 */
void
PDM_dmesh_extract_parent_gnum_get
(
 PDM_dmesh_extract_t     *dme,
 PDM_mesh_entities_t      entity_type,
 int                     *dn_entity,
 PDM_g_num_t            **parent_gnum,
 PDM_ownership_t          ownership
);


/**
 *
 * \brief Get the block_to_part associated to extraction
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   entity_type          Entity type (cell, face, edge, vtx)
 * \param [out]  btp                  block_to_part to transfert data to extract block
 * \param [in]   ownership            KEEP or USER
 *
 */
void
PDM_dmesh_extract_btp_get
(
 PDM_dmesh_extract_t     *dme,
 PDM_mesh_entities_t      entity_type,
 PDM_block_to_part_t    **btp,
 PDM_ownership_t          ownership
);


/**
 *
 * \brief Get the block_to_part associated to extraction for each group
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 * \param [in]   i_group              No of group
 * \param [in]   bound_type           Bound type (cell, face, edge, vtx)
 * \param [out]  btp                  block_to_part to transfert data to extract block
 * \param [in]   ownership            KEEP or USER
 *
 */
void
PDM_dmesh_extract_btp_group_get
(
 PDM_dmesh_extract_t     *dme,
 int                      i_group,
 PDM_bound_type_t         bound_type,
 PDM_block_to_part_t    **btp,
 PDM_ownership_t          ownership
);


/**
 *
 * \brief Free structure
 *
 * \param [in]   dme                  PDM_dmesh_extract_t
 *
 */
void
PDM_dmesh_extract_free
(
  PDM_dmesh_extract_t  *dme
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_EXTRACT_H__ */
