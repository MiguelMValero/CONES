/*
 * \file
 */

#ifndef __PDM_EXTRACT_PART_H__
#define __PDM_EXTRACT_PART_H__

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
#include "pdm_part_to_part.h"
#include "pdm_part_mesh_nodal_elmts.h"

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

typedef struct _pdm_extract_part_t PDM_extract_part_t;


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
 * \param [in]   n_part_in           Number of initial partition
 * \param [in]   n_part_out          Number of final partition
 * \param [in]   extract_kind        Extraction kind : (local/requilibrate/from target)
 * \param [in]   split_dual_method   Split method if requilibrate extract_kind
 * \param [in]   compute_child_gnum  Yes/No computation of a newest global numbering
 * \param [in]   ownership           Tell if you want ownership of resulting
 * \param [in]   comm                MPI communicator
 *
 * \return   Initialized \ref PDM_extract_part_t instance
 */
PDM_extract_part_t*
PDM_extract_part_create
(
 const int                     dim,
 const int                     n_part_in,
 const int                     n_part_out,
       PDM_extract_part_kind_t extract_kind,
       PDM_split_dual_t        split_dual_method,
       PDM_bool_t              compute_child_gnum,
       PDM_ownership_t         ownership,
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
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
);


/**
 *
 * \brief Set the extract number
 *
 * \param [in]   extrp         PDM_extract_part_t
 * \param [in]   i_part        part identifier
 * \param [in]   n_extract     Number of entity to select
 * \param [in]   extract_lnum  List of id to extract (starting at 0)
 *
 */
void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum
);

/**
 *
 * \brief Set the extract target number
 *
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   i_part            part identifier
 * \param [in]   n_target          Number of target to select
 * \param [in]   target_gnum       List of global id to extract
 * \param [in]   target_location   Init location (optional NULL pointer accepted and computed internaly)
 *
 */
void
PDM_extract_part_target_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_target,
  PDM_g_num_t              *target_gnum,
  int                      *target_location
);


/**
 *
 * \brief Keep target_gnum data ownership inside extrp
 *
 * \param [in]   extrp             PDM_extract_part_t
 * 
 */
void
PDM_extract_part_target_gnum_keep_ownnership
(
  PDM_extract_part_t       *extrp
);


/**
 *
 * \brief Set partition
 *
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   i_part            part identifier
 *
 */
void
PDM_extract_part_part_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);



/**
 *
 * \brief Set partition group (optional)
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   bound_type        Kind of group
 * \param [in]   n_group           Number of group of kind bound_type
 */
void
PDM_extract_part_n_group_set
(
  PDM_extract_part_t        *extrp,
  PDM_bound_type_t           bound_type,
  int                        n_group
);

/**
 *
 * \brief Set partition group (optional)
 *
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   i_part            part identifier
 * \param [in]   i_group           group identifier
 * \param [in]   bound_type        Kind of group
 * \param [in]   n_group_entity    Number of entity in current group
 * \param [in]   group_entity      List of entity in group (size = n_group_entity)
 *
 */
void
PDM_extract_part_part_group_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity,
  PDM_g_num_t              *group_entity_ln_to_gn
);

/**
 *
 * \brief Set PDM_part_mesh_nodal_elmts_t
 *
 * \param [in]   extrp            PDM_extract_part_t structure
 * \param [in]   pmne             PDM_part_mesh_nodal_elmts_t corresponding of dimenstion
 *
 */
void
PDM_extract_part_part_nodal_set
(
  PDM_extract_part_t          *extrp,
  PDM_part_mesh_nodal_elmts_t *pmne
);

/**
 *
 * \brief Set entity center (ussefull for equilibrate / hilbert ordering )
 *
 * \param [in]   extrp            PDM_extract_part_t structure
 * \param [in]   i_part           part identifier
 * \param [in]   entity_center    Center of entity (relative of dim throught create)
 *
 */
void
PDM_extract_part_entity_center_set
(
  PDM_extract_part_t          *extrp,
  int                          i_part,
  double                      *entity_center
);


/**
 *
 * \brief Set PDM_part_mesh_nodal_elmts_t
 *
 * \param [in]   extrp               PDM_extract_part_t structure
 * \param [in]   i_part_out          part identifier
 * \param [in]   PDM_mesh_entities_t Kind of entity required \ref PDM_mesh_entities_t
 *
 * \return Number of entity
 */
int
PDM_extract_part_n_entity_get
(
 PDM_extract_part_t       *extrp,
 int                       i_part_out,
 PDM_mesh_entities_t       entity_type
);



/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  extrp               Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part              Id of part
 * \param [in]  connectivity_type   Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect             Connectivity array (size = connect_idx[n_entity] )
 * \param [in]  connect_idx         Connectivity index (size = n_entity+1 )
 * \param [in]  ownership           Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
int
PDM_extract_part_connectivity_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp             Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part            Id of part
 * \param [in]  entity_type       Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_ln_to_gn   Entity local numbering to global numbering (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership         Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t        entity_type,
 PDM_g_num_t              **pentity_ln_to_gn,
 PDM_ownership_t            ownership
);


/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp                  Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part                 Id of part
 * \param [in]  entity_type            Entity kind \ref PDM_mesh_entities_t)
 * \param [out] parent_entity_ln_to_gn Entity local numbering to global numbering of parent entity, correspond to the input mesh (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership              Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_parent_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **parent_entity_ln_to_gn,
 PDM_ownership_t           ownership
);

/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp                  Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part                 Id of part
 * \param [in]  entity_type            Entity kind \ref PDM_mesh_entities_t)
 * \param [out] parent_entity_lnum     Local indexes of parent entity, correspond to the input mesh (size = n_entity)
 * \param [in]  ownership              Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_parent_lnum_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t        entity_type,
 int                      **parent_entity_lnum,
 PDM_ownership_t            ownership
);


/**
 *
 * \brief Get the vertex coordinates on current i_part partition and return number of vertices
 *
 * \param [in]   extrp      Pointer to \ref PDM_extract_part_t object
 * \param [in]   i_part     Id of part
 * \param [out]  vtx_coord  Vertex coordinate (size = 3 * n_vtx)
 * \param [in]   ownership  Ownership for color ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_vtx_coord_get
(
 PDM_extract_part_t         *extrp,
 int                        i_part_out,
 double                   **pvtx_coord,
 PDM_ownership_t            ownership
);


/**
 * \brief Retreive the partitionned mesh
 *
 * \param [in]  extrp             Pointer to \ref PDM_extract_part_t object
 * \param [out] extract_pmne      Partitionned mesh nodal, describe by elements (see \ref PDM_part_mesh_nodal_elmts_t )
 * \param [in]  ownership         Who is responsible to free retreived data ?
 *
 */
void
PDM_extract_part_part_mesh_nodal_get
(
  PDM_extract_part_t           *extrp,
  PDM_part_mesh_nodal_elmts_t **extract_pmne,
  PDM_ownership_t               ownership
);


/**
 *
 * \brief Free the structure
 *
 * \param [in]   extrp      Pointer to \ref PDM_extract_part_t object
 */
void
PDM_extract_part_free
(
  PDM_extract_part_t  *extrp
);

/**
 *
 * \brief Free all resulting array if not owner
 *
 * \param [in]   extrp      Pointer to \ref PDM_extract_part_t object
 */
void
PDM_extract_part_partial_free
(
  PDM_extract_part_t  *extrp
);


/**
 *
 * \brief Get for entity_type (cells/faces/edge/vertices) the associated part_to_part (\ref PDM_part_to_part_t ). The part to part exchange protocol allow user to
 * exchange easily data from input mesh to the extract one.
 *
 * \param [in]   extrp        Pointer to \ref PDM_extract_part_t object
 * \param [in]   entity_type  Bound type \ref PDM_mesh_entities_t
 * \param [out]  ptp          Part to part protocol exchange, to exchange betwenn the input mesh and the output one (\ref PDM_part_to_part_t)
 * \param [in]   ownership    Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_part_to_part_get
(
       PDM_extract_part_t   *extrp,
 const PDM_mesh_entities_t   entity_type,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership

);


/**
 *
 * \brief Get for bound_type the associated part_to_part (\ref PDM_part_to_part_t ). The part to part exchange protocol allow user to
 * exchange easily data from input mesh to the extract one.
 *
 * \param [in]   extrp        Pointer to \ref PDM_extract_part_t object
 * \param [in]   bound_type   Bound type \ref PDM_bound_type_t
 * \param [in]   i_group      Id of group
 * \param [out]  ptp          Part to part protocol exchange, to exchange betwenn the input mesh and the output one (\ref PDM_part_to_part_t)
 * \param [in]   ownership    Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_part_to_part_group_get
(
       PDM_extract_part_t   *extrp,
 const PDM_bound_type_t      bound_type,
       int                   i_group,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership

);


/**
 *
 * \brief Get the bound description for the entity (cell/face/edge/vertices)
 *
 * \param [in]   extrp                              Pointer to \ref PDM_extract_part_t object
 * \param [in]   bound_type                             Bound type \ref PDM_bound_type_t
 * \param [in]   i_part                                 Id of part
 * \param [in]   i_group                                Id of group
 * \param [out]  pn_extract_group_entity                Number of entity in current group
 * \param [out]  pextract_group_entity_ln_to_gn         Entity global id in current partition (size = pn_extract_group_entity)
 * \param [out]  pextract_group_entity_parent_ln_to_gn  Entity global id in parent partition (size = pn_extract_group_entity)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_group_get
(
       PDM_extract_part_t   *extrp,
 const PDM_bound_type_t      bound_type,
       int                   i_part,
       int                   i_group,
       int                  *pn_extract_group_entity,
       int                 **pextract_group_entity,
       PDM_g_num_t         **pextract_group_entity_ln_to_gn,
       PDM_g_num_t         **pextract_group_entity_parent_ln_to_gn,
       PDM_ownership_t       ownership
);




/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXTRACT_PART_H__ */
