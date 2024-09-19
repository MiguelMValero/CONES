/*
 * \file
 */

#ifndef __PDM_MULTIPART_H__
#define __PDM_MULTIPART_H__

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
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_domain_interface.h"

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

typedef struct _pdm_multipart_t PDM_multipart_t;

/**
 * \enum PDM_part_size_t
 * \brief Use homogeneous or heterogeneous partition sizes (only for ParMetis method)
 */
typedef enum {
  PDM_PART_SIZE_HOMOGENEOUS   = 1, /*!< All requested partition have the same size */
  PDM_PART_SIZE_HETEROGENEOUS = 2, /*!< Each requested partition can have a portion (within 0. and 1.) of the mesh */
} PDM_part_size_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure. This method allows to split multiple domains
 *
 * \param [in]   n_domain         Number of domains in the original mesh
 * \param [in]   n_part           Number of partition per rank in each domain
 * \param [in]   merge_domains    Merge or not the domains before splitting
 * \param [in]   split_method     Choice of library used to split the mesh
 * \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
 * \param [in]   part_fraction    Weight (in %) of each partition in heterogeneous case (i.e. if \p part_size_method is set to \p PDM_PART_SIZE_HETEROGENEOUS)
 * \param [in]   comm             PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_multipart_t instance
 */

PDM_multipart_t *
PDM_multipart_create
(
 const int              n_domain,
 const int             *n_part,
 const PDM_bool_t       merge_domains,
 const PDM_split_dual_t split_method,
 const PDM_part_size_t  part_size_method,
 const double          *part_fraction,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
);


/**
 *
 * \brief Set distributed mesh data for the input domain
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain       Domain identifier
 * \param [in]   dmesh          Pointer on \ref PDM_dmesh_t containing all distributed connectivities
 */

void PDM_multipart_dmesh_set
(
 PDM_multipart_t   *multipart,
 const int          domain_id,
       PDM_dmesh_t *dmesh
);

/**
 *
 * \brief Set distributed mesh data for the input domain. The mesh is described by nodal connectivity
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain       Domain identifier
 * \param [in]   dmesh_nodal    Pointer on \ref PDM_dmesh_nodal_t
 */

void PDM_multipart_dmesh_nodal_set
(
 PDM_multipart_t         *multipart,
 const int                domain_id,
       PDM_dmesh_nodal_t *dmesh_nodal
);

/**
 * \brief Set block
 *
 * \param [in]   multipart              Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain               Domain identifier
 * \param [in]   dn_cell                Number of distributed cells
 * \param [in]   dn_face                Number of distributed faces
 * \param [in]   dn_vtx                 Number of distributed vertices
 * \param [in]   n_face_group           Number of face groups
 * \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
 *                                      (size : \p dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face             Distributed cell face connectivity or NULL
 *                                      (size : \p dface_vtx_idx[\p dn_cell], numbering : 1 to n)
 * \param [in]   dface_cell             Distributed face cell connectivity or NULL
 *                                      (size : 2 * \p dn_face, numbering : 1 to n)
 * \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
 *                                      (size : \p dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx              Distributed face to vertex connectivity
 *                                      (size : \p dface_vtx_idx[\p dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord             Distributed vertex coordinates
 *                                      (size : 3 * \p dn_vtx)
 * \param [in]   dface_group_idx        Index of distributed faces list of each group
 *                                      (size = \p n_face_group + 1) or NULL
 * \param [in]   dface_group            Distributed faces list of each group
 *                                      (size = \p dface_group[\p dface_group_idx[\p n_face_group]], numbering : 1 to n)
 *                                      or NULL
 *
 */
void
PDM_multipart_block_set
(
 PDM_multipart_t             *multipart,
 const int                    i_domain,
 const int                    dn_cell,
 const int                    dn_face,
 const int                    dn_vtx,
 const int                    n_face_group,
 const int                   *dcell_face_idx,
 const PDM_g_num_t           *dcell_face,
 const PDM_g_num_t           *dface_cell,
 const int                   *dface_vtx_idx,
 const PDM_g_num_t           *dface_vtx,
 const double                *dvtx_coord,
 const int                   *dface_group_idx,
 const PDM_g_num_t           *dface_group
);

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain              Id of domain which parameters apply (or -1 for all domains)
 * \param [in]   renum_cell_method     Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cache-blocking method :
 *                                     [*n_cell_per_cache_wanted*, *is_asynchronous*, *is_vectorisation*, *n_vect_face*, *split_method*]
 * \param [in]   renum_face_method     Choice of renumbering method for faces
 *
 */
void PDM_multipart_set_reordering_options
(
 PDM_multipart_t *multipart,
 const int        i_domain,
 const char      *renum_cell_method,
 const int       *renum_cell_properties,
 const char      *renum_face_method
);

/**
 *
 * \brief Set the vertex reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain              Id of domain which parameters apply (or -1 for all domains)
 * \param [in]   renum_vtx_method      Choice of renumbering method for vertices
 *
 */

void PDM_multipart_set_reordering_options_vtx
(
 PDM_multipart_t *multipart,
 const int        i_domain,
 const char      *renum_vtx_method
);


/**
 *
 * \brief Construct the partitioned meshes on all domains
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t instance
 */
void
PDM_multipart_compute
(
 PDM_multipart_t *multipart
);


/**
 * \brief Retrieve the partitioned nodal mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [out] pmesh_nodal           Partitioned nodal mesh
 * \param [in]  ownership             Who is responsible to free retrieved data ?
 *
 */

void
PDM_multipart_get_part_mesh_nodal
(
      PDM_multipart_t        *multipart,
const int                     i_domain,
      PDM_part_mesh_nodal_t **pmesh_nodal,
      PDM_ownership_t         ownership
);

/**
 * \brief Retrieve the partitioned mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [out] pmesh                 Partitioned mesh
 *
 */

// void
// PDM_multipart_get_part_mesh
// (
//        PDM_multipart_t  *multipart,
//  const int               i_domain,
//        PDM_part_mesh_t **pmesh
// );

/**
 * \brief Specify interface between domain (see \ref PDM_multidomain_interface)
 *
 * \param [in]  multipart   Pointer to \ref PDM_multipart_t instance
 * \param [in]  ditrf       Pointer to \ref PDM_domain_interface_t instance
 *
 */

void
PDM_multipart_domain_interface_shared_set
(
  PDM_multipart_t        *multipart,
  PDM_domain_interface_t *ditrf
);

/**
 *
 * \brief Returns the dimensions of a given partition
 */
void
PDM_multipart_part_dim_get
(
PDM_multipart_t *multipart,
const int        i_domain,
const int        i_part,
      int       *n_cell,
      int       *n_face,
      int       *n_face_part_bound,
      int       *n_vtx,
      int       *n_proc,
      int       *n_total_part,
      int       *s_cell_face,
      int       *s_face_vtx,
      int       *s_face_bound,
      int       *n_bound_groups
);


/**
 *
 * \brief Get the connection graph between partition for the requested entity type
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Type of mesh entity
 * \param [out] ppart_bound_proc_idx  Partitioning boundary entities index from process (size = n_proc + 1)
 * \param [out] ppart_bound_part_idx  Partitioning boundary entities index from partition (size = n_total_part + 1)
 * \param [out] ppart_bound           Partitioning boundary entities (size = 4 * n_entity_part_bound)
 * \param [in]  ownership             Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_multipart_part_graph_comm_get
(
 PDM_multipart_t      *multipart,
 const int             i_domain,
 const int             i_part,
 PDM_mesh_entities_t   entity_type,
 int                 **ppart_bound_proc_idx,
 int                 **ppart_bound_part_idx,
 int                 **ppart_bound,
 PDM_ownership_t       ownership
);

/**
 *
 * \brief Returns the data arrays of a given partition
 *
 * \deprecated Use \ref PDM_multipart_part_connectivity_get instead
 */
void
PDM_multipart_part_val_get
(
PDM_multipart_t     *multipart,
const int            i_domain,
const int            i_part,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_bound_idx,
      int          **face_bound,
      PDM_g_num_t  **face_bound_ln_to_gn
);


/**
 *
 * \brief Returns the total number of part among all process
 */
int
PDM_multipart_part_tn_part_get
(
PDM_multipart_t                *multipart,
const int                       i_domain
);

/**
 * \brief Get a partitioned connectivity
 *
 * \note If the return \p connect is \p NULL, you may build the missing connectivity using the appropriate **Connectivity transformation** function.
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  connectivity_type     Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect_idx           Connectivity index (size = *n_entity* + 1 )
 * \param [in]  connect               Connectivity array (size = \p connect_idx[*n_entity*] )
 * \param [in]  ownership             Choice of ownership of the resulting arrays \ref PDM_ownership_t
 *
 * \return Number of leading entities
 */
int
PDM_multipart_part_connectivity_get
(
PDM_multipart_t                *multipart,
const int                       i_domain,
const int                       i_part,
      PDM_connectivity_type_t   connectivity_type,
      int                     **connect_idx,
      int                     **connect,
      PDM_ownership_t           ownership
);


/**
 * \brief Get the number of entities with given type.
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t
 *
 * \return Number of entities
 */
int
PDM_multipart_part_n_entity_get
(
PDM_multipart_t            *multipart,
const int                   i_domain,
const int                   i_part,
      PDM_mesh_entities_t   entity_type
);

/**
 *
 * \brief Get the global ids of entities with given type.
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind
 * \param [out] entity_ln_to_gn       Entity local numbering to global numbering (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership             Ownership for \p entity_ln_to_gn
 *
 * \return Number of entities
 */
int
PDM_multipart_part_ln_to_gn_get
(
PDM_multipart_t            *multipart,
const int                   i_domain,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      PDM_g_num_t         **entity_ln_to_gn,
      PDM_ownership_t       ownership
);


/**
 *
 * \brief Get the color of entities with given type.
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind
 * \param [out] entity_color          Entity color (only for specific renumbering option )
 * \param [in]  ownership             Ownership for \p entity_color
 *  *
 * \return Number of entities
 */
int
PDM_multipart_partition_color_get
(
PDM_multipart_t            *multipart,
const int                   i_domain,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      int                 **entity_color,
      PDM_ownership_t       ownership
);

/**
 *
 * \brief Get array containing hyperplane color
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  hyperplane_color      Hyperplane color
 * \param [in]  ownership             Ownership for \p hyperplane_color
 */
void
PDM_multipart_part_hyperplane_color_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **hyperplane_color,
      PDM_ownership_t   ownership
);

/**
 *
 * \brief Get array containing thread color - Only if specific reordering (in ParaDiGMA plugins)
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  thread_color          Thread color
 * \param [in]  ownership             Ownership for \p thread_color
 */
void
PDM_multipart_part_thread_color_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **thread_color,
      PDM_ownership_t   ownership
);


/**
 *
 * \brief Get array containing vtx_ghost_information, useful to have a priority on vertex between multiple partitions
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t instance
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  vtx_ghost_information Integer that gives the current priority of vertices on current partitions
 * \param [in]  ownership             Ownership for \p vtx_ghost_information
 */
void
PDM_multipart_part_ghost_infomation_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **vtx_ghost_information,
      PDM_ownership_t   ownership
);


/**
 *
 * \brief Return times for a given domain
 * (NOT IMPLEMENTED)
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain       Id of current domain
 * \param [out]  elapsed        Elapsed time
 * \param [out]  cpu            CPU time
 * \param [out]  cpu_user       User CPU time
 * \param [out]  cpu_sys        System CPU time
 *
 */
void
PDM_multipart_time_get
(
 PDM_multipart_t *multipart,
 const int        i_domain,
 double         **elapsed,
 double         **cpu,
 double         **cpu_user,
 double         **cpu_sys
);


/**
 *
 * \brief Free the structure
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t instance
 */
void
PDM_multipart_free
(
 PDM_multipart_t *multipart
);

/**
 *
 * \brief Get the vertex coordinates on current i_domain, i_part partition and return number of vertices
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain       Id of current domain
 * \param [in]   i_part         Id of part
 * \param [out]  vtx_coord      Vertex coordinate (size = 3 * n_vtx)
 * \param [in]   ownership      Ownership for \p vtx_coord
 *
 * \return Number of vertices
 *
 */
int
PDM_multipart_part_vtx_coord_get
(
PDM_multipart_t                *multipart,
const int                       i_domain,
const int                       i_part,
      double                  **vtx_coord,
      PDM_ownership_t           ownership
);


/**
 *
 * \brief Get the group description for a given entity
 *
 * \param [in]   multipart              Pointer to \ref PDM_multipart_t instance
 * \param [in]   i_domain               Domain identifier
 * \param [in]   i_part                 Partition identifier
 * \param [in]   entity_type            Type of mesh entity
 * \param [out]  n_group                Number of groups
 * \param [out]  group_entity_idx       Index for group->entity connectivity (size = \p n_group)
 * \param [out]  group_entity           Group->entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [out]  group_entity_ln_to_gn  Group->entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 * \param [in]   ownership              Ownership
 *
 */
void PDM_multipart_group_get
(
 PDM_multipart_t      *multipart,
 const int             i_domain,
 const int             i_part,
 PDM_mesh_entities_t   entity_type,
 int                  *n_group_entity,
 int                 **group_entity_idx,
 int                 **group_entity,
 PDM_g_num_t         **group_entity_ln_to_gn,
 PDM_ownership_t       ownership
);


/**
 *
 * \brief Return statistics
 *
 * \param [in]   multipart                      Pointer to \ref PDM_multipart instance
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells number
 * \param [out]  cells_max                      maximum of cells number
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */
void
PDM_multipart_stat_get
(
 PDM_multipart_t  *multipart,
 int               i_domain,
 int              *cells_average,
 int              *cells_median,
 double           *cells_std_deviation,
 int              *cells_min,
 int              *cells_max,
 int              *bound_part_faces_average,
 int              *bound_part_faces_median,
 double           *bound_part_faces_std_deviation,
 int              *bound_part_faces_min,
 int              *bound_part_faces_max,
 int              *bound_part_faces_sum
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_H__ */
