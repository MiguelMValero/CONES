/*
 * \file
 */

#ifndef __PDM_PARTITIONING_ALGORITHM_H__
#define __PDM_PARTITIONING_ALGORITHM_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_part.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *  \brief Gather the entities split by the partitioner
 *   (usually cells) to their attributed partition, using the array mapping
 *   entities id to their assigned partition number.
 *   Each partition is hold by a (unique) process following the input partition
 *   distribution. The connection between partition members and original entities
 *   is made through the local to global numbering computed by the function.
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   part_distribution     Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution   Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dentity_to_part       Id of assigned partition for each entity (size=dn_entity)
 * \param [in]   dentity_gnum          If not null specified the current gnum in the block
 * \param [in]   dentity_init_location If not null specified the current gnum in the block
 * \param [out]  pn_entities           Number of entities in each partition (size = n_part)
 * \param [out]  pentity_ln_to_gn      Array of local to global entity id for each partition (size = n_part)
 *
 * \return       n_part              Number of partitions managed by this process
*/
int
PDM_part_assemble_partitions
(
 const PDM_MPI_Comm    comm,
       PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
 const PDM_g_num_t    *dentity_gnum,
 const int            *dentity_init_location,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn,
       int          ***pentity_init_location
);

/**
 *  \brief Construct the face->cell connectivity from the cell->face connectivity
 *   Since we assume that a face->cell is an array of two columns, this function
 *   can not be used to reverse other connectivities such as face->vtx.
 *
 *   This function respect the orientation : negative face in cell->face connectivity
 *   indicates that cell is right in the face->cell connectivity.
 *
 * \param [in]   n_part             Number of partitions
 * \param [in]   np_cell            Number of cells in each partition (size=n_part)
 * \param [in]   np_face            Number of faces in each partition (size=n_part)
 * \param [in]   pcell_face_idx     2d array of cell to face connectivity indexes
 *                                  (size = n_part*np_cell[i_part])
 * \param [in]   pcell_face         2d array of cell to face connectivity
 *                                  (size = n_part*pcell_face_idx[i_part][np_cell+1])
 * \param [out]  pface_face         2d array of face to cell connectivity
 *                                  (size = n_part*2*np_face[i_part])
 */
void
PDM_part_reverse_pcellface
(
  const int         n_part,
  const int        *n_cell,
  const int        *n_face,
  const int       **pcell_face_idx,
  const int       **pcell_face,
        int      ***pface_cell
);

/**
 *  \brief Reorient the boundary faces such that they have a outward normal for the boundary cell.
 *   This functions only uses topological information (face->cell connectivity) to determine
 *   if the face must be reoriented. Thus, input face->cell and cell->face connectivities must
 *   contain correct orientation information.
 *
 *
 * \param [in]    n_part             Number of partitions
 * \param [in]    pn_face            Number of faces in each partition (size=n_part)
 * \param [inout] pface_face         On each part, array of face to cell connectivity
 *                                   (size = n_part*2*np_face[i_part])
 * \param [in]    pcell_face_idx     On each part, array of cell to face connectivity indexes
 *                                   (size = n_part*np_cell[i_part])
 * \param [inout] pcell_face         On each part, array of cell to face connectivity
 *                                   (size = n_part*pcell_face_idx[i_part][np_cell+1])
 * \param [in]    pface_vtx_idx      On each part, array of face to vertex connectivity indexes
 *                                   (size = n_part*np_vtx[i_part])
 * \param [inout] pface_vtx          On each part, array of face to vertex connectivity
 *                                  (size = n_part*pface_vtx_idx[i_part][np_vtx+1])
 */
void
PDM_part_reorient_bound_faces
(
  const int         n_part,
  const int        *np_face,
        int       **pface_cell,
  const int       **pcell_face_idx,
        int       **pcell_face,
  const int       **pface_vtx_idx,
        int       **pface_vtx,
        int       **pface_edge_idx,
        int       **pface_edge
);

/**
 *  \brief Recover partitioned entity groups (cell, face, vertex) from distributed
 *   entity groups. Return the list of local element id belonging to each group,
 *   and the position of those entities in the corresponding original (distributed) group.
 *
 *   This function is especially used to retrieve boundary conditions which are defined as
 *   face groups.
 *
 *   n_group is a global data that must be know by each process, even if they
 *   dont hold any group element.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   n_group             Number of groups defined for this entity
 * \param [in]   dgroup_idx          Number of distributed elements in each group (size=n_group+1)
 * \param [in]   dgroup              Global id of entities belonging to the groups (size=dgroup_idx[n_group])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entities         Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pgroup_idx          For each part, number of partitioned elements in each group
 *                                   (size = n_part, each component size = n_group+1)
 * \param [out]  pgroup              For each part, local id of entities belonging to the groups
 *                                   (size = n_part, each component size = pgroup_idx[n_group])
 * \param [out]  pgroup_ln_to_gn     For each part, position of entity in the original groups
 *                                   (size = n_part, each component size = pgroup_idx[n_group])
 */
void
PDM_part_distgroup_to_partgroup
(
 const PDM_MPI_Comm      comm,
 const PDM_g_num_t      *entity_distribution,
 const int               n_group,
 const int              *dgroup_idx,
 const PDM_g_num_t      *dgroup,
 const int               n_part,
 const int              *pn_entity,
 const PDM_g_num_t     **pentity_ln_to_gn,
       int            ***pgroup_idx,
       int            ***pgroup,
       PDM_g_num_t    ***pgroup_ln_to_gn
);


/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   The orientation data (i.e. negative index) present in the distributed connectivity,
 *   if any, are preserved, *meaning that boundary faces can be badly oriented on partitions*.
 *   See PDM_part_reorient_bound_faces function to correct this orientation.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition (size=n_part)
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array
 *                                   (size = n_part, each component size = pn_child_entity[i_part])
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes
 *                                   (size = n_part, each component size = pn_entity[i_part])
 * \param [out]  pconnectivity       For each part, partitioned connectivity (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
 */
void
PDM_part_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
);


/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   The orientation data (i.e. negative index) present in the distributed connectivity,
 *   if any, are preserved, *meaning that boundary faces can be badly oriented on partitions*.
 *   See PDM_part_reorient_bound_faces function to correct this orientation.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   pn_entity           Number of entities in each partition
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition (size = pn_entity)
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array (size = pn_entity)
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes (size = pn_entity+1)
 * \param [out]  pconnectivity       For each part, partitioned connectivity (size = pconnectivity_idx[pn_entity])
 */
void
PDM_part_dconnectivity_to_pconnectivity_sort_single_part
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             pn_entity,
 const PDM_g_num_t    *pentity_ln_to_gn,
       int            *pn_child_entity,
       PDM_g_num_t   **pchild_ln_to_gn,
       int           **pconnectivity_idx,
       int           **pconnectivity
);

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   -- multi-section version
 */
void
PDM_part_multi_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const int             n_part,
 const int             n_section,
 const int            *section_idx,
       PDM_g_num_t   **entity_distribution,
       int            *dconnectivity_idx,
       PDM_g_num_t    *dconnectivity,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int         ****pconnectivity_idx,
       int         ****pconnectivity
);

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   The orientation data (i.e. negative index) present in the distributed connectivity,
 *   if any, are preserved, *meaning that boundary faces can be badly oriented on partitions*.
 *   See PDM_part_reorient_bound_faces function to correct this orientation.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition (size=n_part)
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array
 *                                   (size = n_part, each component size = pn_child_entity[i_part])
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes
 *                                   (size = n_part, each component size = pn_entity[i_part])
 * \param [out]  pconnectivity       For each part, partitioned connectivity (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
 */
void
PDM_part_dconnectivity_to_pconnectivity_hash
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
);

/**
 *  \brief Generates the communication information at the partition interfaces for the
 *   given entity. The communication data associates to
 *   each partitioned entity belonging to an (internal) interface the 4-tuple
 *   (local id, opposite proc number, opposite part number on opposite proc, local id in the
 *   opposite partition).
 *   This list is sorted by opposite proc id, then by part id, and finally with respect
 *   to the entity global_id. Also return the stride indices pproc_bound_idx and ppart_bound_idx
 *   to acces the communication information for a given opposite proc id or global part id.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   part_distribution   Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entity[i_part])
 * \param [in]   pentity_hint        Can be used to indicate whether (1) or not (0) an entity is potentially
 *                                   shared with an other partition in order to minimize exchanged data
 *                                   (size=n_part, each component size = pn_entity[i_part]) or NULL
 * \param [out]  pproc_bound_idx     For each part, indexes of communication information related to the
 *                                   other procs (size=n_part, each component size=n_rank+1)
 * \param [out]  ppart_bound_idx     For each part, indexes of communication information related to the
 *                                   other (global id) parts (size=n_part, each component size=n_part_tot+1)
 * \param [out]  pentity_bound       For each part, communication information (see above) (size=n_part)
 */
void
PDM_part_generate_entity_graph_comm
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t   *part_distribution,
 const PDM_g_num_t   *entity_distribution,
 const int            n_part,
 const int           *pn_entity,
 const PDM_g_num_t  **pentity_ln_to_gn,
 const int          **pentity_hint,
       int         ***pproc_bound_idx,
       int         ***ppart_bound_idx,
       int         ***pentity_bound,
       int         ***pentity_priority
);

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dcoordinates_to_pcoordinates
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  const PDM_g_num_t    *vertex_distribution,
  const double         *dvtx_coord,
  const int            *pn_vtx,
  const PDM_g_num_t   **pvtx_ln_to_gn,
        double       ***pvtx_coord
);

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dfield_to_pfield
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  size_t                s_data,
  const PDM_g_num_t    *field_distribution,
  const unsigned char  *dfield,
  const int            *pn_field,
  const PDM_g_num_t   **pfield_ln_to_gn,
        unsigned char ***pfield
);


/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dfield_to_pfield2
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  const PDM_g_num_t     *field_distribution,
  const int             *dfield_stri,
  const unsigned char   *dfield,
  const int             *pn_field,
  const PDM_g_num_t    **pfield_ln_to_gn,
  int                 ***pfield_stride,
        unsigned char ***pfield
);

/**
 *  \brief Extend an existing ln_to_gn from a connectivity
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   part_distribution   Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dentity_to_part     Id of assigned partition for each entity (size=dn_entity)
 * \param [out]  pn_entities         Number of entities in each partition (size = n_part)
 * \param [out]  pentity_ln_to_gn    Array of local to global entity id for each partition (size = n_part)
 *
 * \return       n_part              Number of partitions managed by this process
 */
void
PDM_extend_mesh
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
 const int             n_part,
 const PDM_g_num_t    *dual_graph_idx,
 const PDM_g_num_t    *dual_graph,
 const int            *pn_entity,
       PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_entity_extented,
       PDM_g_num_t  ***pentity_ln_to_gn_extended
);


/**
 *  \brief Deduce group for each partition from the distributed one (for cells, faces, edges and vtx)
 *
 * \param [in]  comm                PDM_MPI communicator
 * \param [in]  n_part              Number of partitions
 * \param [in]  entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]  dentity_group_idx   Connectivity index between entity and group (size = dn_entity)
 * \param [in]  dentity_group       For each entity the associated group (size = dentity_group_idx[dn_entity])
 * \param [in]  pn_entity           Number of entities in each partition (size = n_part)
 * \param [in]  pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out] pentity_group_idx   Connectivity index between entity and group (size = n_part)
 * \param [out] pentity_group       For each entity the associated group  (size = n_part)
 *
 */
void
PDM_part_dentity_group_to_pentity_group
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  const PDM_g_num_t     *entity_distribution,
  const int             *dentity_group_idx,
  const int             *dentity_group,
  const int             *pn_entity,
  const PDM_g_num_t    **pentity_ln_to_gn,
  int                 ***pentity_group_idx,
  int                 ***pentity_group
);


/**
 *  \brief Compute the explicit distributed connectivity from an implicit one, with a prescribed stride. Use to convert for exemple face_cell or edge_vtx implicit connectivity
 *
 * \param [in]   dn_entity1             Number of entity1
 * \param [in]   stride                 Implicit stride of dentity1_dentity2 connectivity
 * \param [in]   dentity1_dentity2      Implicit connectivity between entity and group (size = stride * dn_entity1)
 * \param [out]  dentity1_dentity2_idx  Connectivity index between entity1 and entity2 (size = dn_entity1)
 * \param [out]  dentity1_dentity2_new  Connectivity index between entity1 and entity2 (size = dentity1_dentity2_idx[dn_entity1])
 *
 */
void
PDM_setup_connectivity_idx
(
  int           dn_entity1,
  int           stride,
  PDM_g_num_t  *dentity1_dentity2,
  int         **dentity1_dentity2_idx,
  PDM_g_num_t **dentity1_dentity2_new
);


/**
 *  \brief Compute the edges for all partitions in an independent of parallelism way. Useful when user only give face_vtx but edges is mandatory for algorithm (ex : Iso-surfaces)
 *
 * \warning This function uses collective MPI communications and must be called simultaneously by all processes in \p comm.
 *
 * \param [in]  comm                PDM_MPI communicator
 * \param [in]  n_part              Number of partitions
 * \param [in]  pn_face             Number of faces for each partition (size = \p n_part)
 * \param [in]  pn_vtx              Number of vertices for each partition (size = \p n_part)
 * \param [in]  pface_vtx_idx       For each part, connectivity index between faces and vertices
 *                                 (size = \p n_part, each component size = \p pn_face[i_part]+1)
 * \param [in]  pface_vtx           For each part, connectivity between faces and vertices
 *                                 (size = \p n_part, each component size = \p pface_vtx_idx[i_part][\p pn_face[i_part]])
 * \param [in]  pface_ln_to_gn      For each part, face global ids (size = \p n_part, each component size = \p pn_face[i_part])
 * \param [in]  pvtx_ln_to_gn       For each part, vertex global ids (size = \p n_part, each component size = \p pn_vtx[i_part])
 * \param [out] pface_edge_idx      For each part, connectivity index between faces and edges
 *                                  (size = \p n_part, each component size = \p pn_face[i_part]+1)
 * \param [in]  pface_edge          For each part, connectivity between faces and edges
 *                                 (size = \p n_part, each component size = \p pface_edge_idx[i_part][\p pn_face[i_part]])
 * \param [out] pn_edge             Number of edges for each partition (size = \p n_part)
 * \param [out] pedge_vtx           For each part, connectivity between edges and vertices
 *                                 (size = n_part, each component size = 2 * \p pn_edge[i_part])
 * \param [out] pedge_ln_to_gn      For each part, edge global ids (size = \p n_part, each component size = \p pn_edge[i_part])
 *
 */
void
PDM_compute_face_edge_from_face_vtx
(
  PDM_MPI_Comm    comm,
  int             n_part,
  int            *pn_face,
  int            *pn_vtx,
  int           **pface_vtx_idx,
  int           **pface_vtx,
  PDM_g_num_t   **pface_ln_to_gn,
  PDM_g_num_t   **pvtx_ln_to_gn,
  int          ***pface_edge_idx,
  int          ***pface_edge,
  int           **pn_edge,
  int          ***pedge_vtx,
  PDM_g_num_t  ***pedge_ln_to_gn
);


/**
 *  \brief Deduce connectivity in a new partition from another one. See \ref PDM_extract_part_t for exemple of use
 *
 * \param [in]  comm                               PDM_MPI communicator
 * \param [in]  n_part1                            Number of partitions in first partitioning
 * \param [in]  n_part1_entity1                    For each part, number of entity1
 * \param [in]  part1_entity1_entity2_idx          For each part, for partition 1, connectivity index between entity1 and entity2
 *                                                  (size = n_part1, each component size = n_part1_entity1[i_part]+1)
 * \param [in]  part1_entity1_entity2              For each part, for partition1, connectivity between entity1 and entity2
 *                                                  (size = n_part1, each component size = part1_entity1_entity2_idx[n_part1_entity1[i_part]])
 * \param [in]  part1_entity1_ln_to_gn             For each part, for partition 1, global id of entity1
 * \param [in]  part1_entity2_ln_to_gn             For each part, for partition 1, global id of entity2
 * \param [in]  n_part2                            Number of partitions in second partitioning
 * \param [in]  n_part2_entity1                    For each part, number of entity1
 * \param [in]  part2_entity1_ln_to_gn             For each part, for partition 2, global id of entity1
 * \param [in]  part2_entity1_to_part1_entity1_idx For each part, for partition 2, connectivity index between part2_entity1 and part1_entity1
 * \param [in]  part2_entity1_to_part1_entity1     For each part, for partition 2, connectivity between part2_entity1 and part1_entity1 (global id)
 * \param [out] n_part2_entity2                    For each part, number of entity2
 * \param [out] part2_entity1_entity2_idx          For each part, for partition 2, connectivity index between entity1 and entity2
 *                                                  (size = n_part1, each component size = n_part2_entity2[i_part]+1)
 * \param [out] part2_entity1_entity2              For each part, for partition 2, connectivity between entity1 and entity2
 *                                                  (size = n_part1, each component size = part2_entity1_entity2_idx[n_part2_entity2[i_part]])
 * \param [out] part2_entity2_child_ln_to_gn       For each part, for partition 2, global id of child entity2
 * \param [out] part2_entity2_ln_to_gn             For each part, for partition 2, global id of entity2
 *
 */
void
PDM_pconnectivity_to_pconnectivity
(
  const PDM_MPI_Comm    comm,
  const int             n_part1,
  const int            *n_part1_entity1,
  const int           **part1_entity1_entity2_idx,
  const int           **part1_entity1_entity2,
  const PDM_g_num_t   **part1_entity1_ln_to_gn,
  const PDM_g_num_t   **part1_entity2_ln_to_gn,
  const int             n_part2,
  const int            *n_part2_entity1,
  const PDM_g_num_t   **part2_entity1_ln_to_gn,
  const int           **part2_entity1_to_part1_entity1_idx,
  const PDM_g_num_t   **part2_entity1_to_part1_entity1,
        int           **n_part2_entity2,
        int          ***part2_entity1_entity2_idx,
        int          ***part2_entity1_entity2,
        PDM_g_num_t  ***part2_entity2_child_ln_to_gn,
        PDM_g_num_t  ***part2_entity2_ln_to_gn
);


/**
 *  \brief Deduce connectivity in a new partition from another one. See \ref PDM_extract_part_t for example of use
 *
 * \param [in]  comm                               PDM_MPI communicator
 * \param [in]  n_part1                            Number of partitions in first partitioning
 * \param [in]  n_part1_entity1                    For each part, number of entity1
 * \param [in]  part1_entity1_entity2_idx          For each part, for partition 1, connectivity index between entity1 and entity2
 *                                                  (size = n_part1, each component size = n_part1_entity1[i_part]+1)
 * \param [in]  part1_entity1_entity2              For each part, for partition1, connectivity between entity1 and entity2
 *                                                  (size = n_part1, each component size = part1_entity1_entity2_idx[n_part1_entity1[i_part]])
 * \param [in]  part1_entity1_ln_to_gn             For each part, for partition 1, global id of entity1
 * \param [in]  part1_entity2_ln_to_gn             For each part, for partition 1, global id of entity2
 * \param [in]  n_part2                            Number of partitions in second partitioning
 * \param [in]  n_part2_entity1                    For each part, number of entity1
 * \param [in]  part2_entity1_ln_to_gn             For each part, for partition 2, global id of entity1
 * \param [in]  part2_entity1_to_part1_entity1_idx For each part, for partition 2, connectivity index between part2_entity1 and part1_entity1
 * \param [in]  part2_entity1_to_part1_entity1     For each part, for partition 2, connectivity between part2_entity1 and part1_entity1 (global id)
 * \param [out] n_part2_entity2                    For each part, number of entity2
 * \param [out] part2_entity1_entity2_idx          For each part, for partition 2, connectivity index between entity1 and entity2
 *                                                  (size = n_part2, each component size = n_part2_entity2[i_part]+1)
 * \param [out] part2_entity1_entity2              For each part, for partition 2, connectivity between entity1 and entity2
 *                                                  (size = n_part2, each component size = part2_entity1_entity2_idx[n_part2_entity2[i_part]])
 * \param [out] part2_entity2_ln_to_gn             For each part, for partition 2, global id of entity2
 * \param [out] part2_entity2_child_ln_to_gn       For each part, for partition 2, global id of child entity2
 * \param [out] ptp                                Part to part exchange protocol (see \ref PDM_part_to_part_t ). Useful to exchange additional data between part1 and part2
 *
 */
void
PDM_pconnectivity_to_pconnectivity_keep
(
  const PDM_MPI_Comm          comm,
  const int                   n_part1,
  const int                  *n_part1_entity1,
  const int                 **part1_entity1_entity2_idx,
  const int                 **part1_entity1_entity2,
  const PDM_g_num_t         **part1_entity1_ln_to_gn,
  const PDM_g_num_t         **part1_entity2_ln_to_gn,
  const int                   n_part2,
  const int                  *n_part2_entity1,
  const PDM_g_num_t         **part2_entity1_ln_to_gn,
  const int                 **part2_entity1_to_part1_entity1_idx,
  const PDM_g_num_t         **part2_entity1_to_part1_entity1,
        int                 **n_part2_entity2,
        int                ***part2_entity1_entity2_idx,
        int                ***part2_entity1_entity2,
        PDM_g_num_t        ***part2_entity2_ln_to_gn,
        PDM_g_num_t        ***part2_entity2_child_ln_to_gn,
        PDM_part_to_part_t  **ptp
);


/**
 *  \brief Deduce connectivity in a new partition from another one. See \ref PDM_extract_part_t for example of use
 *
 * \param [in]  comm                                    PDM_MPI communicator
 * \param [in]  n_part1                                 Number of partitions in first partitioning
 * \param [in]  n_part1_entity1                         For each part, number of entity1
 * \param [in]  part1_entity1_entity2_idx               For each part, for partition 1, connectivity index between entity1 and entity2
 *                                                       (size = n_part1, each component size = n_part1_entity1[i_part]+1)
 * \param [in]  part1_entity1_entity2                   For each part, for partition1, connectivity between entity1 and entity2
 *                                                       (size = n_part1, each component size = part1_entity1_entity2_idx[n_part1_entity1[i_part]])
 * \param [in]  part1_entity1_ln_to_gn                  For each part, for partition 1, global id of entity1
 * \param [in]  part1_entity2_ln_to_gn                  For each part, for partition 1, global id of entity2
 * \param [in]  n_part2                                 Number of partitions in second partitioning
 * \param [in]  n_part2_entity1                         For each part, number of entity1
 * \param [in]  part2_entity1_ln_to_gn                  For each part, for partition 2, global id of entity1
 * \param [in]  part2_entity1_to_part1_entity1_idx      For each part, for partition 2, connectivity index between part2_entity1 and part1_entity1
 * \param [in]  part2_entity1_to_part1_entity1_triplet  For each part, for partition 2, connectivity between part2_entity1 and part1_entity1 by triplet (i_proc, i_part, i_entity)
 * \param [out] n_part2_entity2                         For each part, number of entity2
 * \param [out] part2_entity1_entity2_idx               For each part, for partition 2, connectivity index between entity1 and entity2
 *                                                       (size = n_part2, each component size = n_part2_entity2[i_part]+1)
 * \param [out] part2_entity1_entity2                   For each part, for partition 2, connectivity between entity1 and entity2
 *                                                       (size = n_part2, each component size = part2_entity1_entity2_idx[n_part2_entity2[i_part]])
 * \param [out] part2_entity2_ln_to_gn                  For each part, for partition 2, global id of entity2
 * \param [out] part2_entity2_child_ln_to_gn            For each part, for partition 2, global id of child entity2
 * \param [out] ptp                                     Part to part exchange protocol (see \ref PDM_part_to_part_t ). Useful to exchange additional data between part1 and part2
 *
 */
void
PDM_pconnectivity_to_pconnectivity_from_location_keep
(
  const PDM_MPI_Comm          comm,
  const int                   n_part1,
  const int                  *n_part1_entity1,
  const int                 **part1_entity1_entity2_idx,
  const int                 **part1_entity1_entity2,
  const PDM_g_num_t         **part1_entity2_ln_to_gn,
  const int                   n_part2,
  const int                  *n_part2_entity1,
  const PDM_g_num_t         **part2_entity1_ln_to_gn,
  const int                 **part2_entity1_to_part1_entity1_idx,
  const int                 **part2_entity1_to_part1_entity1_triplet,
        int                 **n_part2_entity2,
        int                ***part2_entity1_entity2_idx,
        int                ***part2_entity1_entity2,
        PDM_g_num_t        ***part2_entity2_ln_to_gn,
        int                ***part2_entity2_to_part1_entity2,
        PDM_part_to_part_t  **ptp_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_ALGORITHM_H__ */
