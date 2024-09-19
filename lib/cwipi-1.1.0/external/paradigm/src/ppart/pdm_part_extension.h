/*
 * \file
 */

#ifndef __PDM_PART_EXTENSION_H__
#define __PDM_PART_EXTENSION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"
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

typedef struct _pdm_part_extension_t PDM_part_extension_t;


typedef enum {

  PDM_EXTEND_FROM_FACE = 0,
  PDM_EXTEND_FROM_EDGE = 1,
  PDM_EXTEND_FROM_VTX  = 2

} PDM_extend_type_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a part extension structure
 *
 * \param [in] n_domain      Number of domains
 * \param [in] n_part        Number of partitions per domain
 * \param [in] extend_type   Extension from which entity ?
 * \param [in] depth         Extension depth
 * \param [in] comm          MPI communicator
 * \param [in] owner         Data ownership
 *
 *  \return   \p PDM_part_extension_t structure instance
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int                n_domain,
 const int               *n_part,
       PDM_extend_type_t  extend_type,
       int                depth,
 const PDM_MPI_Comm       comm,
 const PDM_ownership_t    owner
);

/**
 *
 * \brief Compute extended partitions
 *
 * \param [in]   part_ext          \p PDM_part_extension_t structure instance
 *
 */

void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext
);

/**
 *
 * \brief Set data to perform the partitioned mesh extension
 *
 * \warning Deprecated: use the individual setters instead
 *
 * \param [in]   part_ext                  \p PDM_part_extension_t structure instance
 * \param [in]   i_domain                  Domain identifier
 * \param [in]   i_part                    Partition identifier
 * \param [in]   n_cell                    Number of cells
 * \param [in]   n_face                    Number of faces
 * \param [in]   n_face_part_bound         Number of partition boundary faces
 * \param [in]   n_face_group              Number of face groups
 * \param [in]   n_edge                    Number of edges
 * \param [in]   n_vtx                     Number of vertices
 * \param [in]   cell_face_idx             Cell-face connectivity index (size = \ref n_cell + 1)
 * \param [in]   cell_face                 Cell-face connectivity (size = \ref cell_face_idx(\ref n_cell + 1))
 * \param [in]   face_cell                 Face-cell connectivity (size = 2 * \ref n_face)
 * \param [in]   face_edge_idx             Face-edge connectivity index (size = \ref n_face + 1)
 * \param [in]   face_edge                 Face-edge connectivity (size = \ref face_edge_idx[\ref n_face])
 * \param [in]   face_vtx_idx              Face-vertex connectivity index (size = \ref n_face + 1)
 * \param [in]   face_vtx                  Face-vertex connectivity (size = \ref face_vtx_idx[\ref n_face])
 * \param [in]   edge_vtx                  Edge-vertex connectivity (size = 2 * \ref n_edge)
 * \param [in]   face_bound_idx            Face->group connectivity index (size = \ref n_face_group + 1)
 * \param [in]   face_bound                Face->group connectivity (size = \ref face_edge_idx[\ref n_face_group])
 * \param [in]   face_join_idx             Faces connecting domains connectivity index
 * \param [in]   face_join                 Faces connecting domains connectivity
 * \param [in]   face_part_bound_proc_idx  Partitioning boundary faces index from process (size = n_proc + 1)
 * \param [in]   face_part_bound_part_idx  Partitioning boundary faces index from partition (size = n_total_part + 1)
 * \param [in]   face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
 * \param [in]   vtx_part_bound_proc_idx   Partitioning boundary vertices index from process (size = n_proc + 1)
 * \param [in]   vtx_part_bound_part_idx   Partitioning boundary vertices index from partition (size = n_total_part + 1)
 * \param [in]   vtx_part_bound            Partitioning boundary vertices (size = 4 * n_vertex_part_bound)
 * \param [in]   cell_ln_to_gn             Cell global ids (size = \ref n_cell)
 * \param [in]   face_ln_to_gn             Face global ids (size = \ref n_face)
 * \param [in]   edge_ln_to_gn             Edge global ids (size = \ref n_edge)
 * \param [in]   vtx_ln_to_gn              Vertex global ids (size = \ref n_vtx)
 * \param [in]   face_group_ln_to_gn       Global ids of faces with groups (size = \ref n_face_group)
 * \param [in]   vtx_coord                 Vertex coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                   n_cell,
  int                   n_face,
  int                   n_face_part_bound,
  int                   n_face_group,
  int                   n_edge,
  int                   n_vtx,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  int                  *vtx_part_bound_proc_idx,
  int                  *vtx_part_bound_part_idx,
  int                  *vtx_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn,
  double               *vtx_coord
);

/**
 *
 * \brief Use shared domain interface
 *
 * \param [in]   part_ext                    \p PDM_part_extension_t structure instance
 * \param [in]   PDM_part_domain_interface_t \p PDM_part_domain_interface_t structure instance
 *
 */

void
PDM_part_extension_part_domain_interface_shared_set
(
  PDM_part_extension_t        *part_ext,
  PDM_part_domain_interface_t *pdi
);

/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          \p PDM_part_extension_t structure instance
 *
 */

void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
);


/**
 *
 * \brief Get extended connectivity
 *
 * \param [in]  part_ext            \p PDM_part_extension_t structure instance
 * \param [in]  i_domain            Domain identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  connectivity_type   Connectivity type
 * \param [out] connect_idx         Connectivity index
 * \param [out] connect             Connectivity
 *
 * \return Number of leading entities
 *
 */

int
PDM_part_extension_connectivity_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect_idx,
 int                     **connect
);


/**
 *
 * \brief Get global ids of extended entities
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Entity type
 * \param [out] ln_to_gn     Global ids
 *
 * \return  Number of entities
 *
 */

int
PDM_part_extension_ln_to_gn_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 PDM_g_num_t             **ln_to_gn
);

/**
 *
 * \brief Get interface
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] interface_no Interfaces
 *
 * \return  Number of interfaces
 *
 */

int
PDM_part_extension_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **interface_no
);

/**
 *
 * \brief Get groups for extended entities with given type
 *
 * \param [in]  part_ext                \p PDM_part_extension_t structure instance
 * \param [in]  i_domain                Domain identifier
 * \param [in]  i_part                  Partition identifier
 * \param [in]  mesh_entity             Type of mesh entity
 * \param [out] group_entity_idx        Group->entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [out] group_entity            Index for group->entity connectivity (size = \p n_group)
 * \param [out] group_entity_ln_to_gn   Group->entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 *
 * \return  Number of groups
 *
 */

int
PDM_part_extension_group_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **group_entity_idx,
 int                     **group_entity,
 PDM_g_num_t             **group_ln_to_gn
);


/**
 *
 * \brief Get coordinates of extended vertices
 *
 * \param [in]  part_ext     PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [out] vtx_coord    Vertex coordinates (size = \ref n_vtx * 3)
 *
 * \return  Number of extended vertices
 *
 */

int
PDM_part_extension_vtx_coord_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                  **vtx_coord
);

/**
 *
 * \brief Get composed interface
 *
 * \param [in]  part_ext                 PDM_part_extension_t structure instance
 * \param [out] composed_interface_idx   ??
 * \param [out] composed_interface       ??
 * \param [out] composed_ln_to_gn_sorted ??
 *
 * \return  something ??
 *
 */

int
PDM_part_extension_composed_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                     **composed_interface_idx,
 int                     **composed_interface,
 PDM_g_num_t             **composed_ln_to_gn_sorted
);


/**
 *
 * \brief Create part_to_part from interior and extended elements
 *
 * \param [out]  ptp                             Part to part structure
 * \param [in]   n_part                          Number of partitions
 * \param [in]   n_int_cell                      Number of interior elements
 * \param [in]   int_cell_ln_to_gn               gnum of interior elements
 * \param [in]   n_ext_cell                      Number of extended elements
 * \param [in]   ext_cell_ln_to_gn               gnum of extended elements
 * \param [out]  n_selected_cell_to_send         Number of elements selected for send
 * \param [out]  selected_cell_to_send           Local numbering of elements selected for send
 *
 */

void
PDM_part_to_part_create_from_extension
(
       PDM_part_to_part_t **ptp,
 const int                  n_part,
       int                 *n_int_cell,
 const PDM_g_num_t        **int_cell_ln_to_gn,
       int                 *n_ghost_cell,
 const PDM_g_num_t        **ghost_cell_ln_to_gn,
       int                **n_selected_cell_to_send,
       int               ***selected_cell_to_send,
 const PDM_MPI_Comm         comm
);





/**
 *
 * \brief Set connectivity
 *
 * \param [in]  part_ext           \p PDM_part_extension_t structure instance
 * \param [in]  i_domain           Domain identifier
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  connect            Connectivity
 *
 */

void
PDM_part_extension_connectivity_set
(
 PDM_part_extension_t    *part_ext,
 int                      i_domain,
 int                      i_part,
 PDM_connectivity_type_t  connectivity_type,
 int                     *connect_idx,
 int                     *connect
 );

/**
 *
 * \brief Set global ids
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [in]  n_entity     Local number of entities
 * \param [in]  ln_to_gn     Global ids (size = \p n_entity)
 *
 */

void
PDM_part_extension_ln_to_gn_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                       n_entity,
 PDM_g_num_t              *ln_to_gn
);

/**
 *
 * \brief Set vertex coordinates
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  vtx_coord    Vertex coordinates (size = 3 * *n_vtx*)
 *
 */

void
PDM_part_extension_vtx_coord_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                   *vtx_coord
);

/**
 *
 * \brief Set the connection graph between partitions for the requested entity type
 *
 * \param [in]  multipart             \p PDM_part_extension_t structure instance
 * \param [in]  i_domain              Domain identifier
 * \param [in]  i_part                Partition identifier
 * \param [in]  entity_type           Type of mesh entity
 * \param [in]  part_bound_proc_idx   Partitioning boundary entities index from process (size = *n_rank* + 1)
 * \param [in]  part_bound_part_idx   Partitioning boundary entities index from partition (size = *n_total_part* + 1)
 * \param [in]  part_bound            Partitioning boundary entities (size = 4 * \p part_bound_proc_idx[*n_rank*])
 */

void
PDM_part_extension_part_bound_graph_set
(
 PDM_part_extension_t *part_ext,
 int                   i_domain,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                  *part_bound_proc_idx,
 int                  *part_bound_part_idx,
 int                  *part_bound
);

/**
 *
 * \brief Set group description
 *
 * \param [in]  part_ext               \p PDM_part_extension_t structure instance
 * \param [in]  i_domain               Domain identifier
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 * \param [in]  group_entity_idx       Index for group->entity connectivity (size = \p n_group)
 * \param [in]  group_entity           Group->entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  group_entity_ln_to_gn  Group->entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 *
 */

void
PDM_part_extension_group_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       n_group,
 int                      *group_entity_idx,
 int                      *group_entity,
 PDM_g_num_t              *group_entity_ln_to_gn
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_H__ */
