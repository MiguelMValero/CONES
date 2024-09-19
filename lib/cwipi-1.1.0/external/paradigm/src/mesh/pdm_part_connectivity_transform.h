/*
 * \file
 */

#ifndef PDM_PART_CONNECTIVITY_TRANSFORM_H_
#define PDM_PART_CONNECTIVITY_TRANSFORM_H_

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 *
 * \brief Compress connectivity to be unique - Realloc need to be done extern to this method
 *
 * \param [in]     n_entity1           Number of entity in connectivity
 * \param [inout]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [inout]  entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_compress_connectivity
(
 int  n_entity1,
 int *entity1_entity2_idx,
 int *entity1_entity2
);

/**
 *
 * \brief Combine connectivity between entity1_entity2 and entity2_entity3 to have entity1_entity3
 *
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 * \param [in]  entity2_entity3_idx Connectivity index between entity2 and entity3 (size = n_entity2)
 * \param [in]  entity2_entity3     Connectivity between entity2 and entity3 (size = entity1_entity2_idx[n_entity2] )
 * \param [out] entity1_entity3_idx Connectivity index between entity1 and entity3 (size = n_entity1)
 * \param [out] entity1_entity3     Connectivity between entity1 and entity3 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_combine_connectivity
(
 int   n_entity1,
 int  *entity1_entity2_idx,
 int  *entity1_entity2,
 int  *entity2_entity3_idx,
 int  *entity2_entity3,
 int **entity1_entity3_idx,
 int **entity1_entity3
);


/**
 *
 * \brief Transpose connectivity entity1_entity2 to have entity2_entity1
 *
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  n_entity1           Number of entity2
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 * \param [out] entity2_entity1_idx Connectivity index between entity2 and entity1 (size = n_entity2)
 * \param [out] entity2_entity1     Connectivity between entity2 and entity1 (size = entity1_entity2_idx[n_entity2] )
 *
 */
void
PDM_connectivity_transpose
(
const int   n_entity1,
const int   n_entity2,
      int  *entity1_entity2_idx,
      int  *entity1_entity2,
      int **entity2_entity1_idx,
      int **entity2_entity1
);

/**
 *
 * \brief Transpose connectivity entity1_entity2 to have entity2_entity1 for all partition
 *
 * \param [in]  n_part              Number of partition in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  n_entity1           Number of entity2
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 * \param [out] entity2_entity1_idx Connectivity index between entity2 and entity1 (size = n_entity2)
 * \param [out] entity2_entity1     Connectivity between entity2 and entity1 (size = entity1_entity2_idx[n_entity2] )
 *
 */
void
PDM_part_connectivity_transpose
(
const int    n_part,
const int   *n_entity1,
const int   *n_entity2,
      int  **entity1_entity2_idx,
      int  **entity1_entity2,
      int ***entity2_entity1_idx,
      int ***entity2_entity1
);

/**
 *
 * \brief Combine connectivity between entity1_entity2 and entity2_entity3 to have entity1_entity3 for all partition
 *
 * \param [in]  n_part              Number of partition in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [in]  entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 * \param [in]  entity2_entity3_idx Connectivity index between entity2 and entity3 (size = n_entity2)
 * \param [in]  entity2_entity3     Connectivity between entity2 and entity3 (size = entity1_entity2_idx[n_entity2] )
 * \param [out] entity1_entity3_idx Connectivity index between entity1 and entity3 (size = n_entity1)
 * \param [out] entity1_entity3     Connectivity between entity1 and entity3 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_part_combine_connectivity
(
 const int    n_part,
 int         *n_entity1,
 int        **entity1_entity2_idx,
 int        **entity1_entity2,
 int        **entity2_entity3_idx,
 int        **entity2_entity3,
 int       ***entity1_entity3_idx,
 int       ***entity1_entity3
);

/**
 *
 * \brief Convert implicit pair connectivity, to a connectivity with index. Useful for convert face_cell or edge_vtx.
 *
 * \param [in]  n_part              Number of partition in current process
 * \param [in]  n_entity1           Number of entity1
 * \param [in]  entity1_entity2_in  Implicit connectivity (face_cell for exemple with right cell is boundary face[2*i+1] == 0)
 * \param [out] entity1_entity2_idx Connectivity index between entity1 and entity2 (size = n_entity1)
 * \param [out] entity1_entity2     Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_part_connectivity_to_connectity_idx
(
const int    n_part,
const int   *n_entity1,
      int  **entity1_entity2_in,
      int ***entity1_entity2_idx,
      int ***entity1_entity2
);

/**
 *
 * \brief Generate face_vtx according to the sign connectivity face_edge and edge_vtx
 *
 * \param [in]  n_face        Number of partition in current process
 * \param [in]  face_edge_idx Connectivity index between face and edge (size = n_face+1)
 * \param [in]  face_edge     Connectivity between face and edge (signed) (size = face_edge_idx[n_face])
 * \param [in]  edge_vtx      Connectivity between edge and vtx - Implicit (size = 2 * n_edge)
 * \param [out] face_vtx      Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_compute_face_vtx_from_face_and_edge
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
);

/**
 *
 * \brief Generate face_vtx with face_edge and edge_vtx
 *
 * \param [in]  n_face        Number of partition in current process
 * \param [in]  face_edge_idx Connectivity index between face and edge (size = n_face+1)
 * \param [in]  face_edge     Connectivity between face and edge (size = face_edge_idx[n_face])
 * \param [in]  edge_vtx      Connectivity between edge and vtx - Implicit (size = 2 * n_edge)
 * \param [out] face_vtx      Connectivity between entity1 and entity2 (size = entity1_entity2_idx[n_entity1] )
 *
 */
void
PDM_compute_face_vtx_from_face_and_edge_unsigned
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
);


/**
 *
 * \brief Generate face_vtx with face_edge and edge_vtx
 *
 * \param [in]  comm            PDM_MPI communicator
 * \param [in]  distrib_face    Distribution of faces among process (size = n_rank+1)
 * \param [in]  distrib_edge    Distribution of faces among process (size = n_rank+1)
 * \param [in]  dface_edge_idx  Connectivity index between face and edge (size = dn_face )
 * \param [in]  dface_edge      Connectivity between face and edge (size = dface_edge_idx[dn_face])
 * \param [out] dface_vtx       Connectivity between face and vtx (size = dface_edge_idx[dn_face])
 *
 */
void
PDM_compute_dface_vtx_from_edges_distrib
(
 PDM_MPI_Comm   comm,
 PDM_g_num_t   *distrib_face,
 PDM_g_num_t   *distrib_edge,
 int           *dface_edge_idx,
 PDM_g_num_t   *dface_edge,
 PDM_g_num_t   *dedge_vtx,
 PDM_g_num_t  **dface_vtx
);

/**
 *
 * \brief Generate face_vtx with face_edge and edge_vtx
 *
 * \param [in]  comm            PDM_MPI communicator
 * \param [in]  dn_face         Number of distributed faces
 * \param [in]  dn_edge         Number of distributed edges
 * \param [in]  dface_edge_idx  Connectivity index between face and edge (size = dn_face )
 * \param [in]  dface_edge      Connectivity between face and edge (size = dface_edge_idx[dn_face])
 * \param [out] dface_vtx       Connectivity between entity1 and entity2 (size = dface_edge_idx[dn_face])
 *
 */
void
PDM_compute_dface_vtx_from_edges
(
 PDM_MPI_Comm   comm,
 int            dn_face,
 int            dn_edge,
 int           *dface_edge_idx,
 PDM_g_num_t   *dface_edge,
 PDM_g_num_t   *dedge_vtx,
 PDM_g_num_t  **dface_vtx
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_PART_CONNECTIVITY_TRANSFORM_H_ */
