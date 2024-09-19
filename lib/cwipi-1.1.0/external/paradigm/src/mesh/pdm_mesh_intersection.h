/*
 * \file
 */

#ifndef PDM_MESH_INTERSECTION_H
#define PDM_MESH_INTERSECTION_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_to_part.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_extract_part.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_mesh_intersection_t PDM_mesh_intersection_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a mesh intersection structure
 *
 * \param [in]   dim_mesh_a     Mesh A dimension
 * \param [in]   dim_mesh_b     Mesh B dimension
 * \param [in]   project_coeff  Projection coefficient
 * \param [in]   comm           Associated communicator
 * \param [in]   owner          Results Ownership 
 *
 * \return       mi             Pointer to \ref PDM_mesh_intersection object
 *
 */


PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm,
 const PDM_ownership_t              owner
);


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 * \param [in]   i_mesh         Mesh identifier
 * \param [in]   n_part         Number of partitions
 *
 */

void
PDM_mesh_intersection_n_part_set
(
  PDM_mesh_intersection_t *mi,
  const int                i_mesh,
  const int                n_part
);


/**
 *
 * \brief  Compute meshes intersection
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 *
 */

void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
);

/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 * \param [in]   i_mesh         Mesh identifier
 * \param [in]   mesh           Pointer to \ref PDM_part_mesh_nodal object
 *
 */

void
PDM_mesh_intersection_mesh_nodal_set
(
 PDM_mesh_intersection_t  *mi,
 int                       i_mesh,
 PDM_part_mesh_nodal_t    *mesh
 );


/**
 *
 * \brief Set a mesh partition  
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 * \param [in]   i_mesh         Mesh identifier
 * \param [in]   i_part         Partition identifier
 * \param [in]   n_cell         Number of cells
 * \param [in]   n_face         Number of faces
 * \param [in]   n_edge         Number of edges
 * \param [in]   n_vtx          Number of vertices
 * \param [in]   cell_face_idx  Index in the cell -> face connectivity
 * \param [in]   cell_face      Cell -> face connectivity
 * \param [in]   face_edge_idx  Index in the face -> edge connectivity
 * \param [in]   face_edge      Cell -> face connectivity
 * \param [in]   edge_vtx       Edge -> vertex conectivity 
 * \param [in]   face_vtx_idx   Index in the face -> vertex connectivity
 * \param [in]   face_vtx       Face -> vertex connectivity
 * \param [in]   cell_ln_to_gn  Local cell numbering to global cel numbering
 * \param [in]   face_ln_to_gn  Local face numbering to global face numbering
 * \param [in]   edge_ln_to_gn  Local edge numbering to global edge numbering
 * \param [in]   vtx_ln_to_gn   Local vertex numbering to global vertex numbering
 * \param [in]   vtx_coord      Vertex coordinates
 *
 */

void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  int                       i_mesh,
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

void
PDM_mesh_intersection_tetraisation_pt_set
(
 PDM_mesh_intersection_t* mi,
 int     tetraisation_pt_type,
 double *tetraisation_pt_coord
);

void
PDM_mesh_intersection_stat_get
(
 PDM_mesh_intersection_t* mi,
 double *local_vol_A_B,
 double *global_vol_A_B,
 double *global_vol_A

);


/**
 *
 * \brief A mesh intersection structure  
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 */

void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
);


/**
 * \brief Get part_to_part object to exchange data between the intersected meshes
 *
 * \param [in ] mi         Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_intersection_part_to_part_get
(
 PDM_mesh_intersection_t  *mi,
 PDM_part_to_part_t      **ptp,
 PDM_ownership_t           ownership
 );


/**
 * \brief Get intersection result for the a point of view
 *
 * \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ipart              Partition index
 * \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
 * \param [out] elt_a_elt_b        List of intersected B elements for each A element 
 * \param [out] elt_a_elt_b_volume Volume of each intersection 
 *
 */

void
PDM_mesh_intersection_result_from_a_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       ipart,
       int                     **elt_a_elt_b_idx,
       PDM_g_num_t             **elt_a_elt_b,
       double                  **elt_a_elt_b_volume
);



/**
 * \brief Get intersection result for the b point of view
 * 
 * TODO: Do as PDM_mesh_intersection_result_from_a_get
 *
 * \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ipart              Partition index
 * \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
 * \param [out] elt_a_elt_b        List of intersected B elements for each A element 
 * \param [out] elt_a_elt_b_volume Volume of each intersection 
 *
 */

void
PDM_mesh_intersection_result_from_b_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       ipart,
       double                  **elt_b_elt_a_volume
 );


/**
 * \brief Get intersection result for the b point of view
 * 
 * TODO: Do as PDM_mesh_intersection_result_from_a_get
 *
 * \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ipart              Partition index
 * \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
 * \param [out] elt_a_elt_b        List of intersected B elements for each A element 
 * \param [out] elt_a_elt_b_volume Volume of each intersection 
 *
 */

void
PDM_mesh_intersection_elt_volume_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       imesh,
 const int                       ipart,
       double                  **elt_volume
 );


/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   mi              Pointer to \ref PDM_mesh_intersection object
 * \param [in]   tol             Tolerance
 *
 */
void
PDM_mesh_intersection_tolerance_set
(
       PDM_mesh_intersection_t *mi,
 const double                   tol
);


/**
 * \brief Get preprocessing results 
 *
 * \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] elt_a_elt_b_idx    Index of list of intersected B element candidate for each A element
 *                                 in the extr_mesh distribution 
 * \param [out] elt_a_elt_b        List of intersected B element candidate for each A element in the 
 *                                 extr_mesh distribution 
 * \param [out] extr_mesh_a        Redistributed mesh A with only A element candidate  
 * \param [out] extr_mesh_b        Redistributed mesh B with only B element candidate  
 *
 */

void
PDM_mesh_intersection_preprocessing_get
(
  PDM_mesh_intersection_t  *mi,
  int                     **box_a_box_b_idx,
  int                     **box_a_box_b,
  PDM_extract_part_t      **extr_mesh_a,
  PDM_extract_part_t      **extr_mesh_b
);


/**
 *
 * \brief Get mesh dimension
 *
 * \param [in] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [in] i_mesh             Mesh identifier
 *
 * \return Dimension of mesh \p i_mesh
 */

int
PDM_mesh_intersection_mesh_dimension_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       i_mesh
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_MESH_INTERSECTION_H */
