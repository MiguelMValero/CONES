/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_part_mesh_nodal_tetra_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{
  const int n_face_elt        = 4;
  // const int n_sum_vtx_face    = 12;
  int n_sum_vtx_elt           = 4;

  if(order > 1) {
    n_sum_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRAHO, order);
  }

  int __parent_node[4] = {0, 1, 2, 3};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;
    if (parent_num != NULL) {
      icell = parent_num[ielt];
    }

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];

    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    elmt_face_vtx[idx + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    elmt_face_vtx[idx + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    elmt_face_vtx[idx + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;

}

/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_part_mesh_nodal_prism_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{
  const int n_face_elt     = 5;
  // const int n_sum_vtx_face = 3*4 + 2*3;
  int n_sum_vtx_elt        = 6;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISMHO, order);
  }

  int __parent_node[6] = {0, 1, 2, 3, 4, 5};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;
    if (parent_num != NULL) {
      icell = parent_num[ielt];
    }

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];

    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    elmt_face_vtx[idx + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    elmt_face_vtx[idx + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    elmt_face_vtx[idx + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    elmt_face_vtx[idx + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    elmt_face_vtx[idx + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    elmt_face_vtx[idx + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    elmt_face_vtx[idx + 16] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 17] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
  }

  // *n_elt_current  += n_elt;
  // *n_face_current += n_elt * n_face_elt;

}



/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_part_mesh_nodal_pyra_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{

  const int n_face_elt     = 5;
  // const int n_sum_vtx_face = 1*4 + 4*3;
  int n_sum_vtx_elt        = 5;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PYRAMIDHO, order);
  }

  int __parent_node[5] = {0, 1, 2, 3, 4};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;
    if (parent_num != NULL) {
      icell = parent_num[ielt];
    }

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];

    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    elmt_face_vtx[idx + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    elmt_face_vtx[idx + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    elmt_face_vtx[idx + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    elmt_face_vtx[idx + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
  }

  // *n_elt_current  += n_elt;
  // *n_face_current += n_elt * n_face_elt;

}

void
PDM_part_mesh_nodal_hexa_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{
  const int n_face_elt     = 6;
  // const int n_sum_vtx_face = 24;
  int n_sum_vtx_elt        = 8;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXAHO, order);
  }

  int __parent_node[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */


  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;
    if (parent_num != NULL) {
      icell = parent_num[ielt];
    }

    /* Store the face_cell */
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];
    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    elmt_face_vtx[idx + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    elmt_face_vtx[idx + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    elmt_face_vtx[idx + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    elmt_face_vtx[idx + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    elmt_face_vtx[idx + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    elmt_face_vtx[idx + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    elmt_face_vtx[idx + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    elmt_face_vtx[idx + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    elmt_face_vtx[idx + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    elmt_face_vtx[idx + 16] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 17] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    elmt_face_vtx[idx + 18] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    elmt_face_vtx[idx + 19] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    elmt_face_vtx[idx + 20] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 21] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    elmt_face_vtx[idx + 22] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    elmt_face_vtx[idx + 23] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
  }


  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;
}


void
PDM_part_mesh_nodal_tri_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{

  const int n_face_elt        = 1;
  // const int n_sum_vtx_face    = 3;
  int n_sum_vtx_elt           = 3;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
  }

  int __parent_node[3] = {0, 1, 2};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */


  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;
    if (parent_num != NULL) {
      icell = parent_num[ielt];
    }

    /* Store the face_cell */
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];
    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

  }


  // *n_elt_current  += n_elt;
  // *n_face_current += n_elt * n_face_elt;
}


void
PDM_part_mesh_nodal_quad_decomposes_faces
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_face_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const int         *parent_num,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_face_cell,
       int         *parent_elmt_position
)
{

  const int n_face_elt        = 1;
  // const int n_sum_vtx_face    = 4;
  int n_sum_vtx_elt           = 4;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);
  }

  int __parent_node[4] = {0, 1, 2, 3};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  PDM_UNUSED(n_face_current);
  PDM_UNUSED(n_elt_current);
  // int _n_face_current = *n_face_current;
  // int _n_elt_current  = *n_elt_current;
  // int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  // PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  // int         *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  // int         *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;
  // PDM_g_num_t *_elmt_face_cell            = elmt_face_cell       + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    int icell = ielt;

    /* Store the face_cell */
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
      // _parent_elmt_position     [ielt * n_face_elt + i_face    ] = i_face;
      // _elmt_face_cell           [ielt * n_face_elt + i_face    ] = elmt_ln_to_gn[ielt];
      int idx_face = elmt_cell_face_idx[icell] + i_face;
      parent_elmt_position[idx_face] = i_face;
      elmt_face_cell      [idx_face] = elmt_ln_to_gn[ielt];
    }

    // _elmt_cell_face_idx[ielt+1] = _elmt_cell_face_idx[ielt] + n_face_elt;
    int idx = elmt_face_vtx_idx[elmt_cell_face_idx[icell]];
    elmt_face_vtx[idx + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    elmt_face_vtx[idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    elmt_face_vtx[idx + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    elmt_face_vtx[idx + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

  }


  // *n_elt_current  += n_elt;
  // *n_face_current += n_elt * n_face_elt;
}

void
PDM_part_mesh_nodal_bar_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{

  const int n_edge_elt        = 1;
  const int n_sum_vtx_edge    = 2;
  int n_sum_vtx_elt           = 2;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);
  }

  int __parent_node[2] = {0, 1};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }


  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}



void
PDM_part_mesh_nodal_tri_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{
  const int n_edge_elt        = 3;
  const int n_sum_vtx_edge    = 6;
  int n_sum_vtx_elt           = 3;

  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
  }

  int __parent_node[3] = {0, 1, 2};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}



void
PDM_part_mesh_nodal_quad_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{
  const int n_edge_elt        = 4;
  const int n_sum_vtx_edge    = 8;
  int n_sum_vtx_elt           = 4;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);
  }

  int __parent_node[4] = {0, 1, 2, 3};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}


void
PDM_part_mesh_nodal_tetra_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{
  const int n_edge_elt        = 6;
  const int n_sum_vtx_edge    = 12;
  int n_sum_vtx_elt     = 4;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRAHO, order);
  }

  int __parent_node[4] = {0, 1, 2, 3};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    // E3 = N3 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    // E4 = N1 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    // E5 = N2 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    // E6 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}


void
PDM_part_mesh_nodal_prism_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{

  const int n_edge_elt     = 9;
  const int n_sum_vtx_edge = 18;
  int n_sum_vtx_elt        = 6;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISMHO, order);
  }

  int __parent_node[6] = {0, 1, 2, 3, 4, 5};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    // E3 = N3 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    // E4 = N1 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    // E5 = N2 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E6 = N3 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    // E7 = N4 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E8 = N5 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    // E9 = N6 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 16] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 17] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}


void
PDM_part_mesh_nodal_pyra_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{

  const int n_edge_elt     = 8;
  const int n_sum_vtx_edge = 16;
  int n_sum_vtx_elt        = 5;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PYRAMIDHO, order);
  }

  int __parent_node[5] = {0, 1, 2, 3, 4};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    // E3 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    // E4 = N4 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    // E5 = N1 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E6 = N2 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E7 = N3 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E8 = N4 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15] = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}


void
PDM_part_mesh_nodal_hexa_decomposes_edges
(
       int          n_elt,
       int          order,
       int         *parent_node,
       int         *n_elt_current,
       int         *n_edge_current,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int         *connectivity_elmt_vtx,
 const PDM_g_num_t *elmt_ln_to_gn,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_edge_cell,
       int         *parent_elmt_position
)
{

  const int n_edge_elt     = 12;
  const int n_sum_vtx_edge = 24; // 2 vtx * 12 edge
  int n_sum_vtx_elt        = 8;
  if(order > 1) {
    n_sum_vtx_elt     = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXAHO, order);
  }

  int __parent_node[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = __parent_node;
  } else {
    _parent_node = parent_node;
  }

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _parent_elmt_position     [ielt * n_edge_elt + i_edge    ] = i_edge;
      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];
    }

    _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];

    // E3 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];

    // E4 = N4 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];

    // E5 = N1 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[0]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

    // E6 = N2 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[1]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    // E7 = N3 N7
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[2]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];

    // E8 = N4 N8
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[3]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];

    // E9 = N5 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 16]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 17]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];

    // E10 = N6 N7
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 18]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[5]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 19]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];

    // E11 = N7 N8
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 20]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[6]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 21]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];

    // E12 = N8 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 22]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[7]]-1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 23]  = vtx_ln_to_gn[connectivity_elmt_vtx[n_sum_vtx_elt * ielt + _parent_node[4]]-1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}

void
PDM_part_mesh_nodal_std_decomposes_faces
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_face_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
 const int                  *parent_num,
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_face_vtx_idx,
       PDM_g_num_t          *elmt_face_vtx,
       int                  *elmt_cell_face_idx,
       PDM_g_num_t          *elmt_face_cell,
       int                  *parent_elmt_position
)
{
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     abort();
     break;
   case PDM_MESH_NODAL_BAR2:
   case PDM_MESH_NODAL_BARHO:
   case PDM_MESH_NODAL_BARHO_BEZIER:
     abort();
     break;
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
     PDM_part_mesh_nodal_tri_decomposes_faces(n_elt,
                                              order,
                                              parent_node,
                                              n_elt_current,
                                              n_face_current,
                                              vtx_ln_to_gn,
                                              connectivity_elmt_vtx,
                                              parent_num,
                                              elmt_ln_to_gn,
                                              elmt_face_vtx_idx,
                                              elmt_face_vtx,
                                              elmt_cell_face_idx,
                                              elmt_face_cell,
                                              parent_elmt_position);
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     PDM_part_mesh_nodal_quad_decomposes_faces(n_elt,
                                               order,
                                               parent_node,
                                               n_elt_current,
                                               n_face_current,
                                               vtx_ln_to_gn,
                                               connectivity_elmt_vtx,
                                               parent_num,
                                               elmt_ln_to_gn,
                                               elmt_face_vtx_idx,
                                               elmt_face_vtx,
                                               elmt_cell_face_idx,
                                               elmt_face_cell,
                                               parent_elmt_position);
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     PDM_part_mesh_nodal_tetra_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                parent_num,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     PDM_part_mesh_nodal_pyra_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                parent_num,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     PDM_part_mesh_nodal_prism_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                parent_num,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     PDM_part_mesh_nodal_hexa_decomposes_faces(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_face_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                parent_num,
                                                elmt_ln_to_gn,
                                                elmt_face_vtx_idx,
                                                elmt_face_vtx,
                                                elmt_cell_face_idx,
                                                elmt_face_cell,
                                                parent_elmt_position);
     break;
   default:
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is supported\n");
  }
}


void
PDM_part_mesh_nodal_poly2d_decomposes_edges
(
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
 const int                  *connectivity_elmt_vtx_idx,
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_edge_vtx_idx,
       PDM_g_num_t          *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       PDM_g_num_t          *elmt_edge_cell,
       int                  *parent_elmt_position
)
{

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    // Reminder for poly2d -> Number of vertex = Number of edge
    int n_edge_elt = connectivity_elmt_vtx_idx[ielt+1] - connectivity_elmt_vtx_idx[ielt];
    *n_edge_current += n_edge_elt;
    int idx2 = connectivity_elmt_vtx_idx[ielt];
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[idx + 1] = _current_elmt_edge_vtx_idx[idx] + 2;
      _parent_elmt_position     [idx    ] = i_edge;

      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];

      _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

      int inext = (i_edge + 1) % n_edge_elt;
      _current_elmt_edge_vtx[2 * idx    ]  = vtx_ln_to_gn[connectivity_elmt_vtx[idx2 + i_edge]-1];
      _current_elmt_edge_vtx[2 * idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[idx2 + inext ]-1];

      idx += 1;
    }
  }


  *n_elt_current  += n_elt;

}


void
PDM_part_mesh_nodal_std_decomposes_edges
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_edge_vtx_idx,
       PDM_g_num_t          *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       PDM_g_num_t          *elmt_edge_cell,
       int                  *parent_elmt_position
)
{
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     abort();
     break;
   case PDM_MESH_NODAL_BAR2:
   case PDM_MESH_NODAL_BARHO:
   case PDM_MESH_NODAL_BARHO_BEZIER:
     PDM_part_mesh_nodal_bar_decomposes_edges(n_elt,
                                              order,
                                              parent_node,
                                              n_elt_current,
                                              n_edge_current,
                                              vtx_ln_to_gn,
                                              connectivity_elmt_vtx,
                                              elmt_ln_to_gn,
                                              elmt_edge_vtx_idx,
                                              elmt_edge_vtx,
                                              elmt_cell_edge_idx,
                                              elmt_edge_cell,
                                              parent_elmt_position);
     break;
   case PDM_MESH_NODAL_TRIA3:
   case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER:
     PDM_part_mesh_nodal_tri_decomposes_edges(n_elt,
                                              order,
                                              parent_node,
                                              n_elt_current,
                                              n_edge_current,
                                              vtx_ln_to_gn,
                                              connectivity_elmt_vtx,
                                              elmt_ln_to_gn,
                                              elmt_edge_vtx_idx,
                                              elmt_edge_vtx,
                                              elmt_cell_edge_idx,
                                              elmt_edge_cell,
                                              parent_elmt_position);
     break;
   case PDM_MESH_NODAL_QUAD4:
   case PDM_MESH_NODAL_QUADHO:
     PDM_part_mesh_nodal_quad_decomposes_edges(n_elt,
                                               order,
                                               parent_node,
                                               n_elt_current,
                                               n_edge_current,
                                               vtx_ln_to_gn,
                                               connectivity_elmt_vtx,
                                               elmt_ln_to_gn,
                                               elmt_edge_vtx_idx,
                                               elmt_edge_vtx,
                                               elmt_cell_edge_idx,
                                               elmt_edge_cell,
                                               parent_elmt_position);
     break;
   case PDM_MESH_NODAL_TETRA4:
   case PDM_MESH_NODAL_TETRAHO:
     PDM_part_mesh_nodal_tetra_decomposes_edges(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_edge_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_edge_vtx_idx,
                                                elmt_edge_vtx,
                                                elmt_cell_edge_idx,
                                                elmt_edge_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PYRAMID5:
   case PDM_MESH_NODAL_PYRAMIDHO:
     PDM_part_mesh_nodal_pyra_decomposes_edges(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_edge_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_edge_vtx_idx,
                                                elmt_edge_vtx,
                                                elmt_cell_edge_idx,
                                                elmt_edge_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_PRISM6:
   case PDM_MESH_NODAL_PRISMHO:
     PDM_part_mesh_nodal_prism_decomposes_edges(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_edge_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_edge_vtx_idx,
                                                elmt_edge_vtx,
                                                elmt_cell_edge_idx,
                                                elmt_edge_cell,
                                                parent_elmt_position);
     break;
   case PDM_MESH_NODAL_HEXA8:
   case PDM_MESH_NODAL_HEXAHO:
     PDM_part_mesh_nodal_hexa_decomposes_edges(n_elt,
                                                order,
                                                parent_node,
                                                n_elt_current,
                                                n_edge_current,
                                                vtx_ln_to_gn,
                                                connectivity_elmt_vtx,
                                                elmt_ln_to_gn,
                                                elmt_edge_vtx_idx,
                                                elmt_edge_vtx,
                                                elmt_cell_edge_idx,
                                                elmt_edge_cell,
                                                parent_elmt_position);
     break;
   default:
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is supported\n");
  }
}


void
PDM_part_mesh_nodal_elmts_sections_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_face_vtx_idx,
  PDM_g_num_t                  *elmt_face_vtx,
  int                          *elmt_cell_face_idx,
  PDM_g_num_t                  *elmt_face_cell,
  int                          *parent_elmt_position
)
{

  int  n_section  = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *section_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int n_elt_current  = 0;
  int n_face_current = 0;

  int parent_node[8];

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < n_section; i_section++) {
      int id_section = section_id[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

      switch (t_elt) {
        case PDM_MESH_NODAL_POINT:
        case PDM_MESH_NODAL_BAR2:
        case PDM_MESH_NODAL_TRIA3:
        case PDM_MESH_NODAL_QUAD4:
        case PDM_MESH_NODAL_TETRA4:
        case PDM_MESH_NODAL_PYRAMID5:
        case PDM_MESH_NODAL_PRISM6:
        case PDM_MESH_NODAL_HEXA8:
        {
          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    id_section,
                                                    i_part,
                                                    &connec,
                                                    &numabs,
                                                    &parent_num,
                                                    &parent_entity_g_num,
                                                    PDM_OWNERSHIP_KEEP);


          // Is this legit??
          PDM_g_num_t *g_num = numabs;
          if (numabs == NULL) {
            g_num = parent_entity_g_num;
          }

          assert(g_num != NULL);

          PDM_part_mesh_nodal_std_decomposes_faces(t_elt,
                                                   n_elt,
                                                   1,
                                                   NULL,
                                                   &n_elt_current,
                                                   &n_face_current,
                                                   vtx_ln_to_gn[i_part],
                                                   connec,
                                                   parent_num,
                                                   g_num,
                                                   elmt_face_vtx_idx,
                                                   elmt_face_vtx,
                                                   elmt_cell_face_idx,
                                                   elmt_face_cell,
                                                   parent_elmt_position);

          break;
        }
        case PDM_MESH_NODAL_BARHO:
        case PDM_MESH_NODAL_BARHO_BEZIER:
        case PDM_MESH_NODAL_TRIAHO:
        case PDM_MESH_NODAL_TRIAHO_BEZIER:
        case PDM_MESH_NODAL_QUADHO:
        case PDM_MESH_NODAL_TETRAHO:
        case PDM_MESH_NODAL_PYRAMIDHO:
        case PDM_MESH_NODAL_PRISMHO:
        case PDM_MESH_NODAL_HEXAHO:
        {
          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          int          order               = -1;
          const char  *ho_ordering         = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &connec,
                                                       &numabs,
                                                       &parent_num,
                                                       &parent_entity_g_num,
                                                       &order,
                                                       &ho_ordering,
                                                       PDM_OWNERSHIP_KEEP);
          PDM_Mesh_nodal_ho_parent_node(t_elt,
                                        order,
                                        ho_ordering,
                                        parent_node);

          PDM_part_mesh_nodal_std_decomposes_faces(t_elt,
                                                   n_elt,
                                                   order,
                                                   parent_node,
                                                   &n_elt_current,
                                                   &n_face_current,
                                                   vtx_ln_to_gn[i_part],
                                                   connec,
                                                   parent_num,
                                                   numabs,
                                                   elmt_face_vtx_idx,
                                                   elmt_face_vtx,
                                                   elmt_cell_face_idx,
                                                   elmt_face_cell,
                                                   parent_elmt_position);
        break;
        }
        case PDM_MESH_NODAL_POLY_2D:
        {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_faces : Element type is supported\n");
          break;
        }

        case PDM_MESH_NODAL_POLY_3D:
        {
          int          n_face              = 0;
          PDM_g_num_t *face_ln_to_gn       = NULL;
          int         *face_vtx_idx        = NULL;
          int         *face_vtx            = NULL;
          PDM_g_num_t *cell_ln_to_gn       = NULL;
          int         *cell_face_idx       = NULL;
          int         *cell_face           = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;

          PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &n_face,
                                                       &face_ln_to_gn,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       &cell_ln_to_gn,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       &parent_num,
                                                       &parent_entity_g_num,
                                                       PDM_OWNERSHIP_KEEP);

          // Is this legit??
          PDM_g_num_t *g_num = cell_ln_to_gn;
          if (cell_ln_to_gn == NULL) {
            g_num = parent_entity_g_num;
          }

          assert(g_num != NULL);

          for (int ielt = 0; ielt < n_elt; ielt++) {
            int icell = ielt;
            if (parent_num != NULL) {
              icell = parent_num[ielt];
            }

            int cell_face_n = cell_face_idx[ielt+1] - cell_face_idx[ielt];
            int *_cell_face = cell_face + cell_face_idx[ielt];

            for (int iface = 0; iface < cell_face_n; iface++) {

              int idx_write_face = elmt_cell_face_idx[icell] + iface;

              parent_elmt_position[idx_write_face] = iface;
              elmt_face_cell      [idx_write_face] = g_num[ielt];

              PDM_g_num_t *_elmt_face_vtx = elmt_face_vtx + elmt_face_vtx_idx[idx_write_face];

              int face_id = PDM_ABS(_cell_face[iface]) - 1;

              int face_vtx_n = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
              int *_face_vtx = face_vtx + face_vtx_idx[face_id];

              for (int ivtx = 0; ivtx < face_vtx_n; ivtx++) {
                int vtx_id = -1;
                if (_cell_face[iface] > 0) {
                  vtx_id = _face_vtx[ivtx];
                }
                else {
                  vtx_id = _face_vtx[face_vtx_n-1-ivtx];
                }
                _elmt_face_vtx[ivtx] = vtx_ln_to_gn[i_part][vtx_id-1];
              }

            }
          }
          break;
        }

        default:
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_faces : Element type is supported\n");
      }
    }
  }
}


void
PDM_part_mesh_nodal_elmts_sections_decompose_edges
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_edge_vtx_idx,
  PDM_g_num_t                  *elmt_edge_vtx,
  int                          *elmt_cell_edge_idx,
  PDM_g_num_t                  *elmt_edge_cell,
  int                          *parent_elmt_position
)
{

  int  n_section  = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *section_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int n_elt_current  = 0;
  int n_edge_current = 0;

  int parent_node[8];

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < n_section; i_section++) {
      int id_section = section_id[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      switch (t_elt) {
        case PDM_MESH_NODAL_POINT:
        case PDM_MESH_NODAL_BAR2:
        case PDM_MESH_NODAL_TRIA3:
        case PDM_MESH_NODAL_QUAD4:
        case PDM_MESH_NODAL_TETRA4:
        case PDM_MESH_NODAL_PYRAMID5:
        case PDM_MESH_NODAL_PRISM6:
        case PDM_MESH_NODAL_HEXA8:
        {
          int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    id_section,
                                                    i_part,
                                                    &connec,
                                                    &numabs,
                                                    &parent_num,
                                                    &parent_entity_g_num,
                                                    PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_nodal_std_decomposes_edges(t_elt,
                                                   n_elt,
                                                   1,
                                                   NULL,
                                                   &n_elt_current,
                                                   &n_edge_current,
                                                   vtx_ln_to_gn[i_part],
                                                   connec,
                                                   numabs,
                                                   elmt_edge_vtx_idx,
                                                   elmt_edge_vtx,
                                                   elmt_cell_edge_idx,
                                                   elmt_edge_cell,
                                                   parent_elmt_position);

          break;
        }
        case PDM_MESH_NODAL_BARHO:
        case PDM_MESH_NODAL_BARHO_BEZIER:
        case PDM_MESH_NODAL_TRIAHO:
        case PDM_MESH_NODAL_TRIAHO_BEZIER:
        case PDM_MESH_NODAL_QUADHO:
        case PDM_MESH_NODAL_TETRAHO:
        case PDM_MESH_NODAL_PYRAMIDHO:
        case PDM_MESH_NODAL_PRISMHO:
        case PDM_MESH_NODAL_HEXAHO:
        {
          int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *parent_num          = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          int          order               = -1;
          const char* ho_ordering          = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &connec,
                                                       &numabs,
                                                       &parent_num,
                                                       &parent_entity_g_num,
                                                       &order,
                                                       &ho_ordering,
                                                       PDM_OWNERSHIP_KEEP);
          PDM_Mesh_nodal_ho_parent_node(t_elt,
                                        order,
                                        ho_ordering,
                                        parent_node);

          PDM_part_mesh_nodal_std_decomposes_edges(t_elt,
                                                   n_elt,
                                                   order,
                                                   parent_node,
                                                   &n_elt_current,
                                                   &n_edge_current,
                                                   vtx_ln_to_gn[i_part],
                                                   connec,
                                                   numabs,
                                                   elmt_edge_vtx_idx,
                                                   elmt_edge_vtx,
                                                   elmt_cell_edge_idx,
                                                   elmt_edge_cell,
                                                   parent_elmt_position);
        break;
        }
        case PDM_MESH_NODAL_POLY_2D:
        {
          int         *connec              = NULL;
          int         *connec_idx          = NULL;

          int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);
          // int         *PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
          //                                                       id_section,
          //                                                       i_part);

          PDM_g_num_t *numabs = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                    id_section,
                                                                    i_part,
                                                                    PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &connec_idx,
                                                       &connec,
                                                       PDM_OWNERSHIP_KEEP);

          PDM_part_mesh_nodal_poly2d_decomposes_edges(n_elt,
                                                      &n_elt_current,
                                                      &n_edge_current,
                                                      vtx_ln_to_gn[i_part],
                                                      connec,
                                                      connec_idx,
                                                      numabs,
                                                      elmt_edge_vtx_idx,
                                                      elmt_edge_vtx,
                                                      elmt_cell_edge_idx,
                                                      elmt_edge_cell,
                                                      parent_elmt_position);

          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_edges : Element type is supported\n");
          break;
        }

        case PDM_MESH_NODAL_POLY_3D:
        {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_edges : Element type is supported\n");
          break;
        }

        default:
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_part_mesh_nodal_elmts_sections_decompose_edges : Element type is supported\n");
      }
    }
  }
}



void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t  *pmne,
 int                          *n_elt_tot,
 int                          *n_face_elt_tot,
 int                          *n_sum_vtx_face_tot,
 int                         **elmt_face_vtx_idx,
 int                         **elmt_cell_face_idx
)
{
  /* Get current structure to treat */
  *n_face_elt_tot     = 0;
  *n_sum_vtx_face_tot = 0;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {

      int n_face_elt     = PDM_n_face_elt_per_elmt    (pmne->sections_std[i_section]->t_elt);
      int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(pmne->sections_std[i_section]->t_elt);

      *n_elt_tot          += pmne->sections_std[i_section]->n_elt[i_part];
      *n_face_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_face_elt;
      *n_sum_vtx_face_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_face;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      // int _n_face      = pmne->sections_poly3d[i_section]->n_face[i_part];
      // *n_face_elt_tot += _n_face;
      // int n_face_vtx   = pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      // *n_sum_vtx_face_tot += n_face_vtx;
      int n_elt = pmne->sections_poly3d[i_section]->n_elt[i_part];
      int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
      int *cell_face     = pmne->sections_poly3d[i_section]->_cellfac    [i_part];
      int *face_vtx_idx  = pmne->sections_poly3d[i_section]->_facvtx_idx [i_part];
      for (int i = 0; i < n_elt; i++) {
        for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
          (*n_face_elt_tot)++;
          int face_id = PDM_ABS(cell_face[j]) - 1;
          *n_sum_vtx_face_tot += face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
        }
      }

      *n_elt_tot += n_elt;
    }

    // Not so sure about this...
    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face      = pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_face_elt_tot += _n_face;
      int n_edge_vtx   = pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_edge_vtx;
      *n_elt_tot          += pmne->sections_poly2d[i_section]->n_elt[i_part];
    }

  }


  *elmt_face_vtx_idx  = PDM_array_zeros_int(*n_face_elt_tot + 1);
  *elmt_cell_face_idx = PDM_array_zeros_int(*n_elt_tot      + 1);


  int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              id_section,
                                                              i_part);


      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_error(__FILE__, __LINE__, 0, "TODO\n");
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int i_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
        int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
        for (int i = 0; i < n_elt; i++) {
          int icell = i;
          if (parent_num != NULL) {
            icell = parent_num[i];
          }
          (*elmt_cell_face_idx)[icell+1] = cell_face_idx[i+1] - cell_face_idx[i];
        }
      }
      else {
        int cell_face_n = PDM_n_face_elt_per_elmt(pmne->sections_std[isection]->t_elt);

        for (int i = 0; i < n_elt; i++) {
          (*elmt_cell_face_idx)[i+1] = cell_face_n;
        }
      }

    }

  }

  PDM_array_accumulate_int(*elmt_cell_face_idx + 1, *n_elt_tot);


  int cell_face_vtx_n[6];


  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    for (int isection = 0; isection < n_section; isection++) {


      int id_section = sections_id[isection];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              id_section,
                                                              i_part);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_error(__FILE__, __LINE__, 0, "TODO\n");
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int i_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
        int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
        int *cell_face     = pmne->sections_poly3d[i_section]->_cellfac    [i_part];
        int *face_vtx_idx  = pmne->sections_poly3d[i_section]->_facvtx_idx [i_part];
        for (int i = 0; i < n_elt; i++) {
          int icell = i;
          if (parent_num != NULL) {
            icell = parent_num[i];
          }
          int cell_face_n = cell_face_idx[i+1] - cell_face_idx[i];
          for (int j = 0; j < cell_face_n; j++) {
            int face_id = PDM_ABS(cell_face[cell_face_idx[i]+j]) - 1;
            (*elmt_face_vtx_idx)[(*elmt_cell_face_idx)[icell]+j+1] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
          }
        }
      }
      else {
        int cell_face_n = PDM_n_face_elt_per_elmt(pmne->sections_std[isection]->t_elt);

        switch (t_elt) {
          case PDM_MESH_NODAL_TRIA3:
          case PDM_MESH_NODAL_TRIAHO: {
            cell_face_vtx_n[0] = 3;
            break;
          }
          case PDM_MESH_NODAL_QUAD4:
          case PDM_MESH_NODAL_QUADHO: {
            cell_face_vtx_n[0] = 4;
            break;
          }
          case PDM_MESH_NODAL_TETRA4:
          case PDM_MESH_NODAL_TETRAHO: {
            for (int i = 0; i < 4; i++) {
              cell_face_vtx_n[i] = 3;
            }
            break;
          }
          case PDM_MESH_NODAL_PYRAMID5:
          case PDM_MESH_NODAL_PYRAMIDHO: {
            cell_face_vtx_n[0] = 4;
            cell_face_vtx_n[1] = 3;
            cell_face_vtx_n[2] = 3;
            cell_face_vtx_n[3] = 3;
            cell_face_vtx_n[4] = 3;
            break;
          }
          case PDM_MESH_NODAL_PRISM6:
          case PDM_MESH_NODAL_PRISMHO: {
            cell_face_vtx_n[0] = 3;
            cell_face_vtx_n[1] = 3;
            cell_face_vtx_n[2] = 4;
            cell_face_vtx_n[3] = 4;
            cell_face_vtx_n[4] = 4;
            break;
          }
          case PDM_MESH_NODAL_HEXA8:
          case PDM_MESH_NODAL_HEXAHO: {
            for (int i = 0; i < 6; i++) {
              cell_face_vtx_n[i] = 4;
            }
            break;
          }
          default:
            PDM_error(__FILE__, __LINE__, 0, "Invalid elt type %d\n", (int) t_elt);
        }

        for (int i = 0; i < n_elt; i++) {

          for (int j = 0; j < cell_face_n; j++) {
            (*elmt_face_vtx_idx)[(*elmt_cell_face_idx)[i]+j+1] = cell_face_vtx_n[j];
          }

        }
      }

    }

  }


  PDM_array_accumulate_int(*elmt_face_vtx_idx + 1, *n_face_elt_tot);
  // printf("n_face_elt_tot     ::%i\n", *n_face_elt_tot   );
  // printf("n_sum_vtx_face_tot::%i\n" , *n_sum_vtx_face_tot);
}



void
PDM_part_mesh_nodal_elmts_decompose_edges_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                         *n_elt_tot,
 int                         *n_edge_elt_tot,
 int                         *n_sum_vtx_edge_tot
)
{
  /* Get current structure to treat */
  *n_edge_elt_tot     = 0;
  *n_sum_vtx_edge_tot = 0;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {

      int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (pmne->sections_std[i_section]->t_elt);
      int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(pmne->sections_std[i_section]->t_elt);

      *n_elt_tot          += pmne->sections_std[i_section]->n_elt[i_part];
      *n_edge_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_edge_elt;
      *n_sum_vtx_edge_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_edge;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face          = pmne->sections_poly3d[i_section]->n_face[i_part];
      *n_edge_elt_tot     +=     pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_sum_vtx_edge_tot += 2 * pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_elt_tot          += pmne->sections_poly3d[i_section]->n_elt[i_part];
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face          =     pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_edge_elt_tot     +=     pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_edge_tot += 2 * pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_elt_tot          +=     pmne->sections_poly2d[i_section]->n_elt[i_part];
    }
  }

  // printf("n_edge_elt_tot     ::%i\n", *n_edge_elt_tot   );
  // printf("n_sum_vtx_edge_tot::%i\n" , *n_sum_vtx_edge_tot);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
