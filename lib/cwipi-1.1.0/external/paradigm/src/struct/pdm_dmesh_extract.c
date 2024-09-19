/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_dmesh_extract.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elmts_priv.h"
#include "pdm_dmesh_extract_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_gnum_location.h"
#include "pdm_unique.h"
#include "pdm_dmesh_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_partitioning_nodal_algorithm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_rebuild_group
(
 PDM_dmesh_extract_t *dme,
 PDM_g_num_t         *distrib_entity,
 PDM_mesh_entities_t  entity_type,
 PDM_bound_type_t     bound_type
)
{

  int         *dbound_entity_idx = NULL;
  PDM_g_num_t *dbound_entity     = NULL;
  int n_group_entity = PDM_dmesh_bound_get(dme->dmesh,
                                           bound_type,
                                           &dbound_entity,
                                           &dbound_entity_idx,
                                           PDM_OWNERSHIP_BAD_VALUE);
  if(dbound_entity_idx == NULL) {
    return;
  }

  assert(dme->btp_bound_entity_to_extract_entity[bound_type] == NULL);
  assert(dme->btp_bound_ownership               [bound_type] == NULL);

  dme->btp_bound_entity_to_extract_entity[bound_type] = malloc(n_group_entity * sizeof(PDM_block_to_part_t *));
  dme->btp_bound_ownership               [bound_type] = malloc(n_group_entity * sizeof(PDM_ownership_t      ));

  dme->n_bound[bound_type] = n_group_entity;

  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  int dn_entity = dme->distrib_extract[entity_type][i_rank+1] - dme->distrib_extract[entity_type][i_rank];

  int         **tmp_dextract_bound_entity_idx      = NULL;
  int         **tmp_dextract_bound_entity          = NULL;
  PDM_g_num_t **tmp_dextract_bound_entity_ln_to_gn = NULL;
  PDM_part_distgroup_to_partgroup(dme->comm,
                                  distrib_entity,
                                  n_group_entity,
                                  dbound_entity_idx,
                                  dbound_entity,
                                  1,
                                  &dn_entity,
          (const PDM_g_num_t **)  &dme->parent_extract_gnum[entity_type],
                                  &tmp_dextract_bound_entity_idx,
                                  &tmp_dextract_bound_entity,
                                  &tmp_dextract_bound_entity_ln_to_gn);

  int         *dextract_bound_entity_idx      = tmp_dextract_bound_entity_idx     [0];
  int         *dextract_bound_entity          = tmp_dextract_bound_entity         [0];
  PDM_g_num_t *dextract_bound_entity_ln_to_gn = tmp_dextract_bound_entity_ln_to_gn[0];
  free(tmp_dextract_bound_entity_idx     );
  free(tmp_dextract_bound_entity         );
  free(tmp_dextract_bound_entity_ln_to_gn);

  if(0 == 1) {
    PDM_log_trace_array_int (dextract_bound_entity         , dextract_bound_entity_idx[n_group_entity], "dextract_bound_entity ::");
    PDM_log_trace_array_long(dextract_bound_entity_ln_to_gn, dextract_bound_entity_idx[n_group_entity], "dextract_bound_entity_ln_to_gn ::");
  }

  /*
   * Rebuild all btp for each group
   */
  for(int i_group = 0; i_group < n_group_entity; ++i_group) {

    int dn_group_entity = dbound_entity_idx[i_group+1] - dbound_entity_idx[i_group];
    PDM_g_num_t* distrib_group_entity = PDM_compute_entity_distribution(dme->comm, dn_group_entity);

    int beg = dextract_bound_entity_idx[i_group];
    int dn_extract_group_entity = dextract_bound_entity_idx[i_group+1] - beg;
    PDM_g_num_t *_lextract_bound_entity_ln_to_gn = &dextract_bound_entity_ln_to_gn[beg];

    dme->btp_bound_entity_to_extract_entity[bound_type][i_group] = PDM_block_to_part_create(distrib_group_entity,
                                                                     (const PDM_g_num_t **) &_lextract_bound_entity_ln_to_gn,
                                                                                            &dn_extract_group_entity,
                                                                                            1,
                                                                                            dme->comm);


    // dextract_bound_entity + dme->parent_extract_gnum[entity_type] => dbound_entity_extract (PDM_g_num_t)
    // On peut également faire le transfert du gnum via le btp

    dme->btp_bound_ownership[bound_type][i_group] = PDM_OWNERSHIP_KEEP;
    free(distrib_group_entity);


  }


  /*
   * Reminder : block(view as part) is implicit
   */
  PDM_g_num_t *dextract_bound_entity_gnum = malloc(dextract_bound_entity_idx[n_group_entity] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dextract_bound_entity_idx[n_group_entity]; ++i) {
    dextract_bound_entity_gnum[i] = dme->distrib_extract[entity_type][i_rank] + dextract_bound_entity[i];
  }

  PDM_dmesh_bound_set(dme->dmesh_extract,
                      bound_type,
                      n_group_entity,
                      dextract_bound_entity_gnum,
                      dextract_bound_entity_idx,
                      PDM_OWNERSHIP_KEEP);


  // free(dextract_bound_entity_idx     );
  free(dextract_bound_entity         );
  free(dextract_bound_entity_ln_to_gn);
}




static
void
_dmesh_extract_3d
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  PDM_g_num_t *dcell_face     = NULL;
  int         *dcell_face_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_CELL_FACE,
                             &dcell_face,
                             &dcell_face_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  int from_face_edge = 0;
  int from_face_vtx  = 0;

  PDM_g_num_t *dface_vtx     = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_vtx_idx != NULL){
    from_face_vtx = 1;
  }

  // face edge
  PDM_g_num_t *dface_edge     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_edge_idx != NULL) {
    from_face_edge = 1;
  }

  int dn_cell = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_CELL);
  int dn_face = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_FACE);
  int dn_edge = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_EDGE);
  int dn_vtx  = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_VTX);

  PDM_g_num_t* distrib_cell = PDM_compute_entity_distribution(dme->comm, dn_cell);
  PDM_g_num_t* distrib_face = PDM_compute_entity_distribution(dme->comm, dn_face);
  PDM_g_num_t* distrib_vtx  = PDM_compute_entity_distribution(dme->comm, dn_vtx );
  PDM_g_num_t* distrib_edge = NULL;

  PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                             dme->n_selected,
                                             dme->selected_gnum,
                                             distrib_cell,
                                             dcell_face_idx,
                                             dcell_face,
                                             &dme->distrib_extract                 [PDM_MESH_ENTITY_CELL],
                                             &dme->parent_extract_gnum             [PDM_MESH_ENTITY_CELL],
                                             &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                             &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                             &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_CELL],
                                             &dme->distrib_extract                 [PDM_MESH_ENTITY_FACE],
                                             &dme->parent_extract_gnum             [PDM_MESH_ENTITY_FACE]);

  dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
  dme->dmesh_extract->dn_face = dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank];

  if(from_face_edge == 1) {

    PDM_dconnectivity_to_extract_dconnectivity_block(dme->comm,
                                                     dme->dmesh_extract->dn_face,
                                                     dme->parent_extract_gnum[PDM_MESH_ENTITY_FACE],
                                                     distrib_face,
                                                     dface_edge_idx,
                                                     dface_edge,
                                                     &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                     &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                     &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                                     &dme->distrib_extract                 [PDM_MESH_ENTITY_EDGE],
                                                     &dme->parent_extract_gnum             [PDM_MESH_ENTITY_EDGE]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;

    // edge_vtx
    distrib_edge = PDM_compute_entity_distribution(dme->comm, dn_edge);
    PDM_g_num_t *dedge_vtx     = NULL;
    int         *dedge_vtx_idx = NULL;
    PDM_dmesh_connectivity_get(dme->dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
    int *_dedge_vtx_idx = NULL;
    if(dedge_vtx_idx == NULL)  {
      _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
        _dedge_vtx_idx[i_edge] = 2*i_edge;
      }
    } else {
      _dedge_vtx_idx = dedge_vtx_idx;
    }
    dme->dmesh_extract->dn_edge = dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank];
    PDM_dconnectivity_to_extract_dconnectivity_block(dme->comm,
                                                     dme->dmesh_extract->dn_edge,
                                                     dme->parent_extract_gnum[PDM_MESH_ENTITY_EDGE],
                                                     distrib_edge,
                                                     _dedge_vtx_idx,
                                                     dedge_vtx,
                                                     &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_EDGE],
                                                     &dme->distrib_extract                 [PDM_MESH_ENTITY_VTX],
                                                     &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VTX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
    free(dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
    dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = NULL;
    if(dedge_vtx_idx == NULL)  {
      free(_dedge_vtx_idx);
    }

  } else if(from_face_vtx == 1) {
    PDM_dconnectivity_to_extract_dconnectivity_block(dme->comm,
                                                     dme->dmesh_extract->dn_face,
                                                     dme->parent_extract_gnum[PDM_MESH_ENTITY_FACE],
                                                     distrib_face,
                                                     dface_vtx_idx,
                                                     dface_vtx,
                                                     &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                     &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                     &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                                     &dme->distrib_extract                 [PDM_MESH_ENTITY_VTX],
                                                     &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VTX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_TRUE;
 }

  dme->dmesh_extract->dn_face = dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank];
  dme->dmesh_extract->dn_vtx  = dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank];
  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX] = PDM_block_to_part_create(distrib_vtx,
                                                             (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                    &dme->dmesh_extract->dn_vtx,
                                                                                    1,
                                                                                    dme->comm);


  _rebuild_group(dme, distrib_cell, PDM_MESH_ENTITY_CELL, PDM_BOUND_TYPE_CELL);
  _rebuild_group(dme, distrib_face, PDM_MESH_ENTITY_FACE, PDM_BOUND_TYPE_FACE);
  _rebuild_group(dme, distrib_edge, PDM_MESH_ENTITY_EDGE, PDM_BOUND_TYPE_EDGE);
  _rebuild_group(dme, distrib_vtx , PDM_MESH_ENTITY_VTX,  PDM_BOUND_TYPE_VTX );

  free(distrib_cell);
  free(distrib_face);
  if(distrib_edge != NULL) {
    free(distrib_edge);
  }
  free(distrib_vtx );




}


static
void
_dmesh_extract_2d
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  int from_face_edge = 0;
  int from_face_vtx  = 0;

  PDM_g_num_t *dface_vtx     = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_vtx_idx != NULL){
    from_face_vtx = 1;
  }

  // face edge
  PDM_g_num_t *dface_edge     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_edge_idx != NULL) {
    from_face_edge = 1;
  }

  int dn_face = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_FACE);
  int dn_edge = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_EDGE);
  int dn_vtx  = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_VTX);

  PDM_g_num_t* distrib_face = PDM_compute_entity_distribution(dme->comm, dn_face);
  PDM_g_num_t* distrib_vtx  = PDM_compute_entity_distribution(dme->comm, dn_vtx );
  PDM_g_num_t* distrib_edge = NULL;

  if(from_face_edge == 1) {

    PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                               dme->n_selected,
                                               dme->selected_gnum,
                                               distrib_face,
                                               dface_edge_idx,
                                               dface_edge,
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_FACE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_FACE],
                                               &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                               &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                               &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_EDGE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_EDGE]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;

    // edge_vtx
    distrib_edge = PDM_compute_entity_distribution(dme->comm, dn_edge);
    PDM_g_num_t *dedge_vtx     = NULL;
    int         *dedge_vtx_idx = NULL;
    PDM_dmesh_connectivity_get(dme->dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
    int *_dedge_vtx_idx = NULL;
    if(dedge_vtx_idx == NULL)  {
      _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
        _dedge_vtx_idx[i_edge] = 2*i_edge;
      }
    } else {
      _dedge_vtx_idx = dedge_vtx_idx;
    }
    dme->dmesh_extract->dn_edge = dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank];
    PDM_dconnectivity_to_extract_dconnectivity_block(dme->comm,
                                                     dme->dmesh_extract->dn_edge,
                                                     dme->parent_extract_gnum[PDM_MESH_ENTITY_EDGE],
                                                     distrib_edge,
                                                     _dedge_vtx_idx,
                                                     dedge_vtx,
                                                     &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_EDGE],
                                                     &dme->distrib_extract                 [PDM_MESH_ENTITY_VTX],
                                                     &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VTX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
    free(dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
    dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = NULL;
    if(dedge_vtx_idx == NULL)  {
      free(_dedge_vtx_idx);
    }

  } else if(from_face_vtx == 1) {
    PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                               dme->n_selected,
                                               dme->selected_gnum,
                                               distrib_face,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_FACE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_FACE],
                                               &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                               &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                               &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_VTX],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VTX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_TRUE;
 }

  dme->dmesh_extract->dn_face = dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank];
  dme->dmesh_extract->dn_vtx  = dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank];
  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX] = PDM_block_to_part_create(distrib_vtx,
                                                             (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                    &dme->dmesh_extract->dn_vtx,
                                                                                    1,
                                                                                    dme->comm);


  _rebuild_group(dme, distrib_face, PDM_MESH_ENTITY_FACE, PDM_BOUND_TYPE_FACE);
  _rebuild_group(dme, distrib_edge, PDM_MESH_ENTITY_EDGE, PDM_BOUND_TYPE_EDGE);
  _rebuild_group(dme, distrib_vtx , PDM_MESH_ENTITY_VTX,  PDM_BOUND_TYPE_VTX );

  free(distrib_face);
  if(distrib_edge != NULL) {
    free(distrib_edge);
  }
  free(distrib_vtx );
}


static
void
_dmesh_extract_1d
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  int dn_edge = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_EDGE);
  int dn_vtx  = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_VTX );

  PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(dme->comm, dn_edge);
  PDM_g_num_t* distrib_vtx  = PDM_compute_entity_distribution(dme->comm, dn_vtx );

  // edge_vtx
  PDM_g_num_t *dedge_vtx     = NULL;
  int         *dedge_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             &dedge_vtx,
                             &dedge_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);
  int *_dedge_vtx_idx = NULL;
  if(dedge_vtx_idx == NULL)  {
    _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
    for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
      _dedge_vtx_idx[i_edge] = 2*i_edge;
    }
  } else {
    _dedge_vtx_idx = dedge_vtx_idx;
  }


  PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                             dme->n_selected,
                                             dme->selected_gnum,
                                             distrib_edge,
                                             dedge_vtx_idx,
                                             dedge_vtx,
                                             &dme->distrib_extract                 [PDM_MESH_ENTITY_EDGE],
                                             &dme->parent_extract_gnum             [PDM_MESH_ENTITY_EDGE],
                                             &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                             &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                             &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_EDGE],
                                             &dme->distrib_extract                 [PDM_MESH_ENTITY_VTX],
                                             &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VTX]);

  dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
  free(dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
  dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = NULL;
  if(dedge_vtx_idx == NULL)  {
    free(_dedge_vtx_idx);
  }

  dme->dmesh_extract->dn_edge = dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank];
  dme->dmesh_extract->dn_vtx  = dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VTX ][i_rank];

  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX] = PDM_block_to_part_create(distrib_vtx,
                                                             (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                    &dme->dmesh_extract->dn_vtx,
                                                                                    1,
                                                                                    dme->comm);

  _rebuild_group(dme, distrib_edge, PDM_MESH_ENTITY_EDGE, PDM_BOUND_TYPE_EDGE);
  _rebuild_group(dme, distrib_vtx , PDM_MESH_ENTITY_VTX,  PDM_BOUND_TYPE_VTX );

  free(distrib_edge);
  free(distrib_vtx );

}

static
void
_dmesh_extract_0d
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  int dn_vtx  = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_VTX);

  PDM_g_num_t* distrib_vtx  = PDM_compute_entity_distribution(dme->comm, dn_vtx );

  double *weight = malloc(dme->n_selected * sizeof(double));
  for(int i = 0; i < dme->n_selected; ++i) {
    weight[i] = 1.;
  }
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      &dme->selected_gnum,
                                                      &weight,
                                                      &dme->n_selected,
                                                      1,
                                                      dme->comm);
  free(weight);

  int          dn_extract_vtx      = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t *dextract_gnum_vtx   = PDM_part_to_block_block_gnum_get   (ptb);

  dme->distrib_extract[PDM_MESH_ENTITY_VTX] = PDM_compute_entity_distribution(dme->comm, dn_extract_vtx);

  dme->dmesh_extract->dn_vtx = dme->distrib_extract[PDM_MESH_ENTITY_VTX][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VTX][i_rank];

  dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX] = malloc(dme->dmesh_extract->dn_vtx * sizeof(PDM_g_num_t));
  for(int i = 0; i < dme->dmesh_extract->dn_vtx; ++i) {
    dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX][i] = dextract_gnum_vtx[i];
  }

  PDM_part_to_block_free(ptb);

  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX] = PDM_block_to_part_create(distrib_vtx,
                                                             (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                    &dme->dmesh_extract->dn_vtx,
                                                                                    1,
                                                                                    dme->comm);
  _rebuild_group(dme, distrib_vtx , PDM_MESH_ENTITY_VTX, PDM_BOUND_TYPE_VTX );

  free(distrib_vtx );
}


static
void
_dmesh_extract_nodal
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);
  PDM_MPI_Comm_size(dme->comm, &n_rank);

  /* Equilibrate */
  PDM_dmesh_nodal_elmts_t *dmn_elts = NULL;
  PDM_mesh_entities_t extracted_entity;
  if(dme->dim == 3) {
    dmn_elts = dme->dmesh_nodal->volumic;
    extracted_entity = PDM_MESH_ENTITY_CELL;
  } else if(dme->dim == 2) {
    dmn_elts = dme->dmesh_nodal->surfacic;
    extracted_entity = PDM_MESH_ENTITY_FACE;
  } else if(dme->dim == 1) {
    dmn_elts = dme->dmesh_nodal->ridge;
    extracted_entity = PDM_MESH_ENTITY_EDGE;
  } else {
    dmn_elts = dme->dmesh_nodal->corner;
    extracted_entity = PDM_MESH_ENTITY_VTX;
  }

  // For nodal meshes, parent_extract_gnum is filled for vertices *and* for the selected dim to extract :
  // parent gnum is ordered following *output* sections ordering in the extracted dmesh nodal
  PDM_dmesh_nodal_elmts_t* dmn_elts_extract = PDM_dmesh_nodal_elmts_to_extract_dmesh_nodal_elmts(dmn_elts,
                                                                                                 dme->n_selected,
                                                                                                 dme->selected_gnum,
                                                                                                 &dme->distrib_extract    [PDM_MESH_ENTITY_VTX],
                                                                                                 &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                                 &dme->distrib_extract    [extracted_entity],
                                                                                                 &dme->parent_extract_gnum[extracted_entity]);
  int dn_vtx  = dme->distrib_extract[PDM_MESH_ENTITY_VTX][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VTX][i_rank];

  const PDM_g_num_t *distrib_vtx = PDM_DMesh_nodal_distrib_vtx_get(dme->dmesh_nodal);

  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX] = PDM_block_to_part_create(distrib_vtx,
                                                             (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VTX],
                                                                                    &dn_vtx,
                                                                                    1,
                                                                                    dme->comm);

  PDM_dmesh_nodal_elmts_generate_distribution(dmn_elts_extract);

  double** tmp_dvtx_coord   = NULL;
  int stride_one = 1;
  PDM_block_to_part_exch(dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX],
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
             (void *  )  PDM_DMesh_nodal_vtx_get(dme->dmesh_nodal),
                         NULL,
             (void ***)  &tmp_dvtx_coord);

  PDM_g_num_t n_g_cell = 0;
  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_edge = 0;
  if(dme->dim ==  3) {
    n_g_cell = dmn_elts_extract->section_distribution[dmn_elts_extract->n_section];
  } else if(dme->dim == 2) {
    n_g_face = dmn_elts_extract->section_distribution[dmn_elts_extract->n_section];
  } else if(dme->dim == 1) {
    n_g_edge = dmn_elts_extract->section_distribution[dmn_elts_extract->n_section];
  }

  dme->dmesh_nodal_extract = PDM_DMesh_nodal_create(dme->comm,
                                                    dme->dim,
                                                    dme->distrib_extract[PDM_MESH_ENTITY_VTX][n_rank],
                                                    n_g_cell,
                                                    n_g_face,
                                                    n_g_edge);

  if(dme->dim ==  3) {
    PDM_DMesh_nodal_elmts_free(dme->dmesh_nodal_extract->volumic);
    dme->dmesh_nodal_extract->volumic = dmn_elts_extract;
  } else if(dme->dim == 2) {
    PDM_DMesh_nodal_elmts_free(dme->dmesh_nodal_extract->surfacic);
    dme->dmesh_nodal_extract->surfacic = dmn_elts_extract;
  } else if(dme->dim == 1) {
    PDM_DMesh_nodal_elmts_free(dme->dmesh_nodal_extract->ridge);
    dme->dmesh_nodal_extract->ridge    = dmn_elts_extract;
  }

  PDM_DMesh_nodal_coord_set(dme->dmesh_nodal_extract,
                            dn_vtx,
                            tmp_dvtx_coord[0],
                            PDM_OWNERSHIP_KEEP);

  free(tmp_dvtx_coord);

}



/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_extract_t*
PDM_dmesh_extract_create
(
 const int                     dim,
       PDM_MPI_Comm            comm
)
{
  PDM_dmesh_extract_t *dme = (PDM_dmesh_extract_t *) malloc(sizeof(PDM_dmesh_extract_t));

  dme->dim  = dim;
  dme->comm = comm;

  // Utilisation privé
  dme->dmesh                   = PDM_dmesh_create(PDM_OWNERSHIP_KEEP, 0, 0, 0, 0, comm);
  dme->dmesh_nodal             = NULL;
  dme->dmesh_extract           = NULL;
  dme->dmesh_nodal_extract     = NULL;
  dme->dmesh_extract_ownership = PDM_OWNERSHIP_KEEP;
  dme->dmesh_shared            = NULL;

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    dme->btp_entity_to_extract_entity [i] = NULL;
    dme->btp_ownership                [i] = PDM_OWNERSHIP_KEEP;

    dme->distrib_extract              [i] = NULL;
    dme->parent_extract_gnum          [i] = NULL;

    dme->distrib_extract_ownership    [i] = PDM_OWNERSHIP_KEEP;
    dme->parent_extract_gnum_ownership[i] = PDM_OWNERSHIP_KEEP;

  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    dme->n_bound                           [i] = 0;
    dme->btp_bound_entity_to_extract_entity[i] = NULL;
    dme->btp_bound_ownership               [i] = NULL;

  }

  return dme;
}

void
PDM_dmesh_extract_compute
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);
  if(dme->dmesh_shared == NULL){
    // Synchronize dmesh
    PDM_g_num_t _dn_cell = dme->dmesh->dn_cell;
    PDM_g_num_t _dn_face = dme->dmesh->dn_face;
    PDM_g_num_t _dn_edge = dme->dmesh->dn_edge;
    PDM_g_num_t _dn_vtx  = dme->dmesh->dn_vtx;

    PDM_MPI_Allreduce(&_dn_cell, &dme->dmesh->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_face, &dme->dmesh->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_edge, &dme->dmesh->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_vtx , &dme->dmesh->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);

  } else if (dme->dmesh_nodal == NULL) {
    PDM_dmesh_free(dme->dmesh);
    dme->dmesh = dme->dmesh_shared;
  }

  if(dme->dmesh_nodal == NULL) {
    dme->dmesh_extract = PDM_dmesh_create(PDM_OWNERSHIP_KEEP, 0, 0, 0, 0, dme->comm);
    dme->dmesh_extract->_dvtx_coord = NULL;

    if(dme->dim == 3) {
      _dmesh_extract_3d(dme);
    } else if(dme->dim == 2) {
      _dmesh_extract_2d(dme);
    } else if(dme->dim == 1) {
      _dmesh_extract_1d(dme);
    } else {
      _dmesh_extract_0d(dme);
    }

    /* Exchange coordinates */
    if(dme->dmesh->_dvtx_coord != NULL) {
      double** tmp_dvtx_coord   = NULL;
      int stride_one = 1;
      PDM_block_to_part_exch(dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VTX],
                             3 * sizeof(double),
                             PDM_STRIDE_CST_INTERLACED,
                             &stride_one,
                 (void *  )  dme->dmesh->_dvtx_coord,
                             NULL,
                 (void ***)  &tmp_dvtx_coord);

      PDM_dmesh_vtx_coord_set(dme->dmesh_extract,
                              tmp_dvtx_coord[0],
                              PDM_OWNERSHIP_KEEP);

      free(tmp_dvtx_coord);
    }

    PDM_g_num_t _dn_extract_cell = dme->dmesh_extract->dn_cell;
    PDM_g_num_t _dn_extract_face = dme->dmesh_extract->dn_face;
    PDM_g_num_t _dn_extract_edge = dme->dmesh_extract->dn_edge;
    PDM_g_num_t _dn_extract_vtx  = dme->dmesh_extract->dn_vtx;

    PDM_MPI_Allreduce(&_dn_extract_cell, &dme->dmesh_extract->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_extract_face, &dme->dmesh_extract->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_extract_edge, &dme->dmesh_extract->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
    PDM_MPI_Allreduce(&_dn_extract_vtx , &dme->dmesh_extract->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  } else {
    assert(dme->dmesh_nodal != NULL);
    _dmesh_extract_nodal(dme);
  }
}


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
)
{
  PDM_UNUSED(entity_type);
  dme->n_selected    = n_selected;
  dme->selected_gnum = selected_gnum;

}


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
)
{
  PDM_dmesh_dn_entity_set(dme->dmesh, entity_type, dn_entity);
}


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
)
{
  PDM_dmesh_vtx_coord_set(dme->dmesh, dvtx_coord, PDM_OWNERSHIP_USER);
}

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
)
{
  PDM_dmesh_bound_set(dme->dmesh, bound_type, n_bound, connect, connect_idx, PDM_OWNERSHIP_USER);
}


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
)
{
  PDM_dmesh_connectivity_set(dme->dmesh,
                             connectivity_type,
                             dconnect,
                             dconnect_idx,
                             PDM_OWNERSHIP_USER);
}

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
)
{
  dme->dmesh_shared = dmesh;
}


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
)
{
  dme->dmesh_nodal = dmesh_nodal;
}

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
)
{
  *dmesh_extract = dme->dmesh_extract;
  dme->dmesh_extract_ownership = ownership;
}

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
)
{
  *dmesh_nodal_extract = dme->dmesh_nodal_extract;
  dme->dmesh_extract_ownership = ownership;
}

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
)
{
  if(dme->distrib_extract[entity_type] != NULL){
    int i_rank;
    PDM_MPI_Comm_rank(dme->comm, &i_rank);
    *dn_entity   = dme->distrib_extract[entity_type][i_rank+1] - dme->distrib_extract[entity_type][i_rank];
    *parent_gnum = dme->parent_extract_gnum[entity_type];
    dme->parent_extract_gnum_ownership[entity_type] = ownership;
  } else {
    *dn_entity = 0;
    *parent_gnum = NULL;
  }
}


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
)
{
  *btp = dme->btp_entity_to_extract_entity[entity_type];
  dme->btp_ownership[entity_type] = ownership;
}


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
)
{
  *btp = dme->btp_bound_entity_to_extract_entity[bound_type][i_group];
  dme->btp_bound_ownership[bound_type][i_group] = ownership;
}


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
)
{

  if(dme->dmesh_shared == NULL){
    PDM_dmesh_free(dme->dmesh);
  }

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if (dme->btp_ownership[i] == PDM_OWNERSHIP_KEEP && dme->btp_entity_to_extract_entity[i] != NULL) {
      PDM_block_to_part_free(dme->btp_entity_to_extract_entity[i]);
    }
  }


  for (int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(dme->distrib_extract_ownership[i] == PDM_OWNERSHIP_KEEP && dme->distrib_extract[i] != NULL) {
      free(dme->distrib_extract[i]);
    }

    if(dme->parent_extract_gnum_ownership[i] == PDM_OWNERSHIP_KEEP && dme->parent_extract_gnum[i] != NULL) {
      free(dme->parent_extract_gnum[i]);
    }
  }


  for (int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {

    for(int i_group = 0; i_group < dme->n_bound[i]; ++i_group) {
      if (dme->btp_bound_ownership[i][i_group] == PDM_OWNERSHIP_KEEP && dme->btp_bound_entity_to_extract_entity[i][i_group] != NULL) {
        PDM_block_to_part_free(dme->btp_bound_entity_to_extract_entity[i][i_group]);
      }
    }

    if(dme->btp_bound_entity_to_extract_entity[i] != NULL) {
      free(dme->btp_bound_entity_to_extract_entity[i]);
      free(dme->btp_bound_ownership               [i]);
    }


  }

  if(dme->dmesh_extract_ownership == PDM_OWNERSHIP_KEEP && dme->dmesh_extract != NULL) {
    PDM_dmesh_free(dme->dmesh_extract);
  }


  if(dme->dmesh_extract_ownership == PDM_OWNERSHIP_KEEP && dme->dmesh_nodal_extract != NULL) {
    PDM_DMesh_nodal_free(dme->dmesh_nodal_extract);
  }

  free(dme);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
