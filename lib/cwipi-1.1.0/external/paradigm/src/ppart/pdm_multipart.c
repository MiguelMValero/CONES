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

/*============================================================================
 * TODO : write module description here
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
#include "pdm_timer.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_mesh_nodal.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_unique.h"
#include "pdm_partitioning_nodal_algorithm.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_part_connectivity_transform.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
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
PDM_bound_type_t
_entity_type_to_bound_type
(
 PDM_mesh_entities_t entity_type
 )
{
  PDM_bound_type_t bound_type;
  switch (entity_type) {
    case PDM_MESH_ENTITY_VTX:
      bound_type = PDM_BOUND_TYPE_VTX;
      break;

    case PDM_MESH_ENTITY_EDGE:
      bound_type = PDM_BOUND_TYPE_EDGE;
      break;

    case PDM_MESH_ENTITY_FACE:
      bound_type = PDM_BOUND_TYPE_FACE;
      break;

    default:
      PDM_error(__FILE__, __LINE__, 0, "Entity type %d has no corresponding bound type\n", entity_type);
  }

  return bound_type;
}


/**
 *
 * \brief Translate in _part_t structure (only mapping)
 */
static
_part_t**
_map_part_t_with_part_mesh
(
  PDM_part_mesh_t* pm
)
{
  int n_part = pm->n_part;
  _part_t **pdm_part = (_part_t **) malloc(n_part * sizeof(_part_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    pdm_part[i_part] = _part_create();

    pdm_part[i_part]->n_cell                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_CELL);
    pdm_part[i_part]->n_face                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_FACE);
    pdm_part[i_part]->n_edge                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);
    pdm_part[i_part]->n_vtx                    = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VTX);
    pdm_part[i_part]->n_section                = 0;
    pdm_part[i_part]->n_elt                    = NULL;

    pdm_part[i_part]->n_face_group             = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_FACE);
    pdm_part[i_part]->n_edge_group             = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_EDGE);

    pdm_part[i_part]->n_face_part_bound        = 0;
    pdm_part[i_part]->n_vtx_part_bound         = 0;

    PDM_part_mesh_vtx_coord_get(pm, i_part, &pdm_part[i_part]->vtx, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   &pdm_part[i_part]->face_vtx,
                                   &pdm_part[i_part]->face_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   &pdm_part[i_part]->cell_face,
                                   &pdm_part[i_part]->cell_face_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &pdm_part[i_part]->face_edge,
                                   &pdm_part[i_part]->face_edge_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                   &pdm_part[i_part]->edge_face,
                                   &pdm_part[i_part]->edge_face_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    int* face_cell_idx = NULL;
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                   &pdm_part[i_part]->face_cell,
                                   &face_cell_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    int* edge_vtx_idx = NULL;
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &pdm_part[i_part]->edge_vtx,
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &pdm_part[i_part]->cell_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pdm_part[i_part]->face_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &pdm_part[i_part]->edge_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &pdm_part[i_part]->vtx_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_bound_concat_get(pm,
                                   i_part,
                                   PDM_BOUND_TYPE_FACE,
                                   &pdm_part[i_part]->face_bound_idx,
                                   &pdm_part[i_part]->face_bound,
                                   &pdm_part[i_part]->face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_bound_concat_get(pm,
                                   i_part,
                                   PDM_BOUND_TYPE_EDGE,
                                   &pdm_part[i_part]->edge_bound_idx,
                                   &pdm_part[i_part]->edge_bound,
                                   &pdm_part[i_part]->edge_bound_ln_to_gn,
                                   PDM_OWNERSHIP_KEEP);

    pdm_part[i_part]->face_group_idx           = NULL;
    pdm_part[i_part]->face_group               = NULL;
    pdm_part[i_part]->face_group_ln_to_gn      = NULL;

    pdm_part[i_part]->face_part_bound_proc_idx = NULL;
    pdm_part[i_part]->face_part_bound_part_idx = NULL;
    pdm_part[i_part]->face_part_bound          = NULL;

    pdm_part[i_part]->edge_part_bound_proc_idx = NULL;
    pdm_part[i_part]->edge_part_bound_part_idx = NULL;
    pdm_part[i_part]->edge_part_bound          = NULL;

    pdm_part[i_part]->vtx_part_bound_proc_idx  = NULL;
    pdm_part[i_part]->vtx_part_bound_part_idx  = NULL;
    pdm_part[i_part]->vtx_part_bound           = NULL;

    pdm_part[i_part]->face_join_ln_to_gn       = NULL;
    pdm_part[i_part]->face_join_idx            = NULL;
    pdm_part[i_part]->face_join                = NULL;
    pdm_part[i_part]->edge_join_idx            = NULL;
    pdm_part[i_part]->edge_join                = NULL;

    pdm_part[i_part]->elt_section_ln_to_gn     = NULL;

    pdm_part[i_part]->cell_tag                 = NULL;
    pdm_part[i_part]->face_tag                 = NULL;
    pdm_part[i_part]->edge_tag                 = NULL;
    pdm_part[i_part]->vtx_tag                  = NULL;

    pdm_part[i_part]->cell_weight              = NULL;
    pdm_part[i_part]->face_weight              = NULL;

    pdm_part[i_part]->cell_color               = NULL;
    pdm_part[i_part]->face_color               = NULL;
    pdm_part[i_part]->face_hp_color            = NULL;
    pdm_part[i_part]->edge_color               = NULL;
    pdm_part[i_part]->vtx_color                = NULL;
    pdm_part[i_part]->thread_color             = NULL;
    pdm_part[i_part]->hyperplane_color         = NULL;

    pdm_part[i_part]->vtx_ghost_information    = NULL;

    pdm_part[i_part]->new_to_old_order_cell    = NULL;
    pdm_part[i_part]->new_to_old_order_face    = NULL;
    pdm_part[i_part]->new_to_old_order_edge    = NULL;
    pdm_part[i_part]->new_to_old_order_vtx     = NULL;

    pdm_part[i_part]->subpartlayout            = NULL;
  }

  return pdm_part;
}

/**
 *
 * \brief Free the memory occuped by a partition structure
 *
 * \param [inout]   part          _part_t object
 */
static void
_part_free
(
 _part_t         *part
)
{
  /* Following is not results but internal array */
  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

  if (part->new_to_old_order_edge != NULL)
    free(part->new_to_old_order_edge);
  part->new_to_old_order_edge = NULL;


  if (part->new_to_old_order_vtx != NULL)
    free(part->new_to_old_order_vtx);
  part->new_to_old_order_vtx = NULL;

  if(part->subpartlayout != NULL){
    if(part->subpartlayout->cell_tile_idx!= NULL)
      free(part->subpartlayout->cell_tile_idx);
    if(part->subpartlayout->face_tile_idx!= NULL)
      free(part->subpartlayout->face_tile_idx);
    if(part->subpartlayout->face_bnd_tile_idx!= NULL)
      free(part->subpartlayout->face_bnd_tile_idx);
    if(part->subpartlayout->mask_tile_idx!= NULL)
      free(part->subpartlayout->mask_tile_idx);
    if(part->subpartlayout->cell_vect_tile_idx!= NULL)
      free(part->subpartlayout->cell_vect_tile_idx);
    if(part->subpartlayout->mask_tile_n!= NULL)
      free(part->subpartlayout->mask_tile_n);
    if(part->subpartlayout->cell_vect_tile_n!= NULL)
      free(part->subpartlayout->cell_vect_tile_n);
    if(part->subpartlayout->mask_tile!= NULL)
      free(part->subpartlayout->mask_tile);
    free(part->subpartlayout);
  }


  free(part);
}


static void
_setup_ghost_information
(
const int           i_rank,
const int           n_part,
const int          *pn_vtx,
      int         **pinternal_vtx_priority,
      PDM_g_num_t  *distrib_partition
)
{
  /* 0 : Interior / 1 : owner join (at least one) / 2 : not owner */
  // pinternal_vtx_priority contains value between i_rank + i_part
  for (int ipart = 0; ipart < n_part; ipart++) {
    for(int ivtx = 0; ivtx < pn_vtx[ipart]; ++ivtx) {

      int g_part = pinternal_vtx_priority[ipart][ivtx];
      int t_part = -1;
      if( g_part >= distrib_partition[i_rank] && g_part < distrib_partition[i_rank+1]) {
        t_part = g_part - distrib_partition[i_rank];
      }

      // if(pinternal_vtx_priority[ipart][ivtx] == i_rank){
      if(t_part == ipart){
        pinternal_vtx_priority[ipart][ivtx] = 1;
      } else if(pinternal_vtx_priority[ipart][ivtx] == -1) {
        pinternal_vtx_priority[ipart][ivtx] = 0;
      } else { /* Pas owner / pas interieur donc sur 1 autre proc */
         pinternal_vtx_priority[ipart][ivtx] = 2;
      }
    }
  }
}

static
void
_create_dparent_num_corner
(
  PDM_dmesh_nodal_elmts_t  *dmn_elts,
  PDM_g_num_t             **dparent_gnum_corner,
  PDM_g_num_t             **distrib_corner
)
{
  // Just repart all section in one block
  int n_rank = -1;
  PDM_MPI_Comm_size (dmn_elts->comm, &n_rank);

  int n_section = dmn_elts->n_section;

  int          *dn_corner       = malloc(n_section * sizeof(int          ));
  PDM_g_num_t **corner_ln_to_gn = malloc(n_section * sizeof(PDM_g_num_t *));
  PDM_g_num_t **corner_vtx      = malloc(n_section * sizeof(PDM_g_num_t *));
  int         **corner_vtx_n    = malloc(n_section * sizeof(int         *));

  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = dmn_elts->sections_id[i_section];
    const PDM_g_num_t* distrib = PDM_DMesh_nodal_elmts_distrib_section_get(dmn_elts, id_section);

    PDM_g_num_t beg_elmt_gnum = distrib[dmn_elts->i_rank] + dmn_elts->section_distribution[i_section];

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmn_elts, id_section);

    assert(t_elt == PDM_MESH_NODAL_POINT);

    int n_elt           = PDM_DMesh_nodal_elmts_section_n_elt_get(dmn_elts, id_section);
    PDM_g_num_t* connec = PDM_DMesh_nodal_elmts_section_std_get(dmn_elts, id_section);

    dn_corner[i_section] = n_elt;
    corner_ln_to_gn[i_section] = malloc( n_elt * sizeof(PDM_g_num_t));
    corner_vtx_n   [i_section] = malloc( n_elt * sizeof(int        ));
    corner_vtx     [i_section] = connec;
    for(int i = 0; i < n_elt; ++i) {
      corner_vtx_n   [i_section][i] = 1;
      corner_ln_to_gn[i_section][i] = beg_elmt_gnum + i + 1;
    }

    // PDM_log_trace_array_long(connec, n_elt, "corner_vtx ::");
    // PDM_log_trace_array_long(corner_ln_to_gn[i_section], n_elt, "corner_ln_to_gn ::");

  }

  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      corner_ln_to_gn,
                                                      NULL,
                                                      dn_corner,
                                                      n_section,
                                                      dmn_elts->comm);


  PDM_g_num_t* distrib_ptb = PDM_part_to_block_distrib_index_get(ptb);

  PDM_g_num_t* _distrib_corner = malloc((n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    _distrib_corner[i] = distrib_ptb[i];
  }

  int         *blk_child_n    = NULL;
  PDM_g_num_t *blk_child_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         corner_vtx_n,
             (void **)   corner_vtx,
                        &blk_child_n,
             (void **)  &blk_child_gnum);

  *dparent_gnum_corner = blk_child_gnum;
  free(blk_child_n);

  PDM_part_to_block_free(ptb);
  for (int i_section = 0; i_section < n_section; i_section++) {
    free(corner_ln_to_gn[i_section]);
    free(corner_vtx_n[i_section]);
  }
  free(corner_ln_to_gn);
  free(corner_vtx_n);
  free(corner_vtx);
  free(dn_corner);

  *distrib_corner = _distrib_corner;
}

static PDM_part_mesh_nodal_t*
_compute_part_mesh_nodal_3d
(
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part,
 PDM_ownership_t    ownership
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  PDM_UNUSED(ownership);

  /*
   * Rebuild the volumic part from cell
   */
  int           *pn_cell        = (int          * ) malloc( n_part * sizeof(int          ));
  int           *pn_face        = (int          * ) malloc( n_part * sizeof(int          ));
  int           *pn_vtx         = (int          * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t  **pcell_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  double       **pvtx_coord     = (double      ** ) malloc( n_part * sizeof(double      *));
  for(int i_part = 0; i_part < n_part; ++i_part){

    pn_cell[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_CELL);
    pn_face[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_FACE);
    pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &pcell_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pface_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &pvtx_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_vtx_coord_get(pm->pmesh,
                                i_part,
                                &pvtx_coord[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_vol = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->volumic,
                                                                                        n_part,
                                                                                        pn_vtx,
                                                                                        pvtx_ln_to_gn,
                                                                                        pn_cell,
                                                                                        NULL,
                                                                                        pcell_ln_to_gn,
                                                                                        NULL);

  int          *pn_surf             = NULL;
  int         **psurf_to_entity     = NULL;
  PDM_g_num_t **psurf_gnum          = NULL;
  PDM_g_num_t **psurf_to_face_g_num = NULL;
  PDM_reverse_dparent_gnum(dmn->surfacic->dparent_gnum,
                           NULL, // dparent_sign
                           dmn->surfacic->delmt_child_distrib,
                           n_part,
                           pn_face,
                           pface_ln_to_gn,
                          &pn_surf,
                          &psurf_to_entity,
                          &psurf_gnum,
                          &psurf_to_face_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  if(0 == 1) {
    PDM_log_trace_array_long(dmn->surfacic->dparent_gnum, dmn->surfacic->delmt_child_distrib[i_rank+1] - dmn->surfacic->delmt_child_distrib[i_rank], "dmn->surfacic->dparent_gnum : ");
    for(int i_part = 0; i_part < n_part; ++i_part){
      PDM_log_trace_array_long(pface_ln_to_gn     [i_part], pn_face[i_part], "pface_ln_to_gn      : ");
      PDM_log_trace_array_long(psurf_to_face_g_num[i_part], pn_surf[i_part], "psurf_to_face_g_num : ");
    }
  }

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_surf,
                                                                                         psurf_to_entity,
                                                                                         psurf_gnum,
                                                                                         psurf_to_face_g_num);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(psurf_gnum         [i_part]);
    free(psurf_to_face_g_num[i_part]);
    free(psurf_to_entity    [i_part]);
  }
  free(psurf_to_entity);
  free(pn_surf);
  free(psurf_gnum);
  free(psurf_to_face_g_num);
  free(pface_ln_to_gn);
  free(pn_face);

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = NULL;

  if(dmn->ridge != NULL) {

    PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
    int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_edge[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_EDGE,
                                        &pedge_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);
    }
    int          *pn_ridge             = NULL;
    int         **pridge_to_entity     = NULL;
    PDM_g_num_t **pridge_gnum          = NULL;
    PDM_g_num_t **pridge_to_edge_g_num = NULL;
    PDM_reverse_dparent_gnum(dmn->ridge->dparent_gnum,
                             NULL, // dparent_sign
                             dmn->ridge->delmt_child_distrib,
                             n_part,
                             pn_edge,
                             pedge_ln_to_gn,
                            &pn_ridge,
                            &pridge_to_entity,
                            &pridge_gnum,
                            &pridge_to_edge_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                               n_part,
                                                               pn_vtx,
                                                               pvtx_ln_to_gn,
                                                               pn_ridge,
                                                               pridge_to_entity,
                                                               pridge_gnum,
                                                               pridge_to_edge_g_num);

    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pridge_gnum[i_part]);
      free(pridge_to_edge_g_num[i_part]);
      free(pridge_to_entity[i_part]);
    }
    free(pridge_to_entity);
    free(pn_ridge);
    free(pridge_gnum);
    free(pridge_to_edge_g_num);
    free(pedge_ln_to_gn);
    free(pn_edge);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_corner = NULL;

  if(dmn->corner != NULL) {

    PDM_g_num_t *dparent_gnum_corner = NULL;
    PDM_g_num_t *distrib_corner      = NULL;
    _create_dparent_num_corner(dmn->corner, &dparent_gnum_corner, &distrib_corner);

    int          *pn_corner             = NULL;
    int         **pcorner_to_entity     = NULL;
    PDM_g_num_t **pcorner_gnum          = NULL;
    PDM_g_num_t **pcorner_to_vtx_g_num  = NULL;
    PDM_reverse_dparent_gnum(dparent_gnum_corner,
                             NULL, // dparent_sign
                             distrib_corner,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                            &pn_corner,
                            &pcorner_to_entity,
                            &pcorner_gnum,
                            &pcorner_to_vtx_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_corner = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->corner,
                                                                n_part,
                                                                pn_vtx,
                                                                pvtx_ln_to_gn,
                                                                pn_corner,
                                                                pcorner_to_entity,
                                                                pcorner_gnum,
                                                                pcorner_to_vtx_g_num);
    free(dparent_gnum_corner);
    free(distrib_corner);
    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pcorner_gnum[i_part]);
      free(pcorner_to_vtx_g_num[i_part]);
      free(pcorner_to_entity[i_part]);
    }
    free(pcorner_to_entity);
    free(pn_corner);
    free(pcorner_gnum);
    free(pcorner_to_vtx_g_num);
  }

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_vol);
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf);
  if(pmn_ridge != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge);
  }
  if(pmn_corner != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_corner);
  }
  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Copy coordinates because ownership between part_mesh and part_mesh_nodal is complicated
    double      *lvtx_coords   = (double      *) malloc(3 * pn_vtx[i_part] * sizeof(double     ));
    PDM_g_num_t *lvtx_ln_to_gn = (PDM_g_num_t *) malloc(3 * pn_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < 3 * pn_vtx[i_part]; ++i_vtx) {
      lvtx_coords[i_vtx] = pvtx_coord[i_part][i_vtx];
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      lvtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_part][i_vtx];
    }
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  pn_vtx[i_part],
                                  lvtx_coords,
                                  lvtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  }

  free(pcell_ln_to_gn);
  free(pn_cell);
  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_vtx);

  return pmn;
}

static PDM_part_mesh_nodal_t*
_compute_part_mesh_nodal_2d
(
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part,
 PDM_ownership_t    ownership
)
{
  PDM_UNUSED(ownership);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_face        = (int *  )         malloc( n_part * sizeof(int          ));

  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_vtx         = (int *  )         malloc( n_part * sizeof(int          ));
  double       **pvtx_coord     = (double      ** ) malloc( n_part * sizeof(double      *));
  for(int i_part = 0; i_part < n_part; ++i_part){
    pn_face[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_FACE);
    pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pface_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &pvtx_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_vtx_coord_get(pm->pmesh,
                                i_part,
                                &pvtx_coord[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_face,
                                                                                         NULL,
                                                                                         pface_ln_to_gn,
                                                                                         NULL);

  PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_edge[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_EDGE);
    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &pedge_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }
  int          *pn_ridge;
  int         **pridge_to_entity;
  PDM_g_num_t **pridge_gnum;
  PDM_g_num_t **pridge_to_edge_g_num;
  PDM_reverse_dparent_gnum(dmn->ridge->dparent_gnum,
                           NULL, // dparent_sign
                           dmn->ridge->delmt_child_distrib,
                           n_part,
                           pn_edge,
                           pedge_ln_to_gn,
                          &pn_ridge,
                          &pridge_to_entity,
                          &pridge_gnum,
                          &pridge_to_edge_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                                                          n_part,
                                                                                          pn_vtx,
                                                                                          pvtx_ln_to_gn,
                                                                                          pn_ridge,
                                                                                          pridge_to_entity,
                                                                                          pridge_gnum,
                                                                                          pridge_to_edge_g_num);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(pridge_gnum[i_part]);
    free(pridge_to_edge_g_num[i_part]);
    free(pridge_to_entity[i_part]);
  }
  free(pridge_to_entity);
  free(pn_ridge);
  free(pridge_gnum);
  free(pridge_to_edge_g_num);

  PDM_part_mesh_nodal_elmts_t* pmn_corner = NULL;

  if(dmn->corner != NULL) {

    PDM_g_num_t *dparent_gnum_corner = NULL;
    PDM_g_num_t *distrib_corner      = NULL;
    _create_dparent_num_corner(dmn->corner, &dparent_gnum_corner, &distrib_corner);

    int          *pn_corner             = NULL;
    int         **pcorner_to_entity     = NULL;
    PDM_g_num_t **pcorner_gnum          = NULL;
    PDM_g_num_t **pcorner_to_vtx_g_num  = NULL;
    PDM_reverse_dparent_gnum(dparent_gnum_corner,
                             NULL, // dparent_sign
                             distrib_corner,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                            &pn_corner,
                            &pcorner_to_entity,
                            &pcorner_gnum,
                            &pcorner_to_vtx_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_corner = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->corner,
                                                                n_part,
                                                                pn_vtx,
                                                                pvtx_ln_to_gn,
                                                                pn_corner,
                                                                pcorner_to_entity,
                                                                pcorner_gnum,
                                                                pcorner_to_vtx_g_num);
    free(dparent_gnum_corner);
    free(distrib_corner);
    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pcorner_gnum[i_part]);
      free(pcorner_to_vtx_g_num[i_part]);
      free(pcorner_to_entity[i_part]);
    }
    free(pcorner_to_entity);
    free(pn_corner);
    free(pcorner_gnum);
    free(pcorner_to_vtx_g_num);
  }

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf);
  if(pmn_ridge != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge);
  }
  if(pmn_corner != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_corner);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    // Copy coordinates because ownership between part_mesh and part_mesh_nodal is complicated
    double      *lvtx_coords   = (double      *) malloc(3 * pn_vtx[i_part] * sizeof(double     ));
    PDM_g_num_t *lvtx_ln_to_gn = (PDM_g_num_t *) malloc(    pn_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < 3 * pn_vtx[i_part]; ++i_vtx) {
      lvtx_coords[i_vtx] = pvtx_coord[i_part][i_vtx];
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      lvtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_part][i_vtx];
    }
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  pn_vtx[i_part],
                                  lvtx_coords,
                                  lvtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
    // free(lvtx_ln_to_gn);
  }

  free(pedge_ln_to_gn);
  free(pn_edge);
  free(pface_ln_to_gn);
  free(pn_face);

  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_vtx);
  return pmn;
}

static PDM_part_mesh_nodal_t*
_compute_part_mesh_nodal_1d
(
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part,
 PDM_ownership_t    ownership
)
{
  PDM_UNUSED(ownership);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));

  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_vtx         = (int *  )         malloc( n_part * sizeof(int          ));
  double       **pvtx_coord     = (double      ** ) malloc( n_part * sizeof(double      *));
  for(int i_part = 0; i_part < n_part; ++i_part){
    pn_edge[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_EDGE);
    pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &pedge_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &pvtx_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_vtx_coord_get(pm->pmesh,
                                i_part,
                                &pvtx_coord[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                                                          n_part,
                                                                                          pn_vtx,
                                                                                          pvtx_ln_to_gn,
                                                                                          pn_edge,
                                                                                          NULL,
                                                                                          pedge_ln_to_gn,
                                                                                          NULL);


  PDM_part_mesh_nodal_elmts_t* pmn_corner = NULL;

  if(dmn->corner != NULL) {

    PDM_g_num_t *dparent_gnum_corner = NULL;
    PDM_g_num_t *distrib_corner      = NULL;
    _create_dparent_num_corner(dmn->corner, &dparent_gnum_corner, &distrib_corner);

    int          *pn_corner             = NULL;
    int         **pcorner_to_entity     = NULL;
    PDM_g_num_t **pcorner_gnum          = NULL;
    PDM_g_num_t **pcorner_to_vtx_g_num  = NULL;
    PDM_reverse_dparent_gnum(dparent_gnum_corner,
                             NULL, // dparent_sign
                             distrib_corner,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                            &pn_corner,
                            &pcorner_to_entity,
                            &pcorner_gnum,
                            &pcorner_to_vtx_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_corner = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->corner,
                                                                n_part,
                                                                pn_vtx,
                                                                pvtx_ln_to_gn,
                                                                pn_corner,
                                                                pcorner_to_entity,
                                                                pcorner_gnum,
                                                                pcorner_to_vtx_g_num);
    free(dparent_gnum_corner);
    free(distrib_corner);
    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pcorner_gnum[i_part]);
      free(pcorner_to_vtx_g_num[i_part]);
      free(pcorner_to_entity[i_part]);
    }
    free(pcorner_to_entity);
    free(pn_corner);
    free(pcorner_gnum);
    free(pcorner_to_vtx_g_num);
  }

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge);
  if(pmn_corner != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_corner);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    // Copy coordinates because ownership between part_mesh and part_mesh_nodal is complicated
    double      *lvtx_coords   = (double      *) malloc(3 * pn_vtx[i_part] * sizeof(double     ));
    PDM_g_num_t *lvtx_ln_to_gn = (PDM_g_num_t *) malloc(    pn_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < 3 * pn_vtx[i_part]; ++i_vtx) {
      lvtx_coords[i_vtx] = pvtx_coord[i_part][i_vtx];
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      lvtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_part][i_vtx];
    }
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  pn_vtx[i_part],
                                  lvtx_coords,
                                  lvtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  }

  free(pedge_ln_to_gn);
  free(pn_edge);

  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_vtx);
  return pmn;
}


static
void
_split_graph_hilbert
(
 PDM_MPI_Comm   comm,
 PDM_dmesh_t   *dmesh,
 int            n_part,
 int           *node_part
)
{
  if(dmesh->n_g_cell != 0) {

    int         *dcell_face_idx = NULL;
    PDM_g_num_t *dcell_face     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               &dcell_face,
                               &dcell_face_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    int         *dface_vtx_idx = NULL;
    PDM_g_num_t *dface_vtx     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               &dface_vtx,
                               &dface_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_g_num_t *distrib_face = NULL;
    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);

    int own_distrib_face = 0;
    if(distrib_face == NULL) {
      own_distrib_face = 1;
      distrib_face = PDM_compute_entity_distribution(comm, dmesh->dn_face);
    }

    PDM_g_num_t *distrib_vtx = NULL;
    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_VTX, &distrib_vtx);
    int own_distrib_vtx = 0;
    if(distrib_vtx == NULL) {
      own_distrib_vtx = 1;
      distrib_vtx = PDM_compute_entity_distribution(comm, dmesh->dn_vtx);
    }

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_geom (PDM_PART_GEOM_HILBERT,
                   n_part,
                   comm,
                   dmesh->dn_cell,
                   dcell_face_idx,
                   dcell_face,
                   NULL, //cell_weight
                   dface_vtx_idx,
                   dface_vtx,
                   distrib_face,
                   dvtx_coord,
                   distrib_vtx,
                   node_part);

    if(own_distrib_vtx) {
      free(distrib_vtx);
    }

    if(own_distrib_face) {
      free(distrib_face);
    }

  } else if (dmesh->n_g_face != 0) {

    int         *dface_vtx_idx = NULL;
    PDM_g_num_t *dface_vtx     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               &dface_vtx,
                               &dface_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);


    int         *dedge_vtx_idx = NULL;
    PDM_g_num_t *dedge_vtx     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);


    int         *dface_edge_idx = NULL;
    PDM_g_num_t *dface_edge     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                               &dface_edge,
                               &dface_edge_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_geom_2d(PDM_PART_GEOM_HILBERT,
                     n_part,
                     comm,
                     dmesh->dn_face,
                     dmesh->dn_edge,
                     dmesh->dn_vtx,
                     dface_vtx_idx,
                     dface_vtx,
                     dface_edge_idx,
                     dface_edge,
                     dedge_vtx,
                     dvtx_coord,
                     NULL,
                     node_part);

  } else if (dmesh->n_g_edge != 0) {

    int         *dedge_vtx_idx = NULL;
    PDM_g_num_t *dedge_vtx     = NULL;
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_geom_1d(PDM_PART_GEOM_HILBERT,
                     n_part,
                     comm,
                     dmesh->dn_edge,
                     dmesh->dn_vtx,
                     dedge_vtx,
                     dvtx_coord,
                     NULL,
                     node_part);

  } else if (dmesh->n_g_vtx != 0) {

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_geom_0d(PDM_PART_GEOM_HILBERT,
                     n_part,
                     comm,
                     dmesh->dn_vtx,
                     dvtx_coord,
                     NULL,
                     node_part);
  }
}

static
void
_warm_up_for_split
(
 PDM_MPI_Comm       comm,
 PDM_dmesh_t       *dmesh,
 PDM_g_num_t       *distrib_node,
 int                compute_dual,
 PDM_g_num_t      **out_dual_graph_idx,
 PDM_g_num_t      **out_dual_graph
)
{

  PDM_g_num_t *dual_graph_idx = NULL;
  PDM_g_num_t *dual_graph     = NULL;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int         *darc_to_elmt_idx = NULL; // Donc face_cell OU edge_face
  PDM_g_num_t *darc_to_elmt_tmp = NULL;
  PDM_g_num_t *darc_to_elmt     = NULL;
  int         *delmt_to_arc_idx = NULL; // Donc cell_face OU face_edge
  PDM_g_num_t *delmt_to_arc     = NULL;
  int dn_node = 0;
  int dn_arc  = 0;

  PDM_g_num_t *distrib_arc  = NULL;

  int is1d = 0;

  if(dmesh->n_g_cell != 0) { // Donc 3D
    dn_node = dmesh->dn_cell;
    dn_arc  = dmesh->dn_face;
    distrib_arc  = PDM_compute_entity_distribution(comm, dn_arc );

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               &darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               &delmt_to_arc,
                               &delmt_to_arc_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    if(darc_to_elmt_tmp == NULL) {
      assert(delmt_to_arc_idx != NULL);
      PDM_dcellface_to_dfacecell(distrib_arc,
                                 distrib_node,
                                 delmt_to_arc_idx,
                                 delmt_to_arc,
                                 &darc_to_elmt_tmp,
                                 comm);

      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                 darc_to_elmt_tmp,
                                 darc_to_elmt_idx,
                                 PDM_OWNERSHIP_KEEP);

    }
  } else if(dmesh->n_g_face != 0) { // Donc 2D

    dn_node = dmesh->dn_face;
    dn_arc  = dmesh->dn_edge;
    distrib_arc  = PDM_compute_entity_distribution(comm, dn_arc );

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                               &darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                               &delmt_to_arc,
                               &delmt_to_arc_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    if(darc_to_elmt_tmp == NULL) {
      assert(delmt_to_arc_idx != NULL);
      PDM_dcellface_to_dfacecell(distrib_arc,
                                 distrib_node,
                                 delmt_to_arc_idx,
                                 delmt_to_arc,
                                 &darc_to_elmt_tmp,
                                 comm);

      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                 darc_to_elmt_tmp,
                                 darc_to_elmt_idx,
                                 PDM_OWNERSHIP_KEEP);

    }

  } else if(dmesh->n_g_edge != 0) { // Donc 1D

    dn_node = dmesh->dn_edge;
    dn_arc  = dmesh->dn_vtx;
    is1d    = 1;

    distrib_arc  = PDM_compute_entity_distribution(comm, dn_arc );

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_VTX_EDGE,
                               &darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &delmt_to_arc,
                               &delmt_to_arc_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    if(darc_to_elmt_tmp == NULL) {

      assert(delmt_to_arc_idx == NULL);
      delmt_to_arc_idx = malloc((dn_node+1) * sizeof(int));
      for(int i = 0; i < dn_node+1; ++i) {
        delmt_to_arc_idx[i] = 2*i;
      }

      // PDM_dcellface_to_dfacecell(distrib_arc,
      //                            distrib_node,
      //                            delmt_to_arc_idx,
      //                            delmt_to_arc,
      //                            &darc_to_elmt_tmp,
      //                            comm);
      PDM_dconnectivity_transpose(comm,
                                  distrib_node,
                                  distrib_arc,
                                  delmt_to_arc_idx,
                                  delmt_to_arc,
                                  0,
                                  &darc_to_elmt_idx,
                                  &darc_to_elmt_tmp);

      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_VTX_EDGE,
                                 darc_to_elmt_tmp,
                                 darc_to_elmt_idx,
                                 PDM_OWNERSHIP_KEEP);

    }

  } else if(dmesh->n_g_vtx != 0) { // Donc 0D
    return;
  }

  /*
   * Reminder :
   *  -> 3D : node/elmts = cell | arc = face
   *  -> 2D : node/elmts = face | arc = edge
   *  -> 1D : node/elmts = edge | arc = vtx
   */
  if(is1d == 0) {
    assert(darc_to_elmt_idx == NULL);
    for (int i = 0; i < dn_arc; i++) {
      darc_to_elmt_tmp[2*i+1] = -darc_to_elmt_tmp[2*i+1];
    }

    PDM_setup_connectivity_idx(dn_arc,
                               2,
                               darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               &darc_to_elmt);

    /* Reamke same sign */
    for (int i = 0; i < dn_arc; i++) {
      darc_to_elmt_tmp[2*i+1] = -darc_to_elmt_tmp[2*i+1];
    }
  }

  if(0 == 1) {
    PDM_log_trace_connectivity_long(delmt_to_arc_idx, delmt_to_arc, dn_node, "delmt_to_arc : ");
  }

  if(delmt_to_arc == NULL) {
    if(dmesh->n_g_cell != 0) { // Donc 3D
      PDM_dconnectivity_transpose(comm,
                                  distrib_arc,
                                  distrib_node,
                                  darc_to_elmt_idx,
                                  darc_to_elmt,
                                  1,
                                  &delmt_to_arc_idx,
                                  &delmt_to_arc);

      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 delmt_to_arc,
                                 delmt_to_arc_idx,
                                 PDM_OWNERSHIP_KEEP);

    } else if(dmesh->n_g_face != 0) {
      PDM_dconnectivity_transpose(comm,
                                  distrib_arc,
                                  distrib_node,
                                  darc_to_elmt_idx,
                                  darc_to_elmt,
                                  1,
                                  &delmt_to_arc_idx,
                                  &delmt_to_arc);

      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 delmt_to_arc,
                                 delmt_to_arc_idx,
                                 PDM_OWNERSHIP_KEEP);

    } else if(dmesh->n_g_edge != 0) {
      // PDM_dconnectivity_transpose(comm,
      //                             distrib_node,
      //                             distrib_arc,
      //                             darc_to_elmt_idx,
      //                             darc_to_elmt,
      //                             1,
      //                             &delmt_to_arc_idx,
      //                             &delmt_to_arc);

      // PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_VTX_EDGE,
      //                            delmt_to_arc,
      //                            delmt_to_arc_idx,
      //                            PDM_OWNERSHIP_KEEP);
    }
  }

  if(compute_dual == 1) {
    if(is1d == 1) {
      PDM_deduce_combine_connectivity_dual(comm,
                                           distrib_node,
                                           distrib_arc,
                                           delmt_to_arc_idx,
                                           delmt_to_arc,
                                           darc_to_elmt_idx,
                                           darc_to_elmt_tmp,
                                           1, // is signed
                                           &dual_graph_idx,
                                           &dual_graph);
    } else {
      PDM_deduce_combine_connectivity_dual(comm,
                                           distrib_node,
                                           distrib_arc,
                                           delmt_to_arc_idx,
                                           delmt_to_arc,
                                           darc_to_elmt_idx,
                                           darc_to_elmt,
                                           1, // is signed
                                           &dual_graph_idx,
                                           &dual_graph);
    }

    /* Shift to 0 dual */
    for(int i = 0; i < dual_graph_idx[dn_node]; ++i) {
      dual_graph[i] = dual_graph[i] - 1;
    }
  }
  if(is1d == 0) {
    free(darc_to_elmt_idx);
    free(darc_to_elmt);
  } else {
    free(delmt_to_arc_idx);
  }
  free(distrib_arc);

  *out_dual_graph_idx = dual_graph_idx;
  *out_dual_graph     = dual_graph;

}


static
PDM_g_num_t*
_split_graph
(
      PDM_MPI_Comm       comm,
      PDM_dmesh_t       *dmesh,
      _part_mesh_t      *pmeshes,
      int                n_part,
      PDM_split_dual_t   split_method,
      PDM_part_size_t    part_size_method,
const double            *part_fraction,
      PDM_g_num_t       *distrib_node,
      int              **node_part
)
{
  int verbose = 0;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Split only graph by the most high entity : (cell for 3D / face for 2D)
   *  The most high level entitvy is named elmt
   */
  int  dn_cell = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_CELL);
  int  dn_face = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_FACE);
  int  dn_edge = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_EDGE);
  int  dn_vtx  = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VTX);

  if(verbose && i_rank == 0) {
    printf(" dn_cell = %i \n", dn_cell);
    printf(" dn_face = %i \n", dn_face);
    printf(" dn_edge = %i \n", dn_edge);
    printf(" dn_vtx  = %i \n", dn_vtx );
    printf(" n_g_cell = "PDM_FMT_G_NUM" \n", dmesh->n_g_cell );
    printf(" n_g_face = "PDM_FMT_G_NUM" \n", dmesh->n_g_face );
    printf(" n_g_edge = "PDM_FMT_G_NUM" \n", dmesh->n_g_edge );
    printf(" n_g_vtx  = "PDM_FMT_G_NUM" \n", dmesh->n_g_vtx  );
  }

  int compute_dual = 1;
  if(split_method == PDM_SPLIT_DUAL_WITH_HILBERT ||
     split_method == PDM_SPLIT_DUAL_WITH_IMPLICIT) {
    compute_dual = 0;
  }

  PDM_g_num_t *dual_graph_idx = NULL;
  PDM_g_num_t *dual_graph     = NULL;
  _warm_up_for_split(comm,
                     dmesh,
                     distrib_node,
                     compute_dual,
                     &dual_graph_idx,
                     &dual_graph);

  int dn_node = distrib_node[i_rank+1] - distrib_node[i_rank];
  int *_node_part = malloc(dn_node * sizeof(int));

  // Compute total number of partitions for this domain
  int tn_part;
  PDM_MPI_Allreduce(&n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmeshes->tn_part = tn_part;

  PDM_g_num_t *distrib_partition = PDM_compute_entity_distribution(comm, n_part );
  double *part_fractions = NULL;
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS){
    int *n_part_per_rank = (int    *) malloc( n_rank * sizeof(int   ));
    int *displ           = (int    *) malloc( n_rank * sizeof(int   ));
    part_fractions       = (double *) malloc(tn_part * sizeof(double));
    for (int i =0; i < n_rank; i++){
      n_part_per_rank[i] = distrib_partition[i+1] - distrib_partition[i];
      displ[i] = distrib_partition[i];
    }

    PDM_MPI_Allgatherv((void*) part_fraction,
                       n_part,
                       PDM_MPI_DOUBLE,
                       part_fractions,
                       n_part_per_rank,
                       displ,
                       PDM_MPI_DOUBLE,
                       comm);
    free(n_part_per_rank);
    free(displ);
  }

  if (split_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    _split_graph_hilbert(comm,
                         dmesh,
                         n_part,
                         _node_part);
  } else if(split_method == PDM_SPLIT_DUAL_WITH_IMPLICIT) {
    assert(n_part == 1);
    for(int i = 0; i < dn_node; ++i) {
      _node_part[i] = i_rank;
    }
  } else {
    PDM_para_graph_split (split_method,
                          distrib_node,
                          dual_graph_idx,
                          dual_graph,
                          NULL,
                          NULL,
                          tn_part,
                          part_fractions,
                          _node_part,
                          comm);
  }

  // PDM_log_trace_array_int (_node_part, dn_node, "_node_part :: ");
  if(compute_dual == 1) {
    free(dual_graph_idx);
    free(dual_graph);
  }
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS) {
    free(part_fractions);
  }

  *node_part = _node_part;

  return distrib_partition;
}

static
void
_deduce_part_face_connectivity
(
  PDM_MPI_Comm       comm,
  PDM_dmesh_t       *dmesh,
  _part_mesh_t      *pmeshes,
  int                n_part,
  int               *pn_face,
  PDM_g_num_t      **pface_ln_to_gn,
  PDM_g_num_t       *face_distrib,
  PDM_g_num_t       *edge_distrib,
  int             ***out_pface_vtx_idx,
  int             ***out_pface_vtx,
  int             ***out_pface_edge_idx,
  int             ***out_pface_edge,
  int              **out_pn_edge,
  PDM_g_num_t     ***out_pedge_ln_to_gn,
  int             ***out_pedge_vtx_idx,
  int             ***out_pedge_vtx,
  int              **out_pn_vtx,
  PDM_g_num_t     ***out_pvtx_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int          *pn_edge        = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  int         **pedge_vtx_idx  = NULL;
  int         **pedge_vtx      = NULL;

  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;

  int         **pface_vtx_idx = NULL;
  int         **pface_vtx     = NULL;

  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;


  int from_face_edge = 0;
  int from_face_vtx  = 0;
  // face vtx
  PDM_g_num_t *dface_vtx     = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dmesh,
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
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_edge_idx != NULL) {
    from_face_edge = 1;
  }

  if(from_face_edge == 1) {

    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 face_distrib,
                                                 dface_edge_idx,
                                                 dface_edge,
                                                 n_part,
                                                 pn_face,
                          (const PDM_g_num_t **) pface_ln_to_gn,
                                                 &pn_edge,
                                                 &pedge_ln_to_gn,
                                                 &pface_edge_idx,
                                                 &pface_edge);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_FACE, pn_face[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                     pface_edge    [i_part],
                                     pface_edge_idx[i_part],
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        pface_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
    }

    // edge_vtx
    PDM_g_num_t *dedge_vtx     = NULL;
    int         *dedge_vtx_idx = NULL;
    PDM_dmesh_connectivity_get(dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
    int *_dedge_vtx_idx = NULL;
    if(dedge_vtx_idx == NULL)  {
      int dn_edge = edge_distrib[i_rank+1] - edge_distrib[i_rank];
      _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
        _dedge_vtx_idx[i_edge] = 2*i_edge;
      }
    } else {
      _dedge_vtx_idx = dedge_vtx_idx;
    }

    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 edge_distrib,
                                                 _dedge_vtx_idx,
                                                 dedge_vtx,
                                                 n_part,
                                                 pn_edge,
                           (const PDM_g_num_t **) pedge_ln_to_gn,
                                                 &pn_vtx,
                                                 &pvtx_ln_to_gn,
                                                 &pedge_vtx_idx,
                                                 &pedge_vtx);
    if(dedge_vtx_idx == NULL)  {
      free(_dedge_vtx_idx);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_EDGE, pn_edge[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                     pedge_vtx     [i_part],
                                     NULL,
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_EDGE,
                                        pedge_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
      free(pedge_vtx_idx [i_part]);
    }
  } else if(from_face_vtx == 1) {
    // PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dmesh->dn_face, "dface_vtx ::");
    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 face_distrib,
                                                 dface_vtx_idx,
                                                 dface_vtx,
                                                 n_part,
                                                 pn_face,
                          (const PDM_g_num_t **) pface_ln_to_gn,
                                                 &pn_vtx,
                                                 &pvtx_ln_to_gn,
                                                 &pface_vtx_idx,
                                                 &pface_vtx);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_FACE, pn_face[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     pface_vtx     [i_part],
                                     pface_vtx_idx [i_part],
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        pface_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
    }
  }

  *out_pface_vtx_idx  = pface_vtx_idx;
  *out_pface_vtx      = pface_vtx;
  *out_pface_edge_idx = pface_edge_idx;
  *out_pface_edge     = pface_edge;
  *out_pn_edge        = pn_edge;
  *out_pedge_ln_to_gn = pedge_ln_to_gn;
  *out_pedge_vtx_idx  = pedge_vtx_idx;
  *out_pedge_vtx      = pedge_vtx;
  *out_pn_vtx         = pn_vtx;
  *out_pvtx_ln_to_gn  = pvtx_ln_to_gn;

}

static
void
_rebuild_part_mesh_group
(
 PDM_dmesh_t       *dmesh,
_part_mesh_t       *pmeshes,
 int                n_part,
 PDM_bound_type_t   entity_bound,
 PDM_g_num_t       *entity_distrib,
 int               *pn_entity,
 PDM_g_num_t      **entity_ln_to_gn,
 PDM_MPI_Comm       comm
)
{

  PDM_g_num_t *dentity_bound     = NULL;
  int         *dentity_bound_idx = NULL;
  int n_entity_group = PDM_dmesh_bound_get(dmesh,
                                           entity_bound,
                                           &dentity_bound,
                                           &dentity_bound_idx,
                                           PDM_OWNERSHIP_BAD_VALUE);

  if (n_entity_group == 0) {
    return;
  }


  int         **pentity_bound_idx               = NULL;
  int         **pentity_bound                   = NULL;
  PDM_g_num_t **pentity_bound_ln_to_gn          = NULL;
  PDM_part_distgroup_to_partgroup(comm,
                                  entity_distrib,
                                  n_entity_group,
                                  dentity_bound_idx,
                                  dentity_bound,
                                  n_part,
                                  pn_entity,
           (const PDM_g_num_t **) entity_ln_to_gn,
                                 &pentity_bound_idx,
                                 &pentity_bound,
                                 &pentity_bound_ln_to_gn);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_bound_concat_set(pmeshes->pmesh,
                                   i_part,
                                   entity_bound,
                                   n_entity_group,
                                   pentity_bound_idx     [i_part],
                                   pentity_bound         [i_part],
                                   pentity_bound_ln_to_gn[i_part],
                                   PDM_OWNERSHIP_KEEP);
  }

  free(pentity_bound_idx               );
  free(pentity_bound                   );
  free(pentity_bound_ln_to_gn          );
}

static
void
_deduce_part_connectivity_0d
(
 PDM_MPI_Comm       comm,
 PDM_dmesh_t       *dmesh,
 _part_mesh_t      *pmeshes,
 int                n_part,
 int               *pn_vtx,
 PDM_g_num_t      **pvtx_ln_to_gn,
 PDM_g_num_t       *vtx_distrib
)
{
  // Vertex
  double *dvtx_coord = NULL;
  PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);
  double      **pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        vtx_distrib,
                                        dvtx_coord,
                                        pn_vtx,
                 (const PDM_g_num_t **) pvtx_ln_to_gn,
                                       &pvtx_coord);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_VTX, pn_vtx[i_part]);
    PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      pvtx_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_vtx_coord_set(pmeshes->pmesh,
                                i_part,
                                pvtx_coord   [i_part],
                                PDM_OWNERSHIP_KEEP);
  }

  free(pvtx_coord);
}

static
void
_deduce_part_connectivity_1d
(
 PDM_MPI_Comm    comm,
 PDM_dmesh_t    *dmesh,
 _part_mesh_t   *pmeshes,
 int             n_part,
 int            *pn_edge,
 PDM_g_num_t   **pedge_ln_to_gn,
 PDM_g_num_t    *edge_distrib,
 PDM_g_num_t    *vtx_distrib,
 int          ***out_pedge_vtx_idx,
 int          ***out_pedge_vtx,
 int           **out_pn_vtx,
 PDM_g_num_t  ***out_pvtx_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  /*
   * Reconstruction edge_vtx
   */
  PDM_g_num_t *dedge_vtx     = NULL;
  int         *dedge_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             &dedge_vtx,
                             &dedge_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);
  int *_dedge_vtx_idx = NULL;
  if(dedge_vtx_idx == NULL)  {
    int dn_edge = edge_distrib[i_rank+1] - edge_distrib[i_rank];
    _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
    for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
      _dedge_vtx_idx[i_edge] = 2*i_edge;
    }
  } else {
    _dedge_vtx_idx = dedge_vtx_idx;
  }

  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;

  int         **pedge_vtx_idx = NULL;
  int         **pedge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               edge_distrib,
                                               _dedge_vtx_idx,
                                               dedge_vtx,
                                               n_part,
                                               pn_edge,
                        (const PDM_g_num_t **) pedge_ln_to_gn,
                                               &pn_vtx,
                                               &pvtx_ln_to_gn,
                                               &pedge_vtx_idx,
                                               &pedge_vtx);
  if(dedge_vtx_idx == NULL)  {
    free(_dedge_vtx_idx);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_part_mesh_n_entity_set(pmeshes->pmesh,
                               i_part,
                               PDM_MESH_ENTITY_EDGE,
                               pn_edge[i_part]);

    PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   pedge_vtx    [i_part],
                                   NULL,
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      pedge_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);

    free(pedge_vtx_idx[i_part]);
    pedge_vtx_idx[i_part] = NULL;
  }

  _deduce_part_connectivity_0d(comm,
                               dmesh,
                               pmeshes,
                               n_part,
                               pn_vtx,
                               pvtx_ln_to_gn,
                               vtx_distrib);

  *out_pedge_vtx_idx  = pedge_vtx_idx;
  *out_pedge_vtx      = pedge_vtx;
  *out_pn_vtx         = pn_vtx;
  *out_pvtx_ln_to_gn  = pvtx_ln_to_gn;

}


static
void
_deduce_part_connectivity_3d
(
 PDM_MPI_Comm       comm,
 PDM_dmesh_t       *dmesh,
 _part_mesh_t      *pmeshes,
 int                n_part,
 int               *pn_cell,
 PDM_g_num_t      **pcell_ln_to_gn,
 PDM_g_num_t       *cell_distrib,
 PDM_g_num_t       *face_distrib,
 PDM_g_num_t       *edge_distrib,
 PDM_g_num_t       *vtx_distrib,
 int              **out_pn_face,
 PDM_g_num_t     ***out_pface_ln_to_gn,
 int             ***out_pface_vtx_idx,
 int             ***out_pface_vtx,
 int             ***out_pface_edge_idx,
 int             ***out_pface_edge,
 int              **out_pn_edge,
 PDM_g_num_t     ***out_pedge_ln_to_gn,
 int             ***out_pedge_vtx_idx,
 int             ***out_pedge_vtx,
 int              **out_pn_vtx,
 PDM_g_num_t     ***out_pvtx_ln_to_gn
)
{
  /*
   * Reconstruction cell_face
   */
  PDM_g_num_t *dcell_face     = NULL;
  int         *dcell_face_idx = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_CELL_FACE,
                             &dcell_face,
                             &dcell_face_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  int          *pn_face        = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;

  int         **pcell_face_idx = NULL;
  int         **pcell_face     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               cell_distrib,
                                               dcell_face_idx,
                                               dcell_face,
                                               n_part,
                                               pn_cell,
                        (const PDM_g_num_t **) pcell_ln_to_gn,
                                               &pn_face,
                                               &pface_ln_to_gn,
                                               &pcell_face_idx,
                                               &pcell_face);


  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_part_mesh_n_entity_set(pmeshes->pmesh,
                               i_part,
                               PDM_MESH_ENTITY_CELL,
                               pn_cell[i_part]);

    PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   pcell_face    [i_part],
                                   pcell_face_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      pcell_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);

  }

  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;

  int          *pn_edge        = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  int         **pedge_vtx_idx  = NULL;
  int         **pedge_vtx      = NULL;

  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;

  int         **pface_vtx_idx = NULL;
  int         **pface_vtx     = NULL;
  _deduce_part_face_connectivity(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_face,
                                 pface_ln_to_gn,
                                 face_distrib,
                                 edge_distrib,
                                 &pface_vtx_idx,
                                 &pface_vtx,
                                 &pface_edge_idx,
                                 &pface_edge,
                                 &pn_edge,
                                 &pedge_ln_to_gn,
                                 &pedge_vtx_idx,
                                 &pedge_vtx,
                                 &pn_vtx,
                                 &pvtx_ln_to_gn);

  // Vertex
  _deduce_part_connectivity_0d(comm,
                               dmesh,
                               pmeshes,
                               n_part,
                               pn_vtx,
                               pvtx_ln_to_gn,
                               vtx_distrib);

  /*
   * Setup face_cell
   */
  int **pface_cell = NULL;
  PDM_part_reverse_pcellface(n_part,
                             pn_cell,
                             pn_face,
              (const int **) pcell_face_idx,
              (const int **) pcell_face,
                             &pface_cell);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                   pface_cell     [i_part],
                                   NULL,
                                   PDM_OWNERSHIP_KEEP);
  }
  free(pface_cell);
  free(pcell_face_idx);
  free(pcell_face);


  /*
   * Group (usefull for ordering)
   */
  _rebuild_part_mesh_group(dmesh,
                           pmeshes,
                           n_part,
                           PDM_BOUND_TYPE_FACE,
                           face_distrib,
                           pn_face,
                           pface_ln_to_gn,
                           comm);

  if (edge_distrib != NULL) {
    _rebuild_part_mesh_group(dmesh,
                             pmeshes,
                             n_part,
                             PDM_BOUND_TYPE_EDGE,
                             edge_distrib,
                             pn_edge,
                             pedge_ln_to_gn,
                             comm);
  }

  _rebuild_part_mesh_group(dmesh,
                           pmeshes,
                           n_part,
                           PDM_BOUND_TYPE_VTX,
                           vtx_distrib,
                           pn_vtx,
                           pvtx_ln_to_gn,
                           comm);

  *out_pn_face        = pn_face;
  *out_pface_ln_to_gn = pface_ln_to_gn;
  *out_pface_vtx_idx  = pface_vtx_idx;
  *out_pface_vtx      = pface_vtx;
  *out_pface_edge_idx = pface_edge_idx;
  *out_pface_edge     = pface_edge;
  *out_pn_edge        = pn_edge;
  *out_pedge_ln_to_gn = pedge_ln_to_gn;
  *out_pedge_vtx_idx  = pedge_vtx_idx;
  *out_pedge_vtx      = pedge_vtx;
  *out_pn_vtx         = pn_vtx;
  *out_pvtx_ln_to_gn  = pvtx_ln_to_gn;
}



static
void
_deduce_part_connectivity_2d
(
 PDM_MPI_Comm       comm,
 PDM_dmesh_t       *dmesh,
 _part_mesh_t      *pmeshes,
 int                n_part,
 int               *pn_face,
 PDM_g_num_t      **pface_ln_to_gn,
 PDM_g_num_t       *face_distrib,
 PDM_g_num_t       *edge_distrib,
 PDM_g_num_t       *vtx_distrib,
 int             ***out_pface_vtx_idx,
 int             ***out_pface_vtx,
 int             ***out_pface_edge_idx,
 int             ***out_pface_edge,
 int              **out_pn_edge,
 PDM_g_num_t     ***out_pedge_ln_to_gn,
 int             ***out_pedge_vtx_idx,
 int             ***out_pedge_vtx,
 int              **out_pn_vtx,
 PDM_g_num_t     ***out_pvtx_ln_to_gn
)
{
  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;

  int          *pn_edge        = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  int         **pedge_vtx_idx  = NULL;
  int         **pedge_vtx      = NULL;

  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;

  int         **pface_vtx_idx = NULL;
  int         **pface_vtx     = NULL;

  _deduce_part_face_connectivity(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_face,
                                 pface_ln_to_gn,
                                 face_distrib,
                                 edge_distrib,
                                 &pface_vtx_idx,
                                 &pface_vtx,
                                 &pface_edge_idx,
                                 &pface_edge,
                                 &pn_edge,
                                 &pedge_ln_to_gn,
                                 &pedge_vtx_idx,
                                 &pedge_vtx,
                                 &pn_vtx,
                                 &pvtx_ln_to_gn);

  // Vertex
  _deduce_part_connectivity_0d(comm,
                               dmesh,
                               pmeshes,
                               n_part,
                               pn_vtx,
                               pvtx_ln_to_gn,
                               vtx_distrib);

  int **pedge_face = NULL;
  PDM_part_reverse_pcellface(n_part,
                             pn_face,
                             pn_edge,
              (const int **) pface_edge_idx,
              (const int **) pface_edge,
                             &pedge_face);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                   pedge_face     [i_part],
                                   NULL,
                                   PDM_OWNERSHIP_KEEP);
  }
  free(pedge_face);

  /*
   * Group (usefull for ordering)
   */
  _rebuild_part_mesh_group(dmesh,
                           pmeshes,
                           n_part,
                           PDM_BOUND_TYPE_EDGE,
                           edge_distrib,
                           pn_edge,
                           pedge_ln_to_gn,
                           comm);


  *out_pface_vtx_idx  = pface_vtx_idx;
  *out_pface_vtx      = pface_vtx;
  *out_pface_edge_idx = pface_edge_idx;
  *out_pface_edge     = pface_edge;
  *out_pn_edge        = pn_edge;
  *out_pedge_ln_to_gn = pedge_ln_to_gn;
  *out_pedge_vtx_idx  = pedge_vtx_idx;
  *out_pedge_vtx      = pedge_vtx;
  *out_pn_vtx         = pn_vtx;
  *out_pvtx_ln_to_gn  = pvtx_ln_to_gn;
}



static
void
_run_ppart_domain
(
PDM_dmesh_t       *dmesh,
PDM_dmesh_nodal_t *dmesh_nodal,
_part_mesh_t      *pmeshes,
int                n_part,
PDM_split_dual_t   split_method,
PDM_part_size_t    part_size_method,
const double*      part_fraction,
PDM_MPI_Comm       comm
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // int  dn_cell = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_CELL);
  int  dn_face = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_FACE);
  int  dn_edge = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_EDGE);
  int  dn_vtx  = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VTX);

  int dn_node = 0;
  if(dmesh->n_g_cell != 0) {
    assert(dmesh->dn_cell > 0);
    dn_node = dmesh->dn_cell;
  } else if (dmesh->n_g_face != 0) {
    dn_node = dmesh->dn_face;
  } else if (dmesh->n_g_edge != 0) {
    dn_node = dmesh->dn_edge;
  } else if (dmesh->n_g_vtx != 0) {
    dn_node = dmesh->dn_vtx;
  } else {
    dn_node = 0;
  }

  /*
   *  Split graph (manage 3D/2D automaticaly)
   */
  int *node_part = NULL;
  PDM_g_num_t *distrib_node = PDM_compute_entity_distribution(comm, dn_node);
  PDM_g_num_t* distrib_partition = _split_graph(comm,
                                                dmesh,
                                                pmeshes,
                                                n_part,
                                                split_method,
                                                part_size_method,
                                                part_fraction,
                                                distrib_node,
                                                &node_part);

  /*
   * Deduce node_ln_to_gn
   */
  int          *pn_node        = NULL;
  PDM_g_num_t **pnode_ln_to_gn = NULL;
  PDM_part_assemble_partitions(comm,
                               distrib_partition,
                               distrib_node,
                               node_part,
                               NULL,
                               NULL,
                              &pn_node,
                              &pnode_ln_to_gn,
                               NULL);
  free(node_part);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pnode_ln_to_gn[i_part], pn_node[i_part], "pnode_ln_to_gn :: ");
    }
  }

  PDM_g_num_t *face_distrib = NULL;
  PDM_g_num_t *edge_distrib = NULL;
  PDM_g_num_t *vtx_distrib  = PDM_compute_entity_distribution(comm, dn_vtx);

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &edge_distrib);
  int own_edge_distrib = 0;
  if(edge_distrib == NULL) {
    own_edge_distrib = 1;
    edge_distrib = PDM_compute_entity_distribution(comm, dn_edge);
  }

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &face_distrib);
  int own_face_distrib = 0;
  if(face_distrib == NULL) {
    own_face_distrib = 1;
    face_distrib = PDM_compute_entity_distribution(comm, dn_face);
  }

  /*
   *  Deduce all required connectivity by descending connectivity
   */
  int          *pn_cell        = NULL;
  PDM_g_num_t **pcell_ln_to_gn = NULL;

  int          *pn_face        = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;

  int         **pcell_face_idx = NULL;
  int         **pcell_face     = NULL;

  int          *pn_edge        = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  int         **pedge_vtx_idx  = NULL;
  int         **pedge_vtx      = NULL;

  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;

  int         **pface_vtx_idx = NULL;
  int         **pface_vtx     = NULL;

  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;

  if(dmesh->n_g_cell != 0) {

    _deduce_part_connectivity_3d(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_node,
                                 pnode_ln_to_gn,
                                 distrib_node,
                                 face_distrib,
                                 edge_distrib,
                                 vtx_distrib,
                                 &pn_face,
                                 &pface_ln_to_gn,
                                 &pface_vtx_idx,
                                 &pface_vtx,
                                 &pface_edge_idx,
                                 &pface_edge,
                                 &pn_edge,
                                 &pedge_ln_to_gn,
                                 &pedge_vtx_idx,
                                 &pedge_vtx,
                                 &pn_vtx,
                                 &pvtx_ln_to_gn);

    pn_cell        = pn_node;
    pcell_ln_to_gn = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;

  } else if (dmesh->n_g_face != 0) {
    _deduce_part_connectivity_2d(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_node,
                                 pnode_ln_to_gn,
                                 distrib_node, // = face_distrib
                                 edge_distrib,
                                 vtx_distrib,
                                 &pface_vtx_idx,
                                 &pface_vtx,
                                 &pface_edge_idx,
                                 &pface_edge,
                                 &pn_edge,
                                 &pedge_ln_to_gn,
                                 &pedge_vtx_idx,
                                 &pedge_vtx,
                                 &pn_vtx,
                                 &pvtx_ln_to_gn);

    pn_face        = pn_node;
    pface_ln_to_gn = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;

  } else if (dmesh->n_g_edge != 0) {
    _deduce_part_connectivity_1d(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_node,
                                 pnode_ln_to_gn,
                                 edge_distrib,
                                 vtx_distrib,
                                 &pedge_vtx_idx,
                                 &pedge_vtx,
                                 &pn_vtx,
                                 &pvtx_ln_to_gn);

    pn_edge        = pn_node;
    pedge_ln_to_gn = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;

  } else if (dmesh->n_g_vtx != 0) {
    _deduce_part_connectivity_0d(comm,
                                 dmesh,
                                 pmeshes,
                                 n_part,
                                 pn_node,
                                 pnode_ln_to_gn,
                                 vtx_distrib);

    pn_vtx         = pn_node;
    pvtx_ln_to_gn  = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;

  } else {
    dn_node = 0;
  }

  /*
   * Force // ordering of vertex (needed by other ordering method)
   */
  int         **pinternal_vtx_bound_proc_idx  = NULL;
  int         **pinternal_vtx_bound_part_idx  = NULL;
  int         **pinternal_vtx_bound           = NULL;
  int         **pinternal_vtx_priority        = NULL;
  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      vtx_distrib,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                     &pinternal_vtx_priority);

  _setup_ghost_information(i_rank,
                           n_part,
                           pn_vtx,
                           pinternal_vtx_priority,
                           distrib_partition);

  /* Free in order to be correclty */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pinternal_vtx_bound_proc_idx[ipart]);
    free(pinternal_vtx_bound_part_idx[ipart]);
    free(pinternal_vtx_bound[ipart]);
    // free(pinternal_vtx_priority[ipart]);
  }
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);

  /* MAP on _part_t to reuse ordering */
  _part_t** parts = _map_part_t_with_part_mesh(pmeshes->pmesh);
  for (int i_part = 0; i_part < n_part; i_part++) {
    parts[i_part]->vtx_ghost_information = pinternal_vtx_priority[i_part];
    pmeshes->vtx_ghost_information[i_part] = pinternal_vtx_priority[i_part];
  }

  /*
   * Real re-numebering
   */
  PDM_part_renum_cell(parts, n_part, pmeshes->renum_cell_method, (void *) pmeshes->renum_cell_properties);
  PDM_part_renum_face(parts, n_part, pmeshes->renum_face_method, NULL);
  PDM_part_renum_edge(parts, n_part, pmeshes->renum_edge_method, NULL);
  PDM_part_renum_vtx (parts, n_part, pmeshes->renum_vtx_method , (void *) pinternal_vtx_priority);
  free(pinternal_vtx_priority);

  for (int i_part = 0; i_part < n_part; i_part++) {

    if(parts[i_part]->cell_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_CELL,
                                     parts[i_part]->cell_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->face_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_FACE,
                                     parts[i_part]->face_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->edge_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_EDGE,
                                     parts[i_part]->edge_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->vtx_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_VTX,
                                     parts[i_part]->vtx_color,
                                     PDM_OWNERSHIP_KEEP);
    }

    pmeshes->hyperplane_color[i_part] = parts[i_part]->hyperplane_color;
    pmeshes->thread_color    [i_part] = parts[i_part]->thread_color;

    _part_free(parts[i_part]);
  }
  free(parts);

  /*
   * All entities are reorder - In case of HO mesh we need to append all ho vtx in vtx_ln_to_gn AND pvtx_coord
   */
  int have_ho = 0;
  if(dmesh_nodal != NULL) {
    have_ho =  PDM_dmesh_nodal_have_ho(dmesh_nodal);
  }

  if(have_ho == 1) {

    /* Deduce vtx from the connectivity by elmt */
    int          *pn_vtx_all       = NULL;
    PDM_g_num_t **vtx_all_ln_to_gn = NULL;
    PDM_generate_ho_vtx_ln_to_gn(dmesh_nodal,
                                 n_part,
                                 pn_cell,
                                 pcell_ln_to_gn,
                                 pn_face,
                                 pface_ln_to_gn,
                                 pn_edge,
                                 pedge_ln_to_gn,
                                 pn_vtx,
                                 pvtx_ln_to_gn,
                                 &pn_vtx_all,
                                 &vtx_all_ln_to_gn);
    if(0 == 1) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        log_trace("pn_vtx_all[%i] = %i \n", i_part, pn_vtx_all[i_part]);
      }
    }

    for (int i_part = 0; i_part < n_part; i_part++) {

      PDM_g_num_t *tmp_vtx_ln_to_gn = NULL;
      double      *tmp_vtx_coord    = NULL;
      PDM_part_mesh_entity_ln_to_gn_get(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &tmp_vtx_ln_to_gn,
                                        PDM_OWNERSHIP_USER);
      PDM_part_mesh_vtx_coord_get(pmeshes->pmesh,
                                  i_part,
                                  &tmp_vtx_coord,
                                  PDM_OWNERSHIP_USER);
      free(tmp_vtx_ln_to_gn);
      free(tmp_vtx_coord);

      pvtx_ln_to_gn[i_part] = vtx_all_ln_to_gn[i_part];
    }
    // free(pvtx_ln_to_gn);
    double** pvtx_coord = NULL;

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          n_part,
                                          vtx_distrib,
                                          dvtx_coord,
                                          pn_vtx_all,
                   (const PDM_g_num_t **) vtx_all_ln_to_gn,
                                         &pvtx_coord);


    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_VTX, pn_vtx_all[i_part]);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        pvtx_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_vtx_coord_set(pmeshes->pmesh,
                                  i_part,
                                  pvtx_coord   [i_part],
                                  PDM_OWNERSHIP_KEEP);
    }
    free(pvtx_coord);
    free(vtx_all_ln_to_gn);
    free(pn_vtx_all);
  }

  /*
   *  All data has been reorder, we can now and only now setup desired comm graph
   */
  if(pn_cell != NULL){
    int         **pinternal_face_bound_proc_idx = NULL;
    int         **pinternal_face_bound_part_idx = NULL;
    int         **pinternal_face_bound          = NULL;
    PDM_part_generate_entity_graph_comm(comm,
                                        distrib_partition,
                                        face_distrib,
                                        n_part,
                                        pn_face,
                 (const PDM_g_num_t **) pface_ln_to_gn,
                                        NULL,
                                       &pinternal_face_bound_proc_idx,
                                       &pinternal_face_bound_part_idx,
                                       &pinternal_face_bound,
                                        NULL);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_BOUND_TYPE_FACE,
                                        pinternal_face_bound_proc_idx[i_part],
                                        pinternal_face_bound_part_idx[i_part],
                                        pinternal_face_bound[i_part],
                                        PDM_OWNERSHIP_KEEP);

    }
    free(pinternal_face_bound_proc_idx );
    free(pinternal_face_bound_part_idx );
    free(pinternal_face_bound          );

  } else if(pn_edge  != NULL) {
    int **pinternal_edge_bound_proc_idx = NULL;
    int **pinternal_edge_bound_part_idx = NULL;
    int **pinternal_edge_bound          = NULL;
    PDM_part_generate_entity_graph_comm(comm,
                                        distrib_partition,
                                        edge_distrib,
                                        n_part,
                                        pn_edge,
                 (const PDM_g_num_t **) pedge_ln_to_gn,
                                        NULL,
                                       &pinternal_edge_bound_proc_idx,
                                       &pinternal_edge_bound_part_idx,
                                       &pinternal_edge_bound,
                                        NULL);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_BOUND_TYPE_EDGE,
                                        pinternal_edge_bound_proc_idx[i_part],
                                        pinternal_edge_bound_part_idx[i_part],
                                        pinternal_edge_bound[i_part],
                                        PDM_OWNERSHIP_KEEP);
    }
    free(pinternal_edge_bound_proc_idx );
    free(pinternal_edge_bound_part_idx );
    free(pinternal_edge_bound          );
  }

  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      vtx_distrib,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                      NULL);

  // Finally complete parts structure with internal join data and bounds
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_BOUND_TYPE_VTX,
                                      pinternal_vtx_bound_proc_idx[i_part],
                                      pinternal_vtx_bound_part_idx[i_part],
                                      pinternal_vtx_bound[i_part],
                                      PDM_OWNERSHIP_KEEP);

  }
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);

  if(own_edge_distrib == 1) {
    free(edge_distrib);
  }
  if(own_face_distrib == 1) {
    free(face_distrib);
  }

  if(pn_cell != NULL) {
    free(pn_cell);
    free(pcell_ln_to_gn);
    free(pcell_face_idx);
    free(pcell_face);
  }
  if(pn_face != NULL) {
    free(pn_face);
    free(pface_ln_to_gn);
    free(pface_edge_idx);
    free(pface_edge);
    if(pface_vtx_idx  !=  NULL) {
      free(pface_vtx_idx);
      free(pface_vtx);
    }
  }
  if(pn_edge != NULL) {
    free(pn_edge);
    free(pedge_ln_to_gn);
    free(pedge_vtx_idx);
    free(pedge_vtx);
  }
  if(pn_vtx != NULL) {
    free(pn_vtx);
    free(pvtx_ln_to_gn);
    // free(pvtx_coord);
  }
  free(distrib_node);
  free(vtx_distrib);
  free(distrib_partition);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure. This method allows to split multiple domains
 *
 * \param [in]   n_domain         Number of domains in the original mesh
 * \param [in]   n_part           Number of partition per proc in each domain
 * \param [in]   merge_blocks     Merge or not the domains before splitting
 * \param [in]   split_method     Choice of library used to split the mesh
 * \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
 * \param [in]   part_weight      Weight (in %) of each partition in heterogeneous case if \ref PDM_part_size_t is set at PDM_PART_SIZE_HETEROGENEOUS
 * \param [in]   comm             PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_multipart_t object
 */
PDM_multipart_t *
PDM_multipart_create
(
 const int              n_domain,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_split_dual_t split_method,
 const PDM_part_size_t  part_size_method,
 const double          *part_fraction,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
)
{
  PDM_multipart_t *multipart = (PDM_multipart_t *) malloc(sizeof(PDM_multipart_t));

  multipart->n_domain = n_domain;
  multipart->n_part   = (int * ) malloc( multipart->n_domain * sizeof(int));

  for (int i = 0; i < multipart->n_domain; ++i) {
    multipart->n_part[i] = n_part[i];
  }

  multipart->merge_blocks     = merge_blocks;
  multipart->split_method     = split_method;
  multipart->part_size_method = part_size_method;
  multipart->part_fraction    = part_fraction;
  multipart->comm             = comm;
  multipart->owner            = owner;

  multipart->n_total_joins    = 0;
  multipart->join_to_opposite = NULL;

  // multipart->dmeshes_ids = (int *) malloc(multipart->n_domain * sizeof(int));

  multipart->dmeshes          = (PDM_dmesh_t                **) malloc(multipart->n_domain * sizeof(PDM_dmesh_t                *));
  multipart->dmeshes_nodal    = (PDM_dmesh_nodal_t          **) malloc(multipart->n_domain * sizeof(PDM_dmesh_nodal_t          *));
  multipart->dmn_to_dm        = (PDM_dmesh_nodal_to_dmesh_t **) malloc(multipart->n_domain * sizeof(PDM_dmesh_nodal_to_dmesh_t *));
  multipart->is_owner_dmeshes = (PDM_bool_t                  *) malloc(multipart->n_domain * sizeof(PDM_bool_t                  ));

  for (int idomain = 0; idomain < multipart->n_domain; ++idomain) {
    multipart->dmeshes_nodal   [idomain] = NULL;
    multipart->dmeshes         [idomain] = NULL;
    multipart->dmn_to_dm       [idomain] = NULL;
    multipart->is_owner_dmeshes[idomain] = PDM_FALSE;
  }

  multipart->pmeshes       = (_part_mesh_t *) malloc(multipart->n_domain * sizeof(_part_mesh_t));

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get("PDM_PART_RENUM_CELL_NONE");
  int _renum_face_method = PDM_part_renum_method_face_idx_get("PDM_PART_RENUM_FACE_NONE");
  int _renum_edge_method = PDM_part_renum_method_edge_idx_get("PDM_PART_RENUM_EDGE_NONE");
  int _renum_vtx_method  = PDM_part_renum_method_vtx_idx_get ("PDM_PART_RENUM_VTX_NONE" );
  for (int idomain = 0; idomain < multipart->n_domain; idomain++) {
    multipart->pmeshes[idomain].renum_cell_method = _renum_cell_method;
    multipart->pmeshes[idomain].renum_face_method = _renum_face_method;
    multipart->pmeshes[idomain].renum_edge_method = _renum_edge_method;
    multipart->pmeshes[idomain].renum_vtx_method  = _renum_vtx_method;
    multipart->pmeshes[idomain].renum_cell_properties = NULL;
    multipart->pmeshes[idomain].joins_ids = NULL;
    multipart->pmeshes[idomain].pmesh     = PDM_part_mesh_create(n_part[idomain], comm);
    multipart->pmeshes[idomain].vtx_ghost_information = malloc(n_part[idomain] * sizeof(int *));
    multipart->pmeshes[idomain].hyperplane_color      = malloc(n_part[idomain] * sizeof(int *));
    multipart->pmeshes[idomain].thread_color          = malloc(n_part[idomain] * sizeof(int *));
    for(int i_part = 0; i_part < n_part[idomain]; ++i_part) {
      multipart->pmeshes[idomain].vtx_ghost_information[i_part] = NULL;
      multipart->pmeshes[idomain].hyperplane_color     [i_part] = NULL;
      multipart->pmeshes[idomain].thread_color         [i_part] = NULL;
    }
    multipart->pmeshes[idomain].is_owner_vtx_ghost_information = PDM_TRUE;
    multipart->pmeshes[idomain].is_owner_hyperplane_color      = PDM_TRUE;
    multipart->pmeshes[idomain].is_owner_thread_color          = PDM_TRUE;
  }

  return (PDM_multipart_t *) multipart;
}

/**
 *
 * \brief Set distributed mesh data for the input domain
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   domain_id      Domain identifier
 * \param [in]   dmesh          Pointer on \ref PDM_dmesh_t containaing all distributed connectivities
 */
void PDM_multipart_dmesh_set
(
 PDM_multipart_t   *multipart,
 const int          domain_id,
       PDM_dmesh_t *dmesh
)
{
  assert(domain_id < multipart->n_domain);
  multipart->dmeshes[domain_id] = dmesh;
}

/**
 *
 * \brief Set distributed mesh data for the input domain. The mesh is describe by nodal connectiviy
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   domain_id        Global domain id
 * \param [in]   dmesh_nodal    Pointer on \ref PDM_dmesh_nodal_t
 */
void PDM_multipart_dmesh_nodal_set
(
 PDM_multipart_t         *multipart,
 const int                domain_id,
       PDM_dmesh_nodal_t *dmesh_nodal
)
{
  assert(domain_id < multipart->n_domain);
  assert(multipart->dmeshes_nodal[domain_id] == NULL);
  multipart->dmeshes_nodal[domain_id] = dmesh_nodal;
}


/**
 * \brief Set block
 *
 * \param [in]   multipart              Pointer to \ref PDM_multipart_t object
 * \param [in]   i_domain               Domain identifier
 * \param [in]   dn_cell                Number of distributed cells
 * \param [in]   dn_face                Number of distributed faces
 * \param [in]   dn_vtx                 Number of distributed vertices
 * \param [in]   n_face_group           Number of face groups
 * \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
 *                                      (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face             Distributed cell face connectivity or NULL
 *                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dface_cell             Distributed face cell connectivity or NULL
 *                                      (size : 2 * dn_face, numbering : 1 to n)
 * \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
 *                                      (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx              Distributed face to vertex connectivity
 *                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord             Distributed vertex coordinates
 *                                      (size : 3*dn_vtx)
 * \param [in]   dface_group_idx        Index of distributed faces list of each group
 *                                      (size = n_face_group + 1) or NULL
 * \param [in]   dface_group            Distributed faces list of each group
 *                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
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
)
{

  // Create dmesh
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     0,
                                     dn_vtx,
                                     multipart->comm);

  PDM_dmesh_vtx_coord_set(dm,
             (double  *)  dvtx_coord,
                          PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
             (PDM_g_num_t *) dface_vtx,
             (int         *) dface_vtx_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_CELL,
             (PDM_g_num_t *) dface_cell,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_CELL_FACE,
             (PDM_g_num_t *) dcell_face,
             (int         *) dcell_face_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_bound_set(dm,
                      PDM_BOUND_TYPE_FACE,
                      n_face_group,
      (PDM_g_num_t *) dface_group,
      (int         *) dface_group_idx,
                      PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(multipart, i_domain, dm);
  multipart->is_owner_dmeshes[i_domain] = PDM_TRUE;
}

void
PDM_multipart_domain_interface_shared_set
(
  PDM_multipart_t        *multipart,
  PDM_domain_interface_t *ditrf
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(ditrf);
  abort();
}


/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]   i_domain              Id of domain which parameters apply (or -1 for all domains)
 * \param [in]   renum_cell_method     Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cacheblocking method :
 *                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
                                        n_vect_face, split_method]
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
)
{

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get(renum_cell_method);
  int _renum_face_method = PDM_part_renum_method_face_idx_get(renum_face_method);
  if (_renum_cell_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
  }
  if (_renum_face_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
  }

  if(i_domain < 0) {
    for (int idomain = 0; idomain < multipart->n_domain; idomain++) {
      multipart->pmeshes[idomain].renum_cell_method = _renum_cell_method;
      multipart->pmeshes[idomain].renum_face_method = _renum_face_method;
      multipart->pmeshes[idomain].renum_cell_properties = renum_cell_properties;
    }
  }
  else {
    assert(i_domain < multipart->n_domain);
    multipart->pmeshes[i_domain].renum_cell_method = _renum_cell_method;
    multipart->pmeshes[i_domain].renum_face_method = _renum_face_method;
    multipart->pmeshes[i_domain].renum_cell_properties = renum_cell_properties;
  }
}
void PDM_multipart_set_reordering_options_vtx
(
 PDM_multipart_t *multipart,
 const int        i_domain,
 const char      *renum_vtx_method
)
{

  int _renum_vtx_method = PDM_part_renum_method_vtx_idx_get(renum_vtx_method);
  if (_renum_vtx_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering vtx method\n", renum_vtx_method);
  }

  if(i_domain < 0) {
    for (int idomain = 0; idomain < multipart->n_domain; idomain++) {
      multipart->pmeshes[idomain].renum_vtx_method = _renum_vtx_method;
    }
  }
  else {
    assert(i_domain < multipart->n_domain);
    multipart->pmeshes[i_domain].renum_vtx_method = _renum_vtx_method;
  }
}

/**
 * \brief Construct the partitioned meshes on every domains
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 */
void
PDM_multipart_compute
(
 PDM_multipart_t *multipart
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(multipart->comm, &i_rank);
  PDM_MPI_Comm_size(multipart->comm, &n_rank);


  /*
   * Step 1 : Split the graph (If we have a dmesh_nodal and prepare the dmesh before the treatment for faces and elements are the same)
   * Step 2 : Rebuild all connectivity in coherent manner
   * Step 3 : Apply all reordering
   * Step 4 : Deduce all mesh_nodal connectivity
   */

  if (multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  } else {
    PDM_timer_t *timer = PDM_timer_create();
    double cum_elapsed_time = 0;
    int *starting_part_idx =  PDM_array_new_idx_from_sizes_int(multipart->n_part, multipart->n_domain);

    // int is_by_elt = 0;
    for (int i_domain = 0; i_domain < multipart->n_domain; ++i_domain) {
      PDM_dmesh_nodal_t* dmesh_nodal = multipart->dmeshes_nodal[i_domain];
      if (dmesh_nodal != NULL) { // element representation
        // is_by_elt = 1;
        // PDM_printf("Partitionning elt domain %d/%d \n", i_domain+1, multipart->n_domain);
        PDM_MPI_Comm comm = multipart->comm;
        PDM_split_dual_t split_method = multipart->split_method;
        int n_part = multipart->n_part[i_domain];
        _part_mesh_t* pmesh = &(multipart->pmeshes[i_domain]);

        // _run_ppart_domain_nodal(dmesh_nodal,pmesh,split_method,n_part,comm);

        const double* part_fraction      = &multipart->part_fraction[starting_part_idx[i_domain]];
        PDM_part_size_t part_size_method = multipart->part_size_method;

        //Convert dmesh nodal to dmesh
        PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
        PDM_dmesh_nodal_generate_distribution(dmesh_nodal);
        PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmesh_nodal);

        // printf("dmesh_nodal->n_cell_abs = "PDM_FMT_G_NUM" \n", dmesh_nodal->n_cell_abs );
        if(dmesh_nodal->n_cell_abs != 0) {
          PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
        } else if (dmesh_nodal->n_face_abs != 0) {
          PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
        } else {
          PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_VTX,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX);
        }

        PDM_dmesh_t  *_dmesh = NULL;
        PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &_dmesh);
        _run_ppart_domain(_dmesh, dmesh_nodal, pmesh, n_part, split_method, part_size_method, part_fraction, comm);
        multipart->dmeshes  [i_domain] = _dmesh;
        multipart->dmn_to_dm[i_domain] = dmn_to_dm; /* Store it - We need it for PDM_multipart_get_part_mesh_nodal */
        // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);

      } else { // face representation
        // PDM_printf("Partitionning face domain %d/%d \n", i_domain+1, multipart->n_domain);

        PDM_MPI_Comm comm = multipart->comm;

        PDM_split_dual_t split_method    = multipart->split_method;
        PDM_part_size_t part_size_method = multipart->part_size_method;

        const double* part_fraction = &multipart->part_fraction[starting_part_idx[i_domain]];

        PDM_dmesh_t  *_dmeshes =   multipart->dmeshes[i_domain];
        _part_mesh_t *_pmeshes = &(multipart->pmeshes[i_domain]);

        int n_part = multipart->n_part[i_domain];


        if (0 && i_rank == 0)
          PDM_printf("Running partitioning for block %i...\n", i_domain+1);
        PDM_timer_resume(timer);
        _run_ppart_domain(_dmeshes, NULL, _pmeshes, n_part, split_method, part_size_method, part_fraction, comm);
        PDM_timer_hang_on(timer);
        if (0 && i_rank == 0)
          PDM_printf("...completed (elapsed time : %f)\n", PDM_timer_elapsed(timer) - cum_elapsed_time);
        cum_elapsed_time = PDM_timer_elapsed(timer);
      }
    }
    PDM_timer_free(timer);

    free(starting_part_idx);
  }
}

/**
 * \brief Retreive the partitionned nodal mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Domain identifier
 * \param [out] pmesh_nodal           Nodal partitionned mesh
 * \param [in]  ownership             Who is responsible to free retreived data ?
 *
 */

void
PDM_multipart_get_part_mesh_nodal
(
PDM_multipart_t        *multipart,
const int               i_domain,
PDM_part_mesh_nodal_t **pmesh_nodal,
PDM_ownership_t         ownership
)
{
  assert(i_domain < multipart->n_domain);

  _part_mesh_t      *pmesh       = &(multipart->pmeshes    [i_domain]);
  PDM_dmesh_nodal_t *dmesh_nodal = multipart->dmeshes_nodal[i_domain];
  if (dmesh_nodal == NULL) {
    *pmesh_nodal = NULL;
  }
  else {
    int n_part = multipart->n_part[i_domain];
    if(dmesh_nodal->mesh_dimension == 3){
      *pmesh_nodal = _compute_part_mesh_nodal_3d(dmesh_nodal, pmesh, n_part, ownership);
    } else if(dmesh_nodal->mesh_dimension == 2){
      *pmesh_nodal = _compute_part_mesh_nodal_2d(dmesh_nodal, pmesh, n_part, ownership);
    } else if(dmesh_nodal->mesh_dimension == 1){
      *pmesh_nodal = _compute_part_mesh_nodal_1d(dmesh_nodal, pmesh, n_part, ownership);
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_compute_part_mesh_nodal error : Bad dmesh_nodal dimension \n");
    }
  }
}

/**
 * \brief Retreive the partitionned mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [out] pmesh                 Partitionned mesh
 *
 */


// void
// PDM_multipart_get_part_mesh
// (
//        PDM_multipart_t  *multipart,
//  const int               i_domain,
//        PDM_part_mesh_t **pmesh
// )
// {
//   assert(i_domain < multipart->n_domain);

//   *pmesh = &(multipart->pmeshes    [i_domain]);
// }

/**
 *
 * \brief Returns the dimensions of a given partition
 *
 * \param [in]
 *
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
)
{

  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  *n_cell = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL);
  *n_face = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE);
  *n_vtx  = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VTX);

  PDM_MPI_Comm_size(multipart->comm, n_proc);
  *n_total_part = _pmeshes.tn_part;

  *s_cell_face = 1;

  int *cell_face     = NULL;
  int *cell_face_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 &cell_face,
                                 &cell_face_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  int *face_vtx     = NULL;
  int *face_vtx_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 &face_vtx,
                                 &face_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);
  if(*n_cell > 0) {
    *s_cell_face = cell_face_idx[*n_cell];
  }
  if(face_vtx_idx != NULL) {
    *s_face_vtx  = face_vtx_idx[*n_face];
  } else {
    *s_face_vtx  = 0;
  }

  int                     *face_part_bound_proc_idx;
  int                     *face_part_bound_part_idx;
  int                     *face_part_bound;
  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_BOUND_TYPE_FACE,
                                    &face_part_bound_proc_idx,
                                    &face_part_bound_part_idx,
                                    &face_part_bound,
                                    PDM_OWNERSHIP_BAD_VALUE);


  *n_face_part_bound = 0;
  if(face_part_bound_part_idx != NULL) {
    *n_face_part_bound = face_part_bound_part_idx[*n_total_part];
  }
  *n_bound_groups = PDM_part_mesh_n_bound_get(_pmeshes.pmesh, PDM_BOUND_TYPE_FACE);

  int         *face_bound_idx;
  int         *face_bound;
  PDM_g_num_t *face_bound_ln_to_gn;
  PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_BOUND_TYPE_FACE,
                                 &face_bound_idx,
                                 &face_bound,
                                 &face_bound_ln_to_gn,
                                 PDM_OWNERSHIP_BAD_VALUE);

  *s_face_bound = 0;
  if(face_bound_idx !=NULL) {
    *s_face_bound   = face_bound_idx[*n_bound_groups];
  }
}


/**
 *
 * \brief Returns the connexion graph between partition for the request \ref PDM_bound_type_t
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_domain              Id of part
 * \param [in]  entity_type           Type of mesh entity
 * \param [out] ppart_bound_proc_idx  Partitioning boundary entities block distribution from processus (size = n_proc + 1)
 * \param [out] ppart_bound_part_idx  Partitioning boundary entities block distribution from partition (size = n_total_part + 1)
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
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  PDM_bound_type_t bound_type = _entity_type_to_bound_type(entity_type);

  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    bound_type,
                                    ppart_bound_proc_idx,
                                    ppart_bound_part_idx,
                                    ppart_bound,
                                    ownership);
}

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
)
{

  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  // *cell_tag = NULL;
  // *face_tag = NULL;
  // *vtx_tag  = NULL;

  // *cell_ln_to_gn = _pmeshes.parts[i_part]->cell_ln_to_gn;
  // *face_ln_to_gn = _pmeshes.parts[i_part]->face_ln_to_gn;
  // *vtx_ln_to_gn  = _pmeshes.parts[i_part]->vtx_ln_to_gn;


  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    cell_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    vtx_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  // *cell_face_idx = _pmeshes.parts[i_part]->cell_face_idx;
  // *cell_face     = _pmeshes.parts[i_part]->cell_face;
  // *face_cell     = _pmeshes.parts[i_part]->face_cell;
  // *face_vtx_idx  = _pmeshes.parts[i_part]->face_vtx_idx;
  // *face_vtx      = _pmeshes.parts[i_part]->face_vtx;

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 cell_face,
                                 cell_face_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 face_vtx,
                                 face_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  int *face_cell_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                 face_cell,
                                 &face_cell_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);
  assert(face_cell_idx == NULL);

  PDM_part_mesh_vtx_coord_get(_pmeshes.pmesh, i_part, vtx, PDM_OWNERSHIP_BAD_VALUE);

  // *face_part_bound_proc_idx = _pmeshes.parts[i_part]->face_part_bound_proc_idx;
  // *face_part_bound_part_idx = _pmeshes.parts[i_part]->face_part_bound_part_idx;
  // *face_part_bound          = _pmeshes.parts[i_part]->face_part_bound;

  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_BOUND_TYPE_FACE,
                                    face_part_bound_proc_idx,
                                    face_part_bound_part_idx,
                                    face_part_bound,
                                    PDM_OWNERSHIP_BAD_VALUE);

  int n_cell = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
  if (n_cell > 0) {
    // *face_bound_idx       = _pmeshes.parts[i_part]->face_bound_idx;
    // *face_bound           = _pmeshes.parts[i_part]->face_bound;
    // *face_bound_ln_to_gn  = _pmeshes.parts[i_part]->face_bound_ln_to_gn;
    PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                   i_part,
                                   PDM_BOUND_TYPE_FACE,
                                   face_bound_idx,
                                   face_bound,
                                   face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_BAD_VALUE);
  } else {
    // *face_bound_idx       = _pmeshes.parts[i_part]->edge_bound_idx;
    // *face_bound           = _pmeshes.parts[i_part]->edge_bound;
    // *face_bound_ln_to_gn  = _pmeshes.parts[i_part]->edge_bound_ln_to_gn;
    PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                   i_part,
                                   PDM_BOUND_TYPE_EDGE,
                                   face_bound_idx,
                                   face_bound,
                                   face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_BAD_VALUE);
  }
  // *face_join_idx        = NULL; // _pmeshes.parts[i_part]->face_join_idx;
  // *face_join            = NULL; // _pmeshes.parts[i_part]->face_join;
  // *face_join_ln_to_gn   = NULL; // _pmeshes.parts[i_part]->face_join_ln_to_gn;

  // *elt_vtx_idx          = NULL; // _pmeshes.parts[i_part]->elt_vtx_idx;
  // *elt_vtx              = NULL; // _pmeshes.parts[i_part]->elt_vtx;
  // *elt_section_ln_to_gn = NULL; // _pmeshes.parts[i_part]->elt_section_ln_to_gn;
}


/**
 *
 * \brief Returns the total number of part among all process
 */
int
PDM_multipart_part_tn_part_get
(
PDM_multipart_t *multipart,
const int        i_domain
)
{
  assert(i_domain < multipart->n_domain);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];
  return PDM_part_mesh_tn_part_get(_pmeshes.pmesh);
}

/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  connectivity_type     Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect_idx           Connectivity index (size = n_entity+1 )
 * \param [in]  connect               Connectivity array (size = connect_idx[n_entity] )
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
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];
  int pn_entity = -1;

  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_EDGE  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VTX);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_part_connectivity_get error : Wrong connectivity_type \n");
  }

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 connectivity_type,
                                 connect,
                                 connect_idx,
                                 ownership);

  /* Build the requested connectivity if missing */
  if (*connect == NULL) {
    if (connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX) {
      int *face_edge_idx = NULL;
      int *face_edge     = NULL;
      PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                     &face_edge,
                                     &face_edge_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      if (face_edge_idx != NULL) {
        /* In some case, face_edge connectivity also does not exists (for example 1d meshes)
        so we can not always reconstruct face_vtx connectivity (which should not be
        asked at all, but we preserve old behaviour) */

        assert(face_edge_idx != NULL);
        assert(face_edge     != NULL);

        int *edge_vtx_idx = NULL;
        int *edge_vtx     = NULL;
        PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_BAD_VALUE);

        assert(edge_vtx != NULL);

        PDM_compute_face_vtx_from_face_and_edge(pn_entity,
                                                face_edge_idx,
                                                face_edge,
                                                edge_vtx,
                                                connect);

        // same index as face_edge, do we need a copy?
        *connect_idx = malloc(sizeof(int) * (pn_entity + 1));
        memcpy(*connect_idx, face_edge_idx, sizeof(int) * (pn_entity + 1));
      }
    }
    else {
      // TODO: face_edge/edge_vtx ?
    }

    // Store the connectivity we just built, with the requested ownership
    PDM_part_mesh_connectivity_set(_pmeshes.pmesh,
                                   i_part,
                                   connectivity_type,
                                   *connect,
                                   *connect_idx,
                                   ownership);
  }

  return pn_entity;
}

/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t
 */
int
PDM_multipart_part_n_entity_get
(
PDM_multipart_t            *multipart,
const int                   i_domain,
const int                   i_part,
      PDM_mesh_entities_t   entity_type
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  int pn_entity = 0;
  switch (entity_type) {
    case PDM_MESH_ENTITY_CELL:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL);
      break;
    case PDM_MESH_ENTITY_FACE:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE);
      break;
    case PDM_MESH_ENTITY_EDGE:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_EDGE);
      break;
    case PDM_MESH_ENTITY_VTX:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VTX);
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_part_n_entity_get error : Wrong entity_type \n");
      break;
  }

  return pn_entity;
}

/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_ln_to_gn       Entity local numbering to global numbering (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership             Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
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
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  int pn_entity = PDM_multipart_part_n_entity_get(multipart, i_domain, i_part, entity_type);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    entity_type,
                                    entity_ln_to_gn,
                                    ownership);

  return pn_entity;
}

/**
 *
 * \brief Return number of entity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_color          Entity color (only for specific renumbering option )
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
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
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  int pn_entity = PDM_multipart_part_n_entity_get(multipart, i_domain, i_part, entity_type);
  PDM_part_mesh_entity_color_get(_pmeshes.pmesh,
                                 i_part,
                                 entity_type,
                                 entity_color,
                                 ownership);

  return pn_entity;
}

/**
 *
 * \brief Get array containing hyperplane color
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  hyperplane_color      Hyperplane color
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_hyperplane_color_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **hyperplane_color,
      PDM_ownership_t   ownership
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_domain]);

  *hyperplane_color = _pmeshes->hyperplane_color[i_part];
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_domain].is_owner_hyperplane_color = PDM_FALSE;
  } else {
    multipart->pmeshes[i_domain].is_owner_hyperplane_color = PDM_TRUE;
  }
}

/**
 *
 * \brief Get array containing thread color - Only if specific reordering (in paradigma plugins)
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  thread_color          Thread color
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_thread_color_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **thread_color,
      PDM_ownership_t   ownership
)
{
  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_domain]);

  *thread_color = _pmeshes->thread_color[i_part];
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_domain].is_owner_thread_color = PDM_FALSE;
  } else {
    multipart->pmeshes[i_domain].is_owner_thread_color = PDM_TRUE;
  }
}

/**
 *
 * \brief Get array containing vtx_ghost_information, usefull to have a priority on vertex between 2 partitions
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_domain              Id of domain
 * \param [in]  i_part                Id of part
 * \param [in]  vtx_ghost_information Integer that give the current priority of vertices on current partitions
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_ghost_infomation_get
(
PDM_multipart_t        *multipart,
const int               i_domain,
const int               i_part,
      int             **vtx_ghost_information,
      PDM_ownership_t   ownership
)
{

  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_domain]);

  *vtx_ghost_information = _pmeshes->vtx_ghost_information[i_part];
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_domain].is_owner_vtx_ghost_information = PDM_FALSE;
  } else {
    multipart->pmeshes[i_domain].is_owner_vtx_ghost_information = PDM_TRUE;
  }
}


/**
 *
 * \brief Return times for a given domain
 * (NOT IMPLEMENTED)
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
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
)
{
  assert(i_domain < multipart->n_domain);

  // PDM_printf("PDM_multipart_time_get: Not implemented\n");
  *elapsed  = NULL;
  *cpu      = NULL;
  *cpu_user = NULL;
  *cpu_sys  = NULL;

}

/**
 *
 * \brief Free the structure
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 */

void
PDM_multipart_free
(
 PDM_multipart_t *multipart
)
{
  // free(multipart->dmeshes_ids);
  for (int i_domain = 0; i_domain < multipart->n_domain; i_domain++) {
    if (multipart->pmeshes[i_domain].joins_ids != NULL) {
      free(multipart->pmeshes[i_domain].joins_ids);
    }

    for (int i_part = 0; i_part < multipart->n_part[i_domain]; i_part++) {
      if(multipart->pmeshes[i_domain].vtx_ghost_information[i_part] != NULL) {
        if(multipart->pmeshes[i_domain].is_owner_vtx_ghost_information == PDM_TRUE) {
          free(multipart->pmeshes[i_domain].vtx_ghost_information[i_part]);
        }
      }

      if(multipart->pmeshes[i_domain].hyperplane_color[i_part] != NULL) {
        if(multipart->pmeshes[i_domain].is_owner_hyperplane_color == PDM_TRUE) {
          free(multipart->pmeshes[i_domain].hyperplane_color[i_part]);
        }
      }

      if(multipart->pmeshes[i_domain].thread_color[i_part] != NULL) {
        if(multipart->pmeshes[i_domain].is_owner_thread_color == PDM_TRUE) {
          free(multipart->pmeshes[i_domain].thread_color[i_part]);
        }
      }
    }
    free(multipart->pmeshes[i_domain].vtx_ghost_information);
    free(multipart->pmeshes[i_domain].hyperplane_color);
    free(multipart->pmeshes[i_domain].thread_color);

    PDM_part_mesh_free(multipart->pmeshes[i_domain].pmesh);

    if(multipart->dmeshes[i_domain] != NULL ) {
      if(multipart->is_owner_dmeshes[i_domain] == PDM_TRUE) {
        PDM_dmesh_free(multipart->dmeshes[i_domain]);
      }
    }

    if(multipart->dmn_to_dm[i_domain] != NULL) {
      PDM_dmesh_nodal_to_dmesh_free(multipart->dmn_to_dm[i_domain]);
      multipart->dmn_to_dm[i_domain] = NULL;
    }
  }
  free(multipart->pmeshes);
  free(multipart->dmeshes);
  free(multipart->dmeshes_nodal);
  free(multipart->dmn_to_dm);
  free(multipart->is_owner_dmeshes);
  free(multipart->n_part);

  //PDM_part_renum_method_purge();
  free (multipart);
  multipart = NULL;

  // PDM_printf("Cleaned from PDM_multipart_free\n");
}


/**
 *
 * \brief Get the vertex coordinates on current i_domain, i_part partition and return number of vertices
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_domain       Id of current domain
 * \param [in]   i_part         Id of part
 * \param [out]  vtx_coord      Vertex coordinate (size = 3 * n_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
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
)
{
  PDM_UNUSED(ownership);

  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  PDM_part_mesh_vtx_coord_get(_pmeshes.pmesh,
                              i_part,
                              vtx_coord,
                              ownership);

  return PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VTX);
}


/**
 *
 * \brief Get the group description for a given entity
 *
 * \param [in]   multipart              Pointer to \ref PDM_multipart_t object
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
 int                  *n_group,
 int                 **group_entity_idx,
 int                 **group_entity,
 PDM_g_num_t         **group_entity_ln_to_gn,
 PDM_ownership_t       ownership
)
{

  assert(i_domain < multipart->n_domain && i_part < multipart->n_part[i_domain]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  PDM_bound_type_t bound_type = _entity_type_to_bound_type(entity_type);

  *n_group = PDM_part_mesh_n_bound_get(_pmeshes.pmesh, bound_type);

  PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                 i_part,
                                 bound_type,
                                 group_entity_idx,
                                 group_entity,
                                 group_entity_ln_to_gn,
                                 ownership);
}


/**
 *
 * \brief Return statistics
 *
 * \param [in]   ppart                          Pointer to \ref PDM_part object
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
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
)
{
  int n_rank;
  PDM_MPI_Comm_size(multipart->comm, &n_rank);

  assert(i_domain < multipart->n_domain);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_domain];

  int* dpart_proc = (int *) malloc((n_rank + 1) * sizeof(int));
  PDM_MPI_Allgather((void *) &multipart->n_part[i_domain],
                    1,
                    PDM_MPI_INT,
           (void *) (&dpart_proc[1]),
                    1,
                    PDM_MPI_INT,
                    multipart->comm);

  dpart_proc[0] = 0;
  for (int i = 1; i < n_rank+1; i++) {
    dpart_proc[i] = dpart_proc[i] + dpart_proc[i-1];
  }

  int *n_loc = (int *) malloc(multipart->n_part[i_domain]  * sizeof(int));
  int *n_tot = (int *) malloc(dpart_proc[n_rank]          * sizeof(int));

  int *s_loc = (int *) malloc(multipart->n_part[i_domain]  * sizeof(int));
  int *s_tot = (int *) malloc(dpart_proc[n_rank]          * sizeof(int));

  for (int i = 0; i < multipart->n_part[i_domain]; i++) {
    n_loc[i] = 0;
    s_loc[i] = 0;
  }

  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    n_tot[i] = 0;
    s_tot[i] = 0;
  }


  int tn_part = dpart_proc[n_rank];
  for (int i_part = 0; i_part < multipart->n_part[i_domain]; i_part++) {
    n_loc[i_part] = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );

    int                     *face_part_bound_proc_idx;
    int                     *face_part_bound_part_idx;
    int                     *face_part_bound;
    PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                      i_part,
                                      PDM_BOUND_TYPE_FACE,
                                      &face_part_bound_proc_idx,
                                      &face_part_bound_part_idx,
                                      &face_part_bound,
                                      PDM_OWNERSHIP_BAD_VALUE);
    if(face_part_bound_part_idx != NULL) {
      s_loc[i_part] = face_part_bound_part_idx[tn_part];
    }
  }

  int *n_part_proc = (int *) malloc((n_rank) * sizeof(int));

  for (int i = 0; i < n_rank; i++) {
    n_part_proc[i] = dpart_proc[i+1] - dpart_proc[i];
  }

  PDM_MPI_Allgatherv(n_loc,
                     multipart->n_part[i_domain],
                     PDM_MPI_INT,
                     n_tot,
                     n_part_proc,
                     dpart_proc,
                     PDM_MPI_INT,
                     multipart->comm);

  PDM_MPI_Allgatherv(s_loc,
                     multipart->n_part[i_domain],
                     PDM_MPI_INT,
                     s_tot,
                     n_part_proc,
                     dpart_proc,
                     PDM_MPI_INT,
                     multipart->comm);

  PDM_quick_sort_int(s_tot, 0, dpart_proc[n_rank]-1);
  PDM_quick_sort_int(n_tot, 0, dpart_proc[n_rank]-1);

  double   _cells_average;
  double   _bound_part_faces_average;

  *bound_part_faces_min = -1;
  *bound_part_faces_max = -1;
  *cells_min = -1;
  *cells_max = -1;
  _cells_average = 0;
  _bound_part_faces_average = 0;

  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    if (*bound_part_faces_min < 0)
      *bound_part_faces_min = s_tot[i];
    else
      *bound_part_faces_min = PDM_MIN(*bound_part_faces_min, s_tot[i]);
    if (*bound_part_faces_max < 0)
      *bound_part_faces_max = s_tot[i];
    else
      *bound_part_faces_max = PDM_MAX(*bound_part_faces_max, s_tot[i]);
    if (*cells_min < 0)
      *cells_min = n_tot[i];
    else
      *cells_min = PDM_MIN(*cells_min, n_tot[i]);
    if (*cells_max < 0)
      *cells_max = n_tot[i];
    else
      *cells_max = PDM_MAX(*cells_max, n_tot[i]);

    _cells_average += n_tot[i];
    _bound_part_faces_average += s_tot[i];
  }

  _cells_average = (_cells_average/((double) dpart_proc[n_rank]));
  *bound_part_faces_sum = (int) _bound_part_faces_average;
  _bound_part_faces_average =
    _bound_part_faces_average/((double) dpart_proc[n_rank]);

  *cells_average = (int) round(_cells_average);
  *bound_part_faces_average = (int) round(_bound_part_faces_average);

  *cells_std_deviation = 0.;
  *bound_part_faces_std_deviation = 0.;
  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    *cells_std_deviation += (n_tot[i] - _cells_average) * (n_tot[i] - _cells_average);
    *bound_part_faces_std_deviation += (s_tot[i] - _bound_part_faces_average) *
                                      (s_tot[i] - _bound_part_faces_average);
  }

  *cells_std_deviation = sqrt(*cells_std_deviation/dpart_proc[n_rank]);
  *bound_part_faces_std_deviation =
    sqrt(*bound_part_faces_std_deviation/dpart_proc[n_rank]);

  int mid = dpart_proc[n_rank]/2;
  if (dpart_proc[n_rank] % 2 == 1) {
    *cells_median = n_tot[mid];
    *bound_part_faces_median = s_tot[mid];
  }

  else {
    *cells_median =(int) round((n_tot[mid-1] + n_tot[mid])/2.);
    *bound_part_faces_median = (int) ((s_tot[mid-1] + s_tot[mid])/2.);
  }

  free(n_part_proc);
  free(n_tot);
  free(s_tot);
  free(n_loc);
  free(s_loc);
  free(dpart_proc);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
