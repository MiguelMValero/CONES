
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
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_dmesh_to_dmesh_nodal.h"
#include "pdm_dmesh_to_dmesh_nodal_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_sort.h"
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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
int
_generate_sections
(
 PDM_MPI_Comm           comm,
 PDM_g_num_t           *distrib_entity,
 PDM_g_num_t           *section_n,
 PDM_Mesh_nodal_elt_t  *section_kind,
 int                    n_section_tot,
 PDM_g_num_t          **out_post_section_n,
 int                  **out_post_section_kind,
 int                  **out_local_post_section_n,
 int                  **out_local_post_section_idx
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Rebuild sections
   */
  int *section_kind_n = malloc(n_rank * sizeof(int));
  PDM_MPI_Allgather(&n_section_tot, 1, PDM_MPI_INT, section_kind_n, 1, PDM_MPI_INT, comm);

  int *section_kind_idx = malloc((n_rank+1) * sizeof(int));
  section_kind_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    section_kind_idx[i+1] = section_kind_idx[i] + section_kind_n[i];
  }

  int         *g_section_kind = (int         *) malloc(section_kind_idx[n_rank] * sizeof(int        ));
  PDM_g_num_t *g_section_n    = (PDM_g_num_t *) malloc(section_kind_idx[n_rank] * sizeof(PDM_g_num_t));

  PDM_MPI_Allgatherv(section_kind, n_section_tot, PDM_MPI_INT,
                     g_section_kind, section_kind_n, section_kind_idx, PDM_MPI_INT, comm);
  PDM_MPI_Allgatherv(section_n, n_section_tot, PDM__PDM_MPI_G_NUM,
                     g_section_n, section_kind_n, section_kind_idx, PDM__PDM_MPI_G_NUM, comm);

  if(1 == 0) {
    PDM_log_trace_array_long(section_n, n_section_tot, "section_n ::");
    PDM_log_trace_array_long(distrib_entity, n_rank+1, "distrib_entity ::");
    PDM_log_trace_connectivity_int (section_kind_idx, g_section_kind, n_rank, "g_section_kind : ");
    PDM_log_trace_connectivity_long(section_kind_idx, g_section_n   , n_rank, "g_section_n      : ");
  }

  /*
   * Post-treament
   */
  int n_section_all_rank = section_kind_idx[n_rank];
  int         *post_section_kind = (int         *) malloc(section_kind_idx[n_rank] * sizeof(int        ));
  PDM_g_num_t *post_section_n    = (PDM_g_num_t *) malloc(section_kind_idx[n_rank] * sizeof(PDM_g_num_t));

  int n_section_post = 0;
  PDM_Mesh_nodal_elt_t first = PDM_MESH_NODAL_N_ELEMENT_TYPES;
  for(int i = 0; i < n_section_all_rank; ++i) {
    if(first != (PDM_Mesh_nodal_elt_t) g_section_kind[i]) {
      first = (PDM_Mesh_nodal_elt_t) g_section_kind[i];
      post_section_kind[n_section_post] = first;
      post_section_n   [n_section_post] = g_section_n[i];
      n_section_post++;
    } else {
      post_section_n[n_section_post-1] += g_section_n[i];
    }
  }
  post_section_kind = realloc(post_section_kind, n_section_post * sizeof(int        ));
  post_section_n    = realloc(post_section_n   , n_section_post * sizeof(PDM_g_num_t));

  free(g_section_kind);
  free(g_section_n);
  free(section_kind_idx);
  free(section_kind_n);
  free(section_n);
  free(section_kind);


  if(1 == 0) {
    PDM_log_trace_array_int (post_section_kind, n_section_post, "post_section_kind ::");
    PDM_log_trace_array_long(post_section_n   , n_section_post, "post_section_n    ::");
  }

  /*
   * Compute global shift
   */
  PDM_g_num_t* post_section_idx = malloc((n_section_post+1) * sizeof(PDM_g_num_t));
  post_section_idx[0] = 0;
  for(int i = 0; i < n_section_post; ++i) {
    post_section_idx[i+1] = post_section_idx[i] + post_section_n[i];
  }


  int dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  /*
   * Reanrange in global section the local one
   */
  int* local_post_section_n   = malloc((n_section_post  ) * sizeof(int));
  int* local_post_section_idx = malloc((n_section_post+1) * sizeof(int));
  for(int i = 0; i < n_section_post; ++i) {
    local_post_section_n[i] = 0;
  }

  for(int i_entity = 0; i_entity < dn_entity; ++i_entity) {
    PDM_g_num_t gnum = distrib_entity[i_rank] + i_entity + 1;
    int t_section = PDM_binary_search_gap_long(gnum-1, post_section_idx, n_section_post+1);
    local_post_section_n[t_section]++;
  }

  local_post_section_idx[0] = 0;
  for(int i = 0; i < n_section_post; ++i) {
    local_post_section_idx[i+1] = local_post_section_idx[i] + local_post_section_n[i];
  }

  free(post_section_idx );

  *out_post_section_n    = post_section_n;
  *out_post_section_kind = post_section_kind;

  *out_local_post_section_n   = local_post_section_n;
  *out_local_post_section_idx = local_post_section_idx;

  return n_section_post;
}

static
void
_rebuild_dmesh_nodal_2d
(
  PDM_dmesh_nodal_t  *dmn,
  PDM_MPI_Comm        comm,
  PDM_g_num_t        *distrib_face,
  PDM_g_num_t        *distrib_edge,
  PDM_g_num_t        *distrib_vtx,
  PDM_g_num_t        *dface_edge,
  int                *dface_edge_idx,
  PDM_g_num_t        *dface_vtx,
  int                *dface_vtx_idx,
  PDM_g_num_t        *dedge_vtx,
  int                *n_bound,
  int               **dbound_idx,
  PDM_g_num_t       **dbound,
  int                *n_blk_gnum,
  PDM_g_num_t       **blk_entity_gnum,
  PDM_g_num_t       **blk_elmt_gnum
)
{
  PDM_UNUSED(distrib_vtx);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Create implicit partitionning */
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];
  PDM_g_num_t* pface_ln_to_gn = malloc(dn_face * sizeof(PDM_g_num_t));
  for(int i_face = 0; i_face < dn_face; ++i_face) {
    pface_ln_to_gn[i_face] = distrib_face[i_rank] + i_face + 1;
  }

  int dn_edge = distrib_edge[i_rank+1] - distrib_edge[i_rank];
  int* dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
  for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
    dedge_vtx_idx[i_edge] = 2 * i_edge;
  }

  int pn_vtx = 0;
  int         *pface_vtx_idx = NULL;
  int         *pface_vtx     = NULL;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;

  /* Translate */
  if(dface_vtx == NULL){
    int pn_edge = 0;
    int         *pface_edge_idx = NULL;
    int         *pface_edge     = NULL;
    PDM_g_num_t *pedge_ln_to_gn = NULL;

    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                             distrib_face,
                                                             dface_edge_idx,
                                                             dface_edge,
                                                             dn_face,
                                                             pface_ln_to_gn,
                                                             &pn_edge,
                                                             &pedge_ln_to_gn,
                                                             &pface_edge_idx,
                                                             &pface_edge);

    pn_vtx = 0;
    int         *pedge_vtx_idx = NULL;
    int         *pedge_vtx     = NULL;

    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                             distrib_edge,
                                                             dedge_vtx_idx,
                                                             dedge_vtx,
                                                             pn_edge,
                                                             pedge_ln_to_gn,
                                                             &pn_vtx,
                                                             &pvtx_ln_to_gn,
                                                             &pedge_vtx_idx,
                                                             &pedge_vtx);

    PDM_compute_face_vtx_from_face_and_edge(dn_face,
                                            pface_edge_idx,
                                            pface_edge,
                                            pedge_vtx,
                                            &pface_vtx);

    // PDM_log_trace_connectivity_int(pedge_vtx_idx, pedge_vtx, pn_edge, "pedge_vtx ::");

    pface_vtx_idx = pface_edge_idx;

    free(pface_edge    );
    free(pedge_ln_to_gn);
    free(pedge_vtx_idx);
    free(pedge_vtx    );
  } else {
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                             distrib_face,
                                                             dface_vtx_idx,
                                                             dface_vtx,
                                                             dn_face,
                                                             pface_ln_to_gn,
                                                             &pn_vtx,
                                                             &pvtx_ln_to_gn,
                                                             &pface_vtx_idx,
                                                             &pface_vtx);
  }

  /*
   * Reconstruction of section for surfacic
   */
  PDM_g_num_t          *section_n    = malloc((dn_face+1) * sizeof(PDM_g_num_t         )); // Suralloc
  PDM_Mesh_nodal_elt_t *section_kind = malloc( dn_face    * sizeof(PDM_Mesh_nodal_elt_t)); // Suralloc

  int n_section_tot    = 0;
  int n_section_tri    = 0;
  int n_section_quad   = 0;
  int n_section_poly2d = 0;
  int ln_vtx_old = -1;
  for(int i_face = 0; i_face < dn_face; ++i_face) {
    int ln_vtx = pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face];
    if(ln_vtx_old == ln_vtx){
      section_n[n_section_tot-1]++;
      continue;
    }
    if(ln_vtx == 3) {
      section_kind[n_section_tot] = PDM_MESH_NODAL_TRIA3;
      n_section_tri++;
    } else if(ln_vtx == 4){
      section_kind[n_section_tot] = PDM_MESH_NODAL_QUAD4;
      n_section_quad++;
    } else {
      section_kind[n_section_tot] = PDM_MESH_NODAL_POLY_2D;
      n_section_poly2d++;
    }
    section_n[n_section_tot] = 1;
    n_section_tot++;
    ln_vtx_old = ln_vtx;
  }
  section_n    = realloc(section_n   , (n_section_tot+1) * sizeof(PDM_g_num_t         ));
  section_kind = realloc(section_kind,  n_section_tot    * sizeof(PDM_Mesh_nodal_elt_t));

  PDM_g_num_t *post_section_n         = NULL;
  int         *post_section_kind      = NULL;
  int         *local_post_section_n   = NULL;
  int         *local_post_section_idx = NULL;

  int n_section_post = _generate_sections(comm,
                                          distrib_face,
                                          section_n,
                                          section_kind,
                                          n_section_tot,
                                          &post_section_n,
                                          &post_section_kind,
                                          &local_post_section_n,
                                          &local_post_section_idx);

  /*
   * Requilibrate all block
   */
  for(int i_section = 0; i_section < n_section_post; ++i_section) {

    int beg = local_post_section_idx[i_section];
    int end = local_post_section_idx[i_section+1];
    int nl_elmt = end - beg;

    PDM_g_num_t* distrib_elmt = PDM_compute_uniform_entity_distribution(comm, post_section_n[i_section]);
    PDM_g_num_t* ln_to_gn = malloc(local_post_section_n[i_section] * sizeof(PDM_g_num_t));

    for(int i = 0; i < local_post_section_n[i_section]; ++i) {
      ln_to_gn[i] = distrib_face[i_rank] + local_post_section_idx[i_section] + i + 1;
    }

    PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                     PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                     1.,
                                                                     &ln_to_gn,
                                                                     distrib_elmt,
                                                                     &local_post_section_n[i_section],
                                                                     1,
                                                                     comm);
    int n_face_vtx_tot = pface_vtx_idx[end] - pface_vtx_idx[beg];
    int         *send_face_vtx_n = malloc(nl_elmt        * sizeof(int        ));
    PDM_g_num_t *send_face_vtx   = malloc(n_face_vtx_tot * sizeof(PDM_g_num_t));
    int         *blk_face_vtx_n  = NULL;
    PDM_g_num_t *blk_face_vtx    = NULL;

    int idx_write = 0;
    for(int i = 0; i < nl_elmt; ++i) {
      send_face_vtx_n[i] = pface_vtx_idx[beg+i+1] - pface_vtx_idx[beg+i];
      for(int j = pface_vtx_idx[beg+i]; j < pface_vtx_idx[beg+i+1]; ++j) {
        int i_vtx = pface_vtx[j];
        send_face_vtx[idx_write++] = pvtx_ln_to_gn[i_vtx-1];
      }
    }

    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
             (int  **)     &send_face_vtx_n,
             (void **)     &send_face_vtx,
             (int  **)     &blk_face_vtx_n,
             (void **)     &blk_face_vtx);

    free(send_face_vtx_n);
    free(send_face_vtx);


    PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) post_section_kind[i_section];
    int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                 PDM_GEOMETRY_KIND_SURFACIC,
                                                 t_elt);

    int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
    if(t_elt == PDM_MESH_NODAL_POLY_2D) {

      int *blk_face_vtx_idx = PDM_array_new_idx_from_sizes_int(blk_face_vtx_n, dn_elmt);
      PDM_DMesh_nodal_section_poly2d_set(dmn,
                                         PDM_GEOMETRY_KIND_SURFACIC,
                                         id_section,
                                         dn_elmt,
                                         blk_face_vtx_idx,
                                         blk_face_vtx,
                                         PDM_OWNERSHIP_KEEP);
    } else {
      PDM_DMesh_nodal_section_std_set(dmn,
                                      PDM_GEOMETRY_KIND_SURFACIC,
                                      id_section,
                                      dn_elmt,
                                      blk_face_vtx,
                                      PDM_OWNERSHIP_KEEP);
    }

    free(blk_face_vtx_n);
    // free(blk_face_vtx);

    PDM_part_to_block_free(ptb);
    free(distrib_elmt);
    free(ln_to_gn);

  }

  free(post_section_kind);
  free(post_section_n   );
  free(local_post_section_n  );
  free(local_post_section_idx);

  if(1 == 0) {
    printf("n_section_tri    = %i\n", n_section_tri   );
    printf("n_section_quad   = %i\n", n_section_quad  );
    printf("n_section_poly2d = %i\n", n_section_poly2d);
  }

  /*
   * Recuperation des bords
   */
  int          n_edge_bound   = n_bound   [PDM_BOUND_TYPE_EDGE];
  int         *edge_bound_idx = dbound_idx[PDM_BOUND_TYPE_EDGE];
  PDM_g_num_t *edge_bound     = dbound    [PDM_BOUND_TYPE_EDGE];
  // int n_section_bar = n_edge_bound;

  int n_edge_bnd_tot = edge_bound_idx[n_edge_bound];
  int          pn_vtx_bnd = 0;
  int         *pedge_vtx_bnd_idx = NULL;
  int         *pedge_bnd_vtx     = NULL;
  PDM_g_num_t *pvtx_bnd_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           n_edge_bnd_tot,
                                                           edge_bound,
                                                           &pn_vtx_bnd,
                                                           &pvtx_bnd_ln_to_gn,
                                                           &pedge_vtx_bnd_idx,
                                                           &pedge_bnd_vtx);

  /*
   * On ne veut pas changer la numerotation absolu des bar
   */

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_BAR2;
  int id_section = PDM_DMesh_nodal_section_add(dmn,
                                               PDM_GEOMETRY_KIND_RIDGE,
                                               t_elt);

  PDM_g_num_t *pedge_bnd_vtx_gnum = malloc(2 * n_edge_bnd_tot * sizeof(PDM_g_num_t));
  for(int i = 0; i < 2 * n_edge_bnd_tot; ++i) {
    pedge_bnd_vtx_gnum[i] = pvtx_bnd_ln_to_gn[pedge_bnd_vtx[i]-1];
  }
  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  id_section,
                                  n_edge_bnd_tot,
                                  pedge_bnd_vtx_gnum,
                                  PDM_OWNERSHIP_KEEP);


  PDM_g_num_t* distrib_bar = PDM_compute_entity_distribution(comm, n_edge_bnd_tot);
  int         *out_edge_bound_idx = malloc((n_edge_bound+1) * sizeof(int        ));
  PDM_g_num_t *out_edge_bound     = malloc( n_edge_bnd_tot  * sizeof(PDM_g_num_t));

  out_edge_bound_idx[0] = edge_bound_idx[0];
  for(int i_group = 0; i_group < n_edge_bound; ++i_group) {
    out_edge_bound_idx[i_group+1] = edge_bound_idx[i_group+1];
    for(int idx_edge = edge_bound_idx[i_group]; idx_edge < edge_bound_idx[i_group+1]; ++idx_edge) {
      out_edge_bound[idx_edge] = distrib_bar[i_rank] + idx_edge + 1;
    }
  }

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_edge_bound,
                                         out_edge_bound_idx,
                                         out_edge_bound,
                                         PDM_OWNERSHIP_KEEP);

  // PDM_log_trace_connectivity_long(out_edge_bound_idx, out_edge_bound, n_edge_bound, "out_edge_bound ::");

  /*
   *  Il faut créer un table qui permet d'updater les numero de faces/edges en numero d'element
   *     - dedge_to_elmts_n
   *     - dedge_to_elmts
   *  Donc c'est un ptb partial
   *    -> Si on a une liste de edge_ln_to_gn
   *  PDM_block_to_part_create_from_sparse_block
   *
   */
  PDM_part_to_block_t* ptb_edge = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           &edge_bound,
                                                           NULL,
                                                           &n_edge_bnd_tot,
                                                           1,
                                                           comm);
  PDM_g_num_t *gnum_edge = PDM_part_to_block_block_gnum_get(ptb_edge);
  int n_blk_edge = PDM_part_to_block_n_elt_block_get(ptb_edge);

  PDM_g_num_t *blk_out_edge_bound = NULL;
  PDM_part_to_block_exch(ptb_edge,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
              (void **)  &out_edge_bound,
                         NULL,
              (void **)  &blk_out_edge_bound);


  // PDM_log_trace_array_long(gnum_edge, n_blk_edge, "gnum_edge ::");
  // PDM_log_trace_array_long(blk_out_edge_bound, n_blk_edge, "blk_out_edge_bound ::");

  /* Store it */
  assert(n_blk_gnum     [PDM_BOUND_TYPE_EDGE] == 0);
  assert(blk_entity_gnum[PDM_BOUND_TYPE_EDGE] == NULL);
  assert(blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE] == NULL);

  n_blk_gnum     [PDM_BOUND_TYPE_EDGE] = n_blk_edge;
  blk_entity_gnum[PDM_BOUND_TYPE_EDGE] = malloc(n_blk_edge * sizeof(PDM_g_num_t));
  blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE] = blk_out_edge_bound;

  for(int i = 0; i < n_blk_edge; ++i) {
    blk_elmt_gnum  [PDM_BOUND_TYPE_EDGE][i] = gnum_edge[i];
  }

  PDM_part_to_block_free(ptb_edge);
  /*
   *  Translation face -> elmts
   */
  // PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(gnum_edge,
  //                                                                       n_blk_edge,
  //                                           (const PDM_g_num_t **)      &edge_bound,
  //                                                                       &n_edge_bnd_tot,
  //                                                                       1,
  //                                                                       comm);
  // int stride_one = 1;
  // PDM_g_num_t **tmp_edge_to_elmt = NULL;
  // PDM_block_to_part_exch(btp,
  //                        sizeof(PDM_g_num_t),
  //                        PDM_STRIDE_CST_INTERLACED,
  //                        &stride_one,
  //                        blk_out_edge_bound,
  //                        NULL,
  //             (void ***) &tmp_edge_to_elmt);
  // PDM_g_num_t *edge_to_elmt = tmp_edge_to_elmt[0];
  // free(tmp_edge_to_elmt);

  // PDM_log_trace_array_long(edge_to_elmt, n_edge_bnd_tot, "edge_to_elmt ::");


  // free(edge_to_elmt);
  // // free(blk_out_edge_bound);
  // PDM_block_to_part_free(btp);


  free(distrib_bar);

  free(pedge_vtx_bnd_idx);
  free(pedge_bnd_vtx    );
  free(pvtx_bnd_ln_to_gn);
  free(dedge_vtx_idx);

  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx    );

}

/**
 * \brief Get element type
 *
 *   \param[in] n_cell_face     Number of faces in the current cell
 *   \param[in]  face_cell      Face to cell connectivity
 *   \param[in]  face_cell_idx  Face to vertex connectivity index
 *   \param[in]  face_cell_n    Number of vertices for each face
 *   \param[in]  face_vtx       Face to vertex connectivity
 *   \param[out] tria_vtx       Quadrangles connectivity
 *   \param[out] quad_vtx       Quadrangles connectivity
 *
 *  \return   Cell type
 */

inline static
PDM_Mesh_nodal_elt_t
_type_cell_3D
(
 const int             n_cell_face,
 const PDM_l_num_t    *cell_face,
 const PDM_l_num_t    *face_vtx_idx,
 const PDM_l_num_t    *face_vtx,
 PDM_l_num_t           tria_vtx[],
 PDM_l_num_t           quad_vtx[]
)
{
  PDM_l_num_t  n_trias = 0;
  PDM_l_num_t  n_quads = 0;

  if (n_cell_face > 6) {
    return PDM_MESH_NODAL_POLY_3D;
  }

  for (int i = 0; i < n_cell_face; i++) {

    const int face_id = PDM_ABS(cell_face[i]) - 1;
    const int n_som_face = face_vtx_idx[face_id+1]-face_vtx_idx[face_id];
    PDM_l_num_t idx = face_vtx_idx[face_id];

    if (n_som_face == 3) {
      PDM_l_num_t *cell_som_tria_courant = tria_vtx + 3*n_trias;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_tria_courant[j-idx] = face_vtx[j];
      }
      n_trias += 1;
    }
    else if (n_som_face == 4) {
      PDM_l_num_t *cell_som_quad_courant = quad_vtx + 4*n_quads;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_quad_courant[j-idx] = face_vtx[j];
      }
      n_quads += 1;
    }
    else
      return PDM_MESH_NODAL_POLY_3D;

  }

  PDM_Mesh_nodal_elt_t cell_type;

  if ((n_quads == 0) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_TETRA4;
  else if (n_quads == 6)
    cell_type = PDM_MESH_NODAL_HEXA8;
  else if ((n_quads == 1) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_PYRAMID5;
  else if ((n_quads == 3) && (n_trias == 2)) {
    int trias[6];
    n_trias = 0;
    for (int i = 0; i < n_cell_face; i++) {

      const int face_id = PDM_ABS(cell_face[i]) - 1;
      const int ideb = face_vtx_idx[face_id];

      const int n_som_face = face_vtx_idx[face_id+1]-face_vtx_idx[face_id];

      if (n_som_face == 3) {
        for (int j = 0; j < 3; j++) {
          trias[3*n_trias+j] = face_vtx[ideb+j];
        }
        n_trias += 1;
      }
      if (n_trias >= 2)
        break;
    }

    cell_type = PDM_MESH_NODAL_PRISM6;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (trias[i] == trias[3+j]) {
          cell_type = PDM_MESH_NODAL_POLY_3D;
          break;
        }
      }
      if (cell_type == PDM_MESH_NODAL_POLY_3D)
        break;
    }
  }

  else {
    cell_type = PDM_MESH_NODAL_POLY_3D;
  }

  return cell_type;

}

static inline double
_p_dot
(
 const double a[3],
 const double b[3]
 )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


static inline void
_p_cross
(
 const double a[3],
 const double b[3],
 double c[3]
 )
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = b[0] * a[2] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}


/**
 *
 * \brief Build tetrahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[out] tetra_vtx   Tetrahedron connectivity
 *
 */

static void
_connec_tetra
(
 double      *coords,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t  tetra_vtx[]
 )
{

  /* Initialization */

  tetra_vtx[0] = tria_vtx[0];
  tetra_vtx[1] = tria_vtx[1];
  tetra_vtx[2] = tria_vtx[2];

  for (int i = 3; i < 11; i++) {
    if ((tria_vtx[i] != tetra_vtx[0]) &&
        (tria_vtx[i] != tetra_vtx[1]) &&
        (tria_vtx[i] != tetra_vtx[2]))
      tetra_vtx[3] = tria_vtx[i];
  }

  /* Orientation */
  double v1[3];
  double v2[3];
  double v3[3];
  double n[3];

  for (int i = 0; i < 3; i++) {
    v1[i] = coords[3*(tetra_vtx[1] - 1) + i] - coords[3*(tetra_vtx[0] - 1) + i];
    v2[i] = coords[3*(tetra_vtx[2] - 1) + i] - coords[3*(tetra_vtx[0] - 1) + i];
    v3[i] = coords[3*(tetra_vtx[3] - 1) + i] - coords[3*(tetra_vtx[0] - 1) + i];
  }

  _p_cross(v1, v2, n);
  double orient = _p_dot(v3, n);

  if (orient < 0) {
    tetra_vtx[0] = tria_vtx[2];
    tetra_vtx[1] = tria_vtx[1];
    tetra_vtx[2] = tria_vtx[0];
  }
}

/**
 *
 * \brief Build prism nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] prism_vtx   Prism connectivity
 *
 */

static void
_connec_prism
(
 double      *coords,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t *quad_vtx,
 PDM_l_num_t  prism_vtx[]
 )
{

  /* Initialisation */

  for (int i = 0; i < 6; i++)
    prism_vtx[i] = tria_vtx[i];

  /* Orientation des faces */
  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 3; j++) {
      int isom = prism_vtx[3*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 1.0/3.0;

    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;

    double v1[3];
    double v2[3];
    int isom3 = prism_vtx[3*i+2] - 1 ;
    int isom2 = prism_vtx[3*i+1] - 1;
    int isom1 = prism_vtx[3*i] - 1;

    for (int k = 0; k < 3; k++) {
      v1[k] = coords[3*isom2+k] - coords[3*isom1+k];
      v2[k] = coords[3*isom3+k] - coords[3*isom1+k];
    }
    _p_cross(v1, v2, n + 3*i);
  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = prism_vtx[1];
    prism_vtx[1] = prism_vtx[2];
    prism_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = prism_vtx[4];
    prism_vtx[4] = prism_vtx[5];
    prism_vtx[5] = tmp;
  }

  /* Permutation circulaire */

  int id1 = -1;
  for (int j = 0; j < 12; j++) {
    if (quad_vtx[j] == prism_vtx[0]) {
      id1 = j;
      break;
    }
  }

  int id2 = (id1 / 4) * 4 + (id1 + 1) % 4;
  if ((quad_vtx[id2] == prism_vtx[1]) ||
      (quad_vtx[id2] == prism_vtx[2]))
    id2 =  (id1 / 4) * 4 + (id1 + 3) % 4;

  int id_deb = -1;
  for (int j = 0; j < 3; j++) {
    if (quad_vtx[id2] == prism_vtx[3+j]) {
      id_deb = j;
      break;
    }
  }

  int tmp[3];
  for (int j = 0; j < 3; j++)
    tmp[j] = prism_vtx[3+j];

  for (int j = 0; j < 3; j++) {
    int idx = (id_deb + j) % 3;
    prism_vtx[3+j] = tmp[idx];
  }

}


/**
 *
 * \brief Build pyramid nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] pyramid_vtx Pyramid connectivity
 *
 */

static void
_connec_pyramid
(
 double      *coords,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t *quad_vtx,
 PDM_l_num_t  pyramid_vtx[]
 )
{

  /* Initialisation */

  pyramid_vtx[0] = quad_vtx[0];
  pyramid_vtx[1] = quad_vtx[1];
  pyramid_vtx[2] = quad_vtx[2];
  pyramid_vtx[3] = quad_vtx[3];

  for (int i = 0; i < 9; i++) {
    if ((tria_vtx[i] != pyramid_vtx[0]) &&
        (tria_vtx[i] != pyramid_vtx[1]) &&
        (tria_vtx[i] != pyramid_vtx[2]) &&
        (tria_vtx[i] != pyramid_vtx[3])) {
      pyramid_vtx[4] = tria_vtx[i];
      break;
    }
  }

  /* Orientation */
  double c[3];
  double n[3];

  for (int k = 0; k < 3; k++)
    c[k] = 0.;
  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    for (int k = 0; k < 3; k++)
      c[k] += coords[3*isom+k];
  }
  for (int k = 0; k < 3; k++)
    c[k] *= 0.25;

  for (int k = 0; k < 3; k++)
    n[k] = 0.;

  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    int suiv = (j+1) % 4;
    int isom_suiv = pyramid_vtx[suiv] - 1;

    double v1[3];
    double v2[3];
    for (int k = 0; k < 3; k++) {
      v1[k] = coords[3*isom+k] -  c[k];
      v2[k] = coords[3*isom_suiv+k] -  c[k];
    }

    _p_cross(v1, v2, n);

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = coords[3*(pyramid_vtx[3] - 1) + k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);

  if (orientation < 0) {
    int tmp = pyramid_vtx[0];
    pyramid_vtx[0] = pyramid_vtx[3];
    pyramid_vtx[3] = tmp;
    tmp = pyramid_vtx[1];
    pyramid_vtx[1] = pyramid_vtx[2];
    pyramid_vtx[2] = tmp;
  }

}


/**
 *
 * \brief Build hexahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] hexa_vtx    Hexahedron connectivity
 *
 */

static void
_connec_hexa
(
 double       *coords,
 PDM_l_num_t  *quad_vtx,
 PDM_l_num_t   hexa_vtx[]
)
{

  /* Initialization */

  hexa_vtx[0] = quad_vtx[0];
  hexa_vtx[1] = quad_vtx[1];
  hexa_vtx[2] = quad_vtx[2];
  hexa_vtx[3] = quad_vtx[3];

  PDM_l_num_t face_contact[4] = {-1, -1, -1, -1};

  for (int i = 1; i < 6; i++) {
    int cpt = 0;
    for (int j = 0; j < 4; j++) {
      PDM_l_num_t som_courant = quad_vtx[4*i+j];
      if ((som_courant != hexa_vtx[0]) &&
          (som_courant != hexa_vtx[1]) &&
          (som_courant != hexa_vtx[2]) &&
          (som_courant != hexa_vtx[3]))
        cpt += 1;
    }
    if (cpt == 4) {
      hexa_vtx[4] = quad_vtx[4*i];
      hexa_vtx[5] = quad_vtx[4*i+1];
      hexa_vtx[6] = quad_vtx[4*i+2];
      hexa_vtx[7] = quad_vtx[4*i+3];
    }
    if (cpt == 2) {
      face_contact[0] = quad_vtx[4*i];
      face_contact[1] = quad_vtx[4*i+1];
      face_contact[2] = quad_vtx[4*i+2];
      face_contact[3] = quad_vtx[4*i+3];
    }
  }

  /* Calcul des centres et normales de la base et de la face opposee */
  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 0.25;

    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;

    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      int suiv = (j+1) % 4;
      int isom_suiv = hexa_vtx[4*i+suiv] - 1;

      double v1[3];
      double v2[3];
      for (int k = 0; k < 3; k++) {
        v1[k] = coords[3*isom+k] -  c[3*i+k];
        v2[k] = coords[3*isom_suiv+k] -  c[3*i+k];
      }

      _p_cross(v1, v2, n + 3*i);

    }

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = hexa_vtx[0];
    hexa_vtx[0] = hexa_vtx[3];
    hexa_vtx[3] = tmp;
    tmp = hexa_vtx[1];
    hexa_vtx[1] = hexa_vtx[2];
    hexa_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = hexa_vtx[4];
    hexa_vtx[4] = hexa_vtx[7];
    hexa_vtx[7] = tmp;
    tmp = hexa_vtx[5];
    hexa_vtx[5] = hexa_vtx[6];
    hexa_vtx[6] = tmp;
  }

  /* Permutation circulaire eventuelle de la face sup */

  int id1 = -1;
  int k1 = -1;
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      if (face_contact[j] == hexa_vtx[k]) {
        id1 = j;
        k1 = k;
        break;
      }
      if (id1 != -1)
        break;
    }
  }

  if (k1 == -1) {
    printf("Error connect_hexa : %d %d %d %d %d %d %d %d\n",
               hexa_vtx[0],
               hexa_vtx[1],
               hexa_vtx[2],
               hexa_vtx[3],
               hexa_vtx[4],
               hexa_vtx[5],
               hexa_vtx[6],
               hexa_vtx[7]);

    for (int i10 = 0; i10 < 4; i10++) {
      printf("   face %d : %d %d %d %d\n", i10+1, quad_vtx[4*i10],
                 quad_vtx[4*i10+1],
                 quad_vtx[4*i10+2],
                 quad_vtx[4*i10+3]);
    }
    abort();

  }

  int id2 = (id1 + 1) % 4;
  int k2 = (k1 + 1) % 4;
  int k3 = (k1 + 3) % 4;

  if ((face_contact[id2] == hexa_vtx[k2]) ||
      (face_contact[id2] == hexa_vtx[k3]))
    id2 = (id1 + 3) % 4;

  int id_deb = -1;
  for (int j = 0; j < 4; j++) {
    if (face_contact[id2] == hexa_vtx[4+j]) {
      id_deb = (j - k1);
      if (id_deb < 0)
        id_deb += 4;
      id_deb = id_deb % 4;
      break;
    }
  }

  int tmp[4];
  for (int j = 0; j < 4; j++)
    tmp[j] = hexa_vtx[4+j];

  for (int j = 0; j < 4; j++) {
    int idx = (id_deb + j) % 4;
    hexa_vtx[4+j] = tmp[idx];
  }
}


static
void
_rebuild_dmesh_nodal_by_kind_3d
(
  PDM_dmesh_nodal_t  *dmn,
  PDM_MPI_Comm        comm,
  PDM_g_num_t        *distrib_cell,
  PDM_g_num_t        *distrib_face,
  PDM_g_num_t        *distrib_edge,
  PDM_g_num_t        *distrib_vtx,
  double             *dvtx_coords,
  PDM_g_num_t        *dcell_face,
  int                *dcell_face_idx,
  PDM_g_num_t        *dface_edge,
  int                *dface_edge_idx,
  PDM_g_num_t        *dface_vtx,
  int                *dface_vtx_idx,
  PDM_g_num_t        *dedge_vtx,
  int                *n_bound,
  int               **dbound_idx,
  PDM_g_num_t       **dbound,
  int                *n_blk_gnum,
  PDM_g_num_t       **blk_entity_gnum,
  PDM_g_num_t       **blk_elmt_gnum
)
{
  PDM_UNUSED(dmn);
  PDM_UNUSED(comm);
  PDM_UNUSED(distrib_cell);
  PDM_UNUSED(distrib_face);
  PDM_UNUSED(distrib_edge);
  PDM_UNUSED(distrib_vtx);
  PDM_UNUSED(dcell_face);
  PDM_UNUSED(dcell_face_idx);
  PDM_UNUSED(dface_edge);
  PDM_UNUSED(dface_edge_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dedge_vtx);
  PDM_UNUSED(n_bound);
  PDM_UNUSED(dbound_idx);
  PDM_UNUSED(dbound);
  PDM_UNUSED(n_blk_gnum);
  PDM_UNUSED(blk_entity_gnum);
  PDM_UNUSED(blk_elmt_gnum);


  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_cell = distrib_cell[i_rank+1] - distrib_cell[i_rank];
  PDM_g_num_t* pcell_ln_to_gn = malloc(dn_cell * sizeof(PDM_g_num_t));
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {
    pcell_ln_to_gn[i_cell] = distrib_cell[i_rank] + i_cell + 1;
  }


  int pn_face = 0;
  int         *pcell_face_idx = NULL;
  int         *pcell_face     = NULL;
  PDM_g_num_t *pface_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_cell,
                                                           dcell_face_idx,
                                                           dcell_face,
                                                           dn_cell,
                                                           pcell_ln_to_gn,
                                                           &pn_face,
                                                           &pface_ln_to_gn,
                                                           &pcell_face_idx,
                                                           &pcell_face);

  int pn_vtx = 0;
  int         *pface_vtx_idx = NULL;
  int         *pface_vtx     = NULL;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           pn_face,
                                                           pface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);
  double **tmp_pvtx_coords = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coords,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coords);
  double *pvtx_coords = tmp_pvtx_coords[0];
  free(tmp_pvtx_coords);

  /*
   * Reconstruction of volumic part
   */
  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */

  PDM_g_num_t  *section_n    = malloc(PDM_MESH_NODAL_N_ELEMENT_TYPES * sizeof(PDM_g_num_t         )); // Suralloc
  PDM_g_num_t  *g_section_n  = malloc(PDM_MESH_NODAL_N_ELEMENT_TYPES * sizeof(PDM_g_num_t         )); // Suralloc
  for(int i = 0; i < PDM_MESH_NODAL_N_ELEMENT_TYPES; ++i) {
    section_n  [i] = 0;
    g_section_n[i] = 0;
  }

  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {
    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(pcell_face_idx[i_cell+1]-pcell_face_idx[i_cell],
                                                   pcell_face + pcell_face_idx[i_cell],
                                                   pface_vtx_idx,
                                                   pface_vtx,
                                                   cell_som_tria,
                                                   cell_som_quad);
    section_n[cell_type] += 1;
  }


  int *section_cell_vtx_n   = malloc(PDM_MESH_NODAL_N_ELEMENT_TYPES     * sizeof(int));
  int *section_cell_n       = malloc(PDM_MESH_NODAL_N_ELEMENT_TYPES     * sizeof(int));
  int *section_cell_vtx_idx = malloc((PDM_MESH_NODAL_N_ELEMENT_TYPES+1) * sizeof(int));
  int *section_cell_idx     = malloc((PDM_MESH_NODAL_N_ELEMENT_TYPES+1) * sizeof(int));
  section_cell_vtx_idx[0] = 0;
  section_cell_idx    [0] = 0;
  for(int i_section = 0; i_section < PDM_MESH_NODAL_N_ELEMENT_TYPES; ++i_section) {
    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element((PDM_Mesh_nodal_elt_t) i_section, 1);
    section_cell_vtx_idx[i_section+1] = section_cell_vtx_idx[i_section] + n_vtx_per_elmt * section_n[i_section];
    section_cell_idx    [i_section+1] = section_cell_idx    [i_section] + section_n[i_section];
    section_cell_vtx_n  [i_section] = 0;
    section_cell_n      [i_section] = 0;
  }


  /* Creation des blocks d'envoi */
  int n_tot_cell_vtx = section_cell_vtx_idx[PDM_MESH_NODAL_N_ELEMENT_TYPES];
  PDM_g_num_t *pcell_vtx_gnum = malloc(n_tot_cell_vtx * sizeof(PDM_g_num_t));
  int         *pcell_vtx_n    = malloc(dn_cell        * sizeof(int        ));
  int lconnect[24];
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {

    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(pcell_face_idx[i_cell+1]-pcell_face_idx[i_cell],
                                                   pcell_face + pcell_face_idx[i_cell],
                                                   pface_vtx_idx,
                                                   pface_vtx,
                                                   cell_som_tria,
                                                   cell_som_quad);
    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element((PDM_Mesh_nodal_elt_t) cell_type, 1);

    switch(cell_type) {
      case PDM_MESH_NODAL_TETRA4 :
        _connec_tetra(pvtx_coords,
                      cell_som_tria,
                      lconnect);
        break;
      case PDM_MESH_NODAL_HEXA8 :
        _connec_hexa(pvtx_coords,
                     cell_som_quad,
                     lconnect);
        break;
      case PDM_MESH_NODAL_PRISM6 :
        _connec_prism(pvtx_coords,
                      cell_som_tria,
                      cell_som_quad,
                      lconnect);
        break;
      case PDM_MESH_NODAL_PYRAMID5 :
        _connec_pyramid(pvtx_coords,
                        cell_som_tria,
                        cell_som_quad,
                        lconnect);
        break;
      case PDM_MESH_NODAL_POLY_3D :
        {
          abort();
          break;
        }
      default :
        break;
    }

    pcell_vtx_n[section_cell_idx[cell_type] + section_cell_n[cell_type]++] = n_vtx_per_elmt;
    for(int i = 0; i < n_vtx_per_elmt; ++i) {
      int idx_write = section_cell_vtx_idx[cell_type] + section_cell_vtx_n[cell_type]++;
      pcell_vtx_gnum[idx_write] = pvtx_ln_to_gn[lconnect[i]-1];
    }
  }

  // Gather data
  PDM_MPI_Allreduce(section_n, g_section_n,
                    PDM_MESH_NODAL_N_ELEMENT_TYPES,
                    PDM__PDM_MPI_G_NUM,
                    PDM_MPI_SUM,
                    comm);

  for(int i_section = 0; i_section < PDM_MESH_NODAL_N_ELEMENT_TYPES; ++i_section) {

    if(g_section_n[i_section] == 0) {
      continue;
    }

    PDM_g_num_t* distrib_elmt = PDM_compute_uniform_entity_distribution(comm, section_n[i_section]);
    PDM_g_num_t* ln_to_gn = malloc(section_n[i_section] * sizeof(PDM_g_num_t));

    for(int i = 0; i < section_n[i_section]; ++i) {
      ln_to_gn[i] = distrib_elmt[i_rank] + i + 1;
    }

    int n_elt = (int) section_n[i_section];

    PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                     PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                     1.,
                                                                     &ln_to_gn,
                                                                     distrib_elmt,
                                                                     &n_elt,//&section_n[i_section],
                                                                     1,
                                                                     comm);

    int         *lcell_vtx_n    = &pcell_vtx_n   [section_cell_idx    [i_section]];
    PDM_g_num_t *lcell_vtx_gnum = &pcell_vtx_gnum[section_cell_vtx_idx[i_section]];

    int         *blk_cell_vtx_n = NULL;
    PDM_g_num_t *blk_cell_vtx   = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
             (int  **)     &lcell_vtx_n,
             (void **)     &lcell_vtx_gnum,
             (int  **)     &blk_cell_vtx_n,
             (void **)     &blk_cell_vtx);


    PDM_part_to_block_free(ptb);

    PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) i_section;
    int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                 PDM_GEOMETRY_KIND_VOLUMIC,
                                                 t_elt);

    int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
    if(t_elt == PDM_MESH_NODAL_POLY_3D) {
      abort();
      // int *blk_face_vtx_idx = PDM_array_new_idx_from_sizes_int(blk_face_vtx_n, dn_elmt);
      // PDM_DMesh_nodal_section_poly2d_set(dmn,
      //                                    PDM_GEOMETRY_KIND_SURFACIC,
      //                                    id_section,
      //                                    dn_elmt,
      //                                    blk_face_vtx_idx,
      //                                    blk_face_vtx,
      //                                    PDM_OWNERSHIP_KEEP);
    } else {
      PDM_DMesh_nodal_section_std_set(dmn,
                                      PDM_GEOMETRY_KIND_VOLUMIC,
                                      id_section,
                                      dn_elmt,
                                      blk_cell_vtx,
                                      PDM_OWNERSHIP_KEEP);
    }

    free(blk_cell_vtx_n);
    free(distrib_elmt);
    free(ln_to_gn);
  }

  free(pcell_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx    );
  free(pvtx_ln_to_gn);
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_ln_to_gn);
  free(pvtx_coords);

  free(section_n);
  free(g_section_n);
  free(pcell_vtx_gnum);
  free(pcell_vtx_n);


}

static
void
_rebuild_dmesh_nodal_3d
(
  PDM_dmesh_nodal_t  *dmn,
  PDM_MPI_Comm        comm,
  PDM_g_num_t        *distrib_cell,
  PDM_g_num_t        *distrib_face,
  PDM_g_num_t        *distrib_edge,
  PDM_g_num_t        *distrib_vtx,
  double             *dvtx_coords,
  PDM_g_num_t        *dcell_face,
  int                *dcell_face_idx,
  PDM_g_num_t        *dface_edge,
  int                *dface_edge_idx,
  PDM_g_num_t        *dface_vtx,
  int                *dface_vtx_idx,
  PDM_g_num_t        *dedge_vtx,
  int                *n_bound,
  int               **dbound_idx,
  PDM_g_num_t       **dbound,
  int                *n_blk_gnum,
  PDM_g_num_t       **blk_entity_gnum,
  PDM_g_num_t       **blk_elmt_gnum
)
{
  PDM_UNUSED(dmn);
  PDM_UNUSED(comm);
  PDM_UNUSED(distrib_cell);
  PDM_UNUSED(distrib_face);
  PDM_UNUSED(distrib_edge);
  PDM_UNUSED(distrib_vtx);
  PDM_UNUSED(dcell_face);
  PDM_UNUSED(dcell_face_idx);
  PDM_UNUSED(dface_edge);
  PDM_UNUSED(dface_edge_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dedge_vtx);
  PDM_UNUSED(n_bound);
  PDM_UNUSED(dbound_idx);
  PDM_UNUSED(dbound);
  PDM_UNUSED(n_blk_gnum);
  PDM_UNUSED(blk_entity_gnum);
  PDM_UNUSED(blk_elmt_gnum);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_cell = distrib_cell[i_rank+1] - distrib_cell[i_rank];
  PDM_g_num_t* pcell_ln_to_gn = malloc(dn_cell * sizeof(PDM_g_num_t));
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {
    pcell_ln_to_gn[i_cell] = distrib_cell[i_rank] + i_cell + 1;
  }


  int pn_face = 0;
  int         *pcell_face_idx = NULL;
  int         *pcell_face     = NULL;
  PDM_g_num_t *pface_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_cell,
                                                           dcell_face_idx,
                                                           dcell_face,
                                                           dn_cell,
                                                           pcell_ln_to_gn,
                                                           &pn_face,
                                                           &pface_ln_to_gn,
                                                           &pcell_face_idx,
                                                           &pcell_face);

  int pn_vtx = 0;
  int         *pface_vtx_idx = NULL;
  int         *pface_vtx     = NULL;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           pn_face,
                                                           pface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);
  double **tmp_pvtx_coords = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coords,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coords);
  double *pvtx_coords = tmp_pvtx_coords[0];
  free(tmp_pvtx_coords);

  /*
   * Reconstruction of volumic part
   */
  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */

  PDM_g_num_t          *section_n    = malloc((dn_cell+1) * sizeof(PDM_g_num_t         )); // Suralloc
  PDM_Mesh_nodal_elt_t *section_kind = malloc( dn_cell    * sizeof(PDM_Mesh_nodal_elt_t)); // Suralloc


  int cell_kind_old = -1;
  int n_section_tot    = 0;
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell) {
    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(pcell_face_idx[i_cell+1]-pcell_face_idx[i_cell],
                                                   pcell_face + pcell_face_idx[i_cell],
                                                   pface_vtx_idx,
                                                   pface_vtx,
                                                   cell_som_tria,
                                                   cell_som_quad);
    if(cell_type == (PDM_Mesh_nodal_elt_t) cell_kind_old ){
      section_n[n_section_tot-1]++;
      continue;
    }

    section_kind[n_section_tot] = cell_type;
    section_n[n_section_tot] = 1;
    n_section_tot++;
    cell_kind_old = cell_type;
  }

  section_n    = realloc(section_n   , (n_section_tot+1) * sizeof(PDM_g_num_t         ));
  section_kind = realloc(section_kind,  n_section_tot    * sizeof(PDM_Mesh_nodal_elt_t));

  PDM_g_num_t *post_section_n         = NULL;
  int         *post_section_kind      = NULL;
  int         *local_post_section_n   = NULL;
  int         *local_post_section_idx = NULL;

  int n_section_post = _generate_sections(comm,
                                          distrib_cell,
                                          section_n,
                                          section_kind,
                                          n_section_tot,
                                          &post_section_n,
                                          &post_section_kind,
                                          &local_post_section_n,
                                          &local_post_section_idx);
  // PDM_log_trace_array_int(local_post_section_n  , n_section_post, "local_post_section_n ::");
  // PDM_log_trace_array_int(local_post_section_idx, n_section_post+1, "local_post_section_idx ::");

  /*
   * Requilibrate all block
   */
  for(int i_section = 0; i_section < n_section_post; ++i_section) {

    int beg = local_post_section_idx[i_section];
    int end = local_post_section_idx[i_section+1];
    int nl_elmt = end - beg;

    PDM_g_num_t* distrib_elmt = PDM_compute_uniform_entity_distribution(comm, post_section_n[i_section]);
    PDM_g_num_t* ln_to_gn = malloc(local_post_section_n[i_section] * sizeof(PDM_g_num_t));

    for(int i = 0; i < local_post_section_n[i_section]; ++i) {
      ln_to_gn[i] = distrib_cell[i_rank] + local_post_section_idx[i_section] + i + 1;
    }

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element((PDM_Mesh_nodal_elt_t) post_section_kind[i_section], 1);

    int *pcell_vtx = malloc(nl_elmt * n_vtx_per_elmt * sizeof(int));

    for(int i_cell = 0; i_cell < local_post_section_n[i_section]; ++i_cell) {
      PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(pcell_face_idx[beg+i_cell+1]-pcell_face_idx[beg+i_cell],
                                                     pcell_face + pcell_face_idx[i_cell],
                                                     pface_vtx_idx,
                                                     pface_vtx,
                                                     cell_som_tria,
                                                     cell_som_quad);
      assert(cell_type == (PDM_Mesh_nodal_elt_t) post_section_kind[i_section]);
      int *lconnect = &pcell_vtx[i_cell * n_vtx_per_elmt];
      switch(cell_type) {
        case PDM_MESH_NODAL_TETRA4 :
          _connec_tetra(pvtx_coords,
                        cell_som_tria,
                        lconnect);
          break;
        case PDM_MESH_NODAL_HEXA8 :
          _connec_hexa(pvtx_coords,
                        cell_som_quad,
                        lconnect);
          break;
        case PDM_MESH_NODAL_PRISM6 :
          _connec_prism(pvtx_coords,
                        cell_som_tria,
                        cell_som_quad,
                        lconnect);
          break;
        case PDM_MESH_NODAL_PYRAMID5 :
          _connec_pyramid(pvtx_coords,
                          cell_som_tria,
                          cell_som_quad,
                          lconnect);
          break;
        case PDM_MESH_NODAL_POLY_3D :
          {
            abort();
            break;
          }
        default :
          break;
      }
    }

    PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                     PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                     1.,
                                                                     &ln_to_gn,
                                                                     distrib_elmt,
                                                                     &local_post_section_n[i_section],
                                                                     1,
                                                                     comm);
    int n_cell_vtx_tot = nl_elmt * n_vtx_per_elmt;
    int         *send_cell_vtx_n = malloc(nl_elmt        * sizeof(int        ));
    PDM_g_num_t *send_cell_vtx   = malloc(n_cell_vtx_tot * sizeof(PDM_g_num_t));
    // int         *blk_cell_vtx_n  = NULL;
    PDM_g_num_t *blk_cell_vtx    = NULL;

    for(int i = 0; i < n_cell_vtx_tot; ++i) {
      int i_vtx = pcell_vtx[i];
      send_cell_vtx[i] = pvtx_ln_to_gn[i_vtx-1];
    }
    free(pcell_vtx);

    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           n_vtx_per_elmt,
                           NULL,
             (void **)     &send_cell_vtx,
                           NULL,
             (void **)     &blk_cell_vtx);

    free(send_cell_vtx_n);
    free(send_cell_vtx);


    PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) post_section_kind[i_section];
    int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                 PDM_GEOMETRY_KIND_VOLUMIC,
                                                 t_elt);

    int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
    if(t_elt == PDM_MESH_NODAL_POLY_3D) {
      abort();
      // int *blk_face_vtx_idx = PDM_array_new_idx_from_sizes_int(blk_face_vtx_n, dn_elmt);
      // PDM_DMesh_nodal_section_poly2d_set(dmn,
      //                                    PDM_GEOMETRY_KIND_SURFACIC,
      //                                    id_section,
      //                                    dn_elmt,
      //                                    blk_face_vtx_idx,
      //                                    blk_face_vtx,
      //                                    PDM_OWNERSHIP_KEEP);
    } else {
      PDM_DMesh_nodal_section_std_set(dmn,
                                      PDM_GEOMETRY_KIND_VOLUMIC,
                                      id_section,
                                      dn_elmt,
                                      blk_cell_vtx,
                                      PDM_OWNERSHIP_KEEP);
    }

    // free(blk_face_vtx_n);
    // free(blk_face_vtx);

    PDM_part_to_block_free(ptb);
    free(distrib_elmt);
    free(ln_to_gn);
  }

  free(post_section_n        );
  free(post_section_kind     );
  free(local_post_section_n  );
  free(local_post_section_idx);


  free(pcell_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx    );
  free(pvtx_ln_to_gn);
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_ln_to_gn);
  free(pvtx_coords);

  /*
   * Recuperation des bords
   */
  int          n_face_bound   = n_bound   [PDM_BOUND_TYPE_FACE];
  int         *face_bound_idx = dbound_idx[PDM_BOUND_TYPE_FACE];
  PDM_g_num_t *face_bound     = dbound    [PDM_BOUND_TYPE_FACE];


  int n_face_bnd_tot = face_bound_idx[n_face_bound];
  int          pn_vtx_bnd = 0;
  int         *pface_vtx_bnd_idx = NULL;
  int         *pface_bnd_vtx     = NULL;
  PDM_g_num_t *pvtx_bnd_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           n_face_bnd_tot,
                                                           face_bound,
                                                           &pn_vtx_bnd,
                                                           &pvtx_bnd_ln_to_gn,
                                                           &pface_vtx_bnd_idx,
                                                           &pface_bnd_vtx);

  PDM_g_num_t          *section_face_bnd_n    = malloc((n_face_bnd_tot+1) * sizeof(PDM_g_num_t         )); // Suralloc
  PDM_Mesh_nodal_elt_t *section_face_bnd_kind = malloc( n_face_bnd_tot    * sizeof(PDM_Mesh_nodal_elt_t)); // Suralloc

  int ln_vtx_old = -1;
  int n_section_face_bnd_tot = 0;
  for(int i_face = 0; i_face < n_face_bnd_tot; ++i_face) {

    int ln_vtx = pface_vtx_bnd_idx[i_face+1] - pface_vtx_bnd_idx[i_face];
    if(ln_vtx_old == ln_vtx) {
      section_face_bnd_n[n_section_face_bnd_tot-1]++;
      continue;
    }
    if(ln_vtx == 3) {
      section_face_bnd_kind[n_section_face_bnd_tot] = PDM_MESH_NODAL_TRIA3;
    } else if(ln_vtx == 4){
      section_face_bnd_kind[n_section_face_bnd_tot] = PDM_MESH_NODAL_QUAD4;
    } else {
      section_face_bnd_kind[n_section_face_bnd_tot] = PDM_MESH_NODAL_POLY_2D;
    }
    section_face_bnd_n[n_section_face_bnd_tot] = 1;
    n_section_face_bnd_tot++;
    ln_vtx_old = ln_vtx;
  }
  section_face_bnd_n    = realloc(section_face_bnd_n   ,  n_section_face_bnd_tot * sizeof(PDM_g_num_t         ));
  section_face_bnd_kind = realloc(section_face_bnd_kind,  n_section_face_bnd_tot * sizeof(PDM_Mesh_nodal_elt_t));




  PDM_g_num_t *post_section_face_bnd_n         = NULL;
  int         *post_section_face_bnd_kind      = NULL;
  int         *local_post_section_face_bnd_n   = NULL;
  int         *local_post_section_face_bnd_idx = NULL;

  PDM_g_num_t* distrib_face_bnd = PDM_compute_entity_distribution(comm, n_face_bnd_tot);
  int n_section_face_bnd_post = _generate_sections(comm,
                                                   distrib_face_bnd,
                                                   section_face_bnd_n,
                                                   section_face_bnd_kind,
                                                   n_section_face_bnd_tot,
                                                   &post_section_face_bnd_n,
                                                   &post_section_face_bnd_kind,
                                                   &local_post_section_face_bnd_n,
                                                   &local_post_section_face_bnd_idx);

  // PDM_log_trace_array_long(section_face_bnd_n, n_section_face_bnd_tot, "section_face_bnd_n ::");
  // PDM_log_trace_array_int(local_post_section_face_bnd_n  , n_section_face_bnd_post, "local_post_section_face_bnd_n ::");
  // PDM_log_trace_array_int(local_post_section_face_bnd_idx, n_section_face_bnd_post+1, "local_post_section_face_bnd_idx ::");
  /*
   * Generate and redistribute surfacique
   */
  for(int i_section = 0; i_section < n_section_face_bnd_post; ++i_section) {

    int beg = local_post_section_face_bnd_idx[i_section];
    int end = local_post_section_face_bnd_idx[i_section+1];
    int nl_elmt = end - beg;

    // PDM_g_num_t* distrib_elmt = PDM_compute_uniform_entity_distribution(comm, post_section_face_bnd_n[i_section]);
    PDM_g_num_t* ln_to_gn = malloc(nl_elmt * sizeof(PDM_g_num_t));

    for(int i = 0; i < nl_elmt; ++i) {
      ln_to_gn[i] = local_post_section_face_bnd_idx[i_section] + distrib_face_bnd[i_rank] + i + 1;
    }

    // PDM_log_trace_array_long(ln_to_gn, nl_elmt, "face_bnd_ln_to_gn ::");

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element((PDM_Mesh_nodal_elt_t) post_section_face_bnd_kind[i_section], 1);
    int n_cell_vtx_tot = nl_elmt * n_vtx_per_elmt;
    PDM_g_num_t* face_vtx_gnum = (PDM_g_num_t *) malloc(n_cell_vtx_tot * sizeof(PDM_g_num_t));
    for(int i_face = 0; i_face < nl_elmt; ++i_face) {
      int idx_read = pface_vtx_bnd_idx[beg+i_face];
      for(int j = 0; j < n_vtx_per_elmt; ++j) {
        face_vtx_gnum[n_vtx_per_elmt*i_face+j] = pvtx_bnd_ln_to_gn[pface_bnd_vtx[idx_read+j]-1];
      }
    }

    // PDM_log_trace_array_long(face_vtx_gnum, nl_elmt * n_vtx_per_elmt, "face_vtx_gnum :");
    // PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
    //                                                                  PDM_PART_TO_BLOCK_POST_CLEANUP,
    //                                                                  1.,
    //                                                                  &ln_to_gn,
    //                                                                  distrib_elmt,
    //                                                                  &nl_elmt,
    //                                                                  1,
    //                                                                  comm);
    PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        &ln_to_gn,
                                                        NULL,
                                                        &nl_elmt,
                                                        1,
                                                        comm);

    // PDM_g_num_t* distrib_elmt = PDM_part_to_block_distrib_index_get(ptb);

    printf("n_vtx_per_elmt = %i \n", n_vtx_per_elmt);
    printf("nl_elmt = %i \n", nl_elmt);

    PDM_g_num_t *blk_face_vtx    = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           n_vtx_per_elmt,
                           NULL,
             (void **)     &face_vtx_gnum,
                           NULL,
             (void **)     &blk_face_vtx);

    PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) post_section_face_bnd_kind[i_section];
    int id_section = PDM_DMesh_nodal_section_add(dmn,
                                                 PDM_GEOMETRY_KIND_SURFACIC,
                                                 t_elt);



    int dn_elmt = PDM_part_to_block_n_elt_block_get(ptb); //distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
    if(t_elt == PDM_MESH_NODAL_POLY_2D) {
      abort();
    } else {
      // PDM_log_trace_array_long(blk_face_vtx, dn_elmt * n_vtx_per_elmt, "blk_face_vtx :");
      PDM_DMesh_nodal_section_std_set(dmn,
                                      PDM_GEOMETRY_KIND_SURFACIC,
                                      id_section,
                                      dn_elmt,
                                      blk_face_vtx,
                                      PDM_OWNERSHIP_KEEP);
    }

    PDM_part_to_block_free(ptb);
    free(face_vtx_gnum);
    // free(distrib_elmt);
    free(ln_to_gn);
  }





  free(post_section_face_bnd_n        );
  free(post_section_face_bnd_kind     );
  free(local_post_section_face_bnd_n  );
  free(local_post_section_face_bnd_idx);
  free(distrib_face_bnd);

  free(pface_vtx_bnd_idx);
  free(pface_bnd_vtx    );
  free(pvtx_bnd_ln_to_gn);

}


static
PDM_dmesh_nodal_t*
_dmesh_to_dmesh_nodal
(
 PDM_MPI_Comm    comm,
 PDM_g_num_t    *distrib_cell,
 PDM_g_num_t    *distrib_face,
 PDM_g_num_t    *distrib_edge,
 PDM_g_num_t    *distrib_vtx,
 double         *dvtx_coords,
 PDM_g_num_t    *dcell_face,
 int            *dcell_face_idx,
 PDM_g_num_t    *dface_edge,
 int            *dface_edge_idx,
 PDM_g_num_t    *dface_vtx,
 int            *dface_vtx_idx,
 PDM_g_num_t    *dedge_vtx,
 int            *n_bound,
 int           **dbound_idx,
 PDM_g_num_t   **dbound,
 int            *n_blk_gnum,
 PDM_g_num_t   **blk_entity_gnum,
 PDM_g_num_t   **blk_elmt_gnum
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int is_3d = 0;
  int is_2d = 0;
  int mesh_dimension = 3;
  if(distrib_cell != NULL) {
    is_3d = 1;
    mesh_dimension = 3;
  } else {
    assert(distrib_face != NULL);
    is_2d = 1;
    mesh_dimension = 2;
  }
  printf("_dmesh_to_dmesh_nodal --> is_2d = %i | is_3d = %i \n", is_2d, is_3d);

  PDM_g_num_t  n_vtx  = 0;
  PDM_g_num_t  n_cell = 0;
  PDM_g_num_t  n_face = 0;
  PDM_g_num_t  n_edge = 0;

  if(distrib_cell != NULL) {
    n_cell = distrib_cell[n_rank];
  }

  if(distrib_face != NULL) {
    n_face = distrib_face[n_rank];
  }

  if(distrib_edge != NULL) {
    n_edge = distrib_edge[n_rank];
  }

  if(distrib_vtx != NULL) {
    n_vtx = distrib_vtx[n_rank];
  }

  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(comm,
                                                  mesh_dimension,
                                                  n_vtx,
                                                  n_cell,
                                                  n_face,
                                                  n_edge);

  if(is_3d) {
    if(1 == 1) {
      _rebuild_dmesh_nodal_3d(dmn,
                              comm,
                              distrib_cell,
                              distrib_face,
                              distrib_edge,
                              distrib_vtx,
                              dvtx_coords,
                              dcell_face,
                              dcell_face_idx,
                              dface_edge,
                              dface_edge_idx,
                              dface_vtx,
                              dface_vtx_idx,
                              dedge_vtx,
                              n_bound,
                              dbound_idx,
                              dbound,
                              n_blk_gnum,
                              blk_entity_gnum,
                              blk_elmt_gnum);
    } else {
      _rebuild_dmesh_nodal_by_kind_3d(dmn,
                                      comm,
                                      distrib_cell,
                                      distrib_face,
                                      distrib_edge,
                                      distrib_vtx,
                                      dvtx_coords,
                                      dcell_face,
                                      dcell_face_idx,
                                      dface_edge,
                                      dface_edge_idx,
                                      dface_vtx,
                                      dface_vtx_idx,
                                      dedge_vtx,
                                      n_bound,
                                      dbound_idx,
                                      dbound,
                                      n_blk_gnum,
                                      blk_entity_gnum,
                                      blk_elmt_gnum);
    }
  } else {
    _rebuild_dmesh_nodal_2d(dmn,
                            comm,
                            distrib_face,
                            distrib_edge,
                            distrib_vtx,
                            dface_edge,
                            dface_edge_idx,
                            dface_vtx,
                            dface_vtx_idx,
                            dedge_vtx,
                            n_bound,
                            dbound_idx,
                            dbound,
                            n_blk_gnum,
                            blk_entity_gnum,
                            blk_elmt_gnum);
  }


  return dmn;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_to_dmesh_nodal_t*
PDM_dmesh_to_dmesh_nodal_create
(
 const int             n_mesh,
 const PDM_MPI_Comm    comm
)
{

  PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn = (PDM_dmesh_to_dmesh_nodal_t *) malloc(sizeof(PDM_dmesh_to_dmesh_nodal_t));

  dm_to_dmn->comm              = comm;
  dm_to_dmn->results_is_getted = PDM_FALSE;
  dm_to_dmn->n_mesh            = n_mesh;

  dm_to_dmn->dcell_face            = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dcell_face_idx        = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_edge            = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_edge_idx        = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dface_vtx             = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dface_vtx_idx         = malloc(n_mesh * sizeof(int         *));
  dm_to_dmn->dedge_vtx             = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->dvtx_coords           = malloc(n_mesh * sizeof(double      *));
  dm_to_dmn->dparent_elmt_position = malloc(n_mesh * sizeof(int      *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dcell_face           [i_mesh] = NULL;
    dm_to_dmn->dcell_face_idx       [i_mesh] = NULL;
    dm_to_dmn->dface_edge           [i_mesh] = NULL;
    dm_to_dmn->dface_edge_idx       [i_mesh] = NULL;
    dm_to_dmn->dface_vtx            [i_mesh] = NULL;
    dm_to_dmn->dface_vtx_idx        [i_mesh] = NULL;
    dm_to_dmn->dedge_vtx            [i_mesh] = NULL;
    dm_to_dmn->dvtx_coords          [i_mesh] = NULL;
    dm_to_dmn->dparent_elmt_position[i_mesh] = NULL;
  }

  dm_to_dmn->distrib_cell = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_face = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_edge = malloc(n_mesh * sizeof(PDM_g_num_t *));
  dm_to_dmn->distrib_vtx  = malloc(n_mesh * sizeof(PDM_g_num_t *));

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->distrib_cell[i_mesh] = NULL;
    dm_to_dmn->distrib_face[i_mesh] = NULL;
    dm_to_dmn->distrib_edge[i_mesh] = NULL;
    dm_to_dmn->distrib_vtx [i_mesh] = NULL;
  }

  dm_to_dmn->n_bound         = malloc( n_mesh * sizeof(int          *) );
  dm_to_dmn->dbound          = malloc( n_mesh * sizeof(PDM_g_num_t **) );
  dm_to_dmn->dbound_idx      = malloc( n_mesh * sizeof(int         **) );


  dm_to_dmn->n_blk_gnum      = malloc( n_mesh * sizeof(int          *) );
  dm_to_dmn->blk_entity_gnum = malloc( n_mesh * sizeof(PDM_g_num_t **) );
  dm_to_dmn->blk_elmt_gnum   = malloc( n_mesh * sizeof(PDM_g_num_t **) );

  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {

    dm_to_dmn->n_bound   [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          ) );
    dm_to_dmn->dbound    [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
    dm_to_dmn->dbound_idx[i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );

    dm_to_dmn->n_blk_gnum     [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          ) );
    dm_to_dmn->blk_entity_gnum[i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
    dm_to_dmn->blk_elmt_gnum  [i_mesh] = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
      dm_to_dmn->n_bound   [i_mesh][i] = 0;
      dm_to_dmn->dbound    [i_mesh][i] = NULL;
      dm_to_dmn->dbound_idx[i_mesh][i] = NULL;

      dm_to_dmn->n_blk_gnum     [i_mesh][i] = 0;
      dm_to_dmn->blk_entity_gnum[i_mesh][i] = NULL;
      dm_to_dmn->blk_elmt_gnum  [i_mesh][i] = NULL;
    }
  }

  dm_to_dmn->dmn = malloc( n_mesh * sizeof(PDM_dmesh_nodal_t *) );
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = NULL;
  }

  return dm_to_dmn;
}

void
PDM_dmesh_to_dmesh_nodal_compute
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
    dm_to_dmn->dmn[i_mesh] = _dmesh_to_dmesh_nodal(dm_to_dmn->comm,
                                                   dm_to_dmn->distrib_cell   [i_mesh],
                                                   dm_to_dmn->distrib_face   [i_mesh],
                                                   dm_to_dmn->distrib_edge   [i_mesh],
                                                   dm_to_dmn->distrib_vtx    [i_mesh],
                                                   dm_to_dmn->dvtx_coords    [i_mesh],
                                                   dm_to_dmn->dcell_face     [i_mesh],
                                                   dm_to_dmn->dcell_face_idx [i_mesh],
                                                   dm_to_dmn->dface_edge     [i_mesh],
                                                   dm_to_dmn->dface_edge_idx [i_mesh],
                                                   dm_to_dmn->dface_vtx      [i_mesh],
                                                   dm_to_dmn->dface_vtx_idx  [i_mesh],
                                                   dm_to_dmn->dedge_vtx      [i_mesh],
                                                   dm_to_dmn->n_bound        [i_mesh],
                                                   dm_to_dmn->dbound_idx     [i_mesh],
                                                   dm_to_dmn->dbound         [i_mesh],
                                                   dm_to_dmn->n_blk_gnum     [i_mesh],
                                                   dm_to_dmn->blk_entity_gnum[i_mesh],
                                                   dm_to_dmn->blk_elmt_gnum  [i_mesh]);
  }


}

void
PDM_dmesh_to_dmesh_nodal_set_dmesh
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_t                *dm
)
{
  PDM_UNUSED(dm_to_dmn);
  PDM_UNUSED(i_mesh);
  PDM_UNUSED(dm);
  abort();
}


void
PDM_dmesh_to_dmesh_nodal_dmesh_nodal_get
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_dmesh_nodal_t         **dmn,
        PDM_ownership_t             ownership
)
{
  if(ownership == PDM_OWNERSHIP_USER) {
    dm_to_dmn->results_is_getted = PDM_TRUE;
  }
  *dmn = dm_to_dmn->dmn[i_mesh];
}

void
PDM_dmesh_to_dmesh_nodal_distribution_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_g_num_t                *distrib_cell,
        PDM_g_num_t                *distrib_face,
        PDM_g_num_t                *distrib_edge,
        PDM_g_num_t                *distrib_vtx
)
{
  dm_to_dmn->distrib_cell[i_mesh] = distrib_cell;
  dm_to_dmn->distrib_face[i_mesh] = distrib_face;
  dm_to_dmn->distrib_edge[i_mesh] = distrib_edge;
  dm_to_dmn->distrib_vtx [i_mesh] = distrib_vtx;
}

void
PDM_dmesh_to_dmesh_nodal_connectivity_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        int                        *dcell_face_idx,
        PDM_g_num_t                *dcell_face,
        int                        *dface_edge_idx,
        PDM_g_num_t                *dface_edge,
        PDM_g_num_t                *dedge_vtx,
        int                        *dface_vtx_idx,
        PDM_g_num_t                *dface_vtx,
        double                     *dvtx_coords
)
{
  dm_to_dmn->dcell_face_idx[i_mesh] = dcell_face_idx;
  dm_to_dmn->dcell_face    [i_mesh] = dcell_face;
  dm_to_dmn->dface_edge    [i_mesh] = dface_edge;
  dm_to_dmn->dface_edge_idx[i_mesh] = dface_edge_idx;
  dm_to_dmn->dedge_vtx     [i_mesh] = dedge_vtx;
  dm_to_dmn->dface_vtx_idx [i_mesh] = dface_vtx_idx;
  dm_to_dmn->dface_vtx     [i_mesh] = dface_vtx;
  dm_to_dmn->dvtx_coords   [i_mesh] = dvtx_coords;
}

void
PDM_dmesh_to_dmesh_nodal_group_set
(
        PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
  const int                         i_mesh,
        PDM_bound_type_t            bound_type,
        int                         n_group,
        int                        *dbound_idx,
        PDM_g_num_t                *dbound
)
{
  dm_to_dmn->n_bound   [i_mesh][bound_type] = n_group;
  dm_to_dmn->dbound_idx[i_mesh][bound_type] = dbound_idx;
  dm_to_dmn->dbound    [i_mesh][bound_type] = dbound;
}


void
PDM_dmesh_to_dmesh_nodal_free
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn
)
{
  if(dm_to_dmn == NULL) {
    return;
  }

  free(dm_to_dmn->dcell_face    );
  free(dm_to_dmn->dcell_face_idx);
  free(dm_to_dmn->dface_edge    );
  free(dm_to_dmn->dface_edge_idx);
  free(dm_to_dmn->dface_vtx     );
  free(dm_to_dmn->dface_vtx_idx );
  free(dm_to_dmn->dedge_vtx     );
  free(dm_to_dmn->dvtx_coords   );
  free(dm_to_dmn->dparent_elmt_position);

  free(dm_to_dmn->distrib_cell);
  free(dm_to_dmn->distrib_face);
  free(dm_to_dmn->distrib_edge);
  free(dm_to_dmn->distrib_vtx );

  for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh ) {
    free(dm_to_dmn->n_bound   [i_mesh]);
    free(dm_to_dmn->dbound    [i_mesh]);
    free(dm_to_dmn->dbound_idx[i_mesh]);

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i){
      if(dm_to_dmn->blk_entity_gnum[i_mesh][i] != NULL){
        free(dm_to_dmn->blk_entity_gnum[i_mesh][i]);
      }
      if(dm_to_dmn->blk_elmt_gnum[i_mesh][i] != NULL){
        free(dm_to_dmn->blk_elmt_gnum[i_mesh][i]);
      }
    }

    free(dm_to_dmn->n_blk_gnum     [i_mesh]);
    free(dm_to_dmn->blk_entity_gnum[i_mesh]);
    free(dm_to_dmn->blk_elmt_gnum  [i_mesh]);

  }

  free(dm_to_dmn->n_blk_gnum     );
  free(dm_to_dmn->blk_entity_gnum);
  free(dm_to_dmn->blk_elmt_gnum  );

  free(dm_to_dmn->n_bound      );
  free(dm_to_dmn->dbound       );
  free(dm_to_dmn->dbound_idx   );

  if(dm_to_dmn->results_is_getted == PDM_FALSE) {
    for(int i_mesh = 0; i_mesh < dm_to_dmn->n_mesh; ++i_mesh) {
      PDM_DMesh_nodal_free(dm_to_dmn->dmn[i_mesh]);
    }
  }

  free(dm_to_dmn->dmn);

  free(dm_to_dmn);
  dm_to_dmn = NULL;
}



void
PDM_dmesh_to_dmesh_nodal_update_group
(
 PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn,
 int                         i_mesh,
 PDM_bound_type_t            bound_type,
 int                         dn_elmt_bound,
 PDM_g_num_t                *dentity_bound,
 PDM_g_num_t               **delmt_bound
)
{

  int n_blk_gnum           = dm_to_dmn->n_blk_gnum     [i_mesh][bound_type];
  PDM_g_num_t* gnum_entity = dm_to_dmn->blk_entity_gnum[i_mesh][bound_type];
  PDM_g_num_t* gnum_elmt   = dm_to_dmn->blk_elmt_gnum  [i_mesh][bound_type];

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(gnum_entity,
                                                                        n_blk_gnum,
                                            (const PDM_g_num_t **)      &dentity_bound,
                                                                        &dn_elmt_bound,
                                                                        1,
                                                                        dm_to_dmn->comm);

  int stride_one = 1;
  PDM_g_num_t **tmp_entity_to_elmt = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         gnum_elmt,
                         NULL,
              (void ***) &tmp_entity_to_elmt);
  PDM_g_num_t *entity_to_elmt = tmp_entity_to_elmt[0];
  free(tmp_entity_to_elmt);

  // PDM_log_trace_array_long(entity_to_elmt, dn_elmt_bound, "entity_to_elmt ::");

  *delmt_bound = entity_to_elmt;

  PDM_block_to_part_free(btp);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
