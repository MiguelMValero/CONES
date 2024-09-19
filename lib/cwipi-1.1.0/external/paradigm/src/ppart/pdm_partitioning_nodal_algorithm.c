/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_hash_tab.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"

#include "pdm_partitioning_nodal_algorithm.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_order.h"
#include "pdm_logging.h"
// #include "pdm_para_graph_dual.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

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
void
_delmt_vtx_to_pelmt_vtx
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           n_part,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn,
 int                        ***pelmt_strid_idx_out,
 PDM_g_num_t                ***pelmt_connec_out,
 PDM_Mesh_nodal_elt_t       ***pelmt_types_out
)
{
  int n_section = dmne->n_section;
  PDM_g_num_t          **block_elmts_disbrib_idx = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  PDM_g_num_t          **block_elmts_connec      = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  int                  **block_elmts_n_vtx       = (int                  ** ) malloc( n_section * sizeof(int                  *));
  PDM_Mesh_nodal_elt_t **block_elmts_types       = (PDM_Mesh_nodal_elt_t ** ) malloc( n_section * sizeof(PDM_Mesh_nodal_elt_t *));
  int                  **stride_one              = (int                  ** ) malloc( n_section * sizeof(int                  *));

  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];
    block_elmts_disbrib_idx[i_section] = (PDM_g_num_t *) PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

    stride_one[i_section] = (int * ) malloc( 1 * sizeof(int));
    stride_one[i_section][0] = 1;


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
        int n_elt           = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, 1);
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = n_vtx_per_elmt;
          block_elmts_types[i_section][i] = t_elt;
        }
        block_elmts_connec[i_section] = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);

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
        int order = -1;
        const char *ho_ordering = NULL;
        block_elmts_connec[i_section] = PDM_DMesh_nodal_elmts_section_std_ho_get(dmne,
                                                                                 id_section,
                                                                                 &order,
                                                                                 &ho_ordering);

        int n_elt           = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = n_vtx_per_elmt;
          block_elmts_types[i_section][i] = t_elt;
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      case PDM_MESH_NODAL_POLY_3D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
    }

  }

  // PDM_g_num_t* volumic_distrib =
  PDM_multi_block_to_part_t* mbtp = PDM_multi_block_to_part_create(dmne->section_distribution,
                                                                   n_section,
                                            (const PDM_g_num_t **) block_elmts_disbrib_idx,
                                            (const PDM_g_num_t **) elmt_ln_to_gn,
                                            (const int          *) pn_elmt,
                                                                   n_part,
                                                                   dmne->comm);

  free(block_elmts_disbrib_idx);

  /*
   * Exchange connectivity
   */
  int         **pelmts_stride = NULL;
  PDM_g_num_t **pelmts_connec = NULL;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
                                block_elmts_n_vtx,
                     (void ** ) block_elmts_connec,
                     (int  ***) &pelmts_stride,
                     (void ***) &pelmts_connec);

  int         **pelmts_stride_idx = malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pelmts_stride_idx[i_part] = PDM_array_new_idx_from_sizes_int(pelmts_stride[i_part], pn_elmt[i_part]);
    free(pelmts_stride[i_part]);
  }
  free(pelmts_stride);
  /*
   * Exchange type of elements
   */
  PDM_Mesh_nodal_elt_t **pelmts_types;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_Mesh_nodal_elt_t),
                                PDM_STRIDE_CST_INTERLACED,
                                stride_one,
                     (void ** ) block_elmts_types,
                                NULL,
                     (void ***) &pelmts_types);



  for (int i_section = 0; i_section < n_section; i_section++) {
    free(block_elmts_n_vtx[i_section]);
    free(block_elmts_types[i_section]);
    free(stride_one[i_section]);
  }
  free(block_elmts_n_vtx);
  free(block_elmts_types);
  free(block_elmts_connec);
  free(stride_one);

  /* Set output */
  *pelmt_strid_idx_out = pelmts_stride_idx;
  *pelmt_connec_out    = pelmts_connec;
  *pelmt_types_out     = pelmts_types;

  PDM_multi_block_to_part_free(mbtp);
}

// static
// void
// (
// )
// {

// }


/*=============================================================================
 * Public function definitions
 *============================================================================*/


PDM_part_mesh_nodal_elmts_t*
PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           n_part,
 int                          *pn_vtx,
 PDM_g_num_t                 **vtx_ln_to_gn,
 int                          *pn_elmt,
 int                         **pelmt_to_entity,
 PDM_g_num_t                 **elmt_ln_to_gn,
 PDM_g_num_t                 **pparent_entitity_ln_to_gn
)
{

  int n_section = dmne->n_section;
  PDM_g_num_t          **block_elmts_disbrib_idx = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  PDM_g_num_t          **block_elmts_connec      = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  int                  **block_elmts_n_vtx       = (int                  ** ) malloc( n_section * sizeof(int                  *));
  PDM_Mesh_nodal_elt_t **block_elmts_types       = (PDM_Mesh_nodal_elt_t ** ) malloc( n_section * sizeof(PDM_Mesh_nodal_elt_t *));
  int                  **stride_one              = (int                  ** ) malloc( n_section * sizeof(int                  *));
  int                   *pid_section             = (int                  *  ) malloc( n_section * sizeof(int                   ));
  int                   *section_order           = (int                  *  ) malloc( n_section * sizeof(int                   ));
  const char           **section_ho_ordering     = (const char           ** ) malloc( n_section * sizeof(char                 *));

  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];
    block_elmts_disbrib_idx[i_section] = (PDM_g_num_t *) PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

    stride_one[i_section] = (int * ) malloc( 1 * sizeof(int));
    stride_one[i_section][0] = 1;


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
        section_order      [i_section] = 1;
        section_ho_ordering[i_section] = NULL;
        // pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        int n_elt = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, 1);
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = n_vtx_per_elmt;
          block_elmts_types[i_section][i] = t_elt;
        }
        block_elmts_connec[i_section] = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);

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
        int order = -1;
        const char *ho_ordering = NULL;
        block_elmts_connec[i_section] = PDM_DMesh_nodal_elmts_section_std_ho_get(dmne,
                                                                                 id_section,
                                                                                 &order,
                                                                                 &ho_ordering);

        section_order      [i_section] = order;
        section_ho_ordering[i_section] = ho_ordering;
        // pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        int n_elt = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = n_vtx_per_elmt;
          block_elmts_types[i_section][i] = t_elt;
        }
        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        section_order      [i_section] = 1;
        section_ho_ordering[i_section] = NULL;
        // pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        int n_elt = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);

        int         *connec_idx = NULL;
        PDM_g_num_t *connec     = NULL;
        PDM_DMesh_nodal_elmts_section_poly2d_get(dmne,
                                                 id_section,
                                                 &connec_idx,
                                                 &connec);

        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = connec_idx[i+1] - connec_idx[i];
          block_elmts_types[i_section][i] = t_elt;
        }
        block_elmts_connec[i_section] = connec;
        break;
      }

      case PDM_MESH_NODAL_POLY_3D:
      {
        section_order      [i_section] = 1;
        section_ho_ordering[i_section] = NULL;
        // pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts : Element type %d is not supported\n", (int) t_elt);
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts : Element type %d is not supported\n", (int) t_elt);
    }

  }

  // PDM_g_num_t* volumic_distrib =
  PDM_multi_block_to_part_t* mbtp = PDM_multi_block_to_part_create(dmne->section_distribution,
                                                                   n_section,
                                            (const PDM_g_num_t **) block_elmts_disbrib_idx,
                                            (const PDM_g_num_t **) elmt_ln_to_gn,
                                            (const int          *) pn_elmt,
                                                                   n_part,
                                                                   dmne->comm);

  free(block_elmts_disbrib_idx);

  /*
   * Exchange connectivity
   */
  int         **pelmts_stride = NULL;
  PDM_g_num_t **pelmts_connec = NULL;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
                                block_elmts_n_vtx,
                     (void ** ) block_elmts_connec,
                     (int  ***) &pelmts_stride,
                     (void ***) &pelmts_connec);

  int         **pelmts_stride_idx = malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pelmts_stride_idx[i_part] = PDM_array_new_idx_from_sizes_int(pelmts_stride[i_part], pn_elmt[i_part]);
  }

  /*
   * Exchange type of elements
   */
  PDM_Mesh_nodal_elt_t **pelmts_types;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_Mesh_nodal_elt_t),
                                PDM_STRIDE_CST_INTERLACED,
                                stride_one,
                     (void ** ) block_elmts_types,
                                NULL,
                     (void ***) &pelmts_types);



  for (int i_section = 0; i_section < n_section; i_section++) {
    free(block_elmts_n_vtx[i_section]);
    free(block_elmts_types[i_section]);
    free(stride_one[i_section]);
  }
  free(block_elmts_n_vtx);
  free(block_elmts_types);
  free(block_elmts_connec);
  free(stride_one);

  PDM_multi_block_to_part_free(mbtp);

  PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_elmts_create(dmne->mesh_dimension, n_part, dmne->comm);

  /* Create skeleton */
  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];
    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

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
        pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
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
        pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        break;
      }
      case PDM_MESH_NODAL_POLY_3D:
      {
        pid_section        [i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts : Element type %d is not supported\n", (int) t_elt);
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts : Element type %d is not supported\n", (int) t_elt);
    }

  }


  /*
   * A priori le vtx_ln_to_gn n'est pas trié
   */
  PDM_g_num_t **sorted_vtx_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int         **vtx_order           = (int         ** ) malloc( n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    sorted_vtx_ln_to_gn[i_part] = malloc(pn_vtx[i_part] * sizeof(PDM_g_num_t));
    vtx_order          [i_part] = malloc(pn_vtx[i_part] * sizeof(int        ));

    for(int i = 0; i < pn_vtx[i_part]; ++i) {
      sorted_vtx_ln_to_gn[i_part][i] = vtx_ln_to_gn[i_part][i];
      vtx_order          [i_part][i] = i;
    }
    PDM_sort_long(sorted_vtx_ln_to_gn[i_part], vtx_order[i_part], pn_vtx[i_part]);
    // PDM_log_trace_array_long(sorted_vtx_ln_to_gn[i_part] , pn_vtx[i_part] , "sorted_vtx_ln_to_gn :: ");
  }

  /*
   *  We don't need to exchange the section because we find it with binary search on section_distribution
   */
  int          *pelmt_by_section_n        = (int          *) malloc( (n_section) * sizeof(int          ));
  int         **connec                    = (int         **) malloc( (n_section) * sizeof(int         *));
  int         **connec_idx                = (int         **) malloc( (n_section) * sizeof(int         *));
  PDM_g_num_t **numabs                    = (PDM_g_num_t **) malloc( (n_section) * sizeof(PDM_g_num_t *));
  int         **parent_num                = (int         **) malloc( (n_section) * sizeof(int         *));
  PDM_g_num_t **sparent_entitity_ln_to_gn = (PDM_g_num_t **) malloc( (n_section) * sizeof(PDM_g_num_t *));
  PDM_g_num_t **section_elmts_ln_to_gn    = (PDM_g_num_t **) malloc( (n_part   ) * sizeof(PDM_g_num_t *));

  for(int i_section = 0; i_section < n_section; ++i_section){
    sparent_entitity_ln_to_gn[i_section] = NULL;
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    /* Reset */
    PDM_array_reset_int(pelmt_by_section_n, n_section, 0);

    for(int i_elmt = 0; i_elmt < pn_elmt[i_part]; ++i_elmt) {

      PDM_g_num_t g_num = elmt_ln_to_gn[i_part][i_elmt]-1;
      int i_section = PDM_binary_search_gap_long(g_num, dmne->section_distribution, n_section+1);

      /* We need to sort entry in each section */
      pelmt_by_section_n[i_section]++;
    }

    /* We allocate here and ownership if tranfert to PDM_part_mesh_nodal_elmts_t*/
    for(int i_section = 0; i_section < n_section; ++i_section){
      int n_elmt_in_section = pelmt_by_section_n[i_section];
      // int order      = section_order[i_section];
      int id_section = pid_section  [i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      // int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      // connec    [i_section] = malloc( n_elmt_in_section * n_vtx_per_elmt * sizeof(int        ));
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        connec_idx[i_section] = malloc(sizeof(int) * (n_elmt_in_section + 1));
        connec_idx[i_section][0] = 0;
        pelmt_by_section_n[i_section] = 0;
      }
      numabs    [i_section] = malloc( n_elmt_in_section                  * sizeof(PDM_g_num_t));
      parent_num[i_section] = malloc( n_elmt_in_section                  * sizeof(int        ));

      if(pparent_entitity_ln_to_gn != NULL) {
        sparent_entitity_ln_to_gn[i_section] = malloc( n_elmt_in_section * sizeof(PDM_g_num_t));
      }

    }

    for(int i_elmt = 0; i_elmt < pn_elmt[i_part]; ++i_elmt) {
      PDM_g_num_t g_num = elmt_ln_to_gn[i_part][i_elmt]-1;
      int i_section = PDM_binary_search_gap_long(g_num, dmne->section_distribution, n_section+1);

      int id_section = pid_section  [i_section];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int idx_write = pelmt_by_section_n[i_section]++;
        connec_idx[i_section][idx_write+1] = pelmts_stride_idx[i_part][i_elmt+1] - pelmts_stride_idx[i_part][i_elmt];
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        PDM_error (__FILE__, __LINE__, 0, "Poly3d are not yet supported\n");
      }
    }

    for(int i_section = 0; i_section < n_section; ++i_section){
      int id_section = pid_section[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        for (int i = 0; i < pelmt_by_section_n[i_section]; i++) {
          connec_idx[i_section][i+1] += connec_idx[i_section][i];
        }
        connec[i_section] = malloc(sizeof(int) * connec_idx[i_section][pelmt_by_section_n[i_section]]);
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        PDM_error (__FILE__, __LINE__, 0, "Poly3d are not yet supported\n");
      }
      else {
        int order = section_order[i_section];
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
        connec[i_section] = malloc(pelmt_by_section_n[i_section] * n_vtx_per_elmt * sizeof(int));
      }
    }
    PDM_array_reset_int(pelmt_by_section_n, n_section, 0);

    /* For each section we rebuild the connectivity and the parent_num */
    for(int i_elmt = 0; i_elmt < pn_elmt[i_part]; ++i_elmt) {

      PDM_g_num_t g_num = elmt_ln_to_gn[i_part][i_elmt]-1;
      int i_section = PDM_binary_search_gap_long(g_num, dmne->section_distribution, n_section+1);

      int id_section = pid_section  [i_section];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);
      int idx_write = pelmt_by_section_n[i_section]++;

      int idx_read_connec = pelmts_stride_idx[i_part][i_elmt];

      int idx_write_connec = 0;
      int n_vtx_per_elmt = 0;
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        idx_write_connec = connec_idx[i_section][idx_write];
        n_vtx_per_elmt   = connec_idx[i_section][idx_write+1] - idx_write_connec;
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        abort();
      }
      else {
        int order      = section_order[i_section];
        n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
        idx_write_connec = n_vtx_per_elmt * idx_write;
      }


      for(int i_vtx = 0; i_vtx < n_vtx_per_elmt; ++i_vtx){
        PDM_g_num_t vtx_g_num = pelmts_connec[i_part][idx_read_connec+i_vtx];
        int vtx_l_num = PDM_binary_search_long(vtx_g_num, sorted_vtx_ln_to_gn[i_part], pn_vtx[i_part]);
        assert(vtx_l_num != -1);
        // connec[i_section][n_vtx_per_elmt*idx_write + i_vtx] = vtx_order[i_part][vtx_l_num]+1;
        connec[i_section][idx_write_connec + i_vtx] = vtx_order[i_part][vtx_l_num]+1;
      }

      numabs    [i_section][idx_write] = g_num+1;
      if(pelmt_to_entity != NULL) {
        parent_num[i_section][idx_write] = pelmt_to_entity[i_part][i_elmt];
      }
      else {
        parent_num[i_section][idx_write] = i_elmt;//+1; // 0-based to be coherent with PDM_Mesh_nodal

      }
      if(pparent_entitity_ln_to_gn != NULL) {
        sparent_entitity_ln_to_gn[i_section][idx_write] = pparent_entitity_ln_to_gn[i_part][i_elmt];
      }
    }

    /* Keep a true element gnum */
    section_elmts_ln_to_gn[i_part] = malloc(pn_elmt[i_part] * sizeof(PDM_g_num_t));

    int idx_elmt = 0;
    for(int i_section = 0; i_section < n_section; ++i_section){
      int n_elmt_in_section = pelmt_by_section_n[i_section];
      for(int i_elmt = 0; i_elmt < n_elmt_in_section; ++i_elmt) {
        //int i_parent = parent_num[i_section][i_elmt];
        section_elmts_ln_to_gn[i_part][idx_elmt++] = elmt_ln_to_gn[i_part][i_elmt];
      }
    }
    assert(idx_elmt == pn_elmt[i_part]);


    /* Fill up structure */
    for(int i_section = 0; i_section < n_section; ++i_section){
      int n_elmt_in_section = pelmt_by_section_n[i_section];
      int id_section = pid_section[i_section];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_part_mesh_nodal_elmts_section_poly2d_set(pmne,
                                                   id_section,
                                                   i_part,
                                                   n_elmt_in_section,
                                                   connec_idx[i_section],
                                                   connec    [i_section],
                                                   numabs    [i_section],
                                                   parent_num[i_section],
                                                   PDM_OWNERSHIP_KEEP);
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        abort();
      }
      else {
        PDM_part_mesh_nodal_elmts_std_ho_set(pmne,
                                             id_section,
                                             i_part,
                                             n_elmt_in_section,
                                             connec                   [i_section],
                                             numabs                   [i_section],
                                             parent_num               [i_section],
                                             sparent_entitity_ln_to_gn[i_section],
                                             section_order            [i_section],
                                             section_ho_ordering      [i_section],
                                             PDM_OWNERSHIP_KEEP);
      }

      connec                   [i_section] = NULL;
      connec_idx               [i_section] = NULL;
      numabs                   [i_section] = NULL;
      parent_num               [i_section] = NULL;
      sparent_entitity_ln_to_gn[i_section] = NULL;
    }
  }

  /*
   * Preparation des groupes
   *  On reverse l'information pour avoir pour chaque elemts le group qui le constitue
   */
  int         **pgroup_idx      = NULL;
  int         **pgroup          = NULL;
  PDM_g_num_t **pgroup_ln_to_gn = NULL;
  if(dmne->n_group_elmt > 0) {
    PDM_part_distgroup_to_partgroup(dmne->comm,
                                    NULL,
                                    dmne->n_group_elmt,
                                    dmne->dgroup_elmt_idx,
                                    dmne->dgroup_elmt,
                                    n_part,
                                    pn_elmt,
            (const PDM_g_num_t **)  section_elmts_ln_to_gn,
                                    &pgroup_idx,
                                    &pgroup,
                                    &pgroup_ln_to_gn);

    PDM_part_mesh_nodal_elmts_n_group_set(pmne, dmne->n_group_elmt, PDM_OWNERSHIP_KEEP);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i_group = 0; i_group < dmne->n_group_elmt; ++i_group) {

        int beg          = pgroup_idx[i_part][i_group];
        int n_group_elmt = pgroup_idx[i_part][i_group+1] - beg;
        int         *group_elmt     = malloc(n_group_elmt * sizeof(int        ));
        PDM_g_num_t *group_ln_to_gn = malloc(n_group_elmt * sizeof(PDM_g_num_t));

        for(int i = 0; i < n_group_elmt; ++i) {
          group_elmt    [i] = pgroup         [i_part][beg+i];
          group_ln_to_gn[i] = pgroup_ln_to_gn[i_part][beg+i];
        }

        PDM_part_mesh_nodal_elmts_group_set(pmne,
                                            i_part,
                                            i_group,
                                            n_group_elmt,
                                            group_elmt,
                                            group_ln_to_gn);
      }

      free(pgroup_idx     [i_part]);
      free(pgroup         [i_part]);
      free(pgroup_ln_to_gn[i_part]);

    }

    free(pgroup_idx);
    free(pgroup);
    free(pgroup_ln_to_gn);

  }


  free(pelmt_by_section_n       );
  free(connec                   );
  free(connec_idx               );
  free(parent_num               );
  free(numabs                   );
  free(sparent_entitity_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pelmts_connec         [i_part]);
    free(pelmts_stride         [i_part]);
    free(pelmts_types          [i_part]);
    free(pelmts_stride_idx     [i_part]);
    free(sorted_vtx_ln_to_gn   [i_part]);
    free(vtx_order             [i_part]);
    free(section_elmts_ln_to_gn[i_part]);
  }
  free(sorted_vtx_ln_to_gn);
  free(vtx_order);
  free(pelmts_connec);
  free(pelmts_stride);
  free(pelmts_types );
  free(pid_section  );
  free(section_order);
  free(section_ho_ordering);
  free(pelmts_stride_idx);
  free(section_elmts_ln_to_gn);

  return pmne;
}



void
PDM_reverse_dparent_gnum
(
       PDM_g_num_t    *dparent_gnum,
       int            *dparent_sign,
       PDM_g_num_t    *delmt_child_distrib,
       int             n_part,
       int            *pn_parent,
       PDM_g_num_t   **pparent_gnum,
       int           **pn_child,
       int          ***pelmt_to_entity,
       PDM_g_num_t  ***pchild_gnum,
       PDM_g_num_t  ***pchild_parent_gnum,
       int          ***pchild_parent_sign,
 const PDM_MPI_Comm    comm
)
{
  // PDM_UNUSED(parent_distrib); // TO remove if all test is OK
  // TODO : pchild_parent_gnum Keep or remove ?

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* First part_to_block to map in parent block the child g_num */
  int dn_child_elmt = delmt_child_distrib[i_rank+1] - delmt_child_distrib[i_rank];

  // Take both parent and child into account to compute suitable distribution
  int          *n_elts_both  = malloc((n_part+1) * sizeof(int));
  PDM_g_num_t **lngn_both    = malloc((n_part+1) * sizeof(PDM_g_num_t*));
  double      **weights_both = malloc((n_part+1) * sizeof(double*));
  //Parent
  for (int i_part=0; i_part < n_part; i_part++){
    n_elts_both[i_part] = pn_parent[i_part];
    lngn_both[i_part] = pparent_gnum[i_part];
    weights_both[i_part] = PDM_array_const_double(n_elts_both[i_part], 1.0);
  }
  //Child
  n_elts_both[n_part] = dn_child_elmt;
  lngn_both[n_part] = dparent_gnum;
  weights_both[n_part] = PDM_array_const_double(n_elts_both[n_part], 1.0);

  PDM_g_num_t* distrib = NULL;
  PDM_distrib_weight (2,
                      n_rank,
                      n_part+1,
                      n_elts_both,
(const PDM_g_num_t**) lngn_both,
     (const double**) weights_both,
                      5,
                      0.1,
                      comm,
                      &distrib);
  for (int i_part = 0; i_part < n_part+1; ++i_part) {
    free(weights_both[i_part]);
  }
  free(weights_both);
  free(lngn_both);
  free(n_elts_both);

  PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                                                   1.,
                                                                   &dparent_gnum,
                                                                   distrib,
                                                                   &dn_child_elmt,
                                                                   1,
                                                                   comm);


  int         *pblk_child_n    = (int         *) malloc( dn_child_elmt * sizeof(int        ));
  PDM_g_num_t *pblk_child_gnum = (PDM_g_num_t *) malloc( dn_child_elmt * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_child_elmt; ++i) {
    pblk_child_n   [i] = 1;
    pblk_child_gnum[i] = delmt_child_distrib[i_rank] + i + 1; // Donc le gnum de child ...
  }

  int         *blk_child_n = NULL;
  PDM_g_num_t *blk_child_gnum   = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                        &pblk_child_n,
             (void **)  &pblk_child_gnum,
                        &blk_child_n,
             (void **)  &blk_child_gnum);

  int *blk_dparent_sign_n = NULL;
  int *blk_dparent_sign = NULL;
  if(dparent_sign != NULL) {
    PDM_part_to_block_exch(ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           &pblk_child_n,
                (void **)  &dparent_sign,
                           &blk_dparent_sign_n,
                (void **)  &blk_dparent_sign);
    free(blk_dparent_sign_n);
  }

  // PDM_g_num_t n_g_parent = parent_distrib[n_rank]+1;
  // PDM_g_num_t* block_distrib_tmp_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb, &blk_child_n, n_g_parent);

  int dn_parent                 = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t* blk_dparent_gnum = PDM_part_to_block_block_gnum_get(ptb);

  free(pblk_child_n   );
  free(pblk_child_gnum);

  /*
   * At this stage we have in each block of parent the global number of child
   * We need to resend into part
   */
  // PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_tmp_idx,
  //                             (const PDM_g_num_t **)  pparent_gnum,
  //                                                     pn_parent,
  //                                                     n_part,
  //                                                     comm);
  

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block_and_distrib(distrib,
                                                                                    blk_dparent_gnum,
                                                                                    dn_parent,
                                                            (const PDM_g_num_t **)  pparent_gnum,
                                                                                    pn_parent,
                                                                                    n_part,
                                                                                    comm);



  int         **_pchild_n    = NULL;
  PDM_g_num_t **_pchild_gnum = NULL;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          blk_child_n,
                          blk_child_gnum,
                         &_pchild_n,
              (void ***) &_pchild_gnum);

  int         **_tmp_pchild_n    = NULL;
  PDM_g_num_t **_tmp_pchild_parent_gnum = NULL;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          blk_child_n,
                          blk_dparent_gnum,
                         &_tmp_pchild_n,
              (void ***) &_tmp_pchild_parent_gnum);

  int **_tmp_pchild_parent_sign = NULL;
  if(blk_dparent_sign != NULL) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(_tmp_pchild_n[i_part]);
    }
    free(_tmp_pchild_n);
    PDM_block_to_part_exch(btp,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            blk_child_n,
                            blk_dparent_sign,
                           &_tmp_pchild_n,
                (void ***) &_tmp_pchild_parent_sign);
    free(blk_dparent_sign);
  }

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free(distrib);
  free(blk_child_n);
  free(blk_child_gnum);
  // free(block_distrib_tmp_idx);

  /*
   * At this stage we have for each partition the number AND the gnum of childs inside
   *          -> We sort pchild_gnum but normally it's unique
   *
   */
  *pn_child = (int *) malloc(n_part * sizeof(int));
  int* _pn_child = *pn_child;
  PDM_g_num_t **_pchild_parent_gnum = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  int **_pelmt_to_entity = (int **) malloc(n_part * sizeof(int *));
  int **_pchild_parent_sign = NULL;
  if(blk_dparent_sign != NULL) {
    _pchild_parent_sign = (int **) malloc(n_part * sizeof(int *));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int pn_child_tmp = 0;
    for(int i = 0; i < pn_parent[i_part]; ++i) {
      // printf("_pchild_n[i_part][%i] = %i \n",i, _pchild_n[i_part][i] );
      pn_child_tmp += _pchild_n[i_part][i];
      assert(_pchild_n[i_part][i] <= 1); // DOnc soit 0 soit 1
    }

    int* _tmp_pelmt_to_entity = malloc(pn_child_tmp * sizeof(int));
    pn_child_tmp = 0;
    for(int i = 0; i < pn_parent[i_part]; ++i) {
    if(_pchild_n[i_part][i] == 1) {
        _tmp_pelmt_to_entity[pn_child_tmp++] = i;
      }
    }

    // PDM_log_trace_array_long(_pchild_gnum[i_part], pn_child_tmp, "_pchild_gnum :: ");

    int* unique_order = (int *) malloc( pn_child_tmp * sizeof(int));

    _pn_child   [i_part] = PDM_inplace_unique_long2(_pchild_gnum[i_part], unique_order, 0, pn_child_tmp-1);
    _pchild_gnum[i_part] = realloc(_pchild_gnum[i_part], _pn_child[i_part] * sizeof(PDM_g_num_t));

    _pchild_parent_gnum[i_part] = (PDM_g_num_t *) malloc( pn_child_tmp * sizeof(PDM_g_num_t));
    _pelmt_to_entity[i_part]    = (int         *) malloc( pn_child_tmp * sizeof(int));
    for(int i = 0; i < _pn_child[i_part]; ++i) {
      int idx_order = unique_order[i];
      PDM_g_num_t gnum = _tmp_pchild_parent_gnum[i_part][i];
      _pchild_parent_gnum[i_part][idx_order] = gnum;
      _pelmt_to_entity[i_part][idx_order] = _tmp_pelmt_to_entity[i];
    }

    if(blk_dparent_sign != NULL) {
      _pchild_parent_sign[i_part] = (int         *) malloc( pn_child_tmp * sizeof(int        ));
      for(int i = 0; i < _pn_child[i_part]; ++i) {
        int idx_order = unique_order[i];
        int sgn = _tmp_pchild_parent_sign[i_part][idx_order];
        _pchild_parent_sign[i_part][i] = sgn;
      }
      free(_tmp_pchild_parent_sign[i_part]);
    }
    free(_tmp_pchild_parent_gnum[i_part]);
    free(_tmp_pelmt_to_entity);

    free(unique_order);
    free(_pchild_n[i_part]);
    free(_tmp_pchild_n[i_part]);
  }
  free(_pchild_n);
  free(_tmp_pchild_n);
  free(_tmp_pchild_parent_gnum);

  *pelmt_to_entity   = _pelmt_to_entity;
  *pchild_gnum        = _pchild_gnum;
  *pchild_parent_gnum = _pchild_parent_gnum;
  if(dparent_sign != NULL) {
    *pchild_parent_sign = _pchild_parent_sign;
    free(_tmp_pchild_parent_sign);
  }
}

void
PDM_generate_ho_vtx_ln_to_gn
(
 PDM_dmesh_nodal_t      *dmn,
 int                     n_part,
 int                    *pn_cell,
 PDM_g_num_t           **pcell_ln_to_gn,
 int                    *pn_face,
 PDM_g_num_t           **pface_ln_to_gn,
 int                    *pn_edge,
 PDM_g_num_t           **pedge_ln_to_gn,
 int                    *pn_vtx,
 PDM_g_num_t           **pvtx_ln_to_gn,
 int                   **pn_vtx_all,
 PDM_g_num_t          ***vtx_all_ln_to_gn
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /* Volumic */
  PDM_g_num_t **pelmt_volumic_connec  = NULL;
  PDM_g_num_t **pelmt_surfacic_connec = NULL;
  PDM_g_num_t **pelmt_ridge_connec    = NULL;

  int *s_connec_volumic  = malloc(n_part * sizeof(int));
  int *s_connec_surfacic = malloc(n_part * sizeof(int));
  int *s_connec_ridge    = malloc(n_part * sizeof(int));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    s_connec_volumic [i_part] = 0;
    s_connec_surfacic[i_part] = 0;
    s_connec_ridge   [i_part] = 0;
  }

  // if(dmn->volumic != NULL) {
  if(dmn->volumic->n_section > 0) {

    int                  **pelmt_volumic_strid_idx = NULL;
    PDM_Mesh_nodal_elt_t **pelmt_volumic_type      = NULL;
    _delmt_vtx_to_pelmt_vtx(dmn->volumic,
                            n_part,
                            pn_cell,
                            pcell_ln_to_gn,
                            &pelmt_volumic_strid_idx,
                            &pelmt_volumic_connec,
                            &pelmt_volumic_type);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      s_connec_volumic[i_part] = pelmt_volumic_strid_idx[i_part][pn_cell[i_part]];
      free(pelmt_volumic_strid_idx[i_part]);
      free(pelmt_volumic_type     [i_part]);
    }
    free(pelmt_volumic_strid_idx);
    free(pelmt_volumic_type);
  }

  // if(dmn->surfacic != NULL) {
  if(dmn->surfacic->n_section > 0) {
    /* Translate face information into elmt information first */
    int          *pn_surf             = NULL;
    int         **psurf_to_entity     = NULL;
    PDM_g_num_t **psurf_gnum          = NULL;
    PDM_g_num_t **psurf_to_face_g_num = NULL;

    if (dmn->surfacic->delmt_child_distrib != NULL) {
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

      int                  **pelmt_surfacic_strid_idx = NULL;
      PDM_Mesh_nodal_elt_t **pelmt_surfacic_type      = NULL;
      _delmt_vtx_to_pelmt_vtx(dmn->surfacic,
                              n_part,
                              pn_surf,
                              psurf_gnum,
                              &pelmt_surfacic_strid_idx,
                              &pelmt_surfacic_connec,
                              &pelmt_surfacic_type);
      for(int i_part = 0; i_part < n_part; ++i_part) {
        s_connec_surfacic[i_part] = pelmt_surfacic_strid_idx[i_part][pn_surf[i_part]];
        free(pelmt_surfacic_strid_idx[i_part]);
        free(pelmt_surfacic_type     [i_part]);
        free(psurf_gnum              [i_part]);
        free(psurf_to_face_g_num     [i_part]);
        free(psurf_to_entity         [i_part]);
      }
      free(psurf_to_entity);
      free(pelmt_surfacic_strid_idx);
      free(pelmt_surfacic_type);
      free(psurf_gnum);
      free(psurf_to_face_g_num);
      free(pn_surf);
    }
    else {
      int                  **pelmt_surfacic_strid_idx = NULL;
      PDM_Mesh_nodal_elt_t **pelmt_surfacic_type      = NULL;
      _delmt_vtx_to_pelmt_vtx(dmn->surfacic,
                              n_part,
                              pn_face,
                              pface_ln_to_gn,
                              &pelmt_surfacic_strid_idx,
                              &pelmt_surfacic_connec,
                              &pelmt_surfacic_type);
      for(int i_part = 0; i_part < n_part; ++i_part) {
        s_connec_surfacic[i_part] = pelmt_surfacic_strid_idx[i_part][pn_face[i_part]];
        free(pelmt_surfacic_strid_idx[i_part]);
        free(pelmt_surfacic_type     [i_part]);
      }
      free(pelmt_surfacic_strid_idx);
      free(pelmt_surfacic_type);
    }
  }

  // if(dmn->ridge != NULL) {
  if(dmn->ridge->n_section > 0) {
    /* Translate edge information into elmt information first */
    int          *pn_ridge             = NULL;
    int         **pridge_to_entity     = NULL;
    PDM_g_num_t **pridge_gnum          = NULL;
    PDM_g_num_t **pridge_to_edge_g_num = NULL;

    if (dmn->ridge->delmt_child_distrib != NULL) {
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

      int                  **pelmt_ridge_strid_idx = NULL;
      PDM_Mesh_nodal_elt_t **pelmt_ridge_type      = NULL;
      _delmt_vtx_to_pelmt_vtx(dmn->ridge,
                              n_part,
                              pn_ridge,
                              pridge_gnum,
                              &pelmt_ridge_strid_idx,
                              &pelmt_ridge_connec,
                              &pelmt_ridge_type);
      for(int i_part = 0; i_part < n_part; ++i_part) {
        s_connec_ridge[i_part] = pelmt_ridge_strid_idx[i_part][pn_ridge[i_part]];
        free(pelmt_ridge_strid_idx[i_part]);
        free(pelmt_ridge_type     [i_part]);
        free(pridge_gnum          [i_part]);
        free(pridge_to_edge_g_num [i_part]);
        free(pridge_to_entity     [i_part]);
      }
      free(pelmt_ridge_strid_idx);
      free(pelmt_ridge_type);
      free(pridge_gnum);
      free(pridge_to_edge_g_num);
      free(pridge_to_entity);
      free(pn_ridge);
    }
    else {
      int                  **pelmt_ridge_strid_idx = NULL;
      PDM_Mesh_nodal_elt_t **pelmt_ridge_type      = NULL;
      _delmt_vtx_to_pelmt_vtx(dmn->ridge,
                              n_part,
                              pn_edge,
                              pedge_ln_to_gn,
                              &pelmt_ridge_strid_idx,
                              &pelmt_ridge_connec,
                              &pelmt_ridge_type);
      for(int i_part = 0; i_part < n_part; ++i_part) {
        s_connec_ridge[i_part] = pelmt_ridge_strid_idx[i_part][pn_edge[i_part]];
        free(pelmt_ridge_strid_idx[i_part]);
        free(pelmt_ridge_type     [i_part]);
      }
      free(pelmt_ridge_strid_idx);
      free(pelmt_ridge_type);
      free(pridge_gnum);
    }
  }

  /*
   * We have all connectivity of HO elemnt for all geometric kind
   *  -> We need to unify with current vtx
   */
  PDM_g_num_t **all_vtx_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  int          *nall_vtx         = malloc(n_part * sizeof(int          ));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int s_tot = s_connec_volumic[i_part] + s_connec_surfacic[i_part] + s_connec_ridge[i_part];

    PDM_g_num_t *concat_connec = malloc(s_tot * sizeof(PDM_g_num_t));

    for(int i = 0; i < s_connec_volumic[i_part]; ++i){
      concat_connec[i] = pelmt_volumic_connec[i_part][i];
    }
    int shift = s_connec_volumic[i_part];
    for(int i = 0; i < s_connec_surfacic[i_part]; ++i){
      concat_connec[shift+i] = pelmt_surfacic_connec[i_part][i];
    }
    shift += s_connec_surfacic[i_part];
    for(int i = 0; i < s_connec_ridge[i_part]; ++i){
      concat_connec[shift+i] = pelmt_ridge_connec[i_part][i];
    }

    /*
     * Sort also vtx_ln_to_gn
     */
    PDM_g_num_t *sorted_vtx_ln_to_gn = malloc(pn_vtx[i_part] * sizeof(PDM_g_num_t));

    for(int i = 0; i < pn_vtx[i_part]; ++i) {
      sorted_vtx_ln_to_gn[i] = pvtx_ln_to_gn[i_part][i];
    }
    PDM_sort_long(sorted_vtx_ln_to_gn, NULL, pn_vtx[i_part]);

    int s_tot_unique = PDM_inplace_unique_long(concat_connec, NULL, 0, s_tot-1);
    concat_connec = realloc(concat_connec, sizeof(PDM_g_num_t) * s_tot_unique);

    all_vtx_ln_to_gn[i_part] = malloc(s_tot_unique * sizeof(PDM_g_num_t)); // s_tot_unique is a majorant

    int idx_write = 0;
    for(int i = 0; i < pn_vtx[i_part]; ++i) {
      all_vtx_ln_to_gn[i_part][idx_write++] = pvtx_ln_to_gn[i_part][i];
    }

    for(int i = 0; i < s_tot_unique; ++i) {

      PDM_g_num_t gnum = concat_connec[i];
      int pos = PDM_binary_search_long(gnum, sorted_vtx_ln_to_gn, pn_vtx[i_part]);
      if(pos == -1) {
        all_vtx_ln_to_gn[i_part][idx_write++] = gnum;
      }

    }

    nall_vtx[i_part] = idx_write;

    free(sorted_vtx_ln_to_gn);
    free(concat_connec);

    if(0 == 1) {
      PDM_log_trace_array_long(all_vtx_ln_to_gn[i_part], nall_vtx[i_part], "all_vtx_ln_to_gn ::");
    }


  }

  free(s_connec_volumic    );
  free(s_connec_surfacic   );
  free(s_connec_ridge      );
  if(pelmt_volumic_connec != NULL){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      if (pelmt_volumic_connec [i_part] != NULL) {
        free(pelmt_volumic_connec [i_part]);
      }
    }
    free(pelmt_volumic_connec);
  }
  if(pelmt_surfacic_connec != NULL){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      if (pelmt_surfacic_connec[i_part] != NULL) {
        free(pelmt_surfacic_connec[i_part]);
      }
    }
    free(pelmt_surfacic_connec);
  }
  if(pelmt_ridge_connec != NULL){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      if (pelmt_ridge_connec   [i_part] != NULL) {
        free(pelmt_ridge_connec   [i_part]);
      }
    }
    free(pelmt_ridge_connec);
  }


  *pn_vtx_all       = nall_vtx;
  *vtx_all_ln_to_gn = all_vtx_ln_to_gn;

}



PDM_dmesh_nodal_elmts_t*
PDM_dmesh_nodal_elmts_to_extract_dmesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           dn_elmt,
 PDM_g_num_t                  *delmt_selected,
 PDM_g_num_t                 **extract_vtx_distribution,
 PDM_g_num_t                 **extract_parent_vtx_gnum,
 PDM_g_num_t                 **extract_entity_distribution,
 PDM_g_num_t                 **extract_parent_elt_gnum
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmne->comm, &i_rank);
  PDM_MPI_Comm_size(dmne->comm, &n_rank);

  // Sort query gnum by section - Caution we need to keep the link with the parent in same order
  int n_section = dmne->n_section;
  int *dn_elmt_by_section_n = PDM_array_zeros_int(n_section);

  for(int i = 0; i < dn_elmt; ++i) {
    int i_section = PDM_binary_search_gap_long(delmt_selected[i] - 1,
                                               dmne->section_distribution,
                                               n_section + 1);
    dn_elmt_by_section_n[i_section]++;
  }

  // PDM_log_trace_array_int(dn_elmt_by_section_n, n_section, "dn_elmt_by_section_n :")

  int *dn_elmt_by_section_idx = malloc((n_section+1) * sizeof(int));
  dn_elmt_by_section_idx[0] = 0;
  for (int i = 0; i < n_section; i++) {
    dn_elmt_by_section_idx[i+1] = dn_elmt_by_section_idx[i] + dn_elmt_by_section_n[i];
    dn_elmt_by_section_n[i] = 0;
  }

  PDM_g_num_t *delmt_sorted_by_section = malloc(dn_elmt * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_elmt; ++i) {
    int i_section = PDM_binary_search_gap_long(delmt_selected[i] - 1,
                                               dmne->section_distribution,
                                               n_section + 1);
    int idx_write = dn_elmt_by_section_idx[i_section] + dn_elmt_by_section_n[i_section]++;
    delmt_sorted_by_section[idx_write] = delmt_selected[i];
  }

  // PDM_log_trace_array_int(delmt_sorted_by_section, dn_elmt_by_section_idx[n_section], "delmt_sorted_by_section ::");

  /*
   * For each section we requilibrate the data
   */
  PDM_dmesh_nodal_elmts_t* extract_dmne = PDM_DMesh_nodal_elmts_create(dmne->comm,
                                                                       dmne->mesh_dimension,
                                                                       0);

  int          *dn_extract_elmt           = malloc(n_section * sizeof(int          ));
  PDM_g_num_t **extract_gnum_elmt         = malloc(n_section * sizeof(PDM_g_num_t *));
  PDM_g_num_t **dextract_block_elmts_gnum = malloc(n_section * sizeof(PDM_g_num_t *));
  int dn_extract_elmt_tot = 0;

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3, n_section, PDM_FALSE, 1e-3, dmne->comm, PDM_OWNERSHIP_USER);

  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

    int beg                = dn_elmt_by_section_idx[i_section];
    int dn_elmt_in_section = dn_elmt_by_section_idx[i_section+1] - beg;
    PDM_g_num_t* _ldextract_gnum = &delmt_sorted_by_section[beg];

    double *weight = malloc(dn_elmt_in_section * sizeof(double));
    for(int i = 0; i < dn_elmt_in_section; ++i) {
      weight[i] = 1.;
    }

    PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        &_ldextract_gnum,
                                                        &weight,
                                                        &dn_elmt_in_section,
                                                        1,
                                                        dmne->comm);

    free(weight);

    const PDM_g_num_t* distrib_section = PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);

    dn_extract_elmt[i_section] = PDM_part_to_block_n_elt_block_get(ptb);
    PDM_g_num_t* _extract_gnum_elmt = PDM_part_to_block_block_gnum_get(ptb);

    if(0 == 1) {
      PDM_log_trace_array_long(_extract_gnum_elmt, dn_extract_elmt[i_section], "_extract_gnum_elmt ::");
    }

    // Shift by global section
    extract_gnum_elmt[i_section] = malloc(dn_extract_elmt[i_section]*sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_extract_elmt[i_section]; ++i) {
      extract_gnum_elmt[i_section][i] = _extract_gnum_elmt[i] - dmne->section_distribution[i_section];
    }

    if(0 == 1) {
      PDM_log_trace_array_long(extract_gnum_elmt[i_section], dn_extract_elmt[i_section], "extract_gnum_elmt ::");
      PDM_log_trace_array_long(distrib_section, n_rank+1, "distrib_section ::");
    }

    PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_section,
                                (const PDM_g_num_t **)  &extract_gnum_elmt[i_section],
                                                        &dn_extract_elmt[i_section],
                                                        1,
                                                        dmne->comm);
    PDM_part_to_block_free(ptb);

    // Exchange connectivity
    int n_vtx_per_elmt = 0;
    PDM_g_num_t *block_elmts_connec = NULL;

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
        // int n_elt          = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, 1);
        block_elmts_connec = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);
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
        int order = -1;
        const char *ho_ordering = NULL;
        block_elmts_connec = PDM_DMesh_nodal_elmts_section_std_ho_get(dmne,
                                                                      id_section,
                                                                      &order,
                                                                      &ho_ordering);

        n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      case PDM_MESH_NODAL_POLY_3D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
    }

    /*   */
    PDM_g_num_t** tmp_dextract_block_elmts_gnum = NULL;
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &n_vtx_per_elmt,
               (void **)   block_elmts_connec,
                           NULL,
               (void ***) &tmp_dextract_block_elmts_gnum);
    dextract_block_elmts_gnum[i_section] = tmp_dextract_block_elmts_gnum[0];
    free(tmp_dextract_block_elmts_gnum);

    PDM_block_to_part_free(btp);

    // PDM_log_trace_array_long(dextract_block_elmts_gnum[i_section], n_vtx_per_elmt * dn_extract_elmt[i_section], "dextract_block_elmts_gnum:: ");

    PDM_gnum_set_from_parents (gen_gnum,
                               i_section,
                               n_vtx_per_elmt * dn_extract_elmt[i_section],
                               dextract_block_elmts_gnum[i_section]);

    for(int i = 0; i < dn_extract_elmt[i_section]; ++i) {
      extract_gnum_elmt[i_section][i] += dmne->section_distribution[i_section];
    }
    dn_extract_elmt_tot += dn_extract_elmt[i_section];
  }

  *extract_parent_elt_gnum = malloc(dn_extract_elmt_tot*sizeof(PDM_g_num_t));
  PDM_g_num_t *lextract_parent_elt_gnum = *extract_parent_elt_gnum;
  int idx_write = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {
    for(int i = 0; i < dn_extract_elmt[i_section]; ++i) {
      lextract_parent_elt_gnum[idx_write++] = extract_gnum_elmt[i_section][i];
    }
    free(extract_gnum_elmt[i_section]);
  }
  free(extract_gnum_elmt);

  /* Generate child gnum */
  PDM_gnum_compute (gen_gnum);

  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);
    int n_vtx_per_elmt = 0;

    PDM_g_num_t *child_section_block_elmts_gnum = PDM_gnum_get(gen_gnum, i_section);

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
        int id_extract_section = PDM_DMesh_nodal_elmts_section_add(extract_dmne, t_elt);
        n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, 1);
        PDM_DMesh_nodal_elmts_section_std_set(extract_dmne,
                                              id_extract_section,
                                              dn_extract_elmt[i_section],
                                              child_section_block_elmts_gnum,
                                              PDM_OWNERSHIP_KEEP);
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
        int order = -1;
        const char *ho_ordering = NULL;
        PDM_g_num_t* lblock_elmts_connec = PDM_DMesh_nodal_elmts_section_std_ho_get(dmne,
                                                                      id_section,
                                                                      &order,
                                                                      &ho_ordering);
        PDM_UNUSED(lblock_elmts_connec);
        n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);

        int id_extract_section = PDM_DMesh_nodal_elmts_section_ho_add(extract_dmne, t_elt, order, ho_ordering);
        PDM_DMesh_nodal_elmts_section_std_set(extract_dmne,
                                              id_extract_section,
                                              dn_extract_elmt[i_section],
                                              child_section_block_elmts_gnum,
                                              PDM_OWNERSHIP_KEEP);

        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      case PDM_MESH_NODAL_POLY_3D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not supported\n");
    }

    dn_extract_elmt[i_section] = dn_extract_elmt[i_section] * n_vtx_per_elmt;
  }

  PDM_gnum_free(gen_gnum);

  /*
   * Generate the new child gnum
   */
  PDM_part_to_block_t* ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          dextract_block_elmts_gnum,
                                                          NULL,
                                                          dn_extract_elmt,
                                                          n_section,
                                                          dmne->comm);

  int          dn_extact_vtx         = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  PDM_g_num_t* ptb_dextract_vtx_gnum = PDM_part_to_block_block_gnum_get (ptb_vtx);

  *extract_parent_vtx_gnum = malloc(dn_extact_vtx * sizeof(PDM_g_num_t));
  PDM_g_num_t *lextract_parent_vtx_gnum = *extract_parent_vtx_gnum;
  for(int i = 0; i < dn_extact_vtx; ++i) {
    lextract_parent_vtx_gnum[i] = ptb_dextract_vtx_gnum[i];
  }

  // PDM_log_trace_array_long(lextract_parent_vtx_gnum, dn_extact_vtx, "lextract_parent_vtx_gnum ::");

  *extract_vtx_distribution = PDM_compute_entity_distribution(dmne->comm, dn_extact_vtx);
  *extract_entity_distribution = PDM_compute_entity_distribution(dmne->comm, dn_extract_elmt_tot);

  PDM_part_to_block_free(ptb_vtx);

  for(int i = 0; i < n_section; ++i) {
    free(dextract_block_elmts_gnum[i]);
  }

  free(delmt_sorted_by_section);
  free(dn_elmt_by_section_idx);
  free(dn_elmt_by_section_n);
  free(dextract_block_elmts_gnum);
  free(dn_extract_elmt);

  return extract_dmne;
}
