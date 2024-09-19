#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multi_block_merge.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_domain_interface.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 4;
  double             length  = 1.;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);
  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  int n_block = 2;

  PDM_dcube_nodal_t* dcube1 = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_HEXA8,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube1);

  PDM_dmesh_nodal_t*  dmn1 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube1);
  /*
   * Define distribution of cell
   */
  PDM_dmesh_nodal_dump_vtk(dmn1, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic_dcube1_");

  /*
   * Define distibution of vtx
   */
  PDM_dcube_nodal_t* dcube2 = PDM_dcube_nodal_gen_create(comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         1.,
                                                         0.,
                                                         0.,
                                                         PDM_MESH_NODAL_HEXA8,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube2);

  PDM_dmesh_nodal_t*  dmn2 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube2);

  PDM_dmesh_nodal_dump_vtk(dmn2, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic_dcube2_");

  // Get interface data
  int           *n_group_elt    = malloc(n_block * sizeof(int));
  int         **dgroup_elmt_idx = malloc(n_block * sizeof(int*));
  PDM_g_num_t **dgroup_elmt     = malloc(n_block * sizeof(PDM_g_num_t*));
  int          *dn_vtx          = malloc(n_block * sizeof(int));
  int          *dn_face         = malloc(n_block * sizeof(int));
  int         **dface_vtx_idx   = malloc(n_block * sizeof(int*));
  PDM_g_num_t **dface_vtx       = malloc(n_block * sizeof(PDM_g_num_t*));


  PDM_dmesh_nodal_t* dmn[] = {dmn1, dmn2};
  for (int i_block = 0; i_block < n_block; i_block++) {
    PDM_DMesh_nodal_section_group_elmt_get(dmn[i_block], PDM_GEOMETRY_KIND_SURFACIC,
        &n_group_elt[i_block], &dgroup_elmt_idx[i_block], &dgroup_elmt[i_block]);

    dn_vtx[i_block]  = PDM_DMesh_nodal_n_vtx_get(dmn[i_block]);
    dn_face[i_block] = PDM_DMesh_nodal_section_n_elt_get(dmn[i_block], PDM_GEOMETRY_KIND_SURFACIC, 0);
    dface_vtx_idx[i_block] = (int *) malloc((dn_face[i_block]+1)*sizeof(int));
    for (int i = 0; i < dn_face[i_block]+1; i++)
      dface_vtx_idx[i_block][i] = 4*i;
    dface_vtx[i_block] = PDM_DMesh_nodal_section_std_get(dmn[i_block], PDM_GEOMETRY_KIND_SURFACIC, 0);

  }
  // LAZY SETUP : we assume that we have only 2 blocks of same size to have same distribution :)
  assert (n_block == 2);
  int jn_size  = dgroup_elmt_idx[0][4] - dgroup_elmt_idx[0][3];
  PDM_g_num_t *interface_face_ids = malloc(2*jn_size*sizeof(PDM_g_num_t));
  int         *interface_face_dom = malloc(2*jn_size*sizeof(int        ));
  for (int i=0; i < jn_size; i++) {
    interface_face_ids[2*i]   = dgroup_elmt[0][dgroup_elmt_idx[0][3]+i];
    interface_face_ids[2*i+1] = dgroup_elmt[1][dgroup_elmt_idx[0][2]+i];
    interface_face_dom[2*i]   = 0;
    interface_face_dom[2*i+1] = 1;
  }

  PDM_domain_interface_t *dom_intrf = PDM_domain_interface_create(
      1, n_block, PDM_DOMAIN_INTERFACE_MULT_YES, PDM_OWNERSHIP_KEEP, comm);
  PDM_domain_interface_set(dom_intrf, PDM_BOUND_TYPE_FACE, &jn_size, &interface_face_ids, &interface_face_dom);
  //Apparement bug en 2d !!
  PDM_domain_interface_translate_face2vtx(dom_intrf, dn_vtx, dn_face, dface_vtx_idx, dface_vtx);


  int         *graph_vtx_idx = NULL;
  PDM_g_num_t *graph_vtx_ids = NULL;
  int         *graph_vtx_dom = NULL;
  int graph_vtx_dn = PDM_domain_interface_get_as_graph(dom_intrf, PDM_BOUND_TYPE_VTX,
      &graph_vtx_idx, &graph_vtx_ids, &graph_vtx_dom);
  PDM_log_trace_array_int(graph_vtx_idx, graph_vtx_dn+1, "vtx graph idx");
  PDM_log_trace_array_long(graph_vtx_ids, graph_vtx_idx[graph_vtx_dn], "vtx graph gnums");
  PDM_log_trace_array_int(graph_vtx_dom, graph_vtx_idx[graph_vtx_dn], "vtx graph dom");


  PDM_domain_interface_free(dom_intrf);

  free(n_group_elt);
  free(dgroup_elmt_idx);
  free(dgroup_elmt);
  free(dn_vtx);
  free(dn_face);
  for (int i_block = 0; i_block < n_block; i_block++) {
    free(dface_vtx_idx[i_block]);
  }
  free(dface_vtx_idx);
  free(dface_vtx);

  free(interface_face_ids);
  free(interface_face_dom);

  /*
   * Concatenate blocks
   */
  PDM_g_num_t** block_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_distrib_idx[0] = PDM_dmesh_nodal_vtx_distrib_get(dmn1);
  block_distrib_idx[1] = PDM_dmesh_nodal_vtx_distrib_get(dmn2);

  int* n_selected = malloc(n_block * sizeof(int));
  n_selected[0] = PDM_DMesh_nodal_n_vtx_get(dmn1);
  n_selected[1] = PDM_DMesh_nodal_n_vtx_get(dmn2);

  PDM_g_num_t** selected_g_num = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    selected_g_num[i_block] = malloc(n_selected[i_block] * sizeof(PDM_g_num_t));
    PDM_log_trace_array_long(block_distrib_idx[i_block], n_rank+1, "block_distrib_idx ::");
    for(int i = 0; i < n_selected[i_block]; ++i) {
      selected_g_num[i_block][i] = block_distrib_idx[i_block][i_rank] + i + 1;
    }
  }


  PDM_multi_block_merge_t* mbm = PDM_multi_block_merge_create(block_distrib_idx,
                                                              n_block,
                                                              n_selected,
                                                              selected_g_num,
                                                              graph_vtx_dn,
                                                              graph_vtx_idx,
                                                              graph_vtx_dom,
                                                              graph_vtx_ids,
                                                              comm);


  // Exchange

  double** dvtx_coord = malloc(n_block * sizeof(double *));
  dvtx_coord[0] = PDM_DMesh_nodal_vtx_get(dmn1);
  dvtx_coord[1] = PDM_DMesh_nodal_vtx_get(dmn2);

  int** stride_one = malloc(n_block * sizeof(int *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    stride_one[i_block] = PDM_array_const_int(1, 1);
  }

  double *dmerge_vtx_coord = NULL;
  PDM_multi_block_merge_exch(mbm,
                             3 * sizeof(double),
                             PDM_STRIDE_CST,
                             stride_one,
                 (void * )   dvtx_coord,
                             NULL,
                 (void **)   &dmerge_vtx_coord);


  //
  //
  //
  // PDM_multi_block_merge_exch(mbm_elt,
  //                            3 * sizeof(double),
  //                            PDM_STRIDE_CST,
  //                            stride_one,
  //                (void * )   dcell_vtx,
  //                            NULL,
  //                (void **)   &dmerge_dcell_vtx);


  // origin_cell get_orgin_block (size= n_dmerge_cell)

  // orgin_vtx = 4 * s_orgini_cell

  // origin = 

  // Creer dans PDM_multi_block_merge une fonction qui applique la nouvelle numerotation
  // à un tableau contenant des références à l'ancienne numerotation  
  //
  // Transformer indication numerotation en doublon / numabs origin
  //
  // PDM_multi_block_merge_apply_array(mbm,
  //                            size_dmerge_dcell_vtx,
  //                            dmerge_vtx_origi_block,
  //                            dmerge_dcell_vtx,
  //                            dmerge_dcell_new_vtx);



  free(dvtx_coord);

  /*
   * Same protocol for cells
   */

  PDM_g_num_t** block_elmt_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_elmt_distrib_idx[0] = (PDM_g_num_t *) PDM_DMesh_nodal_distrib_section_get(dmn1, PDM_GEOMETRY_KIND_VOLUMIC, 0);
  block_elmt_distrib_idx[1] = (PDM_g_num_t *) PDM_DMesh_nodal_distrib_section_get(dmn2, PDM_GEOMETRY_KIND_VOLUMIC, 0);

  int* n_elmt_selected = malloc(n_block * sizeof(int));
  n_elmt_selected[0] = PDM_DMesh_nodal_section_n_elt_get(dmn1, PDM_GEOMETRY_KIND_VOLUMIC, 0);
  n_elmt_selected[1] = PDM_DMesh_nodal_section_n_elt_get(dmn2, PDM_GEOMETRY_KIND_VOLUMIC, 0);

  PDM_g_num_t** selected_elmt_g_num = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    selected_elmt_g_num[i_block] = malloc(n_elmt_selected[i_block] * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_elmt_selected[i_block]; ++i) {
      selected_elmt_g_num[i_block][i] = block_elmt_distrib_idx[i_block][i_rank] + i + 1;
    }
  }

  /*
   * Setup graph
   */
  int         *dmerge_elmt_idx      = malloc(1 * sizeof(int        ));
  int         *dmerge_elmt_block_id = malloc(0 * sizeof(int        ));
  PDM_g_num_t *dmerge_elmt_g_num    = malloc(0 * sizeof(PDM_g_num_t));
  dmerge_elmt_idx[0] = 0;

  PDM_multi_block_merge_t* mbm_elmt = PDM_multi_block_merge_create(block_elmt_distrib_idx,
                                                                   n_block,
                                                                   n_elmt_selected,
                                                                   selected_elmt_g_num,
                                                                   0,
                                                                   dmerge_elmt_idx,
                                                                   dmerge_elmt_block_id,
                                                                   dmerge_elmt_g_num,
                                                                   comm);

  /*
   * Exchange + Update
   */
  PDM_g_num_t** block_elmt_vtx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_elmt_vtx[0] = PDM_DMesh_nodal_section_std_get(dmn1, PDM_GEOMETRY_KIND_VOLUMIC, 0);
  block_elmt_vtx[1] = PDM_DMesh_nodal_section_std_get(dmn2, PDM_GEOMETRY_KIND_VOLUMIC, 0);
  int strid_cst = 8; // Because HEXA

  int         *dmerge_elmt_vtx_stride = NULL;
  PDM_g_num_t *dmerge_elmt_vtx = NULL;
  int **stride_cst_ptr = (int **) malloc(n_block*sizeof(int*));
  for (int ib = 0; ib < n_block; ib++) {
    stride_cst_ptr[ib] = &strid_cst;
  }

  PDM_multi_block_merge_exch_and_update(mbm_elmt,
                                        mbm,
                                        PDM_STRIDE_CST,
                                        stride_cst_ptr,
                                        block_elmt_vtx,
                               (int **) NULL,
                                       &dmerge_elmt_vtx_stride,
                                       &dmerge_elmt_vtx); 
  free(dmerge_elmt_vtx_stride);
  free(stride_cst_ptr);

  /*
   * Visualisation
   */
  int dn_merge_vtx               = PDM_multi_block_merge_get_n_block(mbm);
  PDM_g_num_t* distrib_merge_vtx = PDM_multi_block_merge_get_distrib(mbm);

  int dn_merge_elmt               = PDM_multi_block_merge_get_n_block(mbm_elmt);
  PDM_g_num_t* distrib_merge_elmt = PDM_multi_block_merge_get_distrib(mbm_elmt);

  assert(dn_merge_vtx  == distrib_merge_vtx [i_rank+1] - distrib_merge_vtx [i_rank]);
  assert(dn_merge_elmt == distrib_merge_elmt[i_rank+1] - distrib_merge_elmt[i_rank]);

  printf("dn_merge_elmt = %i \n", dn_merge_elmt);

  PDM_g_num_t *merge_elmt_ln_to_gn = malloc( dn_merge_elmt      * sizeof(PDM_g_num_t));
  int         *dconnec_idx         = malloc((dn_merge_elmt + 1) * sizeof(int        ));

  dconnec_idx[0] = 0;
  for(int i = 0; i < dn_merge_elmt; ++i) {
    merge_elmt_ln_to_gn[i] = distrib_merge_elmt[i_rank] + i + 1;
    dconnec_idx[i+1] = dconnec_idx[i] + 8; // Because HEXA
  }

  PDM_log_trace_connectivity_long(dconnec_idx, dmerge_elmt_vtx, dn_merge_elmt, "dmerge_elmt_vtx :: ");


  PDM_g_num_t *pvtx_ln_to_gn;
  int         *pcell_vtx_idx;
  int         *pcell_vtx;
  int          pn_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_merge_elmt,
                                                           dconnec_idx,
                                                           dmerge_elmt_vtx,
                                                           dn_merge_elmt,
                                  (const PDM_g_num_t *)    merge_elmt_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pcell_vtx_idx,
                                                           &pcell_vtx);

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_merge_vtx,
                                        dmerge_vtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord_out = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);


  char filename_elmt[999];
  sprintf(filename_elmt, "merge_mesh_%2.2d.vtk", i_rank);
  PDM_vtk_write_std_elements(filename_elmt,
                             pn_vtx,
                             pvtx_coord_out,
                             pvtx_ln_to_gn,
                             PDM_MESH_NODAL_HEXA8,
                             dn_merge_elmt,
                             pcell_vtx,
                             merge_elmt_ln_to_gn,
                             0,
                             NULL,
                             NULL);

  PDM_g_num_t* merge_vtx_ln_to_gn = malloc(dn_merge_vtx * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_merge_vtx; ++i) {
    merge_vtx_ln_to_gn[i] = distrib_merge_vtx[i_rank] + i + 1;
  }
  char filename[999];
  sprintf(filename, "debug_dvtx_coord_merge_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            dn_merge_vtx,
                            dmerge_vtx_coord,
                            merge_vtx_ln_to_gn,
                            NULL);

  free(dmerge_vtx_coord);
  free(merge_vtx_ln_to_gn);
  free(merge_elmt_ln_to_gn);
  free(dconnec_idx);
  free(pvtx_ln_to_gn);
  free(pcell_vtx_idx);
  free(pcell_vtx);
  free(pvtx_coord_out);

  PDM_multi_block_merge_free(mbm);
  PDM_multi_block_merge_free(mbm_elmt);

  free(dmerge_elmt_vtx);
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    free(selected_g_num[i_block]);
    free(selected_elmt_g_num[i_block]);
  }
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    free(stride_one[i_block]);
  }
  free(stride_one);
  free(block_elmt_vtx);
  free(selected_g_num);
  free(selected_elmt_g_num);
  free(dmerge_elmt_idx     );
  free(dmerge_elmt_block_id);
  free(dmerge_elmt_g_num   );
  free(graph_vtx_idx);
  free(graph_vtx_ids);
  free(graph_vtx_dom);

  free(block_distrib_idx);
  free(block_elmt_distrib_idx);
  free(n_selected);
  free(n_elmt_selected);


  PDM_dcube_nodal_gen_free(dcube1);
  PDM_dcube_nodal_gen_free(dcube2);

  PDM_MPI_Finalize ();
  return 0;
}
