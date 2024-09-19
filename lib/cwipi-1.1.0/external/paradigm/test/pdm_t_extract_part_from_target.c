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
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_gnum.h"
#include "pdm_part_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_unique.h"
#include "pdm_part_geom.h"
#include "pdm_gnum_location.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"

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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
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

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;

  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtxL,
                        &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  /*
   * Create dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     0, // dn_edge
                                     dn_vtx,
                                     comm);

  PDM_dmesh_vtx_coord_set(dm,
                          dvtx_coord,
                          PDM_OWNERSHIP_USER);


  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             dface_vtx,
                             dface_vtx_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_CELL,
                             dface_cell,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_bound_set(dm,
                      PDM_BOUND_TYPE_FACE,
                      n_face_group,
                      dface_group,
                      dface_group_idx,
                      PDM_OWNERSHIP_USER);

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_dmesh_set(mpart, 0, dm);
  PDM_multipart_compute(mpart);

  PDM_g_num_t* distrib_cell = PDM_compute_entity_distribution(comm, dn_cell);

  /*
   * Get the partition domain
   */
  int i_domain = 0;

  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  int          *pn_cell                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_domains * sizeof(int          ));

  int         **pcell_face              = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_domains * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_domains * sizeof(double      *));

  PDM_g_num_t **target_g_num   = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  int          *pn_target_cell = (int          *) malloc( n_part_domains * sizeof(int          ));

  /*
   * Compute gnum location
   */
  PDM_gnum_location_t* gnum_loc = PDM_gnum_location_create(n_part,
                                                           n_part, comm, PDM_OWNERSHIP_KEEP);

  for (int i_part = 0; i_part < n_part_domains; i_part++){

    PDM_g_num_t* cell_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     i_domain,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);

    int *face_vtx     = NULL;
    int *face_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &face_vtx,
                                        PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* face_ln_to_gn = NULL;
    int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                 i_domain,
                                                 i_part,
                                                 PDM_MESH_ENTITY_FACE,
                                                 &face_ln_to_gn,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);

    double *vtx = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     i_domain,
                                     i_part,
                                     &vtx,
                                     PDM_OWNERSHIP_KEEP);

    pn_cell       [i_part] = n_cell;
    pcell_ln_to_gn[i_part] = cell_ln_to_gn;
    pface_ln_to_gn[i_part] = face_ln_to_gn;
    pvtx_ln_to_gn [i_part] = vtx_ln_to_gn;
    pcell_face    [i_part] = cell_face;
    pcell_face_idx[i_part] = cell_face_idx;
    pn_face       [i_part] = n_face;
    pn_vtx        [i_part] = n_vtx;

    pface_vtx    [i_part] = face_vtx;
    pface_vtx_idx[i_part] = face_vtx_idx;
    pvtx_coord   [i_part] = vtx;

    /* On prends tout les gnum impaire dans un ordre reverse */
    int n_target_cell = 0;
    target_g_num[i_part] = malloc(dn_cell * sizeof(PDM_g_num_t));

    for(int i = 0; i < dn_cell; ++i ) {
      PDM_g_num_t gnum = distrib_cell[i_rank+1] - i;
      if(gnum % 2 == 0) {
        target_g_num[i_part][n_target_cell++] = gnum;
      }
    }

    target_g_num  [i_part] = realloc(target_g_num[i_part], n_target_cell * sizeof(PDM_g_num_t));
    pn_target_cell[i_part] = n_target_cell;

    if(0 == 1) {
      PDM_log_trace_array_long(target_g_num  [i_part], pn_target_cell[i_part], "target_g_num :: ");
    }

    PDM_gnum_location_elements_set(gnum_loc,
                                   i_part,
                                   pn_cell[i_part],
                                   pcell_ln_to_gn[i_part]);
    PDM_gnum_location_requested_elements_set(gnum_loc,
                                             i_part,
                                             pn_target_cell[i_part],
                                             target_g_num  [i_part]);
  }
  free(distrib_cell);

  PDM_gnum_location_compute(gnum_loc);

  /*
   * Extract
   */
  int n_part_out = 1;
  PDM_extract_part_t* extrp = PDM_extract_part_create(3,
                                                      n_part,
                                                      n_part_out,
                                                      PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_TRUE,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell[i_part],
                              pn_face[i_part],
                              -1, // pn_edge[i_part],
                              pn_vtx[i_part],
                              pcell_face_idx[i_part],
                              pcell_face[i_part],
                              NULL, //pface_edge_idx[i_part],
                              NULL, //pface_edge[i_part],
                              NULL, //pedge_vtx[i_part],
                              pface_vtx_idx[i_part],
                              pface_vtx[i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              NULL, //pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn[i_part],
                              pvtx_coord[i_part]);

    int *location_idx = NULL;
    int *location     = NULL;
    PDM_gnum_location_get(gnum_loc, i_part, &location_idx, &location);

    PDM_extract_part_target_set(extrp,
                                i_part,
                                pn_target_cell[i_part],
                                target_g_num  [i_part],
                                location);
    // PDM_extract_part_target_set(extrp,
    //                             i_part,
    //                             pn_target_cell[i_part],
    //                             target_g_num  [i_part],
    //                             NULL);

  }


  PDM_extract_part_compute(extrp);

  int          *pn_extract_face        = malloc(n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = malloc(n_part_out * sizeof(int          ));
  int         **pextract_face_vtx      = malloc(n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = malloc(n_part_out * sizeof(int         *));
  double      **pextract_vtx           = malloc(n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_face_ln_to_gn = malloc(n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = malloc(n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < n_part_out; ++i_part) {

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    // PDM_extract_part_ln_to_gn_get(extrp,
    //                               i_part,
    //                               PDM_MESH_ENTITY_FACE,
    //                               &pextract_face_ln_to_gn[i_part],
    //                               PDM_OWNERSHIP_KEEP);

    // PDM_extract_part_ln_to_gn_get(extrp,
    //                               i_part,
    //                               PDM_MESH_ENTITY_VTX,
    //                               &pextract_vtx_ln_to_gn[i_part],
    //                               PDM_OWNERSHIP_KEEP);
    pextract_face_ln_to_gn[i_part] = NULL;
    pextract_vtx_ln_to_gn [i_part] = NULL;

    // PDM_g_num_t* pextract_parent_cell_ln_to_gn = NULL;
    // int n_cell = PDM_extract_part_parent_ln_to_gn_get(extrp,
    //                                                   i_part,
    //                                                   PDM_MESH_ENTITY_CELL,
    //                                                   &pextract_parent_cell_ln_to_gn,
    //                                                   PDM_OWNERSHIP_KEEP);

    // for(int i = 0; i < n_cell; ++i) {
    //   assert(pextract_parent_cell_ln_to_gn[i] == target_g_num[i_part][i]);
    // }
    // PDM_log_trace_array_long(pextract_parent_cell_ln_to_gn, n_cell, "pextract_parent_cell_ln_to_gn ::");

  }
  PDM_gnum_location_free(gnum_loc);

  PDM_part_to_part_t *ptp_face = NULL;
  PDM_extract_part_part_to_part_get(extrp,
                                    PDM_MESH_ENTITY_FACE,
                                    &ptp_face,
                                    PDM_OWNERSHIP_KEEP);

  // 1) Compute field on origin faces
  int    **part1_stride = (int    **) malloc(sizeof(double *) * n_part_domains);
  double **pface_field  = (double **) malloc(sizeof(double *) * n_part_domains);
  for (int i_part = 0; i_part < n_part_domains; i_part++) {
    pface_field [i_part] = malloc(sizeof(double) * pn_face[i_part]);
    part1_stride[i_part] = PDM_array_const_int(pn_face[i_part], 1);

    double *surface_vector = malloc(sizeof(double) * pn_face[i_part] * 3);
    double *center         = malloc(sizeof(double) * pn_face[i_part] * 3);
    PDM_geom_elem_polygon_properties(pn_face      [i_part],
                                     pface_vtx_idx[i_part],
                                     pface_vtx    [i_part],
                                     pvtx_coord   [i_part],
                                     surface_vector,
                                     center,
                                     NULL,
                                     NULL);

    for (int iface = 0; iface < pn_face[i_part]; iface++) {
      pface_field[i_part][iface] = center[3*iface];
    }

    free(surface_vector);
    free(center);
  }

  // // 1.5) Test
  // int **p1_data = malloc(sizeof(int *) * n_part_out);
  // for (int i_part = 0; i_part < n_part_out; ++i_part) {
  //   p1_data[i_part] = malloc(sizeof(int) * pn_extract_face[i_part]);
  //   for (int iface = 0; iface < pn_extract_face[i_part]; iface++) {
  //     p1_data[i_part][iface] = iface;
  //   }
  // }

  // int **p2_data = NULL;

  // int request2 = -1;
  // PDM_part_to_part_iexch(ptp_face,
  //                        PDM_MPI_COMM_KIND_P2P,
  //                        PDM_STRIDE_CST_INTERLACED,
  //                        PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
  //                        1,
  //                        sizeof(int),
  //                        NULL,
  //       (const void  **) p1_data,
  //                        NULL,
  //       (      void ***) &p2_data,
  //                        &request2);
  // PDM_part_to_part_iexch_wait(ptp_face, request2);


  // 2) Exchange to extracted part
  int request = -1;
  // int **part2_stride = NULL;
  double **pextract_face_field = NULL;
  PDM_part_to_part_reverse_iexch(ptp_face,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,//(const int **) part1_stride,
                (const void  **) pface_field,
                                 NULL,//&part2_stride,
                (      void ***) &pextract_face_field,
                                 &request);

  // 3) Shift field on extracted part
  PDM_part_to_part_reverse_iexch_wait(ptp_face, request);
  // sleep(1);

  int          *n_elt1             = NULL;
  int         **part1_to_part2_idx = NULL;
  PDM_g_num_t **part1_to_part2     = NULL;
  PDM_part_to_part_part1_to_part2_get(ptp_face,
                                      &n_elt1,
                                      &part1_to_part2_idx,
                                      &part1_to_part2);

  for (int i_part = 0; i_part < n_part_out; ++i_part) {
    // PDM_log_trace_array_int(part1_to_part2_idx[i_part],
    //                         pn_extract_face[i_part]+1,
    //                         "part1_to_part2_idx : ");
    // PDM_log_trace_connectivity_long(part1_to_part2_idx[i_part],
    //                                 part1_to_part2[i_part],
    //                                 n_elt1[i_part],
    //                                 "part1_to_part2 : ");
    for (int iface = 0; iface < pn_extract_face[i_part]; iface++) {
      pextract_face_field[i_part][iface] *= 10.;
    }
  }

  // 4) Send back to origin frame
  for (int i_part = 0; i_part < n_part_domains; i_part++) {
    free(part1_stride[i_part]);
  }
  free(part1_stride);
  double **pface_field2 = NULL;
  PDM_part_to_part_iexch(ptp_face,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(double),
                         NULL,//(const int **) part2_stride, //NULL,
        (const void  **) pextract_face_field,
                         NULL,//&part1_stride, //NULL,
        (      void ***) &pface_field2,
                         &request);

  PDM_part_to_part_iexch_wait(ptp_face, request);

  int  *n_ref_face = NULL;
  int **ref_face   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_face,
                                 &n_ref_face,
                                 &ref_face);

  for (int i_part = 0; i_part < n_part_domains; ++i_part) {
    for (int iface = 0; iface < n_ref_face[i_part]; iface++) {
      pface_field[i_part][ref_face[i_part][iface] - 1] = pface_field2[i_part][iface];
    }
    free(pface_field2[i_part]);
  }
  free(pface_field2);




  /*
   * Export vtk en légende
   */
  if(post) {
    for(int i_part = 0; i_part < n_part_domains; ++i_part) {

      char filename[999];
      sprintf(filename, "face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata_field(filename,
                                   pn_vtx        [i_part],
                                   pvtx_coord    [i_part],
                                   pvtx_ln_to_gn [i_part],
                                   pn_face       [i_part],
                                   pface_vtx_idx [i_part],
                                   pface_vtx     [i_part],
                                   pface_ln_to_gn[i_part],
                                   "field",
                  (const double *) pface_field[i_part],
                                   NULL,
                                   NULL);
    }

    for(int i_part = 0; i_part < n_part_out; ++i_part) {

      char filename[999];
      sprintf(filename, "extract_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_extract_vtx[i_part],
                                pextract_vtx[i_part],
                                NULL, NULL);

      // PDM_log_trace_connectivity_int(pextract_face_vtx_idx[i_part],
      //                                pextract_face_vtx    [i_part],
      //                                pn_extract_face[i_part], " pextract_face_vtx :: ");

      sprintf(filename, "extract_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      // PDM_vtk_write_polydata(filename,
      //                        pn_extract_vtx[i_part],
      //                        pextract_vtx[i_part],
      //                        pextract_vtx_ln_to_gn[i_part],
      //                        pn_extract_face[i_part],
      //                        pextract_face_vtx_idx[i_part],
      //                        pextract_face_vtx[i_part],
      //                        pextract_face_ln_to_gn[i_part],
      //                        NULL);
      PDM_vtk_write_polydata_field(filename,
                                   pn_extract_vtx[i_part],
                                   pextract_vtx[i_part],
                                   pextract_vtx_ln_to_gn[i_part],
                                   pn_extract_face[i_part],
                                   pextract_face_vtx_idx[i_part],
                                   pextract_face_vtx[i_part],
                                   pextract_face_ln_to_gn[i_part],
                                   "field",
                  (const double *) pextract_face_field[i_part],
                                   NULL,
                                   NULL);
    }
  }

  for (int i_part = 0; i_part < n_part_domains; i_part++) {
    free(pface_field[i_part]);
  }
  free(pface_field);

  for (int i_part = 0; i_part < n_part_out; i_part++) {
    free(pextract_face_field[i_part]);
  }
  free(pextract_face_field);

  free(pn_extract_face);
  free(pn_extract_vtx);
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );



  PDM_extract_part_free(extrp);


  for (int i_part = 0; i_part < n_part_domains; i_part++){
    free(target_g_num    [i_part]);
  }
  free(pn_target_cell);
  free(target_g_num);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);

  free(pcell_ln_to_gn  );
  free(pface_ln_to_gn  );
  free(pvtx_ln_to_gn  );
  free(pcell_face    );
  free(pcell_face_idx);
  free(pface_vtx     );
  free(pface_vtx_idx );
  free(pvtx_coord    );

  PDM_multipart_free(mpart);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);

  if (i_rank == 0) {
   PDM_printf("-- End\n");
   fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
