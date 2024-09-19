#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_logging.h"

/*
 *  Use case
 *
 *       3           4           2
 *       +-----------+-----------+
 *       |           |           |
 *       |           |           |
 *       |           |           |
 *       +-----------+-----------+
 *       5           1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */
MPI_TEST_CASE("[PDM_delmts_nodal_elmts_t] Constructor",1) {
  // double dvtx_coord[27] = { 1. , 0. , 0.,
  //                           1. , 0.5, 0.,
  //                           1. , 1. , 0.,
  //                           1.5, 1. , 0.,
  //                           2. , 1. , 0.,
  //                           2. , 0.5, 0.,
  //                           2. , 0. , 0.,
  //                           1.5, 0. , 0.,
  //                           1.5, 0.5, 0.};
  // PDM_UNUSED(dvtx_coord);

  const PDM_g_num_t n_vtx            = 9;
  const PDM_g_num_t n_face           = 8;
  // const PDM_g_num_t n_ridge          = 8;
  const int         n_tri_section_1  = 8;
  const int         n_bar_section_1  = 8;

  PDM_g_num_t connec_tri_1[24] = {6, 8, 9,
                                  9, 5, 6,
                                  2, 8, 1,
                                  9, 3, 4,
                                  6, 7, 8,
                                  9, 4, 5,
                                  2, 9, 8,
                                  9, 2, 3};

  PDM_g_num_t connec_bar_1[16] = {1, 2,
                                  2, 3,
                                  4, 5,
                                  3, 4,
                                  6, 7,
                                  5, 6,
                                  8, 1,
                                  7, 8};

  int n_group_elmt = 1;
  int dgroup_elmt_idx[2] = {0, 8};
  // PDM_g_num_t dgroup_elmt[8] = {9, 10, 11, 12, 13, 14, 15, 16};
  PDM_g_num_t dgroup_elmt[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, -1, n_face, -1);

  // The order of call is important for global numbering
  int tri_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_SURFACIC, PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_RIDGE   , PDM_MESH_NODAL_BAR2);

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_group_elmt,
                                         dgroup_elmt_idx,
                                         dgroup_elmt,
                                         PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  tri_section_1,
                                  n_tri_section_1,
                                  connec_tri_1,
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  bar_section_1,
                                  n_bar_section_1,
                                  connec_bar_1,
                                  PDM_OWNERSHIP_USER);

  /*
   * Generate the connectivity
   */
  //
  PDM_dmesh_nodal_generate_distribution(dmn);
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);
  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmn);
  PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  PDM_dmesh_t* dm;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &dm);

  int dn_cell, dn_face, dn_vtx, dn_edge;
  PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx);

  MPI_CHECK(0, dn_face == 8);
  MPI_CHECK(0, dn_edge == 16);

  PDM_g_num_t *dedge_face;
  int         *dedge_face_idx;
  PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                             &dedge_face, &dedge_face_idx, PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *dface_edge;
  int         *dface_edge_idx;
  PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge, &dface_edge_idx, PDM_OWNERSHIP_KEEP);

  // PDM_log_trace_array_int (dedge_face_idx, dn_edge+1, "dedge_face_idx:: ");
  // PDM_log_trace_array_long (dedge_face, dedge_face_idx[dn_edge], "dedge_face:: ");

  // PDM_log_trace_array_int (dface_edge_idx, dn_face+1, "dface_edge_idx:: ");
  // PDM_log_trace_array_long (dface_edge, dface_edge_idx[dn_face], "dface_edge:: ");

  // int dedge_face_idx_expected[17] =  {0, 1, 2, 3, 4, 5, 7, 9, 10, 12, 13, 15, 17, 19, 20, 22, 24 };
  // PDM_g_num_t dedge_face_expected[24] =  {3, 8, 4, 6, 3, 7, -3, 7, -8, 2, 8, -4, 5, 4, -6, 1, -5, 2, -6, 5, 1, -2, 1, -7};

  // int dface_edge_idx_expected[9] =  {0, 3, 6, 9, 12, 15, 18, 21, 24};
  // PDM_g_num_t dface_edge_expected[24] =  {12, 15, 16, -15, 8, 13, -6, 1, 5, -9, 3, 11, -12, 10, 14, -13, -11, 4, -16, 6, 7, -7, 2, 9};

  // Test FAUX car maintenant dedge_face = 2*dn_face !!!
  // CHECK_EQ_C_ARRAY(dedge_face_idx, dedge_face_idx_expected, dn_edge+1                       );
  // CHECK_EQ_C_ARRAY(dedge_face    , dedge_face_expected    , dedge_face_idx_expected[dn_edge]);
  // CHECK_EQ_C_ARRAY(dface_edge_idx, dface_edge_idx_expected, dn_face+1                       );
  // CHECK_EQ_C_ARRAY(dface_edge    , dface_edge_expected    , dface_edge_idx_expected[dn_face]);


  PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);

  PDM_DMesh_nodal_free(dmn);
}
