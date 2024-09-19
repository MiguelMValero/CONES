#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_logging.h"

// n_vtx = 18
// double coord_x[n_vtx] = {0.000499536, 0.000499536, -1.65387e-06, -1.65387e-06, 0.000499536, 0.000499536, -1.65387e-06, -1.65387e-06, 0.00100073, 0.00100073, 0.00100073, 0.00100073, 0.000499536, -1.65387e-06, 0.000499536, -1.65387e-06, 0.00100073, 0.00100073};
// double coord_y[n_vtx] = {0.000498807, -2.38663e-06, -2.38663e-06, 0.000498807, 0.000501645, 4.5126e-07, 4.5126e-07, 0.000501645, 0.000498807, -2.38663e-06, 0.000501645, 4.5126e-07, 0.001, 0.001, 0.00100284, 0.00100284, 0.001, 0.00100284};
// double coord_z[n_vtx] = {0.000999549, 0.000999549, 0.000999549, 0.000999549, 0., 0., 0., 0., 0.000999549, 0.000999549, 0., 0., 0.000999549, 0.000999549, 0., 0., 0.000999549, 0.};

// Cas simple EMMA : /stck2/stck2.3/bmaugars/dev/dev-Tools/maia/unit_tests_case/EMMA/cube_simple/Cube_ANSA_hexa_separated.cgns

MPI_TEST_CASE("[PDM_dmesh_nodal_to_dmesh] decomposes hexa ",1) {

  const PDM_g_num_t n_vtx            = 18;
  const PDM_g_num_t n_cell           = 4;
  const PDM_g_num_t n_face           = 16;
  const int         n_hexa_section_1 = 4;
  const int         n_quad_section_1 = 16;
  PDM_g_num_t connec_hexa_1[32] = {1,2,3,4,5,6,7,8,
                                   9,10,2,1,11,12,6,5,
                                   13,1,4,14,15,5,8,16,
                                   17,9,1,13,18,11,5,15};

  PDM_g_num_t connec_quad_1[64] = {6, 5, 8, 7,
                                   12, 11, 5, 6,
                                   5, 15, 16, 8,
                                   11, 18, 15, 5,
                                   15, 13, 14, 16,
                                   13, 15, 18, 17,
                                   1, 2, 3, 4,
                                   9, 10, 2, 1,
                                   13, 1, 4, 14,
                                   17, 9, 1, 13,
                                   2, 6, 7, 3,
                                   6, 2, 10, 12,
                                   8, 4, 3, 7,
                                   4, 8, 16, 14,
                                   9, 11, 12, 10,
                                   11, 9, 17, 18};

  int n_group_elmt = 6;
  int dgroup_elmt_idx[7] = {0, 4,      // Bottom_1
                            6,         // Left_1
                            10,        // Top_1
                            12,          // Right_1
                            14,          // Inlet_1
                            16};         // Outlet_1
  PDM_g_num_t dgroup_elmt[16] = {1, 2, 3, 4,      // Bottom_1
                                 5, 6,           // Left_1
                                 7, 8, 9, 10,  // Top_1
                                 11, 12,          // Right_1
                                 13, 14,          // Inlet_1
                                 15, 16};         // Outlet_1

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, n_cell, n_face, -1);

  // The order of call is important for global numbering
  int hexa_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_VOLUMIC , PDM_MESH_NODAL_HEXA8);
  int quad_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_SURFACIC, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_VOLUMIC,
                                  hexa_section_1,
                                  n_hexa_section_1,
                                  connec_hexa_1,
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  quad_section_1,
                                  n_quad_section_1,
                                  connec_quad_1,
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_SURFACIC,
                                         n_group_elmt,
                                         dgroup_elmt_idx,
                                         dgroup_elmt,
                                         PDM_OWNERSHIP_USER);

  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  SUBCASE("transform to faces ")
  {
    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

    PDM_dmesh_t* dm;
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dm);

    int dn_cell, dn_face, dn_vtx, dn_edge;
    PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx);

    CHECK( dn_cell == 4 );
    CHECK( dn_face == 20);
    // CHECK( dn_vtx  == 18);

    PDM_g_num_t *dface_cell;
    int         *dface_cell_idx;
    PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               &dface_cell, &dface_cell_idx, PDM_OWNERSHIP_KEEP);

    // PDM_log_trace_array_long(dface_cell, 2*dn_face, "dface_cell:: ");

    PDM_g_num_t dface_cell_expected[40] = {1,0,1,2,1,0,1,3,1,0,2,0,1,0,2,4,2,0,3,0,
                                           2,0,3,4,4,0,2,0,3,0,3,0,4,0,4,0,3,0,4,0};
    CHECK_EQ_C_ARRAY(dface_cell   , dface_cell_expected   , 2*dn_face             );

    PDM_g_num_t *dcell_face;
    int         *dcell_face_idx;
    PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               &dcell_face, &dcell_face_idx, PDM_OWNERSHIP_KEEP);

    int         dcell_face_idx_expected[7] = {0, 6, 12, 18, 24};
    PDM_g_num_t dcell_face_expected[24]    = {1, 2, 3, 4, 5, 7, -2, 6, 8, 9, 11, 14, -4, 10, 12, 15, 16, 19, -12, -8, 13, 17, 18, 20};

    if(0 == 1) {
      PDM_log_trace_array_int (dcell_face_idx, dn_cell+1, "dcell_face_idx:: ");
      PDM_log_trace_array_long(dcell_face, dcell_face_idx[dn_cell], "dcell_face:: ");
    }

    CHECK_EQ_C_ARRAY(dcell_face_idx, dcell_face_idx_expected, dn_cell+1              );
    CHECK_EQ_C_ARRAY(dcell_face    , dcell_face_expected    , dcell_face_idx[dn_cell]);

  }

  // SUBCASE("transform to edges  ")
  // {
  //   // Plante car on a pas implementer la decomposition en edge des Hexa
  //   PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
  //                                    PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
  //                                    PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
  // }

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_DMesh_nodal_free(dmn);

}


MPI_TEST_CASE("[PDM_dmesh_nodal_to_dmesh] decomposes tri ",1) {
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
  PDM_g_num_t dgroup_elmt[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, -1, n_face, -1);

  // The order of call is important for global numbering
  int tri_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_SURFACIC, PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_RIDGE   , PDM_MESH_NODAL_BAR2);

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

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_group_elmt,
                                         dgroup_elmt_idx,
                                         dgroup_elmt,
                                         PDM_OWNERSHIP_USER);

  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);


  PDM_dmesh_t* dm;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dm);

  // int         edge_face_idx_expected_p0[9] = {0, 1, 2, 3, 4, 5, 7, 9, 10 };
  // PDM_g_num_t edge_face_expected_p0[32]    = {3, 0, 8, 0, 4, 0, 3, 0, 6, 0, 7, 3, 8, 7, 2, 0,
  //                                             4, 8, 6, 4, 5, 0, 2, 6, 1, 5, 5, 0, 1, 2, 7, 1};

  int dn_cell, dn_face, dn_vtx, dn_edge;
  PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx);

  MPI_CHECK(0, dn_face == 8);

  MPI_CHECK(0, dn_edge == 16);

  PDM_g_num_t *edge_face;
  int         *edge_face_idx;
  PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                             &edge_face, &edge_face_idx, PDM_OWNERSHIP_KEEP);

  // PDM_log_trace_array_long (edge_face, 2 * dn_edge, "edge_face:: ");
  // MPI_CHECK_EQ_C_ARRAY(0, edge_face, edge_face_expected_p0, 2 * dn_edge);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_DMesh_nodal_free(dmn);

}


MPI_TEST_CASE("[PDM_dmesh_nodal_to_dmesh] decomposes tri 2p ",2) {
  // double dvtx_coord_p0[15] = { 1. , 0. , 0.,
  //                              1. , 0.5, 0.,
  //                              1. , 1. , 0.,
  //                              1.5, 1. , 0.,
  //                              2. , 1. , 0.};

  // double dvtx_coord_p1[12] = {2. , 0.5, 0.,
  //                             2. , 0. , 0.,
  //                             1.5, 0. , 0.,
  //                             1.5, 0.5, 0.};

  const PDM_g_num_t n_vtx              = 9;
  const PDM_g_num_t n_face             = 8;
  const int         n_tri_section_1[2] = {4, 4};
  const int         n_bar_section_1[2] = {4, 4};

  PDM_g_num_t connec_tri_1[2][12] = {{6, 8, 9,
                                      9, 5, 6,
                                      2, 8, 1,
                                      9, 3, 4},
                                     {6, 7, 8,
                                      9, 4, 5,
                                      2, 9, 8,
                                      9, 2, 3}};

  PDM_g_num_t connec_bar_1[2][8] = {{1, 2,
                                     2, 3,
                                     4, 5,
                                     3, 4},
                                   { 6, 7,
                                     5, 6,
                                     8, 1,
                                     7, 8}};

  int n_group_elmt = 1;
  int dgroup_elmt_idx[2][2] = {{0, 4}, {0, 4}};
  PDM_g_num_t dgroup_elmt[2][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, -1, n_face, -1);

  // The order of call is important for global numbering
  int tri_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_SURFACIC, PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_RIDGE   , PDM_MESH_NODAL_BAR2);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  tri_section_1,
                                  n_tri_section_1[test_rank],
                                  connec_tri_1[test_rank],
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  bar_section_1,
                                  n_bar_section_1[test_rank],
                                  connec_bar_1[test_rank],
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_group_elmt,
                                         dgroup_elmt_idx[test_rank],
                                         dgroup_elmt[test_rank],
                                         PDM_OWNERSHIP_USER);

  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);


  PDM_dmesh_t* dm;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dm);

  PDM_g_num_t edge_face_expected_p0[16]    = {3,0,8,0,4,0,6,0,3,0,3,7,2,0,7,8};

  PDM_g_num_t edge_face_expected_p1[16]    = {4,8,4,6,5,0,1,5,2,6,1,2,5,0,1,7};

  int dn_cell, dn_face, dn_vtx, dn_edge;
  PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx);

  MPI_CHECK(0, dn_face == 4);
  MPI_CHECK(1, dn_face == 4);

  MPI_CHECK(0, dn_edge == 8);
  MPI_CHECK(1, dn_edge == 8);

  PDM_g_num_t *edge_face;
  int         *edge_face_idx;
  PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                             &edge_face, &edge_face_idx, PDM_OWNERSHIP_KEEP);


  // PDM_log_trace_array_long (edge_face, 2 * dn_edge, "edge_face:: ");
  MPI_CHECK_EQ_C_ARRAY(0, edge_face, edge_face_expected_p0, 2 * dn_edge);
  MPI_CHECK_EQ_C_ARRAY(1, edge_face, edge_face_expected_p1, 2 * dn_edge);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_DMesh_nodal_free(dmn);

}


MPI_TEST_CASE("[PDM_dmesh_nodal_to_dmesh] find missing ridges ",2) {
  const PDM_g_num_t n_vtx  = 9;
  const PDM_g_num_t n_face = 4;
  const int         n_qua_section_1[2] = {2, 2};
  const int         n_bar_section_1[2] = {5, 5};

  PDM_g_num_t connec_qua_1[2][8] = {{1, 2, 5, 4,
                                     2, 3, 6, 5},
                                    {4, 5, 8, 7,
                                     5, 6, 9, 8}};

  PDM_g_num_t connec_bar_1[2][10] = {{1, 2,
                                      2, 3,
                                      8, 7,
                                      9, 8,
                                      5,3},
                                     {3, 6,
                                      7, 4,
                                      6, 9,
                                      1, 5,
                                      4, 1}};

  int n_group_elmt = 1;
  int dgroup_elmt_idx[2][2] = {{0, 2}, {0, 2}};
  PDM_g_num_t dgroup_elmt[2][2] = {{1, 2}, {3, 4}};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, -1, n_face, -1);

  // The order of call is important for global numbering
  int qua_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_SURFACIC, PDM_MESH_NODAL_QUAD4);
  int bar_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_GEOMETRY_KIND_RIDGE   , PDM_MESH_NODAL_BAR2);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  qua_section_1,
                                  n_qua_section_1[test_rank],
                                  connec_qua_1[test_rank],
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  bar_section_1,
                                  n_bar_section_1[test_rank],
                                  connec_bar_1[test_rank],
                                  PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_section_group_elmt_set(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_group_elmt,
                                         dgroup_elmt_idx[test_rank],
                                         dgroup_elmt[test_rank],
                                         PDM_OWNERSHIP_USER);

  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  //PDM_g_num_t edge_face_expected_p0[


  PDM_dmesh_t* dm;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dm);

  int dn_cell, dn_face, dn_vtx, dn_edge;
  PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx);

  MPI_CHECK(0, dn_face == 2);
  MPI_CHECK(1, dn_face == 2);

  /*MPI_CHECK(0, dn_edge == 7);
    MPI_CHECK(1, dn_edge == 7);*/

  int         *dedge_vtx_idx;
  PDM_g_num_t *dedge_vtx;
  PDM_dmesh_connectivity_get(dm,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             &dedge_vtx,
                             &dedge_vtx_idx,
                             PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *distrib_edge = NULL;
  PDM_dmesh_distrib_get (dm,
                         PDM_MESH_ENTITY_EDGE,
                         &distrib_edge);

  // PDM_log_trace_array_long (distrib_edge,
  //                           3,
  //                           "distrib_edge : ");

  // PDM_log_trace_connectivity_long (dedge_vtx_idx,
  //                                  dedge_vtx,
  //                                  dn_edge,
  //                                  "dedge_vtx :");

  PDM_g_num_t *distrib_missing_ridge;
  PDM_g_num_t *dmissing_ridge_parent_g_num;
  PDM_dmesh_nodal_to_dmesh_get_missing (dmntodm,
                                        0,
                                        PDM_GEOMETRY_KIND_RIDGE,
                                        &distrib_missing_ridge,
                                        &dmissing_ridge_parent_g_num);

  // int dn_missing_ridge = (int) (distrib_missing_ridge[test_rank+1] -
  //                               distrib_missing_ridge[test_rank  ]);
  // PDM_log_trace_array_long (distrib_missing_ridge,
  //                           3,
  //                           "distrib_missing_ridge : ");

  // PDM_log_trace_array_long (dmissing_ridge_parent_g_num,
  //                           dn_missing_ridge,
  //                           "dmissing_ridge_parent_g_num : ");

  PDM_dmesh_nodal_to_dmesh_get_missing (dmntodm,
                                        0,
                                        PDM_GEOMETRY_KIND_RIDGE,
                                        &distrib_missing_ridge,
                                        &dmissing_ridge_parent_g_num);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_DMesh_nodal_free(dmn);
}
