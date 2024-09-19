#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_logging.h"

// double coord_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};


MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes hexa",1) {

  // const PDM_g_num_t n_vtx            = 12;
  const PDM_g_num_t n_cell           = 2;
  const int         n_hexa_section_1 = 2;
  PDM_g_num_t connec_hexa_1[16] = {1, 2, 5, 4, 7, 8, 11, 10, // First
                                   2, 3, 6, 5, 8, 9, 12, 11};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 3, n_cell);

  int hexa_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_HEXA8);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        hexa_section_1,
                                        n_hexa_section_1,
                                        connec_hexa_1,
                                        PDM_OWNERSHIP_USER); // Ownership is just a mapping here
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  // printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 12 );
  CHECK( n_sum_vtx_face_tot == 48 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);
  std::vector<int>         dparent_elmt_pos(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                                  dcell_face_vtx_idx.data(),
                                  dcell_face_vtx.data(),
                                  delmt_face_cell.data(),
                                  NULL, NULL,
                                  dparent_elmt_pos.data());
  // dcell_face_vtx_idx[0] = 0;
  // PDM_sections_decompose_faces(dmn,
  //                                 dcell_face_vtx_idx,
  //                                 dcell_face_vtx,
  //                                 delmt_face_cell,
  //                                 NULL, NULL);

  // PDM_log_trace_array_long(delmt_face_cell, n_face_elt_tot, "delmt_face_cell:: ");
  // PDM_log_trace_array_int(dcell_face_vtx_idx, n_face_elt_tot+1, "dcell_face_vtx_idx:: ");
  // PDM_log_trace_array_long(dcell_face_vtx, n_sum_vtx_face_tot, "dcell_face_vtx:: ");
  // PDM_log_trace_array_long(dparent_elmt_pos.data(), n_face_elt_tot, "dparent_elmt_pos:: ");

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {4, 5, 2, 1, 11, 10, 7, 8,
                                                          7, 10, 4, 1, 10, 11, 5, 4,
                                                          5, 11, 8, 2, 2, 8, 7, 1,
                                                          5, 6, 3, 2, 12, 11, 8, 9,
                                                          8, 11, 5, 2, 11, 12, 6, 5,
                                                          6, 12, 9, 3, 3, 9, 8, 2};
  std::vector<int> dparent_elmt_pos_expected           = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  PDM_DMesh_nodal_elmts_free(dmn);
}


// double coord_tetra_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_tetra_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_tetra_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};

MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes tetra",1) {
  // const PDM_g_num_t n_vtx             = 12;
  const PDM_g_num_t n_cell            = 10;
  const int         n_tetra_section_1 = 10;
  PDM_g_num_t connec_tetra_1[40] = {1, 2, 4, 7,
                                    2, 5, 4, 11,
                                    4, 7, 11, 10,
                                    2, 7, 8, 11,
                                    2, 4, 7, 11,
                                    2, 6, 5, 11,
                                    2, 9, 3, 6,
                                    11, 12, 9, 6,
                                    9, 11, 2, 8,
                                    9, 2, 11, 6};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 3, n_cell);

  int tetra_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_TETRA4);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        tetra_section_1,
                                        n_tetra_section_1,
                                        connec_tetra_1,
                                        PDM_OWNERSHIP_USER);
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  // printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 40  );
  CHECK( n_sum_vtx_face_tot == 120 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);
  std::vector<int>         dparent_elmt_pos(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                                  dcell_face_vtx_idx.data(),
                                  dcell_face_vtx.data(),
                                  delmt_face_cell.data(),
                                  NULL, NULL,
                                  dparent_elmt_pos.data());
  // dcell_face_vtx_idx[0] = 0;
  // PDM_sections_decompose_faces(dmn,
  //                                 dcell_face_vtx_idx,
  //                                 dcell_face_vtx,
  //                                 delmt_face_cell,
  //                                 NULL, NULL);

  if( 0 == 1)
  {
    PDM_log_trace_array_long(delmt_face_cell.data(), n_face_elt_tot, "delmt_face_cell:: ");
    PDM_log_trace_array_int(dcell_face_vtx_idx.data(), n_face_elt_tot+1, "dcell_face_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_face_vtx.data(), n_sum_vtx_face_tot, "dcell_face_vtx:: ");
  }

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {1, 4, 2, 1, 2, 7, 1, 7, 4, 2, 4, 7, 2, 4, 5, 2, 5, 11, 2, 11, 4, 5, 4, 11, 4, 11, 7, 4, 7, 10, 4, 10, 11, 7, 11, 10, 2, 8, 7, 2, 7, 11, 2, 11, 8, 7, 8, 11, 2, 7, 4, 2, 4, 11, 2, 11, 7, 4, 7, 11, 2, 5, 6, 2, 6, 11, 2, 11, 5, 6, 5, 11, 2, 3, 9, 2, 9, 6, 2, 6, 3, 9, 3, 6, 11, 9, 12, 11, 12, 6, 11, 6, 9, 12, 9, 6, 9, 2, 11, 9, 11, 8, 9, 8, 2, 11, 2, 8, 9, 11, 2, 9, 2, 6, 9, 6, 11, 2, 11, 6};
  std::vector<int>         dparent_elmt_pos_expected   = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  PDM_DMesh_nodal_elmts_free(dmn);
}


// double coord_pyra_x[9] = {0., 1., 0., 1., 0., 1., 0., 1., 0.5};
// double coord_pyra_y[9] = {0., 0., 1., 1., 0., 0., 1., 1., 0.5};
// double coord_pyra_z[9] = {0., 0., 0., 0., 1., 1., 1., 1., 0.5};

MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes pyra",1) {
  // const PDM_g_num_t n_vtx            = 9;
  const PDM_g_num_t n_cell           = 6;
  const int         n_pyra_section_1 = 6;
  PDM_g_num_t connec_pyra_1[30] = {1, 2, 4, 3, 9,
                                   5, 7, 8, 6, 9,
                                   1, 3, 7, 5, 9,
                                   2, 6, 8, 4, 9,
                                   1, 5, 6, 2, 9,
                                   3, 4, 8, 7, 9};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 3, n_cell);

  int pyra_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_PYRAMID5);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        pyra_section_1,
                                        n_pyra_section_1,
                                        connec_pyra_1,
                                        PDM_OWNERSHIP_USER);
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  // printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 30 );
  CHECK( n_sum_vtx_face_tot == 96 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);
  std::vector<int>         dparent_elmt_pos(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                                  dcell_face_vtx_idx.data(),
                                  dcell_face_vtx.data(),
                                  delmt_face_cell.data(),
                                  NULL, NULL,
                                  dparent_elmt_pos.data());

  if( 0 == 1)
  {
    PDM_log_trace_array_long(delmt_face_cell.data()  , n_face_elt_tot    , "delmt_face_cell:: ");
    PDM_log_trace_array_int(dcell_face_vtx_idx.data(), n_face_elt_tot+1  , "dcell_face_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_face_vtx.data()   , n_sum_vtx_face_tot, "dcell_face_vtx:: ");
  }

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 4, 7, 10, 13, 16, 20, 23, 26, 29, 32, 36, 39, 42, 45, 48, 52, 55, 58, 61, 64, 68, 71, 74, 77, 80, 84, 87, 90, 93, 96};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {3 ,4, 2, 1, 9, 1, 2, 9, 2, 4, 9, 4, 3, 9, 3, 1, 6, 8, 7, 5, 9, 5, 7, 9, 7, 8, 9, 8, 6, 9, 6, 5, 5, 7, 3, 1, 9, 1, 3, 9, 3, 7, 9, 7, 5, 9, 5, 1, 4, 8, 6, 2, 9, 2, 6, 9, 6, 8, 9, 8, 4, 9, 4, 2, 2, 6, 5, 1, 9, 1, 5, 9, 5, 6, 9, 6, 2, 9, 2, 1, 7, 8, 4, 3, 9, 3, 4, 9, 4, 8, 9, 8, 7, 9, 7, 3};
  std::vector<int>         dparent_elmt_pos_expected   = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  PDM_DMesh_nodal_elmts_free(dmn);
}


// double coord_prism_x[12] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2}
// double coord_prism_y[12] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1}
// double coord_prism_z[12] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1}

MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes prism",1) {
  // const PDM_g_num_t n_vtx             = 12;
  const PDM_g_num_t n_cell            = 4;
  const int         n_prism_section_1 = 4;
  PDM_g_num_t connec_prism_1[24] = {1, 5, 4, 7, 11, 10,
                                    1, 2, 5, 7, 8, 11,
                                    2, 6, 5, 8, 12, 11,
                                    2, 3, 6, 8, 9, 12};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 3, n_cell);

  int prism_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_PRISM6);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        prism_section_1,
                                        n_prism_section_1,
                                        connec_prism_1,
                                        PDM_OWNERSHIP_USER);
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  // printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 20 );
  CHECK( n_sum_vtx_face_tot == 72 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);
  std::vector<int>         dparent_elmt_pos(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmn,
                               dcell_face_vtx_idx.data(),
                               dcell_face_vtx.data(),
                               delmt_face_cell.data(),
                               NULL, NULL,
                               dparent_elmt_pos.data());
  if( 0 == 1)
  {
    PDM_log_trace_array_long(delmt_face_cell.data()  , n_face_elt_tot    , "delmt_face_cell:: ");
    PDM_log_trace_array_int(dcell_face_vtx_idx.data(), n_face_elt_tot+1  , "dcell_face_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_face_vtx.data()   , n_sum_vtx_face_tot, "dcell_face_vtx:: ");
  }

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 3, 6, 10, 14, 18, 21, 24, 28, 32, 36, 39, 42, 46, 50, 54, 57, 60, 64, 68, 72};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {4, 5, 1, 11, 10, 7, 10, 11, 5, 4, 11, 7, 1, 5, 7, 10, 4, 1, 5, 2, 1, 8, 11, 7, 11, 8, 2, 5, 8, 7, 1, 2, 7, 11, 5, 1, 5, 6, 2, 12, 11, 8, 11, 12, 6, 5, 12, 8, 2, 6, 8, 11, 5, 2, 6, 3, 2, 9, 12, 8, 12, 9, 3, 6, 9, 8, 2, 3, 8, 12, 6, 2};
  std::vector<int>         dparent_elmt_pos_expected   = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  // printf(" dface_vtx_idx[dn_face]::%i\n",  dface_vtx_idx[dn_face]);
  // PDM_log_trace_array_long(dface_cell, 2*dn_face, "dface_cell:: ");
  // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1, "dface_vtx_idx:: ");
  // PDM_log_trace_array_long(dface_vtx, dface_vtx_idx[dn_face], "dface_vtx:: ");

  PDM_DMesh_nodal_elmts_free(dmn);
}



// double coord_quad_x[9] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
// double coord_quad_y[9] = {0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0};
// double coord_quad_z[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes quad",1) {
  // const PDM_g_num_t n_vtx            = 9;
  const PDM_g_num_t n_face           = 4;
  const int         n_quad_section_1 = 4;
  PDM_g_num_t connec_quad_1[16] = {1, 2, 5, 4,
                                   2, 3, 6, 5,
                                   4, 5, 8, 7,
                                   5, 6, 9, 8,};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 2, n_face);

  int quad_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        quad_section_1,
                                        n_quad_section_1,
                                        connec_quad_1,
                                        PDM_OWNERSHIP_USER);
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_edge_elt_tot     = -1;
  int n_sum_vtx_edge_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmn, &n_edge_elt_tot, &n_sum_vtx_edge_tot);

  // printf("n_edge_elt_tot     = %i\n", n_edge_elt_tot);
  // printf("n_sum_vtx_edge_tot = %i\n", n_sum_vtx_edge_tot);

  CHECK( n_edge_elt_tot     == 16 );
  CHECK( n_sum_vtx_edge_tot == 32 );

  // PDM_g_num_t* delmt_edge_cell    = (PDM_g_num_t*) malloc(  n_edge_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_edge_vtx_idx = (int        *) malloc( (n_edge_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_edge_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_edge_vtx_idx(n_edge_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_edge_vtx(n_sum_vtx_edge_tot);
  std::vector<PDM_g_num_t> delmt_edge_cell(n_edge_elt_tot, -1);
  std::vector<int>         dparent_elmt_pos(n_edge_elt_tot);

  dcell_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmn,
                               dcell_edge_vtx_idx.data(),
                               dcell_edge_vtx.data(),
                               delmt_edge_cell.data(),
                               NULL, NULL,
                               dparent_elmt_pos.data());

  if( 1 == 1)
  {
    PDM_log_trace_array_long(delmt_edge_cell.data()  , n_edge_elt_tot    , "delmt_edge_cell:: ");
    PDM_log_trace_array_int(dcell_edge_vtx_idx.data(), n_edge_elt_tot+1  , "dcell_edge_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_edge_vtx.data()   , n_sum_vtx_edge_tot, "dcell_edge_vtx:: ");
  }

  std::vector<PDM_g_num_t> delmt_edge_cell_expected    = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
  std::vector<int>         dcell_edge_vtx_idx_expected = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32};
  std::vector<PDM_g_num_t> dcell_edge_vtx_expected     = {1, 2, 2, 5, 5, 4, 4, 1, 2, 3, 3, 6, 6, 5, 5, 2, 4, 5, 5, 8, 8, 7, 7, 4, 5, 6, 6, 9, 9, 8, 8, 5};
  std::vector<int>         dparent_elmt_pos_expected   = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};

  CHECK( delmt_edge_cell    == delmt_edge_cell_expected);
  CHECK( dcell_edge_vtx_idx == dcell_edge_vtx_idx_expected);
  CHECK( dcell_edge_vtx     == dcell_edge_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  PDM_DMesh_nodal_elmts_free(dmn);
}




// // double coord_tri_x[9] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
// // double coord_tri_y[9] = {0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0};
// // double coord_tri_z[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes tri",1) {
  // const PDM_g_num_t n_vtx            = 9;
  const PDM_g_num_t n_face           = 8;
  const int         n_tri_section_1  = 8;
  PDM_g_num_t connec_tri_1[24] = {1, 2, 5,
                                  1, 5, 4,
                                  2, 3, 6,
                                  2, 6, 5,
                                  4, 5, 8,
                                  4, 8, 7,
                                  5, 6, 9,
                                  5, 9, 8};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_elmts_t* dmn = PDM_DMesh_nodal_elmts_create(pdm_comm, 2, n_face);

  int tri_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn, PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_elmts_section_std_set(dmn,
                                        tri_section_1,
                                        n_tri_section_1,
                                        connec_tri_1,
                                        PDM_OWNERSHIP_USER);
  PDM_dmesh_nodal_elmts_generate_distribution(dmn);

  int n_edge_elt_tot     = -1;
  int n_sum_vtx_edge_tot = -1;
  PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmn, &n_edge_elt_tot, &n_sum_vtx_edge_tot);

  // printf("n_edge_elt_tot     = %i\n", n_edge_elt_tot);
  // printf("n_sum_vtx_edge_tot = %i\n", n_sum_vtx_edge_tot);

  CHECK( n_edge_elt_tot     == 24 );
  CHECK( n_sum_vtx_edge_tot == 48 );

  // PDM_g_num_t* delmt_edge_cell    = (PDM_g_num_t*) malloc(  n_edge_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_edge_vtx_idx = (int        *) malloc( (n_edge_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_edge_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_edge_vtx_idx(n_edge_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_edge_vtx(n_sum_vtx_edge_tot);
  std::vector<PDM_g_num_t> delmt_edge_cell(n_edge_elt_tot, -1);
  std::vector<int>         dparent_elmt_pos(n_edge_elt_tot, -1);

  dcell_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmn,
                               dcell_edge_vtx_idx.data(),
                               dcell_edge_vtx.data(),
                               delmt_edge_cell.data(),
                               NULL, NULL,
                               dparent_elmt_pos.data());
  // // dcell_edge_vtx_idx[0] = 0;
  // // PDM_dmesh_nodal_decompose_edges(dmn,
  // //                                 dcell_edge_vtx_idx,
  // //                                 dcell_edge_vtx,
  // //                                 delmt_edge_cell,
  // //                                 NULL, NULL);

  if( 0 == 1)
  {
    PDM_log_trace_array_long(delmt_edge_cell.data()  , n_edge_elt_tot    , "delmt_edge_cell:: ");
    PDM_log_trace_array_int(dcell_edge_vtx_idx.data(), n_edge_elt_tot+1  , "dcell_edge_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_edge_vtx.data()   , n_sum_vtx_edge_tot, "dcell_edge_vtx:: ");
  }

  std::vector<PDM_g_num_t> delmt_edge_cell_expected    = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8};
  std::vector<int>         dcell_edge_vtx_idx_expected = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48 };
  std::vector<PDM_g_num_t> dcell_edge_vtx_expected     = {1, 2, 2, 5, 5, 1, 1, 5, 5, 4, 4, 1, 2, 3, 3, 6, 6, 2, 2, 6, 6, 5, 5, 2, 4, 5, 5, 8, 8, 4, 4, 8, 8, 7, 7, 4, 5, 6, 6, 9, 9, 5, 5, 9, 9, 8, 8, 5};
  std::vector<int>         dparent_elmt_pos_expected   = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};

  CHECK( delmt_edge_cell    == delmt_edge_cell_expected);
  CHECK( dcell_edge_vtx_idx == dcell_edge_vtx_idx_expected);
  CHECK( dcell_edge_vtx     == dcell_edge_vtx_expected);
  CHECK( dparent_elmt_pos   == dparent_elmt_pos_expected);

  PDM_DMesh_nodal_elmts_free(dmn);
}




