#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_logging.h"

/*
 *  Use case
 *     ./paradigm/test/pdm_t_partitioning_dcube -n 3 -n_part 1 -pt-scotch (Générateur dcube_gen sur 8 cellules)
 */

MPI_TEST_CASE("[pdm_dconnectivity_transform] - 1p - dcell_face + dface_vtx = dcell_vtx ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_g_num_t cell_distrib[2]   = {0, 8};
  PDM_g_num_t face_distrib[2]   = {0, 36};
  PDM_g_num_t vtx_distrib[2]    = {0, 27};
  int         dcell_face_idx[9] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
  PDM_g_num_t dcell_face[48]    = {1, 5, 13, 17, 25, 29, 2, 6, -17, 21, 27, 31, 3, 7, 14, 18, -29, 33, 4, 8, -18, 22, -31, 35, -5, 9, 15, 19, 26, 30, -6, 10, -19, 23, 28, 32, -7, 11, 16, 20, -30, 34, -8, 12, -20, 24, -32, 36};
  int         dface_vtx_idx[37] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144};
  PDM_g_num_t dface_vtx[144]    = {1, 4, 5, 2, 2, 5, 6, 3, 4, 7, 8, 5, 5, 8, 9, 6, 10, 13, 14, 11, 11, 14, 15, 12, 13, 16, 17, 14, 14, 17, 18, 15, 19, 22, 23, 20, 20, 23, 24, 21, 22, 25, 26, 23, 23, 26, 27, 24, 1, 4, 13, 10, 4, 7, 16, 13, 10, 13, 22, 19, 13, 16, 25, 22, 2, 5, 14, 11, 5, 8, 17, 14, 11, 14, 23, 20, 14, 17, 26, 23, 3, 6, 15, 12, 6, 9, 18, 15, 12, 15, 24, 21, 15, 18, 27, 24, 1, 2, 11, 10, 10, 11, 20, 19, 2, 3, 12, 11, 11, 12, 21, 20, 4, 5, 14, 13, 13, 14, 23, 22, 5, 6, 15, 14, 14, 15, 24, 23, 7, 8, 17, 16, 16, 17, 26, 25, 8, 9, 18, 17, 17, 18, 27, 26};

  int*         dcell_vtx_idx;
  PDM_g_num_t* dcell_vtx;
  PDM_deduce_combine_connectivity(pdm_comm,
                                  cell_distrib,
                                  face_distrib,
                                  dcell_face_idx,
                                  dcell_face,
                                  dface_vtx_idx,
                                  dface_vtx,
                                  0,
                                  &dcell_vtx_idx,
                                  &dcell_vtx);

  int dn_cell = cell_distrib[test_rank+1] - cell_distrib[test_rank];

  int         dcell_vtx_idx_expected[9] = {0, 8, 16, 24, 32, 40, 48, 56, 64};
  PDM_g_num_t dcell_vtx_expected[64]    = {1, 2, 4, 5, 10, 11, 13, 14, 2, 3, 5, 6, 11, 12, 14, 15, 4, 5, 7, 8, 13, 14, 16, 17, 5, 6, 8, 9, 14, 15, 17, 18, 10, 11, 13, 14, 19, 20, 22, 23, 11, 12, 14, 15, 20, 21, 23, 24, 13, 14, 16, 17, 22, 23, 25, 26, 14, 15, 17, 18, 23, 24, 26, 27};

  MPI_CHECK_EQ_C_ARRAY(0, dcell_vtx_idx, dcell_vtx_idx_expected, dn_cell+1                      );
  MPI_CHECK_EQ_C_ARRAY(0, dcell_vtx    , dcell_vtx_expected    , dcell_vtx_idx_expected[dn_cell]);

  int*         dvtx_cell_idx;
  PDM_g_num_t* dvtx_cell;
  PDM_dconnectivity_transpose(pdm_comm,
                              cell_distrib,
                              vtx_distrib,
                              dcell_vtx_idx,
                              dcell_vtx,
                              1,
                              &dvtx_cell_idx,
                              &dvtx_cell);

  // int*         dcell_vtx_idx_2;
  // PDM_g_num_t* dcell_vtx_2;
  // PDM_dconnectivity_transpose(pdm_comm,
  //                              vtx_distrib,
  //                              cell_distrib,
  //                              dvtx_cell_idx,
  //                              dvtx_cell,
  //                              0,
  //                              &dcell_vtx_idx_2,
  //                              &dcell_vtx_2);

  free(dcell_vtx_idx);
  free(dcell_vtx);
  free(dvtx_cell_idx);
  free(dvtx_cell);

}

MPI_TEST_CASE("[pdm_dconnectivity_transform] - 2p - dcell_face + dface_vtx = dcell_vtx ",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  // PDM_g_num_t vtx_distrib[3]    = {1, 14, 28};

  std::vector<PDM_g_num_t> cell_distrib_p0   = {0, 4, 8};
  std::vector<PDM_g_num_t> face_distrib_p0   = {0, 18, 36};
  std::vector<int        > dcell_face_idx_p0 = {0, 6, 12, 18, 24};
  std::vector<PDM_g_num_t> dcell_face_p0     = {1, 5, 13, 17, 25, 29, 2, 6, -17, 21, 27, 31, 3, 7, 14, 18, -29, 33, 4, 8, -18, 22, -31, 35};
  std::vector<int        > dface_vtx_idx_p0  = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72};
  std::vector<PDM_g_num_t> dface_vtx_p0      = {1, 4, 5, 2, 2, 5, 6, 3, 4, 7, 8, 5, 5, 8, 9, 6, 10, 13, 14, 11, 11, 14, 15, 12, 13, 16, 17, 14, 14, 17, 18, 15, 19, 22, 23, 20, 20, 23, 24, 21, 22, 25, 26, 23, 23, 26, 27, 24, 1, 4, 13, 10, 4, 7, 16, 13, 10, 13, 22, 19, 13, 16, 25, 22, 2, 5, 14, 11, 5, 8, 17, 14};

  std::vector<PDM_g_num_t> cell_distrib_p1   = {0, 4, 8};
  std::vector<PDM_g_num_t> face_distrib_p1   = {0, 18, 36};
  std::vector<int        > dcell_face_idx_p1 = {0, 6, 12, 18, 24};
  std::vector<PDM_g_num_t> dcell_face_p1     = {-5, 9, 15, 19, 26, 30, -6, 10, -19, 23, 28, 32, -7, 11, 16, 20, -30, 34, -8, 12, -20, 24, -32, 36};
  std::vector<int        > dface_vtx_idx_p1  = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72};
  std::vector<PDM_g_num_t> dface_vtx_p1      = {11, 14, 23, 20, 14, 17, 26, 23, 3, 6, 15, 12, 6, 9, 18, 15, 12, 15, 24, 21, 15, 18, 27, 24, 1, 2, 11, 10, 10, 11, 20, 19, 2, 3, 12, 11, 11, 12, 21, 20, 4, 5, 14, 13, 13, 14, 23, 22, 5, 6, 15, 14, 14, 15, 24, 23, 7, 8, 17, 16, 16, 17, 26, 25, 8, 9, 18, 17, 17, 18, 27, 26};

  PDM_g_num_t *cell_distrib   = NULL;
  PDM_g_num_t *face_distrib   = NULL;
  int         *dcell_face_idx = NULL;
  PDM_g_num_t *dcell_face     = NULL;
  int         *dface_vtx_idx  = NULL;
  PDM_g_num_t *dface_vtx      = NULL;

  if(test_rank == 0) {
    cell_distrib   = cell_distrib_p0  .data();
    face_distrib   = face_distrib_p0  .data();
    dcell_face_idx = dcell_face_idx_p0.data();
    dcell_face     = dcell_face_p0    .data();
    dface_vtx_idx  = dface_vtx_idx_p0 .data();
    dface_vtx      = dface_vtx_p0     .data();
  } else if(test_rank == 1) {
    cell_distrib   = cell_distrib_p1  .data();
    face_distrib   = face_distrib_p1  .data();
    dcell_face_idx = dcell_face_idx_p1.data();
    dcell_face     = dcell_face_p1    .data();
    dface_vtx_idx  = dface_vtx_idx_p1 .data();
    dface_vtx      = dface_vtx_p1     .data();
  }

  int*         dcell_vtx_idx;
  PDM_g_num_t* dcell_vtx;
  PDM_deduce_combine_connectivity(pdm_comm,
                                  cell_distrib,
                                  face_distrib,
                                  dcell_face_idx,
                                  dcell_face,
                                  dface_vtx_idx,
                                  dface_vtx,
                                  0,
                                  &dcell_vtx_idx,
                                  &dcell_vtx);

  int dn_cell = cell_distrib[test_rank+1] - cell_distrib[test_rank];

  PDM_g_num_t dcell_vtx_idx_expected_p0[5] = {0, 8, 16, 24, 32};
  PDM_g_num_t dcell_vtx_expected_p0[32]    = {1, 2, 4, 5, 10, 11, 13, 14, 2, 3, 5, 6, 11, 12, 14, 15, 4, 5, 7, 8, 13, 14, 16, 17, 5, 6, 8, 9, 14, 15, 17, 18};

  PDM_g_num_t dcell_vtx_idx_expected_p1[5] = {0, 8, 16, 24, 32};
  PDM_g_num_t dcell_vtx_expected_p1[32]    = {10, 11, 13, 14, 19, 20, 22, 23, 11, 12, 14, 15, 20, 21, 23, 24, 13, 14, 16, 17, 22, 23, 25, 26, 14, 15, 17, 18, 23, 24, 26, 27};

  MPI_CHECK_EQ_C_ARRAY(0, dcell_vtx_idx, dcell_vtx_idx_expected_p0, dn_cell+1                         );
  MPI_CHECK_EQ_C_ARRAY(0, dcell_vtx    , dcell_vtx_expected_p0    , dcell_vtx_idx_expected_p0[dn_cell]);

  MPI_CHECK_EQ_C_ARRAY(1, dcell_vtx_idx, dcell_vtx_idx_expected_p1, dn_cell+1                         );
  MPI_CHECK_EQ_C_ARRAY(1, dcell_vtx    , dcell_vtx_expected_p1    , dcell_vtx_idx_expected_p1[dn_cell]);

  // int*         dvtx_cell_idx;
  // PDM_g_num_t* dvtx_cell;
  // PDM_dconnectivity_transpose(pdm_comm,
  //                              cell_distrib,
  //                              vtx_distrib,
  //                              dcell_vtx_idx,
  //                              dcell_vtx,
  //                              0,
  //                              &dvtx_cell_idx,
  //                              &dvtx_cell);
  // free(dvtx_cell_idx);
  // free(dvtx_cell);

  free(dcell_vtx_idx);
  free(dcell_vtx);

}


MPI_TEST_CASE("[pdm_dconnectivity_transform] - 1p - PDM_dorder_reverse",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<PDM_g_num_t> entity_distrib = {0, 5};
  int dn_entity = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  std::vector<PDM_g_num_t> new_to_old = {5, 3, 2, 1, 4};

  PDM_g_num_t* old_to_new;
  PDM_dorder_reverse(pdm_comm,
                     entity_distrib.data(),
                     new_to_old.data(),
                     &old_to_new);

  std::vector<PDM_g_num_t> old_to_new_expected = {4, 3, 2, 5, 1};
  MPI_CHECK_EQ_C_ARRAY(0, old_to_new, old_to_new_expected.data(), dn_entity);

  // PDM_log_trace_array_long(old_to_new, dn_entity, "old_to_new : ");

}


MPI_TEST_CASE("[pdm_dconnectivity_transform] - 2p - PDM_dorder_reverse",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<PDM_g_num_t> entity_distrib = {0, 3, 5};
  std::vector<PDM_g_num_t> new_to_old;
  if(i_rank == 0) {
    new_to_old = {5, 3, 2};
  } else {
    new_to_old = {1, 4};
  }

  int dn_entity = entity_distrib[i_rank+1] - entity_distrib[i_rank];

  PDM_g_num_t* old_to_new;
  PDM_dorder_reverse(pdm_comm,
                     entity_distrib.data(),
                     new_to_old.data(),
                     &old_to_new);

  std::vector<PDM_g_num_t> old_to_new_expected_p0 = {4, 3, 2};
  std::vector<PDM_g_num_t> old_to_new_expected_p1 = {5, 1};

  MPI_CHECK_EQ_C_ARRAY(0, old_to_new, old_to_new_expected_p0.data(), dn_entity);
  MPI_CHECK_EQ_C_ARRAY(1, old_to_new, old_to_new_expected_p1.data(), dn_entity);

  // PDM_log_trace_array_long(old_to_new, dn_entity, "old_to_new : ");

}
