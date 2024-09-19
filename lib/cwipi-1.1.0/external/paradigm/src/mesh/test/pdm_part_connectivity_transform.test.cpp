#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_part_connectivity_transform.h"

/*
 *  Use case
 *     ./paradigm/test/pdm_t_partitioning_dcube -n 3 -n_part 1 -pt-scotch (Générateur dcube_gen sur 8 cellules)
 */

MPI_TEST_CASE("[pdm_part_connectivity_transform] - 1p - pdm_connectivity_transform ",1) {

  const int n_cell = 8;
  const int n_face = 36;

  int face_cell_idx[37] = {0 , 1, 2, 3, 4, 6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26,
                           28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 44, 45, 46, 47, 48};
  int face_cell    [48] = {1, 6, 7, 4, 1, 5, 6, 2, 7, 3, 4, 8, 5, 2, 3, 8, 1, 7, 5, 3, 1, 6, 7,
                           4, 5, 2, 3, 8, 6, 4, 2, 8, 1, 5, 6, 2, 1, 7, 5, 3, 6, 4, 2, 8, 7, 3, 4, 8};

  int cell_face_idx[n_cell+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
  int cell_face    [48]       = { 1,  5,  13, 17, 25, 29,
                                 -6, 10, -19, 23, 28, 32,
                                 -7, 11,  16, 20,-30, 34,
                                  4,  8, -18, 22,-31, 35,
                                 -5,  9,  15, 19, 26, 30,
                                  2,  6, -17, 21, 27, 31,
                                  3,  7, 14, 18, -29, 33,
                                 -8, 12,-20, 24, -32, 36};
  int* cell_cell_idx;
  int* cell_cell;
  PDM_combine_connectivity(n_cell,
                           cell_face_idx,
                           cell_face,
                           face_cell_idx,
                           face_cell,
                           &cell_cell_idx,
                           &cell_cell);

  // PDM_log_trace_array_int(cell_cell_idx, n_cell+1         , "cell_cell_idx::");
  // PDM_log_trace_array_int(cell_cell, cell_cell_idx[n_cell], "cell_cell::");

  int cell_cell_idx_expected[n_cell+1] = {0, 4, 8, 12, 16, 20, 24, 28, 32};
  int cell_cell_expected[32]           = {1, 5, 6, 7, 2, 5, 6, 8, 3, 5, 7, 8, 4, 6, 7, 8, 1, 2, 3, 5, 1, 2, 4, 6, 1, 3, 4, 7, 2, 3, 4, 8};

  CHECK_EQ_C_ARRAY(cell_cell_idx, cell_cell_idx_expected, n_cell+1);
  CHECK_EQ_C_ARRAY(cell_cell    , cell_cell_expected    , cell_cell_idx[n_cell]);

  free(cell_cell_idx);
  free(cell_cell);


  int* cell_face_from_transpose_idx;
  int* cell_face_from_transpose;
  PDM_connectivity_transpose(n_face,
                             n_cell,
                             face_cell_idx,
                             face_cell,
                             &cell_face_from_transpose_idx,
                             &cell_face_from_transpose);

  int cell_face_from_transpose_idx_expected[n_cell+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
  int cell_face_from_transpose_expected    [48]       = {1, 5, 13, 17, 25, 29, 6, 10, 19, 23, 28, 32, 7, 11, 16, 20, 30, 34, 4, 8, 18, 22, 31, 35, 5, 9, 15, 19, 26, 30, 2, 6, 17, 21, 27, 31, 3, 7, 14, 18, 29, 33, 8, 12, 20, 24, 32, 36};

  CHECK_EQ_C_ARRAY(cell_face_from_transpose_idx, cell_face_from_transpose_idx_expected, n_cell+1);
  CHECK_EQ_C_ARRAY(cell_face_from_transpose    , cell_face_from_transpose_expected    , cell_face_from_transpose_idx_expected[n_cell]);

  free(cell_face_from_transpose_idx);
  free(cell_face_from_transpose);

  int* face_cell_from_transpose_idx;
  int* face_cell_from_transpose;
  PDM_connectivity_transpose(n_cell,
                             n_face,
                             cell_face_idx,
                             cell_face,
                             &face_cell_from_transpose_idx,
                             &face_cell_from_transpose);

  int face_cell_from_transpose_idx_expected[n_face+1] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 44, 45, 46, 47, 48};
  int face_cell_from_transpose_expected    [48]       = {1, 6, 7, 4, 1, -5, -2, 6, -3, 7, 4, -8, 5, 2, 3, 8, 1, 7, 5, 3, 1, -6, -4, 7, -2, 5, 3, -8, 6, 4, 2, 8, 1, 5, 6, 2, 1, -7, -3, 5, -4, 6, 2, -8, 7, 3, 4, 8};

  CHECK_EQ_C_ARRAY(face_cell_from_transpose_idx, face_cell_from_transpose_idx_expected, n_face+1);
  CHECK_EQ_C_ARRAY(face_cell_from_transpose    , face_cell_from_transpose_expected    , face_cell_from_transpose_idx_expected[n_face]);

  free(face_cell_from_transpose_idx);
  free(face_cell_from_transpose);

}
