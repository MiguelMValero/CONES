#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"

#include "pdm_array.h"

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_zeros", 1) {
  int *array = PDM_array_zeros_int(5);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == 0);
  free(array);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_const", 1) {
  int *array = PDM_array_const_int(7, 42);
  for (int i = 0; i < 7; i++)
    CHECK(array[i] == 42);
  free(array);

  PDM_g_num_t *array_gnum = PDM_array_const_gnum(10, 33550336);
  for (int i = 0; i < 10; i++)
    CHECK(array_gnum[i] == 33550336);
  free(array_gnum);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_reset", 1) {
  int *array = (int *) malloc(5*sizeof(int));
  for (int i = 0; i < 5; i++)
    array[i] = i;

  PDM_array_reset_int(array, 5, 0);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == 0);
  PDM_array_reset_int(array, 5, -999);
  for (int i = 0; i < 5; i++)
    CHECK(array[i] == -999);

  free(array);

  PDM_g_num_t *array_gnum = (PDM_g_num_t *) malloc(5*sizeof(PDM_g_num_t));
  for (int i = 0; i < 5; i++)
    array_gnum[i] = i;

  PDM_array_reset_gnum(array_gnum, 5, 0);
  for (int i = 0; i < 5; i++)
    CHECK(array_gnum[i] == 0);
  free(array_gnum);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_new_idx_from_sizes", 1) {
  int size_array[] = {5, 5, 2, 5};
  int *idx_array = PDM_array_new_idx_from_sizes_int(size_array, 4);
  int expected_idx_array[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array, expected_idx_array, 4+1);
  free(idx_array);

  PDM_g_num_t *idx_array_gnum = PDM_array_new_idx_from_sizes_gnum(size_array, 4);
  PDM_g_num_t expected_idx_array_gnum[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array_gnum, expected_idx_array_gnum, 4+1);
  free(idx_array_gnum);

  int empty_array[] = {};
  idx_array = PDM_array_new_idx_from_sizes_int(empty_array, 0);
  CHECK(idx_array[0] == 0);
  free(idx_array);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_idx_from_sizes", 1) {
  int size_array[] = {5, 5, 2, 5};
  int idx_array[] = {-1,-1,-1,-1,-1};
  PDM_array_idx_from_sizes_int(size_array, 4, idx_array);
  int expected_idx_array[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array, expected_idx_array, 4+1);

  PDM_g_num_t idx_array_gnum[] = {-1,-1,-1,-1,-1};
  PDM_array_idx_from_sizes_gnum(size_array, 4, idx_array_gnum);
  int expected_idx_array_gnum[] = {0, 5, 10, 12, 17};
  CHECK_EQ_C_ARRAY(idx_array_gnum, expected_idx_array_gnum, 4+1);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_are_equal", 1) {
  int arrayA[] = {0,3,6,4,2};
  int arrayB[] = {0,3,6,4,2};
  int arrayC[] = {0,3,6,4,3};
  CHECK(PDM_array_are_equal_int(arrayA, arrayB, 5) == 1);
  CHECK(PDM_array_are_equal_int(arrayA, arrayC, 5) == 0);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_accumulate", 1) {
  int array[] = {0,3,6,4,2};
  int expected_array[] = {0,3,9,13,15};
  PDM_array_accumulate_int(array, 5);
  CHECK_EQ_C_ARRAY(array, expected_array, 5);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_count_per_col", 1) {
  int color_array[] = {3,2,1,2,4,3,2,1,2,3,2,1,2,3,4};
  int n_per_col[5];
  int expected_n_per_col[] = {0,3,6,4,2};
  PDM_array_count_per_col_int(5, 15, color_array, n_per_col);
  CHECK_EQ_C_ARRAY(n_per_col, expected_n_per_col, 5);
}

MPI_TEST_CASE("[pdm_array] - 1p - PDM_array_repart_per_col", 1) {
  int color_array[] = {3,2,1,2,4,3,2,1,2,3,2,1,2,3,4};
  int ordered_idx[5+1];
  int ordered[15];
  int expected_ordered_idx[] = {0, 0, 3, 9, 13, 15};
  int expected_ordered[] = {2,7,11,1,3,6,8,10,12,0,5,9,13,4,14};
  PDM_array_repart_per_col_int(5, 15, color_array, ordered_idx, ordered);
  CHECK_EQ_C_ARRAY(ordered_idx, expected_ordered_idx, 5+1);
  CHECK_EQ_C_ARRAY(ordered, expected_ordered, 15);
}

