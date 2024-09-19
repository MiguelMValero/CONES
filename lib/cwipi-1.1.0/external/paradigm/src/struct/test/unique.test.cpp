#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_unique.h"
#include "pdm_logging.h"

// PDM_log_trace_array_int(block_pts_per_elt_n , block_n_elt, "block_pts_per_elt_n  AFTER :: ");

TEST_CASE("[pdm_unique] - PDM_inplace_unique") {

  std::vector<int> array = { 2, 4, 5, 4, 2, 1 };
  int size = PDM_inplace_unique(array.data(), 0, static_cast<int>(array.size()-1));

  array.resize(size);
  std::vector<int> array_expected = { 1, 2, 4, 5};

  CHECK( size == 4);
  CHECK( array == array_expected);

}

TEST_CASE("[pdm_unique] - PDM_inplace_unique_long") {

  std::vector<PDM_g_num_t> array = { 2, 4, 5, 4, 2, 1 };
  std::vector<int> order(array.size());
  std::iota(begin(order), end(order), 1 );
  int size = PDM_inplace_unique_long(array.data(), order.data(), 0, static_cast<int>(array.size()-1));

  array.resize(size);
  std::vector<PDM_g_num_t> array_expected = { 1, 2, 4, 5};

  CHECK( size == 4);
  CHECK( array == array_expected);

  // PDM_log_trace_array_int(order.data(), 6, "order:: ");

  std::vector<int> order_expected = {6, 1, 5, 2, 4, 3};
  CHECK( order == order_expected);

  // So : array[order[i]] = idx_old // Ici array[order[0]] = 6 ( Car 1 est bien à la fin de array)

}

TEST_CASE("[pdm_unique] - PDM_inplace_unique_long2") {

  std::vector<PDM_g_num_t> array = { 2, 4, 5, 4, 2, 1 };
  std::vector<int> unique_order(array.size());
  std::iota(begin(unique_order), end(unique_order), 1 );
  int size = PDM_inplace_unique_long2(array.data(), unique_order.data(), 0, static_cast<int>(array.size()-1));

  array.resize(size);
  std::vector<PDM_g_num_t> array_expected = { 1, 2, 4, 5};

  CHECK( size == 4);
  CHECK( array == array_expected);

  // PDM_log_trace_array_int(unique_order.data(), 6, "unique_order:: ");

  std::vector<int> unique_order_expected = {1, 2, 3, 2, 1, 0};
  CHECK( unique_order == unique_order_expected);

  // So : array[order_unique[i]] = no_unique // Ici array[order_unique[0]] = 1 car 2 est en position 1


}
