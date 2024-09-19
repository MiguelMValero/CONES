#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_order.h"
#include "pdm_logging.h"


TEST_CASE("[pdm_order] - PDM_order_inplace_unique_and_order_long - stride = 1") {

  std::vector<int> array = { 2, 4, 5, 4, 2, 1 };

  int n_entity = array.size();
  std::vector<int> order(array.size());

  int n_unique = PDM_order_inplace_unique_and_order_long(n_entity,
                                                         1,
                                                         array.data(),
                                                         order.data());

  // PDM_log_trace_array_int(array.data() , n_unique, "array :: ");
  // PDM_log_trace_array_int(order.data() , n_unique, "order :: ");

  std::vector<int> array_expected = { 1, 2, 4, 5};
  std::vector<int> order_expected = { 5, 0, 3, 2};

  CHECK( n_unique == 4);

  array.resize(n_unique);
  order.resize(n_unique);

  CHECK( array == array_expected);
  CHECK( order == order_expected);

}
