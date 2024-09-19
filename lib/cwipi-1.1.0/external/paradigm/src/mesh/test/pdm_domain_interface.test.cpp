#include "doctest/extensions/doctest_mpi.h"
#include <vector>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"

// #include "pdm_domain_interface.c"


// MPI_TEST_CASE("[2p] interface_to_graph", 2) {

//   PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
//   int i_rank, n_rank;
//   PDM_MPI_Comm_rank (pdm_comm, &i_rank);
//   PDM_MPI_Comm_size (pdm_comm, &n_rank);

//   int n_interface = 2;

//   std::vector<int> interface_dn;
//   std::vector<PDM_g_num_t*> interface_ids_ptr;
//   std::vector<int*> interface_dom_ptr;

//   std::vector<int> expected_graph_idx;
//   std::vector<PDM_g_num_t> expected_graph_ids;
//   std::vector<int> expected_graph_dom;
//   PDM_domain_interface_mult_t multizone_intrf;

//   std::vector<std::vector<int>> interface_dom; //Only to keep memory alive
//   std::vector<std::vector<PDM_g_num_t>> interface_ids;
//   SUBCASE("mult interface") {
//     if (i_rank == 0) {
//       interface_dn = {3,2};
//       interface_ids = {{4,1, 8,5, 12,9}, {13,1, 2,14}};
//       interface_dom = {{0,1, 0,1, 0,1}, {2,1, 1,2}};
//     }
//     else {
//       interface_dn = {1,2};
//       interface_ids = {{13,16}, {15,3, 16,4}};
//       interface_dom = {{1,0}, {2,1, 2,1}};
//     }
//     multizone_intrf = PDM_DOMAIN_INTERFACE_MULT_YES;
//   }
//   SUBCASE("simple interface") {
//     if (i_rank == 0) {
//       interface_dn = {3,2};
//       interface_ids = {{4,1, 8,5, 12,9}, {1,13, 2,14}};
//       interface_dom = {{0,1}, {1,2}};
//     }
//     else {
//       interface_dn = {1,2};
//       interface_ids = {{16,13}, {3,15, 4,16}};
//       interface_dom = {{0,1}, {1,2}};
//     }
//     multizone_intrf = PDM_DOMAIN_INTERFACE_MULT_NO;
//   }

//   if (i_rank == 0) {
//     expected_graph_idx = {0,3,5,7};
//     expected_graph_ids = {4,1,13, 8,5, 12,9};
//     expected_graph_dom = {0,1,2, 0,1, 0,1};
//   }
//   else {
//     expected_graph_idx = {0,2,4,6,8};
//     expected_graph_ids = {16,13, 2,14, 3,15, 4,16};
//     expected_graph_dom = {0,1, 1,2, 1,2, 1,2};
//   }
//   for (int k = 0; k < n_interface; ++k) {
//     interface_dom_ptr.push_back(interface_dom[k].data());
//     interface_ids_ptr.push_back(interface_ids[k].data());
//   }

//   int         *graph_idx = NULL;
//   PDM_g_num_t *graph_ids = NULL;
//   int         *graph_dom = NULL;

//   int graph_dn = _interface_to_graph(n_interface, multizone_intrf,
//       interface_dn.data(), interface_ids_ptr.data(), interface_dom_ptr.data(),
//       &graph_idx, &graph_ids, &graph_dom, pdm_comm);

//   CHECK_EQ_C_ARRAY(graph_idx, expected_graph_idx.data(), graph_dn);
//   CHECK_EQ_C_ARRAY(graph_ids, expected_graph_ids.data(), graph_idx[graph_dn]);
//   CHECK_EQ_C_ARRAY(graph_dom, expected_graph_dom.data(), graph_idx[graph_dn]);

//   free(graph_idx);
//   free(graph_ids);
//   free(graph_dom);
// }


