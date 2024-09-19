#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_distrib.h"

MPI_TEST_CASE("[pdm_distrib] - 1p - PDM_distrib_compute",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t* distrib = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  int dnelmt = 10;

  SUBCASE("No offset ") {
    int offset = 0;
    PDM_distrib_compute(dnelmt, distrib, offset, pdm_comm);

    std::vector<PDM_g_num_t> distrib_expected = { 1, 11};
    MPI_CHECK_EQ_C_ARRAY(0, distrib, distrib_expected.data(), n_rank+1);
  }

  SUBCASE(" With offset") {

    int offset = 20;
    PDM_distrib_compute(dnelmt, distrib, offset, pdm_comm);

    std::vector<PDM_g_num_t> distrib_expected = { 21, 31};
    MPI_CHECK_EQ_C_ARRAY(0, distrib, distrib_expected.data(), n_rank+1);
  }

  free(distrib);
}


MPI_TEST_CASE("[pdm_distrib] - 2p - PDM_distrib_compute",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t* distrib = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  int dnelmt = -1;
  if(i_rank == 0) { dnelmt = 10;}
  if(i_rank == 1) { dnelmt = 5; }

  SUBCASE("No offset ") {
    int offset = 0;
    PDM_distrib_compute(dnelmt, distrib, offset, pdm_comm);

    std::vector<PDM_g_num_t> distrib_expected = { 1, 11, 16};
    CHECK_EQ_C_ARRAY(distrib, distrib_expected.data(), n_rank+1);
  }

  SUBCASE(" With offset") {

    int offset = 20;
    PDM_distrib_compute(dnelmt, distrib, offset, pdm_comm);

    std::vector<PDM_g_num_t> distrib_expected = { 21, 31, 36};
    CHECK_EQ_C_ARRAY(distrib, distrib_expected.data(), n_rank+1);
  }

  free(distrib);
}


MPI_TEST_CASE("[pdm_distrib] - 1p - PDM_compute_entity_distribution",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int dnelmt = 10;

  PDM_g_num_t* distrib = PDM_compute_entity_distribution(pdm_comm, dnelmt);

  std::vector<PDM_g_num_t> distrib_expected = { 0, 10};
  MPI_CHECK_EQ_C_ARRAY(0, distrib, distrib_expected.data(), n_rank+1);

  free(distrib);

}

MPI_TEST_CASE("[pdm_distrib] - 2p - PDM_compute_entity_distribution",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int dnelmt = -1;
  if(i_rank == 0) { dnelmt = 10;}
  if(i_rank == 1) { dnelmt = 5; }

  PDM_g_num_t* distrib = PDM_compute_entity_distribution(pdm_comm, dnelmt);

  std::vector<PDM_g_num_t> distrib_expected = { 0, 10, 15};
  CHECK_EQ_C_ARRAY(distrib, distrib_expected.data(), n_rank+1);

  free(distrib);

}


MPI_TEST_CASE("[pdm_distrib] - 1p - PDM_compute_uniform_entity_distribution",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_g_elmt = 9;
  PDM_g_num_t* distrib = PDM_compute_uniform_entity_distribution(pdm_comm, n_g_elmt);

  int dn_elemt = PDM_compute_uniform_dn_entity(pdm_comm, n_g_elmt);

  MPI_CHECK(0, dn_elemt == 9);

  std::vector<PDM_g_num_t> distrib_expected = { 0, 9};
  MPI_CHECK_EQ_C_ARRAY(0, distrib, distrib_expected.data(), n_rank+1);

  free(distrib);
}

MPI_TEST_CASE("[pdm_distrib] - 2p - PDM_compute_uniform_entity_distribution",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n_g_elmt = 9;
  PDM_g_num_t* distrib = PDM_compute_uniform_entity_distribution(pdm_comm, n_g_elmt);

  int dn_elemt = PDM_compute_uniform_dn_entity(pdm_comm, n_g_elmt);

  MPI_CHECK(0, dn_elemt == 5);
  MPI_CHECK(1, dn_elemt == 4);

  std::vector<PDM_g_num_t> distrib_expected = { 0, 5, 9};
  CHECK_EQ_C_ARRAY(distrib, distrib_expected.data(), n_rank+1);


}
