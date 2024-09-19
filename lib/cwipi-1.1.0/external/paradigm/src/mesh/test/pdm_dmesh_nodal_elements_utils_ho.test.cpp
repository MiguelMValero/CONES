#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_logging.h"

// double coord_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};


MPI_TEST_CASE("[pdm_dmesh_nodal_elements_utils] decomposes hexa ho",1) {

  // const PDM_g_num_t n_vtx            = 12;
  // const PDM_g_num_t n_cell           = 2;
  // const int         n_hexa_section_1 = 2;
  // PDM_g_num_t connec_hexa_1[16] = {1, 2, 5, 4, 7, 8, 11, 10, // First
  //                                  2, 3, 6, 5, 8, 9, 12, 11};

  // PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  // printf("coucou\n");

  // int toto = 10;
  // int toto_expected = 10;
  // CHECK( toto == toto_expected );

}
