#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_logging.h"

/*
 *  Use case
 *
 *       3           4           2
 *       +-----------+-----------+
 *       |           |           |
 *       |           |           |
 *       |           |           |
 *       +-----------+-----------+
 *       5           1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */
MPI_TEST_CASE("[pdm_part_mesh_nodal_elmts] Constructor",1) {
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

  // const PDM_g_num_t n_vtx            = 9;
  // const PDM_g_num_t n_face           = 8;
  // const PDM_g_num_t n_ridge          = 8;
  const int         n_tri_section_1  = 8;
  const int         n_bar_section_1  = 8;

  int connec_tri_1[24] = {6, 8, 9,
                          9, 5, 6,
                          2, 8, 1,
                          9, 3, 4,
                          6, 7, 8,
                          9, 4, 5,
                          2, 9, 8,
                          9, 2, 3};

  int connec_bar_1[16] = {1, 2,
                          2, 3,
                          4, 5,
                          3, 4,
                          6, 7,
                          5, 6,
                          8, 1,
                          7, 8};

  int n_part = 1;
  // int n_group_elmt = 1;
  // int dgroup_elmt_idx[2] = {0, 8};
  // // PDM_g_num_t dgroup_elmt[8] = {9, 10, 11, 12, 13, 14, 15, 16};
  // PDM_g_num_t dgroup_elmt[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_part_mesh_nodal_t* pmn  = PDM_part_mesh_nodal_create(2, n_part, pdm_comm);

  PDM_part_mesh_nodal_elmts_t* pelmts_surf  = PDM_part_mesh_nodal_elmts_create(2, n_part, pdm_comm);
  PDM_part_mesh_nodal_elmts_t* pelmts_ridge = PDM_part_mesh_nodal_elmts_create(1, n_part, pdm_comm);

  int tri_section_1 = PDM_part_mesh_nodal_elmts_add(pelmts_surf , PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_part_mesh_nodal_elmts_add(pelmts_ridge, PDM_MESH_NODAL_BAR2 );

  int i_part = 0;
  PDM_part_mesh_nodal_elmts_std_set(pelmts_surf,
                                    tri_section_1,
                                    i_part,
                                    n_tri_section_1,
                                    connec_tri_1,
                                    NULL,
                                    NULL,
                                    NULL,
                                    PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_elmts_std_set(pelmts_ridge,
                                    bar_section_1,
                                    i_part,
                                    n_bar_section_1,
                                    connec_bar_1,
                                    NULL,
                                    NULL,
                                    NULL,
                                    PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pelmts_surf );
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pelmts_ridge);

  // PDM_part_mesh_nodal_elmts_free(pelmts_surf);
  // PDM_part_mesh_nodal_elmts_free(pelmts_ridge);

  PDM_part_mesh_nodal_free(pmn);
}
