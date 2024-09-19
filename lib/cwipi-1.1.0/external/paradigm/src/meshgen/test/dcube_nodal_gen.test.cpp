#include "doctest/extensions/doctest_mpi.h"
#include <array>
#include <vector>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_dcube_nodal_gen.h"

MPI_TEST_CASE("[dcube_nodal_gen] - 1p - hexahedron",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 2, 5, 4, 10, 11, 14, 13, 2, 3, 6, 5, 11, 12, 15, 14, 4, 5, 8, 7, 13, 14, 17, 16, 5, 6, 9, 8, 14, 15, 18, 17, 10, 11, 14, 13, 19, 20, 23, 22, 11, 12, 15, 14, 20, 21, 24, 23, 13, 14, 17, 16, 22, 23, 26, 25, 14, 15, 18, 17, 23, 24, 27, 26} //HEXA
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {2, 1, 4, 5, 3, 2, 5, 6, 5, 4, 7, 8, 6, 5, 8, 9, 19, 20, 23, 22, 20, 21, 24, 23, 22, 23, 26, 25, 23, 24, 27, 26, 4, 1, 10, 13, 7, 4, 13, 16, 13, 10, 19, 22, 16, 13, 22, 25, 3, 6, 15, 12, 6, 9, 18, 15, 12, 15, 24, 21, 15, 18, 27, 24, 10, 1, 2, 11, 19, 10, 11, 20, 11, 2, 3, 12, 20, 11, 12, 21, 7, 16, 17, 8, 16, 25, 26, 17, 8, 17, 18, 9, 17, 26, 27, 18}
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 8} // HEXA
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 24}, // All QUAD are together
  };

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_HEXA8};

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4,
                                                                   PDM_MESH_NODAL_QUAD4, PDM_MESH_NODAL_QUAD4};

  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_HEXA8,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_HEXA8,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  8);
  CHECK( n_face_abs ==  0);
  CHECK( n_edge_abs ==  0);
  CHECK( n_vtx_abs  == 27);

  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // All QUAD are together

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}

MPI_TEST_CASE("[dcube_nodal_gen] - 1p - prism",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 2, 4, 10, 11, 13, 5, 4, 2, 14, 13, 11, 2, 3, 5, 11, 12, 14, 6, 5, 3, 15, 14, 12, 4, 5, 7, 13, 14, 16, 8, 7, 5, 17, 16, 14, 5, 6, 8, 14, 15, 17, 9, 8, 6, 18, 17, 15, 10, 11, 13, 19, 20, 22, 14, 13, 11, 23, 22, 20, 11, 12, 14, 20, 21, 23, 15, 14, 12, 24, 23, 21, 13, 14, 16, 22, 23, 25, 17, 16, 14, 26, 25, 23, 14, 15, 17, 23, 24, 26, 18, 17, 15, 27, 26, 24}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 4, 2, 5, 2, 4, 2, 5, 3, 6, 3, 5, 4, 7, 5, 8, 5, 7, 5, 8, 6, 9, 6, 8, 19, 20, 22, 23, 22, 20, 20, 21, 23, 24, 23, 21, 22, 23, 25, 26, 25, 23, 23, 24, 26, 27, 26, 24},
    {4, 1, 10, 13, 7, 4, 13, 16, 13, 10, 19, 22, 16, 13, 22, 25, 3, 6, 15, 12, 6, 9, 18, 15, 12, 15, 24, 21, 15, 18, 27, 24, 10, 1, 2, 11, 19, 10, 11, 20, 11, 2, 3, 12, 20, 11, 12, 21, 7, 16, 17, 8, 16, 25, 26, 17, 8, 17, 18, 9, 17, 26, 27, 18},
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 16}, // PRISM
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 16}, //All TRI
    {0, 16}, //All QUAD
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_PRISM6};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_QUAD4};

  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_PRISM6,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_PRISM6,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs == 16);
  CHECK( n_face_abs ==  0);
  CHECK( n_edge_abs ==  0);
  CHECK( n_vtx_abs  == 27);
  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 2); // QUAD + TRI

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}



MPI_TEST_CASE("[dcube_nodal_gen] - 1p - Tetrahedron ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_vol = {
    {1, 2, 4, 10, 11, 10, 14, 2, 13, 14, 10, 4, 5, 4, 2, 14, 2, 14, 4, 10, 2, 3, 6, 12, 6, 15, 14, 12, 5, 2, 6, 14, 2, 12, 14, 11, 2, 12, 6, 14, 4, 5, 8, 14, 8, 17, 16, 14, 7, 4, 8, 16, 4, 14, 16, 13, 4, 14, 8, 16, 5, 6, 8, 14, 15, 14, 18, 6, 17, 18, 14, 8, 9, 8, 6, 18, 6, 18, 8, 14, 10, 11, 14, 20, 14, 23, 22, 20, 13, 10, 14, 22, 10, 20, 22, 19, 10, 20, 14, 22, 11, 12, 14, 20, 21, 20, 24, 12, 23, 24, 20, 14, 15, 14, 12, 24, 12, 24, 14, 20, 13, 14, 16, 22, 23, 22, 26, 14, 25, 26, 22, 16, 17, 16, 14, 26, 14, 26, 16, 22, 14, 15, 18, 24, 18, 27, 26, 24, 17, 14, 18, 26, 14, 24, 26, 23, 14, 24, 18, 26}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 4, 2, 5, 2, 4, 3, 2, 6, 5, 6, 2, 5, 4, 8, 7, 8, 4, 5, 8, 6, 9, 6, 8, 19, 20, 22, 23, 22, 20, 21, 24, 20, 23, 20, 24, 23, 26, 22, 25, 22, 26, 23, 24, 26, 27, 26, 24, 1, 10, 4, 13, 4, 10, 7, 4, 16, 13, 16, 4, 13, 10, 22, 19, 22, 10, 13, 22, 16, 25, 16, 22, 3, 6, 12, 15, 12, 6, 9, 18, 6, 15, 6, 18, 15, 24, 12, 21, 12, 24, 15, 18, 24, 27, 24, 18, 1, 2, 10, 11, 10, 2, 19, 10, 20, 11, 20, 10, 11, 2, 12, 3, 12, 2, 11, 12, 20, 21, 20, 12, 7, 16, 8, 17, 8, 16, 25, 26, 16, 17, 16, 26, 17, 18, 8, 9, 8, 18, 17, 26, 18, 27, 18, 26}
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_vol = {
    {0, 40}, // PRISM
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0,  48}, // All TRI
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_vol = {PDM_MESH_NODAL_TETRA4};

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_TRIA3};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_TETRA4,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_TETRA4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs == 40);
  CHECK( n_face_abs ==  0);
  CHECK( n_edge_abs ==  0);
  CHECK( n_vtx_abs  == 27);

  int n_section_vol = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  CHECK( n_section_vol == 1); // HEXA + 6 * QUAD

  int* sections_id_vol = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC);

  for(int i_section = 0; i_section < n_section_vol; ++i_section) {

    int id_section = sections_id_vol[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_vol[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_vol[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_vol [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // All TRI

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }


  PDM_dcube_nodal_gen_free(dcube);
}


MPI_TEST_CASE("[dcube_nodal_gen] - 1p - quad ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 2, 5, 4, 2, 3, 6, 5, 4, 5, 8, 7, 5, 6, 9, 8}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_ridge = {
    {1, 2, 2, 3, 8, 7, 9, 8, 4, 1, 7, 4, 3, 6, 6, 9}
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 4} // QUAD
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_ridge = {
    {0, 8}, // All BAR
  };

  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf  = {PDM_MESH_NODAL_QUAD4};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_ridge = {PDM_MESH_NODAL_BAR2};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_QUAD4,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  0);
  CHECK( n_face_abs ==  4);
  CHECK( n_edge_abs ==  0);
  CHECK( n_vtx_abs  ==  9);

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_ridge = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  CHECK( n_section_ridge == 1); // All BAR

  int* sections_id_ridge = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  for(int i_section = 0; i_section < n_section_ridge; ++i_section) {

    int id_section = sections_id_ridge[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_ridge[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_ridge[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_ridge [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}


MPI_TEST_CASE("[dcube_nodal_gen] - 1p - Triangle",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected_surf = {
    {1, 2, 4, 5, 4, 2, 2, 3, 5, 6, 5, 3, 4, 5, 7, 8, 7, 5, 5, 6, 8, 9, 8, 6}
  };

  std::vector<std::vector<PDM_g_num_t>> connec_expected_ridge = {
    {1, 2, 2, 3, 8, 7, 9, 8, 4, 1, 7, 4, 3, 6, 6, 9}
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_surf = {
    {0, 8}, // TRI
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected_ridge = {
    {0, 8}, // All BAR
  };
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_surf = {PDM_MESH_NODAL_TRIA3};
  std::vector<PDM_Mesh_nodal_elt_t> section_types_expexted_ridge = {PDM_MESH_NODAL_BAR2};


  int n_vtx_seg = 3;

  /*PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_TRIA3,
                                                      PDM_OWNERSHIP_KEEP);*/
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  0);
  CHECK( n_face_abs ==  8);
  CHECK( n_edge_abs ==  0);
  CHECK( n_vtx_abs  ==  9);

  int n_section_surf = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  CHECK( n_section_surf == 1); // HEXA + 6 * QUAD

  int* sections_id_surf = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);

  for(int i_section = 0; i_section < n_section_surf; ++i_section) {

    int id_section = sections_id_surf[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_surf[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_surf[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_surf [i_section], n_vtx_per_elmt * dn_elmt);
  }

  int n_section_ridge = PDM_DMesh_nodal_n_section_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  CHECK( n_section_ridge == 1); // All BAR

  int* sections_id_ridge = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);

  for(int i_section = 0; i_section < n_section_ridge; ++i_section) {

    int id_section = sections_id_ridge[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, PDM_GEOMETRY_KIND_RIDGE, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK( section_types_expexted_ridge[i_section] == t_elmt);
    CHECK_EQ_C_ARRAY( distrib, distrib_expected_ridge[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected_ridge [i_section], n_vtx_per_elmt * dn_elmt);
  }

  PDM_dcube_nodal_gen_free(dcube);
}
