#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_para_graph_dual.h"



MPI_TEST_CASE("[pdm_para_graph_dual] - 1p - dual from arc2node", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_g_num_t* dual_graph_idx;
  PDM_g_num_t* dual_graph;
  int*         dcell_face_idx = NULL;
  PDM_g_num_t* dcell_face = NULL;

  SUBCASE("Simple graph") {

    /* Simple graph with 5 vertices connected by 6 edges (no vtx to nothing edge) */
    static const int dNnode = 5;
    static const int dNarc  = 6;
    static PDM_g_num_t node_distribution[2] = {0, dNnode};
    static PDM_g_num_t arc_distribution[2]  = {0, dNarc};
    static PDM_g_num_t darc2node[2*6] = {
      1, 2, 2, 5, 5, 4, 4, 3, 4, 2, 3, 1};

    PDM_para_graph_dual_from_arc2node(pdm_comm,
                                      node_distribution,
                                      arc_distribution,
                                      darc2node,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph,
                                      0,
                      (int        **) &dcell_face_idx,
                      (PDM_g_num_t**) &dcell_face);

    // Check
    static PDM_g_num_t expc_graph_idx[dNnode+1] = {0, 2, 5, 7, 10, 12};
    CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, dNnode+1);
    // Nb graph starts at 0
    static PDM_g_num_t expc_graph[12] = {1, 2, 0, 3, 4, 0, 3, 1, 2, 4, 1, 3};
    CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 12);
    CHECK(dcell_face_idx == nullptr);
    CHECK(dcell_face == nullptr);
  }

  SUBCASE("Hexa mesh 8 cells") {

    /* Comes from 8 hexa cells connected by faces. Note that face_cell connect
    some cells to the boundary (0) */
    static const int dcell_size = 8;
    static const int dface_size = 36;

    static PDM_g_num_t cell_distribution[2] = {0, dcell_size};
    static PDM_g_num_t face_distribution[2] = {0, dface_size};
    static PDM_g_num_t dface_cell[2*dface_size] = {
      1, 0, 2, 0, 3, 0, 4, 0, 1, 5, 2, 6, 3, 7, 4, 8, 5, 0, 6, 0, 7, 0, 8, 0,
      1, 0, 3, 0, 5, 0, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 0, 4, 0, 6, 0, 8, 0,
      1, 0, 5, 0, 2, 0, 6, 0, 1, 3, 5, 7, 2, 4, 6, 8, 3, 0, 7, 0, 4, 0, 8, 0};

    PDM_para_graph_dual_from_arc2node(pdm_comm,
                                      cell_distribution,
                                      face_distribution,
                                      dface_cell,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph,
                                      1,
                      (int        **) &dcell_face_idx,
                      (PDM_g_num_t**) &dcell_face);

    // Check
    static PDM_g_num_t expc_graph_idx[dcell_size+1] = {0, 3, 6, 9, 12, 15, 18, 21, 24};
    CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, dcell_size+1);
    // Nb graph starts at 0
    static PDM_g_num_t expc_graph[24] = {
     1, 2, 4, 0, 3, 5, 0, 3, 6, 1, 2, 7, 0, 5, 6, 1, 4, 7, 2, 4, 7, 3, 5, 6};
    CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 24);

    static int expc_cell_face_idx[dcell_size+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
    CHECK_EQ_C_ARRAY(dcell_face_idx, expc_cell_face_idx, dcell_size+1);
    static PDM_g_num_t expc_cell_face[48] = {
     1, 5, 13, 17, 25, 29, 2, 6,-17, 21, 27, 31, 3, 7, 14, 18,
     -29, 33, 4, 8,-18, 22,-31, 35,-5, 9, 15, 19, 26, 30,-6,
     10,-19, 23, 28, 32,-7, 11, 16, 20,-30, 34,-8, 12,-20, 24,-32, 36};
    CHECK_EQ_C_ARRAY(dcell_face, expc_cell_face, 48);
   }

  free(dual_graph_idx);
  free(dual_graph);
  if (dcell_face_idx != NULL)
    free(dcell_face_idx);
  if (dcell_face != NULL)
    free(dcell_face);
}

MPI_TEST_CASE("[pdm_para_graph_dual] - 3p - dual from arc2node", 3) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_g_num_t* dual_graph_idx = NULL;
  PDM_g_num_t* dual_graph = NULL;
  int*         dcell_face_idx = NULL;
  PDM_g_num_t* dcell_face = NULL;

  SUBCASE("Hexa mesh 8 cells") {

    /* Comes from 8 hexa cells connected by faces. Note that face_cell connect
    some cells to the boundary (0) */

    static PDM_g_num_t cell_distribution[4] = {0, 3, 6, 8};
    static PDM_g_num_t face_distribution[4] = {0, 12, 24, 36};

    static PDM_g_num_t dface_cell_p0[2*12] = {
        1, 0, 2, 0, 3, 0, 4, 0, 1, 5, 2, 6, 3, 7, 4, 8, 5, 0, 6, 0, 7, 0, 8, 0};
    static PDM_g_num_t dface_cell_p1[2*12] = {
        1, 0, 3, 0, 5, 0, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 0, 4, 0, 6, 0, 8, 0};
    static PDM_g_num_t dface_cell_p2[2*12] = {
        1, 0, 5, 0, 2, 0, 6, 0, 1, 3, 5, 7, 2, 4, 6, 8, 3, 0, 7, 0, 4, 0, 8, 0};
    PDM_g_num_t* dface_cell = NULL;
    // static PDM_g_num_t dface_cell[2*12]
    if (test_rank == 0) {
        dface_cell = dface_cell_p0;
    }
    else if (test_rank == 1) {
        dface_cell = dface_cell_p1;
    }
    else if (test_rank == 2) {
        dface_cell = dface_cell_p2;
    }

    PDM_para_graph_dual_from_arc2node(pdm_comm,
                                      cell_distribution,
                                      face_distribution,
                                      dface_cell,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph,
                                      1,
                      (int        **) &dcell_face_idx,
                      (PDM_g_num_t**) &dcell_face);

    // Check
    if (test_rank == 0) {
      static PDM_g_num_t expc_graph_idx[3+1] = {0, 3, 6, 9};
      static PDM_g_num_t expc_graph[9] = {1, 2, 4, 0, 3, 5, 0, 3, 6};
      CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, 3+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 9);
      static int expc_dcellface_idx[3+1] = {0, 6, 12, 18};
      static PDM_g_num_t expc_dcellface[18] = {
        1, 5, 13, 17, 25, 29, 2, 6, -17, 21, 27, 31, 3, 7, 14, 18, -29, 33};
      CHECK_EQ_C_ARRAY(dcell_face_idx, expc_dcellface_idx, 3+1);
      CHECK_EQ_C_ARRAY(dcell_face, expc_dcellface, 18);
    }
    else if (test_rank == 1){
      static PDM_g_num_t expc_graph_idx[3+1] = {0, 3, 6, 9};
      static PDM_g_num_t expc_graph[9] = {1, 2, 7, 0, 5, 6, 1, 4, 7};
      MPI_CHECK_EQ_C_ARRAY(1, dual_graph_idx, expc_graph_idx, 3+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 9);
      static int expc_dcellface_idx[3+1] = {0, 6, 12, 18};
      static PDM_g_num_t expc_dcellface[18] = {
        4, 8, -18, 22, -31, 35, -5, 9, 15, 19, 26, 30, -6, 10, -19, 23, 28, 32};
      CHECK_EQ_C_ARRAY(dcell_face_idx, expc_dcellface_idx, 3+1);
      CHECK_EQ_C_ARRAY(dcell_face, expc_dcellface, 18);
    }
    else if (test_rank ==2){
      static PDM_g_num_t expc_graph_idx[2+1] = {0, 3, 6};
      static PDM_g_num_t expc_graph[6] = {2, 4, 7, 3, 5, 6};
      MPI_CHECK_EQ_C_ARRAY(2, dual_graph_idx, expc_graph_idx, 2+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 6);
      static int expc_dcellface_idx[2+1] = {0, 6, 12};
      static PDM_g_num_t expc_dcellface[12] = {
        -7, 11, 16, 20, -30, 34, -8, 12, -20, 24, -32, 36};
      CHECK_EQ_C_ARRAY(dcell_face_idx, expc_dcellface_idx, 2+1);
      CHECK_EQ_C_ARRAY(dcell_face, expc_dcellface, 12);
    }
   }

  if(dual_graph_idx != NULL) free(dual_graph_idx);
  if(dual_graph     != NULL) free(dual_graph    );
  if(dcell_face_idx != NULL) free(dcell_face_idx);
  if(dcell_face     != NULL) free(dcell_face    );
}



MPI_TEST_CASE("[pdm_para_graph_dual] - 1p - dual from node2arc", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_g_num_t* dual_graph_idx = NULL;
  PDM_g_num_t* dual_graph = NULL;

  SUBCASE("Simple graph") {

    /* Simple graph with 5 vertices connected by 6 edges (no vtx to nothing edge) */
    static const int dNnode = 5;
    static const int dNarc  = 6;
    static PDM_g_num_t node_distribution[2] = {0, dNnode};
    static PDM_g_num_t arc_distribution[2]  = {0, dNarc};
    static int         dnode2arc_idx[dNnode+1] = {0, 2, 5, 7, 10, 12};
    static PDM_g_num_t dnode2arc[12] = {1, -6, -1, 2, -5, -4, 6, -3, 4, 5, -2, 3};

    PDM_para_graph_dual_from_node2arc(pdm_comm,
                                      node_distribution,
                                      arc_distribution,
                                      dnode2arc_idx,
                                      dnode2arc,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph);

    // Check
    static PDM_g_num_t expc_graph_idx[dNnode+1] = {0, 2, 5, 7, 10, 12};
    CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, dNnode+1);
    // Nb graph starts at 0
    static PDM_g_num_t expc_graph[12] = {1, 2, 0, 3, 4, 0, 3, 1, 2, 4, 1, 3};
    CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 12);
  }

  SUBCASE("Hexa mesh 8 cells") {

    /* Comes from 8 hexa cells connected by faces. Note that face_cell connect
    some cells to the boundary (0) */
    static const int dcell_size = 8;
    static const int dface_size = 36;

    static PDM_g_num_t cell_distribution[2] = {0, dcell_size};
    static PDM_g_num_t face_distribution[2] = {0, dface_size};
    static int         dnode2arc_idx[dcell_size+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
    static PDM_g_num_t dnode2arc[48] = {
     1, 5, 13, 17, 25, 29, 2, 6,-17, 21, 27, 31, 3, 7, 14, 18,
     -29, 33, 4, 8,-18, 22,-31, 35,-5, 9, 15, 19, 26, 30,-6,
     10,-19, 23, 28, 32,-7, 11, 16, 20, -30, 34,-8, 12,-20, 24,-32, 36};

    PDM_para_graph_dual_from_node2arc(pdm_comm,
                                      cell_distribution,
                                      face_distribution,
                                      dnode2arc_idx,
                                      dnode2arc,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph);


    // Check
    static PDM_g_num_t expc_graph_idx[dcell_size+1] = {0, 3, 6, 9, 12, 15, 18, 21, 24};
    CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, dcell_size+1);
    // Nb graph starts at 0
    static PDM_g_num_t expc_graph[24] = {
     1, 2, 4, 0, 3, 5, 0, 3, 6, 1, 2, 7, 0, 5, 6, 1, 4, 7, 2, 4, 7, 3, 5, 6};
    CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 24);
   }

  free(dual_graph_idx);
  free(dual_graph);
}

MPI_TEST_CASE("[pdm_para_graph_dual] - 3p - dual from node2arc", 3) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_g_num_t* dual_graph_idx = NULL;
  PDM_g_num_t* dual_graph = NULL;

  SUBCASE("Hexa mesh 8 cells") {

    /* Comes from 8 hexa cells connected by faces. Note that face_cell connect
    some cells to the boundary (0) */
    //static const int dcell_size = 8;
    //static const int dface_size = 36;

    static PDM_g_num_t cell_distribution[4] = {0, 3, 6, 8};
    static PDM_g_num_t face_distribution[4] = {0, 12, 24, 36};

    static int dcellface_idx_p0[3+1] = {0, 6, 12, 18};
    static PDM_g_num_t dcellface_p0[18] = {
     1, 5, 13, 17, 25, 29, 2, 6, -17, 21, 27, 31, 3, 7, 14, 18, -29, 33};

    static int dcellface_idx_p1[3+1] = {0, 6, 12, 18};
    static PDM_g_num_t dcellface_p1[18] = {
      4, 8, -18, 22, -31, 35, -5, 9, 15, 19, 26, 30, -6, 10, -19, 23, 28, 32};

    static int dcellface_idx_p2[2+1] = {0, 6, 12};
    static PDM_g_num_t dcellface_p2[12] = {
      -7, 11, 16, 20, -30, 34, -8, 12, -20, 24, -32, 36};

    int         *dcellface_idx = NULL;
    PDM_g_num_t *dcellface = NULL;

    if (test_rank == 0){
      dcellface_idx = dcellface_idx_p0;
      dcellface     = dcellface_p0;
    }
    else if(test_rank == 1){
      dcellface_idx = dcellface_idx_p1;
      dcellface     = dcellface_p1;
    }
    else if(test_rank == 2){
      dcellface_idx = dcellface_idx_p2;
      dcellface     = dcellface_p2;
    }

    PDM_para_graph_dual_from_node2arc(pdm_comm,
                                      cell_distribution,
                                      face_distribution,
                                      dcellface_idx,
                                      dcellface,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph);


    // Check
    if (test_rank == 0) {
      static PDM_g_num_t expc_graph_idx[3+1] = {0, 3, 6, 9};
      static PDM_g_num_t expc_graph[9] = {1, 2, 4, 0, 3, 5, 0, 3, 6};
      CHECK_EQ_C_ARRAY(dual_graph_idx, expc_graph_idx, 3+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 9);
    }
    else if (test_rank == 1){
      static PDM_g_num_t expc_graph_idx[3+1] = {0, 3, 6, 9};
      static PDM_g_num_t expc_graph[9] = {1, 2, 7, 0, 5, 6, 1, 4, 7};
      MPI_CHECK_EQ_C_ARRAY(1, dual_graph_idx, expc_graph_idx, 3+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 9);
    }
    else if (test_rank ==2){
      static PDM_g_num_t expc_graph_idx[2+1] = {0, 3, 6};
      static PDM_g_num_t expc_graph[6] = {2, 4, 7, 3, 5, 6};
      MPI_CHECK_EQ_C_ARRAY(2, dual_graph_idx, expc_graph_idx, 2+1);
      CHECK_EQ_C_ARRAY(dual_graph, expc_graph, 6);
    }
  }
  free(dual_graph_idx);
  free(dual_graph);
}

MPI_TEST_CASE("[pdm_para_graph_dual] - 3p - Split graph with PARMetis", 3) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  SUBCASE("Hexa mesh 8 cells") {

    /* Comes from 8 hexa cells connected by faces. Note that face_cell connect
    some cells to the boundary (0) */

    static PDM_g_num_t cell_distribution[4] = {0, 3, 6, 8};

    static PDM_g_num_t graph_idx_p0[3+1] = {0, 3, 6, 9};
    static PDM_g_num_t graph_p0[9] = {1, 2, 4, 0, 3, 5, 0, 3, 6};
    static PDM_g_num_t graph_idx_p1[3+1] = {0, 3, 6, 9};
    static PDM_g_num_t graph_p1[9] = {1, 2, 7, 0, 5, 6, 1, 4, 7};
    static PDM_g_num_t graph_idx_p2[2+1] = {0, 3, 6};
    static PDM_g_num_t graph_p2[6] = {2, 4, 7, 3, 5, 6};

    int dn_cell = -1;
    PDM_g_num_t *graph_idx = NULL;
    PDM_g_num_t *graph     = NULL;

    if (test_rank == 0){
      graph_idx = graph_idx_p0;
      graph     = graph_p0;
      dn_cell = 3;
    }
    else if(test_rank == 1){
      graph_idx = graph_idx_p1;
      graph     = graph_p1;
      dn_cell = 3;
    }
    else if(test_rank == 2){
      graph_idx = graph_idx_p2;
      graph     = graph_p2;
      dn_cell = 2;
    }
    int *cell_part = (int *) malloc(dn_cell * sizeof(int));

    SUBCASE("n_part == n_proc") {
      int n_part = 3;
      PDM_para_graph_split (PDM_SPLIT_DUAL_WITH_PARMETIS,
                            cell_distribution,
                            graph_idx,
                            graph,
                            NULL,
                            NULL,
                            n_part,
                            NULL,
                            cell_part,
                            pdm_comm);

      // Check
      if (test_rank == 0) {
        static PDM_g_num_t expc_cell_part[3] = {0, 0, 0};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank == 1){
        static PDM_g_num_t expc_cell_part[3] = {2, 1, 1};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank ==2){
        static PDM_g_num_t expc_cell_part[2] = {1, 2};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
    }
    SUBCASE("n_part < n_proc") {
      int n_part = 2;
      PDM_para_graph_split (PDM_SPLIT_DUAL_WITH_PARMETIS,
                            cell_distribution,
                            graph_idx,
                            graph,
                            NULL,
                            NULL,
                            n_part,
                            NULL,
                            cell_part,
                            pdm_comm);

      // Check
      if (test_rank == 0) {
        static PDM_g_num_t expc_cell_part[3] = {1, 1, 1};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank == 1){
        static PDM_g_num_t expc_cell_part[3] = {1, 0, 0};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank ==2){
        static PDM_g_num_t expc_cell_part[2] = {0, 0};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
    }
    SUBCASE("n_part > n_proc") {
      int n_part = 4;
      PDM_para_graph_split (PDM_SPLIT_DUAL_WITH_PARMETIS,
                            cell_distribution,
                            graph_idx,
                            graph,
                            NULL,
                            NULL,
                            n_part,
                            NULL,
                            cell_part,
                            pdm_comm);

      // Check
      if (test_rank == 0) {
        static PDM_g_num_t expc_cell_part[3] = {2, 2, 3};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank == 1){
        static PDM_g_num_t expc_cell_part[3] = {3, 0, 1};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
      else if (test_rank ==2){
        static PDM_g_num_t expc_cell_part[2] = {0, 1};
        CHECK_EQ_C_ARRAY(cell_part, expc_cell_part, dn_cell);
      }
    }
    free(cell_part);
  }
}
