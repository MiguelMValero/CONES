#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_block_to_part.h"



MPI_TEST_CASE("[pdm_block_to_part] - 1p - block_to_part",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  static const int darray_size = 6;
  static PDM_g_num_t darray[darray_size] = {5, 4, 6, 1, 3, 2};

  static PDM_g_num_t block_distrib_idx[2] = {0, 7};

  SUBCASE("map with 1 partition ") {
    static const int n_part         = 1;
    static const int n_elmts_part_0 = 6;
    PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {1, 2, 4, 3, 5, 6};

    // Convenient array
    int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
    PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    ln_to_gn[0] = ln_to_gn_p0;
    n_elmts[0]  = n_elmts_part_0;

    // Setup the exchange protocol
    PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_idx,
                                 (const PDM_g_num_t**)  ln_to_gn,
                                                        n_elmts,
                                                        n_part,
                                                        pdm_comm);

    // Allocate the partition array
    PDM_g_num_t** parray = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    for(int i_part = 0; i_part < n_part; i_part++){
      parray[i_part] = (PDM_g_num_t *) malloc( n_elmts[i_part] * sizeof(PDM_g_num_t));
    }

    // Exchange
    int stride_one = 1;
    PDM_block_to_part_exch_in_place(btp, sizeof(PDM_g_num_t), PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                   (void*) darray,
                           NULL,
                  (void**) parray);

    // Check
    static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {5, 4, 1, 6, 3, 2};

    CHECK_EQ_C_ARRAY(parray[0], parray_expected_p0, n_elmts_part_0);
    MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected_p0, n_elmts_part_0);

    // Free
    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
    }
    free(parray);
    free(ln_to_gn);
    free(n_elmts);
    PDM_block_to_part_free(btp);

  }

  SUBCASE("map with 2 partition ") {

    static const int n_part         = 2;
    static const int n_elmts_part_0 = 3;
    static const int n_elmts_part_1 = 3;
    PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {1, 2, 4};
    PDM_g_num_t ln_to_gn_p1[n_elmts_part_1] = {3, 5, 6};

    // Convenient array
    int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
    PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    ln_to_gn[0] = ln_to_gn_p0;
    ln_to_gn[1] = ln_to_gn_p1;
    n_elmts[0]  = n_elmts_part_0;
    n_elmts[1]  = n_elmts_part_1;

    // Setup the exchange protocol
    PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_idx,
                                 (const PDM_g_num_t**)  ln_to_gn,
                                                        n_elmts,
                                                        n_part,
                                                        pdm_comm);

    // Allocate the partition array
    PDM_g_num_t** parray = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    for(int i_part = 0; i_part < n_part; i_part++){
      parray[i_part] = (PDM_g_num_t *) malloc( n_elmts[i_part] * sizeof(PDM_g_num_t));
    }

    // Exchange
    int stride_one = 1;
    PDM_block_to_part_exch_in_place(btp, sizeof(PDM_g_num_t), PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
                   (void*) darray,
                           NULL,
                  (void**) parray);

    // Check
    static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {5, 4, 1};
    static PDM_g_num_t parray_expected_p1[n_elmts_part_1] = {6, 3, 2};

    MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected_p0, n_elmts_part_0);
    MPI_CHECK_EQ_C_ARRAY(0, parray[1], parray_expected_p1, n_elmts_part_1);

    // Free
    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
    }
    free(parray);
    free(ln_to_gn);
    free(n_elmts);
    PDM_block_to_part_free(btp);
  }

}



MPI_TEST_CASE("[pdm_block_to_part] - 2p - block_to_part",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static const int darray_size = 3;
  static PDM_g_num_t darray_p0[darray_size] = {5, 4, 6};
  static PDM_g_num_t darray_p1[darray_size] = {1, 3, 2};

  static PDM_g_num_t block_distrib_idx[3] = {0, 3, 7};

  static const int n_part         = 1;
  static const int n_elmts_part_0 = 3;
  static const int n_elmts_part_1 = 3;
  PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {1, 2, 4};
  PDM_g_num_t ln_to_gn_p1[n_elmts_part_1] = {3, 5, 6};

  // Convenient array
  int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t* darray = NULL;
  if(test_rank == 0){
    n_elmts [0] = n_elmts_part_0;
    ln_to_gn[0] = ln_to_gn_p0;
    darray      = darray_p0;
  } else {
    n_elmts [0] = n_elmts_part_1;
    ln_to_gn[0] = ln_to_gn_p1;
    darray      = darray_p1;
  }

  // Setup the exchange protocol
  PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_idx,
                               (const PDM_g_num_t**)  ln_to_gn,
                                                      n_elmts,
                                                      n_part,
                                                      pdm_comm);

  // Allocate the partition array
  PDM_g_num_t** parray = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; i_part++){
    parray[i_part] = (PDM_g_num_t *) malloc( n_elmts[i_part] * sizeof(PDM_g_num_t));
  }

  // Exchange
  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp, sizeof(PDM_g_num_t), PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                 (void*) darray,
                         NULL,
                (void**) parray);

  // Check
  static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {5, 4, 1};
  static PDM_g_num_t parray_expected_p1[n_elmts_part_1] = {6, 3, 2};

  MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected_p0, n_elmts_part_0);
  MPI_CHECK_EQ_C_ARRAY(1, parray[0], parray_expected_p1, n_elmts_part_1);

  // Free
  for(int i_part = 0; i_part < n_part; i_part++){
    free(parray[i_part]);
  }
  free(parray);
  free(ln_to_gn);
  free(n_elmts);
  PDM_block_to_part_free(btp);

}
