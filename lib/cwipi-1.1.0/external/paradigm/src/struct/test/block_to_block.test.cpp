#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_block_to_block.h"



MPI_TEST_CASE("[pdm_block_to_block] - 1p - block_to_block",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static PDM_g_num_t block_distrib_ini_idx[2] = {0, 7};
  static PDM_g_num_t block_distrib_end_idx[2] = {0, 7};

  PDM_block_to_block_t* btb = PDM_block_to_block_create(block_distrib_ini_idx,
                                                        block_distrib_end_idx,
                                                        pdm_comm);

  SUBCASE("Exchange stride constante ") {

    int blk_data_ini[7] = {5, 8, 11, 0, 3, 1, 9};

    int *blk_data_end = NULL;
    PDM_block_to_block_exch(btb,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                            blk_data_ini,
                            NULL,
              (void **)    &blk_data_end);

    int blk_data_end_expected[7] = {5, 8, 11, 0, 3, 1, 9};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected, block_distrib_ini_idx[1]);

    free(blk_data_end);

  }

  SUBCASE("Exchange stride constante with mpi_type ") {

    PDM_MPI_Datatype mpi_double_int_type;
    PDM_MPI_Type_create_contiguous(2, PDM_MPI_INT, &mpi_double_int_type);
    PDM_MPI_Type_commit(&mpi_double_int_type);

    int blk_data_ini[14] = {5, 50, 8, 80, 11, 110, 0, 00, 3, 30, 1, 10, 9, 90};

    int *blk_data_end = NULL;
    PDM_block_to_block_exch_with_mpi_type(btb,
                                          PDM_STRIDE_CST_INTERLACED,
                                          mpi_double_int_type,
                                          NULL,
                                          blk_data_ini,
                                          NULL,
                            (void **)    &blk_data_end);

    int blk_data_end_expected[14] = {5, 50, 8, 80, 11, 110, 0, 00, 3, 30, 1, 10, 9, 90};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected, 2 * block_distrib_ini_idx[1]);

    free(blk_data_end);
    PDM_MPI_Type_free(&mpi_double_int_type);
  }

  SUBCASE("Exchange stride variable ") {

    int blk_data_ini[7] = {5, 8, 11, 0, 3, 1, 9};
    int blk_stri_ini[7] = {0, 2,  0, 3, 2, 0, 0};

    int  blk_stri_end[7] = {-1, -1, -1, -1, -1, -1, -1};
    int *blk_data_end = NULL;
    PDM_block_to_block_exch(btb,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            blk_stri_ini,
                            blk_data_ini,
              (int   *)     blk_stri_end,
              (void **)    &blk_data_end);

    int blk_stri_end_expected[7] = {0, 2,  0, 3, 2, 0, 0};
    MPI_CHECK_EQ_C_ARRAY(0, blk_stri_end, blk_stri_end_expected, block_distrib_ini_idx[1]);

    int blk_data_end_expected[7] = {5, 8, 11, 0, 3, 1, 9};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected, block_distrib_ini_idx[1]);


    free(blk_data_end);

  }

  SUBCASE("Exchange stride variable with mpi_type ") {

    PDM_MPI_Datatype mpi_double_int_type;
    PDM_MPI_Type_create_contiguous(2, PDM_MPI_INT, &mpi_double_int_type);
    PDM_MPI_Type_commit(&mpi_double_int_type);

    int blk_stri_ini[7] = {0, 2,  0, 3, 2, 0, 0};
    int blk_data_ini[14] = {5, 50, 8, 80, 11, 110, 0, 00, 3, 30, 1, 10, 9, 90};

    int  blk_stri_end[7] = {-1, -1, -1, -1, -1, -1, -1};
    int *blk_data_end = NULL;
    PDM_block_to_block_exch_with_mpi_type(btb,
                                          PDM_STRIDE_VAR_INTERLACED,
                                          mpi_double_int_type,
                                          blk_stri_ini,
                                          blk_data_ini,
                            (int   *)     blk_stri_end,
                            (void **)    &blk_data_end);

    int blk_stri_end_expected[7] = {0, 2,  0, 3, 2, 0, 0};
    MPI_CHECK_EQ_C_ARRAY(0, blk_stri_end, blk_stri_end_expected, block_distrib_ini_idx[1]);

    int blk_data_end_expected[14] = {5, 50, 8, 80, 11, 110, 0, 00, 3, 30, 1, 10, 9, 90};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected, 2 * block_distrib_ini_idx[1]);

    free(blk_data_end);
    PDM_MPI_Type_free(&mpi_double_int_type);
  }

  PDM_block_to_block_free(btb);

}

MPI_TEST_CASE("[pdm_block_to_block] - 2p - block_to_block",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static PDM_g_num_t block_distrib_ini_idx[3] = {0, 2, 7};
  static PDM_g_num_t block_distrib_end_idx[3] = {0, 4, 7};

  PDM_block_to_block_t* btb = PDM_block_to_block_create(block_distrib_ini_idx,
                                                        block_distrib_end_idx,
                                                        pdm_comm);

  SUBCASE("Exchange stride constante ") {

    int blk_data_ini_p0[2] = {5, 8};
    int blk_data_ini_p1[5] = {11, 0, 3, 1, 9};
    int *blk_data_ini = NULL;
    if(test_rank == 0) {
      blk_data_ini = blk_data_ini_p0;
    } else {
      blk_data_ini = blk_data_ini_p1;
    }

    int *blk_data_end = NULL;
    PDM_block_to_block_exch(btb,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                            blk_data_ini,
                            NULL,
              (void **)    &blk_data_end);

    int blk_data_end_expected_p0[4] = {5, 8, 11, 0};
    int blk_data_end_expected_p1[3] = {3, 1, 9};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected_p0, 4);
    MPI_CHECK_EQ_C_ARRAY(1, blk_data_end, blk_data_end_expected_p1, 3);

    free(blk_data_end);

  }

  SUBCASE("Exchange stride constante with mpi_type ") {

    PDM_MPI_Datatype mpi_double_int_type;
    PDM_MPI_Type_create_contiguous(2, PDM_MPI_INT, &mpi_double_int_type);
    PDM_MPI_Type_commit(&mpi_double_int_type);

    int blk_data_ini_p0[4 ] = {5, 50, 8, 80};
    int blk_data_ini_p1[10] = {11, 110, 0, 0, 3, 30, 1, 10, 9, 90};

    int *blk_data_ini = NULL;
    if(test_rank == 0) {
      blk_data_ini = blk_data_ini_p0;
    } else {
      blk_data_ini = blk_data_ini_p1;
    }

    int *blk_data_end = NULL;
    PDM_block_to_block_exch_with_mpi_type(btb,
                                          PDM_STRIDE_CST_INTERLACED,
                                          mpi_double_int_type,
                                          NULL,
                                          blk_data_ini,
                                          NULL,
                            (void **)    &blk_data_end);

    int blk_data_end_expected_p0[8] = {5, 50, 8, 80, 11, 110, 0, 0};
    int blk_data_end_expected_p1[6] = {3, 30, 1, 10, 9, 90};

    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected_p0, 2 * 4);
    MPI_CHECK_EQ_C_ARRAY(1, blk_data_end, blk_data_end_expected_p1, 2 * 3);

    free(blk_data_end);
    PDM_MPI_Type_free(&mpi_double_int_type);
  }

  SUBCASE("Exchange stride variable ") {

    int blk_stri_ini_p0[2] = {0, 2};
    int blk_stri_ini_p1[5] = {0, 3, 2, 0, 0};

    int blk_data_ini_p0[2] = {5, 8};
    int blk_data_ini_p1[5] = {11, 0, 3, 1, 9};

    int *blk_data_ini = NULL;
    int *blk_stri_ini = NULL;
    if(test_rank == 0) {
      blk_stri_ini = blk_stri_ini_p0;
      blk_data_ini = blk_data_ini_p0;
    } else {
      blk_stri_ini = blk_stri_ini_p1;
      blk_data_ini = blk_data_ini_p1;
    }

    int  blk_stri_end[7] = {-1, -1, -1, -1, -1, -1, -1};
    int *blk_data_end = NULL;
    PDM_block_to_block_exch(btb,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            blk_stri_ini,
                            blk_data_ini,
              (int   *)     blk_stri_end,
              (void **)    &blk_data_end);

    int blk_stri_end_expected_p0[4] = {0, 2, 0, 3};
    int blk_stri_end_expected_p1[3] = {2, 0, 0};
    MPI_CHECK_EQ_C_ARRAY(0, blk_stri_end, blk_stri_end_expected_p0, 4);
    MPI_CHECK_EQ_C_ARRAY(1, blk_stri_end, blk_stri_end_expected_p1, 3);

    int blk_data_end_expected_p0[5] = {5, 8, 11, 0, 3};
    int blk_data_end_expected_p1[2] = {1, 9};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected_p0, 5);
    MPI_CHECK_EQ_C_ARRAY(1, blk_data_end, blk_data_end_expected_p1, 2);


    free(blk_data_end);

  }

  SUBCASE("Exchange stride variable with mpi_type ") {

    PDM_MPI_Datatype mpi_double_int_type;
    PDM_MPI_Type_create_contiguous(2, PDM_MPI_INT, &mpi_double_int_type);
    PDM_MPI_Type_commit(&mpi_double_int_type);

    int blk_stri_ini_p0[2] = {0, 2};
    int blk_stri_ini_p1[5] = {0, 3, 2, 0, 0};

    int blk_data_ini_p0[4] = {5, 50, 8, 80};
    int blk_data_ini_p1[10] = {11, 110, 0, 00, 3, 30, 1, 10, 9, 90};

    int *blk_data_ini = NULL;
    int *blk_stri_ini = NULL;
    if(test_rank == 0) {
      blk_stri_ini = blk_stri_ini_p0;
      blk_data_ini = blk_data_ini_p0;
    } else {
      blk_stri_ini = blk_stri_ini_p1;
      blk_data_ini = blk_data_ini_p1;
    }

    int  blk_stri_end[7] = {-1, -1, -1, -1, -1, -1, -1};
    int *blk_data_end = NULL;
    PDM_block_to_block_exch_with_mpi_type(btb,
                                          PDM_STRIDE_VAR_INTERLACED,
                                          mpi_double_int_type,
                                          blk_stri_ini,
                                          blk_data_ini,
                            (int   *)     blk_stri_end,
                            (void **)    &blk_data_end);

    int blk_stri_end_expected_p0[4] = {0, 2, 0, 3};
    int blk_stri_end_expected_p1[3] = {2, 0, 0};
    MPI_CHECK_EQ_C_ARRAY(0, blk_stri_end, blk_stri_end_expected_p0, 4);
    MPI_CHECK_EQ_C_ARRAY(1, blk_stri_end, blk_stri_end_expected_p1, 3);

    int blk_data_end_expected_p0[10] = {5, 50, 8, 80, 11, 110, 0, 00, 3, 30};
    int blk_data_end_expected_p1[4]  = {1, 10, 9, 90};
    MPI_CHECK_EQ_C_ARRAY(0, blk_data_end, blk_data_end_expected_p0, 2 * 5);
    MPI_CHECK_EQ_C_ARRAY(1, blk_data_end, blk_data_end_expected_p1, 2 * 2);

    free(blk_data_end);
    PDM_MPI_Type_free(&mpi_double_int_type);
  }

  PDM_block_to_block_free(btb);

}
