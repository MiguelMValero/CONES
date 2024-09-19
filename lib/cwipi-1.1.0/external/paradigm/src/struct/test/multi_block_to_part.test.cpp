#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_multi_block_to_part.h"

MPI_TEST_CASE("[pdm_multi_block_to_part] - 1p - Simple ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static int n_block = 2;
  static PDM_g_num_t multi_distrib_idx[3] = {0, 6, 9};
  static PDM_g_num_t block_distrib1_idx[2] = {0, 6};
  static PDM_g_num_t block_distrib2_idx[2] = {0, 3};

  PDM_g_num_t** block_distrib_idx = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
  block_distrib_idx[0] = block_distrib1_idx;
  block_distrib_idx[1] = block_distrib2_idx;

  static const int n_part         = 1;
  static const int n_elmts_part_0 = 6;
  PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {9, 7, 5, 2, 1, 3};

  // Convenient array
  int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  ln_to_gn[0] = ln_to_gn_p0;
  n_elmts[0]  = n_elmts_part_0;

  PDM_multi_block_to_part_t* mbtp =
    PDM_multi_block_to_part_create(multi_distrib_idx,
                                   n_block,
          (const PDM_g_num_t**)    block_distrib_idx,
          (const PDM_g_num_t**)    ln_to_gn,
                                   n_elmts,
                                   n_part,
                                   pdm_comm);

  SUBCASE("constant stride = 1") {
    // Stride cst with 2 blk_data reprensenting a implicite block distribution of 6+3 elemnts
    static const int darray_size1 = 6;
    static const int darray_size2 = 3;
    static PDM_g_num_t darray1[darray_size1] = {5 , 4 , 6, 1, 3, 2};
    static PDM_g_num_t darray2[darray_size2] = {10, 11, 8};
    PDM_g_num_t** darray = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
    darray[0] = darray1;
    darray[1] = darray2;

    /* Exchange */
    int** stride_one = (int ** ) malloc( n_block * sizeof(int *));
    for(int i_block = 0; i_block < n_block; ++i_block){
      stride_one[i_block] = (int * ) malloc( 1 * sizeof(int));
      stride_one[i_block][0] = 1;
    }

    PDM_g_num_t** parray = NULL;
    PDM_multi_block_to_part_exch2(mbtp, sizeof(PDM_g_num_t), PDM_STRIDE_CST_INTERLACED,
                                  stride_one,
                       (void ** ) darray,
                                  NULL,
                       (void ***) &parray);

    static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {8, 10, 3, 4, 5, 6};

    // for(int i_part = 0; i_part < n_part; ++i_part) {
    //   for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
    //     printf(" parray[%i][%i] = %i \n", i_part, ielmt, parray[i_part][ielmt]);
    //   }
    // }
    CHECK_EQ_C_ARRAY(parray[0], parray_expected_p0, n_elmts_part_0);

    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
    }
    free(parray);
    for(int i_block = 0; i_block < n_block; ++i_block){
      free(stride_one[i_block]);
    }
    free(stride_one);
    free(darray);
  }


  SUBCASE("variable stride ") {

    // Stride cst with 2 blk_data reprensenting a implicite block distribution of 6+3 elemnts
    static const int dnelmt1 = 6;
    static const int dnelmt2 = 3;
    static const int darray_size1 = 11;
    static const int darray_size2 = 6;
    static PDM_g_num_t darray1[darray_size1] = {5 , 50, 4 , 6, 60, 600, 1, 3, 2, 20, 200};
    static int         dstrid1[dnelmt1]      = {2     , 1 , 3         , 1, 1, 3         };
    static PDM_g_num_t darray2[darray_size2] = {10, 100, 11, 8, 80, 800};
    static int         dstrid2[dnelmt2]      = {2      , 1 , 3         };
    PDM_g_num_t** darray = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
    darray[0] = darray1;
    darray[1] = darray2;
    int** dstri  = (int **) malloc( n_block * sizeof(int *));
    dstri[0] = dstrid1;
    dstri[1] = dstrid2;

    int**         pstrid = NULL;
    PDM_g_num_t** parray = NULL;
    PDM_multi_block_to_part_exch2(mbtp, sizeof(PDM_g_num_t), PDM_STRIDE_VAR_INTERLACED,
                                  dstri,
                       (void ** ) darray,
                                  &pstrid,
                       (void ***) &parray);

    // for(int i_part = 0; i_part < n_part; ++i_part) {
    //   int idx_data = 0;
    //   for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
    //     printf(" pstrid[%i][%i] = %i \n", i_part, ielmt, pstrid[i_part][ielmt]);
    //     for(int i_data = 0; i_data < pstrid[i_part][ielmt]; ++i_data) {
    //       printf(" parray[%i][%i] = %i \n", i_part, ielmt, parray[i_part][idx_data++]);
    //     }
    //   }
    // }

   static PDM_g_num_t parray_expected_p0[12]             = {8, 80, 800, 10, 100, 3, 4, 5, 50, 6, 60, 600};
   static int         pstrid_expected_p0[n_elmts_part_0] = {3         , 2      , 1, 1, 2    , 3         };

    CHECK_EQ_C_ARRAY(pstrid[0], pstrid_expected_p0, n_elmts_part_0);
    CHECK_EQ_C_ARRAY(parray[0], parray_expected_p0, 12);

    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
      free(pstrid[i_part]);
    }
    free(pstrid);
    free(parray);
    free(dstri);
    free(darray);

  }

  // Free
  PDM_multi_block_to_part_free(mbtp);
  free(block_distrib_idx);
  free(ln_to_gn);
  free(n_elmts);

}


MPI_TEST_CASE("[pdm_multi_block_to_part] - 2p - Simple ",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static int n_block = 2;
  static PDM_g_num_t multi_distrib_idx[3] = {0, 6, 9};
  static PDM_g_num_t block_distrib1_idx[3] = {0, 3, 6};
  static PDM_g_num_t block_distrib2_idx[3] = {0, 2, 3};

  PDM_g_num_t** block_distrib_idx = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
  block_distrib_idx[0] = block_distrib1_idx;
  block_distrib_idx[1] = block_distrib2_idx;

  static const int n_part         = 1;
  static const int n_elmts_part_0 = 3;
  static const int n_elmts_part_1 = 3;
  PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {9, 7, 5};
  PDM_g_num_t ln_to_gn_p1[n_elmts_part_1] = {2, 1, 3};

  // Convenient array
  int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  if(test_rank == 0){
    ln_to_gn[0] = ln_to_gn_p0;
    n_elmts[0]  = n_elmts_part_0;
  } else {
    ln_to_gn[0] = ln_to_gn_p1;
    n_elmts[0]  = n_elmts_part_1;
  }

  PDM_multi_block_to_part_t* mbtp =
    PDM_multi_block_to_part_create(multi_distrib_idx,
                                   n_block,
          (const PDM_g_num_t**)    block_distrib_idx,
          (const PDM_g_num_t**)    ln_to_gn,
                                   n_elmts,
                                   n_part,
                                   pdm_comm);

  SUBCASE("constant stride = 1") {
    // Stride cst with 2 blk_data reprensenting a implicite block distribution of 6+3 elemnts
    static const int darray_size1_p0 = 3;
    static const int darray_size2_p0 = 2;
    static PDM_g_num_t darray1_p0[darray_size1_p0] = {5 , 4 , 6};
    static PDM_g_num_t darray2_p0[darray_size2_p0] = {10, 11};

    static const int darray_size1_p1 = 3;
    static const int darray_size2_p1 = 1;
    static PDM_g_num_t darray1_p1[darray_size1_p1] = {1, 3, 2};
    static PDM_g_num_t darray2_p1[darray_size2_p1] = {8};
    PDM_g_num_t** darray = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));

    if(test_rank == 0){
      darray[0] = darray1_p0;
      darray[1] = darray2_p0;
    } else {
      darray[0] = darray1_p1;
      darray[1] = darray2_p1;
    }

    /* Exchange */
    int** stride_one = (int ** ) malloc( n_block * sizeof(int *));
    for(int i_block = 0; i_block < n_block; ++i_block){
      stride_one[i_block] = (int * ) malloc( 1 * sizeof(int));
      stride_one[i_block][0] = 1;
    }

    PDM_g_num_t** parray = NULL;
    PDM_multi_block_to_part_exch2(mbtp, sizeof(PDM_g_num_t), PDM_STRIDE_CST_INTERLACED,
                                  stride_one,
                       (void ** ) darray,
                                  NULL,
                       (void ***) &parray);

    static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {8, 10, 3};
    static PDM_g_num_t parray_expected_p1[n_elmts_part_1] = {4, 5, 6};

    // for(int i_part = 0; i_part < n_part; ++i_part) {
    //   for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
    //     printf(" parray[%i][%i] = %i \n", i_part, ielmt, parray[i_part][ielmt]);
    //   }
    // }
    MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected_p0, n_elmts_part_0);
    MPI_CHECK_EQ_C_ARRAY(1, parray[0], parray_expected_p1, n_elmts_part_1);

    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
    }
    free(parray);
    for(int i_block = 0; i_block < n_block; ++i_block){
      free(stride_one[i_block]);
    }
    free(stride_one);
    free(darray);
  }

  SUBCASE("variable stride ") {
    // Stride cst with 2 blk_data reprensenting a implicite block distribution of 6+3 elemnts
    static const int dnelmt1_p0 = 3;
    static const int dnelmt2_p0 = 2;

    static const int darray_size1_p0 = 6;
    static const int darray_size2_p0 = 3;

    static PDM_g_num_t darray1_p0[darray_size1_p0] = {5 , 50, 4 , 6, 60, 600};
    static int         dstrid1_p0[dnelmt1_p0]      = {2     , 1 , 3         };
    static PDM_g_num_t darray2_p0[darray_size2_p0] = {10, 100, 11};
    static int         dstrid2_p0[dnelmt2_p0]      = {2      , 1 };

    static const int dnelmt1_p1 = 3;
    static const int dnelmt2_p1 = 1;

    static const int darray_size1_p1 = 5;
    static const int darray_size2_p1 = 3;

    static PDM_g_num_t darray1_p1[darray_size1_p1] = {1, 3, 2, 20, 200};
    static int         dstrid1_p1[dnelmt1_p1]      = {1, 1, 3         };
    static PDM_g_num_t darray2_p1[darray_size2_p1] = {8, 80, 800};
    static int         dstrid2_p1[dnelmt2_p1]      = {3         };


    PDM_g_num_t** darray = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));

    if(test_rank == 0){
      darray[0] = darray1_p0;
      darray[1] = darray2_p0;
    } else {
      darray[0] = darray1_p1;
      darray[1] = darray2_p1;
    }

    int** dstri  = (int **) malloc( n_block * sizeof(int *));

    if(test_rank == 0){
      dstri[0] = dstrid1_p0;
      dstri[1] = dstrid2_p0;
    } else {
      dstri[0] = dstrid1_p1;
      dstri[1] = dstrid2_p1;
    }

    int**         pstrid = NULL;
    PDM_g_num_t** parray = NULL;
    PDM_multi_block_to_part_exch2(mbtp, sizeof(PDM_g_num_t), PDM_STRIDE_VAR_INTERLACED,
                                  dstri,
                       (void ** ) darray,
                                  &pstrid,
                       (void ***) &parray);

    // for(int i_part = 0; i_part < n_part; ++i_part) {
    //   int idx_data = 0;
    //   for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
    //     printf(" pstrid[%i][%i] = %i \n", i_part, ielmt, pstrid[i_part][ielmt]);
    //     for(int i_data = 0; i_data < pstrid[i_part][ielmt]; ++i_data) {
    //       printf(" parray[%i][%i] = %i \n", i_part, ielmt, parray[i_part][idx_data++]);
    //     }
    //   }
    // }

    static PDM_g_num_t parray_expected_p0[6]              = {8, 80, 800, 10, 100, 3};
    static int         pstrid_expected_p0[n_elmts_part_0] = {3         , 2      , 1};

    static PDM_g_num_t parray_expected_p1[6]              = {4, 5, 50, 6, 60, 600};
    static int         pstrid_expected_p1[n_elmts_part_1] = {1, 2    , 3         };

    MPI_CHECK_EQ_C_ARRAY(0, pstrid[0], pstrid_expected_p0, n_elmts_part_0);
    MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected_p0, 6);

    MPI_CHECK_EQ_C_ARRAY(1, pstrid[0], pstrid_expected_p1, n_elmts_part_1);
    MPI_CHECK_EQ_C_ARRAY(1, parray[0], parray_expected_p1, 6);

    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
      free(pstrid[i_part]);
    }
    free(pstrid);
    free(parray);
    free(dstri);
    free(darray);

  }

  // Free
  PDM_multi_block_to_part_free(mbtp);
  free(block_distrib_idx);
  free(ln_to_gn);
  free(n_elmts);

}

#include <vector>
#include <algorithm>

MPI_TEST_CASE("[pdm_multi_block_to_part] - 3p - n_block=1",3) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_block = 1;
  // int n_rank = 3;

  int n_elt = 7;
  std::vector<PDM_g_num_t> multi_distrib_idx = {0,21};
  std::vector<std::vector<PDM_g_num_t>> block_distribs_storer = {{0,7,14,21}};

  std::vector<PDM_g_num_t*> block_distribs(n_block);
  for (int i=0; i<n_block; ++i) {
    block_distribs[i] = block_distribs_storer[i].data();
  }

  const int n_part = 1;

  std::vector<PDM_g_num_t> ln_to_gn_0(n_elt);
  std::iota(begin(ln_to_gn_0),end(ln_to_gn_0),test_rank*n_elt+1);
  std::vector<PDM_g_num_t*> ln_to_gn = {ln_to_gn_0.data()};
  std::vector<int> n_elts = {n_elt};


  PDM_multi_block_to_part_t* mbtp =
    PDM_multi_block_to_part_create(multi_distrib_idx.data(),
                                   n_block,
             (const PDM_g_num_t**) block_distribs.data(),
             (const PDM_g_num_t**) ln_to_gn.data(),
                                   n_elts.data(),
                                   n_part,
                                   pdm_comm);
  int* stride_one = (int*)malloc( 1 * sizeof(int));
  stride_one[0] = 1;

  std::vector<int32_t> darray(n_elt);
  if (test_rank==0) darray = { 0, 1, 2, 3, 4, 5, 6};
  if (test_rank==1) darray = { 7, 8, 9,10,11,12,13};
  if (test_rank==2) darray = {14,15,16,17,18,19,20};
  std::vector<int32_t*> darray_ptr(n_block);
  darray_ptr[0] = darray.data();

  int32_t** parray = nullptr;
  PDM_multi_block_to_part_exch2(mbtp, sizeof(int32_t), PDM_STRIDE_CST_INTERLACED,
                                &stride_one,
                     (void ** ) darray_ptr.data(),
                                nullptr,
                     (void ***) &parray);

  // since there is only one block, the exchange is just the identity
  MPI_CHECK_EQ_C_ARRAY(0, parray[0], darray.data(), n_elt);
  MPI_CHECK_EQ_C_ARRAY(1, parray[0], darray.data(), n_elt);
  MPI_CHECK_EQ_C_ARRAY(2, parray[0], darray.data(), n_elt);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(parray[i_part]);
  }
  free(parray);
  free(stride_one);
  PDM_multi_block_to_part_free(mbtp);
}


MPI_TEST_CASE("[pdm_multi_block_to_part] - 3p - n_block=2",3) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_block = 2;
  // int n_rank = 3;
  int i_rank = test_rank;

  int n_elt = 7;
  std::vector<PDM_g_num_t> multi_distrib_idx = {0,21,42};
  std::vector<std::vector<PDM_g_num_t>> block_distribs_storer = {{0,7,14,21},{0,7,14,21}};

  std::vector<PDM_g_num_t*> block_distribs(n_block);
  for (int i=0; i<n_block; ++i) {
    block_distribs[i] = block_distribs_storer[i].data();
  }

  const int n_part = 1;

  std::vector<PDM_g_num_t> merged_distri = {0,14,28,42};

  int n_elts_0 = merged_distri[i_rank+1]-merged_distri[i_rank];
  std::vector<PDM_g_num_t> ln_to_gn_0(n_elts_0);
  std::iota(begin(ln_to_gn_0),end(ln_to_gn_0),merged_distri[i_rank]+1);
  std::vector<PDM_g_num_t*> ln_to_gn = {ln_to_gn_0.data()};
  std::vector<int> n_elts = {n_elts_0};

  PDM_multi_block_to_part_t* mbtp =
    PDM_multi_block_to_part_create(multi_distrib_idx.data(),
                                   n_block,
             (const PDM_g_num_t**) block_distribs.data(),
             (const PDM_g_num_t**) ln_to_gn.data(),
                                   n_elts.data(),
                                   n_part,
                                   pdm_comm);
  int stride = 1;
  int** stride_one = (int ** ) malloc( n_block * sizeof(int *));
  for(int i_block = 0; i_block < n_block; ++i_block){
    stride_one[i_block] = (int * ) malloc( 1 * sizeof(int));
    stride_one[i_block][0] = stride;
  }

  std::vector<int32_t> darray0;
  std::vector<int32_t> darray1;
  if (i_rank==0) darray0 = { 0, 1, 2, 3, 4, 5, 6};
  if (i_rank==1) darray0 = { 7, 8, 9,10,11,12,13};
  if (i_rank==2) darray0 = {14,15,16,17,18,19,20};
  if (i_rank==0) darray1 = {  0, 10, 20, 30, 40, 50, 60};
  if (i_rank==1) darray1 = { 70, 80, 90,100,110,120,130};
  if (i_rank==2) darray1 = {140,150,160,170,180,190,200};
  std::vector<int32_t*> darray_ptr(n_block);
  darray_ptr[0] = darray0.data();
  darray_ptr[1] = darray1.data();

  int32_t** parray = nullptr;
  PDM_multi_block_to_part_exch2(mbtp, sizeof(int32_t), PDM_STRIDE_CST_INTERLACED,
                                stride_one,
                     (void ** ) darray_ptr.data(),
                                nullptr,
                     (void ***) &parray);

  // since there is only one block, the exchange is just the identity
  std::vector<int32_t> parray_expected;
  if (i_rank==0) parray_expected = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13};
  if (i_rank==1) parray_expected = { 14, 15, 16, 17, 18, 19, 20,  0, 10, 20, 30, 40, 50, 60};
  if (i_rank==2) parray_expected = { 70, 80, 90,100,110,120,130,140,150,160,170,180,190,200};
  MPI_CHECK_EQ_C_ARRAY(0, parray[0], parray_expected.data(), n_elt*2);
  MPI_CHECK_EQ_C_ARRAY(1, parray[0], parray_expected.data(), n_elt*2);
  MPI_CHECK_EQ_C_ARRAY(2, parray[0], parray_expected.data(), n_elt*2);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(parray[i_part]);
  }
  free(parray);
  for(int i_block = 0; i_block < n_block; ++i_block){
    free(stride_one[i_block]);
  }
  free(stride_one);
  PDM_multi_block_to_part_free(mbtp);
}
