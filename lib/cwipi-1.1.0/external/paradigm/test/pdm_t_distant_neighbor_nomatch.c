#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_distant_neighbor.h"
#include "pdm_points_merge.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/*
 *  Use case with no match
 *
 *                   jn1             jn2
 * ------------------ |               | ------------------
 *         |          |               |            |
 *         |   (2)    |               |            |
 *         |          |               |     (1)    |
 * ------------------ |               |            |
 *         |          |               |            |
 *         |   (1)    |               | ------------------
 *         |          |               |            |
 * ------------------ |               |            |
 *         |          |               |     (0)    |
 *         |   (0)    |               |            |
 *         |          |               |            |
 * ------------------ |               | ------------------
 *
 */

/**
 *
 * \brief  Main
 *
 */

int
main
(
int   argc,
char *argv[]
)
{
  int verbose = 0;


  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  /* Connection de join1 avec join2 */
  int connect_idx_j1[4] = {0, 1, 3, 4};

  int oppRank = 0;
  int oppPart = 0;
  if(n_rank == 1){
    oppRank = 0;
    oppPart = 1;
  } else if (n_rank == 2){
    oppRank = 1;
    oppPart = 0;
  }

  int connect_triplet_j1[12] = {// Fisrt
                                oppRank, oppPart, 0,
                                // Second
                                oppRank, oppPart, 1,
                                oppRank, oppPart, 0,
                                // Third
                                oppRank, oppPart, 1};


  /* Connection de join2 avec join1 */
  // const int n_faces_j2 = 3;
  int connect_idx_j2[3] = {0, 2, 4};
  int connect_triplet_j2[12] = {// First
                                0, 0, 0,
                                0, 0, 1,
                                // Second
                                0, 0, 1,
                                0, 0, 2};
  // Cas avec croisement to check order in recv_data
  // int connect_triplet_j2[12] = {// First
  //                               0, 0, 1,
  //                               0, 0, 2,
  //                               // Second
  //                               0, 0, 0,
  //                               0, 0, 1};

  /*
   * Data
   */
  // Test with stride = 1
  int stride = 1;
  int point_list_j1[3] = {101, 201, 301};
  int point_list_j2[2] = {102, 202};

  // Test with stride = 2
  // int stride = 2;
  // int point_list_j1[6] = {101, 1010,
  //                         201, 2010,
  //                         301, 3010};
  // int point_list_j2[4] = {102, 1020,
  //                         202, 2020};

  // Test 1 :
  // int point_list_var_j1[6] = {101, 1010,
  //                             201, 2010,
  //                             301, 3010};
  // int point_list_var_stri_j1[3] = {2, 2, 2};

  // Test 2 :
  int point_list_var_j1[6] = {101,
                              201, 2010,
                              301, 3010};
  int point_list_var_stri_j1[3] = {1, 2, 2};

  // Test 3 :
  // int point_list_var_j1[6] = {101, 1010,
  //                             201,
  //                             301, 3010, 30100};
  // int point_list_var_stri_j1[3] = {2, 1, 3};

  // Test 1
  int point_list_var_j2[6] = {102, 1020, 10200,
                              202, 2020, 20200};
  int point_list_var_stri_j2[3] = {3, 3, 3};

  // Test 2
  // int point_list_var_j2[6] = {102,
  //                             202, 2020, 20200};
  // int point_list_var_stri_j2[2] = {1, 3};

  // Test 3
  // int point_list_var_j2[6] = {102, 1020, 10200,
  //                             202};
  // int point_list_var_stri_j2[2] = {3, 1};

  /*
   * Begin
   */

  int n_cloud = 0;
  int *n_entity = NULL;
  int **candidates_idx = NULL;
  int **candidates_desc = NULL;
  if(n_rank == 1){
    n_cloud = 2;
    candidates_idx  = (int **) malloc( n_cloud * sizeof(int*));
    candidates_desc = (int **) malloc( n_cloud * sizeof(int*));
    n_entity        = (int * ) malloc( n_cloud * sizeof(int ));

    n_entity[0] = 3;
    candidates_idx[0]  = connect_idx_j1;
    candidates_desc[0] = connect_triplet_j1;

    n_entity[1] = 2;
    candidates_idx[1]  = connect_idx_j2;
    candidates_desc[1] = connect_triplet_j2;
  } else if ( n_rank == 2){
    n_cloud = 1;
    candidates_idx  = (int **) malloc( n_cloud * sizeof(int*));
    candidates_desc = (int **) malloc( n_cloud * sizeof(int*));
    n_entity        = (int * ) malloc( n_cloud * sizeof(int ));
    if(i_rank == 0){
      n_entity[0]        = 3;
      candidates_idx[0]  = connect_idx_j1;
      candidates_desc[0] = connect_triplet_j1;
    } else if (i_rank == 1){
      n_entity[0]        = 2;
      candidates_idx[0]  = connect_idx_j2;
      candidates_desc[0] = connect_triplet_j2;
    }

  } else {

    PDM_error(__FILE__, __LINE__, 0, "pdm_t_distant_neighbor error : Bad number of process for test cases \n");
  }

  /*
   *  Setup exchange protocol
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(PDM_MPI_COMM_WORLD,
                                                           n_cloud,
                                                           n_entity,
                                                           candidates_idx,
                                                           candidates_desc);

  /*
   *  SetUp exchange
   */
  int** send_entity_data = (int **) malloc( n_cloud * sizeof(int* ));
  if(n_rank == 1){
    send_entity_data[0] = point_list_j1;
    send_entity_data[1] = point_list_j2;
  } else if(n_rank == 2){
    if(i_rank == 0){
      send_entity_data[0] = point_list_j1;
    } else if(i_rank == 1){
      send_entity_data[0] = point_list_j2;
    }
  }

  /*
   * Constant stride test
   */
  int** recv_entity_data = NULL;
  PDM_distant_neighbor_exch_int(dn,
                                sizeof(int),
                                PDM_STRIDE_CST_INTERLACED,
                                stride,
                                NULL,
                                send_entity_data,
                                NULL,
                                &recv_entity_data);

  if(verbose){
    log_trace(" Constant strid exchange results ---- \n");
    for(int i_part = 0; i_part < n_cloud; i_part++){
      int *_part_neighbor_idx  = candidates_idx[i_part];
      log_trace("recv_entity_data[%d]::", i_part);
      for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[i_part]]; i_entity++){
        for(int idata = 0; idata < stride; idata++){
          log_trace("%d ", recv_entity_data[i_part][stride*i_entity+idata]);
        }
      }
      log_trace("\n");
    }
  }

  /*
   * Variable stride test
   */
  int** send_entity_var_data = (int **) malloc( n_cloud * sizeof(int *));
  int** send_entity_var_stri = (int **) malloc( n_cloud * sizeof(int *));
  if(n_rank == 1){
    send_entity_var_data[0] = point_list_var_j1;
    send_entity_var_stri[0] = point_list_var_stri_j1;
    send_entity_var_data[1] = point_list_var_j2;
    send_entity_var_stri[1] = point_list_var_stri_j2;
  } else if(n_rank == 2){
    if(i_rank == 0){
    send_entity_var_data[0] = point_list_var_j1;
    send_entity_var_stri[0] = point_list_var_stri_j1;
    } else if(i_rank == 1){
    send_entity_var_data[0] = point_list_var_j2;
    send_entity_var_stri[0] = point_list_var_stri_j2;
    }
  }

  int** recv_entity_var_stri = NULL;
  int** recv_entity_var_data = NULL;
  // PDM_distant_neighbor_exch(dn,
  PDM_distant_neighbor_exch_int(dn,
                                sizeof(int),
                                PDM_STRIDE_VAR_INTERLACED,
                                -1,
                                send_entity_var_stri,
                                send_entity_var_data,
                       (int***) &recv_entity_var_stri,
                                &recv_entity_var_data);

  if(verbose){
    log_trace(" Variable strid exchange results ---- \n");
    for(int i_part = 0; i_part < n_cloud; i_part++){
      int *_part_neighbor_idx  = candidates_idx[i_part];
      log_trace(" ---> recv_entity_data[%d]::", i_part);
      int idx = 0;
      for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[i_part]]; i_entity++){
        log_trace("i_entity::%d - strid::%d --> ", i_entity, recv_entity_var_stri[i_part][i_entity]);
        for(int i_data = 0; i_data < recv_entity_var_stri[i_part][i_entity]; i_data++){
          log_trace("%d ", recv_entity_var_data[i_part][idx++]);
        }
        log_trace("\n");
      }
      log_trace("\n");
    }
  }
  /*
   * Free
   */
  PDM_distant_neighbor_free(dn);
  free(candidates_idx);
  free(candidates_desc);
  free(n_entity);
  free(send_entity_data);
  free(send_entity_var_data);
  free(send_entity_var_stri);
  for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
    free(recv_entity_data[i_cloud]);
    free(recv_entity_var_stri[i_cloud]);
    free(recv_entity_var_data[i_cloud]);
  }
  free(recv_entity_data);
  free(recv_entity_var_stri);
  free(recv_entity_var_data);
  PDM_MPI_Finalize();

  if (i_rank == 0) {
    PDM_printf ("End\n");
  }

  return 0;

}
