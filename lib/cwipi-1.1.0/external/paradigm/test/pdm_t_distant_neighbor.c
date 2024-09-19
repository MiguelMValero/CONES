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
#include "pdm_error.h"
#include "pdm_logging.h"

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{
  int verbose = 0;

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  /*
   * The test case represent two conform interface between 2 partitions, we need to find out connection of each faces
   *    I/  Localisation in the frame of the current/opposite boundary condition
   *    II/ Use distant_neighbor to find out the face number in the partition
   */
  int point_list_j1[15] = {145, 153, 155, 163, 172, 180, 201, 206, 207, 211, 212, 215, 217, 222, 227};
  int point_list_var_j1[45] = {145, 1450, 14500,
                               153, 1530, 15300,
                               155, 1550, 15500,
                               163, 1630, 16300,
                               172, 1720, 17200,
                               180, 1800, 18000,
                               201, 2010, 20100,
                               206, 2060, 20600,
                               207, 2070, 20700,
                               211, 2110, 21100,
                               212, 2120, 21200,
                               215, 2150, 21500,
                               217, 2170, 21700,
                               222, 2220, 22200,
                               227, 2270, 22700};

  int point_list_var_stri_j1[15] = {3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3,
                                    3};

  double xyz_j1[45] = {0.5,  2.5,  4. ,  0.5,  1.5,  4. ,  1.5,  2.5,  4. ,  1.5,  1.5,
                       4. ,  2.5,  1.5,  4. ,  3.5,  1.5,  4. ,  1.5,  0.5,  4. ,  2.5,
                       0.5,  4. ,  2.5,  2.5,  4. ,  3.5,  0.5,  4. ,  3.5,  2.5,  4. ,
                       4.5,  0.5,  4. ,  4.5,  1.5,  4. ,  0.5,  0.5,  4. ,  4.5,  2.5,
                       4. };

  double cln_j1[15] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1., 1.};

  int point_list_j2[15] = {136, 142, 150, 151, 160, 161, 170, 171, 188, 191, 194, 199, 205, 214, 221};

  int point_list_var_j2[30] = {136, 1360,
                               142, 1420,
                               150, 1500,
                               151, 1510,
                               160, 1600,
                               161, 1610,
                               170, 1700,
                               171, 1710,
                               188, 1880,
                               191, 1910,
                               194, 1940,
                               199, 1990,
                               205, 2050,
                               214, 2140,
                               221, 2210};

  int point_list_var_stri_j2[15] = {2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2,
                                    2};


  double xyz_j2[45] = {0.5,  1.5,  4. ,  1.5,  1.5,  4. ,  1.5,  2.5,  4. ,  2.5,  1.5,
                       4. ,  2.5,  2.5,  4. ,  3.5,  1.5,  4. ,  3.5,  2.5,  4. ,  4.5,
                       1.5,  4. ,  0.5,  0.5,  4. ,  1.5,  0.5,  4. ,  2.5,  0.5,  4. ,
                       3.5,  0.5,  4. ,  4.5,  0.5,  4. ,  4.5,  2.5,  4. ,  0.5,  2.5,
                       4.};

  double cln_j2[15] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1., 1.};

  /*
   * Define i_cloud
   */
  int      n_cloud  = -1;
  int      n_points = 15;
  double** coords = NULL;
  double** char_lenght = NULL;
  if(n_rank == 1){
    n_cloud = 2;
    coords      = (double **) malloc( n_cloud * sizeof(double *));
    char_lenght = (double **) malloc( n_cloud * sizeof(double *));

    coords     [0] = xyz_j1;
    coords     [1] = xyz_j2;
    char_lenght[0] = cln_j1;
    char_lenght[1] = cln_j2;

  } else if ( n_rank == 2){
    n_cloud = 1;
    coords      = (double **) malloc( n_cloud * sizeof(double*));
    char_lenght = (double **) malloc( n_cloud * sizeof(double*));
    if(i_rank == 0){
      coords     [0] = xyz_j1;
      char_lenght[0] = cln_j1;
    } else {//if (i_rank == 1){
      coords     [0] = xyz_j2;
      char_lenght[0] = cln_j2;
    }

  } else {
    PDM_error(__FILE__, __LINE__, 0, "pdm_t_distant_neighbor error : Bad number of process for test cases \n");
  }


  /*
   * Call pdm_points_merge to find match between the two block
   */
  double tolerance = 0.1;
  PDM_points_merge_t* pts_merge = PDM_points_merge_create(n_cloud, tolerance, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);

  if(n_rank == 1){
    PDM_points_merge_cloud_set(pts_merge, 0, n_points, coords[0], char_lenght[0]);
    PDM_points_merge_cloud_set(pts_merge, 1, n_points, coords[1], char_lenght[1]);
  } else if(n_rank == 2){
    if(i_rank == 0){
      PDM_points_merge_cloud_set(pts_merge, 0, n_points, coords[0], char_lenght[0]);
    } else if(i_rank == 1){
      PDM_points_merge_cloud_set(pts_merge, 0, n_points, coords[0], char_lenght[0]);
    }
  }

  PDM_points_merge_process(pts_merge);

  /*
   * Get resulting points_merge
   */
  int  *n_entity        = (int * ) malloc( n_cloud * sizeof(int ));
  int **candidates_idx  = (int **) malloc( n_cloud * sizeof(int*));
  int **candidates_desc = (int **) malloc( n_cloud * sizeof(int*));
  for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
    int n_cloud_points = -1;
    int n_desc   = -1;
    PDM_points_merge_candidates_get(pts_merge, i_cloud, &candidates_idx[i_cloud], &candidates_desc[i_cloud]);
    PDM_points_merge_candidates_size_get(pts_merge, i_cloud, &n_cloud_points, &n_desc);

    n_entity[i_cloud] = n_cloud_points;

    assert(n_desc == candidates_idx[i_cloud][n_cloud_points]);
    if(verbose){
      printf("-- n_desc:: %d \n ", n_desc);
      printf("-- candidates_idx[i_cloud][n_cloud_points+1]:: %d \n ", candidates_idx[i_cloud][n_cloud_points]);
      for(int i = 0; i < n_cloud_points; i++){
        printf("-- %d %d ", i_cloud, i);
        for (int j = candidates_idx[i_cloud][i];
                 j <  candidates_idx[i_cloud][i+1]; j++) {
          printf(" : %d", candidates_desc[i_cloud][3*j]);
          printf(" %d", candidates_desc[i_cloud][3*j+1]);
          printf(" %d", candidates_desc[i_cloud][3*j+2]);
        }
        printf("\n");
      }
    }
  }

  /*
   *  Now we have the connection between the two cloud (in this case it's a perfect match )
   *  We need to find out the connection in the partition (not in the cloud )
   *  To do that we setup an exchange with distant neighbor with the local face of each partition
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(PDM_MPI_COMM_WORLD,
                                                           n_cloud,
                                                           n_entity,
                                                           candidates_idx,
                                                           candidates_desc);

  /*
   *  SetUp exchange
   */
  int** send_entity_data = (int **) malloc( n_cloud * sizeof(int*));
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
  int stride = 1;
  int** recv_entity_data = NULL;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            stride,
                            NULL,
                 ( void**)  send_entity_data,
                            NULL,
                 (void***) &recv_entity_data);

  /*
   * Variable stride test
   */
  int** send_entity_var_data = (int **) malloc( n_cloud * sizeof(int*));
  int** send_entity_var_stri = (int **) malloc( n_cloud * sizeof(int*));
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
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            send_entity_var_stri,
                 (void**)   send_entity_var_data,
                 (int***) &recv_entity_var_stri,
                (void***) &recv_entity_var_data);

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
  PDM_points_merge_free(pts_merge);
  free(coords);
  free(char_lenght);
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

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize();

  return 0;
}


