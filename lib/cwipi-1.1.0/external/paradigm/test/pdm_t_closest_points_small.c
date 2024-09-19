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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_closest_points.h"
#include "pdm_version.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/
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

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);


  /* Define the numbers of Source/Target points */

  int n_src = 27;
  PDM_g_num_t src_gnum[27] = {1,21,17,18,19,20,22,15,23,24,25,26,16,14,2,7,3,4,5,6,8,13,9,10,11,12,27};
  double src_coords[3*27] = {
    0.166667, 0.166667, 0.166667,
    0.833333, 0.166667, 0.833333,
    0.500000, 0.833333, 0.500000,
    0.833333, 0.833333, 0.500000,
    0.166667, 0.166667, 0.833333,
    0.500000, 0.166667, 0.833333,
    0.166667, 0.500000, 0.833333,
    0.833333, 0.500000, 0.500000,
    0.500000, 0.500000, 0.833333,
    0.833333, 0.500000, 0.833333,
    0.166667, 0.833333, 0.833333,
    0.500000, 0.833333, 0.833333,
    0.166667, 0.833333, 0.500000,
    0.500000, 0.500000, 0.500000,
    0.500000, 0.166667, 0.166667,
    0.166667, 0.833333, 0.166667,
    0.833333, 0.166667, 0.166667,
    0.166667, 0.500000, 0.166667,
    0.500000, 0.500000, 0.166667,
    0.833333, 0.500000, 0.166667,
    0.500000, 0.833333, 0.166667,
    0.166667, 0.500000, 0.500000,
    0.833333, 0.833333, 0.166667,
    0.166667, 0.166667, 0.500000,
    0.500000, 0.166667, 0.500000,
    0.833333, 0.166667, 0.500000,
    0.833333, 0.833333, 0.833333
  };

  int n_tgt = 56;
  PDM_g_num_t tgt_gnum[56] = {36,41,40,39,38,37,35,34,33,32,31,42,43,29,50,55,54,53,52,51,49,44,48,47,46,45,30,28,1,7,12,11,10,9,8,6,14,5,4,3,2,13,27,21,26,25,24,23,22,20,15,19,18,17,16,56};
  double tgt_coords[3*56] = {
    1.225000, 1.125000, 0.625000,
    1.475000, 1.375000, 0.625000,
    1.225000, 1.375000, 0.625000,
    0.975000, 1.375000, 0.625000,
    1.725000, 1.125000, 0.625000,
    1.475000, 1.125000, 0.625000,
    0.975000, 1.125000, 0.625000,
    1.725000, 0.875000, 0.625000,
    1.475000, 0.875000, 0.625000,
    1.225000, 0.875000, 0.625000,
    1.725000, 0.625000, 0.625000,
    1.725000, 1.375000, 0.625000,
    1.225000, 0.625000, 0.875000,
    1.225000, 0.625000, 0.625000,
    1.225000, 1.125000, 0.875000,
    1.475000, 1.375000, 0.875000,
    1.225000, 1.375000, 0.875000,
    0.975000, 1.375000, 0.875000,
    1.725000, 1.125000, 0.875000,
    1.475000, 1.125000, 0.875000,
    0.975000, 1.125000, 0.875000,
    1.475000, 0.625000, 0.875000,
    1.725000, 0.875000, 0.875000,
    1.475000, 0.875000, 0.875000,
    1.225000, 0.875000, 0.875000,
    1.725000, 0.625000, 0.875000,
    1.475000, 0.625000, 0.625000,
    1.725000, 1.375000, 0.375000,
    1.225000, 0.625000, 0.125000,
    0.975000, 1.125000, 0.125000,
    1.225000, 1.375000, 0.125000,
    0.975000, 1.375000, 0.125000,
    1.725000, 1.125000, 0.125000,
    1.475000, 1.125000, 0.125000,
    1.225000, 1.125000, 0.125000,
    1.725000, 0.875000, 0.125000,
    1.725000, 1.375000, 0.125000,
    1.475000, 0.875000, 0.125000,
    1.225000, 0.875000, 0.125000,
    1.725000, 0.625000, 0.125000,
    1.475000, 0.625000, 0.125000,
    1.475000, 1.375000, 0.125000,
    1.475000, 1.375000, 0.375000,
    0.975000, 1.125000, 0.375000,
    1.225000, 1.375000, 0.375000,
    0.975000, 1.375000, 0.375000,
    1.725000, 1.125000, 0.375000,
    1.475000, 1.125000, 0.375000,
    1.225000, 1.125000, 0.375000,
    1.725000, 0.875000, 0.375000,
    1.225000, 0.625000, 0.375000,
    1.475000, 0.875000, 0.375000,
    1.225000, 0.875000, 0.375000,
    1.725000, 0.625000, 0.375000,
    1.475000, 0.625000, 0.375000,
    1.725000, 1.375000, 0.875000
  };


  //Brute force expected results
  PDM_g_num_t *expected_closest_gnum = (PDM_g_num_t *) malloc(n_tgt*sizeof(PDM_g_num_t));
  double      *expected_closest_dist = (double *) malloc(n_tgt*sizeof(double));
  for (int i = 0; i < n_tgt; i++) {
    PDM_g_num_t arg_min = -1;
    double min_dist = 1E12;
    for (int j = 0; j < n_src; j++) {
      double dist = (tgt_coords[3*i+0] - src_coords[3*j+0])*(tgt_coords[3*i+0] - src_coords[3*j+0])
                  + (tgt_coords[3*i+1] - src_coords[3*j+1])*(tgt_coords[3*i+1] - src_coords[3*j+1])
                  + (tgt_coords[3*i+2] - src_coords[3*j+2])*(tgt_coords[3*i+2] - src_coords[3*j+2]);
      if (dist < min_dist) {
        arg_min = src_gnum[j];
        min_dist = dist;
      }
    }
    expected_closest_gnum[i] = arg_min;
    expected_closest_dist[i] = min_dist;
  }


  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         1,
                                                         PDM_OWNERSHIP_USER);

  PDM_closest_points_n_part_cloud_set (clsp, 1, 1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coords,
                                    src_gnum);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coords,
                                    tgt_gnum);


  PDM_closest_points_compute (clsp);



  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  PDM_closest_points_free (clsp);

  /*
  printf("Results\n");
  for (int i = 0; i < n_tgt; i++) {
    printf("Tgt %d : closest is %d with dist %f -- expected : %d with dist %f\n", tgt_gnum[i], closest_src_gnum[i], closest_src_dist[i], expected_closest_gnum[i], expected_closest_dist[i]);
  }
  */

  for (int i = 0; i < n_tgt; i++) {
    assert (PDM_ABS(closest_src_dist[i] - expected_closest_dist[i]) < 1E-6);
    assert (closest_src_gnum[i] == expected_closest_gnum[i]);
  }

  
  free(expected_closest_gnum);
  free(expected_closest_dist);
  free(closest_src_gnum);
  free(closest_src_dist);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  PDM_MPI_Finalize ();

  return 0;
}



