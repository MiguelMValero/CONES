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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_part_to_part.h"
#include "pdm_array.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(n_rank == 2);

  int n_part1 = 3;
  int n_part2 = 3;

  int *n_elt1 = malloc(sizeof(int) * n_part1);
  PDM_g_num_t **gnum_elt1 = malloc(sizeof(PDM_g_num_t *) * n_part1);
  int **part1_to_part2_idx = malloc(sizeof(int *) * n_part1);
  PDM_g_num_t **part1_to_part2 = malloc(sizeof(PDM_g_num_t *) * n_part1);

  int *n_elt2 = n_elt1;
  PDM_g_num_t **gnum_elt2 = gnum_elt1;


  if (i_rank == 0) {
    n_elt1[0] = 1;
    n_elt1[1] = 4;
    n_elt1[2] = 2;
  } else {
    n_elt1[0] = 2;
    n_elt1[1] = 6;
    n_elt1[2] = 2;
  }

  for (int i = 0; i < n_part1; i++) {
    gnum_elt1[i]          = malloc(sizeof(PDM_g_num_t) * n_elt1[i]);
    part1_to_part2_idx[i] = malloc(sizeof(int)         * (n_elt1[i] + 1));
  }


  if (i_rank == 0) {
    gnum_elt1[0][0] = 9;

    gnum_elt1[1][0] = 15;
    gnum_elt1[1][1] = 13;
    gnum_elt1[1][2] = 5;
    gnum_elt1[1][3] = 2;

    gnum_elt1[2][0] = 3;
    gnum_elt1[2][1] = 8;

    part1_to_part2_idx[0][0] = 0;
    part1_to_part2_idx[0][1] = 0;

    part1_to_part2_idx[1][0] = 0;
    part1_to_part2_idx[1][1] = 2;
    part1_to_part2_idx[1][2] = 4;
    part1_to_part2_idx[1][3] = 5;
    part1_to_part2_idx[1][4] = 6;

    part1_to_part2_idx[2][0] = 0;
    part1_to_part2_idx[2][1] = 1;
    part1_to_part2_idx[2][2] = 1;
  }

  else {
    gnum_elt1[0][0] = 12;
    gnum_elt1[0][1] = 16;

    gnum_elt1[1][0] = 11;
    gnum_elt1[1][1] = 10;
    gnum_elt1[1][2] = 17;
    gnum_elt1[1][3] = 14;
    gnum_elt1[1][4] = 1;
    gnum_elt1[1][5] = 4;

    gnum_elt1[2][0] = 7;
    gnum_elt1[2][1] = 6;

    part1_to_part2_idx[0][0] = 0;
    part1_to_part2_idx[0][1] = 0;
    part1_to_part2_idx[0][2] = 0;

    part1_to_part2_idx[1][0] = 0;
    part1_to_part2_idx[1][1] = 0;
    part1_to_part2_idx[1][2] = 0;
    part1_to_part2_idx[1][3] = 2;
    part1_to_part2_idx[1][4] = 4;
    part1_to_part2_idx[1][5] = 6;
    part1_to_part2_idx[1][6] = 8;

    part1_to_part2_idx[2][0] = 0;
    part1_to_part2_idx[2][1] = 1;
    part1_to_part2_idx[2][2] = 2;
  }

  for (int i = 0; i < n_part1; i++) {
    part1_to_part2[i] = malloc(sizeof(PDM_g_num_t) * part1_to_part2_idx[i][n_elt1[i]]);
  }

  if (i_rank == 0) {
    // part1_to_part2[0][0] = ;

    part1_to_part2[1][0] = 9;
    part1_to_part2[1][1] = 11;
    part1_to_part2[1][2] = 9;
    part1_to_part2[1][3] = 10;
    part1_to_part2[1][4] = 11;
    part1_to_part2[1][5] = 10;

    part1_to_part2[2][0] = 9;
  }

  else {
    // part1_to_part2[0][0] = ;

    part1_to_part2[1][0] = 12;
    part1_to_part2[1][1] = 11;
    part1_to_part2[1][2] = 12;
    part1_to_part2[1][3] = 10;
    part1_to_part2[1][4] = 16;
    part1_to_part2[1][5] = 11;
    part1_to_part2[1][6] = 16;
    part1_to_part2[1][7] = 10;

    part1_to_part2[2][0] = 12;
    part1_to_part2[2][1] = 16;
  }

  /*
   *  Create Part-to-part object
   */
  PDM_part_to_part_t *ptp = PDM_part_to_part_create ((const PDM_g_num_t **) gnum_elt1,
                                                     n_elt1,
                                                     n_part1,
                                                     (const PDM_g_num_t **)gnum_elt2,
                                                     n_elt2,
                                                     n_part2,
                                                     (const int **) part1_to_part2_idx,
                                                     (const PDM_g_num_t **) part1_to_part2,
                                                     comm);
  /*
   *  Check Part2 ref/unref/gnum1_come_from
   */
  int  *n_ref_num2 = NULL;
  int **ref_num2   = NULL;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_ref_num2,
                                  &ref_num2);

  int  *n_unref_num2 = NULL;
  int **unref_num2   = NULL;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_unref_num2,
                                  &unref_num2);

  int **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from = NULL;
  PDM_part_to_part_gnum1_come_from_get (ptp,
                                        &gnum1_come_from_idx,
                                        &gnum1_come_from);
  // for (int i = 0; i < n_part2; i++) {

  //   log_trace("\npart2 %d\n", i);
  //   PDM_log_trace_array_int(ref_num2[i], n_ref_num2[i], "referenced (l_num) : ");
  //   log_trace("unreferenced (g_num) : ");
  //   for (int j = 0; j < n_unref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" ", gnum_elt2[i][unref_num2[i][j]-1]);
  //   }
  //   log_trace("\n");
  //   log_trace("referenced (g_num) -> gnum1_come_from :\n");
  //   for (int j = 0; j < n_ref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" -> ", gnum_elt2[i][ref_num2[i][j]-1]);
  //     for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
  //       log_trace(PDM_FMT_G_NUM" ", gnum1_come_from[i][k]);
  //     }
  //     log_trace("\n");
  //   }


  // }


  int         **part1_stride = malloc(sizeof(int         *) * n_part1);
  PDM_g_num_t **part1_data   = malloc(sizeof(PDM_g_num_t *) * n_part1);

  for (int i = 0; i < n_part1; i++) {

    // log_trace("\npart1 %d\n", i);

    part1_stride[i] = malloc(sizeof(int) * n_elt1[i]);

    int s_part1_data = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      part1_stride[i][j] = (int) gnum_elt1[i][j];
      // part1_stride[i][j] = 1;
      s_part1_data += part1_stride[i][j];
    }
    // PDM_log_trace_array_int(part1_stride[i], n_elt1[i], "part1_stride : ");

    // log_trace("g_num -> data:\n");
    part1_data[i] = malloc(sizeof(PDM_g_num_t) * s_part1_data);
    int idx = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      // int idx0 = idx;
      for (int k = 1; k <= part1_stride[i][j]; k++) {
        part1_data[i][idx++] = k;
      }
      // log_trace(PDM_FMT_G_NUM, gnum_elt1[i][j]);
      // PDM_log_trace_array_long(part1_data[i] + idx0,
      //                          part1_stride[i][j],
      //                          " -> ");
    }
    // for (int j = 0; j < n_elt1[i]; j++) {
    //   part1_data[i][j] = gnum_elt1[i][j];
    // }
  }


  int         **part2_stride = NULL;
  PDM_g_num_t **part2_data   = NULL;
  int request;

  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_VAR_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                          0,
                          sizeof(PDM_g_num_t),
           (const int **) part1_stride,
          (const void **) part1_data,
                          &part2_stride,
               (void ***) &part2_data,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);

  // log_trace("\n\n---- Check iexch ----\n");
  // for (int i = 0; i < n_part2; i++) {

  //   log_trace("\npart2 %d\n", i);
  //   PDM_log_trace_array_int(part2_stride[i], gnum1_come_from_idx[i][n_ref_num2[i]], "stride : ");
  //   log_trace("referenced (g_num) -> data :\n");
  //   int idx = 0;
  //   for (int j = 0; j < n_ref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" :\n", gnum_elt2[i][ref_num2[i][j]-1]);
  //     for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
  //       log_trace("    "PDM_FMT_G_NUM" -> ", gnum1_come_from[i][k]);
  //       for (int l = 0; l < part2_stride[i][k]; l++) {
  //         log_trace(PDM_FMT_G_NUM" ", part2_data[i][idx++]);
  //       }
  //       log_trace("\n");
  //     }
  //     log_trace("\n");
  //   }


  // }




  /*
   *  Exchange an interleaved, constant-stride field
   */
  // log_trace("\n\n---- Exchange an interleaved, constant-stride field ----\n");
  PDM_g_num_t **part1_field = malloc(sizeof(PDM_g_num_t *) * n_part1);
  for (int i = 0; i < n_part1; i++) {
    // int n = part1_to_part2_idx[i][n_elt1[i]];
    // part1_field[i] = malloc(sizeof(PDM_g_num_t) * n * 2);

    // for (int j = 0; j < n_elt1[i]; j++) {
    //   for (int k = part1_to_part2_idx[i][j]; k < part1_to_part2_idx[i][j+1]; k++) {
    //     part1_field[i][k  ] = gnum_elt1[i][j];
    //     part1_field[i][k+n] = part1_to_part2[i][k];
    //   }
    // }
    int n = n_elt1[i];
    part1_field[i] = malloc(sizeof(PDM_g_num_t) * n * 2);

    for (int j = 0; j < n_elt1[i]; j++) {
      part1_field[i][j  ] = gnum_elt1[i][j];
      part1_field[i][j+n] = gnum_elt1[i][j]+1;
    }

    // log_trace("\npart1 %d\n", i);
    // PDM_log_trace_array_long(part1_field[i],     n, "  part1_field (1st comp) : ");
    // PDM_log_trace_array_long(part1_field[i] + n, n, "  part1_field (2nd comp) : ");
  }



  PDM_g_num_t **part2_field = NULL;
  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_CST_INTERLEAVED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,//PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                          2,
                          sizeof(PDM_g_num_t),
                          NULL,
          (const void **) part1_field,
                          NULL,
               (void ***) &part2_field,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);


  // for (int i = 0; i < n_part2; i++) {
  //   log_trace("\npart2 %d\n", i);
  //   int n = gnum1_come_from_idx[i][n_ref_num2[i]];
  //   PDM_log_trace_array_long(part2_field[i],     n, "  part2_field (1st comp) : ");
  //   PDM_log_trace_array_long(part2_field[i] + n, n, "  part2_field (2nd comp) : ");
  // }



  /* Reverse */
  for (int i = 0; i < n_part1; i++) {
    free (part1_field[i]);
  }
  free (part1_field);

  PDM_part_to_part_reverse_iexch (ptp,
                                  PDM_MPI_COMM_KIND_P2P,
                                  PDM_STRIDE_CST_INTERLEAVED,
                                  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                  2,
                                  sizeof(PDM_g_num_t),
                                  NULL,
                  (const void **) part2_field,
                                  NULL,
                       (void ***) &part1_field,
                                  &request);

  PDM_part_to_part_reverse_iexch_wait (ptp,
                                       request);

  // log_trace("Reverse\n");
  // for (int i = 0; i < n_part1; i++) {
  //   int n = part1_to_part2_idx[i][n_elt1[i]];
  //   log_trace("\npart1 %d\n", i);
  //   PDM_log_trace_array_long(part1_field[i],     n, "  part1_field (1st comp) : ");
  //   PDM_log_trace_array_long(part1_field[i] + n, n, "  part1_field (2nd comp) : ");
  // }


  // log_trace("==== P1 -> P2 ====\n");
  /* 2 consecutive iexch in stride var with same stride */
  for (int ipart = 0; ipart < n_part1; ipart++) {
    int s_part1_data = 0;
    part1_stride[ipart] = realloc(part1_stride[ipart], sizeof(int) * n_elt1[ipart]);
    for (int i = 0; i < n_elt1[ipart]; i++) {
      part1_stride[ipart][i] = (int) (gnum_elt1[ipart][i] % 2) + 1;
      s_part1_data += part1_stride[ipart][i];
    }

    part1_data[ipart] = realloc(part1_data[ipart], sizeof(PDM_g_num_t) * s_part1_data);
    int idx = 0;
    for (int i = 0; i < n_elt1[ipart]; i++) {
      for (int j = 0; j < part1_stride[ipart][i]; j++) {
        part1_data[ipart][idx++] = gnum_elt1[ipart][i];
      }
    }
  }

  for (int i = 0; i < n_part2; i++) {
    free (part2_stride[i]);
    free (part2_data  [i]);
  }
  free (part2_stride);
  free (part2_data);


  part2_stride = NULL;
  // log_trace("1\n");
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
         (const int  **) part1_stride,
         (const void **) part1_data,
                         &part2_stride,
              (void ***) &part2_data,
                         &request);
  PDM_part_to_part_iexch_wait (ptp, request);

  // log_trace("2\n");
  PDM_g_num_t **part2_data2 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
         (const int  **) part1_stride,
         (const void **) part1_data,
                         &part2_stride,
              (void ***) &part2_data2,
                         &request);
  PDM_part_to_part_iexch_wait (ptp, request);

  for (int i = 0; i < n_part2; i++) {
    int idx = 0;
    for (int j = 0; j < n_ref_num2[i]; j++) {
      // int id2 = ref_num2[i][j] - 1;
      for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
        // log_trace("gnum2 "PDM_FMT_G_NUM", gnum1 "PDM_FMT_G_NUM", expected stride = %d, got %d\n",
        //           gnum_elt2[i][id2], gnum1_come_from[i][k], (int) (gnum1_come_from[i][k]%2) + 1,
        //           part2_stride[i][k]);
        for (int l = 0; l < part2_stride[i][k]; l++) {
          // log_trace("  "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
          //           part2_data2[i][idx], part2_data[i][idx]);
          assert(part2_data2[i][idx] == part2_data[i][idx]);
          idx++;
        }
      }
    }

    free (part2_data2[i]);
  }
  free (part2_data2);



  // log_trace("==== P1 <- P2 ====\n");
  /* 2 consecutive reverse iexch in stride var with same stride */
  for (int ipart = 0; ipart < n_part2; ipart++) {
    int s_part2_data = 0;
    part2_stride[ipart] = realloc(part2_stride[ipart], sizeof(int) * n_elt2[ipart]);
    for (int i = 0; i < n_elt2[ipart]; i++) {
      part2_stride[ipart][i] = (int) (gnum_elt2[ipart][i] % 2) + 1;
      s_part2_data += part2_stride[ipart][i];
    }

    part2_data[ipart] = realloc(part2_data[ipart], sizeof(PDM_g_num_t) * s_part2_data);
    int idx = 0;
    for (int i = 0; i < n_elt2[ipart]; i++) {
      for (int j = 0; j < part2_stride[ipart][i]; j++) {
        part2_data[ipart][idx++] = gnum_elt2[ipart][i];
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    free (part1_stride[i]);
    free (part1_data  [i]);
  }
  free (part1_stride);
  free (part1_data);

  part1_stride = NULL;
  // log_trace("1\n");
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                 (const int  **) part2_stride,
                 (const void **) part2_data,
                                 &part1_stride,
                      (void ***) &part1_data,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait (ptp, request);
  // for (int i = 0; i < n_part1; i++) {
  //   log_trace("part1 %d\n", i);
  //   PDM_log_trace_array_int(part1_stride[i], part1_to_part2_idx[i][n_elt1[i]], "part1_stride : ");
  // }
  // log_trace("2\n");

  // for (int i = 0; i < n_part1; i++) {
  //   free(part1_stride[i]);
  // }
  // free(part1_stride);
  // part1_stride = NULL;

  PDM_g_num_t **part1_data2 = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                 (const int  **) part2_stride,
                 (const void **) part2_data,
                                 &part1_stride,
                      (void ***) &part1_data2,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait (ptp, request);
  // for (int i = 0; i < n_part1; i++) {
  //   log_trace("part1 %d\n", i);
  //   PDM_log_trace_array_int(part1_stride[i], part1_to_part2_idx[i][n_elt1[i]], "part1_stride : ");
  // }

  for (int i = 0; i < n_part1; i++) {
    int idx = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      for (int k = part1_to_part2_idx[i][j]; k < part1_to_part2_idx[i][j+1]; k++) {
        // log_trace("gnum1 "PDM_FMT_G_NUM", gnum2 "PDM_FMT_G_NUM", expected stride = %d, got %d\n",
        //           gnum_elt1[i][j], part1_to_part2[i][k], (int) (part1_to_part2[i][k]%2) + 1,
        //           part1_stride[i][k]);
        for (int l = 0; l < part1_stride[i][k]; l++) {
          // log_trace("  "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
          //           part1_data2[i][idx], part1_data[i][idx]);
          assert(part1_data2[i][idx] == part1_data[i][idx]);
          idx++;
        }
      }
    }

    free (part1_data2[i]);
  }
  free (part1_data2);


  /*
   *  Free memory
   */
  PDM_part_to_part_free (ptp);

  for (int i = 0; i < n_part1; i++) {
    free (gnum_elt1[i]);
    free (part1_to_part2_idx[i]);
    free (part1_to_part2[i]);

    free (part1_stride[i]);
    free (part1_data[i]);
    free (part1_field[i]);
  }

  for (int i = 0; i < n_part2; i++) {
    free (part2_stride[i]);
    free (part2_data[i]);
    free (part2_field[i]);
  }

  free (n_elt1);
  free (gnum_elt1);
  free (part1_to_part2_idx);
  free (part1_to_part2);

  free (part1_stride);
  free (part1_data);
  free (part1_field);
  free (part2_stride);
  free (part2_data);
  free (part2_field);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
