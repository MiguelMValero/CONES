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
#include "pdm_gnum_from_hash_values.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*
 *  Use case
 *
 *       3           4     4           2
 *       +-----------+     +-----------+
 *       |           |     |           |
 *       |           |     |           |
 *       |           |     |           |
 *       +-----------+     +-----------+
 *       5           1     1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */


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

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);



  int edge_vtx_p0[8] = {1, 4, /* Edge 1 */
                        4, 3,
                        3, 5,
                        5, 1};

  int edge_str_p0[4] = {2, 2, 2, 2};

  int edge_vtx_p1[8] = {2, 4, /* Edge 1 */
                        4, 1,
                        1, 6,
                        6, 2};
  int edge_str_p1[4] = {2, 2, 2, 2};

  PDM_g_num_t expected_ln_to_gn_p0[4] = {1, 5, 7, 2};
  PDM_g_num_t expected_ln_to_gn_p1[4] = {3, 1, 4, 6};

  /*
   * SetUp
   */
  int     n_part = 0;
  int    *n_elmts = NULL;
  int    **part_stri = NULL;
  int    **part_data = NULL;
  size_t **part_key = NULL;
  if(n_rank == 1){
    n_part = 2;
    n_elmts   = (int *     ) malloc( n_part * sizeof(int    ));
    part_stri = (int **    ) malloc( n_part * sizeof(int*   ));
    part_data = (int **    ) malloc( n_part * sizeof(int*   ));
    part_key  = (size_t ** ) malloc( n_part * sizeof(size_t*));

    n_elmts[0] = 4;
    n_elmts[1] = 4;

    part_data[0] = edge_vtx_p0;
    part_data[1] = edge_vtx_p1;

    part_stri[0] = edge_str_p0;
    part_stri[1] = edge_str_p1;

  } else if( n_rank == 2) {
    n_part = 1;
    n_elmts   = (int *     ) malloc( n_part * sizeof(int    ));
    part_stri = (int **    ) malloc( n_part * sizeof(int*   ));
    part_data = (int **    ) malloc( n_part * sizeof(int*   ));
    part_key  = (size_t ** ) malloc( n_part * sizeof(size_t*));

    if( i_rank == 0) {
      n_elmts[0] = 4;
      part_data[0] = edge_vtx_p0;
      part_stri[0] = edge_str_p0;

    } else if( i_rank == 1) {

      n_elmts[0] = 4;
      part_data[0] = edge_vtx_p1;
      part_stri[0] = edge_str_p1;
    }
  }

  /*
   * Compute key
   */
  for(int i_part = 0; i_part < n_part; i_part++){
    part_key[i_part] = (size_t *) malloc(n_elmts[i_part] * sizeof(size_t));
    int idx = 0;
    for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
      size_t key = 0;
      for(int idata = 0; idata < part_stri[i_part][ielmt]; ++idata){
        key += part_data[i_part][idx++];
      }
      // printf(" part_key[%d][%d] = %lu \n", i_part, ielmt, key);
      part_key[i_part][ielmt] = key;
    }
  }


  PDM_bool_t equilibrate = PDM_FALSE;


  PDM_gnum_from_hv_t *gnum_fhv_id = PDM_gnum_from_hash_values_create(n_part,
                                                                     equilibrate,
                                                                     sizeof(int),
                                                                     PDM_operator_compare_connectivity,
                                                                     PDM_operator_equal_connectivity,
                                                                     PDM_MPI_COMM_WORLD,
                                                                     PDM_OWNERSHIP_KEEP);

  for(int i_part = 0; i_part < n_part; ++i_part){
    PDM_gnum_set_hash_values(gnum_fhv_id,
                             i_part,
                             n_elmts[i_part],
                             part_key[i_part],
                             part_stri[i_part],
            (unsigned char*) part_data[i_part]);
  }

  PDM_gnum_from_hv_compute(gnum_fhv_id); /* Passage de part --> block */

  /*
   *
   */
  PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part){
    ln_to_gn[i_part] = PDM_gnum_from_hv_get(gnum_fhv_id, i_part);
    if (verbose) {
      PDM_log_trace_array_long(ln_to_gn[i_part], n_elmts[i_part], "ln_to_gn::");
    }
  }

  // Check results : independant of the parallelisme
  if(n_rank == 1){
    for(int ielmt = 0; ielmt < n_elmts[0]; ++ielmt) {
      assert(ln_to_gn[0][ielmt] == expected_ln_to_gn_p0[ielmt]);
    }
    for(int ielmt = 0; ielmt < n_elmts[1]; ++ielmt) {
      assert(ln_to_gn[1][ielmt] == expected_ln_to_gn_p1[ielmt]);
    }
  } else {
    if( i_rank == 0) {
      for(int ielmt = 0; ielmt < n_elmts[0]; ++ielmt) {
        assert(ln_to_gn[0][ielmt] == expected_ln_to_gn_p0[ielmt]);
      }
    } else if (i_rank == 1) {
      for(int ielmt = 0; ielmt < n_elmts[0]; ++ielmt) {
        assert(ln_to_gn[0][ielmt] == expected_ln_to_gn_p1[ielmt]);
      }
    }
  }





  /*
   * Free
   */
  PDM_gnum_from_hv_free(gnum_fhv_id);
  free(part_stri);
  free(part_data);
  free(n_elmts);
  for(int i_part = 0; i_part < n_part; ++i_part){
    free(part_key[i_part]);
  }
  free(part_key);
  free(ln_to_gn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;

}
