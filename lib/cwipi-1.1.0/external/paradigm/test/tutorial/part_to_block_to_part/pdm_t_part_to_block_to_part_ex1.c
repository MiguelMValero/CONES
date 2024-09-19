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
#include "pdm_gnum.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_array.h"

#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */
static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Total number of elements .\n\n"
     "  -f      <level>  Frequency of extract .\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_elmt,
           int           *freq)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_elmt = atol(argv[i]);
        *n_elmt = (PDM_g_num_t) _n_elmt;
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        int _freq = atoi(argv[i]);
        *freq = (PDM_g_num_t) _freq;
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  PDM_g_num_t n_elmt = 10; // Number of elements in global configuration
  int         freq   = 1;
  _read_args(argc,
             argv,
             &n_elmt,
             &freq);
  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   *  Each proc have an equi repartition of data :
   *    - Define a distribution
   *    - Generate data among all pb by random
   */

  PDM_g_num_t *distrib_init_elmt = PDM_compute_uniform_entity_distribution(comm, n_elmt);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib_init_elmt, n_rank+1, "distrib_elmt : ");
  }
  int n_part  = 1;
  int pn_elmt = freq * (distrib_init_elmt[i_rank+1] - distrib_init_elmt[i_rank]) ;

  PDM_g_num_t *pln_to_to_gn = malloc(pn_elmt * sizeof(PDM_g_num_t));
  int         *pfield       = malloc(pn_elmt * sizeof(int        ));
  int         *pstrid       = malloc(pn_elmt * sizeof(int        )); // Added
  for(int i = 0; i < pn_elmt; ++i) {
    unsigned int seed = (unsigned int) (distrib_init_elmt[i_rank] + i);
    srand(seed);
    pln_to_to_gn[i] = (rand() % n_elmt) + 1;
    pfield[i]       = 1;
    pstrid      [i] = 1; // Added
  }

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn, pn_elmt, "pln_to_to_gn : ");
  }

  /*
   * I want to replace pfield values by their global, element-wise sum
   *
   * For example, if I have
   *   on rank 0 : pln_to_gn = [3, 6, 4], pfield = [1, 1, 1]
   *   on rank 1 : pln_to_gn = [4, 1],    pfield = [1, 1]
   *
   * In the end, I want to get
   *   on rank 0 : pfield = [1, 1, 2]
   *   on rank 1 : pfield = [2, 1]
   *
   * Tips : use part_to_block with PDM_PART_TO_BLOCK_POST_MERGE
   * (/!\ only available in PDM_STRIDE_VAR_INTERLACED)
   *
   */

  /* part_to_block to sum up in block vision using PDM_PART_TO_BLOCK_POST_MERGE */

  // Create part_to_block

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &pln_to_to_gn,
                                                      NULL,
                                                      &pn_elmt,
                                                      n_part,
                                                      comm);

  // Exchange part_to_block

  int *dfield = NULL;
  int *block_stride = NULL;

  int s_block_data = PDM_part_to_block_exch(ptb,
                                            sizeof(int),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            1,
                                            &pstrid,
                                  (void **) &pfield,
                                            &block_stride,
                                  (void **) &dfield);

  PDM_UNUSED(s_block_data);

  /* sum up what has been merged */

  // Compute size of block_stride DOES NOT WORK IF MORE THAN 1 PROC

  // int size_block_stride = 0;
  // int sum = 0;
  // while (sum != 10) {
  //   sum += block_stride[size_block_stride];
  //   size_block_stride++;
  // }

  int nelmt_proc = PDM_part_to_block_n_elt_block_get(ptb); // valeur de size_block_stride
  int size_block_stride = nelmt_proc;

  // Do the summation

  int *dfield_summed = malloc(size_block_stride * sizeof(int));

  int idx = 0;
  int tmp = 0;

  for (int i = 0; i < size_block_stride; i++) {

    tmp = 0;

    for (int j = 0; j < block_stride[i]; j++) {

      tmp += dfield[idx++];

    }

    dfield_summed[i] = tmp;

  }

  PDM_g_num_t *block_gnum = PDM_part_to_block_block_gnum_get(ptb);
  PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get(ptb); // bornes des gnum dans la partition

  PDM_UNUSED(distrib);

  // PDM_log_trace_array_int(dfield_summed, size_block_stride, "dfield_summed : ");
  // PDM_log_trace_array_long(block_gnum, nelmt_proc, "block_gnum : ");

  /* block_to_part to get back to the wanted format */

  // Create block_to_part

  PDM_block_to_part_t *btp = NULL;

  btp = PDM_block_to_part_create(distrib_init_elmt,
         (const PDM_g_num_t  **) &pln_to_to_gn,
                                 &pn_elmt,
                                 n_part,
                                 comm);

  // Create stride tab

  int n_gnum_rank = (distrib_init_elmt[i_rank+1] - distrib_init_elmt[i_rank]);

  int *dstrid = malloc(n_gnum_rank * sizeof(int));

  for (int i = 0; i < n_gnum_rank; i++) {
    dstrid[i] = 0;
  }

  for (int i = 0; i < nelmt_proc; i++) {
    dstrid[block_gnum[i] - distrib_init_elmt[i_rank] - 1] = 1;
  }

  // PDM_log_trace_array_int(dstrid, n_gnum_rank, "dstrid : ");
  // PDM_log_trace_array_int(pstrid, pn_elmt, "pstrid : ");

  // Exchange block_to_part

  PDM_block_to_part_exch_in_place(btp,
                                  sizeof(int),
                                  PDM_STRIDE_VAR_INTERLACED,
                                  dstrid,
                         (void *) dfield_summed,
                                  &pstrid,
                        (void **) &pfield);

  /* Check summed values in partitions */
  // log_trace("Part vision :\n");
  // for (int i = 0; i < pn_elmt; i++) {
  //   log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d\n", pln_to_to_gn[i], pfield[i]);
  // }

  /* Free the used entities */

  // Free part_to_block
  PDM_part_to_block_free(ptb);

  // Free block_to_part
  PDM_block_to_part_free(btp);

  // other free
  free(dfield);
  free(block_stride);
  free(dfield_summed);
  free(pstrid);
  free(dstrid);

  free(pln_to_to_gn);
  free(distrib_init_elmt);
  free(pfield);

  PDM_MPI_Finalize ();
  return 0;
}
