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
  PDM_UNUSED(n_part);
  int pn_elmt = freq * (distrib_init_elmt[i_rank+1] - distrib_init_elmt[i_rank]) ;

  PDM_g_num_t *pln_to_to_gn = malloc(pn_elmt * sizeof(PDM_g_num_t));
  int         *pfield       = malloc(pn_elmt * sizeof(int        ));
  for(int i = 0; i < pn_elmt; ++i) {
    unsigned int seed = (unsigned int) (distrib_init_elmt[i_rank] + i);
    srand(seed);
    pln_to_to_gn[i] = (rand() % n_elmt) + 1;
    pfield[i]       = 1;
  }

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn, pn_elmt, "pln_to_to_gn : ");
  }

  /*
   * I want to replace pfield values by their global element-wise sum
   *
   * Tips : use part_to_block with PDM_PART_TO_BLOCK_POST_MERGE
   * (/!\ only available in PDM_STRIDE_VAR_INTERLACED)
   *
   */

  /*
   *  1) Move to block-vision
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &pln_to_to_gn,
                                                      NULL,
                                                      &pn_elmt,
                                                      n_part,
                                                      comm);

  int *part_stride = PDM_array_const_int(pn_elmt, 1); // Create a stride array full of 1s
  int *block_stride = NULL;
  int *block_data   = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &pfield,
                         &block_stride,
               (void **) &block_data);
  free(part_stride);



  /*
   *  2) Compute sum in block-vision
   */
  int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb);
  int *summed_block_data = malloc(sizeof(int) * n_elt_block);
  int idx = 0;
  for (int i = 0; i < n_elt_block; i++) {
    summed_block_data[i] = 0.;
    for (int j = 0; j < block_stride[i]; j++) {
      summed_block_data[i] += block_data[idx];
      idx++;
    }
  }

  PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get(ptb);
  if(1 == 0) {
    log_trace("Block vision :\n");
    for (int i = 0; i < n_elt_block; i++) {
      log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d\n", block_g_num[i], summed_block_data[i]);
    }
  }


  /*
   *  3) Send summed values back to initial partitions
   */
  /* First method (easier) */
  // log_trace("\n\n~~~ Method 1 ~~~\n");
  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_g_num,
                                                                        n_elt_block,
                                                 (const PDM_g_num_t **) &pln_to_to_gn,
                                                                        &pn_elmt,
                                                                        n_part,
                                                                        comm);

  int one = 1;
  PDM_block_to_part_exch_in_place(btp,
                                  sizeof(int),
                                  PDM_STRIDE_CST_INTERLACED,
                                  &one,
                         (void *) summed_block_data,
                                  NULL,
                        (void **) &pfield);

  PDM_block_to_part_free(btp);


  /* Check summed values in partitions */
  if(1 == 0) {
    log_trace("Part vision :\n");
    for (int i = 0; i < pn_elmt; i++) {
      log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d\n", pln_to_to_gn[i], pfield[i]);
    }
  }


  /* Alternative method (more involved) */
  // log_trace("\n\n~~~ Method 2 ~~~\n");
  PDM_array_reset_int(block_stride, n_elt_block, 1);
  PDM_g_num_t *distrib_full = PDM_part_to_block_adapt_partial_block_to_block(ptb,
                                                                             &block_stride,
                                                                             distrib_init_elmt[n_rank]);
  if(1 == 0) {
    PDM_log_trace_array_long(distrib_full, n_rank + 1, "distrib_full : ");
    PDM_log_trace_array_int(block_stride, (int) (distrib_full[i_rank+1] - distrib_full[i_rank]), "block_stride : ");
    PDM_log_trace_array_int(summed_block_data, n_elt_block, "summed_block_data : ");
  }

  btp = PDM_block_to_part_create(distrib_full,
          (const PDM_g_num_t **) &pln_to_to_gn,
                                 &pn_elmt,
                                 n_part,
                                 comm);

  int **part_stride2 = NULL;
  int **part_data2   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_stride,
                (void *) summed_block_data,
                         &part_stride2,
              (void ***) &part_data2);
  PDM_block_to_part_free(btp);
  free(distrib_full);



  /* Check summed values in partitions */
  if(1 == 0) {
    log_trace("Part vision :\n");
    for (int i = 0; i < pn_elmt; i++) {
      log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d\n", pln_to_to_gn[i], part_data2[0][i]);
    }
  }

  PDM_part_to_block_free(ptb);

  free(pln_to_to_gn);
  free(distrib_init_elmt);
  free(pfield);

  free(block_data);
  free(block_stride);
  free(summed_block_data);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(part_stride2[i_part]);
    free(part_data2[i_part]);
  }
  free(part_stride2);
  free(part_data2);

  PDM_MPI_Finalize ();
  return 0;
}
