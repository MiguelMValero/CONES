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
    PDM_log_trace_array_long(distrib_init_elmt, n_rank+1, "distrib_init_elmt : ");
  }
  int n_part  = 1;
  int pn_elmt = freq * (distrib_init_elmt[i_rank+1] - distrib_init_elmt[i_rank]) ;

  PDM_g_num_t *pln_to_to_gn = malloc(pn_elmt * sizeof(PDM_g_num_t));
  int         *pfield       = malloc(pn_elmt * sizeof(int        ));
  int         *pstrid       = malloc(pn_elmt * sizeof(int        ));
  for(int i = 0; i < pn_elmt; ++i) {
    unsigned int seed = (unsigned int) (distrib_init_elmt[i_rank] + i);
    srand(seed);
    pln_to_to_gn[i] = (rand() % n_elmt) + 1;
    pfield[i]       = 1;
    pstrid      [i] = 1;
  }

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn, pn_elmt, "pln_to_to_gn : ");
  }

  /* Set child ln_to_gn from parent */

  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(1, n_part, PDM_FALSE, 1.e-3, comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gen_gnum, 0, pn_elmt, pln_to_to_gn);

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t *pln_to_to_gn_child = PDM_gnum_get(gen_gnum, 0);
  int pn_elmt_child = PDM_gnum_n_elt_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn_child, pn_elmt_child, "pln_to_to_gn_child : ");
  }

  /* Do part_to_block */

  // Create part_to_block

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &pln_to_to_gn_child,
                                                      NULL,
                                                     &pn_elmt_child,
                                                      n_part,
                                                      comm);

  // Exchange part_to_block

  int *dfield = NULL;
  int *dstrid = NULL;

  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                        &pstrid,
              (void **) &pfield,
                        &dstrid,
              (void **) &dfield);

  PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get(ptb);
  PDM_g_num_t *block_gnum = PDM_part_to_block_block_gnum_get(ptb);
  int nelmt_proc = PDM_part_to_block_n_elt_block_get(ptb);

  int dstrid_counter = 0;
  for (int i = 0; i < nelmt_proc; i++) {
    dstrid_counter += dstrid[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_long(block_gnum, nelmt_proc, "block_gnum : ");
    PDM_log_trace_array_int(dfield, dstrid_counter, "dfield : ");
  }

  /* Compute sum */

  int *dfield_summed = malloc(nelmt_proc * sizeof(int));
  int idx = 0;
  int partial_sum = 0;

  for (int i = 0; i < nelmt_proc; i++) {
    partial_sum = 0;
    for (int j = 0; j < dstrid[i]; j++) {
      partial_sum += dfield[idx++];
    }
    dfield_summed[i] = partial_sum;
  }

  if(0 == 1) {
    PDM_log_trace_array_int(dfield_summed, nelmt_proc, "dfield_summed : ");
  }

  /* Do block_to_part */

  // Create block_to_part

  PDM_block_to_part_t *btp = NULL;

  btp = PDM_block_to_part_create(distrib,
        (const PDM_g_num_t  **) &pln_to_to_gn_child,
                                &pn_elmt_child,
                                 n_part,
                                 comm);

  // Exchange in place block_to_part

  PDM_block_to_part_exch_in_place(btp,
                                  sizeof(int),
                                  PDM_STRIDE_CST_INTERLACED,
                                 &freq,
                         (void *) dfield_summed,
                                 &pstrid,
                       (void **) &pfield);

  /* Check summed values in partitions */
  if(1 == 0) {
    log_trace("Part vision :\n");
    for (int i = 0; i < pn_elmt_child; i++) {
      log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d\n", pln_to_to_gn_child[i], pfield[i]);
    }
  }

  /* Free entities */

  PDM_gnum_free(gen_gnum);
  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
  free(dfield);
  free(pstrid);
  free(dstrid);
  free(dfield_summed);

  free(pln_to_to_gn);
  free(distrib_init_elmt);
  free(pfield);

  PDM_MPI_Finalize ();
  return 0;
}
