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

#include "pdm_gnum.h"


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
  }

  if(0 == 1) {
    PDM_log_trace_array_long(pln_to_to_gn, pn_elmt, "pln_to_to_gn : ");
  }

  /*
   * Extract pair gnum
   */
  int pn_elmt_pair = 0;
  for(int i = 0; i < pn_elmt; ++i) {
    if(pln_to_to_gn[i] % 2 == 0) {
      pn_elmt_pair++;
    }
  }
  PDM_g_num_t *ppair_ln_to_gn = malloc(pn_elmt_pair * sizeof(PDM_g_num_t));
  pn_elmt_pair = 0;
  for(int i = 0; i < pn_elmt; ++i) {
    if(pln_to_to_gn[i] % 2 == 0) {
      ppair_ln_to_gn[pn_elmt_pair++] = pln_to_to_gn[i];
    }
  }

  /*
   * I want to have a global numbering for all pair value
   * Tips : use pdm_gnum after extraction of pair value
   *
   */
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_KEEP);
  PDM_gnum_set_from_parents(gen_gnum, 0 /* i_part*/, pn_elmt_pair, ppair_ln_to_gn);
  PDM_gnum_compute         (gen_gnum);

  /*
   * Print the block_g_num and distrib
   */
  PDM_g_num_t *ppair_child_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

  /*
   *  Exchange field and print it / Check !
   */
  if(0 == 1) {
    PDM_log_trace_array_long(ppair_ln_to_gn      , pn_elmt_pair, "ppair_ln_to_gn       : ");
    PDM_log_trace_array_long(ppair_child_ln_to_gn, pn_elmt_pair, "ppair_child_ln_to_gn : ");
  }


  PDM_gnum_free (gen_gnum);

  free(pln_to_to_gn);
  free(ppair_ln_to_gn);
  free(distrib_init_elmt);
  free(pfield);


  PDM_MPI_Finalize ();
  return 0;
}
