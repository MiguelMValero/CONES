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

#include "pdm_block_to_part.h"

/**
 * \example pdm_t_block_to_part_sol1.c
 *
 * This is an example on how to do a block_to_part operation:
 */


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
  PDM_g_num_t *distrib_elmt = PDM_compute_uniform_entity_distribution(comm, n_elmt);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib_elmt, n_rank+1, "distrib_elmt : ");
  }

  int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];

  PDM_g_num_t* drand_number = malloc(dn_elmt * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_elmt; ++i) {
    unsigned int seed = (unsigned int) (distrib_elmt[i_rank] + i);
    srand(seed);
    drand_number[i] = (rand() % n_elmt) + 1;
  }

  if(0 == 1) {
    PDM_log_trace_array_long(drand_number, dn_elmt, "drand_number : ");
  }

  /*
   * Create partition :
   *   - Exch proc need to define a ln_to_gn that correspond to the desired partitioning you need
   *   - For this exemple you decide to take an arbitrary partitioning on 1 partition :
   *        Each proc take ln_to_gn = One over n_rank elemnts but with a selection strategy by freq
   */
  int n_part                    = 1;
  int          *pn_elmts        = malloc( n_part * sizeof(int           ));
  PDM_g_num_t **pelmts_ln_to_gn = malloc( n_part * sizeof(PDM_g_num_t * ));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_elmts       [i_part] = dn_elmt / freq;
    pelmts_ln_to_gn[i_part] = malloc( pn_elmts[i_part] * sizeof(PDM_g_num_t));
    PDM_g_num_t next_gnum = i_rank+1;
    for(int i = 0; i < pn_elmts[i_part]; ++i) {
      pelmts_ln_to_gn[i_part][i] = next_gnum;
      next_gnum += n_rank * freq;
      if(next_gnum > distrib_elmt[n_rank]) {
        next_gnum = next_gnum % distrib_elmt[n_rank];
      }
    }
  }

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pelmts_ln_to_gn[i_part], pn_elmts[i_part], "pelmts_ln_to_gn : ");
    }
  }

  /*
   * We have a mapping between part and block
   *   -> We want to get the associated part of the array drand_number
   */

  /*
   *  I/ Creation of block_to_part
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_elmt,
                          (const PDM_g_num_t **)      pelmts_ln_to_gn,
                                                      pn_elmts,
                                                      n_part,
                                                      comm);


  /*
   *  II/ Echange drand_number to get prand_number
   *        --> Tips : stride is constant
   */
  int stride_one = 1;
  PDM_g_num_t **prand_number = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         drand_number,
                         NULL,
              (void ***) &prand_number);

  /*
   * III/ Check results in partition
   */
  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(prand_number[i_part], pn_elmts[i_part], "prand_number : ");
    }
  }


  /*
   * IV/ Free exchange protocol
   */
  PDM_block_to_part_free(btp);

  free(drand_number   );
  free(distrib_elmt   );
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pelmts_ln_to_gn[i_part]);
    free(prand_number   [i_part]);
  }
  free(pelmts_ln_to_gn);
  free(prand_number);
  free(pn_elmts);

  PDM_MPI_Finalize ();
  return 0;
}
