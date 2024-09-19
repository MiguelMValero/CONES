#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_linear_algebra.h"



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
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 int                   *n_row,
 int                   *n_col,
 int                   *stride,
 int                   *seed,
 double                *tol
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n_row") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_row = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-n_col") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_col = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-stride") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *stride = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *seed = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-tol") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tol = atof(argv[i]);
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static double *
_rand_array
(
 const int siz
)
{
  double *a = malloc(sizeof(double) * siz);

  for (int i = 0; i < siz; i++) {
    a[i] = (double) rand() / (double) RAND_MAX;
    a[i] = 2*a[i] - 1;
  }

  return a;
}


int main
(
 int   argc,
 char *argv[]
 )
{
  int    n_row  = 3;
  int    n_col  = 3;
  int    stride = 2;
  int    seed   = -1;
  double tol    = 0.;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_row,
             &n_col,
             &stride,
             &seed,
             &tol);

  if (seed < 0) {
    srand(time(NULL));
  }
  else {
    srand(seed);
  }


  PDM_MPI_Init(&argc, &argv);


  double *a = _rand_array(n_row * n_col);
  double *b = _rand_array(n_row * stride);

  double *_a = malloc(sizeof(double) * n_row * n_col);
  memcpy(_a, a, sizeof(double) * n_row * n_col);


  double *x = malloc(sizeof(double) * n_col * stride);


  int stat = PDM_linear_algebra_linsolve_svd(n_row,
                                             n_col,
                                             stride,
                                             tol,
                                             a,
                                             b,
                                             x);

  printf("stat = %d\n", stat);


  /* Check Ax = b */
  printf("b - Ax = \n");
  for (int i = 0; i < n_row; i++) {
    for (int k = 0; k < stride; k++) {
      for (int j = 0; j < n_col; j++) {
        b[stride*i+k] -= _a[n_col*i+j] * x[stride*j+k];
      }

      printf("%e ", b[stride*i+k]);
    }
    printf("\n");
  }

  free(a);
  free(b);
  free(x);
  free(_a);

  PDM_MPI_Finalize();

  return 0;
}
