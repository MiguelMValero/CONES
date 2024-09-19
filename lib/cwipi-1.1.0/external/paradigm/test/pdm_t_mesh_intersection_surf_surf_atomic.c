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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_array.h"

#include "pdm_mesh_intersection_surf_surf_atomic.h"

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private functions
 *============================================================================*/

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}

static void
_read_args
(
 int    argc,
 char **argv,
 char **filename
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



static void _read_data
(
 char   *filename,
 double  triaA_coord[9],
 double  edgeB_coord[6],
 double  edgeB_normal[6]
 )
{
  FILE *f = fopen(filename, "r");

  if (f == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
  }

  for (int i = 0; i < 9; i++) {
    fscanf(f, "%lf", &triaA_coord[i]);
  }

  for (int i = 0; i < 6; i++) {
    fscanf(f, "%lf", &edgeB_coord[i]);
  }

  for (int i = 0; i < 6; i++) {
    fscanf(f, "%lf", &edgeB_normal[i]);
  }

  fclose(f);
}

/*============================================================================
 * Main
 *============================================================================*/

int main(int argc, char *argv[])
{
  // Init
  PDM_MPI_Init(&argc, &argv);

  char *filename = NULL;
  _read_args(argc,
             argv,
             &filename);

  double triaA_coord[9] = {
    0., 0., 0.,
    1., 0., 0.,
    0., 1., 0.
  };
  double edgeB_coord[6] = {
    -0.18, 1.30, 0.,
    0.88, 0.73, 0.
  };
  double edgeB_normal[6] = {
    0., 0.1, 1.,
    0., -0.1, 1.
  };


  if (filename != NULL) {
    _read_data(filename,
               triaA_coord,
               edgeB_coord,
               edgeB_normal);
  }


  for (int i = 0; i < 2; i++) {
    double mag = PDM_MODULE(&edgeB_normal[3*i]);
PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
    if (mag == 0) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid normal (null vector)\n");
    }
PDM_GCC_SUPPRESS_WARNING_POP
    for (int j = 0; j < 3; j++) {
      edgeB_normal[3*i+j] /= mag;
    }
  }


  // vtk?

  double area = PDM_mesh_intersection_surf_surf_atomic_compute(triaA_coord,
                                                               edgeB_coord,
                                                               edgeB_normal);

  printf("area = %20.16f\n", area);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}

