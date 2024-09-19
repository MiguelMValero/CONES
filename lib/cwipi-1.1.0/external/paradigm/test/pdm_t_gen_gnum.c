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
#include "pdm_printf.h"
#include "pdm_gnum.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_point_cloud_gen.h"


/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -v               Verbose mode\n",
     "  -n_part <n_part> Number of partitions\n",
     "  -n      <n>      Global number of points per partition\n",
     "  -m               Merge quasi-coincident points\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc     Number of arguments
 * \param [in]    argv     Arguments
 * \param [inout] verbose  Verbose mode
 * \param [inout] n_part   Number of partitions
 * \param [inout] gn_pts   Global number of points per partition
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *verbose,
 int           *n_part,
 PDM_g_num_t   *gn_pts,
 PDM_bool_t    *merge
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }

    else if (strcmp (argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }

    else if (strcmp (argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) n;
      }
    }

    else if (strcmp (argv[i], "-merge") == 0) {
      *merge = (PDM_bool_t) 1;
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
  PDM_MPI_Init(&argc, &argv);

  int         verbose = 0;
  int         n_part  = 1;
  PDM_g_num_t gn_elts = 10;
  PDM_bool_t  merge   = PDM_FALSE;
  _read_args(argc,
             argv,
             &verbose,
             &n_part,
             &gn_elts,
             &merge);

  /*
   * Generate a random point cloud
   */
  int n_rank;
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  double *char_length = malloc(sizeof(double) * gn_elts);
  for (int i = 0; i < gn_elts; i++) {
    char_length[i] = 1e-3;
  }

  int     *n_elts = malloc(sizeof(int     ) * n_part);
  double **coords = malloc(sizeof(double *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_g_num_t *gnum = NULL;
    PDM_point_cloud_gen_random(PDM_MPI_COMM_WORLD,
                               n_rank*i_part + i_rank, // seed
                               0,
                               gn_elts,
                               0., 0., 0.,
                               1., 1., 1.,
                               &n_elts[i_part],
                               &coords[i_part],
                               &gnum);
    free(gnum);
  }

  /*
   * Generate a global numbering for the points from their coordinates
   */

  // First, create a PDM_gen_gnum_t instance and set some parameters
  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,     // dimension
                                             n_part,
                                             merge,
                                             1.e-3, // tolerance
                                             PDM_MPI_COMM_WORLD,
                                             PDM_OWNERSHIP_USER);

  // Then, provide the coordinates array for each partition
  // (`char_length` can be NULL if `merge` is disabled)
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_coords(gen_gnum,
                             i_part,
                             n_elts[i_part],
                             coords[i_part],
                             char_length);
  }

  // Once all partitions have been set, build the global numbering
  PDM_gnum_compute(gen_gnum);

  // Finally, retrieve the computed global id arrays
  PDM_g_num_t **gnum = malloc(sizeof(PDM_g_num_t *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum[i_part] = PDM_gnum_get(gen_gnum,
                                i_part);

  }

  // Deallocate the PDM_gen_gnum_t instance
  PDM_gnum_free(gen_gnum);

  if (verbose) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      printf("part %d: ", i_part);
      for (int i = 0; i < n_elts[i_part]; i++) {
        printf(PDM_FMT_G_NUM " ", gnum[i_part][i]);
      }
      printf("\n");
    }
  }

  /*
   * Free memory
   */
  free(char_length);
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(coords[i_part]);
    free(gnum  [i_part]);
  }
  free(n_elts);
  free(coords);
  free(gnum);




  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
