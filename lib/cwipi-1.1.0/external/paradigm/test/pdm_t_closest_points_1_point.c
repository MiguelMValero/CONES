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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"
#include "pdm_closest_points.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_pts,
 double        *length,
 int           *n_closest_pts,
 int           *post
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) _n;
      }
    }

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_closest_pts = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static void
_random_pts
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_pts,
 double              length,
 int                *n_pts,
 PDM_g_num_t       **pts_g_num,
 double            **pts_coord
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  *n_pts = (int) (gn_pts / n_rank);
  if (i_rank < gn_pts % n_rank) {
    (*n_pts)++;
  }

  PDM_g_num_t* distrib_pts = PDM_compute_entity_distribution(comm, (*n_pts));
  for(int i = 0; i < 3 * distrib_pts[i_rank]; ++i) {
    rand();
  }

  *pts_g_num = malloc (sizeof(PDM_g_num_t) * (*n_pts));
  for (int i = 0; i < *n_pts; i++) {
    (*pts_g_num)[i] = distrib_pts[i_rank] + i + 1;
  }

  *pts_coord = malloc (sizeof(double) * (*n_pts) * 3);
  for (int i = 0; i < *n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      (*pts_coord)[3*i + j] = length * (double) rand() / ((double) RAND_MAX);
    }
  }

  free (distrib_pts);
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main (int argc, char *argv[])
{
  srand(0);

  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);



  PDM_g_num_t gn_pts           = 10;
  double      length           = 3.14;
  int         n_closest_points = 1;
  int         post             = 0;

  _read_args (argc,
              argv,
              &gn_pts,
              &length,
              &n_closest_points,
              &post);



  /*
   *  Random boxes
   */
  int n_tgt;
  PDM_g_num_t *tgt_g_num = NULL;
  double      *tgt_coord = NULL;
  _random_pts (comm,
               gn_pts,
               length,
               &n_tgt,
               &tgt_g_num,
               &tgt_coord);

  if (post) {
    for (int i = 0; i < n_tgt; i++) {
      log_trace("tgt "PDM_FMT_G_NUM" coord = %f %f %f\n",
                tgt_g_num[i],
                tgt_coord[3*i + 0],
                tgt_coord[3*i + 1],
                tgt_coord[3*i + 2]);
    }
  }


  /*
   *  Random point
   */
  int n_src = 0;
  double      *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;
  if (i_rank == 0) {
    n_src = 1;
    src_coord = malloc (sizeof(double) * n_src * 3);
    src_g_num = malloc (sizeof(PDM_g_num_t) * n_src);

    for (int i = 0; i < 3; i++) {
      src_coord[i] = 0.5 * length;
    }

    src_g_num[0] = 1;
  }



  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         n_closest_points,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set (clsp,
                                       1,
                                       1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coord,
                                    src_g_num);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num);


  PDM_closest_points_compute (clsp);


  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;
  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  int *closest_src_idx = malloc (sizeof(int) * (n_tgt + 1));
  for (int i = 0; i <= n_tgt; i++) {
    closest_src_idx[i] = n_closest_points * i;
  }
  if (post) {
    PDM_log_trace_connectivity_long (closest_src_idx,
                                     closest_src_gnum,
                                     n_tgt,
                                     "closest_src_gnum : ");
  }

  for (int i = 0; i < closest_src_idx[n_tgt]; i++) {
    //assert(closest_src_gnum[i] == 1);
  }

  PDM_closest_points_free (clsp);

  free (tgt_coord);
  free (tgt_g_num);
  free (closest_src_idx);

  if (i_rank == 0) {
    printf("-- End\n");
    free (src_coord);
    free (src_g_num);
  }


  PDM_MPI_Finalize();

  return 0;
}
