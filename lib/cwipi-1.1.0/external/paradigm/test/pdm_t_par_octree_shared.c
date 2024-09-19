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
#include "pdm_point_cloud_gen.h"
#include "pdm_para_octree.h"
#include "pdm_logging.h"

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
 PDM_g_num_t   *nPts,
 double        *radius,
 int           *local,
 int           *rand,
 int           *shared
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
        long _nPts = atol(argv[i]);
        *nPts = (PDM_g_num_t) _nPts;
      }
    }

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else if (strcmp(argv[i], "-shared") == 0) {
      *shared = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

/* _random01 from pdm_t_intersect_line_box */

static double
_random01
(
 void
)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS (rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
  // return (double) rand() / (double) RAND_MAX;
}

/* code from pdm_t_intersect_line_box */

static int
_generate_random_boxes
(
PDM_MPI_Comm  comm,
PDM_g_num_t   gn_box,
int           i_rank,
double      **box_extents,
PDM_g_num_t **box_ln_to_gn
)
{
  PDM_g_num_t *distrib_box = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_box);
  int n_box = (int) (distrib_box[i_rank+1] - distrib_box[i_rank]);
  for (PDM_g_num_t i = 0; i < 6*distrib_box[i_rank]; i++) {
    rand();
  }
  free (distrib_box);

  double *box_centers = malloc (sizeof(double) * n_box * 3);
  *box_extents = malloc (sizeof(double) * n_box * 6);
  double *_box_extents = *box_extents;
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = _random01();
      double x2 = _random01();

      box_centers[3*i + j] = 0.5 * (x1 + x2);
      _box_extents[6*i + j]     = PDM_MIN (x1, x2);
      _box_extents[6*i + j + 3] = PDM_MAX (x1, x2);
    }
  }


  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3,
                                              1,
                                              PDM_FALSE,
                                              1.e-3,
                                              comm,
                                              PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_coords (gen_gnum,
                            0,
                            n_box,
                            box_centers,
                            NULL);

  PDM_gnum_compute (gen_gnum);

  *box_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (box_centers);

  return n_box;
}

static
void
end_timer_and_print(const char* msg, PDM_MPI_Comm comm, double t1){

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t,
                     &delta_max,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     comm);

  PDM_MPI_Allreduce (&delta_t,
                     &delta_min,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MIN,
                     comm);

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);
  if(i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e %12.5e \n", n_rank, msg, delta_min, delta_max);
  }
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
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t nPts   = 10;
  double radius = 10.;
  int local = 0;
  int rand = 0;
  int shared = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand,
             &shared);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }

  /* Random point cloud */
  int _n_pts_l;
  double      *coords = NULL;
  PDM_g_num_t *gnum   = NULL;
  PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                              0, // seed
                              0, // geometric_g_num
                              nPts,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &_n_pts_l,
                              &coords,
                              &gnum);

  PDM_g_num_t gn_box = nPts;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;

  int n_boxes = _generate_random_boxes(PDM_MPI_COMM_WORLD,
                                       gn_box,
                                       i_rank,
                                       &box_extents,
                                       &box_ln_to_gn);

  /* Parallel octree */

  const int n_point_cloud = 1;
  const int depth_max = 31;
  const int points_in_leaf_max = 1;

  const int build_leaf_neighbours = 1;
  PDM_para_octree_t *octree = PDM_para_octree_create (n_point_cloud,
                                                      depth_max,
                                                      points_in_leaf_max,
                                                      build_leaf_neighbours,
                                                      PDM_MPI_COMM_WORLD);

  PDM_para_octree_point_cloud_set (octree, 0, _n_pts_l, coords, gnum);

  // PDM_para_octree_build (octree, NULL);

  // Mettre en shared l'octree
  // PDM_para_octree_dump (octree);

  // Solicitations shared
  int         *box_pts_idx = NULL;
  PDM_g_num_t *box_pts     = NULL;
  double      *pts_coords  = NULL;
  if(shared == 0) {
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
    double t1 = PDM_MPI_Wtime();
    PDM_para_octree_build (octree, NULL);
    end_timer_and_print("PDM_para_octree_build ", PDM_MPI_COMM_WORLD, t1 );

    double t2 = PDM_MPI_Wtime();
    PDM_para_octree_points_inside_boxes(octree,
                                        n_boxes,
                                        box_extents,
                                        box_ln_to_gn,
                                        &box_pts_idx,
                                        &box_pts,
                                        &pts_coords);
    end_timer_and_print("PDM_para_octree_points_inside_boxes ", PDM_MPI_COMM_WORLD, t2 );
  } else {
    PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
    double t1 = PDM_MPI_Wtime();
    PDM_para_octree_build_shared(octree, NULL);
    end_timer_and_print("PDM_para_octree_build_shared ", PDM_MPI_COMM_WORLD, t1 );

    double t2 = PDM_MPI_Wtime();
    PDM_para_octree_points_inside_boxes_shared(octree,
                                               n_boxes,
                                               box_extents,
                                               box_ln_to_gn,
                                               &box_pts_idx,
                                               &box_pts,
                                               &pts_coords);
    end_timer_and_print("PDM_para_octree_points_inside_boxes_shared ", PDM_MPI_COMM_WORLD, t2 );
  }

  // PDM_para_octree_build_shared(octree);

  PDM_para_octree_dump_times (octree);

  PDM_para_octree_free (octree);

  // PDM_log_trace_connectivity_long(box_pts_idx, box_pts, n_boxes, "box_pts :: ");

  /* Free */
  free(box_pts_idx);
  free(box_pts);
  free(pts_coords);

  free (coords);
  free (gnum);
  free (box_extents);
  free (box_ln_to_gn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
