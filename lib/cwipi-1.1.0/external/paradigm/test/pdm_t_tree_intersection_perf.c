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
#include "pdm_point_cloud_gen.h"
#include "pdm_box_gen.h"
#include "pdm_point_tree_seq.h"
#include "pdm_box_tree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"


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
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_pts,
 int           *points_in_leaf_max,
 int           *cartesian_pts,
 int           *cartesian_box
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
        long n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-points_in_leaf_max") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *points_in_leaf_max = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-pcart") == 0) {
      *cartesian_pts = 1;
    }

    else if (strcmp(argv[i], "-bcart") == 0) {
      *cartesian_box = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
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
 int   argc,
 char *argv[]
 )
{
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t gn_pts             = 10;
  double      radius             = 1000.;
  int         points_in_leaf_max = 10;
  int         cartesian_pts      = 0;
  int         cartesian_box      = 0;


  _read_args(argc,
             argv,
             &gn_pts,
             &points_in_leaf_max,
             &cartesian_pts,
             &cartesian_box);

  PDM_g_num_t gn_box = gn_pts;

  double t1, t2;

  double _n = PDM_MAX(2, 1 + pow(gn_pts, 1./3.));

  /* Generate point cloud */
  int          n_pts     = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;

  if (cartesian_pts) {
    PDM_point_cloud_gen_cartesian(comm,
                                  (int) _n,
                                  (int) _n,
                                  (int) _n,
                                  -radius, -radius, -radius,
                                  radius, radius, radius,
                                  &n_pts,
                                  &pts_coord,
                                  &pts_g_num);
  }
  else {
    PDM_point_cloud_gen_random(comm,
                               0, // seed
                               0, // geometric_g_num
                               gn_pts,
                               -radius, -radius, -radius,
                               radius, radius, radius,
                               &n_pts,
                               &pts_coord,
                               &pts_g_num);
  }


  /* Generate boxes */
  int          n_box       = 0;
  double      *box_extents = NULL;
  PDM_g_num_t *box_g_num   = NULL;

  if (cartesian_box) {
    PDM_box_gen_cartesian(comm,
                          (int) _n,
                          (int) _n,
                          (int) _n,
                          -radius, -radius, -radius,
                          radius, radius, radius,
                          &n_box,
                          &box_extents,
                          &box_g_num);
  }
  else {
    double avg_size = 2*radius/(_n - 1);
    double min_size = 0.5*avg_size;
    double max_size = 1.5*avg_size;
    PDM_box_gen_random(comm,
                       0,
                       0,
                       gn_box,
                       min_size,
                       max_size,
                       -radius, -radius, -radius,
                       radius, radius, radius,
                       &n_box,
                       &box_extents,
                       &box_g_num);
  }







  double t_build_ptree[2];
  double t_build_btree;
  double t_sollicitation[5];

  /* Build point trees */
  PDM_doctree_local_tree_t tree_type[2] = {
    PDM_DOCTREE_LOCAL_TREE_OCTREE,
    PDM_DOCTREE_LOCAL_TREE_KDTREE
  };

  int depth_max          = 31;
  const double tolerance = 1e-6;

  PDM_point_tree_seq_t *ptree[2];

  for (int iptree = 0; iptree < 2; iptree++) {
    ptree[iptree] = PDM_point_tree_seq_create(tree_type[iptree],
                                              depth_max,
                                              points_in_leaf_max,
                                              tolerance);

    PDM_point_tree_seq_point_cloud_set(ptree[iptree],
                                       n_pts,
                                       pts_coord);

    t1 = PDM_MPI_Wtime();
    PDM_point_tree_seq_build(ptree[iptree]);
    t2 = PDM_MPI_Wtime();

    t_build_ptree[iptree] = t2 - t1;
  }



  /* Build box tree */
  int *init_location_box = malloc(3 * n_box * sizeof(int));
  for(int i = 0; i < n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  PDM_box_set_t *box_set = PDM_box_set_create(3,
                                              1,
                                              0,  // No projection to preserve initial extents
                                              n_box,
                                              box_g_num,
                                              box_extents,
                                              1,
                                              &n_box,
                                              init_location_box,
                                              comm);

  int   max_boxes_leaf = 30;
  int   max_tree_depth = 10;
  float max_box_ratio  = 30;

  PDM_box_tree_t *btree = PDM_box_tree_create (max_tree_depth,
                                               max_boxes_leaf,
                                               max_box_ratio);

  t1 = PDM_MPI_Wtime();
  PDM_box_tree_set_boxes (btree,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);
  t2 = PDM_MPI_Wtime();
  t_build_btree = t2 - t1;
  free(init_location_box);


  /* Sollicitations */
  int *box_pts_idx = NULL;
  int *box_pts     = NULL;
  PDM_g_num_t *box_pts_g_num = NULL;
  double      *box_pts_coord = NULL;

  /* Point tree / box tree */
  for (int iptree = 0; iptree < 2; iptree++) {
    t1 = PDM_MPI_Wtime();
    PDM_tree_intersection_point_box(btree,
                                    ptree[iptree],
                                    &box_pts_idx,
                                    &box_pts);
    t2 = PDM_MPI_Wtime();
    t_sollicitation[iptree] = t2 - t1;

    free(box_pts_idx);
    free(box_pts);
  }


  /* Point tree */
  for (int iptree = 0; iptree < 2; iptree++) {
    t1 = PDM_MPI_Wtime();
    PDM_point_tree_seq_points_inside_boxes(ptree[iptree],
                                           n_box,
                                           box_extents,
                                           &box_pts_idx,
                                           &box_pts);
    t2 = PDM_MPI_Wtime();
    t_sollicitation[2+iptree] = t2 - t1;
    free(box_pts_idx);
    free(box_pts);
  }


  /* Box tree */
  t1 = PDM_MPI_Wtime();
  PDM_box_tree_points_inside_boxes(btree,
                                   n_pts,
                                   pts_g_num,
                                   pts_coord,
                                   &box_pts_idx,
                                   &box_pts_g_num,
                                   &box_pts_coord);
  t2 = PDM_MPI_Wtime();
  t_sollicitation[4] = t2 - t1;

  free(box_pts_idx);
  free(box_pts_g_num);
  free(box_pts_coord);


  printf(PDM_FMT_G_NUM" %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
         gn_pts,
         t_build_ptree[0],
         t_build_ptree[1],
         t_build_btree,
         t_sollicitation[0],
         t_sollicitation[1],
         t_sollicitation[2],
         t_sollicitation[3],
         t_sollicitation[4]);


  /* Free memory */
  PDM_point_tree_seq_free(ptree[0]);
  PDM_point_tree_seq_free(ptree[1]);
  PDM_box_tree_destroy(&btree);
  PDM_box_set_destroy (&box_set);

  free(pts_coord);
  free(pts_g_num);

  free(box_extents);
  free(box_g_num);


  PDM_MPI_Finalize ();

  return 0;
}
