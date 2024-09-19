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
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"

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
 int           *visu
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

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
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

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t nPts   = 10;
  double      radius = 10.;
  int         local  = 0;
  int         rand   = 0;
  int         visu   = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand,
             &visu);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }

  /* Random point cloud */
  /* Generate src and tgt point clouds */
  int          n_src     = 0;
  double      *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              nPts,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_src,
                              &src_coord,
                              &src_g_num);

  int depth_max = 5;
  int points_in_leaf_max = 1;
  const double tolerance = 1e-4;
  PDM_octree_seq_t *oct_orig = PDM_octree_seq_create(1, // n_point_cloud
                                                     depth_max,
                                                     points_in_leaf_max,
                                                     tolerance);

  PDM_octree_seq_point_cloud_set(oct_orig,
                                 0,
                                 n_src,
                                 src_coord);
  PDM_octree_seq_build(oct_orig);

  if(visu) {
    char filename[999];
    sprintf(filename, "octree_orig_%i.vtk", i_rank);
    PDM_octree_seq_write_octants(oct_orig, filename);

  }

  PDM_octree_seq_free(oct_orig);

  double *weight =  malloc( n_src * sizeof(double));
  for(int i = 0; i < n_src; ++i) {
    weight[i] = 1.;
  }
  PDM_MPI_Barrier(comm);
  // double t1 = PDM_MPI_Wtime();
  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           PDM_PART_GEOM_HILBERT,
                                                           &src_coord,
                                                           &src_g_num,
                                                           &weight,
                                                           &n_src,
                                                           1,
                                                           comm);
  free(weight);
  // double t2 = PDM_MPI_Wtime();
  // log_trace("PDM_part_to_block_geom_create = %12.5e \n", t2 -t1);

  double *blk_src_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &src_coord,
                         NULL,
               (void **) &blk_src_coord);

  int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  // PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  // PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  PDM_part_to_block_free(ptb);

  /*
   * Create octree seq
   */
  PDM_octree_seq_t *oct_equi = PDM_octree_seq_create(1, // n_point_cloud
                                                     depth_max,
                                                     points_in_leaf_max,
                                                     tolerance);

  PDM_octree_seq_point_cloud_set(oct_equi,
                                 0,
                                 n_parent,
                                 blk_src_coord);
  PDM_octree_seq_build(oct_equi);

  if(visu) {
    char filename[999];
    sprintf(filename, "octree_equi_%i.vtk", i_rank);
    PDM_octree_seq_write_octants(oct_equi, filename);
  }


  if (1) {
    /*
     *  Check intersection with boxes
     */

    PDM_g_num_t gn_box = 10;

    int          n_box        = 0;
    double      *box_extents  = NULL;
    PDM_g_num_t *box_ln_to_gn = NULL;

    if (0) {
      PDM_box_gen_cartesian(comm,
                            5,                         // nx
                            5,                         // ny
                            5,                         // nz
                            -radius, -radius, -radius, // x,y,z_min
                            radius, radius, radius,    // x,y,z_max
                            &n_box,
                            &box_extents,
                            &box_ln_to_gn);
    }
    else {
      PDM_box_gen_random(comm,
                         0,                         // seed
                         0,                         // geometric_g_num
                         gn_box,                    // gn_box
                         0.5*radius,                // min_size
                         1.5*radius,                // max_size
                         -radius, -radius, -radius, // x,y,z_min
                         radius, radius, radius,    // x,y,z_max
                         &n_box,
                         &box_extents,
                         &box_ln_to_gn);
    }


    if(visu) {
      char filename[999];
      sprintf(filename, "boxes_%i.vtk", i_rank);
      PDM_vtk_write_boxes(filename,
                          n_box,
                          box_extents,
                          box_ln_to_gn);

      sprintf(filename, "points_%i.vtk", i_rank);
      PDM_vtk_write_point_cloud(filename,
                                n_parent,
                                blk_src_coord,
                                NULL,
                                NULL);
    }

    // PDM_log_trace_array_long(box_ln_to_gn,
    //                          n_box,
    //                          "box_ln_to_gn : ");


    int *box_pts_idx = NULL;
    int *box_pts     = NULL;
    PDM_octree_seq_points_inside_boxes(oct_equi,
                                       n_box,
                                       box_extents,
                                       &box_pts_idx,
                                       &box_pts);

    // PDM_log_trace_connectivity_int(box_pts_idx,
    //                                box_pts,
    //                                n_box,
    //                                "box_pts : ");
    free(box_pts_idx);
    free(box_pts);

    free(box_extents);
    free(box_ln_to_gn);
  }

  PDM_octree_seq_free(oct_equi);


  free(blk_src_coord);


  /* Parallel octree */
  // const int n_point_cloud = 1;
  // const int depth_max = 31;
  // const int points_in_leaf_max = 1;

  // const int build_leaf_neighbours = 1;
  // PDM_para_octree_t *octree = PDM_para_octree_create (n_point_cloud,
  //                                                     depth_max,
  //                                                     points_in_leaf_max,
  //                                                     build_leaf_neighbours,
  //                                                     PDM_MPI_COMM_WORLD);

  // PDM_para_octree_point_cloud_set (octree, 0, _n_pts_l, coords, gnum);

  // PDM_para_octree_build (octree, NULL);

  // // PDM_para_octree_dump (octree);

  // PDM_para_octree_dump_times (octree);

  // PDM_para_octree_free (octree);

  /* Free */

  free (src_coord);
  free (src_g_num);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
