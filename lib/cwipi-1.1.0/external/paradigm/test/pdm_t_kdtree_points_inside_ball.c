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
#include "pdm_kdtree_seq.h"
#include "pdm_octree_seq.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_reader_stl.h"



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
 double        *length,
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

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
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
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t nPts      = 10;
  double      length    = 10.;
  int         local     = 0;
  int         randomize = 0;
  int         visu      = 0;

  _read_args(argc,
             argv,
             &nPts,
             &length,
             &local,
             &randomize,
             &visu);

  /* Initialize random */

  if (randomize) {
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
  if (1) {
    PDM_point_cloud_gen_random (comm,
                                0, // seed
                                0, // geometric_g_num
                                nPts,
                                -length, -length, -length,
                                length, length, length,
                                &n_src,
                                &src_coord,
                                &src_g_num);

    // for (int i = 0; i < n_src; i++) {
    //   src_coord[3*i] *= 2;
    // }
  }
  else {
    PDM_dmesh_nodal_t *dmn = PDM_reader_stl_dmesh_nodal(comm,
                                                        PDM_MESH_DIR"sphere.stl");

    const PDM_g_num_t *distrib = PDM_DMesh_nodal_distrib_vtx_get(dmn);

    n_src = (int) (distrib[i_rank+1] - distrib[i_rank]);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

    src_coord = malloc(sizeof(double) * n_src * 3);
    memcpy(src_coord, dvtx_coord, sizeof(double) * n_src * 3);

    src_g_num = malloc(sizeof(PDM_g_num_t) * n_src);
    for (int i = 0; i < n_src; i++) {
      src_g_num[i] = distrib[i_rank] + i + 1;
    }

    PDM_DMesh_nodal_free(dmn);
  }



  int depth_max = 20;
  int points_in_leaf_max = 2;
  const double tolerance = 1e-4;
  PDM_kdtree_seq_t *kdt = PDM_kdtree_seq_create(depth_max,
                                                points_in_leaf_max,
                                                tolerance);

  PDM_kdtree_seq_point_cloud_set(kdt,
                                 n_src,
                                 src_coord);
  PDM_kdtree_seq_build(kdt);





  int          n_tgt     = 0;
  double      *tgt_coord = NULL;
  PDM_g_num_t *tgt_g_num = NULL;
  PDM_point_cloud_gen_random (comm,
                              123456789, // seed
                              0,         // geometric_g_num
                              nPts,
                              -length, -length, -length,
                              length, length, length,
                              &n_tgt,
                              &tgt_coord,
                              &tgt_g_num);

  // for (int i = 0; i < n_tgt; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     tgt_coord[3*i+j] += 0.1*length*(2*((double) rand() / (double) RAND_MAX) - 1);
  //   }
  // }

  double *ball_radius2 = malloc(sizeof(double) * n_tgt);
  for (int i = 0; i < n_tgt; i++) {
    double r = 0.8*length;
    ball_radius2[i] = r*r;
  }

  int    *pts_inside_ball_idx   = NULL;
  int    *pts_inside_ball_l_num = NULL;
  double *pts_inside_ball_dist2 = NULL;
  PDM_kdtree_seq_points_inside_balls(kdt,
                                     n_tgt,
                                     tgt_coord,
                                     ball_radius2,
                                     &pts_inside_ball_idx,
                                     &pts_inside_ball_l_num,
                                     &pts_inside_ball_dist2);

  if(visu) {
    char filename[999];
    sprintf(filename, "kdtree_%i.vtk", i_rank);
    PDM_kdtree_seq_write_nodes(kdt, filename);

    sprintf(filename, "src_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_src,
                              src_coord,
                              src_g_num,
                              NULL);

    sprintf(filename, "tgt_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_tgt,
                              tgt_coord,
                              tgt_g_num,
                              NULL);
  }

  // log_trace("> Kdtree\n");
  // for (int i = 0; i < n_tgt; i++) {
  //   log_trace("point %d (%f %f %f): n_pib = %d\n",
  //             i,
  //             tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2],
  //             pts_inside_ball_idx[i+1] - pts_inside_ball_idx[i]);
  //   for (int j = pts_inside_ball_idx[i]; j < pts_inside_ball_idx[i+1]; j++) {
  //     log_trace("  l_num %d, at dist2 %f / %f\n",
  //               pts_inside_ball_l_num[j],
  //               pts_inside_ball_dist2[j], ball_radius2[i]);
  //   }
  // }
  free(pts_inside_ball_idx);
  free(pts_inside_ball_l_num);
  free(pts_inside_ball_dist2);





  PDM_octree_seq_t *oct = PDM_octree_seq_create(1, // n_point_cloud
                                                depth_max,
                                                points_in_leaf_max,
                                                tolerance);

  PDM_octree_seq_point_cloud_set(oct,
                                 0,
                                 n_src,
                                 src_coord);
  PDM_octree_seq_build(oct);

  PDM_octree_seq_points_inside_balls(oct,
                                     n_tgt,
                                     tgt_coord,
                                     ball_radius2,
                                     &pts_inside_ball_idx,
                                     &pts_inside_ball_l_num,
                                     &pts_inside_ball_dist2);

  // log_trace("> Octree\n");
  // for (int i = 0; i < n_tgt; i++) {
  //   log_trace("point %d (%f %f %f): n_pib = %d\n",
  //             i,
  //             tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2],
  //             pts_inside_ball_idx[i+1] - pts_inside_ball_idx[i]);
  //   for (int j = pts_inside_ball_idx[i]; j < pts_inside_ball_idx[i+1]; j++) {
  //     // log_trace("  cloud %d, id %d, at dist2 %f / %f\n",
  //               // pts_inside_ball_l_num[2*j], pts_inside_ball_l_num[2*j+1],
  //     log_trace("  l_num %d, at dist2 %f / %f\n",
  //               pts_inside_ball_l_num[2*j+1],
  //               pts_inside_ball_dist2[j], ball_radius2[i]);
  //   }
  // }

  // /* Brute force check */
  // for (int itgt = 0; itgt < n_tgt; itgt++) {

  //   for (int isrc = 0; isrc < n_src; isrc++) {
  //     double dist2 = 0.;
  //     for (int j = 0; j < 3; j++) {
  //       double delta = tgt_coord[3*itgt + j] - src_coord[3*isrc + j];
  //       dist2 += delta*delta;
  //     }
  //   }
  // }
  if (1) {
    int    *closest_kdtree_pt_id    = malloc(sizeof(int   ) * n_tgt);
    double *closest_kdtree_pt_dist2 = malloc(sizeof(double) * n_tgt);
    PDM_kdtree_seq_closest_point(kdt,
                                 n_tgt,
                                 tgt_coord,
                                 closest_kdtree_pt_id,
                                 closest_kdtree_pt_dist2);

    int    *closest_octree_pt_id    = malloc(sizeof(int   ) * n_tgt * 2);
    double *closest_octree_pt_dist2 = malloc(sizeof(double) * n_tgt);
    PDM_octree_seq_closest_point(oct,
                                 n_tgt,
                                 tgt_coord,
                                 closest_octree_pt_id,
                                 closest_octree_pt_dist2);


    // for (int i = 0; i < n_tgt; i++) {
    //   log_trace("%5d : %5d %5d %f / %5d %5d %f\n",
    //             i,
    //             // closest_kdtree_pt_id[2*i], closest_kdtree_pt_id[2*i+1], closest_kdtree_pt_dist2[i],
    //             0, closest_kdtree_pt_id[i], closest_kdtree_pt_dist2[i],
    //             closest_octree_pt_id[2*i], closest_octree_pt_id[2*i+1], closest_octree_pt_dist2[i]);
    // }

    free(closest_kdtree_pt_id);
    free(closest_kdtree_pt_dist2);
    free(closest_octree_pt_id);
    free(closest_octree_pt_dist2);
  }




  /* Free */
  PDM_kdtree_seq_free(kdt);
  PDM_octree_seq_free(oct);
  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);
  free(ball_radius2);
  free(pts_inside_ball_idx);
  free(pts_inside_ball_l_num);
  free(pts_inside_ball_dist2);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
