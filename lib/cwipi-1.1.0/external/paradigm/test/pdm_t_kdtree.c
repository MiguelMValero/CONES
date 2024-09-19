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
  if (1) {
    PDM_point_cloud_gen_random (comm,
                                0, // seed
                                0, // geometric_g_num
                                nPts,
                                -radius, -radius, -radius,
                                radius, radius, radius,
                                &n_src,
                                &src_coord,
                                &src_g_num);

    for (int i = 0; i < n_src; i++) {
      src_coord[3*i] *= 2;
    }
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



  int depth_max = 31;
  int points_in_leaf_max = 10;
  const double tolerance = 1e-4;
  PDM_kdtree_seq_t *kdt_orig = PDM_kdtree_seq_create(depth_max,
                                                     points_in_leaf_max,
                                                     tolerance);

  PDM_kdtree_seq_point_cloud_set(kdt_orig,
                                 n_src,
                                 src_coord);
  PDM_kdtree_seq_build(kdt_orig);

  int *pts_order = NULL;
  PDM_kdtree_seq_point_new_to_old_get(kdt_orig,
                                      &pts_order);

  double *pts_coord = NULL;
  PDM_kdtree_seq_sorted_points_get(kdt_orig,
                                   &pts_coord);

  if (visu) {
    char filename[999];
    sprintf(filename, "kdtree_orig_%i.vtk", i_rank);
    PDM_kdtree_seq_write_nodes(kdt_orig, filename);

    sprintf(filename, "points_orig_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_src,
                              pts_coord,//src_coord,
                              NULL,//src_g_num,
                              pts_order);
  }

  PDM_kdtree_seq_free(kdt_orig);











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




  PDM_kdtree_seq_t *kdt_equi = PDM_kdtree_seq_create(depth_max,
                                                     points_in_leaf_max,
                                                     tolerance);

  PDM_kdtree_seq_point_cloud_set(kdt_equi,
                                 n_parent,
                                 blk_src_coord);
  PDM_kdtree_seq_build(kdt_equi);

  PDM_kdtree_seq_point_new_to_old_get(kdt_equi,
                                      &pts_order);

  PDM_kdtree_seq_sorted_points_get(kdt_equi,
                                   &pts_coord);


  if(visu) {
    char filename[999];
    sprintf(filename, "kdtree_equi_%i.vtk", i_rank);
    PDM_kdtree_seq_write_nodes(kdt_equi, filename);

    sprintf(filename, "points_equi_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_parent,
                              pts_coord,//blk_src_coord,
                              NULL,//parent_gnum,
                              pts_order);
  }
  PDM_kdtree_seq_free(kdt_equi);


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
  PDM_octree_seq_free(oct_equi);

  PDM_part_to_block_free(ptb);

  free(blk_src_coord);




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
