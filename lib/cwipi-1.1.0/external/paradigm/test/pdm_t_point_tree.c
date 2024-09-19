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
#include "pdm_point_tree_seq.h"
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
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_pts,
 double        *radius,
 int           *tree_type,
 char         **filename,
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
        long n = atol(argv[i]);
        *gn_pts = (PDM_g_num_t) n;
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

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *tree_type = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename = argv[i];
      }
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

  PDM_g_num_t               gn_pts    = 10;
  double                    radius    = 10.;
  PDM_doctree_local_tree_t  tree_type = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  char                     *filename  = NULL;
  int                       visu      = 0;

  _read_args(argc,
             argv,
             &gn_pts,
             &radius,
     (int *) &tree_type,
             &filename,
             &visu);


  /* Random point cloud */
  /* Generate src and tgt point clouds */
  int          n_pts     = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  if (filename == NULL) {
    PDM_point_cloud_gen_random (comm,
                                0, // seed
                                0, // geometric_g_num
                                gn_pts,
                                -radius, -radius, -radius,
                                radius, radius, radius,
                                &n_pts,
                                &pts_coord,
                                &pts_g_num);
  }
  else {
    PDM_dmesh_nodal_t *dmn = PDM_reader_stl_dmesh_nodal(comm,
                                                        filename);

    const PDM_g_num_t *distrib = PDM_DMesh_nodal_distrib_vtx_get(dmn);

    n_pts = (int) (distrib[i_rank+1] - distrib[i_rank]);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

    pts_coord = malloc(sizeof(double) * n_pts * 3);
    memcpy(pts_coord, dvtx_coord, sizeof(double) * n_pts * 3);

    pts_g_num = malloc(sizeof(PDM_g_num_t) * n_pts);
    for (int i = 0; i < n_pts; i++) {
      pts_g_num[i] = distrib[i_rank] + i + 1;
    }

    PDM_DMesh_nodal_free(dmn);
  }



  /* Rearrange points (Hilbert) */
  double *weight =  malloc( n_pts * sizeof(double));
  for(int i = 0; i < n_pts; ++i) {
    weight[i] = 1.;
  }
  PDM_MPI_Barrier(comm);
  // double t1 = PDM_MPI_Wtime();
  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           PDM_PART_GEOM_HILBERT,
                                                           &pts_coord,
                                                           &pts_g_num,
                                                           &weight,
                                                           &n_pts,
                                                           1,
                                                           comm);
  free(weight);
  // double t2 = PDM_MPI_Wtime();
  // log_trace("PDM_part_to_block_geom_create = %12.5e \n", t2 -t1);

  double *blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &pts_coord,
                         NULL,
               (void **) &blk_pts_coord);

  int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  // PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);





  /* Build point tree */
  int depth_max          = 31;
  int points_in_leaf_max = 10;
  const double tolerance = 1e-4;
  PDM_point_tree_seq_t *ptree = PDM_point_tree_seq_create(tree_type,
                                                          depth_max,
                                                          points_in_leaf_max,
                                                          tolerance);

  PDM_point_tree_seq_point_cloud_set(ptree,
                                     n_parent,
                                     blk_pts_coord);

  PDM_point_tree_seq_build(ptree);

  if (visu) {
    char filename2[999];
    sprintf(filename2, "ptree_%d_%i.vtk", (int) tree_type, i_rank);
    PDM_point_tree_seq_write_nodes(ptree, filename2);

    sprintf(filename2, "points_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename2,
                              n_parent,
                              blk_pts_coord,
                              parent_gnum,
                              NULL);
  }

  /* Free */
  PDM_point_tree_seq_free(ptree);
  PDM_part_to_block_free(ptb);

  free(blk_pts_coord);

  free (pts_coord);
  free (pts_g_num);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
