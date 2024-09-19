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
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_doctree.h"
#include "pdm_box_gen.h"

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
 int           *rand
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
  double radius = 10.;
  int local = 0;
  int rand = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand);

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
                              0., 0., 0.,
                              1., 1., 1.,
                              &n_src,
                              &src_coord,
                              &src_g_num);

  // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
  PDM_doctree_t *doct = PDM_doctree_create(comm,
                                           3,
                                           1,
                                           NULL, // global_extents
                                           local_tree_kind);

  int *init_location_pts = malloc(3 * n_src * sizeof(int));
  for(int i = 0; i < n_src; ++i) {
    init_location_pts[3*i  ] = i_rank;
    init_location_pts[3*i+1] = 0; // i_part
    init_location_pts[3*i+2] = i;
  }

  PDM_doctree_point_set(doct,
                        0,
                        n_src,
                        NULL,
                        src_g_num,
                        src_coord);
  int n_box = 0;
  int n_vtx_x = 10;
  int n_vtx_y = 10;
  int n_vtx_z = 10;
  double      *box_extents = NULL;
  PDM_g_num_t *box_gnum    = NULL;
  PDM_box_gen_cartesian(comm,
                        n_vtx_x,
                        n_vtx_y,
                        n_vtx_z,
                        0., 0., 0.,
                        1., 1., 1.,
                        &n_box,
                        &box_extents,
                        &box_gnum);

  int *init_location_box = malloc(3 * n_box * sizeof(int));
  for(int i = 0; i < n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  PDM_doctree_solicitation_set(doct,
                               PDM_TREE_SOLICITATION_BOXES_POINTS,
                               1,
                               &n_box,
                               &init_location_box,
                               &box_gnum,
                               &box_extents);

  PDM_doctree_build(doct);

  if(1 == 1) {

    int         *box_pts_idx = NULL;
    PDM_g_num_t *box_pts     = NULL;
    double      *pts_coord   = NULL;
    PDM_doctree_results_in_orig_frame_get(doct,
                                          n_box,
                                          box_gnum,
                                          &box_pts_idx,
                                          &box_pts,
                                          &pts_coord);

    free(box_pts_idx);
    free(box_pts    );
    free(pts_coord  );
  }


  PDM_doctree_free(doct);


  free(box_gnum);
  free(box_extents);
  free(init_location_box);
  free(init_location_pts);

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
