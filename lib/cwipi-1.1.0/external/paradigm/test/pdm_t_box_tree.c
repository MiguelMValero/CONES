#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"

#include "pdm_vtk.h"
#include "pdm_array.h"
#include "pdm_box_tree.h"
#include "pdm_box.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

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
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           PDM_g_num_t   *n_vtx_seg_cloud,
           PDM_g_num_t   *n_vtx_seg_tgt,
           double        *length,
           int           *post)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg_cloud = atol(argv[i]);
        *n_vtx_seg_cloud = (PDM_g_num_t) _n_vtx_seg_cloud;
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg_tgt = atol(argv[i]);
        *n_vtx_seg_tgt = (PDM_g_num_t) _n_vtx_seg_tgt;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static
void
_generate_cartesian_cloud
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  double      **pt_coord_out,
  PDM_g_num_t **pt_gnum_out
)
{
  PDM_g_num_t n_pts = n_vtx_x * n_vtx_y * n_vtx_z;
  double      *pt_coord = malloc( 3 * n_pts * sizeof(double));
  PDM_g_num_t *pt_gnum  = malloc(     n_pts * sizeof(PDM_g_num_t));

  double step_x = length / (double) (n_vtx_x - 1);
  double step_y = length / (double) (n_vtx_y - 1);
  double step_z = length / (double) (n_vtx_z - 1);

  for (int i_vtx = 0; i_vtx < n_pts; ++i_vtx) {


    pt_gnum[i_vtx] = i_vtx;

    PDM_g_num_t indi = i_vtx % n_vtx_x;
    PDM_g_num_t indj = ((i_vtx - indi) / n_vtx_x) % n_vtx_y;
    PDM_g_num_t indk = i_vtx / (n_vtx_x * n_vtx_y);

    pt_coord[3 * i_vtx    ] = indi * step_x; //+ dcube->zero_x;
    pt_coord[3 * i_vtx + 1] = indj * step_y; //+ dcube->zero_y;
    pt_coord[3 * i_vtx + 2] = indk * step_z; //+ dcube->zero_z;

  }

  *pt_coord_out = pt_coord;
  *pt_gnum_out  = pt_gnum;

}


static
void
_generate_cartesian_boxes
(
  int           n_vtx_x,
  int           n_vtx_y,
  int           n_vtx_z,
  double        length,
  int          *n_box_out,
  double      **box_coord_out,
  PDM_g_num_t **box_gnum_out
)
{
  PDM_g_num_t n_box = (n_vtx_x - 1) * ( n_vtx_y - 1) * ( n_vtx_z - 1);
  // PDM_g_num_t n_box = (n_vtx_x) * ( n_vtx_y) * ( n_vtx_z);
  *n_box_out = n_box;
  double      *box_coord = malloc( 6 * n_box * sizeof(double));
  PDM_g_num_t *box_gnum  = malloc(     n_box * sizeof(PDM_g_num_t));

  double step_x = length / (double) (n_vtx_x - 1);
  double step_y = length / (double) (n_vtx_y - 1);
  double step_z = length / (double) (n_vtx_z - 1);

  int n_box_x = n_vtx_x - 1;
  int n_box_y = n_vtx_y - 1;
  // int n_box_z = n_vtx_z - 1;

  for (int i_box = 0; i_box < n_box; ++i_box) {


    box_gnum[i_box] = i_box;

    PDM_g_num_t ind_box_i = i_box % n_box_x;
    PDM_g_num_t ind_box_j = ((i_box - ind_box_i) / n_box_x) % n_box_y;
    PDM_g_num_t ind_box_k = i_box / (n_box_x * n_box_y);

    int i_vtx = ind_box_i + ind_box_j * n_vtx_x + ind_box_k * n_vtx_x * n_vtx_y;

    PDM_g_num_t indi = i_vtx % n_vtx_x;
    PDM_g_num_t indj = ((i_vtx - indi) / n_vtx_x) % n_vtx_y;
    PDM_g_num_t indk = i_vtx / (n_vtx_x * n_vtx_y);

    box_coord[6 * i_box    ] = indi * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 1] = indj * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 2] = indk * step_z; //+ dcube->zero_z;

    box_coord[6 * i_box + 3] = (indi+1) * step_x; //+ dcube->zero_x;
    box_coord[6 * i_box + 4] = (indj+1) * step_y; //+ dcube->zero_y;
    box_coord[6 * i_box + 5] = (indk+1) * step_z; //+ dcube->zero_z;

  }

  *box_coord_out = box_coord;
  *box_gnum_out  = box_gnum;

}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t           n_vtx_seg       = 10;
  PDM_g_num_t           n_vtx_seg_cloud = 10;
  PDM_g_num_t           n_vtx_seg_tgt   = 10;
  double                length          = 1.;
  int                   post            = 0;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_vtx_seg_cloud,
             &n_vtx_seg_tgt,
             &length,
             &post);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Generate cloud
   */
  double      *pt_coord = NULL;
  PDM_g_num_t *pt_gnum  = NULL;
  PDM_g_num_t n_pts = n_vtx_seg_cloud * n_vtx_seg_cloud * n_vtx_seg_cloud;
  _generate_cartesian_cloud(n_vtx_seg_cloud,
                            n_vtx_seg_cloud,
                            n_vtx_seg_cloud,
                            length,
                            &pt_coord,
                            &pt_gnum);
  for(int i = 0; i < 3*n_pts; ++i) {
    pt_coord[i] = 0.25 * pt_coord[i];
  }


  if(post) {
    PDM_vtk_write_point_cloud("pts_cloud.vtk", n_pts, pt_coord, pt_gnum, NULL);
  }

  /*
   *  Generate target_boxes
   */
  double      *tgt_box_extents = NULL;
  PDM_g_num_t *tgt_box_gnum    = NULL;
  int          n_tgt_box = 0;
  _generate_cartesian_boxes(n_vtx_seg_tgt,
                            n_vtx_seg_tgt,
                            n_vtx_seg_tgt,
                            length,
                            &n_tgt_box,
                            &tgt_box_extents,
                            &tgt_box_gnum);
  // for(int i = 0; i < 6*n_tgt_box; ++i) {
  //   tgt_box_extents[i] = 0.25 * tgt_box_extents[i];
  // }

  if(post) {
    PDM_vtk_write_boxes("target_boxes.vtk", n_tgt_box, tgt_box_extents, tgt_box_gnum);
  }

  /*
   *  Generate boxes
   */
  double      *box_extents = NULL;
  PDM_g_num_t *box_gnum    = NULL;
  int n_box = 0;
  _generate_cartesian_boxes(n_vtx_seg,
                            n_vtx_seg,
                            n_vtx_seg,
                            length,
                            &n_box,
                            &box_extents,
                            &box_gnum);

  int perturb = 0;
  if(perturb == 1) {
    for(int i = 0; i < n_box; ++i) {
      // box_extents[i] = 0.5 * box_extents[i]+0.5;
      for(int j = 0; j < 3; ++j) {
        double mid  = 0.5 * (box_extents[6*i+j] + box_extents[6*i+j+3]);
        double size = 0.5 * (box_extents[6*i+j+3] - box_extents[6*i+j]);
        box_extents[6*i+j  ] = mid - 0.5 *size;
        box_extents[6*i+j+3] = mid + 0.5 *size;
      }
    }
  }

  if(post) {
    PDM_vtk_write_boxes("boxes.vtk", n_box, box_extents, box_gnum);
  }


  /*
   * Start box tree
   */
  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 4; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
  // int   max_boxes_leaf_shared = 1; // Max number of boxes in a leaf for coarse shared BBTree
  // int   max_tree_depth_shared = 4; // Max tree depth for coarse shared BBTree
  // float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)

  const int n_info_location = 3;
  int *init_location_proc = PDM_array_zeros_int (n_info_location * n_box);

  PDM_box_set_t* box_set = PDM_box_set_create(3,             // dim
                                              1,             // normalize
                                              0,             // allow_projection
                                              n_box,
                                              box_gnum,
                                              box_extents,
                                              1,
                                              &n_box,
                                              init_location_proc,
                                              PDM_MPI_COMM_WORLD);

  PDM_box_tree_t *box_tree = PDM_box_tree_create (max_tree_depth_shared,
                                                  max_boxes_leaf_shared,
                                                  max_box_ratio_shared);

  PDM_box_tree_set_boxes (box_tree,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  if(post) {
    PDM_box_tree_write_vtk("box_tree.vtk", box_tree, -1, 0);
  }

  int *init_location_tgt_box = PDM_array_zeros_int (n_info_location * n_tgt_box);
  PDM_box_set_t* box_target = PDM_box_set_create(3,             // dim
                                                 0,             // normalize
                                                 0,             // allow_projection
                                                 n_tgt_box,
                                                 tgt_box_gnum,
                                                 tgt_box_extents,
                                                 1,
                                                 &n_tgt_box,
                                                 init_location_tgt_box,
                                                 PDM_MPI_COMM_WORLD);

  /*
   * Intersect
   */
  double t1 = PDM_MPI_Wtime();
  int  *shared_to_box_idx = NULL;
  int  *shared_to_box     = NULL;
  // PDM_box_tree_get_boxes_intersects (box_tree,
  //                                    box_target,
  //                                    &shared_to_box_idx,
  //                                    &shared_to_box);
  // PDM_box_tree_intersect_boxes_boxes(box_tree,
  //                                    -1,
  //                                    n_tgt_box,
  //                                    tgt_box_extents,
  //                                    &shared_to_box_idx,
  //                                    &shared_to_box);
  PDM_box_tree_intersect_boxes_boxes2(box_tree,
                                      -1,
                                      n_tgt_box,
                                      tgt_box_extents,
                                      &shared_to_box_idx,
                                      &shared_to_box);

  if(post) {
    PDM_log_trace_connectivity_int(shared_to_box_idx, shared_to_box, n_box, "shared_to_box :");
  }

  double dt = PDM_MPI_Wtime()-t1;
  double complexity = n_tgt_box;
  printf("PDM_box_tree_get_boxes_intersects time : %12.5e -  %12.5e - complexity = %12.5e \n", dt, dt/complexity, complexity);

  free(shared_to_box_idx);
  free(shared_to_box    );

  t1 = PDM_MPI_Wtime();
  int  *box_pts_idx = NULL;
  int  *box_pts     = NULL;
  PDM_box_tree_boxes_containing_points(box_tree,
                                       -1,
                                       n_pts,
                                       pt_coord,
                                       &box_pts_idx,
                                       &box_pts);
  dt = PDM_MPI_Wtime()-t1;
  complexity = (double) n_pts;
  printf("PDM_box_tree_boxes_containing_points time : %12.5e - %12.5e - complexity = %12.5e \n", dt, dt/complexity, complexity);

  free(box_pts_idx);
  free(box_pts);


  PDM_box_set_destroy (&box_set);
  PDM_box_set_destroy (&box_target);
  PDM_box_tree_destroy(&box_tree);

  free(init_location_proc);
  free(init_location_tgt_box);
  free(box_extents);
  free(box_gnum);
  free(tgt_box_extents);
  free(tgt_box_gnum);
  free(pt_coord);
  free(pt_gnum);


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
