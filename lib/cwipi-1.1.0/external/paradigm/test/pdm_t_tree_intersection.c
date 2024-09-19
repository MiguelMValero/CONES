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
#include "pdm_point_tree_seq_priv.h"
#include "pdm_box_tree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_unique.h"

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
 PDM_g_num_t   *gn_box,
 double        *radius,
 int           *tree_type,
 int           *points_in_leaf_max,
 int           *visu,
 int           *cartesian_pts,
 int           *cartesian_box,
 double        *randomization,
 int           *permutation,
 int           *rotation,
 char         **filename
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

    else if (strcmp(argv[i], "-b") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *gn_box = (PDM_g_num_t) n;
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

    else if (strcmp(argv[i], "-pil") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *points_in_leaf_max = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-pcart") == 0) {
      *cartesian_pts = 1;
    }

    else if (strcmp(argv[i], "-bcart") == 0) {
      *cartesian_box = 1;
    }

    else if (strcmp(argv[i], "-noise") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *randomization = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-perm") == 0) {
      *permutation = 1;
    }

    else if (strcmp(argv[i], "-rot") == 0) {
      *rotation = 1;
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


  PDM_g_num_t               gn_pts    = 10;
  PDM_g_num_t               gn_box    = 10;
  double                    radius    = 10.;
  PDM_doctree_local_tree_t  tree_type = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  int                       visu      = 0;
  int                       points_in_leaf_max = 10;
  int                       cartesian_pts      = 0;
  int                       cartesian_box      = 0;
  double                    randomization      = 0.;
  int                       permutation        = 0;
  int                       rotation           = 0;
  char                     *filename           = NULL;

  _read_args(argc,
             argv,
             &gn_pts,
             &gn_box,
             &radius,
     (int *) &tree_type,
             &points_in_leaf_max,
             &visu,
             &cartesian_pts,
             &cartesian_box,
             &randomization,
             &permutation,
             &rotation,
             &filename);


  double t1, t2;




  /* (Random) boxes */
  int n_box = 0;
  double      *box_extents = NULL;
  PDM_g_num_t *box_g_num   = NULL;

  double _n = PDM_MAX(2, 1 + pow(gn_box, 1./3.));

  if (filename != NULL) {

    printf("Read STL mesh...\n");
    PDM_dmesh_nodal_t *dmn = PDM_reader_stl_dmesh_nodal(comm,
                                                        filename);

    PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_SURFACIC;

    int *section_ids = PDM_DMesh_nodal_sections_id_get(dmn,
                                                       geom_kind);
    int id_section = section_ids[0];

    PDM_g_num_t *delt_vtx = PDM_DMesh_nodal_section_std_get(dmn,
                                                            geom_kind,
                                                            id_section);

    PDM_Mesh_nodal_elt_t elt_type = PDM_DMesh_nodal_section_elt_type_get(dmn,
                                                                       geom_kind,
                                                                       id_section);
    int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element(elt_type, 1);

    const PDM_g_num_t *distrib = PDM_DMesh_nodal_distrib_section_get(dmn,
                                                                     geom_kind,
                                                                     id_section);

    n_box = (int) (distrib[i_rank+1] - distrib[i_rank]);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

    box_extents = malloc(sizeof(double) * n_box * 6);
    for (int i = 0; i < n_box; i++) {

      double *e = box_extents + 6*i;

      for (int j = 0; j < 3; j++) {
        e[j  ] =  HUGE_VAL;
        e[j+3] = -HUGE_VAL;
      }

      for (int idx_vtx = n_vtx_elt*i; idx_vtx < n_vtx_elt*(i+1); idx_vtx++) {
        int vtx_id = (int) delt_vtx[idx_vtx] - 1;
        for (int j = 0; j < 3; j++) {
          e[j  ] = PDM_MIN(e[j  ], dvtx_coord[3*vtx_id+j]);
          e[j+3] = PDM_MAX(e[j+3], dvtx_coord[3*vtx_id+j]);
        }
      }

    }

    box_g_num = malloc(sizeof(PDM_g_num_t) * n_box);
    for (int i = 0; i < n_box; i++) {
      box_g_num[i] = distrib[i_rank] + i + 1;
    }

    PDM_DMesh_nodal_free(dmn);

  }
  else {
    if (cartesian_box) {
      int n_vtx_x = (int) _n;
      int n_vtx_y = (int) _n;
      int n_vtx_z = (int) _n;
      PDM_box_gen_cartesian(comm,
                            n_vtx_x,
                            n_vtx_y,
                            n_vtx_z,
                            -4*radius, -radius, -radius,
                            4*radius, radius, radius,
                            &n_box,
                            &box_extents,
                            &box_g_num);

      // Randomize
      double noise = randomization * 2*radius/(double) (_n - 1);
      for (int i = 0; i < n_box; i++) {
        for (int j = 0; j < 3; j++) {
          double r = noise * (2*((double) rand()/(double) RAND_MAX) - 1);
          box_extents[6*i+j  ] += r;
          box_extents[6*i+j+3] += r;
        }
      }

      // Random permutation
      if (permutation) {
        double *_box_extents = malloc(sizeof(double) * n_box * 6);
        int *rand_val = malloc(sizeof(int) * n_box);
        int *order    = malloc(sizeof(int) * n_box);
        for (int i = 0; i < n_box; i++) {
          order[i]    = i;
          rand_val[i] = rand();
        }

        PDM_sort_int(rand_val, order, n_box);
        free(rand_val);

        for (int i = 0; i < n_box; i++) {
          memcpy(_box_extents + 6*i, box_extents + 6*order[i], sizeof(double) * 6);
          box_g_num[order[i]] = i+1;
        }
        free(order);
        free(box_extents);
        box_extents = _box_extents;
      }

    }
    else {
      double avg_size = 2*radius/(double) (_n - 1);
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

    if (rotation) {
      double R[3][3] = {
        {0.9362934,  -0.2896295,  0.1986693},
        {0.3129918,   0.9447025, -0.0978434},
        {-0.1593451,  0.1537920,  0.9751703}};

      for (int i = 0; i < 2*n_box; i++) {
        double x = box_extents[3*i];
        double y = box_extents[3*i+1];
        double z = box_extents[3*i+2];

        for (int j = 0; j < 3; j++) {
          box_extents[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
        }
      }
    }
  }

  printf("n_box = %d\n", n_box);



  double g_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      g_extents[j  ] = PDM_MIN(g_extents[j  ], box_extents[6*i+j  ]);
      g_extents[j+3] = PDM_MIN(g_extents[j+3], box_extents[6*i+j+3]);
    }
  }




  /* Random point cloud */
  int          n_pts     = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;

  _n = PDM_MAX(2, 1 + pow(gn_pts, 1./3.));

  if (cartesian_pts) {
    PDM_point_cloud_gen_cartesian(comm,
                                  (int) _n,
                                  (int) _n,
                                  (int) _n,
                                  // -radius, -radius, -radius,
                                  // radius, radius, radius,
                                  g_extents[0], g_extents[1], g_extents[2],
                                  g_extents[3], g_extents[4], g_extents[5],
                                  &n_pts,
                                  &pts_coord,
                                  &pts_g_num);

    // Random permutation
    if (permutation) {
      double *_pts_coord = malloc(sizeof(double) * n_pts * 3);
      int *rand_val = malloc(sizeof(int) * n_pts);
      int *order    = malloc(sizeof(int) * n_pts);
      for (int i = 0; i < n_pts; i++) {
        order[i]    = i;
        rand_val[i] = rand();
      }

      PDM_sort_int(rand_val, order, n_pts);
      free(rand_val);

      for (int i = 0; i < n_pts; i++) {
        memcpy(_pts_coord + 3*i, pts_coord + 3*order[i], sizeof(double) * 3);
        pts_g_num[order[i]] = i+1;
      }
      free(order);
      free(pts_coord);
      pts_coord = _pts_coord;
    }

  }
  else {
    PDM_point_cloud_gen_random(comm,
                               0, // seed
                               0, // geometric_g_num
                               gn_pts,
                               // -radius, -radius, -radius,
                               // radius, radius, radius,
                               g_extents[0], g_extents[1], g_extents[2],
                               g_extents[3], g_extents[4], g_extents[5],
                               &n_pts,
                               &pts_coord,
                               &pts_g_num);
  }


  /* Build point tree */
  int depth_max          = 31;
  // int points_in_leaf_max = 10;
  const double tolerance = 1e-6;
  PDM_point_tree_seq_t *ptree = PDM_point_tree_seq_create(tree_type,
                                                          depth_max,
                                                          points_in_leaf_max,
                                                          tolerance);

  PDM_point_tree_seq_point_cloud_set(ptree,
                                     n_pts,
                                     pts_coord);

  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build(ptree);
  t2 = PDM_MPI_Wtime();
  double t_point_tree = t2 - t1;
  printf("PDM_point_tree_seq_build        : %12.5es\n", t2 - t1);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "point_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(ptree, filename2);

    sprintf(filename2, "points_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename2,
                              n_pts,
                              pts_coord,
                              pts_g_num,
                              NULL);
  }














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
  int   max_tree_depth = 40;
  float max_box_ratio  = 30;

  PDM_box_tree_t *btree = PDM_box_tree_create (max_tree_depth,
                                               max_boxes_leaf,
                                               max_box_ratio);

  t1 = PDM_MPI_Wtime();
  PDM_box_tree_set_boxes (btree,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);
  t2 = PDM_MPI_Wtime();
  double t_box_tree = t2 - t1;
  printf("PDM_box_tree_set_boxes          : %12.5es\n", t2 - t1);
  free(init_location_box);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "box_tree_%i.vtk", i_rank);
    PDM_box_tree_write_vtk(filename2,
                           btree,
                           -1,
                           0);

    sprintf(filename2, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename2,
                        n_box,
                        box_extents,
                        box_g_num);
  }


  t1 = PDM_MPI_Wtime();
  int *box_pts_idx = NULL;
  int *box_pts     = NULL;
  PDM_tree_intersection_point_box(btree,
                                  ptree,
                                  &box_pts_idx,
                                  &box_pts);
  t2 = PDM_MPI_Wtime();
  double t_intersection = t2 - t1;
  printf("PDM_tree_intersection_point_box : %12.5es\n", t2 - t1);
  // PDM_log_trace_connectivity_int(box_pts_idx,
  //                                box_pts,
  //                                n_box,
  //                                "box_pts0 : ");

  if (visu) {
    PDM_g_num_t *box_pts_g_num = malloc(sizeof(PDM_g_num_t) * box_pts_idx[n_box]);
    for (int i = 0; i < box_pts_idx[n_box]; i++) {
      box_pts_g_num[i] = pts_g_num[box_pts[i]];
    }

    // PDM_log_trace_connectivity_long(box_pts_idx,
    //                                 box_pts_g_num,
    //                                 n_box,
    //                                 "box_pts  : ");
    free(box_pts_g_num);
  }


  t1 = PDM_MPI_Wtime();
  int         *box_pts_idx2   = NULL;
  PDM_g_num_t *box_pts_g_num2 = NULL;
  double      *box_pts_coord2 = NULL;
  PDM_box_tree_points_inside_boxes(btree,
                                   n_pts,
                                   pts_g_num,
                                   pts_coord,
                                   &box_pts_idx2,
                                   &box_pts_g_num2,
                                   &box_pts_coord2);
  t2 = PDM_MPI_Wtime();
  double t_old = t2 - t1;
  printf("PDM_box_tree_points_inside_boxes: %12.5es\n", t2 - t1);

  if (visu) {
    for (int i = 0; i < n_box; i++) {
      PDM_inplace_unique_long(box_pts_g_num2,
                              NULL,
                              box_pts_idx2[i],
                              box_pts_idx2[i+1] - 1);
    }

    // PDM_log_trace_connectivity_long(box_pts_idx2,
    //                                 box_pts_g_num2,
    //                                 n_box,
    //                                 "box_pts2 : ");
  }
  free(box_pts_idx2);
  free(box_pts_g_num2);
  free(box_pts_coord2);







  double *box_center = malloc(sizeof(double) * 3 * n_box);
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      box_center[3*i+j] = 0.5*(box_extents[6*i+j] + box_extents[6*i+j+3]);
    }
  }


  PDM_point_tree_seq_t *pbtree = PDM_point_tree_seq_create(tree_type,
                                                           max_tree_depth,
                                                           max_boxes_leaf,
                                                           tolerance);

  PDM_point_tree_seq_point_cloud_set(pbtree,
                                     n_box,
                                     box_center);

  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build(pbtree);

  /* Fix extents */
  for (int i = 0; i < pbtree->n_nodes; i++) {

    // log_trace("Node %d (leaf? %d):\n", i, pbtree->nodes->is_leaf[i]);

    double *e = pbtree->nodes->extents + 6*i;

    for (int j = 0; j < 3; j++) {
      e[j  ] =  HUGE_VAL;
      e[j+3] = -HUGE_VAL;
    }

    for (int k = pbtree->nodes->range[2*i]; k < pbtree->nodes->range[2*i+1]; k++) {
      int ibox = pbtree->new_to_old[k];
      // log_trace("  box %d (%f %f %f  %f %f %f)\n",
      //           ibox,
      //           box_extents[6*ibox  ], box_extents[6*ibox+1], box_extents[6*ibox+2],
      //           box_extents[6*ibox+3], box_extents[6*ibox+4], box_extents[6*ibox+5]);
      for (int j = 0; j < 3; j++) {
        e[j  ] = PDM_MIN(e[j  ], box_extents[6*ibox+j  ]);
        e[j+3] = PDM_MAX(e[j+3], box_extents[6*ibox+j+3]);
      }
    }

    // log_trace("node extents = %f %f %f  %f %f %f\n\n",
    //             e[0], e[1], e[2],
    //             e[3], e[4], e[5]);
  }


  t2 = PDM_MPI_Wtime();
  double t_point_box_tree = t2 - t1;
  printf("PDM_point_tree_seq_build (boxes): %12.5es\n", t2 - t1);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "point_box_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(pbtree, filename2);

    sprintf(filename2, "points_box_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename2,
                              n_box,
                              box_center,
                              box_g_num,
                              NULL);
  }


  int *box_pts_idx3 = NULL;
  int *box_pts3     = NULL;
  t1 = PDM_MPI_Wtime();
  PDM_tree_intersection_point_box2(pbtree,
                                   ptree,
                                   box_extents,
                                   &box_pts_idx3,
                                   &box_pts3);
  t2 = PDM_MPI_Wtime();
  double t_intersection2 = t2 - t1;
  printf("PDM_tree_intersection_point_box2: %12.5es\n", t2 - t1);
  // PDM_log_trace_connectivity_int(box_pts_idx3,
  //                                box_pts3,
  //                                n_box,
  //                                "box_pts3 : ");


  printf("Total intersection : %12.5es\n", t_point_tree + t_box_tree + t_intersection);
  printf("Total old          : %12.5es\n", t_box_tree + t_old);
  printf("Total intersection2: %12.5es\n", t_point_tree + t_point_box_tree + t_intersection2);

  free(box_center);
  PDM_point_tree_seq_free(pbtree);


  if (0) {
    // Check
    int n_err = 0;
    for (int i = 0; i < n_box; i++) {
      int _n_pts  = box_pts_idx [i+1] - box_pts_idx [i];
      int _n_pts3 = box_pts_idx3[i+1] - box_pts_idx3[i];
      if (_n_pts != _n_pts3) {
        log_trace("error box %d:\n", i);
        PDM_log_trace_array_int(box_pts  + box_pts_idx[i],  _n_pts, "0 : ");
        PDM_log_trace_array_int(box_pts3 + box_pts_idx3[i], _n_pts3, "3 : ");
        n_err++;
      }
    }

    printf("n_err = %d\n", n_err);
  }
  free(box_pts_idx3);
  free(box_pts3);
  free(box_pts_idx);
  free(box_pts);


  /* Free */
  PDM_point_tree_seq_free(ptree);
  PDM_box_tree_destroy(&btree);
  PDM_box_set_destroy (&box_set);

  free(pts_coord);
  free(pts_g_num);

  free(box_extents);
  free(box_g_num);
  PDM_MPI_Finalize ();

  return 0;
}
