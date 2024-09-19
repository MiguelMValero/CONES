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
#include "pdm_order.h"
#include "pdm_hilbert.h"


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
 int           *points_in_leaf_max,
 int           *visu,
 int           *rotation
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

    else if (strcmp(argv[i], "-rot") == 0) {
      *rotation = 1;
    }


    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}





inline static int
_intersect_box_box
(
 const int              dim,
 const double *restrict box_extents_a,
 const double *restrict box_extents_b
 )
{
  for (int i = 0; i < dim; i++) {
    if (box_extents_a[i] > box_extents_b[i+dim] || box_extents_b[i] > box_extents_a[i+dim]) {
      return 0;
    }
  }

  return 1;
}

inline static int
_point_inside_box
(
 const int              dim,
 const double *restrict extents,
 const double *restrict coords
 )
{
  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim] || coords[i] < extents[i]) {
      return 0;
    }
  }

  return 1;
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

  PDM_doctree_local_tree_t  tree_type = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  // PDM_doctree_local_tree_t  tree_type = PDM_DOCTREE_LOCAL_TREE_KDTREE;

  double radius = 10.;

  PDM_g_num_t gn_pts    = 10;
  PDM_g_num_t gn_box    = 10;
  int         visu      = 0;
  int         points_in_leaf_max = 10;
  int         rotation           = 0;

  _read_args(argc,
             argv,
             &gn_pts,
             &gn_box,
             &points_in_leaf_max,
             &visu,
             &rotation);

  /* (Random) boxes */
  int n_box = 0;
  double      *box_extents = NULL;
  PDM_g_num_t *box_g_num   = NULL;

  double _n = PDM_MAX(2, 1 + pow(gn_box, 1./3.));

  int n_vtx_x = (int) _n;
  int n_vtx_y = (int) _n;
  int n_vtx_z = (int) _n;
  PDM_box_gen_cartesian(comm,
                        n_vtx_x,
                        n_vtx_y,
                        n_vtx_z,
                        -radius, -radius, -radius,
                        radius, radius, radius,
                        &n_box,
                        &box_extents,
                        &box_g_num);

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

  /* Build "point-box"_tree */
  double *blk_box_center = malloc(sizeof(double) * 3 * n_box);
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      blk_box_center[3*i+j] = 0.5*(box_extents[6*i+j] + box_extents[6*i+j+3]);
    }
  }


  /* Random point cloud */
  int          n_pts     = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;

  PDM_point_cloud_gen_random(comm,
                               0, // seed
                               0, // geometric_g_num
                               gn_pts,
                               -2*radius, -2*radius, -2*radius,
                               2*radius, 2*radius, 2*radius,
                               &n_pts,
                               &pts_coord,
                               &pts_g_num);


  double t1, t2, t3;

  int reorder = 1;
  if(reorder == 1) {
    PDM_MPI_Barrier(comm);
    t1 = PDM_MPI_Wtime();
    int dim = 3;
    double extents[2*dim]; /** DIM x 2**/
    PDM_hilbert_code_t *hilbert_box_codes     = (PDM_hilbert_code_t *) malloc (n_box * sizeof(PDM_hilbert_code_t));
    int *box_order     = (int *) malloc (n_box * sizeof(int));
    PDM_hilbert_get_coord_extents_par(dim, n_pts, pts_coord, extents, comm);
    // PDM_hilbert_get_coord_extents_par(dim, n_box, blk_box_center, extents, comm);
    PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, n_box, blk_box_center, hilbert_box_codes);
    PDM_hilbert_local_order(n_box, hilbert_box_codes, box_order);
    PDM_order_array(n_box, 6 * sizeof(double), box_order, box_extents);
    PDM_order_array(n_box, sizeof(PDM_g_num_t), box_order, box_g_num);

    for(int i = 0; i < n_box; ++i) {
      box_g_num[i] = i + 1;
    }

    t2 = PDM_MPI_Wtime();
    free(hilbert_box_codes);
    free(box_order);

    printf("Hilbert box times : %12.5e \n", t2 - t1);

    PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (n_pts * sizeof(PDM_hilbert_code_t));
    int *order     = (int *) malloc (n_pts * sizeof(int));
    PDM_hilbert_get_coord_extents_par(dim, n_pts, pts_coord, extents, comm);
    PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, n_pts, pts_coord, hilbert_codes);
    PDM_hilbert_local_order(n_pts, hilbert_codes, order);
    PDM_order_array(n_pts, 3 * sizeof(double)     , order, pts_coord);
    PDM_order_array(n_pts,     sizeof(PDM_g_num_t), order, pts_g_num);
    t2 = PDM_MPI_Wtime();
    printf("Hilbert times : %12.5e \n", t2 - t1);
    free(hilbert_codes);
    free(order);

    for(int i = 0; i < n_pts; ++i) {
      pts_g_num[i] = i + 1;
    }

    t1 = PDM_MPI_Wtime();
    int n_part = n_pts/24;
    PDM_hilbert_code_t *rank_index = malloc((n_part+1) * sizeof(PDM_hilbert_code_t));

    double *expli_box_extents = malloc(3 * 8 * n_box * sizeof(double));
    PDM_hilbert_code_t *expli_box_codes   = malloc( 8 * n_box * sizeof(PDM_hilbert_code_t));

    double *weight     = (double *) malloc (8 * n_box * sizeof(double));

    for(int i = 0; i < 8 * n_box; ++i) {
      weight   [i] = 1.;
    }


    // Explicit box pts
    for(int i = 0; i < n_box; ++i) {

      double *_expli_box_extents = expli_box_extents + 3*8*i;

      double x1 = box_extents[6*i  ];
      double y1 = box_extents[6*i+1];
      double z1 = box_extents[6*i+2];

      double dx = box_extents[6*i+3] - box_extents[6*i  ];
      double dy = box_extents[6*i+4] - box_extents[6*i+1];
      double dz = box_extents[6*i+5] - box_extents[6*i+2];

      int j = 0;
      _expli_box_extents[3*j  ] = x1;
      _expli_box_extents[3*j+1] = y1;
      _expli_box_extents[3*j+2] = z1;

      j = 1;
      _expli_box_extents[3*j  ] = x1;
      _expli_box_extents[3*j+1] = y1 + dy;
      _expli_box_extents[3*j+2] = z1;

      j = 2;
      _expli_box_extents[3*j  ] = x1 + dx;
      _expli_box_extents[3*j+1] = y1;
      _expli_box_extents[3*j+2] = z1;

      j = 3;
      _expli_box_extents[3*j  ] = x1 + dx;
      _expli_box_extents[3*j+1] = y1 + dy;
      _expli_box_extents[3*j+2] = z1;

      j = 4;
      _expli_box_extents[3*j  ] = x1;
      _expli_box_extents[3*j+1] = y1;
      _expli_box_extents[3*j+2] = z1 + dz;

      j = 5;
      _expli_box_extents[3*j  ] = x1;
      _expli_box_extents[3*j+1] = y1 + dy;
      _expli_box_extents[3*j+2] = z1 + dz;

      j = 6;
      _expli_box_extents[3*j  ] = x1 + dx;
      _expli_box_extents[3*j+1] = y1;
      _expli_box_extents[3*j+2] = z1 + dz;

      j = 7;
      _expli_box_extents[3*j  ] = x1 + dx;
      _expli_box_extents[3*j+1] = y1 + dy;
      _expli_box_extents[3*j+2] = z1 + dz;

    }

    // Faire une fonction dedié : PDM_hilbert_build_boxes_rank_index pour eviter de copier
    PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, 8 * n_box, expli_box_extents, expli_box_codes);
    PDM_hilbert_build_rank_index(dim, n_part, 8 * n_box, expli_box_codes, weight, NULL, rank_index, comm);
    t2 = PDM_MPI_Wtime();

    printf("PDM_hilbert_build_rank_index : %12.5e \n", t2 - t1);
    // PDM_log_trace_array_double(rank_index, n_part, "rank_index :: ");

    free(rank_index);
    free(expli_box_codes);
    free(expli_box_extents);
    free(blk_box_center);
    free(weight);
  }


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



  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build_from_boxes(ptree,
                                      n_box,
                                      box_extents);
  t2 = PDM_MPI_Wtime();

  // PDM_log_trace_connectivity_int(ptree->leaf_box_idx, ptree->leaf_box_ids, ptree->n_leaf, "leaf_box :");

  // int *box_pts_idx = malloc(sizeof(int) * (n_box + 1));
  // box_pts_idx[0] = 0;
  // int s_box_pts = 4*n_box;
  // int *box_pts = malloc(sizeof(int) * s_box_pts);

  // for (int ibox = 0; ibox < n_box; ibox++) {
  //   box_pts_idx[ibox+1] = box_pts_idx[ibox];

  //   double *be = box_extents + 6*ibox;

  //   for (int ileaf = 0; ileaf < ptree->n_leaf; ileaf++) {
  //     int leaf_id = ptree->leaf_ids[ileaf];

  //     double *le = ptree->nodes->extents + 6*leaf_id;

  //     if (_intersect_box_box(3, be, le)) {

  //       for (int ipt = ptree->nodes->range[2*leaf_id]; ipt < ptree->nodes->range[2*leaf_id+1]; ipt++) {

  //         int point_id = ptree->new_to_old[ipt];

  //         double *pc = ptree->_pts_coord + 3*ipt;

  //         if (_point_inside_box(3, be, pc)) {

  //           if (box_pts_idx[ibox+1] >= s_box_pts) {
  //             s_box_pts *= 2;
  //             box_pts = realloc(box_pts, sizeof(int) * s_box_pts);
  //           }
  //           box_pts[box_pts_idx[ibox+1]++] = point_id;
  //         }
  //       }
  //     }
  //   }
  // }

  // Tentative Bruno
  // Plus simple dans l'autre sens nan ?
  int *pts_box_idx = malloc(sizeof(int) * (n_pts + 1));
  int *pts_box     = malloc(sizeof(int) * (8 * n_pts + 1));
  for(int i = 0; i < n_pts+1; ++i)  {
    pts_box_idx[i] = 0;
  }

  int idx_write_pts = 0;
  for (int ileaf = 0; ileaf < ptree->n_leaf; ileaf++) {
    int leaf_id = ptree->leaf_ids[ileaf];
    for (int ipt = ptree->nodes->range[2*leaf_id]; ipt < ptree->nodes->range[2*leaf_id+1]; ipt++) {
      // int point_id = ptree->new_to_old[ipt];
      double *pc = ptree->_pts_coord + 3*ipt;

      /* Brute force because already pre-conditionned */
      for(int idx_box = ptree->leaf_box_idx[ileaf]; idx_box < ptree->leaf_box_idx[ileaf+1]; ++idx_box) {
        int boxe_id = ptree->leaf_box_ids[idx_box];
        double *be = box_extents + 6*boxe_id;

        if (_point_inside_box(3, be, pc)) {
          pts_box[idx_write_pts++] = boxe_id;
          pts_box_idx[ipt+1]++;
        }
      }
    }
  }


  // Accumulate
  for(int i = 0; i < n_pts; ++i)  {
    pts_box_idx[i+1] += pts_box_idx[i];
  }

  // // PDM_log_trace_connectivity_int(pts_box_idx, pts_box, n_pts, "pts_box :: ");

  // Revert
  // int *box_pts_idx = NULL;
  // int *box_pts     = NULL;
  int *box_pts_idx = malloc(sizeof(int) * (n_box + 1));
  int *box_pts_n   = malloc(sizeof(int) * (n_box + 1));
  for(int i = 0; i < n_box; ++i)  {
    box_pts_n[i] = 0;
  }

  for(int i = 0; i < pts_box_idx[n_pts]; ++i){
    box_pts_n[pts_box[i]]++;
  }

  box_pts_idx[0] = 0;
  for(int i = 0; i < n_box; ++i)  {
    box_pts_idx[i+1] = box_pts_idx[i] + box_pts_n[i];
    box_pts_n[i] = 0;
  }
  // PDM_log_trace_array_int(box_pts_idx, n_box, "box_pts_idx : ");
  // Fill
  int *box_pts   = malloc(sizeof(int) * (box_pts_idx[n_box]));
  for(int ipt = 0; ipt < n_pts; ++ipt) {
    for(int idx_box = pts_box_idx[ipt]; idx_box < pts_box_idx[ipt+1]; ++idx_box) {
      int box_id = pts_box[idx_box];
      int idx_write = box_pts_idx[box_id] + box_pts_n[box_id]++;
      // printf("idx_write = %i \n", idx_write);
      box_pts[idx_write] = ptree->new_to_old[ipt];
    }
  }
  free(pts_box_idx);
  free(pts_box    );
  free(box_pts_n    );

  t3 = PDM_MPI_Wtime();
  printf("new : build %12.5es, sollicitate %12.5es, total %12.5es\n", t2-t1, t3-t2, t3-t1);

  // PDM_log_trace_connectivity_int(box_pts_idx,
  //                                box_pts,
  //                                n_box,
  //                                "box_pts : ");


  PDM_point_tree_seq_t *ptree2 = PDM_point_tree_seq_create(tree_type,
                                                           depth_max,
                                                           points_in_leaf_max,
                                                           tolerance);

  PDM_point_tree_seq_point_cloud_set(ptree2,
                                     n_pts,
                                     pts_coord);

  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build(ptree2);
  t2 = PDM_MPI_Wtime();

  int *box_pts_idx2 = NULL;
  int *box_pts2     = NULL;
  PDM_point_tree_seq_points_inside_boxes(ptree2,
                                         n_box,
                                         box_extents,
                                         &box_pts_idx2,
                                         &box_pts2);

  t3 = PDM_MPI_Wtime();
  printf("old : build %12.5es, sollicitate %12.5es, total %12.5es\n", t2-t1, t3-t2, t3-t1);



  if (0) {
    // Check
    int n_err = 0;
    for (int i = 0; i < n_box; i++) {
      int _n_pts  = box_pts_idx [i+1] - box_pts_idx [i];
      int _n_pts2 = box_pts_idx2[i+1] - box_pts_idx2[i];
      if (_n_pts != _n_pts2) {
        log_trace("error box %d:\n", i);
        PDM_log_trace_array_int(box_pts  + box_pts_idx[i],  _n_pts, "new : ");
        PDM_log_trace_array_int(box_pts2 + box_pts_idx2[i], _n_pts2, "old : ");
        n_err++;
      }
    }

    printf("n_err = %d\n", n_err);
  }
  free(box_pts_idx2);
  free(box_pts2);
  free(box_pts_idx);
  free(box_pts);

  PDM_point_tree_seq_free(ptree);
  PDM_point_tree_seq_free(ptree2);



  if (visu) {
    char filename2[999];
    sprintf(filename2, "point_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(ptree, filename2);

    sprintf(filename2, "point_tree2_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(ptree2, filename2);

    sprintf(filename2, "points_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename2,
                              n_pts,
                              pts_coord,
                              pts_g_num,
                              NULL);

    sprintf(filename2, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename2,
                        n_box,
                        box_extents,
                        box_g_num);
  }

  free(pts_coord);
  free(pts_g_num);
  free(box_extents);
  free(box_g_num);


  PDM_MPI_Finalize ();

  return 0;
}
