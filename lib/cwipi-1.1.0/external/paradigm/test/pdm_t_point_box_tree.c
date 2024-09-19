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
#include "pdm_multipart.h"
#include "pdm_point_tree_seq.h"
#include "pdm_point_tree_seq_priv.h"
#include "pdm_box_tree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_array.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_reader_stl.h"
#include "pdm_reader_gamma.h"

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
 int           *tree_type,
 int           *rotation,
 char         **filename,
 int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *tree_type = atoi(argv[i]);
      }
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

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


//https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
static const char *
_get_filename_extension
(
 const char *filename
 )
{
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}


static
void
_split_surface_mesh
(
 const PDM_MPI_Comm        comm,
 const PDM_split_dual_t    part_method,
 const int                 n_part,
       PDM_dmesh_nodal_t  *dmn,
       PDM_multipart_t   **_mpart
)
{
  int n_domain = 1;
  int *n_part_domains = (int *) malloc(sizeof(int) * n_domain);
  n_part_domains[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  free(n_part_domains);

  *_mpart = mpart;
}

static void
_read_and_split_distributed_mesh
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int           *n_vtx,
 int           *n_face,
 int          **face_vtx_idx,
 int          **face_vtx,
 double       **vtx_coord,
 PDM_g_num_t  **face_ln_to_gn,
 PDM_g_num_t  **vtx_ln_to_gn
)
{
  int n_part = 1;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  /* Read distributed mesh */
  PDM_dmesh_nodal_t *dmn = NULL;

  const char *file_extension = _get_filename_extension(filename);

  if (strcmp(file_extension, "stl") == 0) {
    dmn = PDM_reader_stl_dmesh_nodal(comm,
                                     filename);
  }
  else if (strcmp(file_extension, "mesh") == 0) {
    dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                       filename,
                                       0,
                                       0);
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Unknown mesh format %s\n", file_extension);
  }

  PDM_multipart_t *mpart = NULL;
  _split_surface_mesh(comm,
                      part_method,
                      n_part,
                      dmn,
                      &mpart);

  int i_part = 0;
  double *_vtx_coord = NULL;
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            i_part,
                                            &_vtx_coord,
                                            PDM_OWNERSHIP_KEEP);
  *vtx_coord = malloc(sizeof(double) * (*n_vtx) * 3);
  memcpy(*vtx_coord, _vtx_coord, sizeof(double) * (*n_vtx) * 3);

  int *_face_edge_idx = NULL;
  int *_face_edge     = NULL;
  *n_face = PDM_multipart_part_connectivity_get(mpart,
                                                0,
                                                i_part,
                                                PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                &_face_edge_idx,
                                                &_face_edge,
                                                PDM_OWNERSHIP_KEEP);
  *face_vtx_idx = malloc(sizeof(int) * (*n_face + 1));
  memcpy(*face_vtx_idx, _face_edge_idx, sizeof(int) * (*n_face + 1));
  // *face_vtx = malloc(sizeof(int) * _face_vtx_idx[*n_face]);
  // memcpy(*face_vtx, _face_vtx, sizeof(int) * _face_vtx_idx[*n_face]);


  int *_edge_vtx_idx = NULL;
  int *_edge_vtx     = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &_edge_vtx_idx,
                                      &_edge_vtx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_compute_face_vtx_from_face_and_edge(*n_face,
                                          _face_edge_idx,
                                          _face_edge,
                                          _edge_vtx,
                                          face_vtx);

  PDM_g_num_t *_vtx_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  *vtx_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_vtx));
  memcpy(*vtx_ln_to_gn, _vtx_ln_to_gn, sizeof(PDM_g_num_t) * (*n_vtx));

  PDM_g_num_t *_face_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &_face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  *face_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_face));
  memcpy(*face_ln_to_gn, _face_ln_to_gn, sizeof(PDM_g_num_t) * (*n_face));

  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
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
  int                       rotation  = 0;
  char                     *filename  = NULL;
  int                       visu      = 0;

  _read_args(argc,
             argv,
     (int *) &tree_type,
             &rotation,
             &filename,
             &visu);


  /*
   *  Read mesh
   */
  if (i_rank == 0) {
    printf("-- Read mesh\n");
    fflush(stdout);
  }

  int          sm_n_vtx  = 0;
  int          sm_n_face = 0;
  int         *sm_face_vtx_idx  = NULL;
  int         *sm_face_vtx      = NULL;
  double      *sm_vtx_coord     = NULL;
  PDM_g_num_t *sm_face_ln_to_gn = NULL;
  PDM_g_num_t *sm_vtx_ln_to_gn  = NULL;

  _read_and_split_distributed_mesh(comm,
                                   filename,
                                   &sm_n_vtx,
                                   &sm_n_face,
                                   &sm_face_vtx_idx,
                                   &sm_face_vtx,
                                   &sm_vtx_coord,
                                   &sm_face_ln_to_gn,
                                   &sm_vtx_ln_to_gn);

  double *box_extents = malloc(sizeof(double) * sm_n_face * 6);
  double *box_coord   = malloc(sizeof(double) * sm_n_face * 3);
  for (int i = 0; i < sm_n_face; i++) {
    double *e = box_extents + 6*i;

    for (int j = 0; j < 3; j++) {
      e[j  ] =  HUGE_VAL;
      e[j+3] = -HUGE_VAL;
    }

    for (int idx_vtx = sm_face_vtx_idx[i]; idx_vtx < sm_face_vtx_idx[i+1]; idx_vtx++) {
      int vtx_id = sm_face_vtx[idx_vtx] - 1;
      for (int j = 0; j < 3; j++) {
        e[j  ] = PDM_MIN(e[j  ], sm_vtx_coord[3*vtx_id+j]);
        e[j+3] = PDM_MAX(e[j+3], sm_vtx_coord[3*vtx_id+j]);
      }
    }

    for (int j = 0; j < 3; j++) {
      box_coord[3*i+j] = 0.5*(e[j] + e[j+3]);
    }
  }

  /*
   *  Compute global extents of surface mesh
   */
  double l_extents[6] = {
    HUGE_VAL, HUGE_VAL, HUGE_VAL,
    -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
  };

  for (int i = 0; i < sm_n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      double x = sm_vtx_coord[3*i + j];
      l_extents[j]   = PDM_MIN(l_extents[j],   x);
      l_extents[3+j] = PDM_MAX(l_extents[3+j], x);
    }
  }




  double *weight = malloc(sm_n_face * sizeof(double));
  for(int i = 0; i < sm_n_face; ++i) {
    weight[i] = 1.;
  }
  PDM_part_to_block_t *ptb_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                               1.,
                                                               PDM_PART_GEOM_HILBERT,
                                                               &box_coord,
                                                               &sm_face_ln_to_gn,
                                                               &weight,
                                                               &sm_n_face,
                                                               1,
                                                               comm);
  free(box_coord);
  free(weight);

  double *blk_box_extents = NULL;
  PDM_part_to_block_exch(ptb_box,
                         6 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &box_extents,
                         NULL,
               (void **) &blk_box_extents);

  int          blk_n_box     = PDM_part_to_block_n_elt_block_get(ptb_box);
  PDM_g_num_t *blk_box_g_num = PDM_part_to_block_block_gnum_get (ptb_box);

  printf("[%4d] blk_n_box = %d\n", i_rank, blk_n_box);

  PDM_MPI_Barrier(comm);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename2,
                        blk_n_box,
                        blk_box_extents,
                        blk_box_g_num);
  }



  /* Build box_tree */
  double t1, t2;

  int   max_boxes_leaf = 30;
  int   max_tree_depth = 40;
  float max_box_ratio  = 30;

  int *init_location_box = malloc(3 * blk_n_box * sizeof(int));
  for(int i = 0; i < blk_n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  PDM_MPI_Comm comm_alone;
  PDM_MPI_Comm_split(comm, i_rank, 0, &(comm_alone));

  PDM_box_set_t *box_set = PDM_box_set_create(3,
                                              1,
                                              0,
                                              blk_n_box,
                                              blk_box_g_num,
                                              blk_box_extents,
                                              1,
                                              &blk_n_box,
                                              init_location_box,
                                              comm_alone);

  PDM_box_tree_t *btree = PDM_box_tree_create (max_tree_depth,
                                               max_boxes_leaf,
                                               max_box_ratio);

  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  PDM_box_tree_set_boxes (btree,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);
  t2 = PDM_MPI_Wtime();
  printf("[%4d] PDM_box_tree_set_boxes  : %12.5es (%12.5es/box)\n", i_rank, t2 - t1, (t2 - t1) / (double) blk_n_box);
  free(init_location_box);

  if (visu) {
    char filename2[999];
    sprintf(filename2, "box_tree_%i.vtk", i_rank);
    PDM_box_tree_write_vtk(filename2,
                           btree,
                           -1,
                           0);
  }


  /* Build "point-box"_tree */
  double *blk_box_center = malloc(sizeof(double) * 3 * blk_n_box);
  for (int i = 0; i < blk_n_box; i++) {
    for (int j = 0; j < 3; j++) {
      blk_box_center[3*i+j] = 0.5*(blk_box_extents[6*i+j] + blk_box_extents[6*i+j+3]);
    }
  }


  const double tolerance = 1e-6;
  PDM_point_tree_seq_t *pbtree = PDM_point_tree_seq_create(tree_type,
                                                           max_tree_depth,
                                                           max_boxes_leaf,
                                                           tolerance);

  PDM_point_tree_seq_point_cloud_set(pbtree,
                                     blk_n_box,
                                     blk_box_center);

  PDM_MPI_Barrier(comm);
  t1 = PDM_MPI_Wtime();
  PDM_point_tree_seq_build(pbtree);
  free(blk_box_center);

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
      for (int j = 0; j < 3; j++) {
        e[j  ] = PDM_MIN(e[j  ], blk_box_extents[6*ibox+j  ]);
        e[j+3] = PDM_MAX(e[j+3], blk_box_extents[6*ibox+j+3]);
      }
    }
  }


  t2 = PDM_MPI_Wtime();
  printf("[%4d] PDM_point_tree_seq_build: %12.5es (%12.5es/box)\n", i_rank, t2 - t1, (t2 - t1) / (double) blk_n_box);


  if (visu) {
    char filename2[999];
    sprintf(filename2, "point_box_tree_%i.vtk", i_rank);
    PDM_point_tree_seq_write_nodes(pbtree, filename2);
  }

  /* Free memory */
  PDM_point_tree_seq_free(pbtree);
  PDM_box_tree_destroy(&btree);
  PDM_box_set_destroy (&box_set);
  PDM_part_to_block_free(ptb_box);
  free(blk_box_extents);

  PDM_MPI_Comm_free(&comm_alone);


  PDM_MPI_Finalize ();

  return 0;
}
