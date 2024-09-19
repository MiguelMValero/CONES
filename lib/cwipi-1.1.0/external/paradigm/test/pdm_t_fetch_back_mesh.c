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
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_dbbtree.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_multipart.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_plane.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_triangle.h"
#include "pdm_part_to_part.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"
#include "pdm_para_octree.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const int verbose = 0;
static const int vtk     = 0;

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
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_back,
 PDM_g_num_t   *n_work,
 int           *method
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nb") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *n_back = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-nw") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *n_work = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *method = atoi(argv[i]);
      }
    }

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


static void
_rotate
(
 const int     n_pts,
       double *coord
 )
{
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}


static void
_deformation
(
 const double  length,
 const int     n_vtx,
       double *vtx_coord
 )
{
  double amplitude = 0.15;
  double frequency = 4.;

  for (int i = 0; i < n_vtx; i++) {
    double x = (vtx_coord[3*i    ] - 0.5) / length;
    double y = (vtx_coord[3*i + 1] - 0.5) / length;
    double z = (vtx_coord[3*i + 2] - 0.5) / length;

    vtx_coord[3*i    ] += amplitude*length*cos(frequency*y);
    vtx_coord[3*i + 1] += amplitude*length*cos(frequency*z);
    vtx_coord[3*i + 2] += amplitude*length*cos(frequency*x);
  }
}



static
void
_fetch_back_mesh_partial_dist_cloud_surf
(
 PDM_MPI_Comm   comm,
 int            pback_n_vtx,
 double        *pback_vtx_coord,
 PDM_g_num_t   *pback_vtx_ln_to_gn,
 // int            pback_n_face,
 // int           *pback_face_vtx_idx,
 // int           *pback_face_vtx,
 // PDM_g_num_t   *pback_face_ln_to_gn,
 PDM_dbbtree_t *dbbt,
 // int            pwork_n_vtx,
 double        *pwork_vtx_coord,
 // PDM_g_num_t   *pwork_vtx_ln_to_gn,
 int            pwork_n_edge,
 int           *pwork_edge_vtx,
 PDM_g_num_t   *pwork_edge_ln_to_gn,
 int          **pwork_edge_back_face_idx,
 PDM_g_num_t  **pwork_edge_back_face
)
{
  /*
   *  For each point inserted on a work edge, find its closest back vtx --> upper bound on distance
   */
  const int depth_max = 31;
  const int points_in_leaf_max = 4;
  PDM_para_octree_t *octree = PDM_para_octree_create(1,
                                                     depth_max,
                                                     points_in_leaf_max,
                                                     0,
                                                     comm);

  PDM_para_octree_point_cloud_set(octree,
                                  0,
                                  pback_n_vtx,
                                  pback_vtx_coord,
                                  pback_vtx_ln_to_gn);

  PDM_para_octree_build(octree, NULL);


  // -->> take this as an input (riemannian midpoint...)
  double *inserted_pt_coord = malloc(sizeof(double) * pwork_n_edge * 3);
  double *pwork_edge_length = malloc(sizeof(double) * pwork_n_edge);
  for (int edge_id = 0; edge_id < pwork_n_edge; edge_id++) {
    int vtx_id1 = pwork_edge_vtx[2*edge_id  ] - 1;
    int vtx_id2 = pwork_edge_vtx[2*edge_id+1] - 1;

    pwork_edge_length[edge_id] = 0.;
    for (int j = 0; j < 3; j++) {
      inserted_pt_coord[3*edge_id+j] = 0.5*(pwork_vtx_coord[3*vtx_id1+j] + pwork_vtx_coord[3*vtx_id2+j]);
      double x = pwork_vtx_coord[3*vtx_id1+j] - pwork_vtx_coord[3*vtx_id2+j];
      pwork_edge_length[edge_id] += x*x;
    }
  }
  // <<--

  PDM_g_num_t *closest_back_vtx_gnum  = malloc(sizeof(PDM_g_num_t) * pwork_n_edge);
  double      *closest_back_vtx_dist2 = malloc(sizeof(double)      * pwork_n_edge);
  PDM_para_octree_single_closest_point(octree,
                                       pwork_n_edge,
                                       inserted_pt_coord,
                                       pwork_edge_ln_to_gn,
                                       closest_back_vtx_gnum,
                                       closest_back_vtx_dist2);
  PDM_para_octree_free(octree);
  free(closest_back_vtx_gnum); // unused

  /* Threshold upper bound distance */
  if (0) {
    for (int edge_id = 0; edge_id < pwork_n_edge; edge_id++) {
      closest_back_vtx_dist2[edge_id] = PDM_MAX(pwork_edge_length[edge_id],
                                                closest_back_vtx_dist2[edge_id]);
    }
  }
  free(pwork_edge_length);


  /*
   *  For each point inserted on a work edge, find all back elements closer than the upper bound distance
   */
  PDM_dbbtree_closest_upper_bound_dist_boxes_get(dbbt,
                                                 pwork_n_edge,
                                                 inserted_pt_coord,
                                                 pwork_edge_ln_to_gn,
                                                 closest_back_vtx_dist2,
                                                 pwork_edge_back_face_idx,
                                                 pwork_edge_back_face);
  free(inserted_pt_coord);
  free(closest_back_vtx_dist2);
}




static
void
_fetch_back_mesh_Bruno
(
 PDM_MPI_Comm   comm,
 int            pback_n_vtx,
 double        *pback_vtx_coord,
 PDM_g_num_t   *pback_vtx_ln_to_gn,
 int            pback_n_face,
 int           *pback_face_vtx_idx,
 int           *pback_face_vtx,
 PDM_g_num_t   *pback_face_ln_to_gn,
 double        *pback_face_extents,
 int            pwork_n_vtx,
 double        *pwork_vtx_coord,
 PDM_g_num_t   *pwork_vtx_ln_to_gn,
 int            pwork_n_edge,
 int           *pwork_edge_vtx,
 PDM_g_num_t   *pwork_edge_ln_to_gn,
 double        *g_extents,
 int          **line_to_back_idx,
 PDM_g_num_t  **line_to_back
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_part = 1;
  const double eps_extents = 1.0e-6;

  double *line_coord = malloc(6 * pwork_n_edge * sizeof(double));

  for(int i_edge = 0; i_edge < pwork_n_edge; ++i_edge) {
    int i_vtx1 = PDM_ABS(pwork_edge_vtx[2*i_edge  ])-1;
    int i_vtx2 = PDM_ABS(pwork_edge_vtx[2*i_edge+1])-1;

    line_coord[6*i_edge  ] = pwork_vtx_coord[3*i_vtx1  ];
    line_coord[6*i_edge+1] = pwork_vtx_coord[3*i_vtx1+1];
    line_coord[6*i_edge+2] = pwork_vtx_coord[3*i_vtx1+2];

    line_coord[6*i_edge+3] = pwork_vtx_coord[3*i_vtx2  ];
    line_coord[6*i_edge+4] = pwork_vtx_coord[3*i_vtx2+1];
    line_coord[6*i_edge+5] = pwork_vtx_coord[3*i_vtx2+2];
  }

  // int*         line_to_back2_idx = NULL;
  // PDM_g_num_t* line_to_back2     = NULL;
  // PDM_dbbtree_lines_last_intersect_boxes(dbbt,
  //                                        pwork_n_edge,
  //                                        pwork_edge_ln_to_gn,
  //                                        line_coord,
  //                                        &line_to_back2_idx,
  //                                        &line_to_back2);


  /*
   *  For each back vtx, find closest edge of work mesh
   */
  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 1);

  int *pwork_edge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, pwork_n_edge);

  PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                          0,//i_part,
                                          pwork_n_edge,
                                          pwork_edge_vtx_idx,
                                          pwork_edge_vtx,
                                          pwork_edge_ln_to_gn,
                                          pwork_n_vtx,
                                          pwork_vtx_coord,
                                          pwork_vtx_ln_to_gn);

  PDM_dist_cloud_surf_n_part_cloud_set (dist, 0, n_part);

  PDM_dist_cloud_surf_cloud_set (dist,
                                 0,
                                 0,//i_part,
                                 pback_n_vtx,
                                 pback_vtx_coord,
                                 pback_vtx_ln_to_gn);

  PDM_dist_cloud_surf_compute (dist);

  double      *distance;
  double      *projected;
  PDM_g_num_t *closest_elt_gnum;
  PDM_dist_cloud_surf_get (dist,
                           0,
                           0,//i_part,
                           &distance,
                           &projected,
                           &closest_elt_gnum);

  if (vtk) {

    char filename[999];

    sprintf(filename, "dist_ridge_line_%3.3d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         pwork_n_edge,
                         line_coord,
                         pwork_edge_ln_to_gn,
                         NULL);
    // free (line_coord);



    sprintf(filename, "dist_ridge_pts_%3.3d.vtk", i_rank);
    PDM_vtk_write_point_cloud (filename,
                               pback_n_vtx,
                               pback_vtx_coord,
                               pback_vtx_ln_to_gn,
                               NULL);


    double *proj_line_coord = malloc (sizeof(double) * pback_n_vtx * 6);
    int idx = 0;
    for (int i = 0; i < pback_n_vtx; i++) {

      for (int k = 0; k < 3; k++) {
        proj_line_coord[idx++] = pback_vtx_coord[3*i + k];
      }

      for (int k = 0; k < 3; k++) {
       proj_line_coord[idx++] = projected[3*i + k];
      }
    }

    sprintf(filename, "dist_ridge_proj_%3.3d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         pback_n_vtx,
                         proj_line_coord,
                         pback_vtx_ln_to_gn,
                         NULL);
    free (proj_line_coord);

  }


  /*
   *  On souhaite gonflé les boites du maillage back -> Revient a envoyé au maillage back le projected
   *   On envoie l'arrete + ses coordonnnées pour gonflé
   *     closest_elt_gnum = vtx_back->working_edge
   *     part1 = vtx_back
   *     part2 = edge
   */
  int *closest_elt_gnum_idx = malloc((pback_n_vtx+1) * sizeof(int));
  for(int i = 0; i < pback_n_vtx+1; ++i) {
    closest_elt_gnum_idx[i] = i;
  }
  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) &pback_vtx_ln_to_gn,
                                                                           &pback_n_vtx,
                                                                           1,
                                                    (const PDM_g_num_t **) &pwork_edge_ln_to_gn,
                                                                           &pwork_n_edge,
                                                                           1,
                                                    (const int         **) &closest_elt_gnum_idx,
                                                    (const PDM_g_num_t **) &closest_elt_gnum,
                                                                           comm);

  free(closest_elt_gnum_idx);


  // int    **part_stride = NULL;
  double **tmp_part_data = NULL;
  int request_return_selected;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 6,
                                 sizeof(double),
                 (const int  **) NULL,
                 (const void **) &line_coord,
                                 NULL,
                      (void ***) &tmp_part_data,
                                 &request_return_selected);


  PDM_part_to_part_reverse_iexch_wait(ptp, request_return_selected);
  double *recv_line_coord = tmp_part_data[0];
  free(tmp_part_data);

  if (verbose) {
    // PDM_log_trace_array_double(recv_line_coord, 6*pback_n_vtx, "recv_line_coord ::");
    PDM_log_trace_array_long(closest_elt_gnum, pback_n_vtx, "closest_elt_gnum :: ");

    int          *n_ref_entity1     = NULL;
    int         **ref_l_num_entity1 = NULL;
    PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_entity1, &ref_l_num_entity1);

    PDM_log_trace_array_int(ref_l_num_entity1[0], n_ref_entity1[0], "ref_l_num_entity1 :: ");
  }

  if(vtk) {
    char filename[999];
    sprintf(filename, "inflate_recv_line_%3.3d.vtk", i_rank);
    // PDM_vtk_write_lines (filename,
    //                      pback_n_vtx,
    //                      recv_line_coord,
    //                      closest_elt_gnum,
    //                      NULL);
    PDM_vtk_write_lines (filename,
                         pback_n_vtx,
                         recv_line_coord,
                         pback_vtx_ln_to_gn,
                         NULL);
  }

  /*
   * Correction extents
   */
  for (int iface = 0; iface < pback_n_face; iface++) {

    double *tmp_extents = pback_face_extents + 6*iface;
    // for (int k = 0; k < 3; k++) {
    //   tmp_extents[k]     =  HUGE_VAL;
    //   tmp_extents[3 + k] = -HUGE_VAL;
    // }

    for (int ivtx = pback_face_vtx_idx[iface]; ivtx < pback_face_vtx_idx[iface+1]; ivtx++) {
      int vtx_id = pback_face_vtx[ivtx]-1;
      // log_trace("vtx_id = %i \n", vtx_id);
      double *tmp_coord = recv_line_coord + 6*vtx_id;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }

      tmp_coord = recv_line_coord + 6*vtx_id + 3;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }
      // double *tmp_coord = recv_line_coord + 3*vtx_id;
      // for (int k = 0; k < 3; k++) {
      //   tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
      //   tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      // }
    } // end loop on vertices of iface

    // Threshold/Inflate extents
    for (int k = 0; k < 3; k++) {
      double size = tmp_extents[3 + k] - tmp_extents[k];

      // tmp_extents[k]     -= size*eps_extents;
      // tmp_extents[3 + k] += size*eps_extents;

      if (size < eps_extents) {
        tmp_extents[k]     -= 0.5*eps_extents;
        tmp_extents[3 + k] += 0.5*eps_extents;
      }
    }

  } // end loop on background faces

  if (vtk) {
    char filename[999];

    sprintf(filename, "back_faces_inflate_%d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pback_n_vtx,
                           pback_vtx_coord,
                           pback_vtx_ln_to_gn,
                           pback_n_face,
                           pback_face_vtx_idx,
                           pback_face_vtx,
                           pback_face_ln_to_gn,
                           NULL);

    sprintf(filename, "back_face_inflate_extents_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        pback_n_face,
                        pback_face_extents,
                        pback_face_ln_to_gn);
  }


  const int dim = 3;
  PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                 n_part,
                                                 &pback_n_face,
                               (const double **) &pback_face_extents,
                          (const PDM_g_num_t **) &pback_face_ln_to_gn);


  PDM_dbbtree_lines_intersect_boxes(dbbt,
                                    pwork_n_edge,
                                    pwork_edge_ln_to_gn,
                                    line_coord,
                                    line_to_back_idx,
                                    line_to_back);

  free(recv_line_coord);


  PDM_part_to_part_free(ptp);

  PDM_dist_cloud_surf_dump_times(dist);
  PDM_dist_cloud_surf_free (dist);
  free(line_coord);
  free(pwork_edge_vtx_idx);

  PDM_dbbtree_free(dbbt);
  PDM_box_set_destroy(&box_set);
}






/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int n_part = 1;
  PDM_g_num_t n_back = 10;
  PDM_g_num_t n_work = 1;
  int deform         = 1;
  int method         = 0;
 _read_args (argc,
             argv,
             &n_back,
             &n_work,
             &method);

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  const double x_center = 0.1;
  const double y_center = 0.2;
  const double z_center = 0.3;
  const double radius   = 3.14;


  /*
   *  Generate a fine background distributed mesh (single surface group)
   */
  double      *dback_vtx_coord    = NULL;
  int         *dback_face_vtx_idx = NULL;
  PDM_g_num_t *dback_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx   = NULL;
  PDM_g_num_t *back_distrib_face  = NULL;

  // PDM_sphere_surf_gen(comm,
  //                     2*n_back,
  //                     n_back,
  //                     x_center,
  //                     y_center,
  //                     z_center,
  //                     radius,
  //                     &dback_vtx_coord,
  //                     &dback_face_vtx_idx,
  //                     &dback_face_vtx,
  //                     &back_distrib_vtx,
  //                     &back_distrib_face);

  PDM_sphere_surf_icosphere_gen(comm,
                                n_back,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dback_vtx_coord,
                                &dback_face_vtx_idx,
                                &dback_face_vtx,
                                &back_distrib_vtx,
                                &back_distrib_face);

  if (deform) {
    int dback_n_vtx = back_distrib_vtx[i_rank+1] - back_distrib_vtx[i_rank];
    _deformation(radius,
                 dback_n_vtx,
                 dback_vtx_coord);
  }

  /*
   *  Build the back mesh dbbtree
   */

  /* Assemble partitions from block-distribution of faces */
  int dback_n_face = back_distrib_face[i_rank+1] - back_distrib_face[i_rank];

  PDM_g_num_t *dback_face_ln_to_gn = malloc(dback_n_face * sizeof(PDM_g_num_t));
  for (int i = 0; i < dback_n_face; ++i) {
    dback_face_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }


  PDM_g_num_t *pback_vtx_ln_to_gn = NULL;
  int         *pback_face_vtx_idx = NULL;
  int         *pback_face_vtx     = NULL;
  int          pback_n_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           dback_face_vtx_idx,
                                                           dback_face_vtx,
                                                           dback_n_face,
                                  (const PDM_g_num_t *)    dback_face_ln_to_gn,
                                                           &pback_n_vtx,
                                                           &pback_vtx_ln_to_gn,
                                                           &pback_face_vtx_idx,
                                                           &pback_face_vtx);

  double **tmp_pback_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        dback_vtx_coord,
                                        &pback_n_vtx,
                 (const PDM_g_num_t **) &pback_vtx_ln_to_gn,
                                        &tmp_pback_vtx_coord);
  double *pback_vtx_coord = tmp_pback_vtx_coord[0];
  free(tmp_pback_vtx_coord);



  /* Compute the bounding boxes of local faces */
  double *back_face_extents = malloc(sizeof(double) * dback_n_face * 6);
  const double eps_extents = 1.0e-6;
  for (int iface = 0; iface < dback_n_face; iface++) {

    double *tmp_extents = back_face_extents + 6*iface;
    for (int k = 0; k < 3; k++) {
      tmp_extents[k]     =  HUGE_VAL;
      tmp_extents[3 + k] = -HUGE_VAL;
    }

    for (int ivtx = pback_face_vtx_idx[iface]; ivtx < pback_face_vtx_idx[iface+1]; ivtx++) {
      int vtx_id = pback_face_vtx[ivtx]-1;
      double *tmp_coord = pback_vtx_coord + 3*vtx_id;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }
    } // end loop on vertices of iface

    // Threshold/Inflate extents
    for (int k = 0; k < 3; k++) {
      double size = tmp_extents[3 + k] - tmp_extents[k];

      // tmp_extents[k]     -= size*eps_extents;
      // tmp_extents[3 + k] += size*eps_extents;

      if (size < eps_extents) {
        tmp_extents[k]     -= 0.5*eps_extents;
        tmp_extents[3 + k] += 0.5*eps_extents;
      }
    }

  } // end loop on background faces


  if (vtk) {
    char filename[999];

    sprintf(filename, "back_faces_%d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pback_n_vtx,
                           pback_vtx_coord,
                           pback_vtx_ln_to_gn,
                           dback_n_face,
                           pback_face_vtx_idx,
                           pback_face_vtx,
                           dback_face_ln_to_gn,
                           NULL);

    sprintf(filename, "back_face_extents_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        dback_n_face,
                        back_face_extents,
                        dback_face_ln_to_gn);
  }


  /* Create and build dbbtree */
  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < dback_n_face; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN(l_extents[k],     back_face_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX(l_extents[k + 3], back_face_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
  }
  // Inflate and break symmetry
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }


  /*
   *  Generate a coarse 'work' partitioned mesh
   */

  /* Generate block-distributed nodal mesh */
  PDM_dmesh_nodal_t *dmn = NULL;

  // PDM_sphere_surf_gen_nodal(comm,
  //                           2*n_work,
  //                           n_work,
  //                           x_center,
  //                           y_center,
  //                           z_center,
  //                           radius,
  //                           &dmn);
  PDM_sphere_surf_icosphere_gen_nodal(comm,
                                      n_work,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn);

  /* Split the mesh */
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int n_domain                 = 1;
  int *n_part_domains          = (int *) malloc(sizeof(int) * n_domain);
  n_part_domains[0]            = n_part;

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



  // Vertices
  int     i_part          = 0;
  double *pwork_vtx_coord = NULL;
  int pwork_n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &pwork_vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
  if (1) {
    _rotate(pwork_n_vtx,
            pwork_vtx_coord);
  }

  if (deform) {
    _deformation(radius,
                 pwork_n_vtx,
                 pwork_vtx_coord);
  }

  PDM_g_num_t *pwork_vtx_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &pwork_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Edges
  int *pwork_edge_vtx     = NULL;
  int *pwork_edge_vtx_idx = NULL;
  int  pwork_n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                          0,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                          &pwork_edge_vtx_idx,
                                                          &pwork_edge_vtx,
                                                          PDM_OWNERSHIP_KEEP);

  if (pwork_edge_vtx_idx != NULL) free(pwork_edge_vtx_idx);

  PDM_g_num_t *pwork_edge_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &pwork_edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);


  if (vtk) {
    char filename[999];

    sprintf(filename, "work_edges_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               pwork_n_vtx,
                               pwork_vtx_coord,
                               pwork_vtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               pwork_n_edge,
                               pwork_edge_vtx,
                               pwork_edge_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }





  int         *line_to_back_idx = NULL;
  PDM_g_num_t *line_to_back     = NULL;
  switch (method) {
    case 0: {

      _fetch_back_mesh_Bruno(comm,
                             pback_n_vtx,
                             pback_vtx_coord,
                             pback_vtx_ln_to_gn,
                             dback_n_face,
                             pback_face_vtx_idx,
                             pback_face_vtx,
                             dback_face_ln_to_gn,
                             back_face_extents,
                             pwork_n_vtx,
                             pwork_vtx_coord,
                             pwork_vtx_ln_to_gn,
                             pwork_n_edge,
                             pwork_edge_vtx,
                             pwork_edge_ln_to_gn,
                             g_extents,
                             &line_to_back_idx,
                             &line_to_back);

      break;
    }
    case 1: {
      PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

      PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                     n_part,
                                                     &dback_n_face,
                              (const double      **) &back_face_extents,
                              (const PDM_g_num_t **) &dback_face_ln_to_gn);

      _fetch_back_mesh_partial_dist_cloud_surf(comm,
                                               pback_n_vtx,
                                               pback_vtx_coord,
                                               pback_vtx_ln_to_gn,
                                               // dback_n_face,
                                               // pback_face_vtx_idx,
                                               // pback_face_vtx,
                                               // dback_face_ln_to_gn,
                                               dbbt,
                                               // pwork_n_vtx,
                                               pwork_vtx_coord,
                                               // pwork_vtx_ln_to_gn,
                                               pwork_n_edge,
                                               pwork_edge_vtx,
                                               pwork_edge_ln_to_gn,
                                               &line_to_back_idx,
                                               &line_to_back);


      PDM_dbbtree_free(dbbt);
      PDM_box_set_destroy(&box_set);

      break;
    }

    default: {

      PDM_error(__FILE__, __LINE__, 0, "Method %d not implemented\n", method);

    }
  }


  if (verbose) {
    PDM_log_trace_connectivity_long(line_to_back_idx, line_to_back, pwork_n_edge, "line_to_back ::");
  }


  PDM_g_num_t *pwork_back_vtx_ln_to_gn = NULL;
  int         *pwork_back_face_vtx_idx = NULL;
  int         *pwork_back_face_vtx     = NULL;
  int          pn_work_back_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           dback_face_vtx_idx,
                                                           dback_face_vtx,
                                                           line_to_back_idx[pwork_n_edge],
                                  (const PDM_g_num_t *)    line_to_back,
                                                           &pn_work_back_vtx,
                                                           &pwork_back_vtx_ln_to_gn,
                                                           &pwork_back_face_vtx_idx,
                                                           &pwork_back_face_vtx);
  double **tmp_pwork_back_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        dback_vtx_coord,
                                        &pn_work_back_vtx,
                 (const PDM_g_num_t **) &pwork_back_vtx_ln_to_gn,
                                        &tmp_pwork_back_vtx_coord);
  double *pwork_back_vtx_coord = tmp_pwork_back_vtx_coord[0];
  free(tmp_pwork_back_vtx_coord);


  if (vtk) {
    char filename[999];

    int *halo_color = (int * ) malloc( line_to_back_idx[pwork_n_edge] * sizeof(int));

    for(int i_line = 0; i_line < pwork_n_edge; ++i_line ){
      for(int idx = line_to_back_idx[i_line]; idx < line_to_back_idx[i_line+1]; ++idx) {
        halo_color[idx] = i_line;
      }
    }

    sprintf(filename, "work_back_faces_%d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_work_back_vtx,
                           pwork_back_vtx_coord,
                           pwork_back_vtx_ln_to_gn,
                           line_to_back_idx[pwork_n_edge],
                           pwork_back_face_vtx_idx,
                           pwork_back_face_vtx,
                           line_to_back,
                           halo_color);

    free(halo_color);
  }



  free(pwork_back_vtx_ln_to_gn);
  free(pwork_back_face_vtx_idx);
  free(pwork_back_face_vtx);
  free(pwork_back_vtx_coord);

  free(line_to_back_idx);
  free(line_to_back);


  /*
   *  Free memory
   */
  PDM_multipart_free(mpart);
  PDM_DMesh_nodal_free(dmn);
  free(dback_vtx_coord);
  free(dback_face_vtx_idx);
  free(dback_face_vtx);
  free(back_distrib_vtx);
  free(back_distrib_face);
  free(dback_face_ln_to_gn);
  free(back_face_extents);

  free(pback_vtx_ln_to_gn);
  free(pback_vtx_coord);
  free(pback_face_vtx_idx);
  free(pback_face_vtx);


  PDM_MPI_Finalize ();

  return 0;
}
