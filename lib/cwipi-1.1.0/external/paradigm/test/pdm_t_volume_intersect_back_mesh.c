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
#include "pdm_triangle.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

// static const int    verbose = 0;
// static const int    vtk     = 0;

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
 int           *verbose,
 int           *vtk
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }

    else if (strcmp(argv[i], "-vtk") == 0) {
      *vtk = 1;
    }

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

/* Create volume planes */

static void
_create_volume_4planes
(
double  *edge,
double  *direction,
double   theta,
double   eps,
double  *n,
double  *pt_plane
)
{
  // Gram-Schmidt orthonormalization
  double u[3], v[3], w[3];

  for (int i = 0; i < 3; i++) {
    w[i] = edge[6+i] - edge[i];
  }
  double mw = PDM_MODULE(w);
  for (int i = 0; i < 3; i++) {
    w[i] /= mw;
  }

  double a = -1 * PDM_DOT_PRODUCT(w, direction) / (PDM_MODULE(w)*PDM_MODULE(w));
  for (int i = 0; i < 3; i++) {
    u[i] = direction[i] + a * w[i];
  }
  double mu = PDM_MODULE(u);
  for (int i = 0; i < 3; i++) {
    u[i] /= mu;
  }

  PDM_CROSS_PRODUCT(v, w, u);

  if (PDM_ABS(PDM_DOT_PRODUCT(u, w)) > 1e-13) {
    PDM_error(__FILE__, __LINE__, 0,"!!! u.w = %e\n", PDM_DOT_PRODUCT(u, w));
  }

  double c = cos(theta);
  double s = sin(theta);


  for (int i = 0; i < 3; i++) {

    pt_plane[i] = edge[i] - mw*eps*w[i];
    n[i] = w[i];

    pt_plane[3+i] = edge[6+i] + mw*eps*w[i];
    n[3+i] = -w[i];

    pt_plane[6+i] = edge[3+i] + c*u[i] + s*v[i];
    n[6+i] = s*u[i] - c*v[i];

    pt_plane[9+i] = edge[3+i] + c*u[i] - s*v[i];
    n[9+i] = s*u[i] + c*v[i];

  }

}

/* Ouput volumes in VTK format */

static void
_vtk_write_volume
(
 const char *filename,
 double     *plane_point,
 double     *plane_normal
 )
{
  double u[3], v[3], w[3];
  for (int i = 0; i < 3; i++) {
    u[i] = plane_point[3+i] - plane_point[i];
  }

  PDM_CROSS_PRODUCT(v, u, plane_normal + 6);
  PDM_CROSS_PRODUCT(w, plane_normal + 9, u);

  double vtx_coord[18];
  memcpy(vtx_coord, plane_point, sizeof(double) * 6);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      vtx_coord[6  + 3*i + j] = plane_point[j] + i*u[j] + v[j];
      vtx_coord[12 + 3*i + j] = plane_point[j] + i*u[j] + w[j];
    }
  }

  int face_vtx[14] = {
    1, 5, 3,
    2, 4, 6,
    1, 3, 4, 2,
    1, 2, 6, 5
  };

  int face_vtx_idx[5] = {0, 3, 6, 10, 14};

  PDM_vtk_write_polydata(filename,
                         6,
                         vtx_coord,
                         NULL,
                         4,
                         face_vtx_idx,
                         face_vtx,
                         NULL,
                         NULL);
}

static void
_projection_on_background_mesh_get
(
 double              *pt_to_project_coord,
 int                  n_back_elt,
 int                 *cavity_back_elt,
 int                 *back_elt_vtx_idx,
 int                 *back_elt_vtx,
 double              *back_vtx_coord,
 int                  back_elt_order,
 PDM_Mesh_nodal_elt_t back_elt_type,
 double             **proj_pt_coord
)
{
  int strid = PDM_Mesh_nodal_n_vtx_elt_get(back_elt_type,
                                           back_elt_order);
  int dim   = PDM_Mesh_nodal_elt_dim_get(back_elt_type);

  double best_min_dist2 = 1.0e15;
  *proj_pt_coord = malloc(sizeof(double) * 3);
  double tmp_pt_to_project_coord[3];

  if (dim == 2) {

    double  closestPoint[3];
    double elt_pt_coord[strid * 3];
    double  min_dist2;
    double  *weights  = NULL;
    for (int ielt = 0; ielt < n_back_elt; ielt++) {
      int elt_id = PDM_ABS(cavity_back_elt[ielt]); // not minus 1 here because already -1 in this tab in this test
      int local_head = 0;
      for (int j = back_elt_vtx_idx[elt_id]; j < back_elt_vtx_idx[elt_id+1]; j++) {
        int vtx_id = back_elt_vtx[j] - 1;
        for (int k = 0; k < 3; k++) {
          elt_pt_coord[3*local_head + k] = back_vtx_coord[3*vtx_id + k];
        }
        local_head++;
      }
      memcpy(tmp_pt_to_project_coord, pt_to_project_coord, sizeof(double) * 3);
      PDM_triangle_evaluate_position(tmp_pt_to_project_coord,
                                      elt_pt_coord,
                                      closestPoint,
                                      &min_dist2,
                                      weights);
      if ( min_dist2 < best_min_dist2) {
        best_min_dist2 = min_dist2;
        memcpy(*proj_pt_coord, closestPoint, sizeof(double) * 3);
      }
    } // end loop on back elements

  } // end 2D element case

  // TO DO: 1D element case for ridges

}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int n_part = 1;

  int verbose = 0;
  int vtk     = 0;

  _read_args (argc,
              argv,
              &verbose,
              &vtk);

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_UNUSED(verbose);
  PDM_UNUSED(vtk);

  /* Create surface meshes */

  const double x_center = 0;
  const double y_center = 0;
  const double z_center = 0;

  // Background spherical mesh

  const int    back_nu     = 100;
  const int    back_nv     = 100;
  const double back_radius = 10;

  double      *d_back_vtx_coord    = NULL;
  int         *d_back_face_vtx_idx = NULL;
  PDM_g_num_t *d_back_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx    = NULL;
  PDM_g_num_t *back_distrib_face   = NULL;

  PDM_sphere_surf_gen(comm,
                      back_nu,
                      back_nv,
                      x_center,
                      y_center,
                      z_center,
                      back_radius,
                      &d_back_vtx_coord,
                      &d_back_face_vtx_idx,
                      &d_back_face_vtx,
                      &back_distrib_vtx,
                      &back_distrib_face);

  // "Volume" spherical mesh

  const int    vol_nu     = 10;
  const int    vol_nv     = 10;
  const double vol_radius = 10;

  PDM_dmesh_nodal_t *vol_dmn;

  PDM_sphere_surf_gen_nodal(comm,
                            vol_nu,
                            vol_nv,
                            x_center,
                            y_center,
                            z_center,
                            vol_radius,
                            &vol_dmn);

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

  PDM_multipart_dmesh_nodal_set(mpart, 0, vol_dmn);

  PDM_multipart_compute(mpart);

  free(n_part_domains);

  // Vertices
  int     i_part          = 0;
  double *p_vol_vtx_coord = NULL;
  int p_vol_n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &p_vol_vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *vol_vtx_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &vol_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Edges
  int *p_vol_edge_vtx     = NULL;
  int *p_vol_edge_vtx_idx = NULL;
  int  p_vol_n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                          0,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                          &p_vol_edge_vtx_idx,
                                                          &p_vol_edge_vtx,
                                                          PDM_OWNERSHIP_KEEP);

  if (p_vol_edge_vtx_idx != NULL) free(p_vol_edge_vtx_idx);

  PDM_g_num_t *vol_edge_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &vol_edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Faces
  int *p_vol_face_edge     = NULL;
  int *p_vol_face_edge_idx = NULL;
  int p_vol_n_face = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                         &p_vol_face_edge_idx,
                                                         &p_vol_face_edge,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *vol_face_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &vol_face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  int *p_vol_face_vtx = NULL;
  PDM_compute_face_vtx_from_face_and_edge(p_vol_n_face,
                                          p_vol_face_edge_idx,
                                          p_vol_face_edge,
                                          p_vol_edge_vtx,
                                          &p_vol_face_vtx);

  // Get coordinates of the background mesh associated to the local face distribution

  int dn_back_face = back_distrib_face[i_rank+1] - back_distrib_face[i_rank];

  PDM_g_num_t *d_back_face_ln_to_gn = (PDM_g_num_t * ) malloc(dn_back_face * sizeof(PDM_g_num_t));
  for (int i = 0; i < dn_back_face; ++i) {
    d_back_face_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }

  int dn_back_vtx = back_distrib_vtx[i_rank+1] - back_distrib_vtx[i_rank];

  PDM_g_num_t *d_back_vtx_ln_to_gn = (PDM_g_num_t * ) malloc(dn_back_vtx * sizeof(PDM_g_num_t));
  for (int i = 0; i < dn_back_vtx; ++i) {
    d_back_vtx_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }
  free(d_back_vtx_ln_to_gn);

  PDM_g_num_t *p_back_vtx_ln_to_gn = NULL;
  int         *p_back_face_vtx_idx = NULL;
  int         *p_back_face_vtx     = NULL;
  int          p_back_n_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           d_back_face_vtx_idx,
                                                           d_back_face_vtx,
                                                           dn_back_face,
                                  (const PDM_g_num_t *)    d_back_face_ln_to_gn,
                                                           &p_back_n_vtx,
                                                           &p_back_vtx_ln_to_gn,
                                                           &p_back_face_vtx_idx,
                                                           &p_back_face_vtx);

  int      *n_part_p_back_n_vtx            = malloc(sizeof(int) * n_part);
  n_part_p_back_n_vtx[0]                   = p_back_n_vtx;
  double  **n_part_p_back_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        d_back_vtx_coord,
                                        n_part_p_back_n_vtx,
                 (const PDM_g_num_t **) &p_back_vtx_ln_to_gn,
                                        &n_part_p_back_vtx_coord);
  double *p_back_vtx_coord = n_part_p_back_vtx_coord[0];
  free(n_part_p_back_vtx_coord);

  // Create the extents faces as a partition and get associated coords
  double      *background_box_extents = malloc(sizeof(double)      * dn_back_face * 6);
  double       eps                    = 1.0e-6;
  for (int iface = 0; iface < dn_back_face; iface++) {
    double *tmp_extents = background_box_extents + 6*iface;
    for (int k = 0; k < 3; k++) {
      tmp_extents[k]     =  1.0e15;
      tmp_extents[3 + k] = -1.0e15;
    }
    for (int ivtx = p_back_face_vtx_idx[iface]; ivtx < p_back_face_vtx_idx[iface+1]; ivtx++) {
      int vtx_id = p_back_face_vtx[ivtx]-1;
      double *tmp_coord = p_back_vtx_coord + 3*vtx_id;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }
    } // end loop on vertices of iface
    for (int k = 0; k < 3; k++) {
      if (PDM_ABS(tmp_extents[k] - tmp_extents[3 + k]) < 1.0e-15) { //
        tmp_extents[k]     -= eps;
        tmp_extents[3 + k] += eps;
      }
    }
  } // end loop on background faces

  if (vtk) {
    char filename3[999];
    sprintf(filename3, "dbbtree_boxes_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename3,
                        dn_back_face,
                        background_box_extents,
                        d_back_face_ln_to_gn);
  }

  // Create dbbtree from surface mesh boxes

  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < dn_back_face; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN(l_extents[k],     background_box_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX(l_extents[k + 3], background_box_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                 n_part,
                                                 &dn_back_face,
                               (const double **) &background_box_extents,
                          (const PDM_g_num_t **) &d_back_face_ln_to_gn);

  free(background_box_extents);
  // free(d_back_face_ln_to_gn);

  // Create volumes using normals (first "volume" mesh then background mesh)

  int    *edge_face_normal_stride = PDM_array_zeros_int(p_vol_n_edge);
  double *edge_face_normal        = malloc(sizeof(double) * 3 * p_vol_n_edge * 2);

  double *face_center             = malloc(sizeof(double) * 3 * p_vol_n_face);
  double *face_normal             = malloc(sizeof(double) * 3 * p_vol_n_face);

  for (int iface = 0; iface < p_vol_n_face; iface++) {
    // compute normal vector
    int *tmp_face_edge = p_vol_face_edge + p_vol_face_edge_idx[iface];
    // with face_edge edge_vtx
    int edge_id1 = PDM_ABS(tmp_face_edge[0]) - 1;
    int vtx_id1, vtx_id2;
    if (tmp_face_edge[0] > 0) {
      vtx_id1 = p_vol_edge_vtx[2*edge_id1]-1;
      vtx_id2 = p_vol_edge_vtx[2*edge_id1+1]-1;
    } else {
      vtx_id1 = p_vol_edge_vtx[2*edge_id1+1]-1;
      vtx_id2 = p_vol_edge_vtx[2*edge_id1]-1;
    }
    int edge_id2 = PDM_ABS(tmp_face_edge[1]) - 1;
    int vtx_id3;
    if (p_vol_edge_vtx[2*edge_id2]-1 != vtx_id1 && p_vol_edge_vtx[2*edge_id2]-1 != vtx_id2) {
      vtx_id3 = p_vol_edge_vtx[2*edge_id2]-1;
    } else {
      vtx_id3 = p_vol_edge_vtx[2*edge_id2+1]-1;
    }

    double vect[6];
    for (int i = 0; i < 3; i++) {
      vect[i]     = p_vol_vtx_coord[3*vtx_id2 + i] - p_vol_vtx_coord[3*vtx_id1 + i];
      vect[3 + i] = p_vol_vtx_coord[3*vtx_id3 + i] - p_vol_vtx_coord[3*vtx_id1 + i];
    }
    double normal[3];
    PDM_CROSS_PRODUCT(normal, vect, vect + 3);

    for (int i = 0; i < 3; i++) {
      face_center[3*iface + i] = 1./3 * (p_vol_vtx_coord[3*vtx_id1 + i] +  p_vol_vtx_coord[3*vtx_id2 + i] + p_vol_vtx_coord[3*vtx_id3 + i]);
      face_normal[3*iface + i] = normal[i];
    }


    // fill in edge_face_normal
    for (int iedge = 0; iedge < 3; iedge++) {
      int edge_id = PDM_ABS(tmp_face_edge[iedge])-1;
      for (int i = 0; i < 3; i++) {
        edge_face_normal[6*edge_id + 3*edge_face_normal_stride[edge_id] + i] = normal[i];
      }
      edge_face_normal_stride[edge_id]++;
    } // end loop on all edges of iface
  } // end loop on volume faces

  // Remove empty cell in edge_face_normal
  int *edge_face_normal_idx = PDM_array_new_idx_from_sizes_int(edge_face_normal_stride, p_vol_n_edge);
  int idx_read  = 0;
  int idx_write = 0;
  for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
    idx_write = edge_face_normal_idx[iedge];
    idx_read  = 2*iedge;
    // if (idx_read != idx_write) {
      for (int inormal = 0; inormal < edge_face_normal_stride[iedge]; inormal++) {
        for (int i = 0; i < 3; i++) {
          edge_face_normal[3*idx_write + 3*inormal + i] = edge_face_normal[3*idx_read + 3*inormal + i];
        }
      }
    // }
  } // end loop on edges

  if (verbose) {
    PDM_log_trace_array_double(edge_face_normal, edge_face_normal_idx[p_vol_n_edge], "edge_face_normal");
  }
  free(edge_face_normal_idx);

  // Get normal contributions for other procs

  // Part to Block

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &vol_edge_ln_to_gn,
                                                      NULL,
                                                      &p_vol_n_edge,
                                                      n_part,
                                                      comm);

  int    *block_stride = NULL;
  double *block_data   = NULL;
  PDM_part_to_block_exch(ptb,
                         3*sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &edge_face_normal_stride,
               (void **) &edge_face_normal,
                         &block_stride,
               (void **) &block_data);

  int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get(ptb);

  // Block to Part

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_g_num,
                                                                        n_elt_block,
                                                 (const PDM_g_num_t **) &vol_edge_ln_to_gn,
                                                                        &p_vol_n_edge,
                                                                        n_part,
                                                                        comm);
  int two = 2;
  PDM_block_to_part_exch_in_place(btp,
                                  3*sizeof(double),
                                  PDM_STRIDE_CST_INTERLACED,
                                  &two,
                         (void *) block_data,
                                  NULL,
                        (void **) &edge_face_normal);

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);

  free(block_data);
  free(block_stride);

  // compute volume angle

  // only for vtk
  // double *middle_pt_coord = malloc(sizeof(double) * 3 * p_vol_n_edge);
  double *direction_vect  = malloc(sizeof(double) * 3 * p_vol_n_edge);

  double  theta_min         = 1.0e-1; // WARNING: randomly chosen value
  double  eps2              = 1.0e-1; // WARNING: randomly chosen value
  double *edge_normal       = malloc(sizeof(double) * p_vol_n_edge * 3 * 4);
  double *edge_pt_plane     = malloc(sizeof(double) * p_vol_n_edge * 3 * 4);
  for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
    double theta;
    double direction[3];
    for (int i = 0; i < 3; i++) {
      direction[i] = 0.5 * (edge_face_normal[6*iedge + i] + edge_face_normal[6*iedge + 3 + i]);
    }
    double dot_prod = PDM_DOT_PRODUCT(edge_face_normal + 6*iedge, edge_face_normal + 6*iedge + 3);
    double module1 = PDM_MODULE(edge_face_normal + 6*iedge);
    double module2 = PDM_MODULE(edge_face_normal + 6*iedge + 3);
    theta = acos(dot_prod / (module1 * module2)); // make sure between -1 and 1
    theta += theta_min;
    if (theta > 3) { // WARNING: shouldn't be bigger than PI
      theta = 3;
    }
    if (theta < 0) {
      log_trace("theta is negative!!!\n");
      theta = PDM_ABS(theta);
    }

    int *tmp_edge_vtx = p_vol_edge_vtx + 2*iedge;
    int vtx_id1 = tmp_edge_vtx[0] - 1;
    int vtx_id2 = tmp_edge_vtx[1] - 1;
    // double edge_vector[3];
    double edge[9];
    for (int i = 0; i < 3; i++) {
      edge[i]     = p_vol_vtx_coord[3*vtx_id1 + i];
      edge[6 + i] = p_vol_vtx_coord[3*vtx_id2 + i];
      edge[3 + i]     = 0.5 * (p_vol_vtx_coord[3*vtx_id1 + i] + p_vol_vtx_coord[3*vtx_id2 + i]);
      // edge_vector[i] = edge[6 + i] - edge[i];
    }

    // enfoncer le volume
    for (int i = 0; i < 3; i++) {
      edge[i]     = edge[i] -0.3*direction[i];
      edge[6 + i] = edge[6 + i] -0.3*direction[i];
      edge[3 + i] = 0.5 * (edge[i] + edge[6 + i]);
    }

    double *normal   = edge_normal + 12*iedge;
    double *pt_plane = edge_pt_plane + 12*iedge;

    _create_volume_4planes(edge,
                           direction,
                           theta,
                           eps2,
                           normal,
                           pt_plane);

    if (vtk) {

      char filename4[999];
      sprintf(filename4, "volume_of_edge_id_"PDM_FMT_G_NUM".vtk", vol_edge_ln_to_gn[iedge]);
      _vtk_write_volume(filename4,
                        pt_plane,
                        normal);

    const char *normal_name = "normal";

    char filename38[999];
    sprintf(filename38, "normal_of_edge_id_"PDM_FMT_G_NUM".vtk", vol_edge_ln_to_gn[iedge]);
    PDM_vtk_write_point_cloud_with_field(filename38,
                                         4,
                                         pt_plane,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &normal,
                                         0,
                                         NULL,
                                         NULL);

    }


  } // end loop on edges

  if (vtk) {

    // const char *direction_name = "direction";

    // char filename5[999];
    // sprintf(filename5, "direction_%d.vtk", i_rank);
    // PDM_vtk_write_point_cloud_with_field(filename5,
    //                                      p_vol_n_edge,
    //                                      middle_pt_coord,
    //                                      vol_edge_ln_to_gn,
    //                                      NULL,
    //                                      0,
    //                                      NULL,
    //                                      NULL,
    //                                      1,
    //                      (const char **) &direction_name,
    //                    (const double **) &direction_vect,
    //                                      0,
    //                                      NULL,
    //                                      NULL);

    const char *normal_name = "face_normal";

    char filename6[999];
    sprintf(filename6, "face_normal_%d.vtk", i_rank);
    PDM_vtk_write_point_cloud_with_field(filename6,
                                         p_vol_n_face,
                                         face_center,
                                         vol_face_ln_to_gn,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &face_normal,
                                         0,
                                         NULL,
                                         NULL);
  }

  // Call PDM_dbbtree_volumes_intersect_boxes

  int *volume_plane_idx = PDM_array_new_idx_from_const_stride_int(4, p_vol_n_edge);

  int         *volume_boxes_idx   = NULL;
  PDM_g_num_t *volume_boxes_g_num = NULL;
  PDM_dbbtree_volumes_intersect_boxes(dbbt,
                                      p_vol_n_edge,
                                      vol_edge_ln_to_gn,
                                      volume_plane_idx,
                                      edge_normal,
                                      edge_pt_plane,
                                      &volume_boxes_idx,
                                      &volume_boxes_g_num);

  if (vtk) {
    char filename1[999];
    sprintf(filename1, "dbbt_%d.vtk", i_rank);
    PDM_dbbtree_box_tree_write_vtk(filename1,
                                   dbbt,
                                   -1,
                                   0);
  }

  if (verbose) {
    log_trace("VOLUME-BOX INTERSECTION\n");
    for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
      log_trace("--> volume "PDM_FMT_G_NUM" is intersected by ", vol_edge_ln_to_gn[iedge]);
      log_trace(" %d boxes ", volume_boxes_idx[iedge+1] - volume_boxes_idx[iedge]);
      for (int i = volume_boxes_idx[iedge]; i < volume_boxes_idx[iedge+1]; i++) {
        log_trace("%d ", volume_boxes_g_num[i]);
      }
      log_trace("\n");
    }
  }
  free(edge_normal);
  free(edge_pt_plane);
  free(direction_vect);
  free(volume_plane_idx);
  free(d_back_face_ln_to_gn);

  // VTK output of surface mesh with tagged elements for each volume mesh edges

  int n_boxes = volume_boxes_idx[p_vol_n_edge];

  PDM_g_num_t *p_box_volume_g_num  = malloc(sizeof(PDM_g_num_t) * n_boxes);

  for (int ivol = 0; ivol < p_vol_n_edge; ivol++) {
    for (int ibox = volume_boxes_idx[ivol]; ibox < volume_boxes_idx[ivol+1]; ibox++) {
      p_box_volume_g_num[ibox] = vol_edge_ln_to_gn[ivol];
    } // end loop on ivol boxes
  } // end loop on volumes

  // Part to Block

  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                                                    1.,
                                                                    &volume_boxes_g_num,
                                                                    back_distrib_face,
                                                                    &n_boxes,
                                                                    n_part,
                                                                    comm);

  int         *d_box_volume_stride = NULL;
  PDM_g_num_t *d_box_volume_g_num  = NULL;
  int         *part_stride         = PDM_array_const_int(n_boxes, 1);
  PDM_part_to_block_exch(ptb2,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &p_box_volume_g_num,
                         &d_box_volume_stride,
               (void **) &d_box_volume_g_num);
  free(p_box_volume_g_num);

  PDM_g_num_t* distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb2,
                                                                        &d_box_volume_stride,
                                                                        back_distrib_face[n_rank]);
  free(distrib);

  // int dn_block_face =  PDM_part_to_block_n_elt_block_get(ptb2);

  free(part_stride);
  PDM_part_to_block_free(ptb2);

  int *d_box_volume_idx = PDM_array_new_idx_from_sizes_int(d_box_volume_stride, dn_back_face);
  free(d_box_volume_stride);

  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, vol_dmn);
  PDM_dmesh_nodal_generate_distribution(vol_dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
  PDM_dmesh_t* vol_dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &vol_dmesh);

  PDM_g_num_t* dedge_distrib = NULL;
  PDM_dmesh_distrib_get(vol_dmesh, PDM_MESH_ENTITY_EDGE, &dedge_distrib);

  int total_n_edges = dedge_distrib[n_rank];

  PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  PDM_DMesh_nodal_free(vol_dmn);

  int  **volume       = malloc(sizeof(int  *) * total_n_edges);
  char **volume_names = malloc(sizeof(char *) * total_n_edges);

  for (int ivol = 0; ivol < total_n_edges; ivol++) {
    volume[ivol] = PDM_array_zeros_int(dn_back_face);
    volume_names[ivol] = malloc(sizeof(char) * 99);
    sprintf(volume_names[ivol], "edge_%d.vtk", ivol+1);

  }

  if (vtk) {

  for (int ibox = 0; ibox < dn_back_face; ibox++) {
    // PDM_g_num_t box_gnum = d_back_face_ln_to_gn[ibox];
    for (int ivol = d_box_volume_idx[ibox]; ivol < d_box_volume_idx[ibox+1]; ivol++) {
      PDM_g_num_t vol_gnum = d_box_volume_g_num[ivol];
      volume[vol_gnum-1][ibox] = 1;
    }
  }


    char filename1[999];
    sprintf(filename1, "background_mesh_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename1,
                               p_back_n_vtx,
                               p_back_vtx_coord,
                               p_back_vtx_ln_to_gn,
                               PDM_MESH_NODAL_TRIA3,
                               dn_back_face,
                               p_back_face_vtx,
                               d_back_face_ln_to_gn,
                               total_n_edges,
               (const char **) volume_names,
                (const int **) volume);

    char filename2[999];
    sprintf(filename2, "volume_mesh_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename2,
                               p_vol_n_vtx,
                               p_vol_vtx_coord,
                               vol_vtx_ln_to_gn,
                               PDM_MESH_NODAL_TRIA3,
                               p_vol_n_face,
                               p_vol_face_vtx,
                               vol_face_ln_to_gn,
                               0,
                               NULL,
                               NULL);

  }

  PDM_dbbtree_free(dbbt);
  PDM_box_set_destroy(&box_set);
  free(d_box_volume_idx);
  free(d_box_volume_g_num);
  for (int ivol = 0; ivol < total_n_edges; ivol++) {
    free(volume[ivol]);
    free(volume_names[ivol]);
  }
  free(volume);
  free(volume_names);
  free(p_vol_face_vtx);



  // Create a fake cavity distribution to mimic the pdm_mesh_adaptation setting

  int d_n_cavity = p_vol_n_edge;
  PDM_g_num_t *distrib_cavity  = malloc (sizeof(PDM_g_num_t) * (d_n_cavity + 1));
  PDM_distrib_compute(d_n_cavity,
                      distrib_cavity,
                      -1,
                      comm);
  free(distrib_cavity);
  // PDM_g_num_t *cavity_ln_to_gn = vol_edge_ln_to_gn;

  // Retreive associated edges (just the edges already in partition of rank i)

  // PDM_g_num_t *pcavity_seed_edge_g_num   = vol_edge_ln_to_gn;
  int         *seed_edge_back_face_idx   = volume_boxes_idx;
  PDM_g_num_t *seed_edge_back_face_g_num = volume_boxes_g_num;
  int         *seed_edge_vtx             = p_vol_edge_vtx;
  double      *seed_edge_vtx_coord       = p_vol_vtx_coord;

  // Create a back_face partition

  int p_n_back_face = 0;

  // Unique sort

  PDM_g_num_t *toto_g_num = malloc(sizeof(PDM_g_num_t) * seed_edge_back_face_idx[p_vol_n_edge]);
  memcpy(toto_g_num, seed_edge_back_face_g_num, sizeof(PDM_g_num_t) * seed_edge_back_face_idx[p_vol_n_edge]);

  int* order = (int *) malloc(seed_edge_back_face_idx[p_vol_n_edge] * sizeof(int));
  for(int i = 0; i < seed_edge_back_face_idx[p_vol_n_edge]; ++i){
    order[i] = i;
  }
  PDM_sort_long(seed_edge_back_face_g_num, order, seed_edge_back_face_idx[p_vol_n_edge]);

  if(verbose) {
    PDM_log_trace_array_int(order, seed_edge_back_face_idx[p_vol_n_edge],"order: ");
  }

  int *seed_edge_back_face_l_num = malloc(sizeof(int) * seed_edge_back_face_idx[p_vol_n_edge]);
  seed_edge_back_face_l_num[order[0]] = 0;
  PDM_g_num_t *p_back_face_ln_to_gn = malloc(sizeof(PDM_g_num_t) * seed_edge_back_face_idx[p_vol_n_edge]);
  p_back_face_ln_to_gn[0] = seed_edge_back_face_g_num[0];

  int read_idx = 0;
  int write_idx = 1;
  PDM_g_num_t last_val = seed_edge_back_face_g_num[0];
  int last_l_num = 0;
  for (int i = 1; i < seed_edge_back_face_idx[p_vol_n_edge]; i++) {
    read_idx = i;
    if (seed_edge_back_face_g_num[read_idx] != last_val) {
      last_val = seed_edge_back_face_g_num[read_idx];
      p_back_face_ln_to_gn[write_idx++] = last_val;
      last_l_num++;
    }
    seed_edge_back_face_l_num[order[read_idx]] = last_l_num;
  }

  p_n_back_face = write_idx;

  int* unique_order = (int *) malloc(seed_edge_back_face_idx[p_vol_n_edge] * sizeof(int));
  int new_size = PDM_inplace_unique_long2(toto_g_num, unique_order, 0, seed_edge_back_face_idx[p_vol_n_edge]-1);

  if(verbose) {
    PDM_log_trace_array_long(toto_g_num, new_size, "unique back_face_ln_to_gn: ");
    PDM_log_trace_array_long(p_back_face_ln_to_gn, p_n_back_face, "hand made back_face_ln_to_gn: ");

    PDM_log_trace_array_int(unique_order, seed_edge_back_face_idx[p_vol_n_edge], "unique_order: ");
    PDM_log_trace_array_int(seed_edge_back_face_l_num, seed_edge_back_face_idx[p_vol_n_edge], "seed_edge_back_face_l_num: ");
  }

  // PDM_part_dconnectivity_to_pconnectivity_sort_single_part to get back_face->vtx and coord

  free(p_back_vtx_ln_to_gn);
  free(p_back_face_vtx_idx);
  free(p_back_face_vtx    );
  free(unique_order    );
  free(order    );
  free(toto_g_num    );

  // PDM_g_num_t *p_back_vtx_ln_to_gn = NULL;
  // int         *p_back_face_vtx_idx = NULL;
  // int         *p_back_face_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           d_back_face_vtx_idx,
                                                           d_back_face_vtx,
                                                           p_n_back_face,
                                  (const PDM_g_num_t *)    p_back_face_ln_to_gn,
                                                           &p_back_n_vtx,
                                                           &p_back_vtx_ln_to_gn,
                                                           &p_back_face_vtx_idx,
                                                           &p_back_face_vtx);

  free(n_part_p_back_n_vtx    );
  free(p_back_face_ln_to_gn);
  free(p_back_vtx_coord);

  n_part_p_back_n_vtx    = malloc(sizeof(int) * n_part);
  n_part_p_back_n_vtx[0] = p_back_n_vtx;
  // double  **n_part_p_back_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        d_back_vtx_coord,
                                        n_part_p_back_n_vtx,
                 (const PDM_g_num_t **) &p_back_vtx_ln_to_gn,
                                        &n_part_p_back_vtx_coord);

  p_back_vtx_coord = n_part_p_back_vtx_coord[0];
  free(n_part_p_back_vtx_coord);
  free(n_part_p_back_n_vtx);

  // Do projections

  // for VTK
  // double *back_face_proj_pts = malloc(sizeof(double) * p_n_back_face * 3);
  double *back_face_proj_pts = malloc(sizeof(double) * p_n_back_face * 3);
  // int     n_back_face_proj_pts;

  for (int icav = 0; icav < d_n_cavity; icav++) {
    int edge_id = icav; // tmp in mesh_adaptation need an idx
    // int edge_id = PDM_ABS(pcavity_seed_edge_g_num[icav])-1;
    int vtx_id1 = seed_edge_vtx[2*edge_id]   - 1;
    int vtx_id2 = seed_edge_vtx[2*edge_id+1] - 1;
    double middle_pt[3];
    for (int j = 0; j < 3; j++) {
      middle_pt[j] = 0.5 * (seed_edge_vtx_coord[3*vtx_id1+j] + seed_edge_vtx_coord[3*vtx_id2+j]);
    }

    int n_back_elt = seed_edge_back_face_idx[edge_id+1] - seed_edge_back_face_idx[edge_id];

    double *proj_pt_coord   = NULL;
    int    *cavity_back_elt = seed_edge_back_face_l_num + seed_edge_back_face_idx[edge_id];
    _projection_on_background_mesh_get(middle_pt,
                                       n_back_elt,
                                       cavity_back_elt,
                                       p_back_face_vtx_idx,
                                       p_back_face_vtx,
                                       p_back_vtx_coord,
                                       2,
                                       PDM_MESH_NODAL_TRIA3,
                                       &proj_pt_coord);
    memcpy(back_face_proj_pts + 3*icav, proj_pt_coord, sizeof(double) * 3);
    free(proj_pt_coord);

  } // end loop on cavities

  if (vtk) {

    char filename42[999];
    sprintf(filename42, "cavity_proj_pts_%d.vtk", i_rank);
    PDM_vtk_write_point_cloud_with_field(filename42,
                                         d_n_cavity,
                                         back_face_proj_pts,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL);
  }

    // int face_head = 0;
    // for (int iface = seed_edge_back_face_idx[edge_id]; iface < seed_edge_back_face_idx[edge_id+1]; iface++) {
    //   int back_face_id = seed_edge_back_face_l_num[iface];

    //   double face_pt_coord[9];
    //   int local_head = 0;
    //   for (int j = p_back_face_vtx_idx[back_face_id]; j < p_back_face_vtx_idx[back_face_id+1]; j++) {
    //     int vtx_id = p_back_face_vtx[j] - 1;
    //     for (int k = 0; k < 3; k++) {
    //       face_pt_coord[3*local_head + k] = p_back_vtx_coord[3*vtx_id + k];
    //     }
    //     local_head++;
    //   }
    //   double  closestPoint[3];
    //   double  min_dist2;
    //   double  *weights      = NULL;
    //   PDM_triangle_evaluate_position(middle_pt,
    //                                  face_pt_coord,
    //                                  closestPoint,
    //                                  &min_dist2,
    //                                  weights);

    //   memcpy(back_face_proj_pts + 3*face_head, closestPoint, sizeof(double) * 3);
    //   face_head++;
  PDM_multipart_free(mpart);
  free(seed_edge_back_face_idx);
  free(edge_face_normal_stride);
  free(edge_face_normal);
  free(back_face_proj_pts);
  free(face_center);
  free(face_normal);
  free(p_back_vtx_ln_to_gn);
  free(p_back_face_vtx_idx);
  free(p_back_face_vtx    );
  free(seed_edge_back_face_l_num    );
  free(seed_edge_back_face_g_num    );
  free(p_back_vtx_coord    );

  free(d_back_vtx_coord   );
  free(d_back_face_vtx_idx);
  free(d_back_face_vtx    );
  free(back_distrib_vtx   );
  free(back_distrib_face  );


  PDM_MPI_Finalize ();
}
