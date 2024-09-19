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
#include "pdm_distrib.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"

#include "pdm_vtk.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
           PDM_g_num_t   *gn_src,
           PDM_g_num_t   *gn_tgt,
           int           *n_part,
           int           *post)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *gn_src = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *gn_tgt = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_gen_lines
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn,
 int                *n_vtx,
 double            **vtx_coord,
 PDM_g_num_t       **vtx_g_num,
 int                *n_line,
 PDM_g_num_t       **line_g_num,
 int               **line_vtx_idx,
 int               **line_vtx
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution (comm,
                                                                  gn);
  int n = (int) (distrib[i_rank+1] - distrib[i_rank]);

  *n_vtx = n + 1;
  *vtx_g_num = malloc (sizeof(PDM_g_num_t) * (*n_vtx));
  *vtx_coord = malloc (sizeof(double     ) * (*n_vtx) * 3);
  for (int i = 0; i < (*n_vtx); i++) {
    PDM_g_num_t g = distrib[i_rank] + i;
    (*vtx_g_num)[i] = 1 + g%gn;

    double t = (double) g / (double) (gn - 1);
    double r = 0.7 + 0.3*cos(10*PDM_PI*t);
    (*vtx_coord)[3*i    ] = r*cos(2*PDM_PI*t) + 1.;
    (*vtx_coord)[3*i + 1] = r*sin(2*PDM_PI*t) + 1.;
    (*vtx_coord)[3*i + 2] = 0.2*sin(8*PDM_PI*t) + 1.;
  }


  *n_line = n;
  *line_g_num   = malloc (sizeof(PDM_g_num_t) * (*n_line));
  *line_vtx_idx = malloc (sizeof(int        ) * (*n_line + 1));
  *line_vtx     = malloc (sizeof(int        ) * (*n_line) * 2);
  (*line_vtx_idx)[0] = 0;
  for (int i = 0; i < *n_line; i++) {
    (*line_g_num)[i] = distrib[i_rank] + i + 1;
    (*line_vtx_idx)[i+1] = (*line_vtx_idx)[i];
    (*line_vtx)[(*line_vtx_idx)[i+1]++] = i + 1;
    (*line_vtx)[(*line_vtx_idx)[i+1]++] = i + 2;
  }

  free (distrib);
}



static void
_gen_point_cloud
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_pts,
 const double        length,
 int                *n_pts,
 double            **pts_coord,
 PDM_g_num_t       **pts_g_num
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  *n_pts = (int) (gn_pts / n_rank);
  if (i_rank < gn_pts % n_rank) {
    (*n_pts)++;
  }

  srand(0);


  PDM_g_num_t* distrib_pts = PDM_compute_entity_distribution(comm, (*n_pts));
  for(int i = 0; i < 3 * distrib_pts[i_rank]; ++i) {
    rand();
  }

  *pts_coord = malloc (sizeof(double) * (*n_pts) * 3);
  for (int i = 0; i < *n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      (*pts_coord)[3*i + j] = length * (double) rand() / ((double) RAND_MAX);
    }
  }


  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);
  double *char_length = malloc (sizeof(double) * (*n_pts));
  for (int i = 0; i < *n_pts; i++) {
    char_length[i] = length * 1e-6;
  }

  PDM_gnum_set_from_coords (gen_gnum, 0, *n_pts, *pts_coord, char_length);

  PDM_gnum_compute (gen_gnum);

  *pts_g_num = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (char_length);

  free(distrib_pts);
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

  PDM_g_num_t gn_src = 10;
  PDM_g_num_t gn_tgt = 10;
  int         n_part = 1;
  int         post   = 0;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &gn_src,
             &gn_tgt,
             &n_part,
             &post);
  n_part = 1;

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   *  Create lines
   */
  int          n_vtx;
  PDM_g_num_t *vtx_g_num;
  double      *vtx_coord;
  int          n_line;
  PDM_g_num_t *line_g_num;
  int         *line_vtx_idx;
  int         *line_vtx;
  _gen_lines (comm,
              gn_src,
              &n_vtx,
              &vtx_coord,
              &vtx_g_num,
              &n_line,
              &line_g_num,
              &line_vtx_idx,
              &line_vtx);




  /*
   *  Create target point cloud
   */
  int          n_pts;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  _gen_point_cloud (comm,
                    gn_tgt,
                    2.,
                    &n_pts,
                    &pts_coord,
                    &pts_g_num);

  /*
   *  Create distance structure
   */
  PDM_g_num_t ln_vtx = n_vtx;
  PDM_g_num_t gn_vtx;
  PDM_MPI_Allreduce (&ln_vtx,  &gn_vtx,  1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  PDM_g_num_t ln_line = n_line;
  PDM_g_num_t gn_line;
  PDM_MPI_Allreduce (&ln_line, &gn_line, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 n_part);

  PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                          0,//i_part,
                                          n_line,
                                          line_vtx_idx,
                                          line_vtx,
                                          line_g_num,
                                          n_vtx,
                                          vtx_coord,
                                          vtx_g_num);

  PDM_dist_cloud_surf_n_part_cloud_set (dist, 0, n_part);

  PDM_dist_cloud_surf_cloud_set (dist,
                                 0,
                                 0,//i_part,
                                 n_pts,
                                 pts_coord,
                                 pts_g_num);

  if (i_rank == 0) {
    printf("-- Dist compute\n");
    fflush(stdout);
  }

  // PDM_dist_cloud_surf_compute (dist);
  PDM_dist_cloud_surf_compute (dist);


  if (post) {

    char filename[999];

    double *line_coord = malloc (sizeof(double) * n_line * 6);
    int idx = 0;
    for (int i = 0; i < n_line; i++) {
      for (int j = line_vtx_idx[i]; j < line_vtx_idx[i+1]; j++) {
        int ivtx = line_vtx[j] - 1;
        for (int k = 0; k < 3; k++) {
          line_coord[idx++] = vtx_coord[3*ivtx + k];
        }
      }
    }

    sprintf(filename, "dist_ridge_line_%3.3d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         n_line,
                         line_coord,
                         line_g_num,
                         NULL);
    free (line_coord);



    sprintf(filename, "dist_ridge_pts_%3.3d.vtk", i_rank);
    PDM_vtk_write_point_cloud (filename,
                               n_pts,
                               pts_coord,
                               pts_g_num,
                               NULL);

    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             0,//i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    line_coord = malloc (sizeof(double) * n_pts * 6);
    idx = 0;
    for (int i = 0; i < n_pts; i++) {

      for (int k = 0; k < 3; k++) {
        line_coord[idx++] = pts_coord[3*i + k];
      }

      for (int k = 0; k < 3; k++) {
       line_coord[idx++] = projected[3*i + k];
      }
    }

    sprintf(filename, "dist_ridge_proj_%3.3d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         n_pts,
                         line_coord,
                         pts_g_num,
                         NULL);
    free (line_coord);

  }




  PDM_dist_cloud_surf_dump_times(dist);
  PDM_dist_cloud_surf_free (dist);


  free (line_vtx_idx);
  free (line_vtx);
  free (line_g_num);
  free (vtx_coord);
  free (vtx_g_num);
  free (pts_coord);
  free (pts_g_num);


  PDM_MPI_Finalize();

  return 0;
}
