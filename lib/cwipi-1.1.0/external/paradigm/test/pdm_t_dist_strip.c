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
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Random
 *
 */

static double
_frand_a_b(double a, double b){
  return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

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
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
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
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
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
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[n] : 30,60
// @@@args[part_kind] : -parmetis
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t  n_vtx_seg  = 3;
  double        length  = 1.;
  int           n_part   = 1;
  int           post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */

  if (n_part != 1) {
    printf ("Erreur n_part != 1\n");
    abort();
  }

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  srand(i_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell      = NULL;
  int          *dface_vtx_idx   = NULL;
  PDM_g_num_t  *dface_vtx       = NULL;
  double       *dvtx_coord      = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group     = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  const double xmax = xmin + length;
  const double ymax = ymin + length;
  const double zmax = zmin + length;

  if (i_rank == 0) {
    printf("-- Build cube\n");
    fflush(stdout);
  }

  PDM_dcube_t* dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD,
                                          n_vtx_seg,
                                          length,
                                          xmin,
                                          ymin,
                                          zmin,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dface_vtxL,
                      &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                       &dface_cell,
                       &dface_vtx_idx,
                       &dface_vtx,
                       &dvtx_coord,
                       &dface_group_idx,
                       &dface_group);
  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (i_rank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_t *ppart = PDM_part_create(PDM_MPI_COMM_WORLD,
                                      method,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      "PDM_PART_RENUM_FACE_NONE",
                                      n_property_cell,
                                      renum_properties_cell,
                                      n_property_face,
                                      renum_properties_face,
                                      n_part,
                                      dn_cell,
                                      dn_face,
                                      dn_vtx,
                                      n_face_group,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      have_dcell_part,
                                      dcell_part,
                                      dface_cell,
                                      dface_vtx_idx,
                                      dface_vtx,
                                      NULL,
                                      dvtx_coord,
                                      NULL,
                                      dface_group_idx,
                                      dface_group);

  free(dcell_part);

  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            PDM_MPI_COMM_WORLD,
                                                            PDM_OWNERSHIP_KEEP);

  int **select_face = malloc (sizeof(int *) * n_part);
  int *n_select_face = malloc (sizeof(int) * n_part);
  int **select_vtx = malloc (sizeof(int *) * n_part);
  int *n_select_vtx = malloc (sizeof(int) * n_part);

  int **surface_face_vtx_idx =  malloc (sizeof(int *) * n_part);
  int **surface_face_vtx =  malloc (sizeof(int *) * n_part);
  double **surface_coords = malloc (sizeof(double *) * n_part);

  PDM_g_num_t **surface_face_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surface_vtx_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  const PDM_g_num_t **surface_face_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  const PDM_g_num_t **surface_vtx_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  PDM_gen_gnum_t* gen_gnum_face = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);
  PDM_gen_gnum_t* gen_gnum_vtx  = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);


  if (i_rank == 0) {
    printf("-- mesh dist set\n");
    fflush(stdout);
  }

  double strip = 0.1;

  double **pts_coords = malloc (sizeof(double *) * n_part);
  double **char_length = malloc (sizeof(double *) * n_part);

  int n_pts   = ((n_vtx_seg * n_vtx_seg * n_vtx_seg) / n_rank) / n_part;
  int n_pts_x = (int) (n_pts/(1.+(1.-4.*strip)+(1.-4.*strip) * (1.-4.*strip)));
  int n_pts_y = (int) ((1.-4.*strip) * n_pts);
  int n_pts_z = (int) ((1.-4.*strip) * (1.-4.*strip) * n_pts);

  double _char_length = (xmax-xmin)/n_vtx_seg;

  if ((n_pts - ( n_pts_y + n_pts_z)) < 0)
    n_pts = n_pts_x + n_pts_y + n_pts_z;
  else
    n_pts_x = n_pts - ( n_pts_y + n_pts_z);

  for (int i_part = 0; i_part < n_part; i_part++) {
    pts_coords[i_part] = malloc (sizeof(double) * 3 * n_pts);
    char_length[i_part] = malloc (sizeof(double) * n_pts);

    int idx = 0;
    int idx2 = 0;
    for (int i = 0; i < n_pts_x/2; i++) {
      double x = _frand_a_b (xmin, xmin+strip);
      double y = _frand_a_b (ymin, ymax);
      double z = _frand_a_b (zmin, zmax);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

    for (int i =  n_pts_x/2; i < n_pts_x; i++) {
      double x = _frand_a_b (xmax-strip, xmax);
      double y = _frand_a_b (ymin, ymax);
      double z = _frand_a_b (zmin, zmax);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

    for (int i = 0; i < n_pts_y/2; i++) {
      double x = _frand_a_b (xmin+strip, xmax-strip);
      double y = _frand_a_b (ymin, ymin+strip);
      double z = _frand_a_b (zmin, zmax);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

    for (int i =  n_pts_y/2; i < n_pts_y; i++) {
      double x = _frand_a_b (xmin+strip, xmax-strip);
      double y = _frand_a_b (ymax-strip, ymax);
      double z = _frand_a_b (zmin, zmax);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

    for (int i = 0; i < n_pts_z/2; i++) {
      double x = _frand_a_b (xmin+strip, xmax-strip);
      double y = _frand_a_b (ymin+strip, ymax-strip);
      double z = _frand_a_b (zmin, zmin+strip);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

    for (int i =  n_pts_z/2; i < n_pts_z; i++) {
      double x = _frand_a_b (xmin+strip, xmax-strip);
      double y = _frand_a_b (ymin+strip, ymax-strip);
      double z = _frand_a_b (zmax-strip, zmax);
      char_length[i_part][idx2++] = _char_length;
      pts_coords[i_part][idx++] = x;
      pts_coords[i_part][idx++] = y;
      pts_coords[i_part][idx++] = z;
    }

  }

  PDM_gen_gnum_t* gen_gnum_pts = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_coords (gen_gnum_pts, i_part, n_pts, pts_coords[i_part], char_length[i_part]);
  }

  PDM_gnum_compute(gen_gnum_pts);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

    n_select_face[i_part] = 0;
    n_select_vtx[i_part] = 0;

    select_face[i_part] = malloc (sizeof(int) * n_face);

    for (int i = 0; i < n_face; i++) {
      select_face[i_part][i] = 0;
    }

    select_vtx[i_part] = malloc (sizeof(int) * n_vtx);

    for (int i = 0; i < n_vtx; i++) {
      select_vtx[i_part][i] = 0;
    }

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_bound_proc_idx,
                           &face_part_bound_part_idx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    int iii = 0;
    for (int i = 0; i < n_face; i++) {
      int icel2 = face_cell[2*i+1];
      if (icel2 == 0) {
        iii++;
        select_face[i_part][i] = 1;
      }
    }

    for (int i = 0; i < face_part_bound_proc_idx[n_rank]; i++) {
      select_face[i_part][face_part_bound[4*i]-1] = 0;
    }

    int idx = 1;
    int s_face_vtx = 0;
    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] == 1) {
        select_face[i_part][i] = idx;
        s_face_vtx += (face_vtx_idx[i+1] - face_vtx_idx[i]);
        idx += 1;
      }
    }
    n_select_face[i_part] = idx - 1;

    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] != 0) {
        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          select_vtx[i_part][face_vtx[j]-1] = 1;
        }
      }
    }

    idx = 1;
    for (int i = 0; i < n_vtx; i++) {
      if (select_vtx[i_part][i] == 1) {
        select_vtx[i_part][i] = idx;
        idx += 1;
      }
    }
    n_select_vtx[i_part] = idx - 1;

    surface_face_vtx_idx[i_part] = malloc (sizeof(int) * (n_select_face[i_part] + 1));
    surface_face_vtx_idx[i_part][0] = 0;
    surface_face_vtx[i_part] = malloc (sizeof(int) * s_face_vtx);

    surface_coords[i_part] = malloc (sizeof(double) * 3 * n_select_vtx[i_part]);

    surface_face_parent_gnum[i_part] =
      malloc (sizeof(PDM_g_num_t) * n_select_face[i_part]);
    surface_vtx_parent_gnum[i_part] =
      malloc (sizeof(PDM_g_num_t) * n_select_vtx[i_part]);

    surface_face_gnum[i_part] = NULL;
    surface_vtx_gnum[i_part] = NULL;

    idx = 0;
    int idx1 = 0;
    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] > 0) {
        surface_face_vtx_idx[i_part][idx+1] =
          surface_face_vtx_idx[i_part][idx] + (face_vtx_idx[i+1] - face_vtx_idx[i]);

        surface_face_parent_gnum[i_part][idx] = face_ln_to_gn[i];

        idx += 1;

        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          surface_face_vtx[i_part][idx1++] = select_vtx[i_part][face_vtx[j]-1];
        }
      }
    }

    idx = 0;
    for (int i = 0; i < n_vtx; i++) {
      if (select_vtx[i_part][i] > 0) {

        surface_vtx_parent_gnum[i_part][idx] = vtx_ln_to_gn[i];
        surface_coords[i_part][3*idx  ] = vtx[3*i];
        surface_coords[i_part][3*idx+1] = vtx[3*i+1];
        surface_coords[i_part][3*idx+2] = vtx[3*i+2];

        idx += 1;

      }
    }

    PDM_gnum_set_from_parents (gen_gnum_face,
                               i_part,
                               n_select_face[i_part],
                               surface_face_parent_gnum[i_part]);


    PDM_gnum_set_from_parents (gen_gnum_vtx,
                               i_part,
                               n_select_vtx[i_part],
                               surface_vtx_parent_gnum[i_part]);

    for (int i = 0; i <  n_select_vtx[i_part]; i++) {

    }

  }

  PDM_gnum_compute (gen_gnum_face);

  PDM_gnum_compute (gen_gnum_vtx);

  PDM_g_num_t n_g_face_loc = 0;
  PDM_g_num_t n_g_vtx_loc = 0;

  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_vtx = 0;

  for (int i_part = 0; i_part < n_part; i_part++) {
    surface_face_gnum[i_part] = PDM_gnum_get (gen_gnum_face, i_part);
    surface_vtx_gnum[i_part] = PDM_gnum_get (gen_gnum_vtx, i_part);

    for (int i = 0; i <  n_select_face[i_part]; i++) {
      n_g_face_loc = PDM_MAX(n_g_face_loc, surface_face_gnum[i_part][i]);
    }

    for (int i = 0; i <  n_select_vtx[i_part]; i++) {
      n_g_vtx_loc = PDM_MAX(n_g_vtx_loc, surface_vtx_gnum[i_part][i]);
    }

  }

  PDM_MPI_Allreduce (&n_g_face_loc, &n_g_face, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_MPI_Allreduce (&n_g_vtx_loc, &n_g_vtx, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 n_part);

  PDM_dist_cloud_surf_n_part_cloud_set (dist, 0, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                            i_part,
                                            n_select_face[i_part],
                                            surface_face_vtx_idx[i_part],
                                            surface_face_vtx[i_part],
                                            surface_face_gnum[i_part],
                                            n_select_vtx[i_part],
                                            surface_coords[i_part],
                                            surface_vtx_gnum[i_part]);

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_bound_proc_idx,
                           &face_part_bound_part_idx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    /* PDM_dist_cloud_surf_cloud_set (dist, */
    /*                          0, */
    /*                          i_part, */
    /*                          n_vtx, */
    /*                          vtx, */
    /*                          vtx_ln_to_gn); */

    PDM_g_num_t *pts_gnum =  PDM_gnum_get (gen_gnum_pts, i_part);

    PDM_dist_cloud_surf_cloud_set (dist,
                                   0,
                                   i_part,
                                   n_pts,
                                   pts_coords[i_part],
                                   pts_gnum);

    if (post) {
      char filename[999];
      sprintf(filename, "point_cloud_%d_%2.2d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                n_pts,
                                pts_coords[i_part],
                                pts_gnum,
                                NULL);

      for (int i = 0; i < n_pts; i++) {
        log_trace("point "PDM_FMT_G_NUM" : %f %f %f\n",
                  pts_gnum[i], pts_coords[i_part][3*i], pts_coords[i_part][3*i+1], pts_coords[i_part][3*i+2]);
      }
    }

  }

  if (i_rank == 0) {
    printf("-- Dist compute\n");
    fflush(stdout);
  }

  // PDM_dist_cloud_surf_compute (dist);
  PDM_dist_cloud_surf_compute (dist);

  if (i_rank == 0) {
    printf("-- Dist check\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_bound_proc_idx,
                           &face_part_bound_part_idx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    PDM_g_num_t *pts_gnum = PDM_gnum_get(gen_gnum_pts, i_part);

    int ierr = 0;
    for (int i = 0; i < n_pts; i++) {
      double d1 = PDM_MIN (PDM_ABS (pts_coords[i_part][3*i] - xmin), PDM_ABS (pts_coords[i_part][3*i] - xmax));
      double d2 = PDM_MIN (PDM_ABS (pts_coords[i_part][3*i+1]- ymin), PDM_ABS (pts_coords[i_part][3*i+1] - ymax));
      double d3 = PDM_MIN (PDM_ABS (pts_coords[i_part][3*i+2] - zmin), PDM_ABS (pts_coords[i_part][3*i+2] - zmax));
      double d = PDM_MIN (PDM_MIN (d1,d2), d3);
      d = d * d;
      if (PDM_ABS(distance[i] - d) > 1e-6) {
        log_trace("!!! pt "PDM_FMT_G_NUM" (%f %f %f) error dist = %e / %e, face "PDM_FMT_G_NUM"\n",
                  pts_gnum[i],
                  pts_coords[i_part][3*i], pts_coords[i_part][3*i+1], pts_coords[i_part][3*i+2],
                  distance[i], d, closest_elt_gnum[i]);
        ierr += 1;
        /* printf ("Erreur distance %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e\n", i, */
        /*         pts_coords[i_part][3*i], pts_coords[i_part][3*i+1], pts_coords[i_part][3*i+2], distance[i], d); */
      }
      /* else { */
      /*   printf ("ok distance %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e\n", i, */
      /*           pts_coords[i_part][3*i], pts_coords[i_part][3*i+1], pts_coords[i_part][3*i+2], distance[i], d); */
      /* } */
    }

    if (ierr > 0) {
      printf ("Erreur distance pour %d points\n", ierr);
      abort();
    }

    if (i_rank == 0) {
      printf ("elements surfaciques : "PDM_FMT_G_NUM"\n", 6*(n_vtx_seg-1)*(n_vtx_seg-1));
      printf ("nombre de points     : %d\n", n_pts*n_part*n_rank);
      fflush(stdout);
    }
  }

  PDM_part_free(ppart);

  PDM_dcube_gen_free(dcube);
  PDM_dist_cloud_surf_dump_times(dist);
  PDM_dist_cloud_surf_free (dist);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (select_face[i_part]);
    free (select_vtx[i_part]);

    free (surface_face_vtx_idx[i_part]);
    free (surface_face_vtx[i_part]);
    free (surface_coords[i_part]);

    free (surface_face_parent_gnum[i_part]);
    free (surface_vtx_parent_gnum[i_part]);

    free (char_length[i_part]);
    free (pts_coords[i_part]);

  }

  free (char_length);
  free (pts_coords);

  free (select_face);
  free (select_vtx);

  free (n_select_face);
  free (n_select_vtx);

  free (surface_face_vtx_idx);
  free (surface_face_vtx);
  free (surface_coords);

  free (surface_face_parent_gnum);
  free (surface_vtx_parent_gnum);

  free (surface_face_gnum);
  free (surface_vtx_gnum);

  PDM_gnum_free(gen_gnum_face);
  PDM_gnum_free(gen_gnum_vtx);
  PDM_gnum_free(gen_gnum_pts);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
