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
#include "pdm_geom_elem.h"

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
  double        length    = 1.;
  int           n_part    = 1;
  int           post      = 0;
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

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank    : %d\n", n_rank);
    PDM_printf ("  - n_vtx_seg : "PDM_FMT_G_NUM"\n", n_vtx_seg);
    PDM_printf ("  - length    : %f\n", length);
    PDM_printf ("  - method    : %d\n", method);
  }

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

  /* free (dface_cell); */
  /* free (dface_vtx_idx); */
  /* free (dface_vtx); */
  /* free (dvtx_coord); */
  /* free (dface_group_idx); */
  /* free (dface_group); */

  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t* dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            PDM_MPI_COMM_WORLD,
                                                            PDM_OWNERSHIP_KEEP);

  int **select_face   = malloc (sizeof(int *) * n_part);
  int  *n_select_face = malloc (sizeof(int  ) * n_part);
  int **select_vtx    = malloc (sizeof(int *) * n_part);
  int  *n_select_vtx  = malloc (sizeof(int  ) * n_part);

  int **surface_face_vtx_idx =  malloc (sizeof(int *) * n_part);
  int **surface_face_vtx =  malloc (sizeof(int *) * n_part);
  double **surface_coords = malloc (sizeof(double *) * n_part);

  PDM_g_num_t **surface_face_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surface_vtx_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  const PDM_g_num_t **surface_face_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  const PDM_g_num_t **surface_vtx_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  PDM_gen_gnum_t* gen_gnum_face = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);
  PDM_gen_gnum_t* gen_gnum_vtx  = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);

  double **cell_volume = malloc (sizeof(double *) * n_part);
  double **cell_center = malloc (sizeof(double *) * n_part);

  if (i_rank == 0) {
    printf("-- mesh dist set\n");
    fflush(stdout);
  }

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

    const int     isOriented = 0;
    cell_volume[i_part] = malloc(sizeof(double) * n_cell);
    cell_center[i_part] = malloc(sizeof(double) * 3 * n_cell);

    PDM_geom_elem_polyhedra_properties (isOriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume[i_part],
                                        cell_center[i_part],
                                        NULL,
                                        NULL);

    PDM_dist_cloud_surf_cloud_set (dist,
                                   0,
                                   i_part,
                                   n_cell,
                                   cell_center[i_part],
                                   cell_ln_to_gn);

    // for (int i = 0; i < n_cell; i++) {
    //   log_trace("point "PDM_FMT_G_NUM" : %f %f %f\n",
    //             cell_ln_to_gn[i],
    //             cell_center[i_part][3*i  ],
    //             cell_center[i_part][3*i+1],
    //             cell_center[i_part][3*i+2]);
    // }

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

    int ierr = 0;
    for (int i = 0; i < n_cell; i++) {
      double d1 = PDM_MIN (PDM_ABS (cell_center[i_part][3*i] - xmin),
                           PDM_ABS (cell_center[i_part][3*i] - xmax));
      double d2 = PDM_MIN (PDM_ABS (cell_center[i_part][3*i+1] - ymin),
                           PDM_ABS (cell_center[i_part][3*i+1] - ymax));
      double d3 = PDM_MIN (PDM_ABS (cell_center[i_part][3*i+2] - zmin),
                           PDM_ABS (cell_center[i_part][3*i+2] - zmax));
      double d = PDM_MIN (PDM_MIN (d1,d2), d3);
      d = d * d;
      if (PDM_ABS(distance[i] - d) > 1e-10) {
        ierr += 1;
        log_trace("!!! pt "PDM_FMT_G_NUM" (%f %f %f) error dist = %e / %e, face "PDM_FMT_G_NUM"\n",
                  cell_ln_to_gn[i],
                  cell_center[i_part][3*i],
                  cell_center[i_part][3*i+1],
                  cell_center[i_part][3*i+2],
                  distance[i], d, closest_elt_gnum[i]);
        printf ("Erreur distance "PDM_FMT_G_NUM" (%12.5e %12.5e %12.5e) : %12.5e %12.5e "PDM_FMT_G_NUM"\n",
                cell_ln_to_gn[i],
                cell_center[i_part][3*i],
                cell_center[i_part][3*i+1],
                cell_center[i_part][3*i+2],
                distance[i],
                d,
                closest_elt_gnum[i]);
      }
      /* else { */
      /*   if ((i+1 == 874) || */
      /*       (i+1 == 1266) || */
      /*       (i+1 == 1069) || */
      /*       (i+1 == 1071) || */
      /*       (i+1 == 1056) || */
      /*       (i+1 == 1084) || */
      /*       (i+1 == 1070)){ */

      /*     printf ("Affiche distance %d (%12.5e %12.5e %12.5e) : %12.5e %12.5e %ld\n", i+1, */
      /*             cell_center[i_part][3*i], */
      /*             cell_center[i_part][3*i+1], */
      /*             cell_center[i_part][3*i+2], */
      /*             distance[i], */
      /*             d, */
      /*             closest_elt_gnum[i]); */
      /*   } */
      /* } */
    }

    if (ierr > 0) {
      printf ("Erreur distance pour %d points\n", ierr);
      abort();
    }

    if (i_rank == 0) {
      printf ("elements surfaciques : "PDM_FMT_G_NUM"\n", 6*(n_vtx_seg-1)*(n_vtx_seg-1));
      printf ("nombre de points     : "PDM_FMT_G_NUM"\n", n_vtx_seg*n_vtx_seg*n_vtx_seg);
      fflush(stdout);
    }
  }


  if (post) {
    /* Prepare writer */
    PDM_writer_t *id_cs = PDM_writer_create ("Ensight",
                                             PDM_WRITER_FMT_ASCII,
                                             PDM_WRITER_TOPO_CST,
                                             PDM_WRITER_OFF,
                                             "test_dist",
                                             "dist",
                                             PDM_MPI_COMM_WORLD,
                                             PDM_IO_KIND_MPI_SIMPLE,
                                             1.,
                                             NULL);

    int id_geom = PDM_writer_geom_create (id_cs,
                                          "mesh",
                                          n_part);

    int id_var_dist = PDM_writer_var_create (id_cs,
                                             PDM_WRITER_OFF,
                                             PDM_WRITER_VAR_SCALAR,
                                             PDM_WRITER_VAR_ELEMENTS,
                                             "wall_dist");

    int id_var_closest = PDM_writer_var_create (id_cs,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "closest_bnd_face");

    PDM_writer_step_beg (id_cs, 0.);

    /* Write geometry */
    int **face_vtx_n  = malloc (sizeof(int *) * n_part);
    int **cell_face_n = malloc (sizeof(int *) * n_part);

    PDM_real_t **val_dist    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_closest = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

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

      int         *cell_tag;
      int         *cell_face_idx;
      int         *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int         *face_tag;
      int         *face_cell;
      int         *face_vtx_idx;
      int         *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int         *face_part_bound_proc_idx;
      int         *face_part_bound_part_idx;
      int         *face_part_bound;
      int         *vtx_tag;
      double      *vtx;
      PDM_g_num_t *vtx_ln_to_gn;
      int         *face_group_idx;
      int         *face_group;
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

      PDM_writer_geom_coord_set (id_cs,
                                 id_geom,
                                 i_part,
                                 n_vtx,
                                 vtx,
                                 vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);

      face_vtx_n[i_part] = malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++) {
        face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
      }

      cell_face_n[i_part] = malloc (sizeof(int) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                           id_geom,
                                           i_part,
                                           n_cell,
                                           n_face,
                                           face_vtx_idx,
                                           face_vtx_n[i_part],
                                           face_vtx,
                                           cell_face_idx,
                                           cell_face_n[i_part],
                                           cell_face,
                                           cell_ln_to_gn);

      val_dist[i_part]    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      val_closest[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        val_dist[i_part][i]    = (PDM_real_t) sqrt(distance[i]);
        val_closest[i_part][i] = (PDM_real_t) closest_elt_gnum[i];
      }
    }

    PDM_writer_geom_write (id_cs,
                           id_geom);

    // write variables
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_writer_var_set (id_cs,
                          id_var_dist,
                          id_geom,
                          i_part,
                          val_dist[i_part]);

      PDM_writer_var_set (id_cs,
                          id_var_closest,
                          id_geom,
                          i_part,
                          val_closest[i_part]);
    }

    PDM_writer_var_write (id_cs,
                          id_var_dist);
    PDM_writer_var_write (id_cs,
                          id_var_closest);

    PDM_writer_var_free (id_cs,
                         id_var_dist);
    PDM_writer_var_free (id_cs,
                         id_var_closest);

    PDM_writer_step_end (id_cs);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free (val_dist[i_part]);
      free (val_closest[i_part]);
      free (cell_face_n[i_part]);
      free (face_vtx_n[i_part]);
    }
    free (val_dist);
    free (val_closest);
    free (cell_face_n);
    free (face_vtx_n);

    PDM_writer_geom_data_free (id_cs, id_geom);
    PDM_writer_geom_free (id_cs, id_geom);
    PDM_writer_free (id_cs);
  }


  for (int i_part = 0; i_part < n_part; i_part++) {
    free (cell_center[i_part]);
    free (cell_volume[i_part]);
  }
  free (cell_center);
  free (cell_volume);

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

  }

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

  PDM_MPI_Finalize();

   if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
