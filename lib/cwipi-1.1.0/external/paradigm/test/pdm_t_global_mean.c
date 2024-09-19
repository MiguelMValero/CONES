#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_global_mean.h"
#include "pdm_priv.h"

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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
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

  PDM_g_num_t        n_vtx_seg = 10;
  double             length  = 1.;
  int                n_part   = 1;
  int                post    = 0;
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
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
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

  if (0 == 1) {

    PDM_printf("[%i] n_face_group    : %i\n", i_rank, n_face_group);
    PDM_printf("[%i] dn_cell        : %i\n", i_rank, dn_cell);
    PDM_printf("[%i] dn_face        : %i\n", i_rank, dn_face);
    PDM_printf("[%i] dn_vtx         : %i\n", i_rank, dn_vtx);

    PDM_printf("[%i] dface_cell     : ", i_rank);
    for (int i = 0; i < 2 * dn_face; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_cell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx_idx   : ", i_rank);
    for (int i = 0; i < dn_face + 1; i++)
      PDM_printf(" %i", dface_vtx_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx      : ", i_rank);
    for (int i = 0; i < dface_vtx_idx[dn_face]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_vtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dvtx_coord     : ", i_rank);
    for (int i = 0; i < 3*dn_vtx; i++)
      PDM_printf(" %12.5e", dvtx_coord[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group_idx : ", i_rank);
    for (int i = 0; i < n_face_group + 1; i++)
      PDM_printf(" %i", dface_group_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group    : ", i_rank);
    for (int i = 0; i < dface_group_idx[n_face_group]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_group[i]);
    PDM_printf("\n");

  }
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

  PDM_global_point_mean_t* gmean = PDM_global_mean_create (n_part, comm);

  PDM_part_t *ppart = PDM_part_create(comm,
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

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_part_time_get(ppart,
                    &elapsed,
                    &cpu,
                    &cpu_user,
                    &cpu_sys);

  PDM_printf("[%i]   - elapsed total                    : %12.5e\n", i_rank, elapsed[0]);
  PDM_printf("[%i]   - elapsed building graph           : %12.5e\n", i_rank, elapsed[1]);
  PDM_printf("[%i]   - elapsed splitting graph          : %12.5e\n", i_rank, elapsed[2]);
  PDM_printf("[%i]   - elapsed building mesh partitions : %12.5e\n", i_rank, elapsed[3]);

  PDM_printf("[%i]   - cpu total                        : %12.5e\n", i_rank, cpu[0]);
  PDM_printf("[%i]   - cpu building graph               : %12.5e\n", i_rank, cpu[1]);
  PDM_printf("[%i]   - cpu splitting graph              : %12.5e\n", i_rank, cpu[2]);
  PDM_printf("[%i]   - cpu building mesh partitions     : %12.5e\n", i_rank, cpu[3]);

  PDM_printf("[%i]   - cpu_user total                   : %12.5e\n", i_rank, cpu_user[0]);
  PDM_printf("[%i]   - cpu_user building graph          : %12.5e\n", i_rank, cpu_user[1]);
  PDM_printf("[%i]   - cpu_user splitting graph         : %12.5e\n", i_rank, cpu_user[2]);
  PDM_printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", i_rank, cpu_user[3]);

  PDM_printf("[%i]   - cpu_sys total                    : %12.5e\n", i_rank, cpu_sys[0]);
  PDM_printf("[%i]   - cpu_sys building graph           : %12.5e\n", i_rank, cpu_sys[1]);
  PDM_printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", i_rank, cpu_sys[2]);
  PDM_printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", i_rank, cpu_sys[3]);

  struct timeval t_elaps_fin;
  gettimeofday(&t_elaps_fin, NULL);

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
                         (t_elaps_debut.tv_usec + 1000000 *
                          t_elaps_debut.tv_sec);
  long tranche_elapsed_max = tranche_elapsed;
  double t_elapsed = (double) tranche_elapsed_max/1000000.;

  PDM_printf("[%i]   - TEMPS DANS PART_CUBE  : %12.5e\n", i_rank,  t_elapsed);

  PDM_g_num_t **cellVtxGN = malloc (sizeof(PDM_g_num_t *) * n_part);

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
    int n_faceGroup2;

    PDM_part_part_dim_get(ppart,
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
                          &n_faceGroup2);

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

    PDM_part_part_val_get(ppart,
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

    int n_cellVtx = 0;
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j] - 1;
        n_cellVtx += face_vtx_idx[face+1] - face_vtx_idx[face];
      }
    }

    cellVtxGN[i_part] = malloc (sizeof(PDM_g_num_t) * n_cellVtx);

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
    int n_faceGroup2;

    PDM_part_part_dim_get(ppart,
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
                          &n_faceGroup2);

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

    PDM_part_part_val_get(ppart,
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

    int n_cellVtx = 0;
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j] - 1;
        for (int k = face_vtx_idx[face]; k < face_vtx_idx[face+1]; k++) {
          int iVtx = face_vtx[k] - 1;
          cellVtxGN[i_part][n_cellVtx++] = vtx_ln_to_gn[iVtx];
        }
      }
    }

    PDM_global_mean_set (gmean, i_part, n_cellVtx, cellVtxGN[i_part]);

  }

  double **local_field = malloc (sizeof(double *) * n_part);
  double **local_weight = malloc (sizeof(double *) * n_part);
  double **global_mean_field_ptr = malloc (sizeof(double *) * n_part);

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
    int n_faceGroup2;

    PDM_part_part_dim_get(ppart,
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
                          &n_faceGroup2);

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

    PDM_part_part_val_get(ppart,
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

    int n_cellVtx = 0;
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j] - 1;
        n_cellVtx += face_vtx_idx[face+1] - face_vtx_idx[face];
      }
    }

    local_field[i_part] = malloc (sizeof(double) * n_cellVtx * 3);
    local_weight[i_part] = malloc (sizeof(double) * n_cellVtx);
    global_mean_field_ptr[i_part] = malloc (sizeof(double) * n_cellVtx * 3);

    n_cellVtx = 0;
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j] - 1;
        for (int k = face_vtx_idx[face]; k < face_vtx_idx[face+1]; k++) {
          int iVtx = face_vtx[k] - 1;
          for (int l = 0; l < 3; l++) {
            local_field[i_part][3*n_cellVtx + l] = vtx[3*iVtx+l];
            global_mean_field_ptr[i_part][3*n_cellVtx + l] = 0;
          }
          local_weight[i_part][n_cellVtx] = 1.;
          n_cellVtx += 1;
        }
      }
    }

    PDM_global_mean_field_set (gmean, i_part, 3,
                               local_field[i_part],
                               local_weight[i_part],
                               global_mean_field_ptr[i_part]);

  }

  PDM_global_mean_field_compute (gmean);

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
    int n_faceGroup2;

    PDM_part_part_dim_get(ppart,
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
                          &n_faceGroup2);

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

    PDM_part_part_val_get(ppart,
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

    int n_cellVtx = 0;
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j] - 1;
        n_cellVtx += face_vtx_idx[face+1] - face_vtx_idx[face];
      }
    }

    for (int i = 0; i < 3 * n_cellVtx; i++) {
      if (PDM_ABS (local_field[i_part][i] - global_mean_field_ptr[i_part][i]) > 1e-5) {
        PDM_error (__FILE__, __LINE__, 0, "Error in global mean\n");
        abort();
      }
    }

    free (local_field[i_part]);
    free (local_weight[i_part]);
    free (global_mean_field_ptr[i_part]);
    free (cellVtxGN[i_part]);
  }

  free (local_field);
  free (local_weight);
  free (global_mean_field_ptr);
  free (cellVtxGN);
  free (dcell_part);

  PDM_part_free(ppart);

  PDM_global_mean_free (gmean);

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
