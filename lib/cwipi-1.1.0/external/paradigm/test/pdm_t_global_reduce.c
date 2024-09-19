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
#include "pdm_global_reduce.h"
#include "pdm_priv.h"
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
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *n_vtx_seg,
           double       *length,
           int          *n_part,
           int          *verbose,
           int          *method)
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
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = 3;
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

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t n_vtx_seg = 10;
  double      length    = 1.;
  int         n_part    = 1;
  int         verbose   = 0;

  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &verbose,
             (int *) &method);

  /*
   *  Init
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  /*
   *  Create distributed cube
   */

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;
  int          dface_vtx_l;
  int          dface_group_l;

  PDM_dcube_t *dcube = PDM_dcube_gen_init(comm,
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
                        &dface_vtx_l,
                        &dface_group_l);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);


  /*
   *  Create mesh partitions
   */
  if (i_rank == 0) {
    printf("-- Mesh partitioning\n");
  }
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_t *ppart = PDM_part_create(comm,
                  (PDM_part_split_t)  method,
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




  int *part_n_vtx = malloc (sizeof(int) * n_part);
  PDM_g_num_t **part_vtx_ln_to_gn  = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **part_local_field   = malloc (sizeof(double *)      * n_part);
  double      **part_reduced_field = malloc (sizeof(double *)      * n_part);

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
    int n_face_group2;

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
                           &n_face_group2);

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
    double      *vtx_coord;
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
                           &vtx_coord,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    part_n_vtx[i_part] = n_vtx;
    part_vtx_ln_to_gn[i_part] = malloc (sizeof(PDM_g_num_t) * n_vtx);
    memcpy (part_vtx_ln_to_gn[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);

    part_local_field[i_part] = malloc (sizeof(double) * n_vtx * 3);
    for (int i = 0; i < n_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        part_local_field[i_part][3*i+j] = (i_rank + 1) * vtx_coord[3*i+j];
      }
    }

    part_reduced_field[i_part] = malloc (sizeof(double) * n_vtx * 3);
  }


  /*
   *  Create global reduction object
   */
  PDM_global_reduce_t *gre = PDM_global_reduce_create (n_part,
                                                       comm);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_global_reduce_g_num_set (gre,
                                 i_part,
                                 part_n_vtx[i_part],
                                 part_vtx_ln_to_gn[i_part]);
  }


  double t_start, t_end, elapsed, elapsed_max;

  PDM_MPI_Barrier (comm);

  /*
   *  Global min
   */
  if (i_rank == 0) {
    printf("-- Global min\n");
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_global_reduce_field_set (gre,
                                 i_part,
                                 3,
                                 part_local_field[i_part],
                                 part_reduced_field[i_part]);
  }
  PDM_global_reduce_operation_set (gre,
                                   PDM_REDUCE_OP_MIN);

  t_start = PDM_MPI_Wtime();
  PDM_global_reduce_field_compute (gre);
  t_end = PDM_MPI_Wtime();

  elapsed = t_end - t_start;
  PDM_MPI_Allreduce (&elapsed, &elapsed_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  if (i_rank == 0) {
    printf("  elapsed: %12.5es\n", elapsed_max);
  }

  if (verbose) {
    log_trace ("-- Global min --\n");
    for (int i_part = 0; i_part < n_part; i_part++) {

      log_trace ("  part %d:\n",i_part);

      for (int i = 0; i < part_n_vtx[i_part]; i++) {
        log_trace("    [%d] ("PDM_FMT_G_NUM"): %f %f %f --> %f %f %f\n",
                  i, part_vtx_ln_to_gn[i_part][i],
                  part_local_field[i_part][3*i  ],
                  part_local_field[i_part][3*i+1],
                  part_local_field[i_part][3*i+2],
                  part_reduced_field[i_part][3*i  ],
                  part_reduced_field[i_part][3*i+1],
                  part_reduced_field[i_part][3*i+2]);
      }
    }
  }



  /*
   *  Global max
   */
  if (i_rank == 0) {
    printf("-- Global max\n");
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_global_reduce_field_set (gre,
                                 i_part,
                                 3,
                                 part_local_field[i_part],
                                 part_reduced_field[i_part]);
  }
  PDM_global_reduce_operation_set (gre,
                                   PDM_REDUCE_OP_MAX);

  t_start = PDM_MPI_Wtime();
  PDM_global_reduce_field_compute (gre);
  t_end = PDM_MPI_Wtime();

  elapsed = t_end - t_start;
  PDM_MPI_Allreduce (&elapsed, &elapsed_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  if (i_rank == 0) {
    printf("  elapsed: %12.5es\n", elapsed_max);
  }

  if (verbose) {
    log_trace ("-- Global max --\n");
    for (int i_part = 0; i_part < n_part; i_part++) {

      log_trace ("  part %d:\n",i_part);

      for (int i = 0; i < part_n_vtx[i_part]; i++) {
        log_trace("    [%d] ("PDM_FMT_G_NUM"): %f %f %f --> %f %f %f\n",
                  i, part_vtx_ln_to_gn[i_part][i],
                  part_local_field[i_part][3*i  ],
                  part_local_field[i_part][3*i+1],
                  part_local_field[i_part][3*i+2],
                  part_reduced_field[i_part][3*i  ],
                  part_reduced_field[i_part][3*i+1],
                  part_reduced_field[i_part][3*i+2]);
      }
    }
  }


  /*
   *  Global sum
   */
  if (i_rank == 0) {
    printf("-- Global sum\n");
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_global_reduce_field_set (gre,
                                 i_part,
                                 3,
                                 part_local_field[i_part],
                                 part_reduced_field[i_part]);
  }
  PDM_global_reduce_operation_set (gre,
                                   PDM_REDUCE_OP_MAX);

  t_start = PDM_MPI_Wtime();
  PDM_global_reduce_field_compute (gre);
  t_end = PDM_MPI_Wtime();

  elapsed = t_end - t_start;
  PDM_MPI_Allreduce (&elapsed, &elapsed_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  if (i_rank == 0) {
    printf("  elapsed: %12.5es\n", elapsed_max);
  }

  if (verbose) {
    log_trace ("-- Global sum --\n");
    for (int i_part = 0; i_part < n_part; i_part++) {

      log_trace ("  part %d:\n",i_part);

      for (int i = 0; i < part_n_vtx[i_part]; i++) {
        log_trace("    [%d] ("PDM_FMT_G_NUM"): %f %f %f --> %f %f %f\n",
                  i, part_vtx_ln_to_gn[i_part][i],
                  part_local_field[i_part][3*i  ],
                  part_local_field[i_part][3*i+1],
                  part_local_field[i_part][3*i+2],
                  part_reduced_field[i_part][3*i  ],
                  part_reduced_field[i_part][3*i+1],
                  part_reduced_field[i_part][3*i+2]);
      }
    }
  }


  for (int i_part = 0; i_part < n_part; i_part++) {
    free (part_vtx_ln_to_gn[i_part]);
    free (part_local_field[i_part]);
    free (part_reduced_field[i_part]);
  }
  free (part_n_vtx);
  free (part_vtx_ln_to_gn);
  free (part_local_field);
  free (part_reduced_field);
  free (dcell_part);

  PDM_part_free (ppart);

  PDM_global_reduce_free (gre);

  PDM_dcube_gen_free (dcube);

  if (i_rank == 0) {
    printf("-- End\n");
  }

  PDM_MPI_Finalize();

  return 0;
}
