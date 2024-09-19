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

#include "pdm_writer.h"
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

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t  n_vtx_seg = 10;
  double        length   = 1.;
  int           n_part   = 1;
  int           post     = 0;

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

  PDM_dcube_t* dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD,
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
  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
//                  "PDM_PART_RENUM_CELL_CUTHILL",
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_t *ppart = PDM_part_create(PDM_MPI_COMM_WORLD,
                                      method,
                                      "PDM_PART_RENUM_CELL_CUTHILL",
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

  PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_ASCII,
                                          PDM_WRITER_TOPO_VARIABLE,
                                          PDM_WRITER_OFF,
                                          "test_3d_ens_mv_mesh",
                                          "chrd3d",
                                          PDM_MPI_COMM_WORLD,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);
//                                          "append = 1");

  // PDM_writer_t *id_cs = PDM_writer_create("Ensight",
  //                                         PDM_WRITER_FMT_ASCII,
  //                                         PDM_WRITER_TOPO_VARIABLE,
  //                                         PDM_WRITER_OFF,
  //                                         "test_3d_ens_mv_mesh",
  //                                         "chrd3d",
  //                                         PDM_MPI_COMM_WORLD,
  //                                         PDM_IO_KIND_MPI_SIMPLE,
  //                                         1.,
  //                                         NULL);

  int id_var = PDM_writer_cst_global_var_create (id_cs, "test_var_constante", -1.2345);

  /* Creation de la geometrie */

  int id_geom = PDM_writer_geom_create(id_cs,
                                       "test3d_geom",
                                       n_part);

  /* Creation des variables */

  int id_var_num_part = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_ON,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_coo_x = PDM_writer_var_create(id_cs,
                                           PDM_WRITER_ON,
                                           PDM_WRITER_VAR_SCALAR,
                                           PDM_WRITER_VAR_VERTICES,
                                           "coo_x");

  int id_var_coo_xyz = PDM_writer_var_create(id_cs,
                                             PDM_WRITER_ON,
                                             PDM_WRITER_VAR_VECTOR,
                                             PDM_WRITER_VAR_VERTICES,
                                             "coo_xyz");

  /* Debut d'ecritures */

  int **face_vtxNb = (int **) malloc(sizeof(int *) * n_part);
  int **cell_faceNb = (int **) malloc(sizeof(int *) * n_part);

  PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_coo_x    = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

  int *n_partProcs = (int *) malloc(sizeof(int) * numProcs);

  PDM_MPI_Allgather((void *) &n_part,     1, PDM_MPI_INT,
                    (void *) n_partProcs, 1, PDM_MPI_INT,
                    PDM_MPI_COMM_WORLD);

  int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

  debPartProcs[0] = 0;
  for (int i = 0; i < numProcs; i++) {
    debPartProcs[i+1] = debPartProcs[i] + n_partProcs[i];
  }

  free(n_partProcs);

  /*
   *  Creation des variables :
   *   - numero de partition
   *   - scalaire
   *   - vecteur
   *   - tenseur
   */

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
                          &nEdgeGroup2);

    val_num_part[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
    val_coo_x[i_part]    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx);
    val_coo_xyz[i_part]  = (PDM_real_t *) malloc(sizeof(PDM_real_t) * 3 * n_vtx);
  }

  for (int nstep = 0; nstep < 10; nstep++) {

    double tstep = nstep * 0.01;

    PDM_writer_step_beg(id_cs, tstep);


    PDM_writer_cst_global_var_set (id_cs, id_var, tstep);

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

      assert(sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

      face_vtxNb[i_part] = (int *) malloc(sizeof(int) * n_face);
      cell_faceNb[i_part] = (int *) malloc(sizeof(int) * n_cell);

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

      for (int i = 0; i < n_cell; i++) {
        cell_faceNb[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      for (int i = 0; i < n_face; i++) {
        face_vtxNb[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
      }

      nsom_part[i_part]    = n_vtx;

      for (int i = 0; i < n_cell; i++) {
        val_num_part[i_part][i] = i_part + 1 + debPartProcs[i_rank];
      }

      for (int i = 0; i < n_vtx; i++) {
        val_coo_x[i_part][i]       = vtx[3*i];
        val_coo_xyz[i_part][3*i  ] = vtx[3*i  ];
        val_coo_xyz[i_part][3*i+1] = vtx[3*i+1];
        val_coo_xyz[i_part][3*i+2] = vtx[3*i+2];
      }

      PDM_writer_geom_coord_set(id_cs,
                                id_geom,
                                i_part,
                                n_vtx,
                                vtx,
                                vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);



      /* Construction de la connectivite pour sortie graphique */

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                           id_geom,
                                           i_part,
                                           n_cell,
                                           n_face,
                                           face_vtx_idx,
                                           face_vtxNb[i_part],
                                           face_vtx,
                                           cell_face_idx,
                                           cell_faceNb[i_part],
                                           cell_face,
                                           cell_ln_to_gn);

    }

    PDM_writer_geom_write(id_cs,
                          id_geom);

    for (int i_part = 0; i_part < n_part; i_part++) {

      free(face_vtxNb[i_part]);
      free(cell_faceNb[i_part]);

      PDM_writer_var_set(id_cs,
                         id_var_num_part,
                         id_geom,
                         i_part,
                         val_num_part[i_part]);
    }

    PDM_writer_var_write(id_cs,
                         id_var_num_part);

    for (int i_part = 0; i_part < n_part; i_part++) {

      for (int i = 0; i < nsom_part[i_part]; i++) {
        val_coo_x[i_part][i]       = val_coo_x[i_part][i]       + val_coo_x[i_part][i]/length;
        val_coo_xyz[i_part][3*i]   = val_coo_xyz[i_part][3*i]   + val_coo_xyz[i_part][3*i]/length;
        val_coo_xyz[i_part][3*i+1] = val_coo_xyz[i_part][3*i+1] + val_coo_xyz[i_part][3*i+1]/length;
        val_coo_xyz[i_part][3*i+2] = val_coo_xyz[i_part][3*i+2] + val_coo_xyz[i_part][3*i+2]/length;
      }

      PDM_writer_var_set(id_cs,
                         id_var_coo_x,
                         id_geom,
                         i_part,
                         val_coo_x[i_part]);

      PDM_writer_var_set(id_cs,
                         id_var_coo_xyz,
                         id_geom,
                         i_part,
                         val_coo_xyz[i_part]);

    }

    PDM_writer_var_write(id_cs,
                         id_var_coo_x);

    PDM_writer_var_write(id_cs,
                         id_var_coo_xyz);

    PDM_writer_var_data_free(id_cs,
                             id_var_coo_x);

    PDM_writer_var_data_free(id_cs,
                             id_var_coo_xyz);

    PDM_writer_var_data_free(id_cs,
                             id_var_num_part);

    PDM_writer_geom_data_reset (id_cs, id_geom);

    PDM_writer_step_end(id_cs);
  }

  free(debPartProcs);

  PDM_writer_var_free(id_cs,
                      id_var_num_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(val_num_part[i_part]);
  }
  free(val_num_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    //    free(cell_faceNb[i_part]);
    //free(face_vtxNb[i_part]);
    free(val_coo_x[i_part]);
    free(val_coo_xyz[i_part]);
  }
  free(val_coo_x);
  free(val_coo_xyz);
  free(cell_faceNb);
  free(face_vtxNb);
  free(nsom_part);

  PDM_writer_var_free(id_cs,
                      id_var_coo_x);

  PDM_writer_var_free(id_cs,
                      id_var_coo_xyz);

  /* Liberation memoire */

  PDM_writer_geom_data_free(id_cs,
                            id_geom);

  PDM_writer_geom_free(id_cs,
                       id_geom);

  PDM_writer_free(id_cs);

  /* Calculs statistiques */

  int    cells_average;
  int    cells_median;
  double cells_std_deviation;
  int    cells_min;
  int    cells_max;
  int    bound_part_faces_average;
  int    bound_part_faces_median;
  double bound_part_faces_std_deviation;
  int    bound_part_faces_min;
  int    bound_part_faces_max;
  int    bound_part_faces_sum;

  PDM_part_stat_get(ppart,
                    &cells_average,
                    &cells_median,
                    &cells_std_deviation,
                    &cells_min,
                    &cells_max,
                    &bound_part_faces_average,
                    &bound_part_faces_median,
                    &bound_part_faces_std_deviation,
                    &bound_part_faces_min,
                    &bound_part_faces_max,
                    &bound_part_faces_sum);

  if (i_rank == 0) {
    PDM_printf("Statistics :\n");
    PDM_printf("  - Number of cells :\n");
    PDM_printf("       * average            : %i\n", cells_average);
    PDM_printf("       * median             : %i\n", cells_median);
    PDM_printf("       * standard deviation : %12.5e\n", cells_std_deviation);
    PDM_printf("       * min                : %i\n", cells_min);
    PDM_printf("       * max                : %i\n", cells_max);
    PDM_printf("  - Number of faces exchanging with another partition :\n");
    PDM_printf("       * average            : %i\n", bound_part_faces_average);
    PDM_printf("       * median             : %i\n", bound_part_faces_median);
    PDM_printf("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
    PDM_printf("       * min                : %i\n", bound_part_faces_min);
    PDM_printf("       * max                : %i\n", bound_part_faces_max);
    PDM_printf("       * total              : %i\n", bound_part_faces_sum);
  }

  PDM_part_free(ppart);

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
