#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_writer.h"
#include "pdm_vtk.h"
#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_part_extension.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_cellface_orient.h"

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
_read_args(int                 argc,
           char              **argv,
           PDM_g_num_t        *n_vtx_seg,
           double             *length,
           int                *n_part,
           int                *post,
           PDM_extend_type_t  *extend_type,
           int                *extension_depth,
           int                *method)
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
    else if (strcmp(argv[i], "-ext_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extend_type = (PDM_extend_type_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-from_face") == 0) {
      *extend_type = PDM_EXTEND_FROM_FACE;
    }
    else if (strcmp(argv[i], "-from_vtx") == 0) {
      *extend_type = PDM_EXTEND_FROM_VTX;
    }
    else if (strcmp(argv[i], "-ext_depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth = atoi(argv[i]);
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


static void
_visu
(
 const char      *name_chr,
 int              n_part,
 int              n_cell[],
 int              n_face[],
 int              n_vtx[],
 int             *cell_face_idx[],
 int             *cell_face[],
 PDM_g_num_t     *cell_ln_to_gn[],
 int             *face_vtx_idx[],
 int             *face_vtx[],
 PDM_g_num_t     *face_ln_to_gn[],
 double          *vtx[],
 PDM_g_num_t     *vtx_ln_to_gn[]
)
{
  PDM_UNUSED(face_ln_to_gn);

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  // for (int i = 0; i < n_part; i++) {
  //   char filename[999];
  //   sprintf(filename, "check_faces_part_%d_rank_%d.vtk", i, i_rank);
  //   PDM_vtk_write_polydata(filename,
  //                          n_vtx[i],
  //                          vtx[i],
  //                          vtx_ln_to_gn[i],
  //                          n_face[i],
  //                          face_vtx_idx[i],
  //                          face_vtx[i],
  //                          face_ln_to_gn[i],
  //                          NULL);
  // }


  PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_ASCII,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "test_part_extension",
                                          name_chr,
                                          PDM_MPI_COMM_WORLD,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

  /* Creation de la geometrie */

  int id_geom = PDM_writer_geom_create(id_cs,
                                       "cube_geom",
                                       n_part);

  int *n_part_procs = (int *) malloc(sizeof(int) * n_rank);

  PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                     (void *) n_part_procs, 1, PDM_MPI_INT,
                     PDM_MPI_COMM_WORLD);

  int *distrib_part = (int *) malloc(sizeof(int) * (n_rank + 1));

  distrib_part[0] = 0;
  for (int i = 0; i < n_rank; i++) {
    distrib_part[i+1] = distrib_part[i] + n_part_procs[i];
  }

  free(n_part_procs);

  /* Creation des variables */

  int id_var_num_part = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");


  int id_var_num_rank = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_rank");

  int id_var_cell_num = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "cell_num");

  int id_var_cell_gnum = PDM_writer_var_create(id_cs,
                                               PDM_WRITER_OFF,
                                               PDM_WRITER_VAR_SCALAR,
                                               PDM_WRITER_VAR_ELEMENTS,
                                               "cell_gnum");


  PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_num_rank = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_cell_num = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_cell_gnum = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
  int **face_vtx_n = (int **) malloc(sizeof(int *) * n_part);
  int **cell_face_n = (int **) malloc(sizeof(int *) * n_part);
  PDM_writer_step_beg (id_cs, 0.);
  for (int i_part = 0; i_part < n_part; i_part++) {

    face_vtx_n[i_part] = (int *) malloc(sizeof(int) * n_face[i_part]);
    cell_face_n[i_part] = (int *) malloc(sizeof(int) * n_cell[i_part]);

    for (int i = 0; i < n_face[i_part]; i++) {
      face_vtx_n[i_part][i] = face_vtx_idx[i_part][i+1] - face_vtx_idx[i_part][i];
    }

    for (int i = 0; i < n_cell[i_part]; i++) {
      cell_face_n[i_part][i] = cell_face_idx[i_part][i+1] - cell_face_idx[i_part][i];
    }

    PDM_writer_geom_coord_set(id_cs,
                              id_geom,
                              i_part,
                              n_vtx[i_part],
                              vtx[i_part],
                              vtx_ln_to_gn[i_part],
                              PDM_OWNERSHIP_USER);

    /* Construction de la connectivite pour sortie graphique */

    PDM_writer_geom_cell3d_cellface_add (id_cs,
                                         id_geom,
                                         i_part,
                                         n_cell       [i_part],
                                         n_face       [i_part],
                                         face_vtx_idx [i_part],
                                         face_vtx_n   [i_part],
                                         face_vtx     [i_part],
                                         cell_face_idx[i_part],
                                         cell_face_n  [i_part],
                                         cell_face    [i_part],
                                         cell_ln_to_gn[i_part]);

  }

  PDM_writer_geom_write(id_cs,
                        id_geom);

  for (int i_part = 0; i_part < n_part; i_part++) {

    val_num_part[i_part] = (double *) malloc(sizeof(double) * n_cell[i_part]);
    val_cell_num[i_part] = (double *) malloc(sizeof(double) * n_cell[i_part]);
    val_num_rank[i_part] = (double *) malloc(sizeof(double) * n_cell[i_part]);
    val_cell_gnum[i_part] = (double *) malloc(sizeof(double) * n_cell[i_part]);
    for (int i = 0; i < n_cell[i_part]; i++) {
      val_num_part[i_part][i] = i_part + distrib_part[i_rank];
      val_cell_num[i_part][i] = i ;
      val_num_rank[i_part][i] = i_rank;
      val_cell_gnum[i_part][i] = (double) cell_ln_to_gn[i_part][i];
    }
    // PDM_log_trace_array_double(val_num_part[i_part], n_cell[i_part], "val_num_part :: ");

    PDM_writer_var_set(id_cs,
                       id_var_num_part,
                       id_geom,
                       i_part,
                       val_num_part[i_part]);

    PDM_writer_var_set(id_cs,
                       id_var_num_rank,
                       id_geom,
                       i_part,
                       val_num_rank[i_part]);

    PDM_writer_var_set(id_cs,
                       id_var_cell_num,
                       id_geom,
                       i_part,
                       val_cell_num[i_part]);

    PDM_writer_var_set(id_cs,
                       id_var_cell_gnum,
                       id_geom,
                       i_part,
                       val_cell_gnum[i_part]);
  }

  PDM_writer_var_write(id_cs,
                       id_var_num_part);

  PDM_writer_var_write(id_cs,
                       id_var_num_rank);

  PDM_writer_var_write(id_cs,
                       id_var_cell_num);

  PDM_writer_var_write(id_cs,
                       id_var_cell_gnum);


  PDM_writer_var_free(id_cs,
                      id_var_num_part);

  PDM_writer_var_free(id_cs,
                      id_var_num_rank);

  PDM_writer_var_free(id_cs,
                      id_var_cell_num);

  PDM_writer_var_free(id_cs,
                      id_var_cell_gnum);

  PDM_writer_step_end(id_cs);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (face_vtx_n[i_part]);
    free (cell_face_n[i_part]);
    free (val_num_part[i_part]);
    free (val_cell_num[i_part]);
    free (val_num_rank[i_part]);
    free (val_cell_gnum[i_part]);
  }
  free (distrib_part);
  free (face_vtx_n);
  free (val_num_part);
  free (val_cell_num);
  free (val_cell_gnum);
  free (val_num_rank);
  free (cell_face_n);

  // PDM_writer_geom_free(id_cs,
  //                      id_geom);

  PDM_writer_free(id_cs);

}

/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[n] : 10,20
// @@@param[n_part] : 1,2
// @@@param[ext_type] : 0,2
// @@@args[part_kind] : -parmetis, -pt-scotch
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg       = 10;
  double             length          = 1.;
  int                n_part          = 1;
  int                post            = 0;
  PDM_extend_type_t  extend_type     = PDM_EXTEND_FROM_FACE;
  int                extension_depth = 1;
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
             &extend_type,
             &extension_depth,
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

  PDM_part_extension_t* part_ext = PDM_part_extension_create(1,
                                                             &n_part,
                                                             PDM_EXTEND_FROM_FACE,
                                                             1,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  int          *pn_cell        = malloc(n_part * sizeof(int          ));
  int          *pn_face        = malloc(n_part * sizeof(int          ));
  int          *pn_vtx         = malloc(n_part * sizeof(int          ));
  int         **pcell_face_idx = malloc(n_part * sizeof(int         *));
  int         **pcell_face     = malloc(n_part * sizeof(int         *));
  int         **pface_vtx_idx  = malloc(n_part * sizeof(int         *));
  int         **pface_vtx      = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **pcell_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));
  double      **pvtx           = malloc(n_part * sizeof(double      *));


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
                          &n_face_group2);

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t  *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t  *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t  *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t  *face_group_ln_to_gn;

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

    char *use_multipart_var;
    int use_multipart = 0;
    if (( use_multipart_var = getenv( "PDM_USE_MULTIPART" )) != NULL ) {
      use_multipart = atoi(use_multipart_var);
    }

    if (!use_multipart) {
      PDM_cellface_orient (n_cell,
                           n_face,
                           n_vtx,
                           vtx,
                           cell_face_idx,
                           cell_face,
                           face_cell,
                           face_vtx_idx,
                           face_vtx);
    }

    pn_cell       [i_part] = n_cell;
    pn_face       [i_part] = n_face;
    pn_vtx        [i_part] = n_vtx;
    pcell_face_idx[i_part] = cell_face_idx;
    pcell_face    [i_part] = cell_face;
    pface_vtx_idx [i_part] = face_vtx_idx;
    pface_vtx     [i_part] = face_vtx;
    pcell_ln_to_gn[i_part] = cell_ln_to_gn;
    pface_ln_to_gn[i_part] = face_ln_to_gn;
    pvtx_ln_to_gn [i_part] = vtx_ln_to_gn;
    pvtx          [i_part] = vtx;


    PDM_part_extension_set_part(part_ext, 0, i_part,
                                n_cell,
                                n_face,
                                n_face_part_bound,
                                n_face_group,
                                0,   // n_edge
                                n_vtx,
                                cell_face_idx,
                                cell_face,
                                face_cell,
                                NULL, // face_edge_idx
                                NULL, // face_edge
                                face_vtx_idx,
                                face_vtx,
                                NULL, //edge_vtx
                                face_group_idx,
                                face_group,
                                NULL, // face_join_idx
                                NULL, // face_join
                                face_part_bound_proc_idx,
                                face_part_bound_part_idx,
                                face_part_bound,
                                NULL, // vtx_part_bound_proc_idx
                                NULL, // vtx_part_bound_part_idx
                                NULL, // vtx_part_bound
                                cell_ln_to_gn,
                                face_ln_to_gn,
                                NULL, // edge_ln_to_gn
                                vtx_ln_to_gn,
                                face_group_ln_to_gn,
                                vtx);

    // PDM_log_trace_array_long (cell_ln_to_gn, n_cell, "cell_ln_to_gn::");

    if( 0 == 1) {
      PDM_printf("[%i] n_face_group     : %i\n", i_rank, n_face_group);
      PDM_printf("[%i] n_cell          : %i\n", i_rank, n_cell);
      PDM_printf("[%i] n_face          : %i\n", i_rank, n_face);
      PDM_printf("[%i] n_vtx           : %i\n", i_rank, n_vtx);
      PDM_printf("[%i] n_face_part_bound : %i\n", i_rank, n_face_part_bound);

      PDM_printf("[%i] cell_face     : ", i_rank);
      for (int i = 0; i < n_cell; i++) {
        for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
          PDM_printf(" %i", cell_face[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("\n");

      PDM_printf("[%i]  face_part_bound    : ", i_rank);
      for (int i = 0; i < 4 * n_face_part_bound; i++)
        PDM_printf(" "PDM_FMT_G_NUM, face_part_bound[i]);
      PDM_printf("\n");

      PDM_printf("[%i]  cell_ln_to_gn    : ", i_rank);
      for (int i = 0; i < n_cell; i++)
        PDM_printf(" "PDM_FMT_G_NUM, cell_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_cell     : ", i_rank);
      for (int i = 0; i < 2 * n_face; i++)
        PDM_printf(" %i", face_cell[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_vtx      : ", i_rank);
      for (int i = 0; i < n_face; i++) {
        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          PDM_printf(" %i", face_vtx[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("[%i]  face_ln_to_gn    : ", i_rank);
      for (int i = 0; i < n_face; i++)
        PDM_printf(" "PDM_FMT_G_NUM, face_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] vtx           : ", i_rank);
      for (int i = 0; i < 3 * n_vtx; i++)
        PDM_printf(" %12.5e", vtx[i]);
      PDM_printf("\n");

      PDM_printf("[%i] vtx_ln_to_gn     : ", i_rank);
      for (int i = 0; i <  n_vtx; i++)
        PDM_printf(" "PDM_FMT_G_NUM, vtx_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_group_idx : ", i_rank);
      for (int i = 0; i < n_face_group + 1; i++)
        PDM_printf(" %i", face_group_idx[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_group    : ", i_rank);
      for (int i = 0; i < n_face_group; i++) {
        for (int j = face_group_idx[i]; j < face_group_idx[i+1]; j++) {
          PDM_printf(" %i", face_group[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("[%i] face_group_ln_to_gn   : ", i_rank);
      for (int i = 0; i < n_face_group; i++) {
        for (int j = face_group_idx[i]; j < face_group_idx[i+1]; j++) {
          PDM_printf(" "PDM_FMT_G_NUM, face_group_ln_to_gn[j]);
        }
        PDM_printf("\n");
      }
    }
  }

  if(1 == 1) {

    _visu ("ini_mesh1",
           n_part,
           pn_cell,
           pn_face,
           pn_vtx,
           pcell_face_idx,
           pcell_face,
           pcell_ln_to_gn,
           pface_vtx_idx,
           pface_vtx,
           pface_ln_to_gn,
           pvtx,
           pvtx_ln_to_gn);
  }

  free(pn_cell       );
  free(pn_face       );
  free(pn_vtx        );
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );
  free(pvtx          );


  PDM_part_extension_compute(part_ext);

  for (int i_part = 0; i_part < n_part; i_part++) {

    double* vtx_coord_extended;
    PDM_g_num_t* border_vtx_ln_to_gn;
    PDM_g_num_t* border_cell_ln_to_gn;
    int n_vtx_extended  = PDM_part_extension_vtx_coord_get(part_ext, 0, i_part, &vtx_coord_extended);
    int n_vtx_extended2 = PDM_part_extension_ln_to_gn_get(part_ext, 0, i_part, PDM_MESH_ENTITY_VTX,  &border_vtx_ln_to_gn);
    int n_cell_extended = PDM_part_extension_ln_to_gn_get(part_ext, 0, i_part, PDM_MESH_ENTITY_CELL, &border_cell_ln_to_gn);
    assert(n_vtx_extended == n_vtx_extended2);
    if(0 == 1) {
      for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
        printf("[%i] vtx_coord_extended[%i] = %12.5e %12.5e %12.5e "PDM_FMT_G_NUM" \n", i_part, i_vtx, vtx_coord_extended[3*i_vtx], vtx_coord_extended[3*i_vtx+1], vtx_coord_extended[3*i_vtx+2], border_vtx_ln_to_gn[i_vtx]);
      }

      for(int i_cell = 0; i_cell < n_cell_extended; ++i_cell) {
        printf("[%i] border_cell_ln_to_gn[%i] = "PDM_FMT_G_NUM" \n", i_part, i_cell, border_cell_ln_to_gn[i_cell]);
      }
    }
  }

  PDM_part_extension_free(part_ext);


  free(dcell_part);
  PDM_part_free(ppart);

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
