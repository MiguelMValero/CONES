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
#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "pdm_part.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_dmesh.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_geom_elem.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_dcube_gen.h"

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
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *nx,
           PDM_g_num_t  *ny,
           PDM_g_num_t  *nz,
           double       *lengthx,
           double       *lengthy,
           double       *lengthz,
           int          *n_part,
           int          *randomize,
           int          *random_seed,
           int          *post,
           int          *method,
           int          *use_multipart)
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
        *lengthy = atof(argv[i]);
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ly") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthy = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthz = atof(argv[i]);
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
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *random_seed = atoi(argv[i]);
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
    else if (strcmp(argv[i], "-multipart") == 0) {
      *use_multipart = 1;
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
  PDM_g_num_t nx = 10;
  PDM_g_num_t ny = 10;
  PDM_g_num_t nz = 10;

  double lengthx = 1.;
  double lengthy = 1.;
  double lengthz = 1.;

  int n_part      = 1;
  int post        = 0;
  int randomize   = 0;
  int random_seed = 0;

  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  int use_multipart = 0;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nx,
              &ny,
              &nz,
              &lengthx,
              &lengthy,
              &lengthz,
              &n_part,
              &randomize,
              &random_seed,
              &post,
              (int *) &method,
              &use_multipart);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  /*
   *  Create distributed mesh
   */
  double xmin = 0.;
  double ymin = 0.;
  double zmin = 0.;

  PDM_g_num_t  ng_cell      = 0;
  PDM_g_num_t  ng_face      = 0;
  PDM_g_num_t  ng_vtx       = 0;
  int          n_face_group = 0;
  int          dn_cell      = 0;
  int          dn_face      = 0;
  int          dn_vtx       = 0;
  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (1) {
    PDM_dcube_t *dcube = PDM_dcube_gen_init(comm,
                                            nx,
                                            lengthx,
                                            xmin,
                                            ymin,
                                            zmin,
                                            PDM_OWNERSHIP_USER);

    int sface_vtx;
    int sface_group;
    PDM_dcube_gen_dim_get(dcube,
                          &n_face_group,
                          &dn_cell,
                          &dn_face,
                          &dn_vtx,
                          &sface_vtx,
                          &sface_group);

    PDM_dcube_gen_data_get(dcube,
                           &dface_cell,
                           &dface_vtx_idx,
                           &dface_vtx,
                           &dvtx_coord,
                           &dface_group_idx,
                           &dface_group);

    PDM_dcube_gen_free(dcube);
  }
  else {
    PDM_poly_vol_gen (comm,
                      xmin,
                      ymin,
                      zmin,
                      lengthx,
                      lengthy,
                      lengthz,
                      nx,
                      ny,
                      nz,
                      randomize,
                      random_seed,
                      &ng_cell,
                      &ng_face,
                      &ng_vtx,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dcell_face_idx,
                      &dcell_face,
                      &dface_cell,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);
  }

  if (n_rank == 1) {
    if (dcell_face != NULL) {
      PDM_log_trace_connectivity_long(dcell_face_idx,
                                      dcell_face,
                                      dn_cell,
                                      "dcell_face : ");
    }

    if (dface_cell != NULL) {
      int *dface_cell_idx = PDM_array_new_idx_from_const_stride_int(2, dn_face);
      PDM_log_trace_connectivity_long(dface_cell_idx,
                                      dface_cell,
                                      dn_cell,
                                      "dface_cell : ");
      free(dface_cell_idx);
    }

    int *_face_vtx = malloc(sizeof(int) * dface_vtx_idx[dn_face]);
    for (int i = 0; i < dface_vtx_idx[dn_face]; i++) {
      _face_vtx[i] = (int) dface_vtx[i];
    }



    PDM_vtk_write_polydata("polyvol_gen_faces.vtk",
                           dn_vtx,
                           dvtx_coord,
                           NULL,
                           dn_face,
                           dface_vtx_idx,
                           _face_vtx,
                           NULL,
                           NULL);
    free(_face_vtx);
  }


  if (i_rank == 0) printf("ng_cell = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM", ng_vtx = "PDM_FMT_G_NUM"\n", ng_cell, ng_face, ng_vtx);

  /*
   *  Create mesh partitions
   */
  PDM_multipart_t *mpart = NULL;
  PDM_part_t      *ppart = NULL;
  PDM_dmesh_t     *dmesh = NULL;

  if (use_multipart) {
    /* Initialize multipart */
    mpart = PDM_multipart_create(1, &n_part, PDM_FALSE,
				 method, PDM_PART_SIZE_HOMOGENEOUS,
				 NULL, comm, PDM_OWNERSHIP_KEEP);

    /* Generate dmesh */
    dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                              dn_cell,
                              dn_face,
                              0, // dn_edge
                              dn_vtx,
                              comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_USER);

    PDM_multipart_dmesh_set (mpart, 0, dmesh);

    /* Run */
    PDM_multipart_compute (mpart);

  }

  else {
    int have_dcell_part = 0;

    int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
    int *renum_properties_cell = NULL;
    int *renum_properties_face = NULL;
    int n_property_cell = 0;
    int n_property_face = 0;

    // id_part = 0;

    ppart = PDM_part_create(PDM_MPI_COMM_WORLD,
                            (PDM_part_split_t) method,
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
  }


  if (post) {
    /* Prepare writer */
    PDM_writer_t *id_cs = PDM_writer_create ("Ensight",
                                   PDM_WRITER_FMT_ASCII,
                                   PDM_WRITER_TOPO_CST,
                                   PDM_WRITER_OFF,
                                   "test_polyvol",
                                   "polyvol",
                                   PDM_MPI_COMM_WORLD,
                                   PDM_IO_KIND_MPI_SIMPLE,
                                   1.,
                                   NULL);

    int id_geom = PDM_writer_geom_create (id_cs,
                                          "mesh",
                                          n_part);

    // Cell local id
    int id_var_cell_g_num = PDM_writer_var_create (id_cs,
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAR,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "cell_g_num");

    int id_var_num_part = PDM_writer_var_create (id_cs,
                                                 PDM_WRITER_OFF,
                                                 PDM_WRITER_VAR_SCALAR,
                                                 PDM_WRITER_VAR_ELEMENTS,
                                                 "num_part");

    int id_var_vtx_g_num = PDM_writer_var_create (id_cs,
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAR,
                                                  PDM_WRITER_VAR_VERTICES,
                                                  "vtx_g_num");

    int id_var_coo_x = PDM_writer_var_create (id_cs,
                                              PDM_WRITER_ON,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_VERTICES,
                                              "coo_x");

    int id_var_coo_xyz = PDM_writer_var_create (id_cs,
                                                PDM_WRITER_ON,
                                                PDM_WRITER_VAR_VECTOR,
                                                PDM_WRITER_VAR_VERTICES,
                                                "coo_xyz");

    PDM_writer_step_beg (id_cs, 0.);

    /* Write geometry */
    int **face_vtx_n  = malloc (sizeof(int *) * n_part);
    int **cell_face_n = malloc (sizeof(int *) * n_part);

    PDM_real_t **val_cell_g_num = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_num_part   = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_vtx_g_num  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_x      = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_xyz    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_t_part;
      int s_cell_face;
      int s_face_vtx;
      int s_face_group;
      int *cell_tag      = NULL;
      int *face_tag      = NULL;
      int *vtx_tag       = NULL;
      int *cell_face     = NULL;
      int *cell_face_idx = NULL;
      int *face_cell     = NULL;
      int *face_cell_idx = NULL;
      int *face_vtx      = NULL;
      int *face_vtx_idx  = NULL;
      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      PDM_g_num_t* cell_ln_to_gn = NULL;
      PDM_g_num_t* face_ln_to_gn = NULL;
      PDM_g_num_t* vtx_ln_to_gn  = NULL;

      int          pn_face_group       = 0;
      int         *face_group_idx      = NULL;
      int         *face_group          = NULL;
      PDM_g_num_t *face_group_ln_to_gn = NULL;
      int *face_part_bound_proc_idx = NULL;
      int *face_part_bound_part_idx = NULL;
      int *face_part_bound          = NULL;

      double  *vtx = NULL;

      if (use_multipart) {
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        &cell_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);


        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                            &face_cell_idx,
                                            &face_cell,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_cell_idx == NULL);

        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &face_vtx_idx,
                                            &face_vtx,
                                            PDM_OWNERSHIP_KEEP);

        n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge_idx,
                                                     &face_edge,
                                                     PDM_OWNERSHIP_KEEP);

        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                0,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);
        PDM_multipart_group_get(mpart,
                                0,
                                i_part,
                                PDM_MESH_ENTITY_FACE,
                                &pn_face_group,
                                &face_group_idx,
                                &face_group,
                                &face_group_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);

        PDM_multipart_part_graph_comm_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_FACE,
                                          &face_part_bound_proc_idx,
                                          &face_part_bound_part_idx,
                                          &face_part_bound,
                                          PDM_OWNERSHIP_KEEP);
        PDM_multipart_part_vtx_coord_get(mpart,
                                         0,
                                         i_part,
                                         &vtx,
                                         PDM_OWNERSHIP_KEEP);

      }

      else {
        int n_edge_group2 = 0;
        PDM_part_part_dim_get (ppart,
                               i_part,
                               &n_cell,
                               &n_face,
                               &n_face_part_bound,
                               &n_vtx,
                               &n_proc,
                               &n_t_part,
                               &s_cell_face,
                               &s_face_vtx,
                               &s_face_group,
                               &n_edge_group2);

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
      }

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

      val_cell_g_num[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      val_num_part[i_part]   = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        val_cell_g_num[i_part][i] = (PDM_real_t) cell_ln_to_gn[i];
        val_num_part[i_part][i] = (PDM_real_t) (i_rank*n_part + i_part);
      }

      val_vtx_g_num[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx);
      val_coo_x[i_part]     = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx);
      val_coo_xyz[i_part]   = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_vtx * 3);
      for (int i = 0; i < n_vtx; i++) {
        val_vtx_g_num[i_part][i]   = (PDM_real_t) vtx_ln_to_gn[i];
        val_coo_x[i_part][i]       = vtx[3*i];
        val_coo_xyz[i_part][3*i  ] = vtx[3*i  ];
        val_coo_xyz[i_part][3*i+1] = vtx[3*i+1];
        val_coo_xyz[i_part][3*i+2] = vtx[3*i+2];
      }
    }

    PDM_writer_geom_write (id_cs,
                           id_geom);

    // write variables
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_writer_var_set (id_cs,
                          id_var_cell_g_num,
                          id_geom,
                          i_part,
                          val_cell_g_num[i_part]);

      PDM_writer_var_set (id_cs,
                          id_var_num_part,
                          id_geom,
                          i_part,
                          val_num_part[i_part]);

      PDM_writer_var_set (id_cs,
                          id_var_vtx_g_num,
                          id_geom,
                          i_part,
                          val_vtx_g_num[i_part]);
    }

    PDM_writer_var_write (id_cs,
                          id_var_cell_g_num);
    PDM_writer_var_write (id_cs,
                          id_var_num_part);
    PDM_writer_var_write (id_cs,
                          id_var_vtx_g_num);

    PDM_writer_var_free (id_cs,
                         id_var_cell_g_num);
    PDM_writer_var_free (id_cs,
                         id_var_num_part);
    PDM_writer_var_free (id_cs,
                         id_var_vtx_g_num);

    for (int nstep = 0; nstep < 10; nstep++) {

      double tstep = nstep * 0.01;

      if (nstep > 0) {
        PDM_writer_step_beg(id_cs, tstep);
      }

      for (int i_part = 0; i_part < n_part; i_part++) {
        PDM_writer_var_set (id_cs,
                            id_var_coo_x,
                            id_geom,
                            i_part,
                            val_coo_x[i_part]);
        PDM_writer_var_set (id_cs,
                            id_var_coo_xyz,
                            id_geom,
                            i_part,
                            val_coo_xyz[i_part]);
      }

      PDM_writer_var_write (id_cs,
                            id_var_coo_x);
      PDM_writer_var_write (id_cs,
                            id_var_coo_xyz);

      PDM_writer_var_data_free (id_cs,
                                id_var_coo_x);
      PDM_writer_var_data_free (id_cs,
                                id_var_coo_xyz);

      PDM_writer_step_end (id_cs);
    }


    for (int i_part = 0; i_part < n_part; i_part++) {
      free (val_cell_g_num[i_part]);
      free (val_num_part[i_part]);
      free (val_vtx_g_num[i_part]);
      free (val_coo_x[i_part]);
      free (val_coo_xyz[i_part]);
      free (cell_face_n[i_part]);
      free (face_vtx_n[i_part]);
    }
    free (val_cell_g_num);
    free (val_num_part);
    free (val_vtx_g_num);
    free (val_coo_x);
    free (val_coo_xyz);
    free (cell_face_n);
    free (face_vtx_n);

    PDM_writer_var_free (id_cs,
                         id_var_coo_x);
    PDM_writer_var_free (id_cs,
                         id_var_coo_xyz);

    PDM_writer_geom_data_free (id_cs, id_geom);
    PDM_writer_geom_free (id_cs, id_geom);
    PDM_writer_free (id_cs);
  }

  /*
   *  Finalize
   */
  free (dcell_face_idx);
  free (dcell_face);
  free (dface_cell);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dvtx_coord);
  free (dface_group_idx);
  free (dface_group);

  if (use_multipart) {
    PDM_multipart_free (mpart);
    PDM_dmesh_free (dmesh);
  }
  else {
    PDM_part_free (ppart);
  }

  PDM_MPI_Finalize();

  return 0;
}
