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
#include "pdm_reader_gamma.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
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
           char         **filename,
           int           *visu)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



int main(int argc, char *argv[])
{
  /*
   *  Read args
   */
  char *filename = NULL;
  int   visu     = 0;
  int   n_part   = 1;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  _read_args(argc,
             argv,
             &filename,
             &visu);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"box.mesh";
    printf("No file specified -> exit \n");
    return 0;
  }

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                        filename,
                                                        0,
                                                        0);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
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
  PDM_DMesh_nodal_free(dmn);

  if (visu) {
    // PDM_dmesh_nodal_dump_vtk(dmn,
    //                          PDM_GEOMETRY_KIND_VOLUMIC,
    //                          "dmn_reader_gamma_vol_");





    PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                            PDM_WRITER_FMT_BIN,
                                            PDM_WRITER_TOPO_CST,
                                            PDM_WRITER_OFF,
                                            "reader_gamma",
                                            "reader_gamma",
                                            PDM_MPI_COMM_WORLD,
                                            PDM_IO_KIND_MPI_SIMPLE,
                                            1.,
                                            NULL);

    int id_geom = PDM_writer_geom_create(id_cs,
                                         "reader_gamma",
                                         n_part);

    int id_var_num_part = PDM_writer_var_create(id_cs,
                                                PDM_WRITER_ON,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "num_part");

    PDM_writer_step_beg(id_cs, 0.);

    int **face_vtx_n  = malloc(sizeof(int *) * n_part);
    int **cell_face_n = malloc(sizeof(int *) * n_part);

    int **pface_vtx_idx = malloc(sizeof(int *) * n_part);
    int **pface_vtx     = malloc(sizeof(int *) * n_part);

    PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t *) * n_part);

    int use_edge = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *cell_face_idx;
      int *cell_face;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_KEEP);

      int *face_vtx_idx;
      int *face_vtx;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       PDM_OWNERSHIP_KEEP);



      if (face_vtx == NULL) {
        use_edge = 1;

        int *face_edge;
        int *face_edge_idx;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_edge_idx,
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);

        int *edge_vtx;
        int *edge_vtx_idx;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &edge_vtx_idx,
                                            &edge_vtx,
                                            PDM_OWNERSHIP_KEEP);

        pface_vtx_idx[i_part] = face_edge_idx;
        PDM_compute_face_vtx_from_face_and_edge(n_face,
                                                face_edge_idx,
                                                face_edge,
                                                edge_vtx,
                                                &pface_vtx[i_part]);
      }
      else {
        pface_vtx_idx[i_part] = face_vtx_idx;
        pface_vtx    [i_part] = face_vtx;
      }

      double *vtx_coord;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      PDM_g_num_t *cell_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *face_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *vtx_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      face_vtx_n [i_part] = (int *) malloc(sizeof(int) * n_face);
      cell_face_n[i_part] = (int *) malloc(sizeof(int) * n_cell);

      for (int i = 0; i < n_cell; i++) {
        cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      for (int i = 0; i < n_face; i++) {
        face_vtx_n[i_part][i] = pface_vtx_idx[i_part][i+1] - pface_vtx_idx[i_part][i];
      }

      PDM_writer_geom_coord_set(id_cs,
                                id_geom,
                                i_part,
                                n_vtx,
                                vtx_coord,
                                vtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                           id_geom,
                                           i_part,
                                           n_cell,
                                           n_face,
                                           pface_vtx_idx[i_part],
                                           face_vtx_n   [i_part],
                                           pface_vtx    [i_part],
                                           cell_face_idx,
                                           cell_face_n  [i_part],
                                           cell_face,
                                           cell_ln_to_gn);

      val_num_part[i_part] = malloc(sizeof(PDM_real_t) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        val_num_part[i_part][i] = n_part*i_rank + i_part;
      }

      PDM_writer_var_set(id_cs,
                         id_var_num_part,
                         id_geom,
                         i_part,
                         val_num_part[i_part]);
    }

    PDM_writer_geom_write(id_cs,
                          id_geom);

    PDM_writer_var_write(id_cs,
                         id_var_num_part);

    PDM_writer_var_free(id_cs,
                        id_var_num_part);

    PDM_writer_step_end(id_cs);

    for (int i = 0; i < n_part; i++) {
      free(face_vtx_n[i]);
      free(cell_face_n[i]);
      if (use_edge) {
        free(pface_vtx[i]);
      }
      free(val_num_part[i]);
    }
    free(face_vtx_n);
    free(cell_face_n);
    free(pface_vtx_idx);
    free(pface_vtx);
    free(val_num_part);

    PDM_writer_free(id_cs);
  }


  PDM_multipart_free(mpart);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;
}
