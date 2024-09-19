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
#include "pdm_vtk.h"
#include "pdm_multipart.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_block_to_part.h"
#include "pdm_writer.h"

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
    filename = (char *) PDM_MESH_DIR"bunny1k.vtk";
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


  int            n_vtx_field;
  char         **vtx_field_name;
  int           *vtx_field_stride;
  PDM_data_t    *vtx_field_type;
  void         **dvtx_field_value;
  int            n_elt_field;
  char         **elt_field_name;
  int           *elt_field_stride;
  PDM_data_t    *elt_field_type;
  void         **delt_field_value;

  PDM_dmesh_nodal_t *dmn = PDM_vtk_read_to_dmesh_nodal(comm,
                                                       filename,
                                                       &n_vtx_field,
                                                       &vtx_field_name,
                                                       &vtx_field_type,
                                                       &vtx_field_stride,
                                                       &dvtx_field_value,
                                                       &n_elt_field,
                                                       &elt_field_name,
                                                       &elt_field_type,
                                                       &elt_field_stride,
                                                       &delt_field_value);



  if (visu) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_RIDGE,    "reader_vtk_dmn_ridge_");
    if (dmn->mesh_dimension >= 2) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "reader_vtk_dmn_surface_");
    }
    if (dmn->mesh_dimension == 3) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC,  "reader_vtk_dmn_volume_");
    }
  }


  PDM_multipart_t *mpart = PDM_multipart_create(1,
                                                &n_part,
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


  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);


  if (visu) {
    // PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE,    "reader_vtk_pmn_ridge_");
    if (dmn->mesh_dimension >= 2) {
      PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "reader_vtk_pmn_surface_");
    }
    if (dmn->mesh_dimension == 3) {
      PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_VOLUMIC,  "reader_vtk_pmn_volume_");
    }
  }



  if (n_vtx_field > 0) {
    int          *pn_vtx        = malloc(sizeof(int        ) * n_part);
    PDM_g_num_t **pvtx_ln_to_gn = malloc(sizeof(PDM_g_num_t) * n_part);
    for (int i = 0; i < n_part; i++) {
      pn_vtx[i] = PDM_part_mesh_nodal_n_vtx_get(pmn, i);
      pvtx_ln_to_gn[i] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i);
    }

    PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_vtx,
                                 (const PDM_g_num_t **) pvtx_ln_to_gn,
                                                        pn_vtx,
                                                        n_part,
                                                        comm);

    void ***pvtx_field = malloc(sizeof(void *) * n_vtx_field);

    for (int i = 0; i < n_vtx_field; i++) {

      size_t s_data = 0;
      if (vtx_field_type[i] == PDM_INT) {
        s_data = sizeof(int);
      }
      else if (vtx_field_type[i] == PDM_DOUBLE) {
        s_data = sizeof(double);
      }

      PDM_block_to_part_exch(btp,
                             s_data,
                             PDM_STRIDE_CST_INTERLACED,
                             &vtx_field_stride[i],
                  (void   *) dvtx_field_value[i],
                             NULL,
                  (void ***) &pvtx_field[i]);


      free(dvtx_field_value[i]);
    }
    PDM_block_to_part_free(btp);






    if (visu) {
      PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                            PDM_WRITER_FMT_BIN,
                                            PDM_WRITER_TOPO_CST,
                                            PDM_WRITER_OFF,
                                            "reader_vtk",
                                            "reader_vtk",
                                            comm,
                                            PDM_IO_KIND_MPI_SIMPLE,
                                            1.,
                                            NULL);

      int id_geom = PDM_writer_geom_create_from_mesh_nodal(wrt,
                                                           "reader_vtk",
                                                           pmn);

      int *id_var_vtx = malloc(sizeof(int) * n_vtx_field);
      for (int i = 0; i < n_vtx_field; i++) {
        id_var_vtx[i] = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              (PDM_writer_var_dim_t) vtx_field_stride[i],
                                              PDM_WRITER_VAR_VERTICES,
                                              vtx_field_name[i]);
      }

      PDM_writer_step_beg(wrt, 0.);

      PDM_writer_geom_write(wrt,
                            id_geom);

      for (int i = 0; i < n_vtx_field; i++) {
        PDM_real_t **val_vtx = malloc(sizeof(PDM_real_t *) * n_part);
        for (int j = 0; j < n_part; j++) {

          val_vtx[j] = malloc(sizeof(PDM_real_t) * pn_vtx[j] * vtx_field_stride[i]);

          if (vtx_field_type[i] == PDM_DOUBLE) {
            double *pvf = (double *) pvtx_field[i][j];
            for (int k = 0; k < pn_vtx[j] * vtx_field_stride[i]; k++) {
              val_vtx[j][k] = pvf[k];
            }
          }
          if (vtx_field_type[i] == PDM_INT) {
            int *pvf = (int *) pvtx_field[i][j];
            for (int k = 0; k < pn_vtx[j] * vtx_field_stride[i]; k++) {
              val_vtx[j][k] = pvf[k];
            }
          }

          PDM_writer_var_set(wrt,
                             id_var_vtx[i],
                             id_geom,
                             j,
                             val_vtx[j]);
        }

        PDM_writer_var_write(wrt,
                             id_var_vtx[i]);
        PDM_writer_var_free(wrt,
                            id_var_vtx[i]);

        for (int j = 0; j < n_part; j++) {
          free(val_vtx[j]);
        }
        free(val_vtx);
      }


      PDM_writer_step_end(wrt);

      PDM_writer_free(wrt);

      free(id_var_vtx);
    }








    // if (visu) {
    //   PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(pmn, 0);

    //   for (int j = 0; j < n_part; j++) {

    //     double *val = malloc(sizeof(double) * pn_vtx[j]);

    //     for (int i = 0; i < n_vtx_field; i++) {
    //       for (int k = 0; k < vtx_field_stride[i]; k++) {
    //         char filename2[999];
    //         sprintf(filename2, "vtx_field_%d_comp_%d_part_%d_rank_%d.vtk", i, k, j, i_rank);

    //         double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, j);

    //         if (vtx_field_type[i] == PDM_DOUBLE) {
    //           double *pvf = (double *) pvtx_field[i][j];
    //           for (int l = 0; l < pn_vtx[j]; l++) {
    //             val[l] = pvf[vtx_field_stride[i]*l + k];
    //           }
    //         }
    //         if (vtx_field_type[i] == PDM_INT) {
    //           int *pvf = (int *) pvtx_field[i][j];
    //           for (int l = 0; l < pn_vtx[j]; l++) {
    //             val[l] = pvf[vtx_field_stride[i]*l + k];
    //           }
    //         }

    //         int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, j);


    //         if (t_elt == PDM_MESH_NODAL_POLY_2D) {
    //           int *connec_idx;
    //           int *connec;
    //           PDM_part_mesh_nodal_section_poly2d_get(pmn,
    //                                                  0,
    //                                                  j,
    //                                                  &connec_idx,
    //                                                  &connec);
    //           PDM_vtk_write_polydata_field(filename2,
    //                                        pn_vtx[j],
    //                                        vtx_coord,
    //                                        pvtx_ln_to_gn[j],
    //                                        n_elt,
    //                                        connec_idx,
    //                                        connec,
    //                                        NULL,
    //                                        NULL, NULL,
    //                                        vtx_field_name[i],
    //                                        (const double *) val);
    //         }
    //         else if (t_elt != PDM_MESH_NODAL_POLY_3D) {
    //           int          *connec;
    //           PDM_g_num_t  *numabs;
    //           int          *parent_num;
    //           PDM_g_num_t  *parent_entity_g_num;
    //           int           order;
    //           const char   *ho_ordering;
    //           PDM_part_mesh_nodal_section_std_ho_get(pmn,
    //                                                  0,
    //                                                  j,
    //                                                  &connec,
    //                                                  &numabs,
    //                                                  &parent_num,
    //                                                  &parent_entity_g_num,
    //                                                  &order,
    //                                                  &ho_ordering);

    //           PDM_vtk_write_std_elements_ho_with_vtx_field(filename2,
    //                                                        order,
    //                                                        pn_vtx[j],
    //                                                        vtx_coord,
    //                                                        NULL,//pvtx_ln_to_gn[j],
    //                                                        t_elt,
    //                                                        n_elt,
    //                                                        connec,
    //                                                        NULL,
    //                                                        0,
    //                                                        NULL,
    //                                                        NULL,
    //                                                        1,
    //                                        (const char **) &vtx_field_name[i],
    //                                      (const double **) &val);
    //         }
    //       }
    //     }

    //     free(val);
    //   }
    // }


    free(pn_vtx);
    free(pvtx_ln_to_gn);

    for (int i = 0; i < n_vtx_field; i++) {
      for (int j = 0; j < n_part; j++) {
        free(pvtx_field[i][j]);
      }
      free(pvtx_field[i]);
      free(vtx_field_name  [i]);
    }
    free(pvtx_field);

    free(vtx_field_name  );
    free(vtx_field_type  );
    free(vtx_field_stride);
    free(dvtx_field_value);
  }



  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
  PDM_part_mesh_nodal_free(pmn);




  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;
}
