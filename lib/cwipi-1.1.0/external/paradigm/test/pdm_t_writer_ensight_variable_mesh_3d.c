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
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_generate_mesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);


  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_ASCII,
                                        PDM_WRITER_TOPO_VARIABLE,
                                        PDM_WRITER_OFF,
                                        "writer_ensight_variable_mesh_3d",
                                        "writer_ensight_variable_mesh_3d",
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       "writer_ensight_variable_mesh_3d",
                                       1);


  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;

  PDM_g_num_t *vtx_ln_to_gn = NULL;
  PDM_g_num_t *elt_ln_to_gn = NULL;

  double time = 0.;

  for (int it = 0; it < 5; it++) {

    if (i_rank == 0) {
      printf("it %d\n", it);
    }

    if (it > 0) {
      free(coords);
      coords = NULL;
      free(elt_vtx_idx);
      elt_vtx_idx = NULL;
      free(elt_vtx);
      elt_vtx = NULL;
      free(vtx_ln_to_gn);
      vtx_ln_to_gn = NULL;
      free(elt_ln_to_gn);
      elt_ln_to_gn = NULL;
    }

    PDM_generate_mesh_rectangle_simplified(comm,
                                           it+2,
                                           &n_vtx,
                                           &n_elt,
                                           &coords,
                                           &elt_vtx_idx,
                                           &elt_vtx);

    for (int i = 0; i < n_vtx; i++) {
      coords[3*i+2] = cos(4*coords[3*i]);
    }

    PDM_gen_gnum_t *gnum_vtx = PDM_gnum_create(3,
                                               1,
                                               PDM_FALSE,
                                               1.,
                                               comm,
                                               PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_coords(gnum_vtx,
                             0,
                             n_vtx,
                             coords,
                             NULL);

    PDM_gnum_compute(gnum_vtx);

    vtx_ln_to_gn = PDM_gnum_get(gnum_vtx, 0);

    PDM_gnum_free(gnum_vtx);


    double *elt_center = malloc(sizeof(double) * n_elt * 3);
    for (int i = 0; i < n_elt; i++) {
      for (int j = 0; j < 3; j++) {
        elt_center[3*i+j] = 0;
      }

      for (int k = elt_vtx_idx[i]; k < elt_vtx_idx[i+1]; k++) {
        int vtx_id = elt_vtx[k] - 1;
        for (int j = 0; j < 3; j++) {
          elt_center[3*i+j] += coords[3*vtx_id+j];
        }
      }

      for (int j = 0; j < 3; j++) {
        elt_center[3*i+j] /= elt_vtx_idx[i+1] - elt_vtx_idx[i];
      }
    }

    PDM_gen_gnum_t *gnum_elt = PDM_gnum_create(3,
                                               1,
                                               PDM_FALSE,
                                               1.,
                                               comm,
                                               PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_coords(gnum_elt,
                             0,
                             n_elt,
                             elt_center,
                             NULL);

    PDM_gnum_compute(gnum_elt);

    elt_ln_to_gn = PDM_gnum_get(gnum_elt, 0);

    PDM_gnum_free(gnum_elt);
    free(elt_center);


    PDM_writer_step_beg(wrt, time);

    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              0,
                              n_vtx,
                              coords,
                              vtx_ln_to_gn,
                              PDM_OWNERSHIP_USER);

    int id_block = PDM_writer_geom_bloc_add(wrt,
                                            id_geom,
                                            PDM_WRITER_TRIA3,
                                            PDM_OWNERSHIP_USER);

    PDM_writer_geom_bloc_std_set(wrt,
                                 id_geom,
                                 id_block,
                                 0,
                                 n_elt,
                                 elt_vtx,
                                 elt_ln_to_gn);


    PDM_writer_geom_write(wrt,
                          id_geom);

    PDM_writer_step_end(wrt);

    PDM_writer_geom_data_reset(wrt,
                               id_geom);

    time += 1.;

  }

  PDM_writer_geom_free(wrt,
                       id_geom);

  PDM_writer_free(wrt);

  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(vtx_ln_to_gn);
  free(elt_ln_to_gn);


  PDM_MPI_Finalize();

  return 0;
}
