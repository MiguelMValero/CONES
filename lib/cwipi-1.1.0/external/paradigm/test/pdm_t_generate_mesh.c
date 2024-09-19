#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_generate_mesh.h"
#include "pdm_part_mesh_nodal.h"

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
     "  -visu  Set meshes vtk output.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   visu     Set meshes vtk output
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           int           *visu)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
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
  // Set default values
  int         visu     = 0;

  _read_args(argc,
             argv,
             &visu);
  // Initialize MPI
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Generate sphere mesh
  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_sphere(comm,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        NULL,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        10,
                                                        10,
                                                        1,
                                                        PDM_SPLIT_DUAL_WITH_HILBERT);

  if (visu) {

    char filename[999];
    sprintf(filename, "sphere_mesh_%2.2d.vtk", i_rank);

    int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                               0);

    double* pvtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                                           0);

    PDM_g_num_t *pvtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn,
                                                                   0);

    int pn_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                       0,
                                                       0);

    int         *pelt_vtx            = NULL;
    PDM_g_num_t *pelt_ln_to_gn       = NULL;
    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    PDM_part_mesh_nodal_section_std_get(pmn,
                                        0,
                                        0,
                                        &pelt_vtx,
                                        &pelt_ln_to_gn,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_KEEP);

    PDM_vtk_write_std_elements(filename,
                               pn_vtx,
                               pvtx_coord,
                               pvtx_ln_to_gn,
                               PDM_MESH_NODAL_TRIA3,
                               pn_elt,
                               pelt_vtx,
                               pelt_ln_to_gn,
                               0,
                               NULL,
                               NULL);

  }

  // free
  PDM_part_mesh_nodal_free(pmn);

  // Generate ball mesh
  pmn = PDM_generate_mesh_ball(comm,
                               PDM_MESH_NODAL_TETRA4,
                               1,
                               NULL,
                               1.,
                               0.,
                               0.,
                               0.,
                               0.,
                               10,
                               10,
                               10,
                               0,
                               0.,
                               1,
                               PDM_SPLIT_DUAL_WITH_HILBERT);

    if (visu) {

      char filename[999];
      sprintf(filename, "ball_mesh_%2.2d.vtk", i_rank);

      int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                                 0);

      double* pvtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                                             0);

      PDM_g_num_t *pvtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn,
                                                                     0);

      int pn_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                         0,
                                                         0);

      int         *pelt_vtx            = NULL;
      PDM_g_num_t *pelt_ln_to_gn       = NULL;
      int         *parent_num          = NULL;
      PDM_g_num_t *parent_entity_g_num = NULL;
      PDM_part_mesh_nodal_section_std_get(pmn,
                                          0,
                                          0,
                                          &pelt_vtx,
                                          &pelt_ln_to_gn,
                                          &parent_num,
                                          &parent_entity_g_num,
                                          PDM_OWNERSHIP_KEEP);

      PDM_vtk_write_std_elements(filename,
                                 pn_vtx,
                                 pvtx_coord,
                                 pvtx_ln_to_gn,
                                 PDM_MESH_NODAL_TETRA4,
                                 pn_elt,
                                 pelt_vtx,
                                 pelt_ln_to_gn,
                                 0,
                                 NULL,
                                 NULL);

  }

  // free
  PDM_part_mesh_nodal_free(pmn);

  // Generate ball mesh simplified
  int       n_vtx = 0;
  int       n_elt = 0;
  double   *coords      = NULL;
  int      *elt_vtx_idx = NULL;
  int      *elt_vtx     = NULL;
  PDM_generate_mesh_ball_simplified(comm,
                                    &n_vtx,
                                    &n_elt,
                                    &coords,
                                    &elt_vtx_idx,
                                    &elt_vtx);

  if (visu) {
    log_trace("n_vtx : %d\n", n_vtx);
    log_trace("n_elt : %d\n", n_elt);
    PDM_log_trace_array_double(coords, 3*n_vtx, "coords : ");
    PDM_log_trace_array_int(elt_vtx_idx, n_elt + 1, "elt_vtx_idx : ");
    PDM_log_trace_array_int(elt_vtx, elt_vtx_idx[n_elt], "elt_vtx : ");

    char filename[999];
    sprintf(filename, "ball_mesh_simplified_%2.2d.vtk", i_rank);

    PDM_vtk_write_std_elements(filename,
                               n_vtx,
                               coords,
                               NULL,
                               PDM_MESH_NODAL_TETRA4,
                               n_elt,
                               elt_vtx,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }

  // free
  free(coords     );
  free(elt_vtx_idx);
  free(elt_vtx    );

  // Generate rectangle mesh
  pmn = PDM_generate_mesh_rectangle(comm,
                                    PDM_MESH_NODAL_POLY_2D,
                                    1,
                                    NULL,
                                    0.,
                                    0.,
                                    0.,
                                    10.,
                                    5.,
                                    10,
                                    10,
                                    1,
                                    PDM_SPLIT_DUAL_WITH_HILBERT);

  if (visu) {

    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "rectangle_mesh");
  }

  // free
  PDM_part_mesh_nodal_free(pmn);

  // Generate simplified rectangle mesh
  n_vtx = 0;
  n_elt = 0;
  coords      = NULL;
  elt_vtx_idx = NULL;
  elt_vtx     = NULL;
  PDM_generate_mesh_rectangle_simplified(comm,
                                         10,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx);

  if (visu) {
    log_trace("n_vtx : %d\n", n_vtx);
    log_trace("n_elt : %d\n", n_elt);
    PDM_log_trace_array_double(coords, 3*n_vtx, "coords : ");
    PDM_log_trace_array_int(elt_vtx_idx, n_elt + 1, "elt_vtx_idx : ");
    PDM_log_trace_array_int(elt_vtx, elt_vtx_idx[n_elt], "elt_vtx : ");
  }

  // free
  free(coords     );
  free(elt_vtx_idx);
  free(elt_vtx    );

  // Generate parallelepiped mesh
  pmn = PDM_generate_mesh_parallelepiped(comm,
                                         PDM_MESH_NODAL_PYRAMID5,
                                         1,
                                         NULL,
                                         0.,
                                         0.,
                                         0.,
                                         10.,
                                         1.8,
                                         99.42,
                                         10,
                                         20,
                                         5,
                                         1,
                                         PDM_SPLIT_DUAL_WITH_HILBERT);

  if (visu) {

    char filename[999];
    sprintf(filename, "parallelepiped_mesh_%2.2d.vtk", i_rank);

    int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                               0);

    double* pvtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                                           0);

    PDM_g_num_t *pvtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn,
                                                                   0);

    int pn_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                       0,
                                                       0);

    int         *pelt_vtx            = NULL;
    PDM_g_num_t *pelt_ln_to_gn       = NULL;
    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    PDM_part_mesh_nodal_section_std_get(pmn,
                                        0,
                                        0,
                                        &pelt_vtx,
                                        &pelt_ln_to_gn,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_KEEP);

    PDM_vtk_write_std_elements(filename,
                               pn_vtx,
                               pvtx_coord,
                               pvtx_ln_to_gn,
                               PDM_MESH_NODAL_PYRAMID5,
                               pn_elt,
                               pelt_vtx,
                               pelt_ln_to_gn,
                               0,
                               NULL,
                               NULL);

  }


  // Generate simplified parallelepiped mesh
  n_vtx = 0;
  n_elt = 0;
  coords      = NULL;
  elt_vtx_idx = NULL;
  elt_vtx     = NULL;
  PDM_generate_mesh_parallelepiped_simplified(comm,
                                              10,
                                              &n_vtx,
                                              &n_elt,
                                              &coords,
                                              &elt_vtx_idx,
                                              &elt_vtx);

  if (visu) {
    log_trace("n_vtx : %d\n", n_vtx);
    log_trace("n_elt : %d\n", n_elt);
    PDM_log_trace_array_double(coords, 3*n_vtx, "coords : ");
    PDM_log_trace_array_int(elt_vtx_idx, n_elt + 1, "elt_vtx_idx : ");
    PDM_log_trace_array_int(elt_vtx, elt_vtx_idx[n_elt], "elt_vtx : ");

    char filename[999];
    sprintf(filename, "parallelepiped_mesh_simplified_%2.2d.vtk", i_rank);

    PDM_vtk_write_std_elements(filename,
                               n_vtx,
                               coords,
                               NULL,
                               PDM_MESH_NODAL_TETRA4,
                               n_elt,
                               elt_vtx,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }

  // free
  free(coords     );
  free(elt_vtx_idx);
  free(elt_vtx    );

  // free
  PDM_part_mesh_nodal_free(pmn);

  PDM_MPI_Finalize();

  return 0;
}
