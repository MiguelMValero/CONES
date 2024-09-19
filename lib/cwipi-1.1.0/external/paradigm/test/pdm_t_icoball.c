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
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_geom_elem.h"
#include "pdm_dmesh_nodal.h"

#include "pdm_sphere_vol_gen.h"


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
           PDM_g_num_t   *n,
           double        *x_center,
           double        *y_center,
           double        *z_center,
           double        *radius,
           int           *visu)
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
        *n = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-cx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *x_center = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-cy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *y_center = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-cz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *z_center = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-r") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *radius = atof(argv[i]);
    }
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
  /*
   *  Set default values
   */
  PDM_g_num_t n        = 5;
  double      x_center = 0;
  double      y_center = 0;
  double      z_center = 0;
  double      radius   = 1;
  int         visu     = 0;

  _read_args(argc,
             argv,
             &n,
             &x_center,
             &y_center,
             &z_center,
             &radius,
             &visu);



  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Generate distributed Icoball
   */
  int          dn_vtx        = 0;
  double      *dvtx_coord    = NULL;
  int          dn_cell       = 0;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *dcell_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  PDM_g_num_t *distrib_face  = NULL;
  PDM_g_num_t *distrib_cell  = NULL;


  PDM_sphere_vol_icosphere_gen(comm,
                               n,
                               x_center,
                               y_center,
                               z_center,
                               radius,
                               &dvtx_coord,
                               &dface_vtx,
                               &dcell_vtx,
                               &distrib_vtx,
                               &distrib_face,
                               &distrib_cell);
  dn_vtx  = (int) (distrib_vtx[i_rank+1]  - distrib_vtx[i_rank]);


  // PDM_log_trace_array_long(distrib_vtx,
  //                          n_rank+1,
  //                          "distrib_vtx : ");

  /*
   *  Visu VTK
   */
  PDM_g_num_t *dvtx_ln_to_gn = malloc(sizeof(PDM_g_num_t) * dn_vtx);
  for (int i = 0; i < dn_vtx; i++) {
    dvtx_ln_to_gn[i] = distrib_vtx[i_rank] + i + 1;
  }


  char filename[999];
  // sprintf(filename, "icoball_dvtx_coord_%2.2d.vtk", i_rank);
  // PDM_vtk_write_point_cloud(filename,
  //                           dn_vtx,
  //                           dvtx_coord,
  //                           dvtx_ln_to_gn,
  //                           NULL);


  dn_cell = (int) (distrib_cell[i_rank+1] - distrib_cell[i_rank]);

  int *dcell_vtx_idx = PDM_array_new_idx_from_const_stride_int(4,
                                                               dn_cell);

  int          pn_vtx        = 0;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;
  int         *pcell_vtx_idx = NULL;
  int         *pcell_vtx     = NULL;
  double      *pvtx_coord    = NULL;

  PDM_g_num_t *pcell_ln_to_gn = malloc(sizeof(PDM_g_num_t) * dn_cell);
  for (int i = 0; i < dn_cell; i++) {
    pcell_ln_to_gn[i] = distrib_cell[i_rank] + i + 1;
  }

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_cell,
                                                           dcell_vtx_idx,
                                                           dcell_vtx,
                                                           dn_cell,
                                                           pcell_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pcell_vtx_idx,
                                                           &pcell_vtx);
  /* Coordinates */
  double **tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);



  if (visu) {
    sprintf(filename, "icoball_dcell_%2.2d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               pn_vtx,
                               pvtx_coord,
                               pvtx_ln_to_gn,
                               PDM_MESH_NODAL_TETRA4,
                               dn_cell,
                               pcell_vtx,
                               pcell_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }
  double *volume = malloc(sizeof(double) * dn_cell);
  // PDM_geom_elem_tetra_oriented_volume(dn_cell,
  //                                     pcell_vtx,
  //                                     pvtx_coord,
  //                                     volume,
  //                                     NULL,
  //                                     NULL);
  int count = 0;
  for (int i = 0; i < dn_cell; i++) {

    int *tv = pcell_vtx + 4*i;

    double u[3], v[3], w[3];
    for (int j = 0; j < 3; j++) {
      u[j] = pvtx_coord[3*(tv[1]-1)+j] - pvtx_coord[3*(tv[0]-1)+j];
      v[j] = pvtx_coord[3*(tv[2]-1)+j] - pvtx_coord[3*(tv[0]-1)+j];
      w[j] = pvtx_coord[3*(tv[3]-1)+j] - pvtx_coord[3*(tv[0]-1)+j];
    }

    double uv[3];
    PDM_CROSS_PRODUCT(uv, u, v);
    volume[i] = PDM_DOT_PRODUCT(uv, w) / 6.;

    if (volume[i] < 0) {
      count++;
    }
  }
  // log_trace("%d cells with negative volume / %d\n", count, dn_cell);
  free(volume);

  free(pvtx_ln_to_gn);
  free(pcell_vtx_idx);
  free(pcell_vtx);
  free(pvtx_coord);
  free(dcell_vtx_idx);


  free(dvtx_ln_to_gn);
  free(distrib_vtx);
  free(distrib_face);
  free(distrib_cell);


  free(dvtx_coord);
  free(dface_vtx);
  free(dcell_vtx);



  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_vol_icosphere_gen_nodal(comm,
                                     n,
                                     x_center,
                                     y_center,
                                     z_center,
                                     radius,
                                     &dmn);

  if (visu) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "icoball_volume_");

    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "icoball_surface_");
  }
  PDM_DMesh_nodal_free(dmn);
  free(pcell_ln_to_gn);

  PDM_MPI_Finalize();

  return 0;
}

