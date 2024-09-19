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
           PDM_g_num_t   *n_layer,
           double        *x_center,
           double        *y_center,
           double        *z_center,
           double        *radius,
           double        *geometric_ratio,
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
    else if (strcmp(argv[i], "-n_layer") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *n_layer = (PDM_g_num_t) _n;
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
    else if (strcmp(argv[i], "-ratio") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *geometric_ratio = atof(argv[i]);
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
  PDM_g_num_t n               = 4;
  PDM_g_num_t n_layer         = 3;
  double      x_center        = 0;
  double      y_center        = 0;
  double      z_center        = 0;
  double      radius_interior = 1;
  double      radius_exterior = 2;
  double      geometric_ratio = 1.;
  int         visu            = 0;

  _read_args(argc,
             argv,
             &n,
             &n_layer,
             &x_center,
             &y_center,
             &z_center,
             &radius_interior,
             &geometric_ratio,
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
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_vol_hollow_gen_nodal(comm,
                                  n,
                                  n_layer,
                                  x_center,
                                  y_center,
                                  z_center,
                                  radius_interior,
                                  radius_exterior,
                                  geometric_ratio,
                                  &dmn);

  if (visu) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "hollow_volume_");

    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "hollow_surface_");
  }
  PDM_DMesh_nodal_free(dmn);

  PDM_MPI_Finalize();

  return 0;
}

