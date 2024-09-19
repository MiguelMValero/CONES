#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_ho_basis.h"
#include "pdm_ho_ordering.h"
#include "pdm_ho_bezier.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"

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
_read_args
(
 int      argc,
 char   **argv,
 int     *order,
 double  *noise,
 int     *visu
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-noise") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *noise = atof(argv[i]);
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
  int    order = 3;
  double noise = 0.;
  int    visu  = 0;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &order,
             &noise,
             &visu);


  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);

  int n_node = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  if (order > 3) {
    int *ijk = NULL;
    ijk = PDM_vtk_lagrange_to_ijk(PDM_MESH_NODAL_BARHO, order);
    PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                     PDM_MESH_NODAL_BARHO,
                                     order,
                                     PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order),
                                     ijk);
    free (ijk);
  }


  /*
   *  Build a Bezier P^order curve
   */
  int *elt_node = malloc(sizeof(int) * n_node);
  for (int i = 0; i < n_node; i++) {
    elt_node[i] = i + 1;
  }

  double step = 1. / (double) order;
  double *node_coord = malloc(sizeof(double) * n_node * 3);
  for (int i = 0; i <= order; i++) {
    node_coord[3*i  ] = i*step + noise * (2*(double) rand() / (double) RAND_MAX - 1);
    node_coord[3*i+1] =          noise * (2*(double) rand() / (double) RAND_MAX - 1);
    node_coord[3*i+2] =          noise * (2*(double) rand() / (double) RAND_MAX - 1);
  }



  /*
   *  Evaluate and subdivide
   */
  double t = (double) rand() / (double) RAND_MAX;
  double p[3];
  double *sub_node_coord[2] = {NULL};
  for (int i = 0; i < 2; i++) {
    sub_node_coord[i] = malloc(sizeof(double) * n_node * 3);
  }

  PDM_ho_bezier_de_casteljau_curve(3,
                                   order,
                                   t,
                                   node_coord,
                                   p,
                                   sub_node_coord[0],
                                   sub_node_coord[1]);


  // check weights
  if (order <= 3) {
    double *weight = malloc(sizeof(double) * n_node);
    PDM_ho_bezier_basis(PDM_MESH_NODAL_BARHO_BEZIER,
                        order,
                        1,
                        &t,
                        weight);


    double *weight2 = malloc(sizeof(double) * n_node);
    PDM_ho_bezier_de_casteljau_curve(n_node,
                                     order,
                                     t,
                                     NULL,
                                     weight2,
                                     NULL,
                                     NULL);

    double max_diff = 0;
    for (int i = 0; i < n_node; i++) {
      double diff = PDM_ABS(weight[i] - weight2[i]);
      max_diff = PDM_MAX(max_diff, diff);
    }

    printf("max_diff = %e\n", max_diff);
    free(weight);
    free(weight2);
  }


  /*
   *  Derivative
   */
  double *deriv_node_coord = malloc(sizeof(double) * order * 3);
  PDM_ho_bezier_curve_derivative(3,
                                 order,
                                 node_coord,
                                 deriv_node_coord);

  double dp_dt[3];
  PDM_ho_bezier_de_casteljau_curve(3,
                                   order-1,
                                   t,
                                   deriv_node_coord,
                                   dp_dt,
                                   NULL,
                                   NULL);



  /*
   *  Visu VTK
   */
  PDM_Mesh_nodal_reorder_elt_vtx(PDM_MESH_NODAL_BARHO_BEZIER,
                                 order,
                                 NULL,
                                 "PDM_HO_ORDERING_VTK",
                                 1,
                                 elt_node,
                                 elt_node);

  if (visu) {
    PDM_vtk_write_std_elements_ho("bezier_curve.vtk",
                                  order,
                                  n_node,
                                  node_coord,
                                  NULL,
                                  PDM_MESH_NODAL_BARHO_BEZIER,
                                  1,
                                  elt_node,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);

    char filename[999];
    for (int i = 0; i < 2; i++) {
      sprintf(filename, "bezier_subcurve_%d.vtk", i);
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    n_node,
                                    sub_node_coord[i],
                                    NULL,
                                    PDM_MESH_NODAL_BARHO_BEZIER,
                                    1,
                                    elt_node,
                                    NULL,
                                    0,
                                    NULL,
                                    NULL);
    }


    FILE *f = fopen("bezier_deriv_curve.vtk", "w");
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "bezier_deriv\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS 1 double\n");
    fprintf(f, "%f %f %f\n", p[0], p[1], p[2]);
    fprintf(f, "CELLS 1 2\n1 0\n");
    fprintf(f, "CELL_TYPES 1\n1\n");
    fprintf(f, "POINT_DATA 1\n");
    fprintf(f, "VECTORS vector double\n");
    fprintf(f, "%f %f %f\n", dp_dt[0], dp_dt[1], dp_dt[2]);
    fclose(f);
  }

  /*
   *  Free memory
   */
  free(node_coord);
  for (int i = 0; i < 2; i++) {
    free(sub_node_coord[i]);
  }
  free(deriv_node_coord);
  free(elt_node);

  PDM_MPI_Finalize();

  return 0;
}



