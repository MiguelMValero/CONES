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
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_triangle.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage (EXIT_FAILURE);
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
  PDM_MPI_Init (&argc, &argv);

  int           i_rank;
  int           numProcs;
  int           visu = 0;

  _read_args (argc,
              argv,
              &visu);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  int order = 3;
  int n_vtx_triangle;

  // Set-up test
  double tria_coord[330];

  // order 1
  if (order == 1) {
    n_vtx_triangle = 3;
    // (i,j,k)=(0,0,1)
    tria_coord[0] = 0;
    tria_coord[1] = 0;
    tria_coord[2] = 10;

    // (i,j,k)=(1,0,0)
    tria_coord[0] = 10;
    tria_coord[1] = 0;
    tria_coord[2] = 10;

    // (i,j,k)=(0,1,0)
    tria_coord[0] = 0;
    tria_coord[1] = 10;
    tria_coord[2] = 10;
  }

  // order 2
  if (order == 2) {
    n_vtx_triangle = 6;
    // (i,j,k)=(0,0,2)
    tria_coord[0] = 0;
    tria_coord[1] = 0;
    tria_coord[2] = 10;

    // (i,j,k)=(1,0,1)
    tria_coord[3] = 10;
    tria_coord[4] = 0;
    tria_coord[5] = 20;

    // (i,j,k)=(2,0,0)
    tria_coord[6] = 20;
    tria_coord[7] = 0;
    tria_coord[8] = 10;

    // (i,j,k)=(0,1,1)
    tria_coord[9] = 0;
    tria_coord[10] = 10;
    tria_coord[11] = 20;

    // (i,j,k)=(1,1,0)
    tria_coord[12] = 10;
    tria_coord[13] = 10;
    tria_coord[14] = 20;

    // (i,j,k)=(0,2,0)
    tria_coord[15] = 0;
    tria_coord[16] = 20;
    tria_coord[17] = 10;
  }

  // order 3
  if (order == 3) {
    n_vtx_triangle = 10;
    // (i,j,k)=(0,0,3)
    tria_coord[0] = 0;
    tria_coord[1] = 0;
    tria_coord[2] = 10;

    // (i,j,k)=(1,0,2)
    tria_coord[3] = 10;
    tria_coord[4] = 0;
    tria_coord[5] = 20;

    // (i,j,k)=(2,0,1)
    tria_coord[6] = 20;
    tria_coord[7] = 0;
    tria_coord[8] = 30;

    // (i,j,k)=(3,0,0)
    tria_coord[9] = 30;
    tria_coord[10] = 0;
    tria_coord[11] = 10;

    // (i,j,k)=(0,1,2)
    tria_coord[12] = 0;
    tria_coord[13] = 10;
    tria_coord[14] = 20;

    // (i,j,k)=(1,1,1)
    tria_coord[15] = 10;
    tria_coord[16] = 10;
    tria_coord[17] = 40;

    // (i,j,k)=(2,1,0)
    tria_coord[18] = 20;
    tria_coord[19] = 20;
    tria_coord[20] = 30;

    // (i,j,k)=(0,2,1)
    tria_coord[21] = 0;
    tria_coord[22] = 30;
    tria_coord[23] = 30;

    // (i,j,k)=(1,2,0)
    tria_coord[24] = 10;
    tria_coord[25] = 20;
    tria_coord[26] = 30;

    // (i,j,k)=(0,3,0)
    tria_coord[27] = 0;
    tria_coord[28] = 30;
    tria_coord[29] = 10;

  }

  double vector_du[330];
  double vector_dv[330];
  double normal[330];

  for (int i = 0; i < n_vtx_triangle; i++) {
    for (int j = 0; j < 3; j++) {
      vector_du[3*i + j] = 0;
      vector_dv[3*i + j] = 0;
      normal[3*i + j] = 0;
    }
  }

  // Get xyz coordinates from uvw coordinates

  double uvw[3];
  double weights[n_vtx_triangle];
  double dw_du[n_vtx_triangle];
  double dw_dv[n_vtx_triangle];
  int n_vtx = n_vtx_triangle;
  double pt_du[3];
  double pt_dv[3];
  double pt_normal[3];
  double pt_normal_module;


  for (int i = 0; i < 11; i++) {

    // double rd1 = (double)rand() / (double)RAND_MAX;
    // uvw[0] = rd1;

    // double rd2 = (double)rand() / (double)RAND_MAX;
    // rd2 = rd2*(1-rd1); // normalize between [0,1-rd1]
    // uvw[1] = rd2;

    // uvw[2] = 1 - uvw[0] - uvw[1];

    uvw[0] = 0.1 * i;

    for (int k = 0; k < (10-i+1); k++) {

      uvw[1] = 0.1 * k;
      uvw[2] = PDM_ABS(1 - uvw[0] - uvw[1]);

      PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIAHO,
                          order,
                          1,
                          uvw,
                          weights);

      double xyz_Pn[3] = {0, 0, 0};

      for (int j = 0; j < 3; j++) {
        for (int l = 0; l < n_vtx_triangle; l++) {
          xyz_Pn[j] += weights[l] * tria_coord[3*l+j];
        } // end loop on triangle points
      }

      tria_coord[3*n_vtx    ] = xyz_Pn[0];
      tria_coord[3*n_vtx + 1] = xyz_Pn[1];
      tria_coord[3*n_vtx + 2] = xyz_Pn[2];

      PDM_ho_bezier_basis_derivative(PDM_MESH_NODAL_TRIAHO,
                                     order,
                                     1,
                                     uvw,
                                     dw_du,
                                     dw_dv,
                                     NULL);

      for (int j = 0; j < 3; j++) {
        pt_du[j] = 0;
        pt_dv[j] = 0;
        for (int l = 0; l < n_vtx_triangle; l++) {
          pt_du[j] += dw_du[l] * tria_coord[3*l+j];
          pt_dv[j] += dw_dv[l] * tria_coord[3*l+j];
        } // end loop on triangle points
        vector_du[3*n_vtx + j] = pt_du[j];
        vector_dv[3*n_vtx + j] = pt_dv[j];
      }
      PDM_CROSS_PRODUCT(pt_normal, pt_du, pt_dv);

      // PDM_log_trace_array_double(pt_normal, 3, "pt_normal:");

      pt_normal_module = PDM_MODULE(pt_normal);
      normal[3*n_vtx    ] = pt_normal[0] / pt_normal_module;
      normal[3*n_vtx + 1] = pt_normal[1] / pt_normal_module;
      normal[3*n_vtx + 2] = pt_normal[2] / pt_normal_module;

      // log_trace("n0 %lf n1 %lf n2 %lf\n", pt_normal[0] / pt_normal_module, pt_normal[1] / pt_normal_module, pt_normal[2] / pt_normal_module);

      n_vtx++;

    } // end loop on i
  } // end loop on k

  // vtk ouput of ho element
  char filename1[999];
  sprintf(filename1, "P%d_triangle_ho.vtk", order);

  PDM_g_num_t *vtx_g_num = malloc(sizeof(PDM_g_num_t) * n_vtx);
  int *face_vtx = malloc(sizeof(int) * n_vtx_triangle);
  PDM_g_num_t *face_g_num = malloc(sizeof(PDM_g_num_t) * 1);
  const char  *vtx_field_name = "ho_bezier_basis";
  double      *vtx_field = malloc(sizeof(double) * n_vtx);

  for (int j = 0; j < n_vtx_triangle; j++) {
    vtx_g_num[j] = j + 1;
    vtx_field[j] = 0;
  }

  if (order == 1) {

    face_vtx[0] = 1;
    face_vtx[1] = 2;
    face_vtx[2] = 3;

  }

  if (order == 2) {

    face_vtx[0] = 1;
    face_vtx[1] = 3;
    face_vtx[2] = 6;
    face_vtx[3] = 2;
    face_vtx[4] = 5;
    face_vtx[5] = 4;

  }

  if (order == 3) {

    face_vtx[0] = 1;
    face_vtx[1] = 4;
    face_vtx[2] = 10;
    face_vtx[3] = 2;
    face_vtx[4] = 3;
    face_vtx[5] = 7;
    face_vtx[6] = 9;
    face_vtx[7] = 8;
    face_vtx[8] = 5;
    face_vtx[9] = 6;

  }

  face_g_num[0] = 1;

  if (visu) {
    PDM_vtk_write_std_elements_ho(filename1,
                                  order,
                                  n_vtx_triangle,
                                  tria_coord,
                                  vtx_g_num,
                                  PDM_MESH_NODAL_TRIAHO_BEZIER,
                                  1,
                                  face_vtx,
                                  face_g_num,
                                  0,
                                  NULL,
                                  NULL);
  }

  // vtk output of points
  if (visu) {
    char filename2[999];
    sprintf(filename2, "P%d_triangle_vtx.vtk", order);

    for (int j = n_vtx_triangle; j < n_vtx; j++) {
      vtx_g_num[j] = j + 1;
      vtx_field[j] = 1;
    }

    double *vector_derivatives[2] = {vector_du, vector_dv};

    const char* vector_derivatives_names[] = {"dw_du", "dw_dv", 0};

    double *vector_normal[1] = {normal};

    const char* normal_name[] = {"n", 0};

    PDM_vtk_write_point_cloud_with_field(filename2,
                                         n_vtx,
                                         tria_coord,
                                         vtx_g_num,
                                         NULL,
                                         1,
                         (const char **) &vtx_field_name,
                       (const double **) &vtx_field,
                                         2,
                         (const char **) &vector_derivatives_names,
                       (const double **) &vector_derivatives,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &vector_normal);

  }
  free(face_vtx);
  free(face_g_num);
  free(vtx_g_num);
  free(vtx_field);

  PDM_MPI_Finalize ();

}
