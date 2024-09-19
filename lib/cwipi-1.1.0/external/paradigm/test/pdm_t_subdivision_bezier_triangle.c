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
#include "pdm_array.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

#define ij2idx(i, j, n) ((i) + (j)*((n) + 1) - ((j)-1)*(j)/2)




static void
_subdivide_bezier_triangle
(
 const int  dim,
 const int  order,
 double    *b,
 double    *b0,
 double    *b1,
 double    *b2,
 double    *b3
 )
{
  const int n = (order+1)*(order+2)/2;

  double atr[n*dim], ars[n*dim];
  PDM_ho_bezier_de_casteljau_triangle(dim, order, 0.0, 0.5, b,   NULL,
                                      atr, ars, NULL);

  // double bat[n*dim];
  double bra[n*dim];
  PDM_ho_bezier_de_casteljau_triangle(dim, order, 0.5, 0.5, atr, NULL,
                                      b0, NULL, bra);
                                      //bat, NULL, bra);

  //double csa[n*dim];
  PDM_ho_bezier_de_casteljau_triangle(dim, order, 0.5, 0.5, ars, NULL,
                                      NULL, NULL, b2);
                                      //NULL, NULL, csa);

  // double cbr[n*dim];
  // double cab[n*dim];
  PDM_ho_bezier_de_casteljau_triangle(dim, order, 1.0, 1.0, bra, NULL,
                                      b1, NULL, b3);
                                      //cbr, NULL, cab);

  // TODO: check/transform the local uv frame of each subtriangle...
}



static void
_subdivide_bezier_triangle2
(
 const int     dim,
 const int     order,
 const double  umin,
 const double  vmin,
 const double  wmin,
 double       *b,
 double       *bsub
 )
{
  const int n = (order+1)*(order+2)/2;

  double ctr[n*dim];
  double uc = umin;
  double vc = 1 - umin - wmin;
  PDM_ho_bezier_de_casteljau_triangle(dim, order, uc, vc, b, NULL,
                                      ctr, NULL, NULL);

  double ub = 1 - vmin - wmin;
  double vb = vmin;
  double u1b = 1 - ub + (uc - 1) * vb/vc; // /!\ if vc == 0
  double v1b = ub - uc*vb/vc;             // /!\ if vc == 0
  double bct[n*dim];
  PDM_ho_bezier_de_casteljau_triangle(dim, order, u1b, v1b, ctr, NULL,
                                      bct, NULL, NULL);

  double ua = umin;
  double va = vmin;
  double B = -vb/(vc * v1b);
  double C = (1 - B*uc)/vc;
  double E = (vb/vc - 1) / v1b;
  double F = -(1 + E*uc)/vc;
  double u2a = B*ua + C*va;
  double v2a = 1 + E*ua + F*va;
  PDM_ho_bezier_de_casteljau_triangle(dim, order, u2a, v2a, bct, NULL,
                                      bsub, NULL, NULL);
}



static void
_tri_mesh
(
 const int   n,
 int        *n_node,
 double    **node_uv,
 int        *n_tria,
 int       **tria_node
 )
{
  *n_node = (n+1) * (n+2) / 2;
  *n_tria = n*n;

  *node_uv   = malloc(sizeof(double) * (*n_node) * 2);
  *tria_node = malloc(sizeof(int   ) * (*n_tria) * 3);

  double step = 1. / (double) n;

  int idx = 0;
  for (int j = 0; j <= n; j++) {
    for (int i = 0; i <= n - j; i++) {
      (*node_uv)[idx++] = i * step;
      (*node_uv)[idx++] = j * step;
    }
  }

  idx = 0;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n - j; i++) {
      (*tria_node)[idx++] = 1 + ij2idx(i,   j,   n);
      (*tria_node)[idx++] = 1 + ij2idx(i+1, j,   n);
      (*tria_node)[idx++] = 1 + ij2idx(i,   j+1, n);

      if (i < n-j-1) {
        (*tria_node)[idx++] = 1 + ij2idx(i+1, j,   n);
        (*tria_node)[idx++] = 1 + ij2idx(i+1, j+1, n);
        (*tria_node)[idx++] = 1 + ij2idx(i,   j+1, n);
      }
    }
  }
}




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


  if (order > 3) {
    int *ijk = NULL;
    ijk = PDM_vtk_lagrange_to_ijk(PDM_MESH_NODAL_TRIAHO, order);
    PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                     PDM_MESH_NODAL_TRIAHO,
                                     order,
                                     PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order),
                                     ijk);
    free (ijk);
  }


  /*
   *  Build a Bezier P^order triangle
   */
  int     n_node    = 0;
  int     n_tria    = 0;
  double *node_uv   = NULL;
  int    *tria_node = NULL;
  _tri_mesh(order,
            &n_node,
            &node_uv,
            &n_tria,
            &tria_node);

  double ctr[3] = {0, 0, -1};
  double rad = 2;
  double *node_xyz = malloc(sizeof(double) * n_node * 3);
  for (int i = 0; i < n_node; i++) {
    double u = node_uv[2*i  ];
    double v = node_uv[2*i+1];

    double x = u - 0.5*(1 - v);
    double y = 0.5*sqrt(3)*v;
    double z = 0;

    double vec[3] = {x - ctr[0], y - ctr[1], z - ctr[2]};
    double r = PDM_MODULE(vec);

    vec[2] *= 4.;

    for (int j = 0; j < 3; j++) {
      node_xyz[3*i+j] = ctr[j] + vec[j] * rad/r;
      node_xyz[3*i+j] += noise * (double) rand() / (double) RAND_MAX;
    }

  }


  if (0) {
    log_trace("node_uv =\n");
    for (int i = 0; i < n_node; i++) {
      log_trace("%f %f\n", node_uv[2*i], node_uv[2*i+1]);
    }
    log_trace("node_xyz =\n");
    for (int i = 0; i < n_node; i++) {
      log_trace("%f %f %f\n", node_xyz[3*i], node_xyz[3*i+1], node_xyz[3*i+2]);
    }
  }


  /*
   *  Derivatives as P^(order-1) Bézier triangles
   */
  double *deriv[2];
  for (int i = 0; i < 2; i++) {
    deriv[i] = malloc(sizeof(double) * 3 * order*(order+1)/2);
  }
  PDM_ho_bezier_triangle_derivatives(3,
                                     order,
                                     node_xyz,
                                     deriv[0],
                                     deriv[1]);


  /*
   *  de Casteljau algorithm (evaluate xyz and build 3 subtriangles)
   */
  double uv[2] = {0.17, 0.12};

  double p[3];
  double *sub3_node_xyz[3];
  for (int i = 0; i < 3; i++) {
    sub3_node_xyz[i] = malloc(sizeof(double) * n_node * 3);
  }
  PDM_ho_bezier_de_casteljau_triangle(3,
                                      order,
                                      uv[0],
                                      uv[1],
                                      node_xyz,
                                      p,
                                      sub3_node_xyz[0],//NULL,//
                                      sub3_node_xyz[1],//NULL,//
                                      sub3_node_xyz[2]);//NULL);//

  // log_trace("u = %f, v = %f\n", uv[0], uv[1]);
  // log_trace("p = %f %f %f\n", p[0], p[1], p[2]);


  // evaluate tangent and normal vectors
  double dpdu[3], dpdv[3], normal[3];
  PDM_ho_bezier_de_casteljau_triangle(3, order-1, uv[0], uv[1], deriv[0], dpdu, NULL, NULL, NULL);
  PDM_ho_bezier_de_casteljau_triangle(3, order-1, uv[0], uv[1], deriv[1], dpdv, NULL, NULL, NULL);
  PDM_CROSS_PRODUCT(normal, dpdu, dpdv);

  for (int i = 0; i < 2; i++) {
    free(deriv[i]);
  }

  // check weights
  if (order <= 3) {
    double *weight = malloc(sizeof(double) * n_node);
    PDM_ho_bezier_basis(PDM_MESH_NODAL_TRIAHO_BEZIER,
                        order,
                        1,
                        uv,
                        weight);

    double q[3] = {-p[0], -p[1], -p[2]};
    for (int i = 0; i < n_node; i++) {
      for (int j = 0; j < 3; j++) {
        q[j] += weight[i] * node_xyz[3*i+j];
      }
    }
    // log_trace("diff p = %f %f %f\n", q[0], q[1], q[2]);

    double *weight2 = malloc(sizeof(double) * n_node);
    // de Casteljau to compute weights
    PDM_ho_bezier_de_casteljau_triangle(n_node, order, uv[0], uv[1],
                                        NULL, weight2, NULL, NULL, NULL);
    // log_trace("diff weights :\n");
    // int idx = 0;
    // for (int j = 0; j <= order; j++) {
    //   for (int i = 0; i <= order-j; i++) {
    //     log_trace("  i = %d, j = %d, diff = %e\n",
    //               i, j, PDM_ABS(weight[idx] - weight2[idx]));
    //     idx++;
    //   }
    // }

    free(weight);
    free(weight2);
  }


  /*
   *  Regular subdivision into 4 subtriangles
   */
  double *sub_node_xyz[4];
  for (int i = 0; i < 4; i++) {
    sub_node_xyz[i] = malloc(sizeof(double) * n_node * 3);
  }
  _subdivide_bezier_triangle(3,
                             order,
                             node_xyz,
                             sub_node_xyz[0],
                             sub_node_xyz[1],
                             sub_node_xyz[2],
                             sub_node_xyz[3]);



  /*
   *  Sub-triangle u >= umin, v >= vmin, w >= wmin
   */
  double umin = 0.1;
  double vmin = 0.2;
  double wmin = 0.3;
  double *sub1_node_xyz = malloc(sizeof(double) * n_node * 3);
  _subdivide_bezier_triangle2(3,
                              order,
                              umin,
                              vmin,
                              wmin,
                              node_xyz,
                              sub1_node_xyz);





  /*
   *  Visu VTK
   */
  if (visu) {
    int *ijk_to_vtk = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                      PDM_MESH_NODAL_TRIAHO_BEZIER,
                                                      order);

    int *connec = malloc(sizeof(int) * n_node);
    for (int i = 0; i < n_node; i++) {
      connec[ijk_to_vtk[i]] = i+1;
    }

    PDM_vtk_write_std_elements_ho("bezier_tria.vtk",
                                  order,
                                  n_node,
                                  node_xyz,
                                  NULL,
                                  PDM_MESH_NODAL_TRIAHO_BEZIER,
                                  1,
                                  connec,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);

    PDM_vtk_write_std_elements("bezierCP_tria.vtk",
                               n_node,
                               node_xyz,
                               NULL,
                               PDM_MESH_NODAL_TRIA3,
                               n_tria,
                               tria_node,
                               NULL,
                               0,
                               NULL,
                               NULL);



    char filename[999];
    for (int i = 0; i < 4; i++) {
      sprintf(filename, "bezier_subtria%d.vtk", i);
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    n_node,
                                    sub_node_xyz[i],
                                    NULL,
                                    PDM_MESH_NODAL_TRIAHO_BEZIER,
                                    1,
                                    connec,
                                    NULL,
                                    0,
                                    NULL,
                                    NULL);

      sprintf(filename, "bezierCP_subtria%d.vtk", i);
      PDM_vtk_write_std_elements(filename,
                                 n_node,
                                 sub_node_xyz[i],
                                 NULL,
                                 PDM_MESH_NODAL_TRIA3,
                                 n_tria,
                                 tria_node,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }

    for (int i = 0; i < 3; i++) {
      sprintf(filename, "bezier_sub3tria%d.vtk", i);
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    n_node,
                                    sub3_node_xyz[i],
                                    NULL,
                                    PDM_MESH_NODAL_TRIAHO_BEZIER,
                                    1,
                                    connec,
                                    NULL,
                                    0,
                                    NULL,
                                    NULL);

      sprintf(filename, "bezierCP_sub3tria%d.vtk", i);
      PDM_vtk_write_std_elements(filename,
                                 n_node,
                                 sub3_node_xyz[i],
                                 NULL,
                                 PDM_MESH_NODAL_TRIA3,
                                 n_tria,
                                 tria_node,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }


    double pts[9];
    int color[3] = {0, 1, 2};
    PDM_ho_bezier_de_casteljau_triangle(3, order, umin, vmin, node_xyz, &pts[0], NULL, NULL, NULL);
    PDM_ho_bezier_de_casteljau_triangle(3, order, 1 - vmin - wmin, vmin, node_xyz, &pts[3], NULL, NULL, NULL);
    PDM_ho_bezier_de_casteljau_triangle(3, order, umin, 1 - umin - wmin, node_xyz, &pts[6], NULL, NULL, NULL);
    PDM_vtk_write_point_cloud("bezier_sub1tria_corners.vtk",
                              3,
                              pts,
                              NULL,
                              color);


    PDM_vtk_write_std_elements_ho("bezier_sub1tria.vtk",
                                  order,
                                  n_node,
                                  sub1_node_xyz,
                                  NULL,
                                  PDM_MESH_NODAL_TRIAHO_BEZIER,
                                  1,
                                  connec,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);

    PDM_vtk_write_std_elements("bezierCP_sub1tria.vtk",
                               n_node,
                               sub1_node_xyz,
                               NULL,
                               PDM_MESH_NODAL_TRIA3,
                               n_tria,
                               tria_node,
                               NULL,
                               0,
                               NULL,
                               NULL);
    free(connec);



    FILE *f = fopen("bezier_deriv.vtk", "w");
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "bezier_deriv\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS 3 double\n");
    for (int i = 0; i < 3; i++) {
      fprintf(f, "%f %f %f\n", p[0], p[1], p[2]);
    }
    fprintf(f, "CELLS 3 6\n1 0\n1 1\n1 2\n");
    fprintf(f, "CELL_TYPES 3\n1\n1\n1\n");
    fprintf(f, "POINT_DATA 3\n");
    fprintf(f, "VECTORS vector double\n");
    fprintf(f, "%f %f %f\n", dpdu[0], dpdu[1], dpdu[2]);
    fprintf(f, "%f %f %f\n", dpdv[0], dpdv[1], dpdv[2]);
    fprintf(f, "%f %f %f\n", normal[0], normal[1], normal[2]);
    fclose(f);
  }



  /*
   *  Free memory
   */
  free(node_uv);
  free(node_xyz);
  free(tria_node);
  for (int i = 0; i < 3; i++) {
    free(sub3_node_xyz[i]);
  }
  for (int i = 0; i < 4; i++) {
    free(sub_node_xyz[i]);
  }
  free(sub1_node_xyz);

  PDM_MPI_Finalize();

  return 0;
}


#undef ij2idx
