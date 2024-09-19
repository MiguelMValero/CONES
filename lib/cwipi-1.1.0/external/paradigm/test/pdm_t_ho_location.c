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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_ho_location.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_array.h"
#include "pdm_ho_basis.h"

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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_pts,
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *order,
           int           *t_elt,
           int           *use_newton,
           int           *visu)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *n_pts = (PDM_g_num_t) _n;
      }
    }

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
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-newton") == 0) {
      *use_newton = 1;
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static double
_eval_field
(
 const double x,
 const double y,
 const double z,
 const int    order
 )
{
  int _order = PDM_MAX(1, (int) floor(sqrt(order)));

  double f = x - y + z - 1;

  if (_order == 1) {
    return f;
  }

  f += x*x - y*y + z*z - x*y + y*z - z*x;

  if (_order == 2) {
    return f;
  }

  f +=
  x*x*x - y*y*y + z*z*z -
  x*x*y + x*x*z -
  y*y*x + y*y*z -
  z*z*x + z*z*y;

  return f;
}



static void
_eval_deformation
(
 const int     order,
 const double  x,
 const double  y,
 const double  z,
       double *dx,
       double *dy,
       double *dz
 )
{
  int _order = PDM_MAX(1, (int) ceil(sqrt(order)));

  *dx = 0.4*y;
  *dy = 0.4*z;
  *dz = 0.4*x;

  if (_order == 1) {
    return;
  }

  *dx += 0.2*y*y;
  *dy += 0.2*z*z;
  *dz += 0.2*x*x;

  if (_order == 2) {
    return;
  }

  *dx += 0.1*y*y*y;
  *dy += 0.1*z*z*z;
  *dz += 0.1*x*x*x;

  return;
}



static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

/*=============================================================================
 * Private function definition
 *============================================================================*/

static inline double
_determinant_3x3
(
 const double a[3],
 const double b[3],
 const double c[3]
 )
{
  return a[0] * (b[1]*c[2] - b[2]*c[1])
    +    a[1] * (b[2]*c[0] - b[0]*c[2])
    +    a[2] * (b[0]*c[1] - b[1]*c[0]);
}

/*---------------------------------------------------------------------------
 * Solve the equation "A.x = b" with Cramer's rule.
 *
 * parameters:
 *   A[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   PDM_FALSE if matrix is singular, PDM_TRUE otherwise
 *----------------------------------------------------------------------------*/

static int
_solve_3x3(double  A[3][3],
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
    -   A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    +   A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
          - b[1]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
          + b[2]*(A[0][1]*A[1][2] - A[1][1]*A[0][2])) * det_inv;

  x1 = (  A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2])
          - A[1][0]*(b[0]*A[2][2] - b[2]*A[0][2])
          + A[2][0]*(b[0]*A[1][2] - b[1]*A[0][2])) * det_inv;

  x2 = (  A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1])
          - A[1][0]*(A[0][1]*b[2] - A[2][1]*b[0])
          + A[2][0]*(A[0][1]*b[1] - A[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return PDM_TRUE;
}


static int
_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
 double                     uvw[3]
 )
{
  double duvw     = 1e-8;
  double inv_duvw = 1./duvw;

  int i, j, n_node, iter;
  const int max_iter = 30;
  const double tolerance2 = tolerance * tolerance;
  double dist;
  double a[3][3], b[3], x[3];

  // log_trace("\ntarget: %f %f %f\n", point_coords[0], point_coords[1], point_coords[2]);


  /* Get number of nodes */
  n_node = PDM_Mesh_nodal_n_vertices_element (elt_type, order);
  double *weight     = malloc(sizeof(double) * n_node);
  double *dweight_du = malloc(sizeof(double) * n_node);
  double *dweight_dv = malloc(sizeof(double) * n_node);
  double *dweight_dw = malloc(sizeof(double) * n_node);

  double *dw[3] = {dweight_du, dweight_dv, dweight_dw};

  // assert (elt_type == PDM_MESH_NODAL_PYRAMID5 ||
  //         elt_type == PDM_MESH_NODAL_PRISM6   ||
  //         elt_type == PDM_MESH_NODAL_HEXA8);

  /* Use Newton-method to determine parametric coordinates and shape function */
  for (i = 0; i < 3; i++) {
    uvw[i] = 0.5;
  }

  double _uvw[12];
  double *_weight = malloc(sizeof(double) * n_node * 4);

  for (iter = 0; iter < max_iter; iter++) {

    memcpy(_uvw,   uvw, sizeof(double) * 3);
    memcpy(_uvw+3, uvw, sizeof(double) * 3);
    memcpy(_uvw+6, uvw, sizeof(double) * 3);
    memcpy(_uvw+9, uvw, sizeof(double) * 3);
    _uvw[3]  += duvw;
    _uvw[7]  += duvw;
    _uvw[11] += duvw;

    PDM_ho_basis(elt_type,
                 order,
                 n_node,
                 4,
                 _uvw,
                 _weight);

    for (i = 0; i < n_node; i++) {
      weight[i] = _weight[i];
      for (j = 0; j < 3; j++) {
        dw[j][i] = (_weight[(j+1)*n_node + i] - _weight[i]) * inv_duvw;
      }
    }


    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];


    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        a[i][j] = 0.0;
      }
    }

    for (i = 0; i < n_node; i++) {

      b[0] += (weight[i] * vertex_coords[3*i]);
      b[1] += (weight[i] * vertex_coords[3*i+1]);
      b[2] += (weight[i] * vertex_coords[3*i+2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[j][i] * vertex_coords[3*i]);
        a[1][j] -= (dw[j][i] * vertex_coords[3*i+1]);
        a[2][j] -= (dw[j][i] * vertex_coords[3*i+2]);
      }

    }
    // log_trace("current: %f %f %f,  res_xyz = %e\n",
    //           b[0] + point_coords[0], b[1] + point_coords[1], b[2] + point_coords[2],
    //           PDM_MODULE(b));

    if (_solve_3x3(a, b, x) == PDM_FALSE) {
      printf("_compute_uvw: singular matrix\n");
      return PDM_FALSE;
    }

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    // log_trace("iter %d, dist = %e\n", iter, dist);

    if (dist <= tolerance2) {
      // log_trace("converged :) %e %e %e\n", uvw[0], uvw[1], uvw[2]);

      free(weight    );
      free(dweight_du);
      free(dweight_dv);
      free(dweight_dw);
      free(_weight);
      return 1;
    }

  }

  free(weight    );
  free(dweight_du);
  free(dweight_dv);
  free(dweight_dw);
  free(_weight);
  return 0;
}




static int _accept_uvw
(
 const double               *uvw,
 const PDM_Mesh_nodal_elt_t  t_elt
 )
{
  // PDM_UNUSED(uvw);
  // PDM_UNUSED(t_elt);
  // return 1;

  switch (t_elt) {

    case PDM_MESH_NODAL_TRIAHO: {
      return (uvw[0] >= 0. &&
              uvw[1] >= 0. &&
              uvw[0] + uvw[1] <= 1.);
      break;
    }
    case PDM_MESH_NODAL_QUADHO: {
      return (uvw[0] >= 0. &&
              uvw[0] <= 1. &&
              uvw[1] >= 0. &&
              uvw[1] <= 1.);

      break;
    }
    case PDM_MESH_NODAL_TETRAHO: {
      return (uvw[0] >= 0. &&
              uvw[1] >= 0. &&
              uvw[2] >= 0. &&
              uvw[0] + uvw[1] + uvw[2] <= 1.);
      break;
    }
    case PDM_MESH_NODAL_PYRAMIDHO: {
      return (uvw[0] >= 0. &&
              uvw[0] <= 1. - uvw[2] &&
              uvw[1] >= 0. &&
              uvw[1] <= 1. - uvw[2] &&
              uvw[2] >= 0. &&
              uvw[2] <= 1.);
      break;
    }
    case PDM_MESH_NODAL_PRISMHO: {
      return (uvw[0] >= 0. &&
              uvw[1] >= 0. &&
              uvw[2] >= 0. &&
              uvw[2] <= 1. &&
              uvw[0] + uvw[1] <= 1.);
      break;
    }
    case PDM_MESH_NODAL_HEXAHO: {
      return (uvw[0] >= 0. &&
              uvw[0] <= 1. &&
              uvw[1] >= 0. &&
              uvw[1] <= 1. &&
              uvw[2] >= 0. &&
              uvw[2] <= 1.);

      break;
    }
    default: {
      return 1;
    }

  }

  return 1;
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
  PDM_g_num_t          gn_pts     = 1;
  PDM_g_num_t          nx         = 2;
  PDM_g_num_t          ny         = 2;
  PDM_g_num_t          nz         = 2;
  int                  order      = 3;
  double               length     = 1.;
  int                  use_newton = 0;
  int                  visu       = 0;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXAHO;
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

  assert(t_elt > PDM_MESH_NODAL_POLY_3D);


  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &gn_pts,
             &nx,
             &ny,
             &nz,
             &order,
             (int *) &t_elt,
             &use_newton,
             &visu);

  /*
   *  Init
   */
  // struct timeval t_elaps_debut;

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(n_rank == 1);

  printf("order = %d, field: %d, deformation: %d\n",
         order,
         PDM_MAX(1, (int) floor(sqrt(order))),
         PDM_MAX(1, (int) ceil(sqrt(order))));


  // int dim = 3;
  // if (t_elt == PDM_MESH_NODAL_TRIAHO ||
  //     t_elt == PDM_MESH_NODAL_QUADHO) {
  //   dim = 2;
  // }
  int elt_dim = PDM_Mesh_nodal_elt_dim_get(t_elt);


  if (order > 3) {
    int *ijk = NULL;

    for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BARHO;
         type <= PDM_MESH_NODAL_HEXAHO;
         type++) {

      if (type == PDM_MESH_NODAL_PYRAMIDHO) continue;

      ijk = PDM_vtk_lagrange_to_ijk(type, order);
      PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                       type,
                                       order,
                                       PDM_Mesh_nodal_n_vtx_elt_get(type, order),
                                       ijk);
      free (ijk);
    }
  }


  /*
   *  Create distributed cube
   */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        nx,
                                                        ny,
                                                        nz,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);

  // PDM_dcube_nodal_gen_ordering_set (dcube,
  //                                   "PDM_HO_ORDERING_VTK");

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);



  /* Deform */
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  // double amplitude = 0.1;//0.07;
  // double frequence = 4.;

  if (1) {
    /* Polynomial deformation */
    for (int i = 0; i < dn_vtx; i++) {
      double x = (dvtx_coord[3*i    ] - 0.5) / length;
      double y = (dvtx_coord[3*i + 1] - 0.5) / length;
      double z = (dvtx_coord[3*i + 2] - 0.5) / length;

      double dx, dy, dz;
      _eval_deformation(order,
                        x, y, z,
                        &dx, &dy, &dz);

      dvtx_coord[3*i    ] += length * dx;
      dvtx_coord[3*i + 1] += length * dy;
      dvtx_coord[3*i + 2] += length * dz;
    }

    if (1) {
      /* Rotation */
      for (int i = 0; i < dn_vtx; i++) {
        double x = dvtx_coord[3*i  ];
        double y = dvtx_coord[3*i+1];
        double z = dvtx_coord[3*i+2];

        for (int j = 0; j < 3; j++) {
          dvtx_coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
        }
      }
    }
  }




  /* Isolate a single element */
  PDM_geometry_kind_t geom_kind;
  if (elt_dim == 3) {
    geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  } else if (elt_dim == 2) {
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
  } else {
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
  }

  int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  // int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, geom_kind);

  int id_section = sections_id[0];
  // const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
  // int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
  PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
  assert(PDM_DMesh_nodal_section_type_get(dmn, geom_kind, id_section) == t_elt);


  int n_node = PDM_Mesh_nodal_n_vertices_element(t_elt, order);

  int ielt = 0;
  double *node_coord = malloc(sizeof(double) * n_node * 3);
  for (int i = 0; i < n_node; i++) {
    int ivtx = (int) (dconnec[n_node*ielt + i] - 1);

    memcpy(node_coord + 3*i, dvtx_coord + 3*ivtx, sizeof(double) * 3);
  }


  /* Random point cloud */
  // int n_pts;
  // double      *pts_coord    = NULL;
  // PDM_g_num_t *pts_ln_to_gn = NULL;
  // PDM_point_cloud_gen_random(comm,
  //                            gn_pts,
  //                            0.3,
  //                            0.3,
  //                            0.3,
  //                            0.7,
  //                            0.7,
  //                            0.7,
  //                            &n_pts,
  //                            &pts_coord,
  //                            &pts_ln_to_gn);


  int n_pts = (int) gn_pts;
  double *pts_uvw_init = malloc(sizeof(double) * n_pts * elt_dim);
  double *pts_coord    = malloc(sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    // pts_coord[i] = 0.3 + 0.4*((double) rand() / (double) RAND_MAX);
    while (1) {
      for (int j = 0; j < elt_dim; j++) {
        pts_uvw_init[elt_dim*i+j] = (double) rand() / (double) RAND_MAX;
      }
      if (_accept_uvw(&pts_uvw_init[elt_dim*i], t_elt)) {
        break;
      }
    }
  }

  double *pts_weight = malloc(sizeof(double) * n_pts * n_node);
  PDM_ho_basis(t_elt,
               order,
               n_node,
               n_pts,
               pts_uvw_init,
               pts_weight);

  for (int i = 0; i < n_pts; i++) {

    double *p = pts_coord  + 3*i;
    double *w = pts_weight + n_node*i;

    p[0] = p[1] = p[2] = 0.;
    for (int j = 0; j < n_node; j++) {
      for (int k = 0; k < 3; k++) {
        p[k] += w[j] * node_coord[3*j+k];
      }
    }
  }



  /* HO location */
  double *proj_coord = malloc(sizeof(double) * n_pts * 3);
  double *pts_uvw    = malloc(sizeof(double) * n_pts * elt_dim);
  for (int i = 0; i < n_pts; i++) {
    int stat = 0;

    if (use_newton) {
      stat = _compute_uvw(t_elt,
                          order,
                          pts_coord + 3*i,
                          node_coord,
                          1e-8,
                          pts_uvw + elt_dim*i);
    }

    if (!stat) {
      double dist = PDM_ho_location(t_elt,
                                    order,
                                    n_node,
                                    node_coord,
                                    pts_coord  + 3*i,
                                    proj_coord + 3*i,
                                    pts_uvw    + elt_dim*i);
      PDM_UNUSED(dist);
      // log_trace("i = %d, dist = %e\n", i, dist);
    }
  }

  PDM_ho_basis(t_elt,
               order,
               n_node,
               n_pts,
               pts_uvw,
               pts_weight);



  double *node_uvw = malloc(sizeof(double) * n_node * elt_dim);
  PDM_ho_location_uvw_nodes(t_elt,
                            order,
                            0., 1.,
                            0., 1.,
                            0., 1.,
                            node_uvw);


  double *node_field = malloc(sizeof(double) * n_node);
  for (int i = 0; i < n_node; i++) {
    // node_field[i] = _eval_field(node_uvw[3*i  ],
    //                             node_uvw[3*i+1],
    //                             node_uvw[3*i+2],
    //                             order);
    node_field[i] = _eval_field(node_coord[3*i  ],
                                node_coord[3*i+1],
                                node_coord[3*i+2],
                                order);
  }

  double *interp_field = malloc(sizeof(double) * n_pts);
  double *proj_field   = malloc(sizeof(double) * n_pts);
  double field_max = 0.;
  double err_max   = 0.;
  for (int i = 0; i < n_pts; i++) {

    double *q = pts_coord  + 3*i;
    double *p = proj_coord + 3*i;
    // double *u = pts_uvw_init + 3*i;
    double *w = pts_weight + n_node*i;


    p[0] = p[1] = p[2] = 0.;
    for (int j = 0; j < n_node; j++) {
      for (int k = 0; k < 3; k++) {
        p[k] += w[j] * node_coord[3*j+k];
      }
    }

    double dist = 0;
    for (int k = 0; k < 3; k++) {
      dist += (p[k] - q[k]) * (p[k] - q[k]);
    }
    dist = sqrt(dist);

    // proj_field[i] = _eval_field(u[0],
    //                             u[1],
    //                             u[2],
    //                             order);
    proj_field[i] = _eval_field(p[0],
                                p[1],
                                p[2],
                                order);

    interp_field[i] = 0.;
    for (int j = 0; j < n_node; j++) {
      interp_field[i] += w[j] * node_field[j];
    }

    field_max  = PDM_MAX(field_max, PDM_ABS(proj_field[i]));
    double err = PDM_ABS(proj_field[i] - interp_field[i]);
    err_max = PDM_MAX(err_max, err);
    // log_trace("i = %d, %f %f %f, dist = %e, err = %e\n", i, p[0], p[1], p[2], dist, err);
  }
  // log_trace("max field = %e\n", field_max);

  printf("err_max = %e, relative = %e\n", err_max, err_max/field_max);


  int *connec = malloc(sizeof(int) * n_node);
  for (int i = 0; i < n_node; i++) {
    connec[i] = i + 1;
  }

  // FILE *f = fopen("ho_location_elt.mesh", "w");

  // fprintf(f, "MeshVersionFormatted 3\nDimension 3\n");

  // fprintf(f, "Vertices\n%d\n", n_node);
  // for (int i = 0; i < n_node; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     fprintf(f, "%f ", node_coord[3*i + j]);
  //   }
  //   fprintf(f, "%d\n", i+1);
  // }

  // switch (t_elt) {
  //   case PDM_MESH_NODAL_TETRAHO: {
  //     fprintf(f, "TetrahedraP%dOrdering\n%d\n", order, n_node);
  //     for (int k = 0; k <= order; k++) {
  //       for (int j = 0; j <= order-k; j++) {
  //         for (int i = 0; i <= order-k-j; i++) {
  //           fprintf(f, "%d %d %d\n", i, j, k);
  //         }
  //       }
  //     }

  //     fprintf(f, "TetrahedraP%d\n1\n", order);
  //     for (int i = 0; i < n_node; i++) {
  //       fprintf(f, "%d ", connec[i]);
  //     }
  //     fprintf(f, "1\n");
  //     break;
  //   }

  //   case PDM_MESH_NODAL_HEXAHO: {
  //     fprintf(f, "HexahedraQ%dOrdering\n%d\n", order, n_node);
  //     for (int k = 0; k <= order; k++) {
  //       for (int j = 0; j <= order; j++) {
  //         for (int i = 0; i <= order; i++) {
  //           fprintf(f, "%d %d %d\n", i, j, k);
  //         }
  //       }
  //     }

  //     fprintf(f, "HexahedraQ%d\n1\n", order);
  //     for (int i = 0; i < n_node; i++) {
  //       fprintf(f, "%d ", connec[i]);
  //     }
  //     fprintf(f, "1\n");
  //     break;
  //   }

  //   default: {

  //   }
  // }

  // fprintf(f, "End\n");
  // fclose(f);


  // f = fopen("ho_location_elt.sol", "w");
  // fprintf(f, "MeshVersionFormatted 3\nDimension 3\n");

  // fprintf(f, "SolAtVertices\n%d\n1 1\n", n_node);
  // for (int i = 0; i < n_node; i++) {
  //   fprintf(f, "%f\n", node_field[i]);
  // }

  // fprintf(f, "End\n");
  // fclose(f);





  /* Reorder */
  int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                     t_elt,
                                                     order);
  // PDM_log_trace_array_int(ijk_to_user, n_node, "ijk_to_user : ");

  for (int i = 0; i < n_node; i++) {
    connec[ijk_to_user[i]] = i + 1;
  }

  if (visu) {
    const char *field_name[] = {"field", "field0"};
    PDM_vtk_write_std_elements_ho_with_vtx_field("ho_location_elt.vtk",
                                                 order,
                                                 n_node,
                                                 node_coord,
                                                 NULL,
                                                 t_elt,
                                                 1,
                                                 connec,
                                                 NULL,
                                                 0,
                                                 NULL,
                                                 NULL,
                                                 1,
                                                 (const char   **) field_name,
                                                 (const double **) &node_field);

    const double *field[2] = {interp_field, proj_field};
    PDM_vtk_write_std_elements_double("ho_location_pts.vtk",
                                      n_pts,
                                      pts_coord,
                                      NULL,
                                      PDM_MESH_NODAL_POINT,
                                      n_pts,
                                      NULL,
                                      NULL,
                                      2,
                                      (const char   **) field_name,
                                      (const double **) field);


    PDM_vtk_write_std_elements_double("ho_location_proj.vtk",
                                      n_pts,
                                      proj_coord,
                                      NULL,
                                      PDM_MESH_NODAL_POINT,
                                      n_pts,
                                      NULL,
                                      NULL,
                                      2,
                                      (const char   **) field_name,
                                      (const double **) field);
  }
  free(connec);

  free(node_coord);
  free(pts_coord);
  // free(pts_ln_to_gn);
  free(proj_coord);
  free(pts_uvw);
  free(pts_weight);
  free(pts_uvw_init);

  free(node_uvw);
  free(node_field);
  free(interp_field);
  free(proj_field);









  PDM_dcube_nodal_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;

}

