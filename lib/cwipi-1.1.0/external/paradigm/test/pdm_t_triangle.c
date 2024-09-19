#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_triangle.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private functions
 *============================================================================*/

static inline double
_rand01
(
 void
 )
{
  return (double) rand() / (double) RAND_MAX;
}

static void
_dump_triangle
(
 const char   *filename,
       double  vtx_coord[9]
 )
{
  int connec[3] = {1, 2, 3};

  PDM_g_num_t vtx_id[3] = {1, 2, 3};

  PDM_vtk_write_std_elements(filename,
                             3,
                             vtx_coord,
                             vtx_id,
                             PDM_MESH_NODAL_TRIA3,
                             1,
                             connec,
                             NULL,
                             0, NULL, NULL);
}


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
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -post            Export vtk (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *seed,
 int           *randomize,
 double        *v0,
 double        *v1,
 double        *v2,
 int           *n_pts,
 int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *seed = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-v0") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          v0[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-v1") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          v1[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-v2") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          v2[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_pts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}






static void
_gen_config
(
 double                 tri_coord[9],
 PDM_triangle_status_t  stat,
 double                 dist2,
 double                 coord[3],
 double                 closest_point[3],
 double                 weight[3]
 )
{
  double xu[3] = {tri_coord[3] - tri_coord[0], tri_coord[4] - tri_coord[1], tri_coord[5] - tri_coord[2]};
  double xv[3] = {tri_coord[6] - tri_coord[0], tri_coord[7] - tri_coord[1], tri_coord[8] - tri_coord[2]};
  double normal[3];
  PDM_CROSS_PRODUCT(normal, xu, xv);

  double imag_normal = 1./PDM_MODULE(normal);
  for (int i = 0; i < 3; i++) {
    normal[i] *= imag_normal;
  }

  double u = _rand01();
  double v = _rand01();
  double w = _rand01();
  double isum = 1./(u + v + w);
  u *= isum;
  v *= isum;
  w *= isum;



  if (stat == PDM_TRIANGLE_INSIDE) {
    // inside
    double scale = sqrt(dist2 / PDM_DOT_PRODUCT(normal, normal));
    double r = _rand01();
    if (r < 0.5) {
      scale = -scale;
    }

    for (int i = 0; i < 3; i++) {
      closest_point[i] = w*tri_coord[i] + u*tri_coord[3+i] + v*tri_coord[6+i];
      coord[i] = closest_point[i] + scale*normal[i];
    }

    weight[0] = w;
    weight[1] = u;
    weight[2] = v;
  }

  else {
    // outside
    double r = _rand01();
    if (r < 0.5) {
      // closest point is a vertex
      int ivtx1 = (int) (6*r);
      int ivtx2 = (ivtx1+1)%3;
      int ivtx3 = (ivtx1+2)%3;

      double u1[3], u2[3];
      for (int i = 0; i < 3; i++) {
        u1[i] = tri_coord[3*ivtx1+i] - tri_coord[3*ivtx3+i];
        u2[i] = tri_coord[3*ivtx2+i] - tri_coord[3*ivtx1+i];
      }
      double imag_u1 = 1./PDM_MODULE(u1);
      double imag_u2 = 1./PDM_MODULE(u2);
      for (int i = 0; i < 3; i++) {
        u1[i] *= imag_u1;
        u2[i] *= imag_u2;
      }

      double v1[3];
      PDM_CROSS_PRODUCT(v1, u1, normal);
      double v2[3];
      PDM_CROSS_PRODUCT(v2, u2, normal);

      double theta_max = acos(PDM_MAX(-1, PDM_MIN(1, PDM_DOT_PRODUCT(v1, v2))));
      double theta = theta_max * _rand01();
      double phi   = asin(1 - 2*_rand01());

      for (int i = 0; i < 3; i++) {
        closest_point[i] = tri_coord[3*ivtx1+i];
        coord[i] = closest_point[i] + sqrt(dist2) * (cos(phi)*(cos(theta)*v1[i] + sin(theta)*u1[i]) + sin(phi)*normal[i]);
      }


      weight[ivtx1] = 1;
      weight[ivtx2] = 0;
      weight[ivtx3] = 0;

    }
    else {
      // closest point is on an edge
      double t = _rand01();
      int iedge = (int) (6*(r - 0.5));
      int ivtx1 = iedge;
      int ivtx2 = (iedge+1)%3;
      int ivtx3 = (iedge+2)%3;

      double vec[3];
      for (int i = 0; i < 3; i++) {
        vec[i] = tri_coord[3*ivtx2+i] - tri_coord[3*ivtx1+i];
      }
      double imag_vec = 1./PDM_MODULE(vec);
      for (int i = 0; i < 3; i++) {
        vec[i] *= imag_vec;
      }

      double out[3];
      PDM_CROSS_PRODUCT(out, vec, normal);

      double a = PDM_PI*(_rand01() - 0.5);

      for (int i = 0; i < 3; i++) {
        closest_point[i] = (1-t)*tri_coord[3*ivtx1+i] + t*tri_coord[3*ivtx2+i];
        // coord[i] = closest_point[i] + scale*out[i];
        coord[i] = closest_point[i] + sqrt(dist2)*(cos(a)*out[i] + sin(a)*normal[i]);
      }

      weight[ivtx1] = 1 - t;
      weight[ivtx2] = t;
      weight[ivtx3] = 0;
    }
  }


}





/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  double eps = DBL_EPSILON;

  // printf("(1. + 1.0*eps) - 1. = %20.16e\n", (1. + 1.0*eps) - 1.);
  // printf("(1. + 0.6*eps) - 1. = %20.16e\n", (1. + 0.6*eps) - 1.);
  // printf("(1. + 0.5*eps) - 1. = %20.16e\n", (1. + 0.5*eps) - 1.);

  double scale = 1.;
  for (int i = 0; i < 17; i++) {
    // printf("(%20.16f + %20.16f) - %20.16f = %20.16f\n",
    //        scale, scale*eps, scale*eps,
    //        (scale + scale*eps) - scale*eps);
    // printf("(%20.16f + %20.16f) - %20.16f = %20.16f\n",
    //        scale, scale*eps, scale,
    //        (scale + scale*eps) - scale);
    printf("(%20.16f + %20.16f) = %20.16f\n",
           scale, scale*eps,
           (scale + scale*eps));
    scale *= 10.;
  }

  // printf("(1.e16 +1. ) - 1.e16 = %20.16e\n", (1.e16 +1. ) - 1.e16 );
  // printf("(1.e16 +2. ) - 1.e16 = %20.16e\n", (1.e16 +2. ) - 1.e16 );

  return 0;
  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);

  int    seed      = 0;
  int    randomize = 0;
  double v0[3]     = {0, 0, 0};
  double v1[3]     = {1, 0, 0};
  double v2[3]     = {0, 1, 0};
  int    n_pts     = 1;
  int    visu      = 0;
  _read_args(argc,
             argv,
             &seed,
             &randomize,
             v0,
             v1,
             v2,
             &n_pts,
             &visu);

  if (seed == -1) {
    seed = time(NULL);
  }

  srand(seed);

  double vtx_coord[9] = {
    v0[0], v0[1], v0[2],
    v1[0], v1[1], v1[2],
    v2[0], v2[1], v2[2]
  };

  if (randomize) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        vtx_coord[3*i+j] = _rand01();
      }
    }
  }


  double *pts_coord = malloc(sizeof(double) * n_pts * 3);
  double *pts_proj  = malloc(sizeof(double) * n_pts * 3);
  double *pts_error = malloc(sizeof(double) * n_pts);
  double *pts_stat  = malloc(sizeof(double) * n_pts);
  double *pts_dist  = malloc(sizeof(double) * n_pts);
  PDM_g_num_t *pts_g_num = malloc(sizeof(PDM_g_num_t) * n_pts);

  int count_ok = 0;

  for (int i = 0; i < n_pts; i++) {

    pts_g_num[i] = i+1;
    double *p = pts_coord + 3*i;


    PDM_triangle_status_t true_stat;
    double r = _rand01();
    if (r < 0.5) {
      true_stat = PDM_TRIANGLE_OUTSIDE;
    }
    else {
      true_stat = PDM_TRIANGLE_INSIDE;
    }

    double true_min_dist2 = _rand01();
    // true_min_dist2 *= true_min_dist2;

    double true_closest_point[3], true_weight[3];

    _gen_config(vtx_coord,
                true_stat,
                true_min_dist2,
                p,
                true_closest_point,
                true_weight);

    double s = _rand01();
    if (s < 0.02 && true_stat == PDM_TRIANGLE_OUTSIDE) {
      // log_trace("snap point %d\n", i+1);
      true_min_dist2 = 0.;
      true_stat = PDM_TRIANGLE_INSIDE;
      memcpy(p, true_closest_point, sizeof(double)*3);
    }

    // double closest_point[3];
    double *closest_point = pts_proj + 3*i;
    double min_dist2;
    double weight[3];
    PDM_triangle_status_t stat = PDM_triangle_evaluate_position(p,
                                                                vtx_coord,
                                                                closest_point,
                                                                &min_dist2,
                                                                weight);
    pts_stat[i] = stat;
    pts_dist[i] = sqrt(min_dist2);

    double err_weight[3] = {
      PDM_ABS(weight[0] - true_weight[0]),
      PDM_ABS(weight[1] - true_weight[1]),
      PDM_ABS(weight[2] - true_weight[2])
    };

    if ((stat != true_stat && min_dist2 > 1e-15) ||
        min_dist2 - true_min_dist2 > 1e-15) {
      log_trace("!!! point %d: %f %f %f, status = %d / %d, min_dist2 = %e / %e, err_weight = %e %e %e\n",
                i+1, p[0], p[1], p[2],
                // true_weight[0],
                // true_weight[1],
                // true_weight[2],
                stat,
                true_stat,
                min_dist2,
                true_min_dist2,
                err_weight[0],
                err_weight[1],
                err_weight[2]);


      double err = 0.;
      for (int j = 0; j < 3; j++) {
        err += (closest_point[j] - true_closest_point[j])*(closest_point[j] - true_closest_point[j]);
      }
      err = sqrt(err);
      log_trace("  true_closest_point = %f %f %f, error = %e\n",
                true_closest_point[0], true_closest_point[1], true_closest_point[2],
                err);
      pts_error[i] = 1;
    }
    else {
      pts_error[i] = 0;
      count_ok++;
    }

  }

  if (visu) {
    const char   *field_name[3] = {"error", "status", "distance"};
    const double *field_value[3] = {pts_error, pts_stat, pts_dist};

    PDM_vtk_write_std_elements_double("triangle_points.vtk",
                                      n_pts,
                                      pts_coord,
                                      pts_g_num,
                                      PDM_MESH_NODAL_POINT,
                                      n_pts,
                                      NULL,
                                      pts_g_num,
                                      3,
                                      field_name,
                    (const double **) field_value);

    PDM_vtk_write_std_elements_double("triangle_proj.vtk",
                                      n_pts,
                                      pts_proj,
                                      pts_g_num,
                                      PDM_MESH_NODAL_POINT,
                                      n_pts,
                                      NULL,
                                      pts_g_num,
                                      3,
                                      field_name,
                    (const double **) field_value);

    _dump_triangle("triangle_triangle.vtk", vtx_coord);
  }


  free(pts_coord);
  free(pts_proj );
  free(pts_error);
  free(pts_stat );
  free(pts_dist );
  free(pts_g_num);

  printf("%d OK / %d\n", count_ok, n_pts);



  PDM_MPI_Finalize();

  return n_pts - count_ok;
}
