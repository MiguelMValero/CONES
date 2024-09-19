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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_predicate.h"
#include "pdm_sort.h"
#include "pdm_mesh_nodal.h"
#include "pdm_vtk.h"

static void _compute_orthogonal_complement
(
 double *w,
 double *u,
 double *v
 )
{
  if (PDM_ABS(w[0]) > PDM_ABS(w[1])) {
    double imag = 1. / sqrt(w[0]*w[0] + w[2]*w[2]);
    u[0] = -w[2] * imag;
    u[1] =  0.;
    u[2] =  w[0] * imag;
  }
  else {
    double imag = 1. / sqrt(w[1]*w[1] + w[2]*w[2]);
    u[0] =  0;
    u[1] =  w[2] * imag;
    u[2] = -w[1] * imag;
  }

  PDM_CROSS_PRODUCT (v, w, u);

  // printf("%e %e %e    %e %e %e\n",
  //        PDM_DOT_PRODUCT(u, u) - 1.,
  //        PDM_DOT_PRODUCT(v, v) - 1.,
  //        PDM_DOT_PRODUCT(w, w) - 1.,
  //        PDM_DOT_PRODUCT(u, v),
  //        PDM_DOT_PRODUCT(v, w),
  //        PDM_DOT_PRODUCT(w, u));


  // double error = -0.5*PDM_DOT_PRODUCT (u, v);
  // double u_ort[3], v_ort[3], w_ort[3];
  // for (int i = 0; i < 3; i++) {
  //   u_ort[i] = fma(error, v[i], u[i]);
  //   v_ort[i] = fma(error, u[i], v[i]);
  // }

  // PDM_CROSS_PRODUCT (w_ort, u_ort, v_ort);

  // double uu = PDM_DOT_PRODUCT (u, u);
  // double vv = PDM_DOT_PRODUCT (v, v);
  // double ww = PDM_DOT_PRODUCT (w, w);
  // for (int i = 0; i < 3; i++) {
  //   u[i] = 0.5*(3. - uu)*u_ort[i];
  //   v[i] = 0.5*(3. - vv)*v_ort[i];
  //   w[i] = 0.5*(3. - ww)*w_ort[i];
  // }


  // printf("%e %e %e    %e %e %e\n",
  //        PDM_DOT_PRODUCT(u, u) - 1.,
  //        PDM_DOT_PRODUCT(v, v) - 1.,
  //        PDM_DOT_PRODUCT(w, w) - 1.,
  //        PDM_DOT_PRODUCT(u, v),
  //        PDM_DOT_PRODUCT(v, w),
  //        PDM_DOT_PRODUCT(w, u));
}


static void _random_rotation
(
const int    n,
double       coord[]
)
{
  /*// Random unit quaternion
  double q[4], lq = 0.;
  for (int i = 0; i < 4; i++) {
    q[i] = 2. * (double) rand() / (double) RAND_MAX - 1.;
    lq += q[i]*q[i];
  }

  lq = 1. / sqrt(lq);
  for (int i = 0; i < 4; i++) {
    q[i] *= lq;
  }

  double a = q[0], b = q[1], c = q[2], d = q[3];

  // Quaternion to rotation matrix
  double r[3][3] = {
    {a*a + b*b - c*c - d*d, 2*(b*c - a*d),         2*(a*c + b*d)},
    {2*(a*d + b*c),         a*a - b*b + c*c - d*d, 2*(c*d - a*b)},
    {2*(b*d - a*c),         2*(a*b + c*d),         a*a - b*b - c*c + d*d}
  };*/
  double r[3][3];
  for (int i = 0; i < 3; i++) {
    r[2][i] = 2. * (double) rand() / (double) RAND_MAX - 1.;
  }
  double imag = 1. / PDM_MODULE (r[2]);
  for (int i = 0; i < 3; i++) {
    r[2][i] *= imag;
  }
  _compute_orthogonal_complement (r[2], r[0], r[1]);


  // Rotate
  for (int k = 0; k < n; k++) {

    double x[3] = {coord[3*k], coord[3*k+1], coord[3*k+2]};

    for (int i = 0; i < 3; i++) {
      coord[3*k+i] = 0.;
      for (int j = 0; j < 3; j++) {
        coord[3*k+i] += r[i][j] * x[j];
      }
    }
  }
}


// static void robust_normal
// (
//  const double *a,
//  const double *b,
//  const double *c,
//  double       *normal
//  )
// {
//   double u[2], v[2], w[2];

//   for (int i = 0; i < 3; i++) {
//     int j = (i+1)%3;
//     int k = (i+2)%3;

//     u[0] = a[j]; u[1] = a[k];
//     v[0] = b[j]; v[1] = b[k];
//     w[0] = c[j]; w[1] = c[k];

//     normal[i] = PDM_predicate_orient2d (u, v, w);
//   }
// }

static double
difference_of_products(double a, double b, double c, double d)
{
  double cd  = c * d;
  double err = fma(-c, d,  cd);
  double dop = fma( a, b, -cd);
  return dop + err;
}

static void robust_normal
(
 const double *a,
 const double *b,
 const double *c,
 double       *normal
 )
{
  double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};

  normal[0] = difference_of_products(u[1], v[2], u[2], v[1]);
  normal[1] = difference_of_products(u[2], v[0], u[0], v[2]);
  normal[2] = difference_of_products(u[0], v[1], u[1], v[0]);
}

static void
_triangle_normal
(
 const double *a,
 const double *b,
 const double *c,
 const int     method,
 double       *normal
 )
{
  double vec[9] = {
    b[0] - a[0], b[1] - a[1], b[2] - a[2],
    c[0] - b[0], c[1] - b[1], c[2] - b[2],
    a[0] - c[0], a[1] - c[1], a[2] - c[2]
  };

  switch (method) {
    case 1: {

      double lengths[3] = {
        PDM_MODULE (vec+3),
        PDM_MODULE (vec+6),
        PDM_MODULE (vec  )
      };

      int order[3] = {0, 1, 2};

      PDM_sort_double (lengths, order, 3);

      /*int i1 = order[1];
      int i2 = order[2];

      if (i1 == (i2+1)%3) {
        int tmp = i1;
        i1 = i2;
        i2 = tmp;
      }*/
      // int i3 = (i2+1)%3;
      // printf("  i1: %e %e %e (%e)\n", vec[3*i1], vec[3*i1+1], vec[3*i1+2], PDM_MODULE(vec + 3*i1));
      // printf("  i2: %e %e %e (%e)\n", vec[3*i2], vec[3*i2+1], vec[3*i2+2], PDM_MODULE(vec + 3*i2));
      // printf("  i3: %e %e %e (%e)\n", vec[3*i3], vec[3*i3+1], vec[3*i3+2], PDM_MODULE(vec + 3*i3));

      //PDM_CROSS_PRODUCT (normal, vec + 3*i1, vec + 3*i2);
      int imax = order[2];
      const double *p, *q, *r;
      if (imax == 0) {
        p = a; q = b; r = c;
      } else if (imax == 1) {
        p = b; q = c; r = a;
      } else {
        p = c; q = a; r = b;
      }
      robust_normal (p, q, r, normal);
      break;}

    default:
    PDM_CROSS_PRODUCT (normal, vec, vec+3);
  }
}




static void _kahan_sum_vector
(
 const int     n,
 const double *input,
 double       *sum
 )
{
  for (int i = 0; i < 3; i++) {
    sum[i] = 0.;
  }

  double c[3] = {0.};

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 3; j++) {
      double y = input[3*i + j] - c[j];
      double t = sum[j] + y;
      c[j] = (t - sum[j]) - y;
      sum[j] = t;
    }
  }
}




/**
 * Read command arguments
 **/

static void
_read_args (int      argc,
            char   **argv,
            double  *length,
            int     *sort_faces,
            int     *scale,
            int     *normalize)
{
  int i = 1;

  /* Parse and check command line */
  while (i < argc) {

    if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        exit(1);
      }
      else {
        *length = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-sort") == 0) {
      *sort_faces = 1;
    }

    else if (strcmp(argv[i], "-scale") == 0) {
      *scale = 1;
    }

    else if (strcmp(argv[i], "-norm") == 0) {
      *normalize = 1;
    }

    else {
      exit(1);
    }

    i++;
  }
}



#define FMT "%.16lf" //"%20.16e"


/**
 * Main
 **/

int main (int argc, char *argv[])
{
  PDM_predicate_exactinit();

  double h          = 1e3;
  int    sort_faces = 0;
  int    scale_sum  = 0;
  int    normalize  = 0;
  _read_args (argc,
              argv,
              &h,
              &sort_faces,
              &scale_sum,
              &normalize);

  // Axis-aligned tetrahedron
  double vtx_coord[12] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, h
  };

  double face_normal[12] = {
    0,  0, -1,
    0, -h,  0,
    h,  h,  1,
   -h,  0,  0
  };


  double _coord[24];
  memcpy (_coord,      vtx_coord,   12*sizeof(double));
  memcpy (_coord + 12, face_normal, 12*sizeof(double));


  if (normalize) {
    //double factor = 1. / h;
    double factor = 1. / sqrt(1 * h);
    for (int i = 0; i < 24; i++) {
      _coord[i] *= factor;
    }
  }

  // Apply random rotation
  _random_rotation (8, _coord);

  memcpy (vtx_coord,   _coord,      12*sizeof(double));
  memcpy (face_normal, _coord + 12, 12*sizeof(double));


  // for (int i = 0; i < 4; i++) {
  //   printf("%f %f %f\n", vtx_coord[3*i], vtx_coord[3*i+1], vtx_coord[3*i+2]);
  // }

  // Faces
  const int n_face = 4;
  int face_vtx[12] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4
  };



  printf("h = %e\n", h);

  for (int method = 0; method < 3; method++) {
    printf("method %d:\n", method);

    double face_vec[12];
    double mag_face[4];

    for (int i = 0; i < n_face; i++) {
      int *_face_vtx = face_vtx + 3*i;

      if (method < 2) {
        _triangle_normal (vtx_coord + 3*(_face_vtx[0] - 1),
                          vtx_coord + 3*(_face_vtx[1] - 1),
                          vtx_coord + 3*(_face_vtx[2] - 1),
                          method,
                          face_vec + 3*i);
      } else {
        for (int j = 0; j < 3; j++) {
          face_vec[3*i+j] = face_normal[3*i+j];
        }
      }

      mag_face[i] = PDM_MODULE (face_vec + 3*i);

      printf("  face_vec[%d] = "FMT" "FMT" "FMT"   mag = "FMT"\n", i,
             face_vec[3*i], face_vec[3*i+1], face_vec[3*i+2], mag_face[i]);
    }

    int order[4] = {0, 1, 2, 3};
    if (sort_faces) {
      PDM_sort_double (mag_face, order, n_face);
    }

    double cell_vec[3] = {0., 0., 0.};
    if (scale_sum) {
      double scale = 0.;
      for (int i = 0; i < n_face; i++) {
        scale += mag_face[i];
      }

      scale = 1. / scale;

      for (int i = 0; i < n_face; i++) {
        int iface = order[i];
        for (int j = 0; j < 3; j++) {
        //cell_vec[j] += scale * face_vec[3*iface + j];
          cell_vec[j] = fma(scale, face_vec[3*iface + j], cell_vec[j]);
        }
      }
    }

    else {
      _kahan_sum_vector (n_face, face_vec, cell_vec);
    }

    printf("  cell_vec = "FMT" "FMT" "FMT"\n", cell_vec[0], cell_vec[1], cell_vec[2]);

    double mag = PDM_MODULE (cell_vec);

    printf("  mag(cell_vec) = %e (= %e * h)\n", mag, mag / h);
  }






  // Visu
  if (0) {
    int elt_vtx[4] = {1, 2, 3, 4};
    PDM_vtk_write_std_elements ("test_aniso_tetra.vtk",
                                4,
                                vtx_coord,
                                NULL,
                                PDM_MESH_NODAL_TETRA4,
                                1,
                                elt_vtx,
                                NULL,
                                0,
                                NULL,
                                NULL);
  }

  return 0;
}
