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











/**
 * Read command arguments
 **/

static void
_read_args (int      argc,
            char   **argv,
            int     *type)
{
  int i = 1;

  /* Parse and check command line */
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      exit(1);
    }

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        exit(1);
      }
      else {
        *type = atoi(argv[i]);
      }
    }

    else {
      exit(1);
    }

    i++;
  }
}



static void robust_normal
(
 const double *a,
 const double *b,
 const double *c,
 double       *normal
 )
{
  double u[2], v[2], w[2];

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (i+2)%3;

    u[0] = a[j]; u[1] = a[k];
    v[0] = b[j]; v[1] = b[k];
    w[0] = c[j]; w[1] = c[k];

    normal[i] = PDM_predicate_orient2d (u, v, w);
  }
}


static void print_diff_vector
(
 const double *a,
 const double *b
 )
{
  double err[3];
  for (int i = 0; i < 3; i++) {
    err[i] = a[i] - b[i];
  }

  double ma = PDM_MODULE(a);
  double mb = PDM_MODULE(b);

  printf("err = %2.2e %2.2e %2.2e \t (rel: %2.2e)\n",
         err[0], err[1], err[2], PDM_ABS(ma - mb)/mb);
}


/**
 * Main
 **/

int main (int argc, char *argv[])
{
  PDM_predicate_exactinit();

  int type = 0;
  _read_args (argc,
              argv,
              &type);

  const double _c = 0.866025403784438646763723170753;
  const double _s = 0.5;


  double a[3], b[3], c[3];
  for (int i = 0; i < 3; i++) {
    a[i] = 0.;//(double) rand() / (double) RAND_MAX;
  }


  for (int mag = 0; mag < 11; mag++) {

    double r = pow(10., mag);

    printf("\nr = %e\n", r);

    if (type == 0) {
  // Needle
      b[0] = a[0] + _c;
      b[1] = a[1] + _s;
      b[2] = a[2];

      c[0] = a[0] - r*_s;
      c[1] = a[1] + r*_c;
      c[2] = a[2];
    } else {
  // Sliver
      b[0] = a[0] + _c;
      b[0] += 0.5*r*_s;
      b[1] = a[1] + _s;
      b[1] -= 0.5*r*_c;
      b[2] = a[2];

      c[0] = a[0] + _c;
      c[0] -= 0.5*r*_s;
      c[1] = a[1] + _s;
      c[1] += 0.5*r*_c;
      c[2] = a[2];
    }


    double normal[3] = {0., 0., r};
    double u[3], v[3], w[3];


    for (int i = 0; i < 3; i++) {
      u[i] = b[i] - a[i];
      v[i] = c[i] - a[i];
    }

    PDM_CROSS_PRODUCT (w, u, v);

    print_diff_vector (w, normal);






    for (int i = 0; i < 3; i++) {
      u[i] = c[i] - b[i];
      v[i] = a[i] - b[i];
    }

    PDM_CROSS_PRODUCT (w, u, v);

    print_diff_vector (w, normal);





    for (int i = 0; i < 3; i++) {
      u[i] = a[i] - c[i];
      v[i] = b[i] - c[i];
    }

    PDM_CROSS_PRODUCT (w, u, v);

    print_diff_vector (w, normal);





    robust_normal (a, b, c, w);
    print_diff_vector (w, normal);

    robust_normal (b, c, a, w);
    print_diff_vector (w, normal);

    robust_normal (c, a, b, w);
    print_diff_vector (w, normal);


  // double vtx[9] = {
  //   a[0], a[1], a[2],
  //   b[0], b[1], b[2],
  //   c[0], c[1], c[2]};

  // double len[3];
  // for (int i = 0; i < 3; i++) {
  //   int j = (i+1)%3;
  //   for (int k = 0; k < 3; k++) {
  //     u[k] = vtx[3*j + k] - vtx[3*i + k];
  //   }

  //   len[i] = PDM_MODULE(u);
  // }

  // PDM_sort_double (len, NULL, 3);

  // double area = 0.25*sqrt((len[2] + (len[1] + len[0])) *
  //                         (len[0] - (len[2] - len[1])) *
  //                         (len[0] + (len[2] - len[1])) *
  //                         (len[2] + (len[1] - len[0])));

  // printf("err area = %e\n", area - 0.5*r);
  }

  return 0;
}
