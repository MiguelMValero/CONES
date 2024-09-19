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
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_array.h"
#include "pdm_predicate.h"

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

/**
 * Box-Muller algorithm (polar method) for generating
 * a random number with unit normal distribution
 */
static double
_random_number_normal
(
 void
 )
{
  for (int itry = 0; itry < 100; itry++) {
    double x = 2*_rand01() - 1.;
    double y = 2*_rand01() - 1.;
    double r2 = x*x + y*y;

    if (r2 >= 0 && r2 <= 1) {
      return x*sqrt(-2*log(r2)/r2);
    }
  }

  return 0.;
}


static inline void
_normalize
(
 double a[3]
 )
 {
  double mag = PDM_MODULE(a);
  if (mag >= 0) {
    mag = 1./mag;
    for (int i = 0; i < 3; i++) {
      a[i] *= mag;
    }
  }
 }

static void
_random_orthonormal_basis
(
 double axis[3][3]
 )
{
  double *u = axis[0];
  double *v = axis[1];
  double *w = axis[2];

  int    imin = -1;
  double umin = HUGE_VAL;
  for (int i = 0; i < 3; i++) {
    u[i] = 2*_rand01() - 1.;
    v[i] = 0.;
    if (PDM_ABS(u[i]) < umin) {
      imin = i;
      umin = PDM_ABS(u[i]);
    }
  }
  _normalize(u);

  v[imin] = 1.;
  for (int i = 0; i < 3; i++) {
    v[i] -= u[imin]*u[i];
  }

  _normalize(v);

  PDM_CROSS_PRODUCT(w, u, v);
}


static void
_gen_random_tetrahedron
(
 double    radii[3],
 double    axes[3][3],
 int      *n_vtx,
 double  **vtx_coord,
 int      *n_face,
 int     **face_vtx_idx,
 int     **face_vtx,
 int      *n_cell,
 int     **cell_face_idx,
 int     **cell_face
 )
{
  *n_vtx  = 4;
  *n_face = 4;
  *n_cell = 1;

  /* Four random points uniformly distributed on the surface of an ellipsoid */
  *vtx_coord = malloc(sizeof(double) * 12);
  for (int ivtx = 0; ivtx < 4; ivtx++) {
    double *p = *vtx_coord + 3*ivtx;

    double x = _random_number_normal();
    double y = _random_number_normal();
    double z = _random_number_normal();

    double id = 1./sqrt(x*x + y*y + z*z);

    x *= radii[0] * id;
    y *= radii[1] * id;
    z *= radii[2] * id;

    for (int i = 0; i < 3; i++) {
      p[i] = axes[0][i]*x + axes[1][i]*y + axes[2][i]*z;
    }
  }

  /* Reorient if negative volume */
  if (PDM_predicate_orient3d(*vtx_coord + 3*0,
                             *vtx_coord + 3*1,
                             *vtx_coord + 3*2,
                             *vtx_coord + 3*3) < 0) {
    for (int i = 0; i < 3; i++) {
      double x = (*vtx_coord)[i];
      (*vtx_coord)[  i] = (*vtx_coord)[3+i];
      (*vtx_coord)[3+i] = x;
    }
  }

  /* Faces */
  *face_vtx_idx = PDM_array_new_idx_from_const_stride_int(3, *n_face);
  *face_vtx = malloc(sizeof(int) * 12);
  int *_face_vtx = *face_vtx;
  _face_vtx[ 0] = 2; _face_vtx[ 1] = 3; _face_vtx[ 2] = 4;
  _face_vtx[ 3] = 1; _face_vtx[ 4] = 4; _face_vtx[ 5] = 3;
  _face_vtx[ 6] = 1; _face_vtx[ 7] = 2; _face_vtx[ 8] = 4;
  _face_vtx[ 9] = 1; _face_vtx[10] = 3; _face_vtx[11] = 2;

  *cell_face_idx = PDM_array_new_idx_from_const_stride_int(4, *n_cell);
  *cell_face = malloc(sizeof(int) * 4);
  int *_cell_face = *cell_face;
  _cell_face[0] = 1; _cell_face[1] = 2; _cell_face[2] = 3; _cell_face[3] = 4;
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
 double        *eps,
 int           *n_repeat,
 double        *aniso_max,
 double        *scale,
 int           *visu
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-eps") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *eps = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *seed = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rep") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_repeat = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-aniso") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *aniso_max = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-scale") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *scale = atof(argv[i]);
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


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);

  PDM_predicate_exactinit();

  int    seed      = 0;
  double eps       = 1e-6;
  int    n_repeat  = 1;
  double aniso_max = 1.;
  double scale     = 1.;
  int    visu      = 0;
  _read_args(argc,
             argv,
             &seed,
             &eps,
             &n_repeat,
             &aniso_max,
             &scale,
             &visu);

  if (seed < 0) {
    srand(time(NULL));
  }
  else {
    srand(seed);
  }

  double radii[3];
  radii[0] = scale;
  radii[1] = scale * (1 + _rand01()*aniso_max);
  radii[2] = scale * aniso_max;

  double axes[3][3];
  _random_orthonormal_basis(axes);


  int     n_vtx         = 0;
  double *vtx_coord     = NULL;
  int     n_face        = 0;
  int    *face_vtx_idx  = NULL;
  int    *face_vtx      = NULL;
  int     n_cell        = 0;
  int    *cell_face_idx = NULL;
  int    *cell_face     = NULL;

  _gen_random_tetrahedron(radii,
                          axes,
                          &n_vtx,
                          &vtx_coord,
                          &n_face,
                          &face_vtx_idx,
                          &face_vtx,
                          &n_cell,
                          &cell_face_idx,
                          &cell_face);

  if (visu) {
    int connec[4] = {1, 2, 3, 4};
    PDM_vtk_write_std_elements("tetrahedron.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_TETRA4,
                               n_cell,
                               connec,
                               NULL,
                               0, NULL, NULL);
  }


  free(vtx_coord    );
  free(face_vtx_idx );
  free(face_vtx     );
  free(cell_face_idx);
  free(cell_face    );


  PDM_MPI_Finalize();

  return 0;
}
