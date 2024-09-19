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
#include "pdm_vtk.h"
#include "pdm_array.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

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
     "  -n      <level>  Total number of elements .\n\n"
     "  -f      <level>  Frequency of extract .\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

static void
_read_args(int    argc,
           char **argv,
           int   *n_vtx_seg,
           int   *n_steps)
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
        *n_vtx_seg = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_steps") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_steps = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static inline double _rand (void)
{
  return 2 * (double) rand() / (double) RAND_MAX - 1;
}


/*
 *  Just a function to generate a simple, 2D mesh
 *  with slightly randomized coordinates
 *  (only contains triangles)
 */
static void
_generate_mesh
(
const int      n_vtx_seg,
      int     *n_face,
      int     *n_edge,
      int     *n_vtx,
      int    **face_edge_idx,
      int    **face_edge,
      int    **face_vtx,
      int    **edge_vtx,
      double **vtx_coord
 )
{
  #define ij2idx(i, j) (1 + (i) + n_vtx_seg*(j))

  *n_vtx  = n_vtx_seg * n_vtx_seg;
  *n_edge = 2*(n_vtx_seg-1)*n_vtx_seg + (n_vtx_seg-1)*(n_vtx_seg-1);
  *n_face = 2*(n_vtx_seg-1)*(n_vtx_seg-1);


  /* Vertices */
  double step = 1. / (double) (n_vtx_seg - 1);
  double noise = 0.49 * step;

  *vtx_coord = (double *) malloc(sizeof(double) * (*n_vtx) * 3);
  int idx = 0;
  for (int j = 0; j < n_vtx_seg; j++) {
    double y = step * j;
    int fy = (j > 0 && j < n_vtx_seg-1);

    for (int i = 0; i < n_vtx_seg; i++) {
      double x = step * i;
      int fx = (i > 0 && i < n_vtx_seg-1);

      int f = fx*fy;
      (*vtx_coord)[idx++] = x + f*noise * _rand();
      (*vtx_coord)[idx++] = y + f*noise * _rand();
      (*vtx_coord)[idx++] = 0.;
    }
  }

  /* Edges */
  *edge_vtx = (int *) malloc(sizeof(int) * (*n_edge) * 2);
  idx = 0;
  /* Horizontal edges */
  for (int j = 0; j < n_vtx_seg; j++) {
    for (int i = 0; i < n_vtx_seg-1; i++) {
      if (j < n_vtx_seg-1) {
        (*edge_vtx)[idx++] = ij2idx(i,  j);
        (*edge_vtx)[idx++] = ij2idx(i+1,j);
      } else {
        (*edge_vtx)[idx++] = ij2idx(i+1,j);
        (*edge_vtx)[idx++] = ij2idx(i,  j);
      }
    }
  }

  /* Vertical edges */
  for (int i = 0; i < n_vtx_seg; i++) {
    for (int j = 0; j < n_vtx_seg-1; j++) {
      if (i > 0) {
        (*edge_vtx)[idx++] = ij2idx(i,j);
        (*edge_vtx)[idx++] = ij2idx(i,j+1);
      } else {
        (*edge_vtx)[idx++] = ij2idx(i,j+1);
        (*edge_vtx)[idx++] = ij2idx(i,j);
      }
    }
  }

  /* Diagonal edges */
  for (int j = 0; j < n_vtx_seg-1; j++) {
    for (int i = 0; i < n_vtx_seg-1; i++) {
      (*edge_vtx)[idx++] = ij2idx(i,  j+1);
      (*edge_vtx)[idx++] = ij2idx(i+1,j);
    }
  }

  int idx_vertical = n_vtx_seg*(n_vtx_seg - 1);
  int idx_diagonal = 2*n_vtx_seg*(n_vtx_seg - 1);

  /* Faces */
  *face_edge = (int *) malloc(sizeof(int) * (*n_face) * 3);
  *face_vtx  = (int *) malloc(sizeof(int) * (*n_face) * 3);
  int idx2 = 0;
  idx = 0;
  for (int j = 0; j < n_vtx_seg-1; j++) {
    for (int i = 0; i < n_vtx_seg-1; i++) {
      (*face_vtx)[idx++] = ij2idx(i,  j);
      (*face_vtx)[idx++] = ij2idx(i+1,j);
      (*face_vtx)[idx++] = ij2idx(i,  j+1);

      (*face_edge)[idx2++] = 1 + i + (n_vtx_seg-1)*j;
      (*face_edge)[idx2++] = -(1 + idx_diagonal + i + (n_vtx_seg-1)*j);
      (*face_edge)[idx2++] = 1 + idx_vertical + (n_vtx_seg-1)*i + j;
      if (i > 0) {
        (*face_edge)[idx2-1] = -(*face_edge)[idx2-1];
      }


      (*face_vtx)[idx++] = ij2idx(i,  j+1);
      (*face_vtx)[idx++] = ij2idx(i+1,j);
      (*face_vtx)[idx++] = ij2idx(i+1,j+1);

      (*face_edge)[idx2++] = 1 + idx_diagonal + i + (n_vtx_seg-1)*j;
      (*face_edge)[idx2++] = 1 + idx_vertical + (n_vtx_seg-1)*(i+1) + j;
      (*face_edge)[idx2++] = 1 + i + (n_vtx_seg-1)*(j+1);
      if (j < n_vtx_seg-2) {
        (*face_edge)[idx2-1] = -(*face_edge)[idx2-1];
      }
    }
  }


  *face_edge_idx = PDM_array_new_idx_from_const_stride_int(3, *n_face);

  #undef ij2idx
}




/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{
  int n_vtx_seg = 10; // Number of vtx on each side of the square mesh
  int n_steps   = 10; // Number of smoothing steps
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_steps);

  PDM_MPI_Init (&argc, &argv);

  /*
   *  Generate a simple polygonal mesh (all faces are triangles)
   */
  int     n_face        = 0;
  int     n_edge        = 0;
  int     n_vtx         = 0;
  int    *face_edge_idx = NULL;
  int    *face_edge     = NULL;
  int    *face_vtx      = NULL;
  int    *edge_vtx      = NULL;
  double *vtx_coord     = NULL;

  _generate_mesh(n_vtx_seg,
                 &n_face,
                 &n_edge,
                 &n_vtx,
                 &face_edge_idx,
                 &face_edge,
                 &face_vtx,
                 &edge_vtx,
                 &vtx_coord);


  /*
   *  Goal apply iterative Laplacian smmoothing to the vertices
   *  i.e. the coordinates of a vtx at step 'i_step' are the average of the
   *  coordinates of its neighbor vtx (linked by an edge) at step 'i_step-1'
   *  (https://en.wikipedia.org/wiki/Laplacian_smoothing)
   *
   *  You get a bonus chocolate if you manage to freeze the boundary vertices ;)
   *
   */

  /* Identify boundary edges */
  int *edge_face_n = PDM_array_zeros_int(n_edge);
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    int iedge = PDM_ABS(face_edge[i]) - 1;
    edge_face_n[iedge]++;
  }


  /* Identify boundary vertices */
  int *is_boundary_vtx = PDM_array_zeros_int(n_vtx);
  for (int iedge = 0; iedge < n_edge; iedge++) {
    if (edge_face_n[iedge] < 2) {
      for (int i = 2*iedge; i < 2*(iedge+1); i++) {
        int ivtx = edge_vtx[i] - 1;
        is_boundary_vtx[ivtx] = 1;
      }
    }
  }
  free(edge_face_n);

  /* Compute normalization factor for each vtx (constant over time) */
  double *normalization = malloc(sizeof(double) * n_vtx);

  for (int i = 0; i < n_vtx; i++) {
    normalization[i] = 0;
  }
  for (int i = 0; i < 2*n_edge; i++) {
    int ivtx = edge_vtx[i] - 1;
    normalization[ivtx] += 1;
  }
  for (int i = 0; i < n_vtx; i++) {
    assert(normalization[i] > 0.);
    normalization[i] = 1. / normalization[i];
  }



  double *new_vtx_coord = malloc(sizeof(double) * n_vtx * 3);

  for (int i_step = 0; i_step <= n_steps; i_step++) {

    /*
     *  Export current mesh to vtk format to get a nice video :D
     */
    if(1 == 0) {
      char filename[999];
      sprintf(filename, "laplacian_smoothing_seq_%3.3d.vtk", i_step);

      PDM_vtk_write_polydata(filename,
                             n_vtx,
                             vtx_coord,
                             NULL,
                             n_face,
                             face_edge_idx,
                             face_vtx,
                             NULL,
                             NULL);
    }

    /*
     *  Apply one smoothing step
     */
    for (int i = 0; i < 3*n_vtx; i++) {
      new_vtx_coord[i] = 0;
    }

    for (int iedge = 0; iedge < n_edge; iedge++) {
      int *ev = edge_vtx + 2*iedge;

      for (int i = 0; i < 2; i++) {
        int ivtx1 = ev[i]       - 1;
        int ivtx2 = ev[(i+1)%2] - 1;

        for (int j = 0; j < 3; j++) {
          new_vtx_coord[3*ivtx1+j] += vtx_coord[3*ivtx2+j];
        }
      }
    }

    /* Update coordinates of interior vtx */
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      if (!is_boundary_vtx[ivtx]) {
        for (int i = 0; i < 3; i++) {
          vtx_coord[3*ivtx+i] = normalization[ivtx] * new_vtx_coord[3*ivtx+i];
        }
      }
    }

    if (i_step == n_steps-1) {
      break;
    }
  }
  free(normalization);
  free(new_vtx_coord);
  free(is_boundary_vtx);

  free(face_edge_idx);
  free(face_edge    );
  free(face_vtx     );
  free(edge_vtx     );
  free(vtx_coord    );

  PDM_MPI_Finalize();
  return 0;
}
