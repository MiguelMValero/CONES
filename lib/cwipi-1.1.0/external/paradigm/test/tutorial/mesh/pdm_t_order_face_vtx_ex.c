
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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_dmesh.h"
#include "pdm_part_connectivity_transform.h"

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

/*
 * Auxiliary function
 */

int where_next_couple(int* tab, int target, int start, int end) {
  int idx;
  for (int j = start; j < end; j+=2) {
    if(tab[j] == target) {
      idx = j;
      break;
    }
  }
  return idx;
}

/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{
  int n_vtx_seg = 5; // Number of vtx on each side of the square mesh
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
   *  Goal : get ordered face->vtx for any element type
   *
   */

  /* Output for each face the couples of edge vertices */

  int edge_idx;
  int count = 0;

  for (int i = 0; i < n_face; i++) {
    log_trace("face %d has vertices", i);
    for (int j = face_edge_idx[i]; j < face_edge_idx[i+1]; j++) {
      edge_idx = PDM_ABS(face_edge[j])-1;
      log_trace(" (%d,%d) ", edge_vtx[2*edge_idx], edge_vtx[2*edge_idx+1]);
      count ++;
    }
    log_trace("\n");
  }

  /* Create face edge->(vtx1, vtx2) table with inverted edges and output it */

  int idx = 0;
  int *face_vtx_long = malloc( 2* count * sizeof(int));
  int *face_vtx_long_idx = malloc( (n_face+1) * sizeof(int));

  for (int i = 0; i < n_face; i++) {
    face_vtx_long_idx[i] = idx;
    log_trace("face %d has vertices", i);
    for (int j = face_edge_idx[i]; j < face_edge_idx[i+1]; j++) {
      // if negative edge index
      edge_idx = PDM_ABS(face_edge[j])-1;

      if (face_edge[j] < 0) {
        // make it positive
        face_edge[j] = edge_idx;
        // invert its vertices
        face_vtx_long[idx++] = edge_vtx[2*edge_idx+1];
        face_vtx_long[idx++] = edge_vtx[2*edge_idx];
        log_trace(" (%d,%d) ", face_vtx_long[idx-2], face_vtx_long[idx-1]);
      } else {
        face_vtx_long[idx++] = edge_vtx[2*edge_idx];
        face_vtx_long[idx++] = edge_vtx[2*edge_idx+1];
        log_trace(" (%d,%d) ", face_vtx_long[idx-2], face_vtx_long[idx-1]);
      }

    } // end loop on face edges
    log_trace("\n");
  } // end loop on faces

  face_vtx_long_idx[n_face] = idx;

  /* Create ordered face_vtx */

  idx = 0;
  int idx_get;
  int length;
  int target;
  int *face_vtx_ordered = malloc(count * sizeof(int));

  for (int i = 0; i < n_face; i++) {

    length = (face_vtx_long_idx[i+1] - face_vtx_long_idx[i] - 2) / 2;

    face_vtx_ordered[idx++] = face_vtx_long[face_vtx_long_idx[i]];
    face_vtx_ordered[idx++] = face_vtx_long[face_vtx_long_idx[i]+1];

    target = face_vtx_long[face_vtx_long_idx[i]+1];

    while(length > 1) {

      printf("target %d\n", target);

      idx_get = where_next_couple(face_vtx_long, target, face_vtx_long_idx[i]+2, face_vtx_long_idx[i+1]);

      face_vtx_ordered[idx++] = face_vtx_long[idx_get+1];

      target = face_vtx_long[idx_get+1];

      length--;
    } // end while

  } // end loop on faces

  /* Output ordered face_vtx */

  for (int i = 0; i < n_face; i++) {
    log_trace("face %d has vertices", i);
    for (int j = face_vtx_long_idx[i]/2; j < face_vtx_long_idx[i+1]/2; j++) {
      log_trace(" %d ", face_vtx_ordered[j]);
    }
    log_trace("\n");
  }


  /* Free memory */
  free(face_edge_idx );
  free(face_edge     );
  free(edge_vtx      );
  free(vtx_coord     );

  PDM_MPI_Finalize();
  return 0;
}
