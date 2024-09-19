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

#include "pdm_part_connectivity_transform.h" // Added


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

// Is in table
static int 
is_in_table
(
  int *table, 
  int value,
  int size
) 
{
  int boolean = 0;
  for (int i = 0; i < size; i++) {
    if (table[i] == value) {
      boolean = 1;
    }
  }
  return boolean;
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
   *  You get a bonus chocolate if you manage to freeze the boundary vertices ;) // ici bord si N == 2 ou 4
   *
   */

  // New coordinate data structure
  double *vtx_new_coord = malloc(n_vtx*3 * sizeof(double));

  // Create edge_vtx_idx
  int *edge_vtx_idx = malloc((n_edge+1) * sizeof(int)); // cf pdm_array
  for (int i = 0; i < n_edge+1; i++) {
    edge_vtx_idx[i] = 2*i;
  }

  // Compute vtx_edge connectivity
  int *vtx_edge_idx = NULL;
  int *vtx_edge = NULL;

  PDM_connectivity_transpose(n_edge,
                             n_vtx,
                             edge_vtx_idx,
                             edge_vtx,
                            &vtx_edge_idx,
                            &vtx_edge);

  // Compute edge_face connectivity
  int *edge_face_idx = NULL;
  int *edge_face = NULL;

  PDM_connectivity_transpose(n_face,
                             n_edge,
                             face_edge_idx,
                             face_edge,
                            &edge_face_idx,
                            &edge_face);
  // Entities for vtk output
  char filename[999];

  // Entities for laplace smoothing
  int idx = 0;
  int N = 1; // to avoid 0 zero division
  int idx_vtx_of_edge = 0;
  int idx_vtx_other = 0;
  int idx_edge = 0;
  int n_face_edge = 0;
  double *laplace_smoothing_vtx_coord = malloc(2 * sizeof(double));

  // Create face_vtx_idx
  int face_vtx_idx[n_face+1];
  for (int i = 0; i < n_face+1; i++) {
    face_vtx_idx[i] = 3*i;
  }

  // Create group
  int *bdr_group = malloc(10 * n_vtx * sizeof(double)); // il y a plus élégant...
  int size = 0;

  for (int i = 0; i < n_vtx; i++) {
    N = (vtx_edge_idx[i+1]- vtx_edge_idx[i]);
    for (int j = 0; j < N; j++){
      idx_edge = vtx_edge[idx]-1;
      n_face_edge = (edge_face_idx[idx_edge+1]- edge_face_idx[idx_edge]);
      if (n_face_edge < 2) {
        idx_vtx_of_edge = edge_vtx[edge_vtx_idx[vtx_edge[idx]-1]]-1;
        idx_vtx_other = edge_vtx[edge_vtx_idx[vtx_edge[idx]-1]+1]-1;
        bdr_group[size++] = idx_vtx_of_edge;
        bdr_group[size++] = idx_vtx_other;
      }
      idx++;
    }
  }

  // Do smoothing stepping
  for (int i_step = 0; i_step <= n_steps; i_step++) {

    /*
     *  Export current mesh to vtk format to get a nice video :D
     */
    if(1 == 0) {
      sprintf(filename, "out_mesh_%2.2d.vtk", i_step); // TO DO %2.2d change if more steps
      PDM_vtk_write_polydata(filename,
                             n_vtx,
                             vtx_coord,
                             NULL,
                             n_face,
                             face_vtx_idx,
                             face_vtx,
                             NULL,
                             NULL);
    }

    /*
     *  Apply one smoothing step
     */

    // 2D version
    idx = 0;
    laplace_smoothing_vtx_coord[0] = 0;
    laplace_smoothing_vtx_coord[1] = 0;

    // For each vertex
    for (int i = 0; i < n_vtx; i++) {
      // For each edge coming out of that vertex
      N = (vtx_edge_idx[i+1]- vtx_edge_idx[i]);
      for (int j = 0; j < N; j++){

        if (!(is_in_table(bdr_group, i, size))) {

        idx_vtx_of_edge = edge_vtx[edge_vtx_idx[vtx_edge[idx]-1]]-1;
        idx_vtx_other = edge_vtx[edge_vtx_idx[vtx_edge[idx]-1]+1]-1;

          if (i == idx_vtx_of_edge) {
            laplace_smoothing_vtx_coord[0] += vtx_coord[3*idx_vtx_other];
            laplace_smoothing_vtx_coord[1] += vtx_coord[3*idx_vtx_other+1];
          } else {
            laplace_smoothing_vtx_coord[0] += vtx_coord[3*idx_vtx_of_edge];
            laplace_smoothing_vtx_coord[1] += vtx_coord[3*idx_vtx_of_edge+1];
          }
        }

        idx++;

      }
      vtx_new_coord[3*i] = laplace_smoothing_vtx_coord[0] / N;
      vtx_new_coord[3*i+1] = laplace_smoothing_vtx_coord[1] / N;
      laplace_smoothing_vtx_coord[0] = 0;
      laplace_smoothing_vtx_coord[1] = 0;
    }

    // put back in old table
    for (int i = 0; i < n_vtx; i++) {
      if (!(is_in_table(bdr_group, i, size))) {
          vtx_coord[3*i] = vtx_new_coord[3*i];
          vtx_coord[3*i+1] = vtx_new_coord[3*i+1];
      }
    }

    if (i_step == n_steps-1) {
      break;
    }
  }

  // Additional free
  free(bdr_group);
  free(vtx_new_coord);
  free(vtx_edge_idx);
  free(vtx_edge);
  free(laplace_smoothing_vtx_coord);

  free(face_edge_idx);
  free(face_edge    );
  free(face_vtx     );
  free(edge_vtx     );
  free(vtx_coord    );

  free(edge_face_idx);
  free(edge_face);
  free(edge_vtx_idx);
  PDM_MPI_Finalize();
  return 0;
}
