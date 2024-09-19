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

static int where_next_couple(int* tab, int target, int start, int end) {
  int idx = -1;
  for (int j = start; j < end; j+=2) {
    if(tab[j] == target) {
      idx = j;
      break;
    }
  }
  return idx;
}

/*
 * Get ordered face->vtx
 */

static
void get_ordered_face_vtx
(
  int** face_vtx_ordered,
  int** face_vtx_ordered_idx,
  int*  face_edge_idx,
  int*  face_edge,
  int*  edge_vtx,
  int   n_face
)
{
  /* Create face edge->(vtx1, vtx2) table with inverted edges and output it */

  int edge_idx;
  int count = 0;

  for (int i = 0; i < n_face; i++) {
    for (int j = face_edge_idx[i]; j < face_edge_idx[i+1]; j++) {
      edge_idx = PDM_ABS(face_edge[j])-1;
      count ++;
    }
  }

  int idx = 0;
  int *face_vtx_long     = malloc( 2* count * sizeof(int));
  int *face_vtx_long_idx = malloc( (n_face+1) * sizeof(int));
  *face_vtx_ordered_idx = (int * ) malloc( (n_face+1) * sizeof(int));

  for (int i = 0; i < n_face; i++) {
    face_vtx_long_idx[i] = idx;
    (*face_vtx_ordered_idx)[i] = idx/2;
    for (int j = face_edge_idx[i]; j < face_edge_idx[i+1]; j++) {
      // if negative edge index
      edge_idx = PDM_ABS(face_edge[j])-1;

      if (face_edge[j] < 0) {
        // make it positive
        face_edge[j] = edge_idx;
        // invert its vertices
        face_vtx_long[idx++] = edge_vtx[2*edge_idx+1];
        face_vtx_long[idx++] = edge_vtx[2*edge_idx];
      } else {
        face_vtx_long[idx++] = edge_vtx[2*edge_idx];
        face_vtx_long[idx++] = edge_vtx[2*edge_idx+1];
      }

    } // end loop on face edges
  } // end loop on faces

  face_vtx_long_idx[n_face] = idx;
  (*face_vtx_ordered_idx)[n_face] = idx/2;

  /* Create ordered face_vtx */

  idx = 0;
  int idx_get;
  int length;
  int target;
  *face_vtx_ordered     = (int * ) malloc(count * sizeof(int));

  for (int i = 0; i < n_face; i++) {

    length = (face_vtx_long_idx[i+1] - face_vtx_long_idx[i] - 2) / 2;

    (*face_vtx_ordered)[idx++] = face_vtx_long[face_vtx_long_idx[i]];
    (*face_vtx_ordered)[idx++] = face_vtx_long[face_vtx_long_idx[i]+1];

    target = face_vtx_long[face_vtx_long_idx[i]+1];

    while(length > 1) {

      idx_get = where_next_couple(face_vtx_long, target, face_vtx_long_idx[i]+2, face_vtx_long_idx[i+1]);

      (*face_vtx_ordered)[idx++] = face_vtx_long[idx_get+1];

      target = face_vtx_long[idx_get+1];

      length--;
    } // end while

  } // end loop on faces

  // PDM_log_trace_array_int((*face_vtx_ordered_idx), n_face, "(*face_vtx_ordered_idx): ");
  // PDM_log_trace_array_int((*face_vtx_ordered), count, "(*face_vtx_ordered): ");
  free(face_vtx_long     );
  free(face_vtx_long_idx );

}

/*
 * Get vtx_vtx_idx
 */

static void 
get_vtx_vtx_edge
(
  int** vtx_vtx_edge,
  int*   edge_vtx,
  int n_edge,
  int n_vtx
)
{
  int idx_vtx1;
  int idx_vtx2;
  *vtx_vtx_edge = (int * ) PDM_array_const_int (n_vtx * n_vtx, 0);

  for (int i = 0; i < n_edge; i++) {
    idx_vtx1 = edge_vtx[2*i]-1;
    idx_vtx2 = edge_vtx[2*i+1]-1;

    (*vtx_vtx_edge)[idx_vtx1 * n_vtx + idx_vtx2] = i;
    (*vtx_vtx_edge)[idx_vtx2 * n_vtx + idx_vtx1] = i;

  } // end loop on edges

}

/*
 * Get next edge
 */

static int 
get_next_edge
(
  int  i,
  int  iter_face_idx,
  int  iter_vtx_idx,
  int  n_vtx,
  int* face_vtx_ordered,
  int* face_vtx_ordered_idx,
  int* vtx_vtx_edge
)
{
  int next_edge_idx = -1;
  int other_vtx_idx;

  for (int j = face_vtx_ordered_idx[iter_face_idx]; j < face_vtx_ordered_idx[iter_face_idx+1]; j++) {
    other_vtx_idx = face_vtx_ordered[j]-1;
    if ((vtx_vtx_edge[i * n_vtx + other_vtx_idx] != 0) && (other_vtx_idx != iter_vtx_idx)) {
      next_edge_idx = vtx_vtx_edge[i*n_vtx + other_vtx_idx];
      break;
    }
  }
  return next_edge_idx;
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
   *  Goal : vtx->ordered_vtx_neighbours and vtx->ordered_face_neighbours with same order
   *
   */

  /* Get connectivities */

  int *vtx_edge_idx = NULL;
  int *vtx_edge = NULL;

  int *edge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, n_edge);

  PDM_connectivity_transpose(n_edge,
                             n_vtx,
                             edge_vtx_idx,
                             edge_vtx,
                             &vtx_edge_idx,
                             &vtx_edge);

  int *edge_face_idx = NULL;
  int *edge_face = NULL;

  PDM_connectivity_transpose(n_face,
                             n_edge,
                             face_edge_idx,
                             face_edge,
                             &edge_face_idx,
                             &edge_face);

  int *vtx_vtx_idx = NULL;
  int *vtx_vtx = NULL;

  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_vtx_idx,
                           edge_vtx,
                           &vtx_vtx_idx,
                           &vtx_vtx);
  free(edge_vtx_idx);
  int *vtx_face_idx = NULL;
  int *vtx_face = NULL;

  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_face_idx,
                           edge_face,
                           &vtx_face_idx,
                           &vtx_face);
  int *face_vtx_ordered     = NULL;
  int *face_vtx_ordered_idx = NULL;

  get_ordered_face_vtx(&face_vtx_ordered,
                       &face_vtx_ordered_idx,
                       face_edge_idx,
                       face_edge,
                       edge_vtx,
                       n_face);

  int *vtx_vtx_edge = NULL;

  get_vtx_vtx_edge(&vtx_vtx_edge,
                   edge_vtx,
                   n_edge,
                   n_vtx);

  /* Create ordered vtx neighbour connectivities */

  int iter_edge_idx;
  int iter_face_idx;
  int iter_vtx_idx;
  int count;
  int idx_vtx_ordered_face_neighbours = 0;
  int idx_vtx_ordered_vtx_neighbours  = 0;

  int *vtx_ordered_face_neighbours     = malloc(vtx_face_idx[n_vtx] * sizeof(int));
  int *vtx_ordered_face_neighbours_idx = malloc((n_vtx+1) * sizeof(int));
  int *vtx_ordered_vtx_neighbours      = malloc(vtx_vtx_idx[n_vtx] * sizeof(int));
  int *vtx_ordered_vtx_neighbours_idx  = malloc((n_vtx+1) * sizeof(int));

  vtx_ordered_face_neighbours_idx[0] = idx_vtx_ordered_face_neighbours;
  vtx_ordered_vtx_neighbours_idx[0]  = idx_vtx_ordered_vtx_neighbours;

  for (int i = 0; i < n_vtx; i++) {

    // Initialisation: Choose an edge, a vertex and an adjacent face

    /// Edge
    iter_edge_idx = PDM_ABS(vtx_edge[vtx_edge_idx[i]])-1;

    /// Face and add it
    iter_face_idx = PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]])-1;
    vtx_ordered_face_neighbours[idx_vtx_ordered_face_neighbours++] = iter_face_idx;

    /// Vertex and add it
    if (edge_vtx[2*iter_edge_idx]-1 == i) {
      iter_vtx_idx  = edge_vtx[2*iter_edge_idx+1]-1;
    } else {
      iter_vtx_idx  = edge_vtx[2*iter_edge_idx]-1;
    }
    vtx_ordered_vtx_neighbours[idx_vtx_ordered_vtx_neighbours++] = iter_vtx_idx;

    /// Counter
    count = vtx_vtx_idx[i+1] - vtx_vtx_idx[i]-1; // not to count me

    // Iteration: Loop over the face and get next edge, vertex and face

    while(count > 1) {

      /// Complete face_vtx and get next edge (already ABS - 1)
      iter_edge_idx = get_next_edge(i,
                                    iter_face_idx,
                                    iter_vtx_idx,
                                    n_vtx,
                                    face_vtx_ordered,
                                    face_vtx_ordered_idx,
                                    vtx_vtx_edge);

      /// Face and add it
      if (edge_face_idx[iter_edge_idx+1] - edge_face_idx[iter_edge_idx] == 2) { // if not boundary
        if (PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]])-1 == iter_face_idx) {
          iter_face_idx = PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]+1])-1;
        } else {
          iter_face_idx = PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]])-1;
        }
        vtx_ordered_face_neighbours[idx_vtx_ordered_face_neighbours++] = iter_face_idx;

      }

      /// Vertex and add it
      if (edge_vtx[2*iter_edge_idx]-1 == i) {
        iter_vtx_idx  = edge_vtx[2*iter_edge_idx+1]-1;
      } else {
        iter_vtx_idx  = edge_vtx[2*iter_edge_idx]-1;
      }
      vtx_ordered_vtx_neighbours[idx_vtx_ordered_vtx_neighbours++] = iter_vtx_idx;

      count--;
    } // end while

    vtx_ordered_face_neighbours_idx[i+1] = idx_vtx_ordered_face_neighbours;
    vtx_ordered_vtx_neighbours_idx[i+1]  = idx_vtx_ordered_vtx_neighbours;

  } // end loop on faces

  /* Output */
  if(0 == 1) {
    for (int i = 0; i < n_vtx; i++) {
      log_trace("Sommet: %d, Vertices: ", i);
      for (int j = vtx_ordered_vtx_neighbours_idx[i]; j < vtx_ordered_vtx_neighbours_idx[i+1]; j++) {
        log_trace(" %d ", vtx_ordered_vtx_neighbours[j]);
      }

      log_trace(", Faces: ");
      for (int j = vtx_ordered_face_neighbours_idx[i]; j < vtx_ordered_face_neighbours_idx[i+1]; j++) {
        log_trace(" %d ", vtx_ordered_face_neighbours[j]);
      }
      log_trace("\n");
    }
  }

  /* Output mesh in vtk format TO DO change if not triangles anymore */
  // Create face_vtx_idx
  int face_vtx_idx[n_face+1];
  for (int i = 0; i < n_face+1; i++) {
    face_vtx_idx[i] = 3*i;
  }

  if(1 == 0) {
    char filename[999];
    sprintf(filename, "mesh.vtk");
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

  /* Free memory */
  free(face_edge_idx );
  free(face_edge     );
  free(edge_vtx      );
  free(vtx_coord     );
  free(vtx_edge_idx  );
  free(vtx_edge      );
  free(edge_face_idx );
  free(edge_face     );
  free(vtx_vtx_idx   );
  free(vtx_vtx       );
  free(vtx_vtx_edge);
  free(vtx_face_idx);
  free(vtx_face);
  free(face_vtx);

  free(vtx_ordered_face_neighbours    );
  free(vtx_ordered_face_neighbours_idx);
  free(vtx_ordered_vtx_neighbours     );
  free(vtx_ordered_vtx_neighbours_idx );

  free(face_vtx_ordered);
  free(face_vtx_ordered_idx);

  PDM_MPI_Finalize();
  return 0;
}
