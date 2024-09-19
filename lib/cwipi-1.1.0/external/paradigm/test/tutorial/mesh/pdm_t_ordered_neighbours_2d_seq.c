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
#include "pdm_multipart.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_timer.h"

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
     "  -n      <n_vtx_seg> Number of vertices on each side of the cube mesh.\n\n"
     "  -t      <elt_type>  Type of cells.\n\n"
     "  -h                  This message.\n\n");

  exit(exit_code);
}

static void
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *n_vtx_seg)
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
        *n_vtx_seg = (PDM_g_num_t) atol(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/*
 * Modulo operator (C only has a remainder operator) : Copyed from pdm_poly_clipp.c
 */

static int
_modulo
(
 int val,
 int mod
)
{
  if (val >= 0) {
    return val % mod;
  }
  else {
    return val + mod * ((mod - val - 1)/mod);
  }
}

/*
 * Get vtx_vtx_edge: output the table which gives for two
 * vertices the associated edge_idx or 0 if no edge between those vertices
 */

static 
void 
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
 * Get initial edge: to avoid getting stuck in rotation by boundaries of the mesh
 */

static int 
get_initial_edge
(
  int* vtx_edge,
  int* vtx_edge_idx,
  int  i,
  int* edge_face_idx
)
{
  int idx_edge = -1;
  int idx_edge_out = -1;
  for (int j = vtx_edge_idx[i]; j < vtx_edge_idx[i+1]; j++) {
    idx_edge = PDM_ABS(vtx_edge[j])-1;
    if (edge_face_idx[idx_edge+1] - edge_face_idx[idx_edge] == 1) {
      idx_edge_out = idx_edge;
      break;
    }
  }

  if (idx_edge_out == -1) {
    idx_edge_out = idx_edge;
  }

  return idx_edge_out;
}

/*
 * Get next edge: given the current face and edge, get the next neighbour edge in the face
 */

static int 
get_next_edge
(
  int count,
  int  i,
  int  iter_face_idx,
  int  iter_vtx_idx,
  int  n_vtx,
  int* face_vtx_ordered,
  int* face_vtx_ordered_idx,
  int* vtx_vtx_edge,
  int* vtx_ordered_vtx_neighbours,
  int* idx_vtx_ordered_vtx_neighbours
)
{
  int next_edge_idx = -1;
  int other_vtx_idx = -1;
  int j_i = -1;
  int j_iter = -1;
  int direction;

  int length = face_vtx_ordered_idx[iter_face_idx+1] - face_vtx_ordered_idx[iter_face_idx];
  int *fv = &face_vtx_ordered[face_vtx_ordered_idx[iter_face_idx]];

  // Get information for next step
  for (int j = 0; j < length; j++) {
    if (fv[j]-1 == iter_vtx_idx) {
      j_iter = j;
    }
    if (fv[j]-1 == i) {
      j_i = j;
    }
  }


  // Set turn direction
  if (_modulo(j_iter+1, length)== j_i) {
    direction = -1;
  } else {
    direction = 1;
  }

  for (int j = 0; j < length; j++) {
    other_vtx_idx = fv[_modulo((j_iter + direction * j), length)] - 1;

    if ((other_vtx_idx != i) && (other_vtx_idx != iter_vtx_idx)) {
      if (vtx_vtx_edge[i * n_vtx + other_vtx_idx] != 0) {
        next_edge_idx = vtx_vtx_edge[i*n_vtx + other_vtx_idx];
        if (count > 1) {
          vtx_ordered_vtx_neighbours[(*idx_vtx_ordered_vtx_neighbours)++] = other_vtx_idx;
        }
      } else {
        vtx_ordered_vtx_neighbours[(*idx_vtx_ordered_vtx_neighbours)++] = other_vtx_idx;
      }
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
  /*
   *
   * elt_type :
   *  2 -> tria
   *  3 -> quad
   */

  PDM_g_num_t          n_vtx_seg = 10;  // Number of vtx on each side of the cube mesh

  _read_args(argc,
             argv,
             &n_vtx_seg);


  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Set up timer */

  PDM_timer_t *timer = PDM_timer_create();
  PDM_timer_init(timer);
  PDM_timer_resume(timer);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  /* Maillage polysurfacique */

  double       xmin           = 0;
  double       ymin           = 0;
  double       xmax           = xmin + 1;
  double       ymax           = ymin + 1;
  PDM_g_num_t  nx             = n_vtx_seg;
  PDM_g_num_t  ny             = n_vtx_seg;
  int          n_face;
  int          n_vtx;
  int          n_edge;
  PDM_g_num_t  face_gnum;
  PDM_g_num_t  vtx_gnum;
  PDM_g_num_t  edge_gnum;
  int         *face_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx      = NULL;
  double      *vtx_coord      = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  int          n_edge_group;
  int         *edge_group_idx = NULL;
  PDM_g_num_t *edge_group     = NULL;
  int          have_random    = 0;
  int          init_random    = 0;

  PDM_poly_surf_gen(comm,
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    have_random,
                    init_random,
                    nx,
                    ny,
                    &face_gnum,
                    &vtx_gnum,
                    &edge_gnum,
                    &n_vtx,
                    &vtx_coord,
                    &n_face,
                    &face_vtx_idx,
                    &dface_vtx,
                    &dface_edge,
                    &n_edge,
                    &dedge_vtx,
                    &dedge_face,
                    &n_edge_group,
                    &edge_group_idx,
                    &edge_group);

  /*
   *  Goal : vtx->ordered_vtx_neighbours and vtx->ordered_face_neighbours with same order
   *
   */

  /* Cast because n_rank == 1 and get connectivities */

  assert(n_rank == 1);

  int *face_edge = malloc(face_vtx_idx[n_face] * sizeof(int));
  int *face_vtx  = malloc(face_vtx_idx[n_face] * sizeof(int));

  for (int i = 0; i < face_vtx_idx[n_face]; i ++) {
    face_edge[i] = dface_edge[i];
    face_vtx[i]  = dface_vtx[i];
  }

  int *edge_face_idx = NULL;
  int *edge_face = NULL;

  PDM_connectivity_transpose(n_face,
                             n_edge,
                             face_vtx_idx, // face_edge_idx == face_vtx_idx
                             face_edge,
                             &edge_face_idx,
                             &edge_face);


  int *vtx_edge_idx = NULL;
  int *vtx_edge = NULL;

  int *edge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, n_edge);
  int *edge_vtx     = malloc(edge_vtx_idx[n_edge] * sizeof(int));


  for (int i = 0; i < edge_vtx_idx[n_edge]; i ++) {
    edge_vtx[i] = dedge_vtx[i];
  }

  PDM_connectivity_transpose(n_edge,
                             n_vtx,
                             edge_vtx_idx,
                             edge_vtx,
                             &vtx_edge_idx,
                             &vtx_edge);

  int *vtx_vtx_idx = NULL;
  int *vtx_vtx = NULL;

  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_vtx_idx,
                           edge_vtx,
                           &vtx_vtx_idx,
                           &vtx_vtx);

  int *vtx_face_idx = NULL;
  int *vtx_face = NULL;

  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_face_idx,
                           edge_face,
                           &vtx_face_idx,
                           &vtx_face);

  int *vtx_vtx_edge = NULL;

  get_vtx_vtx_edge(&vtx_vtx_edge,
                   edge_vtx,
                   n_edge,
                   n_vtx);

  /* Create ordered vtx neighbour connectivities */

  int iter_edge_idx = -1;
  int iter_face_idx = -1;
  int iter_face_prec_idx;
  int iter_vtx_idx;
  int count;
  int idx_vtx_ordered_face_neighbours = 0;
  int idx_vtx_ordered_vtx_neighbours  = 0;

  int *vtx_ordered_face_neighbours     = malloc(vtx_face_idx[n_vtx] * sizeof(int));
  int *vtx_ordered_face_neighbours_idx = malloc((n_vtx+1) * sizeof(int));
  int *vtx_ordered_vtx_neighbours      = malloc(vtx_vtx_idx[n_vtx]*n_vtx * sizeof(int)); // TO DO mettre taille raisonnable
  int *vtx_ordered_vtx_neighbours_idx  = malloc((n_vtx+1) * sizeof(int));

  vtx_ordered_face_neighbours_idx[0] = idx_vtx_ordered_face_neighbours;
  vtx_ordered_vtx_neighbours_idx[0]  = idx_vtx_ordered_vtx_neighbours;

  /* Begin measure time */

  PDM_timer_hang_on(timer);
  b_t_elapsed = PDM_timer_elapsed(timer);
  b_t_cpu     = PDM_timer_cpu(timer);
  b_t_cpu_u   = PDM_timer_cpu_user(timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(timer);
  PDM_timer_resume(timer);

  for (int i = 0; i < n_vtx; i++) {

    // Initialisation: Choose an edge, a vertex and an adjacent face

    /// Edge
    iter_edge_idx = get_initial_edge(vtx_edge,
                                     vtx_edge_idx,
                                     i,
                                     edge_face_idx);

    /// Face and add it
    iter_face_prec_idx = -1; // valeur impossible de face
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

    while(count > 0) {

      if (iter_face_prec_idx != iter_face_idx) {

        /// Complete face_vtx and get next edge (already ABS - 1)
        iter_edge_idx = get_next_edge(count,
                                      i,
                                      iter_face_idx,
                                      iter_vtx_idx,
                                      n_vtx,
                                      face_vtx,
                                      face_vtx_idx,
                                      vtx_vtx_edge,
                                      vtx_ordered_vtx_neighbours,
                                      &idx_vtx_ordered_vtx_neighbours);

        /// Face and add it
        iter_face_prec_idx = iter_face_idx;
        if (edge_face_idx[iter_edge_idx+1] - edge_face_idx[iter_edge_idx] == 2) { // if not boundary

          if (count > 1) {

            if (PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]])-1 == iter_face_idx) {
              iter_face_idx = PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]+1])-1;
            } else {
              iter_face_idx = PDM_ABS(edge_face[edge_face_idx[iter_edge_idx]])-1;
            }

            vtx_ordered_face_neighbours[idx_vtx_ordered_face_neighbours++] = iter_face_idx;

          }

        }

        /// Vertex and add it
        if (edge_vtx[2*iter_edge_idx]-1 == i) {
          iter_vtx_idx  = edge_vtx[2*iter_edge_idx+1]-1;
        } else {
          iter_vtx_idx  = edge_vtx[2*iter_edge_idx]-1;
        }

      } // end if not same face

      count--;
    } // end while

    vtx_ordered_face_neighbours_idx[i+1] = idx_vtx_ordered_face_neighbours;
    vtx_ordered_vtx_neighbours_idx[i+1]  = idx_vtx_ordered_vtx_neighbours;

  } // end loop on faces

  PDM_timer_hang_on(timer);
  e_t_elapsed = PDM_timer_elapsed(timer);
  e_t_cpu     = PDM_timer_cpu(timer);
  e_t_cpu_u   = PDM_timer_cpu_user(timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(timer);
  PDM_timer_resume(timer);

  double dt_elapsed = e_t_elapsed - b_t_elapsed;
  double dt_cpu     = e_t_cpu     - b_t_cpu;
  double dt_cpu_u   = e_t_cpu_u   - b_t_cpu_u;
  double dt_cpu_s   = e_t_cpu_s   - b_t_cpu_s;

  printf("nb_vtx %d timer %f\n", n_vtx, dt_cpu_u);

  PDM_UNUSED(dt_elapsed);
  PDM_UNUSED(dt_cpu);
  PDM_UNUSED(dt_cpu_s);

  /* Output */
  if(1 == 0) {
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

  if (0 == 1) {
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
  PDM_timer_free(timer                );
  free(face_vtx_idx                   );
  free(dface_vtx                      );
  free(vtx_coord                      );
  free(dface_edge                     );
  free(dedge_vtx                      );
  free(dedge_face                     );
  free(edge_group_idx                 );
  free(edge_group                     );
  free(face_edge                      );
  free(face_vtx                       );
  free(edge_face_idx                  );
  free(edge_face                      );
  free(vtx_edge_idx                   );
  free(vtx_edge                       );
  free(edge_vtx_idx                   );
  free(edge_vtx                       );
  free(vtx_vtx_idx                    );
  free(vtx_vtx                        );
  free(vtx_face_idx                   );
  free(vtx_face                       );
  free(vtx_vtx_edge                   );
  free(vtx_ordered_face_neighbours    );
  free(vtx_ordered_face_neighbours_idx);
  free(vtx_ordered_vtx_neighbours     );
  free(vtx_ordered_vtx_neighbours_idx );

  PDM_MPI_Finalize();
  return 0;
}
