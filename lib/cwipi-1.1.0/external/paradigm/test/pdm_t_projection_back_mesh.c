#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_ho_basis.h"
#include "pdm_ho_ordering.h"
#include "pdm_array.h"
#include "pdm_triangle.h"
#include "pdm_line.h"
#include "pdm_ho_bezier.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_ho_location.h"
#include "pdm_sort.h"

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
_read_args
(
 int      argc,
 char   **argv,
 char   **filename,
 double  *pt_coord,
 int     *random_seed,
 int     *visu
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-p") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc)
          _usage(EXIT_FAILURE);
        else {
          pt_coord[j] = atof(argv[i]);
        }
      }
    }

    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *random_seed = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_read_back_mesh
(
 const char            *filename,
 PDM_Mesh_nodal_elt_t  *elt_type,
 int                   *elt_order,
 int                   *n_vtx,
 double               **vtx_coord,
 int                   *n_elt,
 int                  **elt_vtx
 )
{
  FILE *f = fopen(filename, "r");

  assert(f != NULL);

  int elt_dim = 0;

  fscanf(f, "%d %d %d %d\n", n_vtx, n_elt, &elt_dim, elt_order);

  if (elt_dim == 1) {
    if (*elt_order == 1) {
      *elt_type = PDM_MESH_NODAL_BAR2;
    }
    else {
      *elt_type = PDM_MESH_NODAL_BARHO_BEZIER;
    }
  }
  else if (elt_dim == 2) {
    if (*elt_order == 1) {
      *elt_type = PDM_MESH_NODAL_TRIA3;
    }
    else {
      *elt_type = PDM_MESH_NODAL_TRIAHO_BEZIER;
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Wrong elt_dim %d\n", elt_dim);
  }


  *vtx_coord = malloc(sizeof(double) * (*n_vtx) * 3);
  for (int i = 0; i < *n_vtx; i++) {
    fscanf(f, "%lf %lf %lf\n",
           *vtx_coord + 3*i, *vtx_coord + 3*i+1, *vtx_coord + 3*i+2);
  }

  int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(*elt_type,
                                               *elt_order);

  *elt_vtx = malloc(sizeof(int) * elt_vtx_n * (*n_elt));
  for (int i = 0; i < *n_elt; i++) {
    for (int j = 0; j < elt_vtx_n; j++) {
      fscanf(f, "%d", *elt_vtx + elt_vtx_n*i + j);
    }
  }

  fclose(f);
}


/**
 * \brief Build the edge-edge connectivity of a manifold curve mesh
 * (i.e. a vertex is connected to at most two edges)
 */
static void
_build_edge_edge_connectivity
(
 const int                   elt_order,
 const int                   n_vtx,
 const int                   n_edge,
       int                  *edge_vtx,
       int                  *parent_node,
       int                  *edge_edge
 )
{
  int stride_ho = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, elt_order);
  int stride    = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, 1);

  PDM_array_reset_int(edge_edge, stride*n_edge, 0);

  int *vtx_edge = PDM_array_zeros_int(n_vtx * 2);

  for (int iedge = 0; iedge < n_edge; iedge++) {
    int *ev = edge_vtx + stride_ho*iedge;

    for (int i = 0; i < stride; i++) {
      int id_vtx;
      if (parent_node == NULL) {
        id_vtx = ev[i] - 1;
      }
      else {
        id_vtx = ev[parent_node[i]] - 1;
      }

      int jedge = vtx_edge[2*id_vtx] - 1;
      if (jedge < 0) {
        vtx_edge[2*id_vtx  ] = iedge + 1;
        vtx_edge[2*id_vtx+1] = i;
      }
      else {
        int j = vtx_edge[2*id_vtx+1];

        edge_edge[stride*iedge + i] = jedge + 1;
        if (edge_edge[stride*jedge + j] != 0) {
          PDM_error(__FILE__, __LINE__, 0,
                    "Vertex %d shared by more than two edges\n",
                    id_vtx+1);
        }
        edge_edge[stride*jedge + j] = iedge + 1;
      }
    }
  }

  free(vtx_edge);
}


/**
 * \brief Build the face-face connectivity of a manifold surface mesh
 * (i.e. an edge is connected to at most two faces)
 */
static void
_build_face_face_connectivity
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   elt_order,
 const int                   n_vtx,
 const int                   n_face,
       int                  *face_vtx,
       int                  *parent_node,
       int                  *face_face
 )
{
  int stride_ho = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, elt_order);
  int stride    = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);

  PDM_array_reset_int(face_face, stride*n_face, 0);

  int key_max = PDM_MAX(0, 2*(n_vtx - 1));
  int n_key = key_max + 1;
  int *key_edge_n = PDM_array_zeros_int(n_key);

  int fvi[2];

  for (int iface = 0; iface < n_face; iface++) {
    int *fv = face_vtx + stride_ho * iface;
    for (int i = 0; i < stride; i++) {
      if (parent_node == NULL) {
        fvi[0] = fv[i] - 1;
        fvi[1] = fv[(i+1)%stride_ho] - 1;
      }
      else {
        fvi[0] = fv[parent_node[i]] - 1;
        fvi[1] = fv[parent_node[(i+1)%stride]] - 1;;
      }
      int key = fvi[0] + fvi[1];
      key_edge_n[key]++;
    } // End of loop on current face's edges
  } // End of loop on faces

  int *key_edge_idx = PDM_array_new_idx_from_sizes_int(key_edge_n, n_key);
  PDM_array_reset_int(key_edge_n, n_key, 0);

  int *key_edge = (int *) malloc(sizeof(int) * key_edge_idx[n_key] * 2);


  for (int iface = 0; iface < n_face; iface++) {
    int *fv = face_vtx + stride_ho * iface;
    for (int i = 0; i < stride; i++) {
      if (parent_node == NULL) {
        fvi[0] = fv[i] - 1;
        fvi[1] = fv[(i+1)%stride] - 1;
      }
      else {
        fvi[0] = fv[parent_node[i]] - 1;
        fvi[1] = fv[parent_node[(i+1)%stride]] - 1;;
      }
      int key = fvi[0] + fvi[1];

      int *ke = key_edge + 2*key_edge_idx[key];

      int found = 0;
      for (int j = 0; j < key_edge_n[key]; j++) {
        int jface = ke[2*j  ];
        int jedge = ke[2*j+1];
        int jvtx;
        if (parent_node == NULL) {
          jvtx = face_vtx[stride_ho*jface + jedge] - 1;
        }
        else {
          jvtx = face_vtx[stride_ho*jface + parent_node[jedge]] - 1;
        }

        if (jvtx == fvi[0] || jvtx == fvi[1]) {
          found = 1;
          face_face[stride*iface + i    ] = jface + 1;
          if (face_face[stride*jface + jedge] != 0) {
            PDM_error(__FILE__, __LINE__, 0,
                      "Edge %d %d shared by more than two faces\n",
                      fvi[0]+1, fvi[1]+1);
          }
          face_face[stride*jface + jedge] = iface + 1;
          break;
        }
      }

      if (!found) {
        ke[2*key_edge_n[key]  ] = iface;
        ke[2*key_edge_n[key]+1] = i;
        key_edge_n[key]++;
      }
    } // End of loop on current face's edges
  } // End of loop on faces
  free(key_edge_n);
  free(key_edge_idx);
  free(key_edge);

  if (0) {
    int *idx = PDM_array_new_idx_from_const_stride_int(stride, n_face);
    PDM_log_trace_connectivity_int(idx, face_face, n_face, "face_face : ");
    free(idx);
  }
}




static int
_projection_on_background_mesh_get2
(
 double               *pt_to_project_coord,
 int                   n_back_elt,
 int                   start_elt,
 // int                  *cavity_back_elt,
 int                  *back_elt_vtx_idx,
 int                  *back_elt_vtx,
 int                  *back_elt_elt,
 double               *back_vtx_coord,
 int                   back_elt_order,
 PDM_Mesh_nodal_elt_t  back_elt_type,
 double               *proj_pt_coord,
 int                  *closest_back_elt,
 int                  *history_elt,
 double               *history_proj
)
{
  const int vb = 0;

  int stride = PDM_Mesh_nodal_n_vtx_elt_get(back_elt_type, 1);

  int *is_visited = PDM_array_zeros_int(n_back_elt);
  int *stack      = malloc(sizeof(int) * n_back_elt);

  int pos_stack = 0;
  stack[pos_stack++] = start_elt;
  is_visited[start_elt] = 1;

  double min_distance = HUGE_VAL;
  int    closest_elt  = 0;

  int iter = -1;


  while (pos_stack > 0) {
    iter++;

    int id_elt = stack[--pos_stack];

    if (history_elt != NULL) {
      history_elt[iter] = id_elt + 1;
    }

    int *ev = back_elt_vtx + back_elt_vtx_idx[id_elt];
    int *ee = back_elt_elt + stride*id_elt;

    double cp[3], distance = 0, uvw[3];

    if (back_elt_type == PDM_MESH_NODAL_BAR2) {

      distance = PDM_line_distance(pt_to_project_coord,
                                   back_vtx_coord + 3*(ev[0] - 1),
                                   back_vtx_coord + 3*(ev[1] - 1),
                                   uvw,
                                   cp);
      uvw[1] = 1 - uvw[0];
      uvw[2] = 0.5;//

    } // End if PDM_MESH_NODAL_BAR2


    else if (back_elt_type == PDM_MESH_NODAL_BARHO) {

      const int n_node = back_elt_order+1;
      double bar_coord[n_node*3];
      for (int i = 0; i < n_node; i++) {
        int id_vtx = ev[i] - 1;
        memcpy(bar_coord + 3*i,
               back_vtx_coord + 3*id_vtx,
               sizeof(double) * 3);
      }

      distance = PDM_ho_location(back_elt_type,
                                 back_elt_order,
                                 n_node,
                                 bar_coord,
                                 pt_to_project_coord,
                                 cp,
                                 uvw);
      uvw[1] = 1 - uvw[0];
      uvw[2] = 0.5;//

    } // End if PDM_MESH_NODAL_BARHO


    else if (back_elt_type == PDM_MESH_NODAL_TRIA3) {

      double tria_coord[9];
      for (int i = 0; i < 3; i++) {
        int id_vtx = ev[i] - 1;
        memcpy(tria_coord + 3*i,
               back_vtx_coord + 3*id_vtx,
               sizeof(double) * 3);
      }

      if (0) {
        // PDM_triangle_evaluate_position_old(pt_to_project_coord,
        //                                    tria_coord,
        //                                    cp,
        //                                    &distance,
        //                                    uvw);
        double weight[3];
        PDM_triangle_evaluate_position(pt_to_project_coord,
                                       tria_coord,
                                       cp,
                                       &distance,
                                       uvw);
        uvw[0] = weight[2];
        uvw[1] = weight[0];
        uvw[2] = weight[1];
        // PERMUTATION??
      } else {
        double weight[3];
        PDM_triangle_closest_point(pt_to_project_coord,
                                   tria_coord,
                                   cp,
                                   &distance,
                                   weight);
        uvw[0] = weight[2];
        uvw[1] = weight[0];
        uvw[2] = weight[1];
      }
    } // End if PDM_MESH_NODAL_TRIA3


    else if (back_elt_type == PDM_MESH_NODAL_TRIAHO) {

      const int n_node = (back_elt_order+1)*(back_elt_order+2)/2;
      double tria_coord[n_node*3];
      for (int i = 0; i < n_node; i++) {
        int id_vtx = ev[i] - 1;
        memcpy(tria_coord + 3*i,
               back_vtx_coord + 3*id_vtx,
               sizeof(double) * 3);
      }

      distance = PDM_ho_location(back_elt_type,
                                 back_elt_order,
                                 n_node,
                                 tria_coord,
                                 pt_to_project_coord,
                                 cp,
                                 uvw);
      double u = uvw[0];
      double v = uvw[1];
      double w = uvw[2];
      uvw[0] = v;
      uvw[1] = w;
      uvw[2] = u;

    } // End if PDM_MESH_NODAL_TRIAHO


    else if (back_elt_type == PDM_MESH_NODAL_TRIAHO_BEZIER) {

      const int n_node = (back_elt_order+1)*(back_elt_order+2)/2;
      double tria_coord[n_node*3];
      for (int i = 0; i < n_node; i++) {
        int id_vtx = ev[i] - 1;
        memcpy(tria_coord + 3*i,
               back_vtx_coord + 3*id_vtx,
               sizeof(double) * 3);
      }

      distance = PDM_ho_bezier_triangle_location(back_elt_order,
                                                 n_node,
                                                 tria_coord,
                                                 pt_to_project_coord,
                                                 cp,
                                                 uvw);
      double u = uvw[0];
      double v = uvw[1];
      double w = uvw[2];
      uvw[0] = v;
      uvw[1] = w;
      uvw[2] = u;

    } // End if PDM_MESH_NODAL_TRIAHO_BEZIER


    else {
      PDM_error(__FILE__, __LINE__, 0, "Projection not yet implemented for type %d\n", (int) back_elt_type);
    }

    if (vb) {
      log_trace("iter %d, id_elt %d, dist = %f",iter, id_elt, distance);
      PDM_log_trace_array_double(uvw, stride-1, ", uvw = ");
      log_trace("cp = %f %f %f\n", cp[0], cp[1], cp[2]);
    }

    if (distance < min_distance) {
      memcpy(proj_pt_coord, cp,  sizeof(double) * 3);
      closest_elt  = id_elt + 1;
      min_distance = distance;
    }

    if (history_proj != NULL) {
      memcpy(history_proj + 3*iter, proj_pt_coord, sizeof(double) * 3);
    }


    int    n_out = 0;
    int    i_out[3];
    double w_out[3];
    for (int i = 0; i < stride; i++) {
      if (uvw[i] < 1e-6) {//<= 0) {
        int id_ngb = ee[i] - 1;
        if (vb) {
          log_trace("  i = %d, uvw[i] = %f, id_ngb = %d\n", i, uvw[i], id_ngb);
        }

        if (id_ngb >= 0) {
          if (!is_visited[id_ngb]) {
            w_out[n_out] = -uvw[i];
            i_out[n_out] = id_ngb;
            n_out++;
          }
        }
      }
    } // End of loop on neighbor elements

    if (vb) {
      PDM_log_trace_array_int(i_out, n_out, "  i_out : ");
    }

    if (n_out == 0) {
      // We found a local minimum of distance :)
      break;
    }
    else {
      PDM_sort_double(w_out, i_out, n_out);
      for (int i = 0; i < n_out; i++) {
        int id_ngb = i_out[i];
        stack[pos_stack++] = id_ngb;
        is_visited[id_ngb] = 1;
      }
    }

  } // End of while loop
  free(is_visited);
  free(stack);

  *closest_back_elt = closest_elt;

  if (vb) {
    printf("closest_elt = %d (0-based)\n", closest_elt-1);
    log_trace("proj_pt_coord = %f %f %f\n",
              proj_pt_coord[0], proj_pt_coord[1], proj_pt_coord[2]);
  }

  return iter+1;
}





static void
_bezier_to_lagrange
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   elt_order,
 const int                   n_vtx,
       double               *vtx_coord,
 const int                   n_elt,
       int                  *elt_vtx
 )
{
  double *lag = malloc(sizeof(double) * n_vtx * 3);

  int *is_set = PDM_array_zeros_int(n_vtx);

  int stride = PDM_Mesh_nodal_n_vtx_elt_get(elt_type,
                                            elt_order);
  double *ec = malloc(sizeof(double) * stride * 3);

  // double step = 1. / (double) elt_order;

  int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

  double *weight   = malloc(sizeof(double) * stride);
  double *uvw_node = malloc(sizeof(double) * stride * elt_dim);

  if (elt_type != PDM_MESH_NODAL_BARHO_BEZIER &&
      elt_type != PDM_MESH_NODAL_TRIAHO_BEZIER) {
    PDM_error(__FILE__, __LINE__, 0,
              "Béizer -> Lagrange not yet implemented for type %d\n", (int) elt_type);
  }

  PDM_ho_location_uvw_nodes(elt_type,
                            elt_order,
                            0, 1,
                            0, 1,
                            0, 1,
                            uvw_node);


  for (int ielt = 0; ielt < n_elt; ielt++) {

    int *ev = elt_vtx + stride*ielt;

    for (int i = 0; i < stride; i++) {
      int id_vtx = ev[i] - 1;
      memcpy(ec + 3*i, vtx_coord + 3*id_vtx, sizeof(double) * 3);
    }

    // int k = 0;
    // for (int j = 0; j <= elt_order; j++) {
    //   for (int i = 0; i <= elt_order - j; i++) {
    //     int id_vtx = ev[k++] - 1;

    //     if (!is_set[id_vtx]) {
    //       PDM_ho_bezier_de_casteljau_triangle(3,
    //                                           elt_order,
    //                                           i*step,
    //                                           j*step,
    //                                           ec,
    //                                           lag + 3*id_vtx,
    //                                           NULL, NULL, NULL);
    //       is_set[id_vtx] = 1;
    //     }
    //   }
    // }
    for (int i = 0; i < stride; i++) {
      int id_vtx = ev[i] - 1;

      if (!is_set[id_vtx]) {

        PDM_ho_bezier_basis(elt_type,
                            elt_order,
                            1,
                            uvw_node + elt_dim*i,
                            weight);

        double *x = lag + 3*id_vtx;
        x[0] = x[1] = x[2] = 0;
        for (int j = 0; j < stride; j++) {
          int jvtx = ev[j] - 1;
          for (int k = 0; k < 3; k++) {
            x[k] += weight[j] * vtx_coord[3*jvtx + k];
          }
        }

        is_set[id_vtx] = 1;
      }
    }

  }
  free(ec);
  free(is_set);
  free(weight);
  free(uvw_node);

  memcpy(vtx_coord, lag, sizeof(double) * n_vtx * 3);
  free(lag);
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
  char *filename = NULL;
  double pt_coord[3] = {0.23062829, 0.39873915, -0.02889847};
  int random_seed = -1;
  int visu        = 0;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &filename,
             pt_coord,
             &random_seed,
             &visu);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"back_faces_P1.dat";
  }

  srand(random_seed);
  // if (random_seed >= 0) {
  //   for (int i = 0; i < 3; i++) {
  //     pt_coord[i] = 2 * ((double) rand() / (double) RAND_MAX) - 1;
  //   }
  // }

  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);



  /*
   *  Read back mesh
   */
  PDM_Mesh_nodal_elt_t  elt_type;
  int                   elt_order;
  int                   n_vtx     = 0;
  double               *vtx_coord = NULL;
  int                   n_elt     = 0;
  int                  *elt_vtx   = NULL;
  _read_back_mesh(filename,
                  &elt_type,
                  &elt_order,
                  &n_vtx,
                  &vtx_coord,
                  &n_elt,
                  &elt_vtx);

  if (random_seed >= 0) {

    int id_vtx = rand() % n_vtx;
    double noise = 0.15;

    for (int i = 0; i < 3; i++) {
      pt_coord[i] = vtx_coord[3*id_vtx + i] + noise*(2*((double) rand()/(double) RAND_MAX) - 1);
    }
  }

  // log_trace("elt_type = %d, elt_order = %d\n", (int) elt_type, elt_order);

  int stride = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);
  int *parent_node = NULL;
  if (elt_order > 1) {
    parent_node = malloc(sizeof(int) * stride);
    PDM_Mesh_nodal_ho_parent_node(elt_type,
                                  elt_order,
                                  NULL,
                                  parent_node);

    if (0) {
      int stride_ho = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, elt_order);
      for (int i = 0; i < n_elt; i++) {
        log_trace("elt %d, principal nodes : ", i);
        for (int j = 0; j < stride; j++) {
          log_trace("%d ", elt_vtx[stride_ho*i + parent_node[j]] - 1);
        }
        log_trace("\n");
      }
    }
  }

  int *elt_elt = malloc(sizeof(int) * stride * n_elt);
  int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
  if (elt_dim == 1) {
    _build_edge_edge_connectivity(elt_order,
                                  n_vtx,
                                  n_elt,
                                  elt_vtx,
                                  parent_node,
                                  elt_elt);
  }
  else if (elt_dim == 2) {
    _build_face_face_connectivity(elt_type,
                                  elt_order,
                                  n_vtx,
                                  n_elt,
                                  elt_vtx,
                                  parent_node,
                                  elt_elt);
  }
  if (parent_node != NULL) free(parent_node);


  if (elt_type == PDM_MESH_NODAL_BARHO_BEZIER) {
    if (0) {
      _bezier_to_lagrange(elt_type,
                          elt_order,
                          n_vtx,
                          vtx_coord,
                          n_elt,
                          elt_vtx);
    }
    elt_type = PDM_MESH_NODAL_BARHO;
  }






  /*
   *  Projection
   */
  // PDM_log_trace_array_double(pt_coord, 3, "pt_coord : ");
  int stride_ho = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, elt_order);
  int *elt_vtx_idx = PDM_array_new_idx_from_const_stride_int(stride_ho, n_elt);

  double  proj_pt_coord[3];
  int     closest_back_elt;
  int    *history_elt  = malloc(sizeof(int) * n_elt);
  double *history_proj = malloc(sizeof(int) * n_elt * 3);

  // Pick a random elt to start from
  int start = rand() % n_elt;

  int n_step = _projection_on_background_mesh_get2(pt_coord,
                                                   n_elt,
                                                   start,
                                                   elt_vtx_idx,
                                                   elt_vtx,
                                                   elt_elt,
                                                   vtx_coord,
                                                   elt_order,
                                                   elt_type,
                                                   proj_pt_coord,
                                                   &closest_back_elt,
                                                   history_elt,
                                                   history_proj);


  /*
   *  Visu VTK
   */
  if (elt_order > 1) {
    PDM_Mesh_nodal_reorder_elt_vtx(elt_type,
                                   elt_order,
                                   NULL,
                                   "PDM_HO_ORDERING_VTK",
                                   n_elt,
                                   elt_vtx,
                                   elt_vtx);
  }


  double *visiting_order = malloc(sizeof(double) * n_elt);
  for (int i = 0; i < n_elt; i++) {
    visiting_order[i] = -n_step;
  }
  for (int i = 0; i < n_step; i++) {
    visiting_order[history_elt[i]-1] = i+1;
  }
  const char* field_name[] = {"visiting_order"};
  const double *field[1] = {visiting_order};

  if (visu) {
    PDM_vtk_write_std_elements_ho("back_mesh.vtk",
                                  elt_order,
                                  n_vtx,
                                  vtx_coord,
                                  NULL,
                                  elt_type,
                                  n_elt,
                                  elt_vtx,
                                  NULL,
                                  1,
                                  field_name,
                                  field);
  }
  free(visiting_order);



  double proj_p[6] = {
    pt_coord[0], pt_coord[1], pt_coord[2],
    proj_pt_coord[0], proj_pt_coord[1], proj_pt_coord[2]
  };

  int proj_e[2] = {1, 2};

  if (visu) {
    PDM_vtk_write_std_elements("proj.vtk",
                               2,
                               proj_p,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               1,
                               proj_e,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }

  int *history_bar = malloc(sizeof(int) * (n_step-1) * 2);
  for (int i = 0; i < n_step-1; i++) {
    history_bar[2*i  ] = i+1;
    history_bar[2*i+1] = i+2;
  }

  if (visu) {
    PDM_vtk_write_std_elements("history.vtk",
                               n_step,
                               history_proj,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               n_step-1,
                               history_bar,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }
  free(history_bar);



  free(vtx_coord);
  free(elt_vtx);
  free(elt_vtx_idx);
  free(elt_elt);
  free(history_elt);
  free(history_proj);

  PDM_MPI_Finalize();

  return 0;
}
