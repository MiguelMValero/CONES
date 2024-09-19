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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_writer.h"
#include "pdm_vtk.h"
#include "pdm_mesh_nodal.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_unique.h"
#include "pdm_sort.h"


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
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                   argc,
 char                **argv,
 PDM_g_num_t          *n,
 int                  *n_layer,
 double               *extrusion_vector,
 PDM_Mesh_nodal_elt_t *elt_type,
 int                  *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        *n = (PDM_g_num_t) _n;
      }
    }

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_layer = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-e") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          extrusion_vector[j] = atof(argv[i]);
        }
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



// static void
// _gen_polyline
// (
//  const PDM_g_num_t   n,
//  int                *base_n_edge,
//  int               **base_edge_vtx,
//  int                *base_n_vtx,
//  double            **base_vtx_coord
//  )
// {
//   *base_n_edge = (int) n;
//   *base_n_vtx  = (int) (n + 1);

//   int *rand_edge = malloc(sizeof(int) * (*base_n_edge));
//   int *perm_edge = malloc(sizeof(int) * (*base_n_edge));
//   for (int i = 0; i < *base_n_edge; i++) {
//     rand_edge[i] = rand();
//     perm_edge[i] = i;
//   }

//   int *rand_vtx = malloc(sizeof(int) * (*base_n_vtx));
//   int *perm_vtx = malloc(sizeof(int) * (*base_n_vtx));
//   for (int i = 0; i < *base_n_vtx; i++) {
//     rand_vtx[i] = rand();
//     perm_vtx[i] = i;
//   }

//   if (1) {
//     PDM_sort_int(rand_edge, perm_edge, *base_n_edge);
//     PDM_sort_int(rand_vtx,  perm_vtx,  *base_n_vtx);
//   }
//   free(rand_edge);
//   free(rand_vtx);

//   *base_edge_vtx = malloc(sizeof(int) * 2 * (*base_n_edge));
//   for (int i = 0; i < *base_n_edge; i++) {
//     int iedge = perm_edge[i];
//     (*base_edge_vtx)[2*iedge  ] = perm_vtx[i  ]+1;
//     (*base_edge_vtx)[2*iedge+1] = perm_vtx[i+1]+1;
//   }


//   *base_vtx_coord = malloc(sizeof(double) * 3 * (*base_n_vtx));
//   double step = 2*PDM_PI / (double) n;
//   for (int i = 0; i < *base_n_vtx; i++) {
//     double x = -PDM_PI + i*step;
//     (*base_vtx_coord)[3*perm_vtx[i]  ] = x;
//     (*base_vtx_coord)[3*perm_vtx[i]+1] = 0;
//     (*base_vtx_coord)[3*perm_vtx[i]+2] = atan(x);
//   }

//   free(perm_edge);
//   free(perm_vtx);
// }


static void
_gen_circle
(
 const PDM_g_num_t   n,
 int                *base_n_edge,
 int               **base_edge_vtx,
 int                *base_n_vtx,
 double            **base_vtx_coord
 )
{
  *base_n_edge = (int) n;
  *base_n_vtx  = (int) (n + 1);

  int *rand_edge = malloc(sizeof(int) * (*base_n_edge));
  int *perm_edge = malloc(sizeof(int) * (*base_n_edge));
  for (int i = 0; i < *base_n_edge; i++) {
    rand_edge[i] = rand();
    perm_edge[i] = i;
  }

  int *rand_vtx = malloc(sizeof(int) * (*base_n_vtx));
  int *perm_vtx = malloc(sizeof(int) * (*base_n_vtx));
  for (int i = 0; i < *base_n_vtx; i++) {
    rand_vtx[i] = rand();
    perm_vtx[i] = i;
  }

  if (1) {
    PDM_sort_int(rand_edge, perm_edge, *base_n_edge);
    PDM_sort_int(rand_vtx,  perm_vtx,  *base_n_vtx);
  }
  free(rand_edge);
  free(rand_vtx);

  *base_edge_vtx = malloc(sizeof(int) * 2 * (*base_n_edge));
  for (int i = 0; i < *base_n_edge; i++) {
    int iedge = perm_edge[i];
    (*base_edge_vtx)[2*iedge  ] = perm_vtx[i  ]+1;
    (*base_edge_vtx)[2*iedge+1] = perm_vtx[i+1]+1;
  }


  *base_vtx_coord = malloc(sizeof(double) * 3 * (*base_n_vtx));
  double step = 2*PDM_PI / (double) n;
  for (int i = 0; i < *base_n_vtx; i++) {
    double t = i*step;
    (*base_vtx_coord)[3*perm_vtx[i]  ] = cos(t);
    (*base_vtx_coord)[3*perm_vtx[i]+1] = sin(t);
    (*base_vtx_coord)[3*perm_vtx[i]+2] = 0;
  }

  free(perm_edge);
  free(perm_vtx);
}


static void
_extrude_polyline
(
PDM_MPI_Comm                 comm,
const int                    base_n_edge,
const int                   *base_edge_vtx,
const int                    base_n_vtx,
const double                *base_vtx_coord,
const int                    n_layer,
      double                *extrusion_vector,
const PDM_Mesh_nodal_elt_t   elt_type,
      int                   *n_face,
      int                  **face_vtx_idx,
      int                  **face_vtx,
      PDM_g_num_t          **face_ln_to_gn,
      int                   *n_vtx,
      double               **vtx_coord,
      PDM_g_num_t          **vtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t gn_face = base_n_edge * n_layer;
  if (elt_type == PDM_MESH_NODAL_TRIA3) {
    gn_face *= 2;
  }


  PDM_g_num_t *distrib_face = PDM_compute_uniform_entity_distribution(comm,
                                                                      gn_face);

  *n_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  *face_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_face));

  int face_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);
  *face_vtx_idx = PDM_array_new_idx_from_const_stride_int(face_vtx_n, *n_face);

  int s_face_vtx = face_vtx_n * (*n_face);

  PDM_g_num_t *face_vtx_gnum = malloc(sizeof(PDM_g_num_t) * s_face_vtx);

  int layer_face_n = base_n_edge;
  if (elt_type == PDM_MESH_NODAL_TRIA3) {
    layer_face_n *= 2;
  }


  for (int iface = 0; iface < *n_face; iface++) {
    PDM_g_num_t g = distrib_face[i_rank] + iface;

    PDM_g_num_t *fv = face_vtx_gnum + (*face_vtx_idx)[iface];

    (*face_ln_to_gn)[iface] = g + 1;

    int ilayer = (int) (g / layer_face_n);
    int ibase_edge;
    if (elt_type == PDM_MESH_NODAL_TRIA3) {
      ibase_edge = (int) ((g/2) % base_n_edge);
    }
    else {
      ibase_edge = (int) (g % base_n_edge);
    }

    int ibase_vtx1 = base_edge_vtx[2*ibase_edge  ];
    int ibase_vtx2 = base_edge_vtx[2*ibase_edge+1];

    if (elt_type == PDM_MESH_NODAL_TRIA3) {

      if (g%2 == 0) {
        fv[0] = ilayer*base_n_vtx + ibase_vtx1;
        fv[1] = ilayer*base_n_vtx + ibase_vtx2;
        fv[2] = fv[1] + base_n_vtx;
      }
      else {
        fv[0] = ilayer*base_n_vtx + ibase_vtx1;
        fv[1] = (ilayer+1)*base_n_vtx + ibase_vtx2;
        fv[2] = fv[0] + base_n_vtx;
      }

    }
    else if (elt_type == PDM_MESH_NODAL_QUAD4) {

      fv[0] = ilayer*base_n_vtx + ibase_vtx1;
      fv[1] = ilayer*base_n_vtx + ibase_vtx2;
      fv[2] = fv[1] + base_n_vtx;
      fv[3] = fv[0] + base_n_vtx;

    }
    else {
      abort();
    }
  }
  free(distrib_face);


  /* Unique vertices */
  *face_vtx = malloc(sizeof(int) * s_face_vtx);
  *n_vtx = PDM_inplace_unique_long2(face_vtx_gnum,
                                    *face_vtx,
                                    0,
                                    s_face_vtx-1);

  for (int i = 0; i < s_face_vtx; i++) {
    (*face_vtx)[i]++;
  }


  *vtx_ln_to_gn = realloc(face_vtx_gnum, sizeof(PDM_g_num_t) * (*n_vtx));
  *vtx_coord = malloc(sizeof(double) * (*n_vtx) * 3);

  double step = 1. / (double) n_layer;

  for (int ivtx = 0; ivtx < (*n_vtx); ivtx++) {

    PDM_g_num_t g = (*vtx_ln_to_gn)[ivtx] - 1;

    int ilayer    = (int) (g / base_n_vtx);
    int ibase_vtx = (int) (g % base_n_vtx);

    for (int i = 0; i < 3; i++) {
      (*vtx_coord)[3*ivtx+i] = base_vtx_coord[3*ibase_vtx+i] + ilayer*step*extrusion_vector[i];
    }

  }

}


/**
 *
 * \brief  Main
 *
 */

int
main
(
 int   argc,
 char *argv[]
 )
{
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t          n                   = 10;
  int                  n_layer             = 5;
  double               extrusion_vector[3] = {0., 1., 0.};
  PDM_Mesh_nodal_elt_t elt_type            = PDM_MESH_NODAL_QUAD4;
  int                  visu                = 0;
  _read_args(argc,
             argv,
             &n,
             &n_layer,
             extrusion_vector,
             &elt_type,
             &visu);



  /*
   *  Generate a globally-shared base polyline
   */
  int     base_n_edge;
  int     base_n_vtx;
  int    *base_edge_vtx  = NULL;
  double *base_vtx_coord = NULL;
  // _gen_polyline(n,
  //               &base_n_edge,
  //               &base_edge_vtx,
  //               &base_n_vtx,
  //               &base_vtx_coord);
  _gen_circle(n,
              &base_n_edge,
              &base_edge_vtx,
              &base_n_vtx,
              &base_vtx_coord);

  int     base_n_edge2 = 0;
  int     base_n_vtx2  = 0;
  int    *base_edge_vtx2  = NULL;
  double *base_vtx_coord2 = NULL;
  // _gen_polyline(n/2,
  //               &base_n_edge2,
  //               &base_edge_vtx2,
  //               &base_n_vtx2,
  //               &base_vtx_coord2);

  base_edge_vtx = realloc(base_edge_vtx, sizeof(int) * (base_n_edge + base_n_edge2) * 2);
  for (int i = 0; i < 2*base_n_edge2; i++) {
    base_edge_vtx[2*base_n_edge + i] = base_edge_vtx2[i] + base_n_vtx;
  }
  base_n_edge += base_n_edge2;

  base_vtx_coord = realloc(base_vtx_coord, sizeof(double) * (base_n_vtx + base_n_vtx2) * 3);
  for (int i = 0; i < base_n_vtx2; i++) {
    base_vtx_coord[3*base_n_vtx + 3*i  ] = base_vtx_coord2[3*i  ];
    base_vtx_coord[3*base_n_vtx + 3*i+1] = base_vtx_coord2[3*i+1];
    base_vtx_coord[3*base_n_vtx + 3*i+2] = 1.5 + 0.5*base_vtx_coord2[3*i+2] + 0.25*base_vtx_coord2[3*i  ];
  }
  base_n_vtx += base_n_vtx2;

  free(base_edge_vtx2);
  free(base_vtx_coord2);

  if (visu && i_rank == 0) {
    PDM_vtk_write_std_elements("base_polyline.vtk",
                               base_n_vtx,
                               base_vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               base_n_edge,
                               base_edge_vtx,
                               NULL,
                               0,
                               NULL, NULL);
  }




  /*
   *  Extrusion
   */

  int          n_face        = 0;
  int         *face_vtx_idx  = NULL;
  int         *face_vtx      = NULL;
  PDM_g_num_t *face_ln_to_gn = NULL;
  int          n_vtx         = 0;
  double      *vtx_coord     = NULL;
  PDM_g_num_t *vtx_ln_to_gn  = NULL;
  _extrude_polyline(comm,
                    base_n_edge,
                    base_edge_vtx,
                    base_n_vtx,
                    base_vtx_coord,
                    n_layer,
                    extrusion_vector,
                    elt_type,
                    &n_face,
                    &face_vtx_idx,
                    &face_vtx,
                    &face_ln_to_gn,
                    &n_vtx,
                    &vtx_coord,
                    &vtx_ln_to_gn);


  if (visu) {
    PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_BIN,
                                          PDM_WRITER_TOPO_CST,
                                          PDM_WRITER_OFF,
                                          "extruded_polyline",
                                          "extruded_polyline",
                                          PDM_MPI_COMM_WORLD,
                                          PDM_IO_KIND_MPI_SIMPLE,
                                          1.,
                                          NULL);

    int id_geom = PDM_writer_geom_create(wrt,
                                         "extruded_polyline",
                                         1);

    int id_var_part = PDM_writer_var_create(wrt,
                                            PDM_WRITER_ON,
                                            PDM_WRITER_VAR_SCALAR,
                                            PDM_WRITER_VAR_ELEMENTS,
                                            "num_part");

    PDM_writer_step_beg(wrt, 0.);

    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              0,
                              n_vtx,
                              vtx_coord,
                              vtx_ln_to_gn,
                              PDM_OWNERSHIP_USER);

    int id_block = PDM_writer_geom_bloc_add(wrt,
                                            id_geom,
                                            PDM_WRITER_POLY_2D,
                                            PDM_OWNERSHIP_USER);

    PDM_writer_geom_bloc_poly2d_set(wrt,
                                    id_geom,
                                    id_block,
                                    0,
                                    n_face,
                                    face_vtx_idx,
                                    face_vtx,
                                    face_ln_to_gn);

    PDM_writer_geom_write(wrt,
                          id_geom);

    PDM_real_t *val_part = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_face);
    for (int i = 0; i < n_face; i++) {
      val_part[i] = (PDM_real_t) i_rank;
    }

    PDM_writer_var_set(wrt,
                       id_var_part,
                       id_geom,
                       0,
                       (const PDM_real_t *) val_part);
    PDM_writer_var_write(wrt,
                         id_var_part);
    PDM_writer_var_free(wrt,
                        id_var_part);
    free(val_part);

    PDM_writer_step_end(wrt);

    PDM_writer_free(wrt);
  }



  /* Free memory */
  free(base_edge_vtx);
  free(base_vtx_coord);

  free(face_vtx_idx );
  free(face_vtx     );
  free(face_ln_to_gn);
  free(vtx_coord    );
  free(vtx_ln_to_gn );

  PDM_MPI_Finalize();

  return 0;
}
