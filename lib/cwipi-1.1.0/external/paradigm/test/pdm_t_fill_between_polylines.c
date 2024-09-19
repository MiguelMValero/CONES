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
#include "pdm_dmesh_nodal.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_unique.h"
#include "pdm_dmesh_nodal_priv.h"


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
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n,
 PDM_g_num_t           *n_layer,
 PDM_Mesh_nodal_elt_t  *elt_type,
 double                *geometric_ratio,
 int                   *visu
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

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-ratio") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *geometric_ratio = atof(argv[i]);
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


static void
_gen_circle
(
 const PDM_g_num_t   n,
       int          *base_n_edge,
       int         **base_edge_vtx,
       int          *base_n_vtx,
       double      **base_vtx_coord
 )
{
  *base_n_edge = (int) n;
  *base_n_vtx  = (int) n;

  *base_edge_vtx = malloc(sizeof(int) * 2 * (*base_n_edge));
  for (int iedge = 0; iedge < *base_n_edge; iedge++) {
    (*base_edge_vtx)[2*iedge  ] = iedge       + 1;
    (*base_edge_vtx)[2*iedge+1] = (iedge+1)%n + 1;
  }

  *base_vtx_coord = malloc(sizeof(double) * 3 * (*base_n_vtx));
  double step = 2*PDM_PI / (double) n;
  for (int ivtx = 0; ivtx < *base_n_vtx; ivtx++) {
    double t = ivtx*step;
    (*base_vtx_coord)[3*ivtx  ] = cos(t);
    (*base_vtx_coord)[3*ivtx+1] = sin(t);
    (*base_vtx_coord)[3*ivtx+2] = 0;
  }
}


static void
_gen_naca_00xx
(
 const PDM_g_num_t   n,
 const double        xx,
       int          *base_n_edge,
       int         **base_edge_vtx,
       int          *base_n_vtx,
       double      **base_vtx_coord
 )
{
  const double a =  0.2969;
  const double b = -0.1260;
  const double c = -0.3516;
  const double d =  0.2843;
  const double e = -0.1036;

  *base_n_edge = (int) n;
  *base_n_vtx  = (int) n;

  *base_edge_vtx = malloc(sizeof(int) * 2 * (*base_n_edge));
  for (int iedge = 0; iedge < *base_n_edge; iedge++) {
    (*base_edge_vtx)[2*iedge  ] = iedge       + 1;
    (*base_edge_vtx)[2*iedge+1] = (iedge+1)%n + 1;
  }

  *base_vtx_coord = malloc(sizeof(double) * 3 * (*base_n_vtx));
  // double step = 2*PDM_PI / (double) (n+1);
  PDM_g_num_t half = n / 2;
  double step = 2*PDM_PI / (double) n;
  for (int ivtx = 0; ivtx < *base_n_vtx; ivtx++) {
    double t;
    if (ivtx < half) {
      t = ivtx*step;
    }
    else {
      t = (ivtx - half)*step + PDM_PI;
    }

    // log_trace("i = %d/%d : t = %f * PI\n", ivtx, *base_n_vtx, t/PDM_PI);
    double x = 0.5*(1 + cos(t));
    double y = 5. * xx * (a*sqrt(x) + x*(b + x*(c + x*(d + x*e))));
    if (ivtx >= half) {
      y = -y;
    }

    (*base_vtx_coord)[3*ivtx  ] = x - 0.5;
    (*base_vtx_coord)[3*ivtx+1] = y;
    (*base_vtx_coord)[3*ivtx+2] = 0;
  }
}




static PDM_dmesh_nodal_t *
_set_dmesh_nodal
(
 const PDM_MPI_Comm          comm,
       double               *dvtx_coord,
       PDM_g_num_t          *dface_vtx,
       PDM_g_num_t          *distrib_vtx,
       PDM_g_num_t          *distrib_face,
       PDM_Mesh_nodal_elt_t  face_type
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                     face_type);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int n_group = 1;
  int *dgroup_elt_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  dgroup_elt_idx[1] = dn_face;

  PDM_g_num_t *dgroup_elt = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  for (int i = 0; i < dn_face; i++) {
    dgroup_elt[i] = distrib_face[i_rank] + i + 1;
  }
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}


/* TO DO: boundary edges from extruded boundry vertices */
static void
_fill_between_polylines
(
      PDM_MPI_Comm           comm,
const int                    base_n_edge,
const int                   *base_edge_vtx,
const int                    base_n_vtx,
const double                *base_vtx_coord[2],
const int                    n_layer,
const PDM_Mesh_nodal_elt_t   elt_type,
      PDM_g_num_t          **distrib_face,
      int                  **dface_vtx_idx,
      PDM_g_num_t          **dface_vtx,
      PDM_g_num_t          **distrib_edge,
      PDM_g_num_t          **dedge_vtx,
      PDM_g_num_t          **distrib_vtx,
      double               **dvtx_coord
 )
{
  int dbg = 0;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);


  double step = 1. / (double) n_layer;


  /* Compute vertex normals */
  double *base_vtx_normal[2] = {NULL, NULL};
  for (int iline = 0; iline < 2; iline++) {
    base_vtx_normal[iline] = malloc(sizeof(double) * base_n_vtx * 3);

    for (int i = 0; i < base_n_vtx * 3; i++) {
      base_vtx_normal[iline][i] = 0.;
    }

    for (int iedge = 0; iedge < base_n_edge; iedge++) {
      int ivtx1 = base_edge_vtx[2*iedge  ] - 1;
      int ivtx2 = base_edge_vtx[2*iedge+1] - 1;
      double dx = base_vtx_coord[iline][3*ivtx2  ] - base_vtx_coord[iline][3*ivtx1  ];
      double dy = base_vtx_coord[iline][3*ivtx2+1] - base_vtx_coord[iline][3*ivtx1+1];

      base_vtx_normal[iline][3*ivtx1  ] +=  dy;
      base_vtx_normal[iline][3*ivtx2  ] +=  dy;
      base_vtx_normal[iline][3*ivtx1+1] += -dx;
      base_vtx_normal[iline][3*ivtx2+1] += -dx;
    }

    /* Scale normals appropriately */
    for (int i = 0; i < base_n_vtx; i++) {
      double mag_normal = PDM_DOT_PRODUCT(base_vtx_normal[iline] + 3*i,
                                          base_vtx_normal[iline] + 3*i);
      assert(mag_normal >= 0);

      double dist = 0.;
      for (int j = 0; j < 3; j++) {
        double delta = base_vtx_coord[0][3*i+j] - base_vtx_coord[1][3*i+j];
        dist += delta*delta;
      }

      double scale = 2*step * sqrt(dist / mag_normal);
      for (int j = 0; j < 3; j++) {
        base_vtx_normal[iline][3*i+j] *= scale;
      }
    }
  }

  if (dbg) {
    char filename[999];
    for (int i = 0; i < 2; i++) {
      sprintf(filename, "base_vtx_normal_%d.vtk", i);

      const char   *field_name [1] = {"base_vtx_normal"};
      const double *field_value[1] = {base_vtx_normal[i]};
      PDM_vtk_write_point_cloud_with_field(filename,
                                           base_n_vtx,
                                           base_vtx_coord[i],
                                           NULL,
                                           NULL,
                                           0, NULL, NULL,
                                           1,
                                           field_name,
                                           field_value,
                                           0, NULL, NULL);
    }
  }


  /* Build vertices */
  PDM_g_num_t gn_vtx = base_n_vtx * (n_layer + 1);

  *distrib_vtx = PDM_compute_uniform_entity_distribution(comm,
                                                         gn_vtx);

  int dn_vtx = (int) ((*distrib_vtx)[i_rank+1] - (*distrib_vtx)[i_rank]);

  *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);

  for (int ivtx = 0; ivtx < dn_vtx; ivtx++) {

    double *vc = *dvtx_coord + 3*ivtx;

    PDM_g_num_t g = (*distrib_vtx)[i_rank] + ivtx;

    PDM_g_num_t ilayer = g/base_n_vtx;
    int ibase_vtx = g%base_n_vtx;

    double t = ilayer * step;

    /* Hermite cubic interpolation */
    const double *vc1 = base_vtx_coord[0] + 3*ibase_vtx;
    const double *vc2 = base_vtx_coord[1] + 3*ibase_vtx;

    double *vn1 = base_vtx_normal[0] + 3*ibase_vtx;
    double *vn2 = base_vtx_normal[1] + 3*ibase_vtx;

    double h00 = 1 + t*(t*(-3 + 2*t));
    double h10 = t*(1 + t*(-2 + t));
    double h01 = t*t*(3 - 2*t);
    double h11 = t*t*(-1 + t);

    for (int i = 0; i < 3; i++) {
      vc[i] = h00*vc1[i] + h10*vn1[i] + h01*vc2[i] + h11*vn2[i];
    }

  } // End of loop on vertices
  free(base_vtx_normal[0]);
  free(base_vtx_normal[1]);



  /* Build faces */
  PDM_g_num_t gn_face = base_n_edge * n_layer;

  if (elt_type == PDM_MESH_NODAL_TRIA3) {
    gn_face *= 2;
  }


  *distrib_face = PDM_compute_uniform_entity_distribution(comm,
                                                          gn_face);

  int dn_face = (int) ((*distrib_face)[i_rank+1] - (*distrib_face)[i_rank]);

  int face_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);
  *dface_vtx_idx = PDM_array_new_idx_from_const_stride_int(face_vtx_n, dn_face);

  int s_face_vtx = face_vtx_n * dn_face;

  *dface_vtx = malloc(sizeof(PDM_g_num_t) * s_face_vtx);

  int layer_face_n = base_n_edge;
  if (elt_type == PDM_MESH_NODAL_TRIA3) {
    layer_face_n *= 2;
  }


  for (int iface = 0; iface < dn_face; iface++) {

    PDM_g_num_t *fv = *dface_vtx + (*dface_vtx_idx)[iface];

    PDM_g_num_t g = (*distrib_face)[i_rank] + iface;

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

  } // End of loop on faces


  *distrib_edge = NULL;
  *dedge_vtx = NULL;


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


  /*
   *  Set default values
   */
  PDM_g_num_t          n               = 4;
  PDM_g_num_t          n_layer         = 3;
  PDM_Mesh_nodal_elt_t elt_type        = PDM_MESH_NODAL_QUAD4;
  double               geometric_ratio = 1.;
  int                  visu            = 0;

  _read_args(argc,
             argv,
             &n,
             &n_layer,
             &elt_type,
             &geometric_ratio,
             &visu);






  int     base_n_edge;
  int     base_n_vtx;
  int    *base_edge_vtx   = NULL;
  double *base_vtx_coord1 = NULL;
  // _gen_circle(n,
  //             &base_n_edge,
  //             &base_edge_vtx,
  //             &base_n_vtx,
  //             &base_vtx_coord1);
  _gen_naca_00xx(n,
                 0.12,
                 &base_n_edge,
                 &base_edge_vtx,
                 &base_n_vtx,
                 &base_vtx_coord1);
  free(base_edge_vtx);

  double *base_vtx_coord2 = NULL;
  _gen_circle(n,
              &base_n_edge,
              &base_edge_vtx,
              &base_n_vtx,
              &base_vtx_coord2);

  for (int i = 0; i < base_n_vtx; i++) {
    base_vtx_coord2[3*i  ] *= 150;
    base_vtx_coord2[3*i+1] *= 150;
  }


  const double *base_vtx_coord[2] = {base_vtx_coord1, base_vtx_coord2};

  if (i_rank == 0 && visu) {
    char filename[999];
    for (int i = 0; i < 2; i++) {
      sprintf(filename, "base_polyline_%d.vtk", i);

      PDM_vtk_write_std_elements(filename,
                                 base_n_vtx,
                                 base_vtx_coord[i],
                                 NULL,
                                 PDM_MESH_NODAL_BAR2,
                                 base_n_edge,
                                 base_edge_vtx,
                                 NULL,
                                 0, NULL, NULL);
    }
  }

  PDM_g_num_t *distrib_face  = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  PDM_g_num_t *distrib_edge  = NULL;
  PDM_g_num_t *dedge_vtx     = NULL;
  PDM_g_num_t *distrib_vtx   = NULL;
  double      *dvtx_coord    = NULL;

  _fill_between_polylines(comm,
                          base_n_edge,
                          base_edge_vtx,
                          base_n_vtx,
                          base_vtx_coord,
                          n_layer,
                          elt_type,
                          &distrib_face,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &distrib_edge,
                          &dedge_vtx,
                          &distrib_vtx,
                          &dvtx_coord);

  free(base_edge_vtx);
  free(base_vtx_coord1);
  free(base_vtx_coord2);



  PDM_dmesh_nodal_t *dmn = _set_dmesh_nodal(comm,
                                            dvtx_coord,
                                            dface_vtx,
                                            distrib_vtx,
                                            distrib_face,
                                            elt_type);



  if (visu) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_RIDGE,
                             "fill_ridge_");

    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "fill_surface_");
  }
  PDM_DMesh_nodal_free(dmn);

  free(distrib_face );
  free(dface_vtx_idx);
  free(distrib_vtx  );

  PDM_MPI_Finalize();

  return 0;
}

