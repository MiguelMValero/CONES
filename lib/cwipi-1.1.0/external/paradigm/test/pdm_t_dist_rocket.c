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

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_timer.h"
#include "pdm_part.h"
#include "pdm_geom_elem.h"
#include "pdm_dcube_gen.h"
#include "pdm_mpi_node_first_rank.h"

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
     "  -n      <level>  Number of vertices.\n\n"
     "  -l      <level>  Rocket length.\n\n"
     "  -t      <level>  Number of Target points (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_faceSeg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   nTgt       Number of Target points
 * \param [inout]   n_part     Number of partitions par process
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nv,
           double        *length,
           double        *domain_size,
           PDM_g_num_t   *nTgt,
           int           *grid,
           int           *n_max_per_leaf,
           int           *on_ground,
           int           *n_proc_data_src,
           int           *post)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nv") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nv = atol(argv[i]);
        *nv = (PDM_g_num_t) _nv;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
      }
    }
    else if (strcmp(argv[i], "-grid") == 0) {
      *grid = 1;
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *domain_size = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-mpl") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_max_per_leaf = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-ground") == 0) {
      *on_ground = 1;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_proc_data_src = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double _rand01(void) {
  return (double) rand() / (double) RAND_MAX;
}

static void
_gen_cloud_random
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_pts,
 const double        origin[3],
 const double        length,
 int                *_n_pts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  // Define distribution
  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  distrib[0] = 0;
  PDM_g_num_t step = n_pts / n_rank;
  PDM_g_num_t remainder = n_pts % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] = step;
    const int i1 = i - 1;
    if (i1 < remainder) {
      distrib[i]++;
    }
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] += distrib[i-1];
  }

  PDM_g_num_t dn_pts = distrib[i_rank+1] - distrib[i_rank];
  *_n_pts = (int) dn_pts;

  *g_num = malloc (sizeof(PDM_g_num_t) * dn_pts);
  *coord = malloc (sizeof(double)      * dn_pts * 3);
  for (int i = 0; i < *_n_pts; i++) {
    (*g_num)[i] = 1 + i + distrib[i_rank];
    for (int j = 0; j < 3; j++) {
      (*coord)[3*i+j] = origin[j] + length * _rand01();
    }
  }

  free (distrib);
}



static void
_gen_cloud_grid
(
 PDM_MPI_Comm   comm,
 int            n_part,
 PDM_g_num_t    n_vtx_seg,
 double         origin[3],
 double         length,
 int          **n_pts,
 PDM_g_num_t ***pts_g_num,
 double      ***pts_coord
 )
{
  int n_rank, i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_t *dcube = PDM_dcube_gen_init (comm,
                                           n_vtx_seg,
                                           length,
                                           origin[0],
                                           origin[1],
                                           origin[2],
                                           PDM_OWNERSHIP_KEEP);

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  int          dface_vtx_l;
  int          dface_group_l;

  PDM_dcube_gen_dim_get (dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);


  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  PDM_dcube_gen_data_get (dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  /*
   *  Create mesh partitions
   */

  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  // int ppart_id = 0;
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc (sizeof(int) * dn_cell);

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;
  PDM_part_t *ppart = PDM_part_create (PDM_MPI_COMM_WORLD,
                                       part_method,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       "PDM_PART_RENUM_FACE_NONE",
                                       n_property_cell,
                                       renum_properties_cell,
                                       n_property_face,
                                       renum_properties_face,
                                       n_part,
                                       dn_cell,
                                       dn_face,
                                       dn_vtx,
                                       n_face_group,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dface_cell,
                                       dface_vtx_idx,
                                       dface_vtx,
                                       NULL,
                                       dvtx_coord,
                                       NULL,
                                       dface_group_idx,
                                       dface_group);

  free(dcell_part);


  *n_pts     = malloc (sizeof(int)           * n_part);
  *pts_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  *pts_coord = malloc (sizeof(double *)      * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_t_part,
                           &s_cell_face,
                           &s_face_vtx,
                           &s_face_group,
                           &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    (*n_pts)[i_part] = n_cell;
    (*pts_g_num)[i_part] = malloc (sizeof(PDM_g_num_t) * n_cell);
    for (int i = 0; i < n_cell; i++) {
      (*pts_g_num)[i_part][i] = cell_ln_to_gn[i];
    }


    const int is_oriented = 0;
    (*pts_coord)[i_part] = malloc (sizeof(double) * n_cell * 3);
    double *cell_volume = malloc (sizeof(double) * n_cell);
    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume,
                                        (*pts_coord)[i_part],
                                        NULL,
                                        NULL);
    free (cell_volume);
  }

  PDM_part_free (ppart);
}








static void
_gen_rocket
(
 PDM_MPI_Comm     comm,
 PDM_g_num_t      nv,
 double           length,
 double           radius,
 PDM_g_num_t     *ng_face,
 PDM_g_num_t     *ng_vtx,
 PDM_g_num_t     *ng_edge,
 int             *dn_vtx,
 double         **dvtx_coord,
 int             *dn_face,
 int            **dface_vtx_idx,
 PDM_g_num_t    **dface_vtx,
 PDM_g_num_t    **dface_edge,
 int             *dn_edge,
 PDM_g_num_t    **dedge_vtx,
 PDM_g_num_t    **dedge_face,
 int             *n_edge_group,
 int            **dedge_group_idx,
 PDM_g_num_t    **dedge_group
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  double length_cap = length / 14.;
  double length_cyl = length - length_cap;

  double C = sqrt(length_cap*length_cap + radius*radius);
  PDM_g_num_t _nu = (int) (0.5 * (nv - 1) * (PDM_PI*radius/3. + C) / length_cyl);
  _nu = PDM_MAX (_nu, 1);

  PDM_g_num_t nu = 6*_nu;

  PDM_g_num_t ng_vtx_cyl    = nu * (nv - 1);
  PDM_g_num_t ng_edge_cyl_u = nu * (nv - 1);
  PDM_g_num_t ng_edge_cyl_v = nu * (nv - 1);
  PDM_g_num_t ng_face_cyl   = nu * (nv - 1);

  PDM_g_num_t nv_cap      = _nu;
  PDM_g_num_t ng_vtx_cap  = 1 + 3*nv_cap*(nv_cap + 1);
  PDM_g_num_t ng_face_cap = 6 * nv_cap * nv_cap;
  PDM_g_num_t ng_edge_cap = ng_face_cap + ng_vtx_cap - 1;

  PDM_g_num_t ng_edge_lim = nu;

  *ng_vtx  = ng_vtx_cyl + ng_vtx_cap;
  *ng_edge = ng_edge_cyl_u + ng_edge_cyl_v + ng_edge_cap;
  *ng_face = ng_face_cyl + ng_face_cap;
  *n_edge_group = 1;

  if (i_rank == 0) {
    printf("ng_vtx = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM"\n", *ng_vtx, *ng_face);
  }

  /* Define distributions */
  PDM_g_num_t *distrib_vtx  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge_lim = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  distrib_vtx[0]      = 0;
  distrib_edge[0]     = 0;
  distrib_face[0]     = 0;
  distrib_edge_lim[0] = 0;

  PDM_g_num_t step_vtx       = *ng_vtx / n_rank;
  PDM_g_num_t remainder_vtx  = *ng_vtx % n_rank;

  PDM_g_num_t step_edge      = *ng_edge / n_rank;
  PDM_g_num_t remainder_edge = *ng_edge % n_rank;

  PDM_g_num_t step_face      = *ng_face / n_rank;
  PDM_g_num_t remainder_face = *ng_face % n_rank;

  PDM_g_num_t step_edge_lim      = ng_edge_lim / n_rank;
  PDM_g_num_t remainder_edge_lim = ng_edge_lim % n_rank;

  for (int i = 0; i < n_rank; i++) {
    distrib_vtx[i+1] = distrib_vtx[i] + step_vtx;
    if (i < remainder_vtx) {
      distrib_vtx[i+1]++;
    }

    distrib_edge[i+1] = distrib_edge[i] + step_edge;
    if (i < remainder_edge) {
      distrib_edge[i+1]++;
    }

    distrib_face[i+1] = distrib_face[i] + step_face;
    if (i < remainder_face) {
      distrib_face[i+1]++;
    }

    distrib_edge_lim[i+1] = distrib_edge_lim[i] + step_edge_lim;
    if (i < remainder_edge_lim) {
      distrib_edge_lim[i+1]++;
    }
  }
  *dn_vtx  = (int) distrib_vtx[i_rank+1]  - distrib_vtx[i_rank];
  *dn_edge = (int) distrib_edge[i_rank+1] - distrib_edge[i_rank];
  *dn_face = (int) distrib_face[i_rank+1] - distrib_face[i_rank];
  int dn_edge_lim = (int) distrib_edge_lim[i_rank+1] - distrib_edge_lim[i_rank];


  /*
   *  Vertices
   */
  *dvtx_coord = malloc (sizeof(double) * (*dn_vtx) * 3);
  double *_dvtx_coord = *dvtx_coord;

  double step_u = 2.*PDM_PI / (double) nu;
  double step_v = length_cyl / (double) (nv - 1);

  PDM_g_num_t b_vtx_v;
  PDM_g_num_t r_vtx_v;
  PDM_g_num_t ivtx = 0;

  /* Cylinder */
  if (distrib_vtx[i_rank] < ng_vtx_cyl) {
    b_vtx_v = distrib_vtx[i_rank] / nu;
    r_vtx_v = distrib_vtx[i_rank] % nu;

    for (PDM_g_num_t j = b_vtx_v; j < nv-1; j++) {

      PDM_g_num_t _b_vtx_u = 0;
      if (j == b_vtx_v) {
        _b_vtx_u = r_vtx_v;
      }

      double v = j * step_v;

      for (PDM_g_num_t i = _b_vtx_u; i < nu; i++) {
        double u = i * step_u;
        _dvtx_coord[3*ivtx    ] = radius * cos(u);
        _dvtx_coord[3*ivtx + 1] = radius * sin(u);
        _dvtx_coord[3*ivtx + 2] = v;
        ivtx++;
        if (ivtx == *dn_vtx) break;
      }
      if (ivtx == *dn_vtx) break;
    }
  }

  /* Cap */
  if (ivtx < *dn_vtx) {
    PDM_g_num_t gvtx = ng_vtx_cyl;

    double d = length_cap / 3.;
    double ta = tan(PDM_PI * 60./180.);
    double z2 = length_cap - sqrt(d*d / (1. + ta*ta));
    double r2 = (length_cap - z2) * ta;

    for (PDM_g_num_t j = nv_cap; j > 0; j--) {
      float t = (float) j / (float) nv_cap;

      double b3 = (1. - t) * (1. - t) * (1. - t);
      double b2 = 3. * (1. - t) * (1. - t) * t;
      double b1 = 3. * (1. - t) * t * t;
      double b0 = t * t * t;

      double r = (b0 + b1)*radius + b2*r2;
      double z = length_cyl + b1*d + b2*z2 + b3*length_cap;

      PDM_g_num_t __nu = 6*j;
      step_u = 2.*PDM_PI / (double) __nu;

      for (PDM_g_num_t i = 0; i < __nu; i++) {
        if (gvtx >= distrib_vtx[i_rank]) {
          double u = i * step_u;
          _dvtx_coord[3*ivtx    ] = r * cos(u);
          _dvtx_coord[3*ivtx + 1] = r * sin(u);
          _dvtx_coord[3*ivtx + 2] = z;
          ivtx++;
          if (ivtx == *dn_vtx) break;
        }
        gvtx++;

      }
      if (ivtx == *dn_vtx) break;
    }

    if (ivtx < *dn_vtx) {
      _dvtx_coord[3*ivtx    ] = 0.;
      _dvtx_coord[3*ivtx + 1] = 0.;
      _dvtx_coord[3*ivtx + 2] = length;
      ivtx++;
    }
  }
  free (distrib_vtx);


  /*
   *  Edges
   */
  *dedge_vtx  = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  *dedge_face = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  PDM_g_num_t  *_dedge_vtx = *dedge_vtx;
  PDM_g_num_t  *_dedge_face = *dedge_face;

  PDM_g_num_t iedg = 0;
  PDM_g_num_t ifac;

  /* Cylinder - horizontal */
  if (distrib_edge[i_rank] < ng_edge_cyl_u) {
    const PDM_g_num_t b_edge_uv = distrib_edge[i_rank] / nu;
    const PDM_g_num_t r_edge_uv = distrib_edge[i_rank] % nu;

    for (PDM_g_num_t j = b_edge_uv; j < nv-1; j++) {

      PDM_g_num_t _b_edge_uu = 0;
      if (j == b_edge_uv) {
        _b_edge_uu = r_edge_uv;
      }

      for (PDM_g_num_t i = _b_edge_uu; i < nu; i++) {
        _dedge_vtx[2*iedg    ] = 1 + i        + nu*j;
        _dedge_vtx[2*iedg + 1] = 1 + (i+1)%nu + nu*j;

        _dedge_face[2*iedg] = 1 + i + nu*j;
        if (j == 0) {
          _dedge_face[2*iedg + 1] = 0;
        } else {
          _dedge_face[2*iedg + 1] = 1 + i + nu*(j-1);
        }
        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }

  /* Cylinder - vertical */
  if (iedg < *dn_edge && distrib_edge[i_rank] <= ng_edge_cyl_u + ng_edge_cyl_v) {
    const PDM_g_num_t b_edge_vv = (distrib_edge[i_rank] + iedg - ng_edge_cyl_u) / nu;
    const PDM_g_num_t r_edge_vv = (distrib_edge[i_rank] + iedg - ng_edge_cyl_u) % nu;

    for (PDM_g_num_t j = b_edge_vv; j < nv-1; j++) {

      PDM_g_num_t _b_edge_vu = 0;
      if (j == b_edge_vv) {
        _b_edge_vu = r_edge_vv;
      }

      for (PDM_g_num_t i = _b_edge_vu; i < nu; i++) {
        _dedge_vtx[2*iedg    ] = 1 + i + nu*j;
        _dedge_vtx[2*iedg + 1] = 1 + i + nu*(j+1);

        _dedge_face[2*iedg    ] = 1 + (i+nu-1)%nu + nu*j;
        _dedge_face[2*iedg + 1] = 1 + i           + nu*j;
        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }

  /* Cap */
  if (iedg < *dn_edge) {
    ivtx = ng_vtx_cyl + 1;
    ifac = ng_face_cyl + 1;
    PDM_g_num_t gedg = ng_edge_cyl_u + ng_edge_cyl_v;
    for (PDM_g_num_t k = nv_cap; k > 0; k--) {
      PDM_g_num_t n_vtx_k = 6*k;
      PDM_g_num_t n_vtx_km = PDM_MAX (1, 6*(k-1));

      PDM_g_num_t n_tri_k = 6*(2*k - 1);
      PDM_g_num_t n_tri_kp = 6*(2*k + 1);

      for (int j = 0; j < 6; j++) {
        for (PDM_g_num_t i = 0; i < k; i++) {
          PDM_g_num_t v0 = ivtx + k*j + i;
          PDM_g_num_t q = 2*i + (2*k-1)*j;

          v0 = k*j + i;
          if (gedg >= distrib_edge[i_rank]) {
            _dedge_vtx[2*iedg    ] = ivtx + v0;
            _dedge_vtx[2*iedg + 1] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

            _dedge_face[2*iedg    ] = ifac + (q + n_tri_k - 1)%n_tri_k;
            _dedge_face[2*iedg + 1] = ifac + q;
            iedg++;
            if (iedg == *dn_edge) break;
          }
          gedg++;

          if (gedg >= distrib_edge[i_rank]) {
            _dedge_vtx[2*iedg    ] = ivtx + v0;
            _dedge_vtx[2*iedg + 1] = ivtx + (v0+1)%n_vtx_k;

            _dedge_face[2*iedg] = ifac + q;
            if (k == nv_cap) {
              _dedge_face[2*iedg + 1] = ng_face_cyl - nu + 1 + i + k*j;
            } else {
              _dedge_face[2*iedg + 1] = ifac - n_tri_kp + 1 + 2*i + (2*k+1)*j;
            }
            iedg++;
            if (iedg == *dn_edge) break;
          }
          gedg++;

          if (i < k-1) {
            if (gedg >= distrib_edge[i_rank]) {
              _dedge_vtx[2*iedg    ] = ivtx + (v0+1)%n_vtx_k;
              _dedge_vtx[2*iedg + 1] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

              _dedge_face[2*iedg    ] = ifac + q;
              _dedge_face[2*iedg + 1] = ifac + q + 1;
              iedg++;
              if (iedg == *dn_edge) break;
            }
            gedg++;
          }
        }
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;

      ivtx += n_vtx_k;
      ifac += n_tri_k;
    }
  }
  free (distrib_edge);


  /* Edge groups */
  *dedge_group_idx = malloc (sizeof(int) * (*n_edge_group + 1));
  int *_dedge_group_idx = *dedge_group_idx;
  _dedge_group_idx[0] = 0;
  _dedge_group_idx[1] = 0;

  *dedge_group = malloc (sizeof(PDM_g_num_t) * dn_edge_lim);
  PDM_g_num_t *_dedge_group = *dedge_group;

  for (PDM_g_num_t i = distrib_edge_lim[i_rank]; i < distrib_edge_lim[i_rank+1]; i++) {
    _dedge_group[_dedge_group_idx[1]++] = 1 + i;
  }
  free (distrib_edge_lim);

  /*
   *  Faces
   */
  *dface_vtx_idx = malloc (sizeof(int) * (*dn_face + 1));
  int *_dface_vtx_idx = *dface_vtx_idx;
  _dface_vtx_idx[0] = 0;

  int n_quad;
  if (distrib_face[i_rank+1] <= ng_face_cyl) {
    n_quad = *dn_face;
  } else {
    n_quad = (int) PDM_MAX (0, ng_face_cyl - distrib_face[i_rank]);
  }

  for (int i = 0; i < n_quad; i++) {
    _dface_vtx_idx[i+1] = 4 + _dface_vtx_idx[i];
  }
  for (int i = n_quad; i < *dn_face; i++) {
    _dface_vtx_idx[i+1] = 3 + _dface_vtx_idx[i];
  }


  *dface_vtx  = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  *dface_edge = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  PDM_g_num_t *_dface_vtx = *dface_vtx;
  PDM_g_num_t *_dface_edge = *dface_edge;


  ifac = 0;
  /* Cylinder */
  if (n_quad > 0) {
    const PDM_g_num_t b_face_v = distrib_face[i_rank] / nu;
    const PDM_g_num_t r_face_v = distrib_face[i_rank] % nu;

    for (PDM_g_num_t j = b_face_v; j < nv-1; j++) {

      PDM_g_num_t _b_face_u = 0;
      if (j == b_face_v) {
        _b_face_u = r_face_v;
      }

      for (PDM_g_num_t i = _b_face_u; i < nu; i++) {

        PDM_g_num_t ip = (i+1)%nu;

        _dface_vtx[4*ifac    ] = 1 + i  + nu*j;
        _dface_vtx[4*ifac + 1] = 1 + ip + nu*j;
        _dface_vtx[4*ifac + 2] = 1 + ip + nu*(j+1);
        _dface_vtx[4*ifac + 3] = 1 + i  + nu*(j+1);

        _dface_edge[4*ifac    ] = 1 + i  + nu*j;
        _dface_edge[4*ifac + 1] = 1 + ng_edge_cyl_u + ip + nu*j;
        if (j == nv-2) {
          PDM_g_num_t _j = i / nv_cap;
          PDM_g_num_t _i = i % nv_cap;
          _dface_edge[4*ifac + 2] = -(ng_edge_cyl_u + ng_edge_cyl_v + 3*_i + (3*nv_cap-1)*_j + 2);
        } else {
          _dface_edge[4*ifac + 2] = -(1 + i  + nu*(j+1));
        }
        _dface_edge[4*ifac + 3] = -(1 + ng_edge_cyl_u + i + nu*j);

        ifac++;
        if (ifac == n_quad) break;
      }
      if (ifac == n_quad) break;
    }
  }

  /* Cap */
  if (ifac < *dn_face) {
    ivtx = ng_vtx_cyl + 1;
    iedg = ng_edge_cyl_u + ng_edge_cyl_v + 1;
    PDM_g_num_t gfac = ng_face_cyl;

    for (PDM_g_num_t k = nv_cap; k > 0; k--) {
      PDM_g_num_t n_vtx_k = 6*k;
      PDM_g_num_t n_vtx_km = PDM_MAX (1, n_vtx_k - 6);
      PDM_g_num_t n_edge_k = 18*k - 6;

      for (int j = 0; j < 6; j++) {
        for (PDM_g_num_t i = 0; i < k; i++) {
          if (gfac >= distrib_face[i_rank]) {
            _dface_vtx[n_quad + 3*ifac    ] = ivtx + k*j + i;
            _dface_vtx[n_quad + 3*ifac + 1] = ivtx + (k*j + i + 1)%n_vtx_k;
            _dface_vtx[n_quad + 3*ifac + 2] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

            PDM_g_num_t e = 3*i + (3*k-1)*j;
            _dface_edge[n_quad + 3*ifac    ] = iedg + e + 1;
            _dface_edge[n_quad + 3*ifac + 1] = iedg + (e+2)%n_edge_k;
            _dface_edge[n_quad + 3*ifac + 2] = iedg + e;
            ifac++;
            if (ifac == *dn_face) break;
          }
          gfac++;

          if (i < k-1) {
            if (gfac >= distrib_face[i_rank]) {
              _dface_vtx[n_quad + 3*ifac    ] = ivtx + k*j + i + 1;
              _dface_vtx[n_quad + 3*ifac + 1] = ivtx + n_vtx_k + ((k-1)*j + i + 1)%n_vtx_km;
              _dface_vtx[n_quad + 3*ifac + 2] = ivtx + n_vtx_k + (k-1)*j + i;

              PDM_g_num_t e = 3*i + (3*k-1)*j;
              _dface_edge[n_quad + 3*ifac    ] = iedg + (e+3)%n_edge_k;
              _dface_edge[n_quad + 3*ifac + 1] = iedg + n_edge_k + 1 + 3*i + (3*k-4)*j;
              _dface_edge[n_quad + 3*ifac + 2] = iedg + (e+2)%n_edge_k;
              ifac++;
              if (ifac == *dn_face) break;
            }
            gfac++;
          }
        }
        if (ifac == *dn_face) break;
      }
      if (ifac == *dn_face) break;
      ivtx += n_vtx_k;
      iedg += n_edge_k;
    }
  }
  free (distrib_face);
}


static void
_get_connectivity
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **nFace,
 int         ***faceEdgeIdx,
 int         ***faceEdge,
 int         ***faceVtxIdx,
 int         ***faceVtx,
 PDM_g_num_t ***faceLNToGN,
 int          **nEdge,
 int         ***edgeVtxIdx,
 int         ***edgeVtx,
 int          **nVtx,
 double      ***vtxCoord,
 PDM_g_num_t ***vtxLNToGN
 )
{
  *nFace = (int *) malloc(sizeof(int) * n_part);
  *faceEdgeIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceEdge = (int **) malloc(sizeof(int *) * n_part);
  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *nEdge = (int *) malloc(sizeof(int) * n_part);
  *edgeVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *edgeVtx = (int **) malloc(sizeof(int *) * n_part);

  *nVtx = (int *) malloc(sizeof(int) * n_part);
  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);


  for (int ipart = 0; ipart < n_part; ipart++) {

    int _nFace;
    int _nEdge;
    int _nEdgePartBound;
    int _nVtx;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           ipart,
                           &_nFace,
                           &_nEdge,
                           &_nEdgePartBound,
                           &_nVtx,
                           &_nProc,
                           &_nTPart,
                           &_sFaceEdge,
                           &_sEdgeVtx,
                           &_sEdgeGroup,
                           &_nEdgeGroup2);

    int         *_faceTag;
    int         *_faceEdgeIdx;
    int         *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int         *_edgeTag;
    int         *_edgeFace;
    int         *_edgeVtxIdx;
    int         *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int         *_edgePartBoundProcIdx;
    int         *_edgePartBoundPartIdx;
    int         *_edgePartBound;
    int         *_vtxTag;
    double      *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int         *_edgeGroupIdx;
    int         *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &_faceTag,
                           &_faceEdgeIdx,
                           &_faceEdge,
                           &_faceLNToGN,
                           &_edgeTag,
                           &_edgeFace,
                           &_edgeVtxIdx,
                           &_edgeVtx,
                           &_edgeLNToGN,
                           &_edgePartBoundProcIdx,
                           &_edgePartBoundPartIdx,
                           &_edgePartBound,
                           &_vtxTag,
                           &_vtx,
                           &_vtxLNToGN,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    /*for (int i = 0; i < _nFace; i++) {
      printf("face ("PDM_FMT_G_NUM"), edges =", _faceLNToGN[i]);
      for (int j = _faceEdgeIdx[i]; j < _faceEdgeIdx[i+1]; j++) {
      printf(" ("PDM_FMT_G_NUM")", _edgeLNToGN[PDM_ABS(_faceEdge[j])-1]);
      }
      printf("\n");
      }*/


    /* Faces */
    (*nFace)[ipart] = _nFace;
    (*faceEdgeIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceEdge)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

    memcpy ((*faceEdgeIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceEdge)[ipart], _faceEdge, _sFaceEdge * sizeof(int));
    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    /* Edges */
    (*nEdge)[ipart] = _nEdge;
    (*edgeVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
    (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

    memcpy ((*edgeVtxIdx)[ipart], _edgeVtxIdx, (_nEdge + 1) * sizeof(int));
    memcpy ((*edgeVtx)[ipart], _edgeVtx, _sEdgeVtx * sizeof(int));

    /* Vertices */
    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));


    /* Compute face-vtx connectivity */
    int *_faceVtx = (*faceVtx)[ipart];

    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));

    for (int i = 0; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i];
      int ivtx2 = _edgeVtx[2*i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
    for (int i = 0; i < _nVtx; i++) {
      vtxEdgeN[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i] - 1;
      int ivtx2 = _edgeVtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
      vtxEdgeN[ivtx1] += 1;
      vtxEdgeN[ivtx2] += 1;
    }
    free(vtxEdgeN);

    for (int i = 0; i < _nFace; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i+1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _faceVtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edgeVtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edgeVtx[2*(edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtxEdgeIdx[vtx_cur - 1]; j <  vtxEdgeIdx[vtx_cur]; j++) {
          for (int k = 0; k < __nEdge; k++) {
            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edgeVtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edgeVtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edgeVtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          printf("error face ("PDM_FMT_G_NUM"), vtx tmp:\n", _faceLNToGN[i]);
          for (int l = 0; l < idxVtx; l++) {
            printf("  %d ("PDM_FMT_G_NUM")\n", _vertices[l], _vtxLNToGN[_vertices[l]-1]);
          }
          printf("\n");
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtxedge !!!!\n");
          abort();
        }
      }
      /*printf("face ("PDM_FMT_G_NUM"), vtx :", _faceLNToGN[i]);
        for (int l = 0; l < __nEdge; l++) {
        printf(" ("PDM_FMT_G_NUM")", _vtxLNToGN[_vertices[l]-1]);
        }
        printf("\n");*/
    }

    free (vtxEdge);
    free (vtxEdgeIdx);

  }
}




static void
_gen_src_mesh
(
 const int            active_rank,
 PDM_MPI_Comm         comm,
 int                  n_part,
 const PDM_g_num_t    nv,
 const double         length,
 int                **n_vtx,
 PDM_g_num_t       ***vtx_g_num,
 double            ***vtx_coord,
 int                **n_face,
 PDM_g_num_t       ***face_g_num,
 int               ***face_vtx_idx,
 int               ***face_vtx
 )
{
  if (active_rank) {
    /*
     *  Distributed mesh generation
     */
    PDM_g_num_t     ng_face;
    PDM_g_num_t     ng_vtx;
    PDM_g_num_t     ng_edge;
    int             dn_vtx;
    double         *dvtx_coord      = NULL;
    int             dn_face;
    int            *dface_vtx_idx   = NULL;
    PDM_g_num_t    *dface_vtx       = NULL;
    PDM_g_num_t    *dface_edge      = NULL;
    int             dn_edge;
    PDM_g_num_t    *dedge_vtx       = NULL;
    PDM_g_num_t    *dedge_face      = NULL;
    int             n_edge_group;
    int            *dedge_group_idx = NULL;
    PDM_g_num_t    *dedge_group     = NULL;

    double radius = 0.05 * length;

    _gen_rocket (comm,
                 nv,
                 length,
                 radius,
                 &ng_face,
                 &ng_vtx,
                 &ng_edge,
                 &dn_vtx,
                 &dvtx_coord,
                 &dn_face,
                 &dface_vtx_idx,
                 &dface_vtx,
                 &dface_edge,
                 &dn_edge,
                 &dedge_vtx,
                 &dedge_face,
                 &n_edge_group,
                 &dedge_group_idx,
                 &dedge_group);


    /*
     *  Create mesh partitions
     */
    int have_dface_part = 0;

    int *dface_part    = (int *) malloc (sizeof(int) * dn_face);
    int *dedge_vtx_idx = (int *) malloc (sizeof(int) * (dn_edge + 1));

    dedge_vtx_idx[0] = 0;
    for (int i = 0; i < dn_edge; i++) {
      dedge_vtx_idx[i+1] = 2 + dedge_vtx_idx[i];
    }


    /*
     *  Split mesh
     */
    PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

    int n_property_face = 0;
    int *renum_properties_face = NULL;
    int n_property_edge = 0;
    int *renum_properties_edge = NULL;

    PDM_part_t *ppart = PDM_part_create (comm,
                                         part_method,
                                         "PDM_PART_RENUM_CELL_NONE",
                                         "PDM_PART_RENUM_FACE_NONE",
                                         n_property_face,
                                         renum_properties_face,
                                         n_property_edge,
                                         renum_properties_edge,
                                         n_part,
                                         dn_face,
                                         dn_edge,
                                         dn_vtx,
                                         n_edge_group,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL,
                                         have_dface_part,
                                         dface_part,
                                         dedge_face,
                                         dedge_vtx_idx,
                                         dedge_vtx,
                                         NULL,
                                         dvtx_coord,
                                         NULL,
                                         dedge_group_idx,
                                         dedge_group);

    free (dface_part);
    free (dvtx_coord);
    free (dface_vtx_idx);
    free (dface_vtx);
    free (dface_edge);
    free (dedge_vtx_idx);
    free (dedge_vtx);
    free (dedge_face);
    free (dedge_group_idx);
    free (dedge_group);


    int **face_edge_idx = NULL;
    int **face_edge     = NULL;
    int  *n_edge        = NULL;
    int **edge_vtx_idx  = NULL;
    int **edge_vtx      = NULL;
    _get_connectivity (ppart,
                       n_part,
                       n_face,
                       &face_edge_idx,
                       &face_edge,
                       face_vtx_idx,
                       face_vtx,
                       face_g_num,
                       &n_edge,
                       &edge_vtx_idx,
                       &edge_vtx,
                       n_vtx,
                       vtx_coord,
                       vtx_g_num);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free (face_edge_idx[i_part]);
      free (face_edge[i_part]);
      free (edge_vtx_idx[i_part]);
      free (edge_vtx[i_part]);
    }
    free (n_edge);
    free (face_edge_idx);
    free (face_edge);
    free (edge_vtx_idx);
    free (edge_vtx);

    PDM_part_free (ppart);
  }

  else {
    *n_vtx        = malloc (sizeof(int)           * n_part);
    *vtx_g_num    = malloc (sizeof(PDM_g_num_t *) * n_part);
    *vtx_coord    = malloc (sizeof(double *)      * n_part);
    *n_face       = malloc (sizeof(int)           * n_part);
    *face_vtx_idx = malloc (sizeof(int *)         * n_part);
    *face_vtx     = malloc (sizeof(int *)         * n_part);
    *face_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      (*n_vtx)[i_part] = 0;
      (*vtx_g_num)[i_part] = malloc (sizeof(PDM_g_num_t) * (*n_vtx)[i_part]);
      (*vtx_coord)[i_part] = malloc (sizeof(double)      * (*n_vtx)[i_part] * 3);

      (*n_face)[i_part] = 0;
      (*face_g_num)[i_part] = malloc (sizeof(PDM_g_num_t) * (*n_face)[i_part]);
      (*face_vtx_idx)[i_part] = malloc (sizeof(int)       * ((*n_face)[i_part] + 1));
      (*face_vtx_idx)[i_part][0] = 0;
      (*face_vtx)[i_part] = malloc (sizeof(int) * (*face_vtx_idx)[i_part][(*n_face)[i_part]]);
    }
  }
}



static void
_write_point_cloud
(
 const char        *filename,
 const char        *header,
 const int          n_pts,
 const double       coord[],
 const PDM_g_num_t  g_num[],
 const double       field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (field != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS field double\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, "%f\n", field[i]);
    }
  }

  fclose(f);
}



static void
_write_polydata
(
 const char        *filename,
 const char        *header,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
    }
  }


  fclose(f);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  srand (time(NULL) + i_rank);


  /*
   *  Set default values
   */
  PDM_g_num_t nv              = 10;
  double      length          = 1.;
  double      domain_size     = 3.;
  PDM_g_num_t ng_tgt          = 10;
  int         grid            = 0;
  int         n_max_per_leaf  = 10;
  int         on_ground       = 0;
  int         n_proc_data_src = -1;
  int         post            = 0;
  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nv,
              &length,
              &domain_size,
              &ng_tgt,
              &grid,
              &n_max_per_leaf,
              &on_ground,
              &n_proc_data_src,
              &post);


  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  int          *n_tgt     = NULL;
  double      **tgt_coord = NULL;
  PDM_g_num_t **tgt_g_num = NULL;

  double xmin = -0.5*domain_size;
  double origin[3] = {xmin, xmin, xmin};
  if (grid) {
    _gen_cloud_grid (PDM_MPI_COMM_WORLD,
                     n_part_tgt,
                     ng_tgt,
                     origin,
                     domain_size,
                     &n_tgt,
                     &tgt_g_num,
                     &tgt_coord);
  }

  else {
    n_part_tgt = 1;
    n_tgt     = malloc (sizeof(int)           * n_part_tgt);
    tgt_g_num = malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
    tgt_coord = malloc (sizeof(double *)      * n_part_tgt);

    _gen_cloud_random (PDM_MPI_COMM_WORLD,
                       ng_tgt,
                       origin,
                       domain_size,
                       n_tgt,
                       tgt_g_num,
                       tgt_coord);
  }

  /*
   *  Define the source point cloud
   */
  PDM_MPI_Comm src_comm = PDM_MPI_COMM_WORLD;
  int active_rank_src = 1;
  if (n_proc_data_src > 0 && n_proc_data_src < n_rank) {
    int rank_in_node = PDM_io_mpi_node_rank (PDM_MPI_COMM_WORLD);

    int n_nodes = 0;
    int i_node = -1;
    int master_rank = 0;
    if (rank_in_node == 0) {
      master_rank = 1;
    }

    int *rank_in_nodes = malloc(sizeof(int) * n_rank);

    PDM_MPI_Allreduce (&master_rank, &n_nodes, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
    PDM_MPI_Allgather (&rank_in_node, 1, PDM_MPI_INT, rank_in_nodes, 1, PDM_MPI_INT, PDM_MPI_COMM_WORLD);

    active_rank_src = 0;

    for (int i = 0; i < i_rank; i++) {
      if (rank_in_nodes[i] == 0) {
        i_node += 1;
      }
    }

    if (n_proc_data_src <= n_nodes) {
      if (i_node < n_proc_data_src && rank_in_node == 0) {
        active_rank_src = 1;
      }
    }

    else {

      if (rank_in_node < (n_proc_data_src / n_nodes)) {
        active_rank_src = 1;
      }
      if ((rank_in_node == (n_proc_data_src / n_nodes)) &&
          (i_node < (n_proc_data_src % n_nodes))) {
        active_rank_src = 1;
      }

    }

    PDM_MPI_Comm_split (PDM_MPI_COMM_WORLD, active_rank_src, i_rank, &src_comm);

    free (rank_in_nodes);
  }


  int n_part_src = 1;
  int          *n_vtx     = NULL;
  double      **vtx_coord = NULL;
  PDM_g_num_t **vtx_g_num = NULL;
  int          *n_face       = NULL;
  int         **face_vtx_idx = NULL;
  int         **face_vtx     = NULL;
  PDM_g_num_t **face_g_num   = NULL;

  _gen_src_mesh (active_rank_src,
                 src_comm,
                 n_part_src,
                 nv,
                 length,
                 &n_vtx,
                 &vtx_g_num,
                 &vtx_coord,
                 &n_face,
                 &face_g_num,
                 &face_vtx_idx,
                 &face_vtx);


  PDM_g_num_t n_g_face_loc = 0;
  PDM_g_num_t n_g_vtx_loc = 0;

  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_vtx = 0;

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    for (int i = 0; i < n_face[i_part]; i++) {
      n_g_face_loc = PDM_MAX (n_g_face_loc, face_g_num[i_part][i]);
    }

    for (int i = 0; i < n_vtx[i_part]; i++) {
      n_g_vtx_loc = PDM_MAX (n_g_vtx_loc, vtx_g_num[i_part][i]);
    }
  }
  PDM_MPI_Allreduce (&n_g_face_loc, &n_g_face, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_MPI_Allreduce (&n_g_vtx_loc, &n_g_vtx, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);


  double z_shift;
  if (on_ground) {
    z_shift = 0.5*domain_size;
  } else {
    z_shift = 0.5*length;
  }

  for (int i_part = 0; i_part <  n_part_src; i_part++) {
    for (int i = 0; i < n_vtx[i_part]; i++) {
      vtx_coord[i_part][3*i+2] -= z_shift;
    }
  }


  if (post) {
    char filename[999];

    for (int i_part = 0; i_part < n_part_tgt; i_part++) {
      sprintf(filename, "tgt_%3.3d.vtk", n_part_tgt*i_rank + i_part);
      _write_point_cloud (filename,
                          "tgt",
                          n_tgt[i_part],
                          tgt_coord[i_part],
                          tgt_g_num[i_part],
                          NULL);
    }

    for (int i_part = 0; i_part < n_part_src; i_part++) {
      sprintf(filename, "src_mesh_%3.3d.vtk", n_part_src*i_rank + i_part);
      _write_polydata (filename,
                       "src_mesh",
                       n_vtx[i_part],
                       vtx_coord[i_part],
                       vtx_g_num[i_part],
                       n_face[i_part],
                       face_vtx_idx[i_part],
                       face_vtx[i_part],
                       face_g_num[i_part]);
    }
  }



  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t *id_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                               n_point_cloud,
                                                               PDM_MPI_COMM_WORLD,
                                                               PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (id_dist,
                                                 n_part_src);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_dist_cloud_surf_surf_mesh_part_set (id_dist,
                                            i_part,
                                            n_face[i_part],
                                            face_vtx_idx[i_part],
                                            face_vtx[i_part],
                                            face_g_num[i_part],
                                            n_vtx[i_part],
                                            vtx_coord[i_part],
                                            vtx_g_num[i_part]);
  }


  PDM_dist_cloud_surf_n_part_cloud_set (id_dist,
                                        0,
                                        n_part_tgt);

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    PDM_dist_cloud_surf_cloud_set (id_dist,
                                   0,
                                   i_part,
                                   n_tgt[i_part],
                                   tgt_coord[i_part],
                                   tgt_g_num[i_part]);
  }

  /* Compute distance */
  // PDM_dist_cloud_surf_compute (id_dist);
  PDM_dist_cloud_surf_compute (id_dist);

  PDM_dist_cloud_surf_dump_times(id_dist);



  PDM_dist_cloud_surf_free (id_dist);

  /*
   *  Finalize
   */
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    free (tgt_coord[i_part]);
    free (tgt_g_num[i_part]);
  }
  free (n_tgt);
  free (tgt_coord);
  free (tgt_g_num);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    free (vtx_coord[i_part]);
    free (vtx_g_num[i_part]);
    free (face_vtx_idx[i_part]);
    free (face_vtx[i_part]);
    free (face_g_num[i_part]);
  }
  free (n_vtx);
  free (vtx_coord);
  free (vtx_g_num);
  free (n_face);
  free (face_vtx_idx);
  free (face_vtx);
  free (face_g_num);


  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
