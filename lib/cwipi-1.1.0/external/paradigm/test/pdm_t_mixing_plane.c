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
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_mesh_intersection.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_part_connectivity_transform.h"

#include "pdm_triangulate.h"
#include "pdm_binary_search.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA     <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -nA     <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -n_part <level>  Number vtx in side of mesh B (default : 10).\n\n"
     "  -t               Element kind .\n\n"
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
 PDM_g_num_t           *n_vtx_seg,
 int                   *n_part,
 PDM_Mesh_nodal_elt_t  *elt_type,
 int                   *n_band,
 double                *ratio,
 double                *length,
 int                   *verbose
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-n_band") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_band = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ratio") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *ratio = atof(argv[i]);
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
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static void
_add_depth
(
 double *coord,
 double  length
 )
{
  double r0 = 0.3;
  double r1 = 0.6;

  double x = coord[0]/length;
  double y = coord[1]/length;
  // double z = coord[2];
  //y = 0.5*(1 - cos(PDM_PI*y));
  y = y*y*(3 - 2*y);

  if (0) {
    // angular sector
    double t = PDM_PI * (0.5 + (x - 0.5 + 0.3*cos(3*y)) / 6.);
    double r = r0 + (r1 - r0) * y;

    coord[0] = length*0.02*sin(5*y);
    coord[1] = length*r * cos(t);
    coord[2] = length*r * sin(t);
  }
  else {
    coord[0] = length*x;
    coord[1] = length*y;
  }
}


static void
_generate_mesh
(
 const PDM_MPI_Comm           comm,
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    rotate,
 const double                 xmin,
 const double                 ymin,
 const double                 zmin,
 const double                 length,
 const PDM_split_dual_t       part_method,
 const int                    n_part,
       PDM_multipart_t      **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  double* vtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

  // randomize
  if (1) {
    double noise = 0.2*length/(double) (n_vtx_seg - 1);
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      if (PDM_ABS(vtx_coord[3*i_vtx  ] - xmin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx  ] - xmin - length) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin         ) > 1.e-9 &&
          PDM_ABS(vtx_coord[3*i_vtx+1] - ymin - length) > 1.e-9) {
        srand(distrib_vtx[i_rank] + i_vtx);
        for (int i = 0; i < 2; i++) {
          vtx_coord[3*i_vtx+i] += noise*0.5*(2*rand()/(double) RAND_MAX - 1);
        }
      }
    }
  }

  if (rotate) {
    // Do something...
  }

  if (1) {
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      _add_depth(&vtx_coord[3*i_vtx], length);
    }
  }

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "mixing_plane_mesh_");
  }




  int n_domain = 1;

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  *_mpart = mpart;

  PDM_DMesh_nodal_free(dmn);
}


static inline double
_eval_field
(
 const double x,
 const double y,
 const double z
 )
{
  PDM_UNUSED(x);
  PDM_UNUSED(z);
  return y;
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static void
_define_levels
(
 const int     n_level,
       double  val_min,
       double  val_max,
       double  ratio,
       double *levels
)
{
  levels[0] = -HUGE_VAL;

  double h0;
  double idenom;
  if (ratio == 1) {
    h0 = 1./(double) (n_level - 1);
  }
  else {
    int m = (n_level - 1) / 2;
    if (n_level%2 == 0) {
      h0 = (1 - ratio)/(2 - pow(ratio, m)*(ratio + 1));
    }
    else {
      h0 = 0.5*(1 - ratio)/(1 - pow(ratio, m));
    }
    idenom = h0/(1 - ratio);
  }

  h0 *= val_max - val_min;

  for (int i = 1; i < n_level-1; i++) {
    if (ratio == 1) {
      levels[i] = val_min + i*h0;
    }
    else {
      if (i <= (n_level-1)/2) {
        levels[i] = val_min + (1 - pow(ratio, i)) * idenom;
      }
      else {
        levels[i] = val_max - (1 - pow(ratio, n_level-i-1)) * idenom;
      }
    }
  }

  levels[n_level-1] = HUGE_VAL;
}
PDM_GCC_SUPPRESS_WARNING_POP


/* Linked list node for convex polygon clipping */

struct _ll_node_t {

  struct _ll_node_t *next;
  double             u;
  double             v;
  double             f;

};

typedef struct _ll_node_t _ll_node_t;


static inline double
_polygon_area
(
 _ll_node_t *origin
 )
{
  double area = 0.;

  _ll_node_t *i = origin->next;
  _ll_node_t *j = i->next;

  while (j != origin) {
    area +=
    (i->u - origin->u)*(j->v - origin->v) -
    (i->v - origin->v)*(j->u - origin->u);
    i = i->next;
    j = j->next;
  }

  return area;
}


static _ll_node_t *
_clip_convex_polygon
(
       _ll_node_t *origin,
       _ll_node_t *link_below0,
       _ll_node_t *link_below1,
       _ll_node_t *link_above0,
       _ll_node_t *link_above1,
 const double      value
)
{
  _ll_node_t *prev_below0 = NULL;
  _ll_node_t *prev_above1 = NULL;

  _ll_node_t *link_below[2] = {link_below0, link_below1};
  _ll_node_t *link_above[2] = {link_above0, link_above1};

  for (int i = 0; i < 2; i++) {
    link_below[i]->next = NULL;
    link_above[i]->next = NULL;
  }

  int idx = 0;

  _ll_node_t *current = origin;
  while (1) {
    _ll_node_t *next = current->next;

    double f_current = current->f - value;
    double f_next    = next->f    - value;

    if (PDM_SIGN(f_current) != PDM_SIGN(f_next)) {

      double t = f_current/(f_current - f_next);
      double u = (1-t)*current->u + t*next->u;
      double v = (1-t)*current->v + t*next->v;

      if (idx == 0) {
        prev_below0 = current;
        link_below[0]->next = link_below[1];
        link_above[0]->next = next;
      }
      else {
        link_below[1]->next = next;
        link_above[1]->next = link_above[0];
        prev_above1 = current;
      }

      link_below[idx]->u = u;
      link_below[idx]->v = v;
      link_below[idx]->f = value;

      link_above[idx]->u = u;
      link_above[idx]->v = v;
      link_above[idx]->f = value;

      idx++;
      if (idx == 2) {
        break;
      }
    }

    // Move to next node
    current = current->next;
    if (current == origin) {
      break;
    }
  }
  assert(idx == 2);

  prev_below0->next = link_below[0];
  prev_above1->next = link_above[1];

  return link_above[0];
}



static void
_get_polygon
(
 double      tri0[3],
 double      tri1[3],
 double      tri2[3],
 _ll_node_t *origin,
 int        *n_vtx,
 double     *vtx_coord
 )
 {
  int idx = 0;
  _ll_node_t *current = origin;
  while (1) {
    double u = current->u;
    double v = current->v;
    double w = 1 - u - v;
    for (int i = 0; i < 3; i++) {
      vtx_coord[3*idx+i] = w*tri0[i] + u*tri1[i] + v*tri2[i];
    }
    idx++;

    current = current->next;
    if (current == origin) {
      break;
    }
  }

  *n_vtx = idx;
 }

 static void
 _export_triangle_band
 (
  const int         i_rank,
  const int         ipart,
  const int         face_id,
  const int         tria_id,
  const int         band_id,
        double      tri0[3],
        double      tri1[3],
        double      tri2[3],
        double     *dbg_coord,
        _ll_node_t *origin
  )
 {
  int _n_vtx = 0;
  _get_polygon(tri0,
               tri1,
               tri2,
               origin,
               &_n_vtx,
               dbg_coord);
  // log_trace("between levels %d and %d:\n", level_id, level_id+1);
  for (int i = 0; i < _n_vtx; i++) {
    log_trace("  %3.3f %3.3f %3.3f\n",
              dbg_coord[3*i  ],
              dbg_coord[3*i+1],
              dbg_coord[3*i+2]);
  }

  int connec_idx[2] = {0, _n_vtx};
  int *connec = malloc(sizeof(int) * _n_vtx);
  for (int i = 0; i < _n_vtx; i++) {
    connec[i] = i+1;
  }

  const char *field_name[]   = {"face_id", "band_id"};
  const int  *field_value[2] = {&face_id, &band_id};

  char filename[999];
  sprintf(filename, "triangles/rank_%d_part%d_face%d_tria%d_band%d.vtk",
          i_rank, ipart, face_id, tria_id, band_id);
  PDM_vtk_write_polydata_with_field(filename,
                                    _n_vtx,
                                    dbg_coord,
                                    NULL,
                                    1,
                                    connec_idx,
                                    connec,
                                    NULL,
                                    NULL,
                                    2,
                                    field_name,
                                    field_value);
  free(connec);
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static void
_mixing_plane
(
       PDM_MPI_Comm    comm,          // only for debug
 const int             n_band,
       double         *levels,
 const int             n_part,
       int            *n_vtx,
       double        **vtx_coord,
       double        **vtx_field,
       PDM_g_num_t   **vtx_ln_to_gn,  // only for debug
       int            *n_face,
       int           **face_vtx_idx,
       int           **face_vtx,
       PDM_g_num_t   **face_ln_to_gn, // only for debug
       int          ***face_band_idx,
       int          ***face_band,
       double       ***face_band_area
 )
{
  int dbg  = 0;
  int visu = 0;

  double *dbg_coord = NULL;
  if (visu) {
    dbg_coord = malloc(sizeof(double) * 120); // size?
  }

  /* Prepare for triangulation */
  int max_face_vtx_n = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int iface = 0; iface < n_face[ipart]; iface++) {
      int face_vtx_n = face_vtx_idx[ipart][iface+1] - face_vtx_idx[ipart][iface];
      max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    }
  }

  PDM_triangulate_state_t *tri_state = NULL;
  int *tri_vtx = NULL;

  if (max_face_vtx_n > 4) {
    tri_state = PDM_triangulate_state_create(max_face_vtx_n);
  }

  tri_vtx = malloc(sizeof(int) * (max_face_vtx_n - 2) * 3);


  *face_band_idx  = malloc(sizeof(int    *) * n_part);
  *face_band      = malloc(sizeof(int    *) * n_part);
  *face_band_area = malloc(sizeof(double *) * n_part);

  int max_face_band_n = 3;
  _ll_node_t *nodes = malloc(sizeof(_ll_node_t) * (3+4*(max_face_band_n-1)));

  /* Intial nodes : triangle vertices */
  nodes[0].u = 0.;
  nodes[0].v = 0.;

  nodes[1].u = 1.;
  nodes[1].v = 0.;

  nodes[2].u = 0.;
  nodes[2].v = 1.;

  for (int ipart = 0; ipart < n_part; ipart++) {

    if (dbg) {
      int i_rank;
      PDM_MPI_Comm_rank(comm, &i_rank);

      char filename[999];
      sprintf(filename, "mixing_plane_part_%d_%2.2d.vtk", ipart, i_rank);
      PDM_vtk_write_polydata_field(filename,
                                   n_vtx        [ipart],
                                   vtx_coord    [ipart],
                                   vtx_ln_to_gn [ipart],
                                   n_face       [ipart],
                                   face_vtx_idx [ipart],
                                   face_vtx     [ipart],
                                   face_ln_to_gn[ipart],
                                   NULL,
                                   NULL,
                                   "field",
                                   vtx_field[ipart]);
    }

    int s_face_band = 4*n_face[ipart];
    (*face_band_idx )[ipart]  = malloc(sizeof(int) * (n_face[ipart] + 1));
    (*face_band     )[ipart] = malloc(sizeof(int   ) * s_face_band);
    (*face_band_area)[ipart] = malloc(sizeof(double) * s_face_band);

    int    *_face_band_idx  = (*face_band_idx )[ipart];
    int    *_face_band      = (*face_band     )[ipart];
    double *_face_band_area = (*face_band_area)[ipart];

    _face_band_idx[0] = 0;

    /* Loop on faces */
    for (int face_id = 0; face_id < n_face[ipart]; face_id++) {

      if (dbg) {
        log_trace("face_id = %d, idx = %d\n",
                  face_id,
                  _face_band_idx[face_id]);
      }

      int *_face_vtx = face_vtx[ipart] + face_vtx_idx[ipart][face_id];
      int face_vtx_n = face_vtx_idx[ipart][face_id+1] - face_vtx_idx[ipart][face_id];

      /* Compute min/max field values in current face */
      double face_min_field =  HUGE_VAL;
      double face_max_field = -HUGE_VAL;
      for (int ivtx = 0; ivtx < face_vtx_n; ivtx++) {
        int vtx_id = _face_vtx[ivtx] - 1;
        face_min_field = PDM_MIN(face_min_field, vtx_field[ipart][vtx_id]);
        face_max_field = PDM_MAX(face_max_field, vtx_field[ipart][vtx_id]);
      }

      if (dbg) {
        log_trace("  %f <= field <= %f\n", face_min_field, face_max_field);
      }

      /* Get ids of bands crossed by current face */
      int face_band_min = PDM_binary_search_gap_double(face_min_field,
                                                  levels,
                                                  n_band + 1);
      assert(face_band_min >= 0);

      int face_band_max = PDM_binary_search_gap_double(face_max_field,
                                                  levels + face_band_min,
                                                  n_band + 1 - face_band_min);
      assert(face_band_max >= 0);
      face_band_max += face_band_min;

      if (face_max_field == levels[face_band_max]) {
        face_band_max--;
      }

      if (dbg) {
        log_trace("  %d <= band <= %d\n", face_band_min, face_band_max);
      }

      /* Number of bands crossed by current face */
      int face_band_n = face_band_max - face_band_min + 1;
      _face_band_idx[face_id+1] = _face_band_idx[face_id] + face_band_n;

      /* Realloc if necessary */
      if (_face_band_idx[face_id+1] >= s_face_band) {
        s_face_band = PDM_MAX(2*s_face_band,
                              _face_band_idx[face_id+1]);
        (*face_band     )[ipart] = realloc((*face_band     )[ipart], sizeof(int   ) * s_face_band);
        (*face_band_area)[ipart] = realloc((*face_band_area)[ipart], sizeof(double) * s_face_band);
        _face_band      = (*face_band     )[ipart];
        _face_band_area = (*face_band_area)[ipart];
      }

      int    *__face_band      = _face_band      + _face_band_idx[face_id];
      double *__face_band_area = _face_band_area + _face_band_idx[face_id];

      if (face_band_n > max_face_band_n) {
        max_face_band_n = PDM_MAX(2*max_face_band_n,
                                  face_band_n);
        nodes = realloc(nodes,
                        sizeof(_ll_node_t) * (3+4*(max_face_band_n-1)));
      }

      /* Initialize face-band areas */
      for (int iband = 0; iband < face_band_n; iband++) {
        __face_band     [iband] = face_band_min + iband;
        __face_band_area[iband] = 0.;
      }



      /* Triangulate current face */
      int n_tri;
      if (face_vtx_n == 3) {
        /* Triangular face */
        n_tri = 1;
        memcpy(tri_vtx, _face_vtx, sizeof(int) * 3);
      }
      else if (face_vtx_n == 4) {
        /* Quadrilateral face */
        n_tri = PDM_triangulate_quadrangle(3,
                                           vtx_coord[ipart],
                                           NULL,
                                           _face_vtx,
                                           tri_vtx);
      }
      else {
        /* Polygonal face */
        n_tri = PDM_triangulate_polygon(3,
                                        face_vtx_n,
                                        vtx_coord[ipart],
                                        NULL,
                                        _face_vtx,
                                        PDM_TRIANGULATE_MESH_DEF,
                                        tri_vtx,
                                        tri_state);
      }


      /* Loop on subtriangles */
      for (int itri = 0; itri < n_tri; itri++) {

        int *_tri_vtx = tri_vtx + 3*itri;
        int vtx_id0 = _tri_vtx[0] - 1;
        int vtx_id1 = _tri_vtx[1] - 1;
        int vtx_id2 = _tri_vtx[2] - 1;

        if (dbg) {
          log_trace("    tri %d: %d %d %d\n",
                    itri, vtx_id0, vtx_id1, vtx_id2);
        }

        /* Reset linked list */
        for (int i = 0; i < 3; i++) {
          nodes[i].next = &nodes[(i+1)%3];
          nodes[i].f    = vtx_field[ipart][_tri_vtx[i]-1];
        }

        /* Compute surface area of current triangle */
        double vec01[3];
        double vec02[3];
        for (int i = 0; i < 3; i++) {
          vec01[i] = vtx_coord[ipart][3*vtx_id1+i] - vtx_coord[ipart][3*vtx_id0+i];
          vec02[i] = vtx_coord[ipart][3*vtx_id2+i] - vtx_coord[ipart][3*vtx_id0+i];
        }
        double normal[3];
        PDM_CROSS_PRODUCT(normal, vec01, vec02);

        double tria_area = 0.5*PDM_MODULE(normal);

        /* Get min/max field values in current triangle */
        double tria_min_field =  HUGE_VAL;
        double tria_max_field = -HUGE_VAL;
        int imin = -1;
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          if (nodes[ivtx].f < tria_min_field) {
            imin = ivtx;
            tria_min_field = nodes[ivtx].f;
          }
          if (nodes[ivtx].f > tria_max_field) {
            tria_max_field = nodes[ivtx].f;
          }
        }

        if (dbg) {
          log_trace("    %f <= field <= %f\n", tria_min_field, tria_max_field);
        }

        /* Get ids of bands crossed by current triangle */
        /* (tria_band ids must be offset by face_band_min) */
        int tria_band_min = PDM_binary_search_gap_double(tria_min_field,
                                                         levels + face_band_min,
                                                         face_band_n + 1);
        assert(tria_band_min >= 0);

        int tria_band_max = PDM_binary_search_gap_double(tria_max_field,
                                                         levels + face_band_min + tria_band_min,
                                                         face_band_n + 1 - tria_band_min);
        assert(tria_band_max >= 0);

        tria_band_max += tria_band_min;

        if (tria_max_field == levels[tria_band_max]) {
          tria_band_max--;
        }

        if (dbg) {
          log_trace("    %d <= band <= %d\n",
                    face_band_min + tria_band_min,
                    face_band_min + tria_band_max);
        }

        /* Number of bands crossed by current triangle */
        int tria_band_n = tria_band_max - tria_band_min + 1;

        if (tria_band_n == 1) {

          /* Current triangle crosses a single band */
          __face_band_area[tria_band_min] += tria_area;

        }
        else {

          /* Current triangle crosses multiple bands */
          _ll_node_t *below = &nodes[imin];

          int idx = 0;
          for (int iband = tria_band_min; iband < tria_band_max; iband++) {
            int level_id = face_band_min + iband + 1;

            _ll_node_t *above = _clip_convex_polygon(below,
                                                     &nodes[3+4*idx],
                                                     &nodes[3+4*idx+1],
                                                     &nodes[3+4*idx+2],
                                                     &nodes[3+4*idx+3],
                                                     levels[level_id]);
            idx++;

            if (visu) {
              int i_rank;
              PDM_MPI_Comm_rank(comm, &i_rank);
              _export_triangle_band(i_rank,
                                    ipart,
                                    face_id,
                                    itri,
                                    face_band_min + iband,
                                    vtx_coord[ipart] + 3*vtx_id0,
                                    vtx_coord[ipart] + 3*vtx_id1,
                                    vtx_coord[ipart] + 3*vtx_id2,
                                    dbg_coord,
                                    below);
            }

            __face_band_area[iband] += _polygon_area(below)*tria_area;

            below = above;
          }

          /* Remaining polygon goes to last triangle band */
          __face_band_area[tria_band_max] += _polygon_area(below)*tria_area;

          if (visu) {
            int i_rank;
            PDM_MPI_Comm_rank(comm, &i_rank);
            _export_triangle_band(i_rank,
                                  ipart,
                                  face_id,
                                  itri,
                                  face_band_min + tria_band_max,
                                  vtx_coord[ipart] + 3*vtx_id0,
                                  vtx_coord[ipart] + 3*vtx_id1,
                                  vtx_coord[ipart] + 3*vtx_id2,
                                  dbg_coord,
                                  below);
          }

        }

      } // End of loop on current face's subtriangles


      if (dbg) {
        for (int iband = _face_band_idx[face_id]; iband < _face_band_idx[face_id+1]; iband++) {
          log_trace("face %d, band %d, area = %f\n",
                    face_id,
                    _face_band[iband],
                    _face_band_area[iband]);
        }
      }


    } // End of loop on current part's faces

    s_face_band = _face_band_idx[n_face[ipart]];
    (*face_band     )[ipart] = realloc((*face_band     )[ipart], sizeof(int   ) * s_face_band);
    (*face_band_area)[ipart] = realloc((*face_band_area)[ipart], sizeof(double) * s_face_band);

  } // End of loop on parts


  /* Free memory */
  if (tri_state != NULL) {
    PDM_triangulate_state_destroy(tri_state);
  }
  free(tri_vtx);
  free(nodes);

  if (visu) {
    free(dbg_coord);
  }
}
PDM_GCC_SUPPRESS_WARNING_POP



/*============================================================================
 * Public function definitions
 *============================================================================*/

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
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);


  /* Read command line arguments */
  PDM_g_num_t          n_vtx_seg = 10;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TRIA3;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  int    n_part = 1;
  int    n_band = 5;
  double ratio  = 1.;
  double length = 1.;
  int    verbose = 0;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_part,
             &elt_type,
             &n_band,
             &ratio,
             &length,
             &verbose);

  /* Generate mesh */
  PDM_multipart_t *mpart = NULL;
  _generate_mesh(comm,
                 n_vtx_seg,
                 elt_type,
                 0,
                 0., 0., 0.,
                 length,
                 part_method,
                 n_part,
                 &mpart);

  /* Define bands */
  double *levels = malloc(sizeof(double) * (n_band + 1));
  _define_levels(n_band+1,
                 0,     //val_min,
                 length,//val_max,
                 ratio,
                 levels);

  if (verbose) {
    PDM_log_trace_array_double(levels, n_band+1, "levels : ");
  }



  /* Mixing plane */
  int          *n_vtx         = malloc(sizeof(int          ) * n_part);
  double      **vtx_coord     = malloc(sizeof(double      *) * n_part);
  double      **vtx_field     = malloc(sizeof(double      *) * n_part);
  int          *n_face        = malloc(sizeof(int          ) * n_part);
  int         **face_vtx_idx  = malloc(sizeof(int         *) * n_part);
  int         **face_vtx      = malloc(sizeof(int         *) * n_part);
  PDM_g_num_t **vtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **face_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int *face_edge_idx;
    int *face_edge;
    n_face[ipart] = PDM_multipart_part_connectivity_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge_idx,
                                                        &face_edge,
                                                        PDM_OWNERSHIP_KEEP);
    face_vtx_idx[ipart] = face_edge_idx;

    int *edge_vtx_idx;
    int *edge_vtx;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &edge_vtx,
                                        PDM_OWNERSHIP_KEEP);

    PDM_compute_face_vtx_from_face_and_edge(n_face[ipart],
                                            face_edge_idx,
                                            face_edge,
                                            edge_vtx,
                                            &(face_vtx[ipart]));

    n_vtx[ipart] = PDM_multipart_part_vtx_coord_get(mpart,
                                                    0,
                                                    ipart,
                                                    &(vtx_coord[ipart]),
                                                    PDM_OWNERSHIP_KEEP);


    vtx_field[ipart] = malloc(sizeof(double) * n_vtx[ipart]);
    for (int i = 0; i < n_vtx[ipart]; i++) {
      vtx_field[ipart][i] = _eval_field(vtx_coord[ipart][3*i  ],
                                        vtx_coord[ipart][3*i+1],
                                        vtx_coord[ipart][3*i+2]);
    }

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    ipart,
                                    PDM_MESH_ENTITY_VTX,
                                    &(vtx_ln_to_gn[ipart]),
                                    PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    ipart,
                                    PDM_MESH_ENTITY_FACE,
                                    &(face_ln_to_gn[ipart]),
                                    PDM_OWNERSHIP_KEEP);
  }


  int    **face_band_idx  = NULL;
  int    **face_band      = NULL;
  double **face_band_area = NULL;
  _mixing_plane(comm,
                n_band,
                levels,
                n_part,
                n_vtx,
                vtx_coord,
                vtx_field,
                vtx_ln_to_gn,
                n_face,
                face_vtx_idx,
                face_vtx,
                face_ln_to_gn,
                &face_band_idx,
                &face_band,
                &face_band_area);




  double *local_band_area = malloc(sizeof(double) * n_band);
  for (int band_id = 0; band_id < n_band; band_id++) {
    local_band_area[band_id] = 0;
  }


  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int face_id = 0; face_id < n_face[ipart]; face_id++) {
      for (int iband = face_band_idx[ipart][face_id]; iband < face_band_idx[ipart][face_id+1]; iband++) {
        int band_id = face_band[ipart][iband];
        local_band_area[band_id] += face_band_area[ipart][iband];
      }
    }
  }


  double *global_band_area = malloc(sizeof(double) * n_band);
  PDM_MPI_Allreduce(local_band_area, global_band_area, n_band,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);
  free(local_band_area);

  for (int band_id = 0; band_id < n_band; band_id++) {
    double hi = PDM_MIN(levels[band_id+1], length);
    double lo = PDM_MAX(levels[band_id], 0);
    double exact_area = (hi - lo) * length;
    double err = PDM_ABS(global_band_area[band_id] - exact_area);
    if (verbose) {
      log_trace("band %d [%f, %f]: %f / %f, err abs = %e, rel = %e\n",
                band_id, lo, hi,
                global_band_area[band_id], exact_area,
                err, err/exact_area);
    }
  }
  free(global_band_area);

  /* Free memory */
  free(levels);
  PDM_multipart_free(mpart);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(face_band_idx [ipart]);
    free(face_band     [ipart]);
    free(face_band_area[ipart]);

    free(vtx_field[ipart]);
    free(face_vtx [ipart]);
  }
  free(face_band_idx );
  free(face_band     );
  free(face_band_area);

  free(n_vtx);
  free(vtx_coord);
  free(vtx_field);
  free(n_face);
  free(face_vtx_idx);
  free(face_vtx);
  free(vtx_ln_to_gn);
  free(face_ln_to_gn);

  PDM_MPI_Barrier(comm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize ();

  return 0;
}
