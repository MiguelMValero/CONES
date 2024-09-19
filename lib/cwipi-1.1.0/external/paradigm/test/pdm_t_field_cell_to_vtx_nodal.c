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
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_extension.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_vtk.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_logging.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"
#include "pdm_domain_interface_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_field_cell_to_vtx.h"

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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *n_dom_i,
           int           *n_dom_j,
           int           *n_dom_k,
           int           *periodic_i,
           int           *periodic_j,
           int           *periodic_k,
           int           *t_elt,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ni") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_i = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nj") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_j = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nk") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_k = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else if (strcmp(argv[i], "-pj") == 0) {
      *periodic_j = 1;
    }
    else if (strcmp(argv[i], "-pk") == 0) {
      *periodic_k = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}
static void
_compute_face_vtx
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  int dbg = 0;

  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  for(int i_face = 0; i_face <  face_edge_idx[n_face]; ++i_face) {
    (*face_vtx)[i_face] = 10000;
  }

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *edge_tag = PDM_array_zeros_int(n_edge);

  for (int iface = 0; iface < n_face; iface++) {
    int *_face_vtx  = *face_vtx  + face_edge_idx[iface];
    int *_face_edge =  face_edge + face_edge_idx[iface];

    if (dbg) {
      log_trace("\nFace %d\n", iface);
      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }

    int _n_edge = face_edge_idx[iface+1] - face_edge_idx[iface];
    // first edge
    int iedge = PDM_ABS(_face_edge[0]) - 1;
    edge_tag[iedge] = 1;
    _face_vtx[0] = PDM_ABS(edge_vtx[2*iedge  ]);
    _face_vtx[1] = PDM_ABS(edge_vtx[2*iedge+1]);

    for (int i = 2; i < _n_edge; i++) {

      for (int j = 1; j < _n_edge; j++) {
        iedge = PDM_ABS(_face_edge[j]) - 1;

        if (edge_tag[iedge]) {
          continue;
        }

        if (edge_vtx[2*iedge] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge+1]);
          edge_tag[iedge] = 1;
          break;
        }
        else if (edge_vtx[2*iedge+1] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge]);
          edge_tag[iedge] = 1;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[iface+1] - face_edge_idx[iface]; ivtx++) {
        log_trace("%d ", _face_vtx[ivtx]);
      }
      log_trace("\n");
    }

    // reset tags
    for (int i = 0; i < _n_edge; i++) {
      iedge = PDM_ABS(_face_edge[i]) - 1;
      edge_tag[iedge] = 0;
    }
  }
  free(edge_tag);
}


static
void
_cell_center_3d
(
  int      pn_cell,
  int     *pcell_face_idx,
  int     *pcell_face,
  int     *pface_edge_idx,
  int     *pface_edge,
  int     *pface_vtx_idx,
  int     *pface_vtx,
  int     *pedge_vtx,
  double  *pvtx_coord,
  double **cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  if(pface_edge     != NULL) {
    from_edge = 1;
  }
  if(pface_vtx     != NULL) {
    from_face = 1;
  }
  assert(pvtx_coord     != NULL);

  double* entity_center = malloc(3 * pn_cell * sizeof(double ));

  if(from_face == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double) pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
          int i_vtx = pface_vtx[idx_vtx]-1;
          fcx += pvtx_coord[3*i_vtx  ];
          fcy += pvtx_coord[3*i_vtx+1];
          fcz += pvtx_coord[3*i_vtx+2];
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  } else if( from_edge == 1) {
    for(int i_cell = 0; i_cell < pn_cell; ++i_cell) {

      entity_center[3*i_cell  ] = 0.;
      entity_center[3*i_cell+1] = 0.;
      entity_center[3*i_cell+2] = 0.;

      double inv = 1./((double)  pcell_face_idx[i_cell+1] - pcell_face_idx[i_cell]);

      double fcx = 0;
      double fcy = 0;
      double fcz = 0;
      for(int idx_face = pcell_face_idx[i_cell]; idx_face < pcell_face_idx[i_cell+1]; ++idx_face) {
        int i_face = PDM_ABS(pcell_face[idx_face])-1;

        double inv2 = 1./((double)  pface_edge_idx[i_face+1] - pface_edge_idx[i_face]);

        for(int idx_edge = pface_edge_idx[i_face]; idx_edge < pface_edge_idx[i_face+1]; ++idx_edge) {
          int i_edge = PDM_ABS(pface_edge[idx_edge])-1;
          int i_vtx1 = pedge_vtx[2*i_edge  ] - 1;
          int i_vtx2 = pedge_vtx[2*i_edge+1] - 1;
          fcx += 0.5 * (pvtx_coord[3*i_vtx1  ] + pvtx_coord[3*i_vtx2  ]);
          fcy += 0.5 * (pvtx_coord[3*i_vtx1+1] + pvtx_coord[3*i_vtx2+1]);
          fcz += 0.5 * (pvtx_coord[3*i_vtx1+2] + pvtx_coord[3*i_vtx2+2]);
        }
        fcx = fcx * inv2;
        fcy = fcy * inv2;
        fcz = fcz * inv2;

        entity_center[3*i_cell  ] += fcx;
        entity_center[3*i_cell+1] += fcy;
        entity_center[3*i_cell+2] += fcz;
      }

      entity_center[3*i_cell  ] = entity_center[3*i_cell  ] * inv;
      entity_center[3*i_cell+1] = entity_center[3*i_cell+1] * inv;
      entity_center[3*i_cell+2] = entity_center[3*i_cell+2] * inv;
    } /* End cell */
  }

  *cell_center = entity_center;
}

/**
 *
 * \brief  Main
 *
 */
int main
(
 int   argc,
 char *argv[]
)
{

  /*
   *  Set default values
   */
  PDM_g_num_t          nx         = 10;
  PDM_g_num_t          ny         = 10;
  PDM_g_num_t          nz         = 10;
  int                  n_dom_i    = 1;
  int                  n_dom_j    = 1;
  int                  n_dom_k    = 1;
  int                  periodic_i = 0;
  int                  periodic_j = 0;
  int                  periodic_k = 0;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  post       = 0;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXA8;

  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &n_dom_i,
             &n_dom_j,
             &n_dom_k,
             &periodic_i,
             &periodic_j,
             &periodic_k,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Initialize structs
   */
  int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

  if (dim == 2) {
    n_dom_k    = 1;
    nz         = 1;
    periodic_k = 0;
  }

  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  // PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *) * n_domain);
  PDM_dcube_nodal_t **dcube = NULL;
  PDM_domain_interface_t *dom_intrf = NULL;

  const int order = 1;

  PDM_dcube_nodal_cart_topo(comm,
                            n_dom_i,
                            n_dom_j,
                            n_dom_k,
                            periodic_i,
                            periodic_j,
                            periodic_k,
                            nx,
                            ny,
                            nz,
                            length,
                            0.,
                            0.,
                            0.,
                            t_elt,
                            order,
                            &dcube,
                            &dom_intrf,
                            PDM_OWNERSHIP_KEEP);

  /*
   * Partitionnement
   */
  int* n_part_by_domain = (int *) malloc(n_domain * sizeof(int));
  for(int i = 0; i < n_domain; ++i) {
    n_part_by_domain[i] = n_part;
  }
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                   n_part_by_domain,
                                                   PDM_FALSE,
                                                   method,
                                                   PDM_PART_SIZE_HOMOGENEOUS,
                                                   NULL,
                                                   comm,
                                                   PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < n_domain; i++) {
    dmn[i] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i]);
    PDM_multipart_dmesh_nodal_set(mpart, i, dmn[i]);
  }

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t **pmn = malloc(n_domain * sizeof(PDM_part_mesh_nodal_t *));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      i_domain,
                                      &pmn[i_domain],
                                      PDM_OWNERSHIP_KEEP);
    if (post) {
      PDM_part_mesh_nodal_dump_vtk(pmn[i_domain],
                                   PDM_GEOMETRY_KIND_SURFACIC,
                                   "part_mesh_surfacic");
    }
  }

  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part      = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_cell        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_face        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pcell_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int         ***pcell_face     = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pcell_face_idx = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pface_vtx_idx  = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pface_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  double      ***pvtx_coord     = (double      ***) malloc( n_domain * sizeof(double      **));
  double      ***cell_center    = (double      ***) malloc( n_domain * sizeof(double      **));

  int* n_group_by_domain = (int *) malloc(n_domain * sizeof(int));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part     [i_dom] = n_part;
    pn_cell       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_face       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_edge       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_vtx        [i_dom] = (int          *) malloc( n_part * sizeof(int          ));

    pcell_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pface_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pvtx_ln_to_gn [i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));

    pcell_face_idx[i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pcell_face    [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pface_vtx_idx [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pface_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pvtx_coord    [i_dom] = (double      **) malloc( n_part * sizeof(double      *));
    cell_center   [i_dom] = (double      **) malloc( n_part * sizeof(double      *));
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {

      pn_cell[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_CELL,
                                                               &pcell_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);

      pn_face[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_FACE,
                                                               &pface_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_edge[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_EDGE,
                                                               &pedge_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VTX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &pcell_face_idx[i_dom][i_part],
                                          &pcell_face[i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);
      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                       &face_edge_idx,
                                                       &face_edge,
                                                       PDM_OWNERSHIP_KEEP);

      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx_idx,
                                          &edge_vtx,
                                          PDM_OWNERSHIP_KEEP);
      pface_vtx    [i_dom][i_part] = NULL;
      pface_vtx_idx[i_dom][i_part] = face_edge_idx;
      _compute_face_vtx(n_face,
                        face_edge_idx,
                        face_edge,
                        edge_vtx,
                        &pface_vtx[i_dom][i_part]);


      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pvtx_coord[i_dom][i_part],
                                       PDM_OWNERSHIP_KEEP);
      int  n_bound = 0;
      int* group_face_idx      = NULL;
      int* group_face          = NULL;
      PDM_g_num_t* face_group_ln_to_gn = NULL;

      PDM_multipart_group_get(mpart,
                              i_dom,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &n_bound,
                              &group_face_idx,
                              &group_face,
                              &face_group_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);
      n_group_by_domain[i_dom] = 1;
      // printf("n_bound = %i \n", n_bound);

      _cell_center_3d(pn_cell[i_dom][i_part],
                      pcell_face_idx[i_dom][i_part],
                      pcell_face[i_dom][i_part],
                      NULL,
                      NULL,
                      pface_vtx_idx[i_dom][i_part],
                      pface_vtx[i_dom][i_part],
                      NULL,
                      pvtx_coord[i_dom][i_part],
                      &cell_center[i_dom][i_part]);

    }
  }

  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   pn_face,
                                                                                   pn_edge,
                                                                                   pn_vtx,
                                                                                   pface_ln_to_gn,
                                                                                   pedge_ln_to_gn,
                                                                                   pvtx_ln_to_gn);


  /*
   * Mesh interpolation
   */
  PDM_field_cell_to_vtx_t* mi = PDM_field_cell_to_vtx_create(n_domain,
                                                             n_part_by_domain,
                                                             n_group_by_domain,
                                                             PDM_CELL_TO_VTX_INTERP_KIND_IDW,
                                                             1,
                                                             comm);
  PDM_field_cell_to_vtx_part_domain_interface_shared_set(mi, pdi);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part_by_domain[i_domain]; ++i_part) {

      PDM_field_cell_to_vtx_part_set(mi,
                                     i_domain,
                                     i_part,
                                     pn_cell       [i_domain][i_part],
                                     0,
                                     0,
                                     pn_vtx        [i_domain][i_part],
                                     NULL,
                                     NULL,
                                     NULL, // pface_edge_idx[i_domain][i_part],
                                     NULL, // pface_edge    [i_domain][i_part],
                                     NULL, // pedge_vtx     [i_domain][i_part],
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL,
                                     NULL,
                                     pvtx_ln_to_gn [i_domain][i_part],
                                     pvtx_coord    [i_domain][i_part]);

      PDM_field_cell_to_vtx_part_mesh_nodal_set(mi, i_domain, pmn[i_domain]);
      // int  n_bound = 0;
      // int* group_face_idx      = 0;
      // int* group_face          = 0;
      // int* face_group_ln_to_gn = 0;

      // PDM_multipart_group_get(mpart,
      //                         i_domain,
      //                         i_part,
      //                         PDM_MESH_ENTITY_FACE,
      //                         &n_bound,
      //                         &group_face_idx,
      //                         &group_face,
      //                         &face_group_ln_to_gn);

      // PDM_field_cell_to_vtx_part_group_set(mi,
      //                                     i_domain,
      //                                     i_part,
      //                                     0,
      //                                     PDM_BOUND_TYPE_FACE,
      //                                     group_face_idx[1]-group_face_idx[0],
      //                                     &group_face[group_face_idx[0]]);


    }
  }
  PDM_field_cell_to_vtx_compute(mi);

  /*
   * Exchange
   */
  // int n_face_group_field = 1;
  // double  ***pfield       = malloc(n_domain * sizeof(double  **));
  // double ****pfield_bound = malloc(n_domain * sizeof(double ***));
  // for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
  //   pfield      [i_domain] = malloc(n_part_by_domain[i_domain] * sizeof(double  *));
  //   pfield_bound[i_domain] = malloc(n_part_by_domain[i_domain] * sizeof(double **));

  //   for(int i_part = 0; i_part < n_part_by_domain[i_domain]; ++i_part) {
  //     pfield      [i_domain][i_part] = malloc(3 * pn_cell[i_domain][i_part] * sizeof(double  ));

  //     for(int i_cell = 0; i_cell < pn_cell[i_domain][i_part]; ++i_cell) {
  //       // pfield      [i_domain][i_part][i_cell] = 1.;
  //       pfield      [i_domain][i_part][3*i_cell  ] = cell_center[i_domain][i_part][3*i_cell  ];
  //       pfield      [i_domain][i_part][3*i_cell+1] = cell_center[i_domain][i_part][3*i_cell+1];
  //       pfield      [i_domain][i_part][3*i_cell+2] = cell_center[i_domain][i_part][3*i_cell+2];
  //     }

  //     int  n_bound = 0;
  //     int* group_face_idx      = 0;
  //     int* group_face          = 0;
  //     int* face_group_ln_to_gn = 0;

  //     PDM_multipart_group_get(mpart,
  //                             i_domain,
  //                             i_part,
  //                             PDM_MESH_ENTITY_FACE,
  //                             &n_bound,
  //                             &group_face_idx,
  //                             &group_face,
  //                             &face_group_ln_to_gn);

  //     pfield_bound[i_domain][i_part] = malloc(n_face_group_field  * sizeof(double *));
  //     for(int i_group = 0; i_group < n_face_group_field; ++i_group) {
  //       int n_face_in_group = group_face_idx[1]-group_face_idx[0];
  //       pfield_bound[i_domain][i_part][i_group] = malloc(3 * n_face_in_group * sizeof(double));

  //       for(int idx_face = 0; idx_face < n_face_in_group; ++idx_face) {

  //         int i_face = group_face[group_face_idx[0]+idx_face]-1;

  //         double center_face[3] = {0., 0., 0.};
  //         for(int idx_vtx = pface_vtx_idx[i_domain][i_part][i_face]; idx_vtx < pface_vtx_idx[i_domain][i_part][i_face+1]; ++idx_vtx){
  //           int i_vtx = pface_vtx[i_domain][i_part][idx_vtx]-1;
  //           for(int k = 0; k < 3; ++k) {
  //             center_face[k] += pvtx_coord[i_domain][i_part][3*i_vtx+k];
  //           }
  //         }
  //         double ipond = 1./((double)pface_vtx_idx[i_domain][i_part][i_face+1] - pface_vtx_idx[i_domain][i_part][i_face] );
  //         for(int k = 0; k < 3; ++k) {
  //           pfield_bound[i_domain][i_part][i_group][3*idx_face+k] = center_face[k] * ipond;
  //         }

  //       }
  //     }
  //   }
  // }

  // // int stride = 1;
  // double ***result_field = NULL;
  // PDM_field_cell_to_vtx_exch(mi,
  //                           PDM_FIELD_KIND_COORDS,
  //                           pfield,
  //                           pfield_bound,
  //                           &result_field);

  /*
   * Visu
   */
  // if(0 == 1) {

  //   for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
  //     for(int i_part = 0; i_part < n_part_by_domain[i_domain]; ++i_part) {

  //       char filename[999];
  //       sprintf(filename, "out_interp_%i_%i.vtk", i_part, i_rank);

  //       int n_vtx = pn_vtx[i_domain][i_part];

  //       int* elmt_vtx = malloc(pn_vtx[i_domain][i_part] * sizeof(int));
  //       double *field_transpose = malloc(3* pn_vtx[i_domain][i_part] * sizeof(double));
  //       for(int i_vtx = 0; i_vtx < pn_vtx[i_domain][i_part]; ++i_vtx) {
  //         elmt_vtx[i_vtx] = i_vtx+1;
  //         field_transpose[        i_vtx] = result_field[i_domain][i_part][3*i_vtx];
  //         field_transpose[  n_vtx+i_vtx] = result_field[i_domain][i_part][3*i_vtx+1];
  //         field_transpose[2*n_vtx+i_vtx] = result_field[i_domain][i_part][3*i_vtx+2];
  //       }

  //       const char* field_name[] = {"result_field_x", "result_field_y", "result_field_z", 0 };
  //       const double * field  [] = {field_transpose, &field_transpose[n_vtx], &field_transpose[2*n_vtx]};
  //       PDM_vtk_write_std_elements_double(filename,
  //                                         n_vtx,
  //                                         pvtx_coord    [i_domain][i_part],
  //                                         pvtx_ln_to_gn [i_domain][i_part],
  //                                         PDM_MESH_NODAL_POINT,
  //                                         pn_vtx        [i_domain][i_part],
  //                                         elmt_vtx,
  //                                         pvtx_ln_to_gn [i_domain][i_part],
  //                                         3,
  //                                         field_name,
  //                      (const double **)  field);

  //       sprintf(filename, "idw_coords_%i_%i.vtk", i_rank, i_part);
  //       PDM_vtk_write_point_cloud(filename,
  //                                 n_vtx,
  //                                 result_field[i_domain][i_part],
  //                                 NULL,
  //                                 NULL);
  //       free(elmt_vtx);
  //       free(field_transpose);

  //     }
  //   }
  // }

  /*
   * Free
   */
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part_by_domain[i_domain]; ++i_part) {
    //   for(int i_group = 0; i_group < n_face_group_field; ++i_group) {
    //     free(pfield_bound[i_domain][i_part][i_group]);
    //   }
    //   free(pfield      [i_domain][i_part]);
    //   free(pfield_bound[i_domain][i_part]);
    //   free(result_field[i_domain][i_part]);
      free(cell_center [i_domain][i_part]);
    }
    // free(pfield      [i_domain]);
    // free(pfield_bound[i_domain]);
    // free(result_field[i_domain]);
    free(cell_center [i_domain]);
  }
  // free(pfield      );
  // free(pfield_bound);
  // free(result_field);
  free(cell_center);


  PDM_field_cell_to_vtx_free(mi);

  free(n_group_by_domain);
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    PDM_part_mesh_nodal_free(pmn[i_domain]);
  }

  free(pmn);

  /*
   *  Free memory
   */
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for(int i_part = 0; i_part < n_part_by_domain[i_dom]; ++i_part) {
      free(pface_vtx    [i_dom][i_part]);
    }
    free(pn_cell       [i_dom]);
    free(pn_face       [i_dom]);
    free(pn_edge       [i_dom]);
    free(pn_vtx        [i_dom]);
    free(pcell_ln_to_gn[i_dom]);
    free(pface_ln_to_gn[i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pcell_face    [i_dom]);
    free(pcell_face_idx[i_dom]);
    free(pface_vtx     [i_dom]);
    free(pface_vtx_idx [i_dom]);
    free(pvtx_coord    [i_dom]);
  }
  free(pn_cell       );
  free(pn_face       );
  free(pn_edge       );
  free(pn_vtx        );
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pedge_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pcell_face);
  free(pcell_face_idx);
  free(pface_vtx);
  free(pface_vtx_idx);
  free(pn_n_part);
  free(pvtx_coord);

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }

  PDM_multipart_free(mpart);

  if(pdi != NULL) {
    PDM_part_domain_interface_free(pdi);
  }
  PDM_domain_interface_free(dom_intrf);
  free(dcube);
  free(dmn);
  free(n_part_by_domain);

  PDM_MPI_Finalize();

  return 0;
}
