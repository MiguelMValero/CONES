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
#include "pdm_dmesh_nodal_to_dmesh.h"

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



static void
_compute_face_vtx2
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx_out
 )
{

  *face_vtx_out = malloc (sizeof(int) * face_edge_idx[n_face]);
  int *face_vtx = *face_vtx_out;

  int *face_vtx_idx = malloc (sizeof(int) * (n_face+1));

  for(int i_face = 0; i_face <  face_edge_idx[n_face]; ++i_face) {
    face_vtx[i_face] = 10000;
  }

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *is_treated = PDM_array_const_int(n_face, -1);
  face_vtx_idx[0] = 0;
  for (int i_face = 0; i_face < n_face; i_face++) {

    int n_vtx_on_face = (face_edge_idx[i_face+1] - face_edge_idx[i_face]);
    face_vtx_idx[i_face+1] = face_vtx_idx[i_face] + n_vtx_on_face;

    int idx_write = face_vtx_idx[i_face];

    // printf(" --------------------------------------- i_face = %i \n", i_face);
    // for(int idx_edge = face_edge_idx[i_face]; idx_edge < face_edge_idx[i_face+1]; ++idx_edge) {
    //   int i_edge = PDM_ABS(face_edge[idx_edge]) - 1;
    //   int i_vtx1 = edge_vtx[2*i_edge  ];
    //   int i_vtx2 = edge_vtx[2*i_edge+1];
    //   printf("i_edge = %i (%i)-> %i /%i (%i/%i) \n", i_edge, face_edge[idx_edge], i_vtx1, i_vtx2, i_vtx1-1, i_vtx2-1);
    // }

    // Find first vtx of face
    int idx_first_edge = face_edge_idx[i_face];
    int i_edge_first = PDM_ABS (face_edge[idx_first_edge]) - 1;
    int sgn          = PDM_SIGN(face_edge[idx_first_edge]);
    int next_vtx = -1;
    if(sgn == 1) {
      face_vtx[idx_write++] = edge_vtx[2*i_edge_first  ];
      next_vtx = PDM_ABS(edge_vtx[2*i_edge_first+1]);
    } else {
      face_vtx[idx_write++] = edge_vtx[2*i_edge_first+1];
      next_vtx = PDM_ABS(edge_vtx[2*i_edge_first]);
    }

    is_treated[0] = 1;
    int found = 1;
    int i_step = 0;
    while(found < n_vtx_on_face){
      for(int idx_edge = face_edge_idx[i_face]+1; idx_edge < face_edge_idx[i_face+1]; ++idx_edge) {

        int i_edge = PDM_ABS(face_edge[idx_edge]) - 1;
        // printf("idx_edge :: %i | i_edge ::  %i | treated ::  %i (  %i )\n", idx_edge, i_edge, is_treated[idx_edge-face_edge_idx[i_face]], idx_edge-face_edge_idx[i_face]);

        // Find the next edge that refers to next_vtx
        if(is_treated[idx_edge-face_edge_idx[i_face]] == -1) {
          int sgn2   = PDM_SIGN(face_edge[idx_edge]);
          int i_vtx1 = PDM_ABS (edge_vtx[2*i_edge  ]);
          int i_vtx2 = PDM_ABS (edge_vtx[2*i_edge+1]);

          if(i_vtx1 == next_vtx) {
            if(sgn2 == 1) {
              face_vtx[idx_write++] = i_vtx1;
              next_vtx = i_vtx2;
            } else {
              face_vtx[idx_write++] = i_vtx2;
              next_vtx = i_vtx1;
            }

            is_treated[idx_edge-face_edge_idx[i_face]] = 1;
            found++;
          } else if(i_vtx2 == next_vtx) {

            if(sgn2 == 1) {
              face_vtx[idx_write++] = i_vtx1;
              next_vtx = i_vtx2;
            } else {
              face_vtx[idx_write++] = i_vtx2;
              next_vtx = i_vtx1;
            }

            is_treated[idx_edge-face_edge_idx[i_face]] = 1;
            found++;
          }
        }
      } /* End for */

      if(i_step > 10) {
        printf(" abort for i_face = %i with foud = %i \n", i_face, found);
        abort();
      }
      i_step++;

    } /* End while */

    /*
     * Reset
     */
    for(int i = 0; i < found; ++i) {
      is_treated[i] = -1;
    }
  }
  // PDM_log_trace_connectivity_int(face_vtx_idx, face_vtx, n_face, "face_vtx ::");

  free(face_vtx_idx);
  free(is_treated);

}

/**
 *
 * \brief  Main
 *
 */
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -t 8
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -ni 2 -pi 1
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi -pj
// Repro LS89 : mpirun -np 2 ./test/pdm_t_join_domains -nx 5 -ny 4 -nz 2 -pt-scotch -pj
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -pt-scotch -n_part 2
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
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TETRA4;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_PRISM6;
  // 2 -> tria
  // 3 -> quad
  // 5 -> tetra
  // 6 -> pyramid
  // 7 -> prism
  // 8 -> hexa

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

  // int n_interface =
  // n_dom_j*n_dom_k*(n_dom_i - 1 + periodic_i) +
  // n_dom_k*n_dom_i*(n_dom_j - 1 + periodic_j) +
  // n_dom_i*n_dom_j*(n_dom_k - 1 + periodic_k);

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
  // /* Check */
  // if (n_rank == 1) {

  //   double vect[4][3] = {
  //     {0.,             0.,             0.},
  //     {length*n_dom_i, 0.,             0.},
  //     {0,              length*n_dom_j, 0.},
  //     {0.,             0.,             length*n_dom_k}
  //   };

  //   double **all_coord = (double **) malloc(sizeof(double *) * n_domain);
  //   for (int i = 0; i < n_domain; i++) {
  //     all_coord[i] = PDM_DMesh_nodal_vtx_get(dmn[i]);
  //   }

  //   for (i_interface = 0; i_interface < n_interface; i_interface++) {

  //     int i_domain1 = interface_dom[i_interface][0];// - 1?
  //     int i_domain2 = interface_dom[i_interface][1];// - 1?

  //     for (int i = 0; i < interface_dn[i_interface]; i++) {

  //       PDM_g_num_t ivtx_1 = interface_ids[i_interface][2*i]   - 1;
  //       PDM_g_num_t ivtx_2 = interface_ids[i_interface][2*i+1] - 1;

  //       // if(ivtx_1 < 0 || ivtx_1 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1])) {
  //       //   printf("error : interface %d, domain 1:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain1, i, ivtx_1+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1]));
  //       // }
  //       // if(ivtx_2 < 0 || ivtx_2 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2])) {
  //       //   printf("error : interface %d, domain 2:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain2, i, ivtx_2+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2]));
  //       // }

  //       double *vtx_1 = all_coord[i_domain1] + 3*ivtx_1;
  //       double *vtx_2 = all_coord[i_domain2] + 3*ivtx_2;

  //       double dist2 = 0;
  //       for (int j = 0; j < 3; j++) {
  //         double delta = vtx_1[j] - (vtx_2[j] + vect[i_period[i_interface]][j]);
  //         dist2 += delta * delta;
  //       }

  //       // assert(dist2 < 1e-9);
  //       if (dist2 > 1e-9) {
  //         printf("!! interface %d, domains %d %d, vtx ("PDM_FMT_G_NUM": %f %f %f) ("PDM_FMT_G_NUM": %f %f %f), dist = %f\n",
  //                i_interface, i_domain1, i_domain2,
  //                ivtx_1+1,
  //                vtx_1[0], vtx_1[1], vtx_1[2],
  //                ivtx_2+1,
  //                vtx_2[0], vtx_2[1], vtx_2[2],
  //                sqrt(dist2));
  //       }
  //     }
  //   }
  //   free (all_coord);
  // }

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
    PDM_dmesh_nodal_generate_distribution(dmn[i]);
    PDM_multipart_dmesh_nodal_set(mpart, i, dmn[i]);
  }

  // const int renum_properties_cell[6] = {1024, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {12, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {256, 0, 1, 64, 1, 1};
  // const int renum_properties_cell[6] = {16, 0, 1, 64, 1, 1};
  // PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_HPC",
  //                                           (void * ) renum_properties_cell,
  //                                                    "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_compute(mpart);



  /*
   * Create dmesh_nodal_to_dmesh to setup face and edge
   */
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_domain, comm, PDM_OWNERSHIP_KEEP);
  for (int i = 0; i < n_domain; i++) {
    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, i, dmn[i]);
  }

  PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t **dm = malloc(n_domain * sizeof(PDM_dmesh_t *));

  int          *dn_vtx        = malloc( n_domain * sizeof(int          ));
  int          *dn_face       = malloc( n_domain * sizeof(int          ));
  int          *dn_edge       = malloc( n_domain * sizeof(int          ));
  int         **dface_vtx_idx = malloc( n_domain * sizeof(int         *));
  PDM_g_num_t **dface_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **dedge_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));
  int         **dedge_vtx_idx = malloc( n_domain * sizeof(int *        ));
 
  for (int i = 0; i < n_domain; i++) {
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, i, &dm[i]);

    dn_vtx [i] = PDM_dmesh_dn_entity_get(dm[i], PDM_MESH_ENTITY_VTX);
    dn_face[i] = PDM_dmesh_connectivity_get(dm[i],
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &dface_vtx[i],
                                            &dface_vtx_idx[i],
                                            PDM_OWNERSHIP_KEEP);
    assert(dface_vtx_idx[i] != NULL);
    assert(dface_vtx[i] != NULL);

    dedge_vtx    [i] = NULL;
    dedge_vtx_idx[i] = NULL;
    dn_edge[i] = PDM_dmesh_connectivity_get(dm[i],
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &dedge_vtx[i],
                                            &dedge_vtx_idx[i],
                                            PDM_OWNERSHIP_KEEP);
    assert(dedge_vtx_idx[i] != NULL);
    assert(dedge_vtx[i] != NULL);

    // PDM_log_trace_array_int(dedge_vtx_idx[i], dn_edge[i], "dedge_vtx_idx[i] :: ");
    // PDM_log_trace_array_long(dedge_vtx[i], dedge_vtx_idx[i][dn_edge[i]], "dedge_vtx[i] :: ");
    // PDM_log_trace_connectivity_long(dedge_vtx_idx[i], dedge_vtx[i], dn_edge[i], "dedge_vtx ::");

  }

  /*
   * Transform interface by vtx by interface by faces
   */
  // PDM_domain_interface_translate_vtx2face(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_face,
  //                                         dface_vtx_idx,
  //                                         dface_vtx);

  // PDM_domain_interface_translate_vtx2edge(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_edge,
  //                                         dedge_vtx_idx,
  //                                         dedge_vtx);

  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part      = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_face        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge        = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int         ***pface_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part     [i_dom] = n_part;
    pn_face       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_edge       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pface_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pn_vtx        [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pvtx_ln_to_gn [i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pface_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
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
    }
  }


  // int         *graph_vtx_idx = NULL;
  // PDM_g_num_t *graph_vtx_ids = NULL;
  // int         *graph_vtx_dom = NULL;
  // int graph_vtx_dn = PDM_domain_interface_get_as_graph(dom_intrf, PDM_BOUND_TYPE_VTX,
  //                                                      &graph_vtx_idx, &graph_vtx_ids, &graph_vtx_dom);
  // PDM_log_trace_array_int(graph_vtx_idx, graph_vtx_dn+1, "vtx graph idx");
  // PDM_log_trace_array_long(graph_vtx_ids, graph_vtx_idx[graph_vtx_dn], "vtx graph gnums");
  // PDM_log_trace_array_int(graph_vtx_dom, graph_vtx_idx[graph_vtx_dn], "vtx graph dom");


  /*
   * Deduction en partition du graphe entre domaine
   */
  // PDM_part_domain_interface_t* pdi = NULL;
  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   pn_face,
                                                                                   pn_edge,
                                                                                   pn_vtx,
                                                                                   pface_ln_to_gn,
                                                                                   pedge_ln_to_gn,
                                                                                   pvtx_ln_to_gn);

  /*
   * Extension
   */
  int n_depth = 1;
  PDM_extend_type_t  extend_type = PDM_EXTEND_FROM_VTX;
  // PDM_extend_type_t  extend_type = PDM_EXTEND_FROM_FACE;
  // int n_depth = 2;
  PDM_part_extension_t* part_ext = PDM_part_extension_create(n_domain,
                                                             n_part_by_domain,
                                                             extend_type,
                                                             n_depth,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);
  PDM_part_extension_part_domain_interface_shared_set(part_ext, pdi);
  int shift_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){

      PDM_g_num_t* cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_dom,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      int *cell_face     = NULL;
      int *cell_face_idx = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_KEEP);

      int *face_vtx     = NULL;
      int *face_vtx_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx_idx,
                                          &face_vtx,
                                          PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* face_ln_to_gn = NULL;
      int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   i_dom,
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE,
                                                   &face_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                  i_dom,
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

      double *vtx = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &vtx,
                                       PDM_OWNERSHIP_KEEP);


      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &vtx_part_bound_proc_idx,
                                        &vtx_part_bound_part_idx,
                                        &vtx_part_bound,
                                        PDM_OWNERSHIP_KEEP);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face2 = PDM_multipart_part_connectivity_get(mpart,
                                                        i_dom,
                                                        i_part,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge_idx,
                                                        &face_edge,
                                                        PDM_OWNERSHIP_KEEP);
      // assert(n_face == n_face2);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx_idx,
                                                       &edge_vtx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);

      int *face_cell     = NULL;
      int *face_cell_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                          &face_cell_idx,
                                          &face_cell,
                                          PDM_OWNERSHIP_KEEP);
      assert(face_cell_idx == NULL);

      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    i_dom,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

      PDM_UNUSED(n_face2);
      PDM_UNUSED(n_edge2);

      int          pn_face_group       = 0;
      int*         group_face_idx      = NULL;
      int*         group_face          = NULL;
      PDM_g_num_t* group_face_ln_to_gn = NULL;
      PDM_multipart_group_get(mpart,
                              i_dom,
                              i_part,
                              PDM_MESH_ENTITY_FACE,
                              &pn_face_group,
                              &group_face_idx,
                              &group_face,
                              &group_face_ln_to_gn,
                              PDM_OWNERSHIP_KEEP);



      // Concatenate face_vtx
      pface_vtx[i_dom][i_part] = NULL;
      _compute_face_vtx(n_face2,
                        face_edge_idx,
                        face_edge,
                        edge_vtx,
                        &pface_vtx[i_dom][i_part]);
      // _compute_face_vtx2(n_face2,
      //                   face_edge_idx,
      //                   face_edge,
      //                   edge_vtx,
      //                   &pface_vtx[i_dom][i_part]);

      int *face_part_bound_proc_idx = NULL;
      int *face_part_bound_part_idx = NULL;
      int *face_part_bound          = NULL;
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_part_bound_proc_idx,
                                        &face_part_bound_part_idx,
                                        &face_part_bound,
                                        PDM_OWNERSHIP_KEEP);

      // assert(n_edge2 == n_edge);
      PDM_part_extension_set_part(part_ext, i_dom, i_part,
                                  n_cell,
                                  n_face,
                                  0, // n_part_joins,
                                  pn_face_group,
                                  n_edge,
                                  n_vtx,
                                  cell_face_idx,
                                  cell_face,
                                  face_cell,
                                  face_edge_idx,
                                  face_edge,
                                  face_edge_idx,
                                  pface_vtx[i_dom][i_part],
                                  edge_vtx,
                                  group_face_idx,
                                  group_face,
                                  NULL, // face_join_idx
                                  NULL, // face_join
                                  face_part_bound_proc_idx,
                                  face_part_bound_part_idx,
                                  face_part_bound,
                                  vtx_part_bound_proc_idx,
                                  vtx_part_bound_part_idx,
                                  vtx_part_bound,
                                  cell_ln_to_gn,
                                  face_ln_to_gn,
                                  edge_ln_to_gn, // edge_ln_to_gn
                                  vtx_ln_to_gn,
                                  group_face_ln_to_gn,
                                  vtx);

      // for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
      //   int i_vtx1 = PDM_ABS(edge_vtx[2*i_edge  ])-1;
      //   int i_vtx2 = PDM_ABS(edge_vtx[2*i_edge+1])-1;
      //   log_trace("i_edge(%i) = %i / %i | %i = %i / %i\n", i_edge,
      //             i_vtx1, i_vtx2, edge_ln_to_gn[i_edge],
      //             vtx_ln_to_gn[i_vtx1], vtx_ln_to_gn[i_vtx2]);
      //   // log_trace("i_edge(%i) = %i / %i \n", edge_ln_to_gn[i_edge],
      //   //           vtx_ln_to_gn[i_vtx1], vtx_ln_to_gn[i_vtx2]);
      // }

      /*
       *  Mini-Bricoloage
       */
      if(0 == 1) {
        PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn[i_dom],
                                                                                         PDM_GEOMETRY_KIND_VOLUMIC,
                                                                                         1, // n_part
                                                                                         &n_vtx,
                                                                                         &vtx_ln_to_gn,
                                                                                         &n_cell,
                                                                                         &cell_ln_to_gn,
                                                                                         NULL);
        int id_section = 0;
        PDM_Mesh_nodal_elt_t t_elt_loc = PDM_part_mesh_nodal_elmts_section_type_get(pmne_vol, id_section);
        int         *elmt_vtx                 = NULL;
        int         *parent_num               = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_section_std_get(pmne_vol,
                                                  id_section,
                                                  0,
                                                  &elmt_vtx,
                                                  &numabs,
                                                  &parent_num,
                                                  &parent_entitity_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

        int* cell_num = (int *) malloc(n_cell * sizeof(int));
        for(int i = 0; i < n_cell; ++i) {
          cell_num[i] = i;
        }

        if (post) {
          char filename[999];
          sprintf(filename, "out_volumic_%i_%i_%i.vtk", i_dom, i_part, i_rank);
          const char* field_name[] = {"cell_num", 0 };
          const int * field     [] = {cell_num};
          PDM_vtk_write_std_elements(filename,
                                     n_vtx,
                                     vtx,
                                     vtx_ln_to_gn,
                                     t_elt_loc,
                                     n_cell,
                                     elmt_vtx,
                                     cell_ln_to_gn,
                                     1,
                                     field_name,
                                     (const int **)  field);
        }
        free(cell_num);

        PDM_part_mesh_nodal_elmts_free(pmne_vol);
      }
    }
    shift_part += pn_n_part[i_dom];
  }
  PDM_UNUSED(shift_part);


  double t1 = PDM_MPI_Wtime();
  PDM_part_extension_compute(part_ext);
  double dt = PDM_MPI_Wtime() - t1;
  printf("PDM_part_extension_compute : %12.5e \n", dt);
  t1 = PDM_MPI_Wtime();

  int         *composed_interface_idx   = NULL;
  int         *composed_interface       = NULL;
  PDM_g_num_t *composed_ln_to_gn_sorted = NULL;
  int n_interf_composed = PDM_part_extension_composed_interface_get(part_ext,
                                                                    &composed_interface_idx,
                                                                    &composed_interface,
                                                                    &composed_ln_to_gn_sorted);

  if(0 == 1) {
    PDM_log_trace_connectivity_int(composed_interface_idx, composed_interface, n_interf_composed, "composed_interface ::");
    PDM_log_trace_array_long(composed_ln_to_gn_sorted, n_interf_composed, "composed_ln_to_gn_sorted ::");
  }

  /*
   * Export current domain with elements
   *   - Interface with dmesh_nodal
   *   - Export vtk
   */
  shift_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++){
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {

      PDM_g_num_t* cell_ln_to_gn = NULL;
      int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   i_dom,
                                                   i_part,
                                                   PDM_MESH_ENTITY_CELL,
                                                   &cell_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                   i_dom,
                                                   i_part,
                                                   PDM_MESH_ENTITY_VTX,
                                                   &vtx_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

      double *vtx = NULL;
      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &vtx,
                                       PDM_OWNERSHIP_KEEP);

      double* vtx_coord_extended;
      PDM_g_num_t* border_vtx_ln_to_gn;
      PDM_g_num_t* border_cell_ln_to_gn;
      PDM_g_num_t* border_face_ln_to_gn;
      int n_vtx_extended  = PDM_part_extension_vtx_coord_get(part_ext, i_dom, i_part, &vtx_coord_extended);
      int n_vtx_extended2 = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_VTX,  &border_vtx_ln_to_gn);
      int n_cell_extended = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_CELL, &border_cell_ln_to_gn);
      int n_face_extended = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_FACE, &border_face_ln_to_gn);
      assert(n_vtx_extended == n_vtx_extended2);


      int *border_vtx_interface  = NULL;
      int *border_face_interface = NULL;
      int *border_cell_interface = NULL;
      PDM_part_extension_interface_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_VTX,  &border_vtx_interface);
      PDM_part_extension_interface_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_FACE, &border_face_interface);
      PDM_part_extension_interface_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_CELL, &border_cell_interface);

      if(0 == 1) {
        PDM_log_trace_array_long(border_vtx_ln_to_gn , n_vtx_extended , "border_vtx_ln_to_gn :: ");
        PDM_log_trace_array_long(border_face_ln_to_gn, n_face_extended, "border_face_ln_to_gn :: ");
        PDM_log_trace_array_long(border_cell_ln_to_gn, n_cell_extended, "border_cell_ln_to_gn :: ");

        PDM_log_trace_array_int (border_vtx_interface , n_vtx_extended , "border_vtx_interface :: ");
        PDM_log_trace_array_int (border_face_interface, n_face_extended, "border_face_interface :: ");
        PDM_log_trace_array_int (border_cell_interface, n_cell_extended, "border_cell_interface :: ");
        PDM_log_trace_array_long(cell_ln_to_gn, n_cell, "cell_ln_to_gn :: ");
      }

      /* Apply periodicity */
      int n_interface = PDM_part_domain_interface_n_interface_get(pdi);
      double **translation_vector = malloc(n_interface * sizeof(double * ));
      for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
        translation_vector[i_interf] = NULL;
        PDM_part_domain_interface_translation_get(pdi, i_interf, &translation_vector[i_interf]);
      }

      char filename_pts[999];
      if(0 == 1) {
        sprintf(filename_pts, "out_vtx_extended_raw_%i_%i_%i.vtk", i_dom, i_part, i_rank);
        PDM_vtk_write_point_cloud(filename_pts,
                                  n_vtx_extended,
                                  vtx_coord_extended,
                                  border_vtx_ln_to_gn,
                                  border_vtx_interface);
      }

      int *i_num_vtx = malloc(n_vtx_extended * sizeof(int));
      for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
        i_num_vtx[i_vtx] = i_vtx;
        // i_num_vtx[i_vtx] = border_vtx_ln_to_gn[i_vtx];
        if(border_vtx_interface[i_vtx] == -40000) {
          continue;
        }
        // int i_interf = PDM_ABS(border_vtx_interface[i_vtx])-1;
        // if(i_interf < n_interface){
        //   for(int k = 0; k < 3; ++k) {
        //     vtx_coord_extended[3*i_vtx+k] += PDM_SIGN(border_vtx_interface[i_vtx]) * translation_vector[i_interf][k];
        //   }
        // } else {
        //   // int l_interf = PDM_ABS(border_vtx_interface[i_vtx]) - n_interface - 1;
        //   int l_interf = PDM_binary_search_long(border_vtx_interface[i_vtx], composed_ln_to_gn_sorted, n_interf_composed);
        //   for(int idx_comp = composed_interface_idx[l_interf]; idx_comp < composed_interface_idx[l_interf+1]; ++idx_comp) {
        //     int i_tr = composed_interface[idx_comp];
        //     for(int k = 0; k < 3; ++k) {
        //       vtx_coord_extended[3*i_vtx+k] += PDM_SIGN(i_tr) * translation_vector[PDM_ABS(i_tr)-1][k];
        //     }
        //   }
        // }
      }

      if(0 == 1) {
        sprintf(filename_pts, "out_vtx_extended_%i_%i_%i.vtk", i_dom, i_part, i_rank);
        // PDM_vtk_write_point_cloud(filename_pts,
        //                           n_vtx_extended,
        //                           vtx_coord_extended,
        //                           NULL,
        //                           border_vtx_interface);
        PDM_vtk_write_point_cloud(filename_pts,
                                  n_vtx_extended,
                                  vtx_coord_extended,
                                  NULL,
                                  i_num_vtx);
      }
      free(i_num_vtx);

      // PDM_g_num_t* border_face_ln_to_gn;
      // PDM_g_num_t* border_edge_ln_to_gn;
      // int n_face_extended  = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_FACE, &border_face_ln_to_gn);
      // int n_edge_extended  = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_EDGE, &border_edge_ln_to_gn);

      int *extend_face_edge_idx = NULL;
      int *extend_face_edge     = NULL;
      int n_face_extended2      = 0;
      int n_edge_extended2      = 0;
      int *extend_edge_vtx_idx = NULL;
      int *extend_edge_vtx     = NULL;
      int *extend_face_vtx_idx = NULL;
      int *extend_face_vtx     = NULL;
      if(extend_type == PDM_EXTEND_FROM_VTX) {
        n_face_extended2 = PDM_part_extension_connectivity_get(part_ext, i_dom, i_part,
                                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                               &extend_face_edge_idx,
                                                               &extend_face_edge);

        n_edge_extended2 = PDM_part_extension_connectivity_get(part_ext, i_dom, i_part,
                                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                               &extend_edge_vtx_idx,
                                                               &extend_edge_vtx);
      } else {
        n_face_extended2 = PDM_part_extension_connectivity_get(part_ext, i_dom, i_part,
                                                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                               &extend_face_vtx_idx,
                                                               &extend_face_vtx);
      }
      PDM_UNUSED(n_face_extended2);

      // PDM_g_num_t* edge_ln_to_gn = NULL;
      // int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
      //                                               i_dom,
      //                                               i_part,
      //                                               PDM_MESH_ENTITY_EDGE,
      //                                               &edge_ln_to_gn,
      //                                               PDM_OWNERSHIP_KEEP);

      // int n_edge_tot  = n_edge + n_edge_extended;
      // PDM_g_num_t *concat_edge_ln_to_gn  = malloc(    n_edge_tot  * sizeof(PDM_g_num_t));
      // for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
      //   concat_edge_ln_to_gn[i_edge] = edge_ln_to_gn[i_edge];
      // }
      // for(int i_edge = 0; i_edge < n_edge_extended; ++i_edge) {
      //   concat_edge_ln_to_gn[n_edge+i_edge] = border_edge_ln_to_gn[i_edge];
      // }

      // PDM_log_trace_part_connectivity_gnum(pface_edge_idx,
      //                                      pface_edge,
      //                                      border_face_ln_to_gn,
      //                                      concat_edge_ln_to_gn,
      //                                      n_face_extended,
      //                                      "face_edge_extented");
      // free(concat_edge_ln_to_gn);

      int n_vtx_tot  = n_vtx + n_vtx_extended;
      int n_cell_tot = n_cell + n_cell_extended;

      PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    n_vtx_tot  * sizeof(PDM_g_num_t));
      PDM_g_num_t *concat_cell_ln_to_gn = malloc(    n_cell_tot * sizeof(PDM_g_num_t));
      int         *is_extend            = malloc(    n_cell_tot * sizeof(int        ));
      double      *concat_vtx_coord     = malloc(3 * n_vtx_tot  * sizeof(double     ));
      for(int i_vtx = 0; i_vtx < 3 * n_vtx_tot; ++i_vtx) {
        concat_vtx_coord[i_vtx] = -10000.;
      }

      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
        concat_vtx_ln_to_gn[i_vtx] = vtx_ln_to_gn[i_vtx];
        concat_vtx_coord[3*i_vtx  ] = vtx[3*i_vtx  ];
        concat_vtx_coord[3*i_vtx+1] = vtx[3*i_vtx+1];
        concat_vtx_coord[3*i_vtx+2] = vtx[3*i_vtx+2];
      }

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        concat_cell_ln_to_gn[i_cell] = cell_ln_to_gn[i_cell];
        is_extend           [i_cell] = 0;
      }

      for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
        concat_vtx_ln_to_gn[n_vtx+i_vtx] = border_vtx_ln_to_gn[i_vtx];
        concat_vtx_coord[3*(n_vtx+i_vtx)  ] = vtx_coord_extended[3*i_vtx  ];
        concat_vtx_coord[3*(n_vtx+i_vtx)+1] = vtx_coord_extended[3*i_vtx+1];
        concat_vtx_coord[3*(n_vtx+i_vtx)+2] = vtx_coord_extended[3*i_vtx+2];
      }
      for(int i_cell = 0; i_cell < n_cell_extended; ++i_cell) {
        concat_cell_ln_to_gn[n_cell+i_cell] = border_cell_ln_to_gn[i_cell];
        is_extend           [n_cell+i_cell] = 1;
      }

      if(0 == 1) {
        // PDM_log_trace_part_connectivity_gnum(pedge_vtx_idx,
        //                                      pedge_vtx,
        //                                      border_edge_ln_to_gn,
        //                                      concat_vtx_ln_to_gn,
        //                                      n_edge_extended,
        //                                      "edge_vtx_extented");

        PDM_log_trace_array_long(concat_vtx_ln_to_gn , n_vtx_tot , "concat_vtx_ln_to_gn :: ");
        PDM_log_trace_array_long(concat_cell_ln_to_gn, n_cell_tot, "concat_cell_ln_to_gn :: ");
      }

      PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn[i_dom],
                                                                                       PDM_GEOMETRY_KIND_VOLUMIC,
                                                                                       1, // n_part
                                                                                       &n_vtx_tot,
                                                                                       &concat_vtx_ln_to_gn,
                                                                                       &n_cell_tot,
                                                                                       &concat_cell_ln_to_gn,
                                                                                       NULL);

      int id_section = 0;
      PDM_Mesh_nodal_elt_t t_elt_loc = PDM_part_mesh_nodal_elmts_section_type_get(pmne_vol, id_section);
      int         *elmt_vtx                 = NULL;
      int         *parent_num               = NULL;
      PDM_g_num_t *numabs                   = NULL;
      PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_section_std_get(pmne_vol,
                                                id_section,
                                                0,
                                                &elmt_vtx,
                                                &numabs,
                                                &parent_num,
                                                &parent_entitity_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);

      if (post) {
        char filename[999];
        sprintf(filename, "out_volumic_extended_%i_%i_%i.vtk", i_dom, i_part, i_rank);

        const char* field_name[] = {"is_extend", 0 };
        const int * field     [] = {is_extend};
        PDM_vtk_write_std_elements(filename,
                                   n_vtx_tot,
                                   concat_vtx_coord,
                                   concat_vtx_ln_to_gn,
                                   t_elt_loc,
                                   n_cell_tot,
                                   elmt_vtx,
                                   concat_cell_ln_to_gn,
                                   1,
                                   field_name,
                                   (const int **)   field);
      }
      PDM_part_mesh_nodal_elmts_free(pmne_vol);

      // Rebuild poly connectivity to check mesh
      if(extend_type == PDM_EXTEND_FROM_VTX && post) {
        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        int n_face2 = PDM_multipart_part_connectivity_get(mpart,
                                                          i_dom,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                          &face_edge_idx,
                                                          &face_edge,
                                                          PDM_OWNERSHIP_KEEP);
        // assert(n_face == n_face2);
        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                         i_dom,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                         &edge_vtx_idx,
                                                         &edge_vtx,
                                                         PDM_OWNERSHIP_KEEP);


        int n_face_tot = n_face2+n_face_extended;

        int *concat_face_edge_idx = malloc( (n_face_tot+1) * sizeof(int));
        int *concat_face_edge     = malloc( (face_edge_idx[n_face2] + extend_face_edge_idx[n_face_extended])* sizeof(int));
        int *concat_edge_vtx      = malloc( 2 * (n_edge + n_edge_extended2 ) * sizeof(int));


        // for(int i = 0; )

        PDM_g_num_t *color_face = (PDM_g_num_t *) malloc(n_face_tot * sizeof(PDM_g_num_t));
        concat_face_edge_idx[0] = 0;
        for(int i_face = 0; i_face < n_face2; ++i_face) {
          concat_face_edge_idx[i_face+1] = concat_face_edge_idx[i_face];
          for(int idx_edge = face_edge_idx[i_face]; idx_edge < face_edge_idx[i_face+1]; ++idx_edge) {
            concat_face_edge[concat_face_edge_idx[i_face+1]++] = face_edge[idx_edge];
          }
          color_face[i_face] = 0;
        }

        for(int i_face = 0; i_face < n_face_extended; ++i_face) {
          int l_face = i_face + n_face2;
          concat_face_edge_idx[l_face+1] = concat_face_edge_idx[l_face];
          for(int idx_edge = extend_face_edge_idx[i_face]; idx_edge < extend_face_edge_idx[i_face+1]; ++idx_edge) {
            concat_face_edge[concat_face_edge_idx[l_face+1]++] = extend_face_edge[idx_edge];
          }
          color_face[l_face] = 1;
        }


        for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
          concat_edge_vtx[2*i_edge  ] = edge_vtx[2*i_edge  ];
          concat_edge_vtx[2*i_edge+1] = edge_vtx[2*i_edge+1];
        }

        for(int i_edge = 0; i_edge < n_edge_extended2; ++i_edge) {
          int l_edge = i_edge + n_edge;
          concat_edge_vtx[2*l_edge  ] = extend_edge_vtx[2*i_edge  ];
          concat_edge_vtx[2*l_edge+1] = extend_edge_vtx[2*i_edge+1];
        }

        // PDM_log_trace_connectivity_int(concat_face_edge_idx,
        //                                concat_face_edge,
        //                                n_face_tot,
        //                                "concat_face_edge ::");

        // Concatenate face_vtx
        int *concat_face_vtx = NULL;
        // _compute_face_vtx(n_face_tot,
        //                   concat_face_edge_idx,
        //                   concat_face_edge,
        //                   concat_edge_vtx,
        //                   &concat_face_vtx);
        _compute_face_vtx2(n_face_tot,
                           concat_face_edge_idx,
                           concat_face_edge,
                           concat_edge_vtx,
                           &concat_face_vtx);
        PDM_log_trace_connectivity_int(concat_face_edge_idx,
                                       concat_face_vtx,
                                       n_face_tot,
                                       "concat_face_vtx ::");

        sprintf(filename_pts, "out_face_vtx_%i_%i_%i.vtk", i_dom, i_part, i_rank);
        PDM_vtk_write_polydata(filename_pts,
                               n_vtx_tot,
                               concat_vtx_coord,
                               concat_vtx_ln_to_gn,
                               n_face_tot,
                               concat_face_edge_idx,
                               concat_face_vtx,
                               color_face,
                               NULL);
        free(color_face);

        free(concat_face_vtx);
        free(concat_face_edge_idx);
        free(concat_face_edge    );
        free(concat_edge_vtx     );
      } else if (post){

        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        int n_face2 = PDM_multipart_part_connectivity_get(mpart,
                                                          i_dom,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                          &face_edge_idx,
                                                          &face_edge,
                                                          PDM_OWNERSHIP_KEEP);

        assert(n_face2 == pn_face[i_dom][i_part]);
        int n_face_tot = pn_face[i_dom][i_part]+n_face_extended;

        int *concat_face_vtx_idx = malloc( (n_face_tot+1) * sizeof(int));
        int *concat_face_vtx     = malloc( (face_edge_idx[n_face2] + extend_face_vtx_idx[n_face_extended])* sizeof(int));


        PDM_g_num_t *color_face = (PDM_g_num_t *) malloc(n_face_tot * sizeof(PDM_g_num_t));
        concat_face_vtx_idx[0] = 0;
        for(int i_face = 0; i_face < n_face2; ++i_face) {
          concat_face_vtx_idx[i_face+1] = concat_face_vtx_idx[i_face];
          for(int idx_vtx = face_edge_idx[i_face]; idx_vtx < face_edge_idx[i_face+1]; ++idx_vtx) {
            concat_face_vtx[concat_face_vtx_idx[i_face+1]++] = pface_vtx[i_dom][i_part][idx_vtx];
          }
          color_face[i_face] = 0;
        }

        for(int i_face = 0; i_face < n_face_extended; ++i_face) {
          int l_face = i_face + n_face2;
          concat_face_vtx_idx[l_face+1] = concat_face_vtx_idx[l_face];
          // printf("i_face = %i | l_face = %i | idx = %i \n", i_face, l_face, concat_face_vtx_idx[l_face+1]);
          for(int idx_vtx = extend_face_vtx_idx[i_face]; idx_vtx < extend_face_vtx_idx[i_face+1]; ++idx_vtx) {
            concat_face_vtx[concat_face_vtx_idx[l_face+1]++] = PDM_ABS(extend_face_vtx[idx_vtx]);
          }
          color_face[l_face] = 1;
        }

        // PDM_log_trace_connectivity_int(concat_face_vtx_idx,
        //                                concat_face_vtx,
        //                                n_face_tot,
        //                                "concat_face_vtx ::");

        sprintf(filename_pts, "out_face_vtx_%i_%i_%i.vtk", i_dom, i_part, i_rank);
        PDM_vtk_write_polydata(filename_pts,
                               n_vtx_tot,
                               concat_vtx_coord,
                               concat_vtx_ln_to_gn,
                               n_face_tot,
                               concat_face_vtx_idx,
                               concat_face_vtx,
                               color_face,
                               NULL);
        free(color_face);
        free(concat_face_vtx_idx);
        free(concat_face_vtx);


      }


      // PDM_part_mesh_nodal_elmts_free(pmne_vol);
      free(concat_vtx_ln_to_gn );
      free(concat_cell_ln_to_gn);
      free(concat_vtx_coord    );


      // pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn[i_dom],
      //                                                     PDM_GEOMETRY_KIND_VOLUMIC,
      //                                                     1, // n_part
      //                                                     &n_vtx_extended,
      //                                                     &border_vtx_ln_to_gn,
      //                                                     &n_cell_extended,
      //                                                     &border_cell_ln_to_gn,
      //                                                     NULL);

      // id_section = 0;
      // t_elt_loc = PDM_part_mesh_nodal_elmts_section_type_get(pmne_vol, id_section);
      // PDM_part_mesh_nodal_elmts_section_std_get(pmne_vol, id_section, 0, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

      // sprintf(filename, "out_volumic_only_extended_%i_%i_%i.vtk", i_dom, i_part, i_rank);

      // PDM_vtk_write_std_elements(filename,
      //                            n_vtx_extended,
      //                            vtx_coord_extended,
      //                            border_vtx_ln_to_gn,
      //                            t_elt_loc,
      //                            n_cell_extended,
      //                            elmt_vtx,
      //                            border_cell_ln_to_gn,
      //                            1,
      //                            field_name,
      //           (const int **)   field);

      // PDM_part_mesh_nodal_elmts_free(pmne_vol);


      free(is_extend    );

      for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
        free(translation_vector[i_interf]);
      }
      free(translation_vector);
      free(pface_vtx[i_dom][i_part]);

    }
    shift_part += pn_n_part[i_dom];
  }
  dt = PDM_MPI_Wtime() - t1;
  printf("post-treatment : %12.5e \n", dt);

  PDM_part_extension_free(part_ext);


  /*
   *  Pour le debug : extration des faces avec le extract part + vtk
   */

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    free(pn_face       [i_dom]);
    free(pn_edge       [i_dom]);
    free(pface_ln_to_gn[i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pface_vtx     [i_dom]);
  }
  free(pn_face       );
  free(pn_edge       );
  free(pface_ln_to_gn);
  free(pedge_ln_to_gn);
  free(pn_vtx       );
  free(pvtx_ln_to_gn);
  free(pface_vtx);
  free(pn_n_part);

  /*
   *  Free memory
   */
  // free(i_period);

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }
  PDM_multipart_free(mpart);

  if(pdi != NULL) {
    PDM_part_domain_interface_free(pdi);
  }
  PDM_UNUSED(pdi);

  free(dm);
  free(dn_vtx);
  free(dn_face);
  free(dn_edge);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dedge_vtx_idx);
  free(dedge_vtx);

  PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  PDM_domain_interface_free(dom_intrf);
  free(n_part_by_domain);

  // for (int i = 0; i < n_interface; i++) {
  //   free(interface_ids[i]);
  //   free(interface_dom[i]);
  // }
  // free(interface_dn);
  // free(interface_ids);
  // free(interface_dom);
  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
