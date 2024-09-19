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

/**
 *
 * \brief  Main
 *
 */
// mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -nx 3 -ny 3 -nz 2 -pi -t 8
// mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -nx 3 -ny 3 -nz 2 -t 8 -ni 2 -pi 1
// mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -nx 3 -ny 3 -nz 2 -t 8 -pi
// mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -nx 3 -ny 3 -nz 2 -t 8 -pi -pj
// Debug this one : mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -n 3 -pi -pj -ni 2
// Debug this one : mpirun -np 1 ./test/pdm_t_domain_interface_vtx_translate -nx 4 -ny 3 -nz 2 -post -ni 2
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
  int                  post       = 0;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_QUAD4;
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

  int n_interface =
  n_dom_j*n_dom_k*(n_dom_i - 1 + periodic_i) +
  n_dom_k*n_dom_i*(n_dom_j - 1 + periodic_j) +
  n_dom_i*n_dom_j*(n_dom_k - 1 + periodic_k);

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
   * Create dmesh_nodal_to_dmesh to setup face and edge
   */
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_domain, comm, PDM_OWNERSHIP_KEEP);
  for (int i = 0; i < n_domain; i++) {
    dmn[i] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i]);
    PDM_dmesh_nodal_generate_distribution(dmn[i]);
    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, i, dmn[i]);
  }

  if(dim == 2) {
    PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  } else {
    PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  }

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

    if(dim == 3) {
      dn_face[i] = PDM_dmesh_connectivity_get(dm[i],
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              &dface_vtx[i],
                                              &dface_vtx_idx[i],
                                              PDM_OWNERSHIP_KEEP);
      assert(dface_vtx_idx[i] != NULL);
      assert(dface_vtx[i] != NULL);
    } else {
      dn_face[i] = 0;
    }

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

  }

  /*
   * Transform interface by vtx by interface by faces
   */
  if(dim == 3) {
    // Il faut rajouter un is_implicit_order
    printf("dn_vtx %i \n", dn_vtx[0]);
    printf("dn_face %i \n", dn_face[0]);
    double t1 = PDM_MPI_Wtime();
    // TODO : Extract dface_vtx concerns by interface before
    PDM_domain_interface_translate_vtx2face(dom_intrf,
                                            dn_vtx,
                                            dn_face,
                                            dface_vtx_idx,
                                            dface_vtx);
    double dt = PDM_MPI_Wtime() - t1;
    printf("PDM_domain_interface_translate_vtx2face %12.5e \n", dt);

  }

  PDM_domain_interface_translate_vtx2edge(dom_intrf,
                                          dn_vtx,
                                          dn_edge,
                                          dedge_vtx_idx,
                                          dedge_vtx);

  if(post) {

    if(dim == 3) {
      int          *interface_face_dn  = NULL;
      PDM_g_num_t **interface_face_ids = NULL;
      int         **interface_face_dom = NULL;
      PDM_domain_interface_get(dom_intrf,
                               PDM_BOUND_TYPE_FACE,
                               &interface_face_dn,
                               &interface_face_ids,
                               &interface_face_dom);

      for(int i = 0; i < n_interface; ++i) {
        PDM_log_trace_array_long(interface_face_ids[i], 2 * interface_face_dn[i], "interface_face_ids");
        PDM_log_trace_array_int (interface_face_dom[i], 2 * interface_face_dn[i], "interface_face_dom");
      }
    }

    int          *interface_edge_dn  = NULL;
    PDM_g_num_t **interface_edge_ids = NULL;
    int         **interface_edge_dom = NULL;
    PDM_domain_interface_get(dom_intrf,
                             PDM_BOUND_TYPE_EDGE,
                             &interface_edge_dn,
                             &interface_edge_ids,
                             &interface_edge_dom);

    for(int i = 0; i < n_interface; ++i) {
      PDM_log_trace_array_long(interface_edge_ids[i], 2 * interface_edge_dn[i], "interface_edge_ids");
      PDM_log_trace_array_int (interface_edge_dom[i], 2 * interface_edge_dn[i], "interface_edge_dom");
    }

  }

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
  }

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
  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
