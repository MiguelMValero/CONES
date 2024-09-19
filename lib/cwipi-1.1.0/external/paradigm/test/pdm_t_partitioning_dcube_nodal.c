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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
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
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
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

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;

  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  if (post) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  /*
   *  Warm up dmesh from dmesh_nodal
   */
  // PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
  // PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmn);

  // PDM_dmesh_nodal_generate_distribution(dmn);
  // PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  // PDM_dmesh_t* dm = NULL;
  // PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &dm);

  // /*
  //  * Test tranpose
  //  */
  // PDM_g_num_t *dedge_distrib;
  // PDM_dmesh_distrib_get(dm, PDM_MESH_ENTITY_EDGE, &dedge_distrib);
  // int         *dedge_group_idx;
  // PDM_g_num_t *dedge_group;
  // int n_group = PDM_dmesh_bound_get(dm,
  //                                   PDM_BOUND_TYPE_EDGE,
  //                                   &dedge_group,
  //                                   &dedge_group_idx,
  //                                   PDM_OWNERSHIP_KEEP);
  // int *dedge_to_group;
  // int *dedge_to_group_idx;
  // PDM_dgroup_entity_transpose(n_group,
  //                             dedge_group_idx,
  //                             dedge_group,
  //                             dedge_distrib,
  //                             &dedge_to_group_idx,
  //                             &dedge_to_group,
  //                             comm);

  // int dn_edge = dedge_distrib[i_rank+1] - dedge_distrib[i_rank];
  // if(0 == 1) {
  //   PDM_log_trace_array_int (dedge_group_idx, n_group+1, "dedge_group_idx ::");
  //   PDM_log_trace_array_long(dedge_group, dedge_group_idx[n_group], "dedge_group ::");
  //   PDM_log_trace_connectivity_int(dedge_to_group_idx, dedge_to_group, dn_edge, "dedge_to_group ::");
  // }

  // /* Reverse */

  // int         *dgroup_edge_check_idx;
  // PDM_g_num_t *dgroup_edge_check;
  // PDM_dentity_group_transpose(n_group,
  //                             dedge_to_group_idx,
  //                             dedge_to_group,
  //                             dedge_distrib,
  //                             &dgroup_edge_check_idx,
  //                             &dgroup_edge_check,
  //                             comm);
  // if(1 == 1) {
  //   PDM_log_trace_array_int (dedge_group_idx, n_group+1               , "dedge_group_idx ::");
  //   PDM_log_trace_array_long(dedge_group    , dedge_group_idx[n_group], "dedge_group     ::");
  //   PDM_log_trace_array_int (dgroup_edge_check_idx, n_group+1                     , "dgroup_edge_check_idx ::");
  //   PDM_log_trace_array_long(dgroup_edge_check    , dgroup_edge_check_idx[n_group], "dgroup_edge_check     ::");
  // }
  // free(dgroup_edge_check_idx);
  // free(dgroup_edge_check);
  // free(dedge_to_group);
  // free(dedge_to_group_idx);


  /*
   *  Pour le 3D, il faut calculer les edge + les faces (suffit de branché dans dmesh_nodal to dmesh)
   *   --> Appel du partitionnement et rebuild de toutes les connectivités disponibles
   *   --> Transort de tout les elements à la fin
   */


  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_domain,
                                                   &n_part_domains,
                                                   PDM_FALSE,
                                                   part_method,
                                                   PDM_PART_SIZE_HOMOGENEOUS,
                                                   NULL,
                                                   comm,
                                                   PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart_id, 0, dmn);
  // PDM_multipart_dmesh_set(mpart_id, 0, dm);

  PDM_multipart_compute(mpart_id);

  PDM_multipart_free(mpart_id);


  // A faire : Extraction d'un maillage surfacique à partir du volumique (en dmesh et dmesh_nodal )
  //           Et également en part_mesh et part_mesh_nodal


  PDM_dcube_nodal_gen_free(dcube);
  // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);

  PDM_MPI_Finalize();

  return 0;
}
