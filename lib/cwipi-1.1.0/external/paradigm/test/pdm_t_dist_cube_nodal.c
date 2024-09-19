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
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_gnum.h"
#include "pdm_error.h"
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_gnum_location.h"
#include "pdm_distrib.h"

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
     "  -t               Element type.\n\n"
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
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t          *n_vtx_seg,
           double                *length,
           int                   *n_part,
           int                   *post,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *part_method)
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
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
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_HEXA8;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             &elt_type,
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

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         0.,
                                                         0.,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if(post) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

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

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE",
                                                     NULL,
                                                     "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart_id, 0, dmn);
  PDM_multipart_compute(mpart_id);

  PDM_part_mesh_nodal_t *pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart_id, 0, &pmesh_nodal, PDM_OWNERSHIP_KEEP);

  /*
   * Extract part to select the bnd to compute distance
   */
  // int n_part_out = 1;
  // int dim_extract = PDM_Mesh_nodal_elt_dim_get(elt_type)-1;
  // PDM_extract_part_t* extrp = PDM_extract_part_create(dim_extract,
  //                                                     n_part,
  //                                                     n_part_out,
  //                                                     PDM_EXTRACT_PART_KIND_LOCAL,
  //                                                     PDM_SPLIT_DUAL_WITH_PTSCOTCH,
  //                                                     PDM_TRUE, // compute_child_gnum
  //                                                     PDM_OWNERSHIP_KEEP,
  //                                                     comm);
  // PDM_part_mesh_nodal_elmts_t *pmne = NULL;
  // if(dim_extract == 2) {
  //   pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC);
  // } else {
  //   pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmesh_nodal, PDM_GEOMETRY_KIND_RIDGE);
  // }
  // PDM_extract_part_part_nodal_set(extrp, pmne);


  int  *pn_selected    = malloc(n_part * sizeof(int  ));
  int **selected_l_num = malloc(n_part * sizeof(int *));
  // for(int i_part = 0; i_part < n_part; ++i_part) {

  //   PDM_g_num_t *vtx_ln_to_gn  = NULL;
  //   PDM_g_num_t *cell_ln_to_gn = NULL;
  //   double      *vtx_coord     = NULL;

  //   int n_vtx  = PDM_multipart_part_ln_to_gn_get(mpart_id, 0, i_part, PDM_MESH_ENTITY_VTX,  &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);
  //   int n_cell = PDM_multipart_part_ln_to_gn_get(mpart_id, 0, i_part, PDM_MESH_ENTITY_CELL, &cell_ln_to_gn, PDM_OWNERSHIP_KEEP);
  //   PDM_multipart_part_vtx_coord_get(mpart_id, 0, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

  //   PDM_extract_part_part_set(extrp,
  //                             i_part,
  //                             n_cell,
  //                             0,
  //                             0, // pn_edge[i_part],
  //                             n_vtx,
  //                             NULL, // pcell_face_idx[i_part],
  //                             NULL, // pcell_face[i_part],
  //                             NULL, // pface_edge_idx[i_part],
  //                             NULL, // pface_edge[i_part],
  //                             NULL, // pedge_vtx[i_part],
  //                             NULL, // pface_vtx_idx[i_part],
  //                             NULL, // pface_vtx[i_part],
  //                             cell_ln_to_gn,
  //                             NULL,
  //                             NULL, //pedge_ln_to_gn[i_part],
  //                             vtx_ln_to_gn,
  //                             vtx_coord);

  //   /*
  //    * Extract group
  //    */
  //   int n_group = PDM_part_mesh_nodal_elmts_n_group_get(pmne);
  //   assert(n_group > 0);
  //   int i_group_extract = 5;

  //   int  n_group_elmt = 0;
  //   int         *group_elmt     = 0;
  //   PDM_g_num_t *group_ln_to_gn = 0;
  //   PDM_part_mesh_nodal_elmts_group_get(pmne,
  //                                       i_part,
  //                                       i_group_extract,
  //                                       &n_group_elmt,
  //                                       &group_elmt,
  //                                       &group_ln_to_gn);
  //   pn_selected   [i_part] = n_group_elmt;
  //   selected_l_num[i_part] = group_elmt;

  //   PDM_log_trace_array_int(group_elmt, n_group_elmt, "group_elmt :");

  //   PDM_extract_part_selected_lnum_set(extrp,
  //                                      i_part,
  //                                      pn_selected[i_part],
  //                                      selected_l_num[i_part]);


  // }


  // PDM_extract_part_compute(extrp);


  // // PDM_part_mesh_nodal_elmts_t* extract_pmne = NULL;
  // // PDM_extract_part_part_mesh_nodal_get(extrp, &extract_pmne, PDM_OWNERSHIP_KEEP);

  // PDM_extract_part_free(extrp);

  PDM_multipart_free(mpart_id);

  PDM_dcube_nodal_gen_free(dcube);
  PDM_part_mesh_nodal_free(pmesh_nodal);

  free(pn_selected);
  free(selected_l_num);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}


