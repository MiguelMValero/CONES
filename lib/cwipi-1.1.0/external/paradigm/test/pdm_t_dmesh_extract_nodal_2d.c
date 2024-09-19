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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_extract.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"

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
_read_args(int                     argc,
           char                  **argv,
           PDM_g_num_t           *n_vtx_seg,
           PDM_Mesh_nodal_elt_t  *t_elt,
           double                *length,
           int                   *post)
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
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
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
  // PDM_Mesh_nodal_elt_t t_elt   = PDM_MESH_NODAL_HEXA8;
  PDM_Mesh_nodal_elt_t t_elt   = PDM_MESH_NODAL_PRISM6;
  int post                     = 0;
  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &t_elt,
             &length,
             &post);

  int order = 1;
  if(PDM_Mesh_nodal_elmt_is_ho(t_elt) == 1) {
    order = 2;
  }

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if(0 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /* Selection on border*/
  int          n_group_face    = 0;
  int         *dgroup_face_idx = NULL;
  PDM_g_num_t *dgroup_face     = NULL;
  PDM_DMesh_nodal_section_group_elmt_get(dmn,
                                         PDM_GEOMETRY_KIND_SURFACIC,
                                         &n_group_face,
                                         &dgroup_face_idx,
                                         &dgroup_face);
  int i_group = 3;
  int dn_extract = dgroup_face_idx[i_group+1] - dgroup_face_idx[i_group];
  PDM_g_num_t *selected = &dgroup_face[dgroup_face_idx[i_group]];

  // PDM_log_trace_array_long(selected, dn_extract, "selected :");

  PDM_dmesh_extract_t *dme = PDM_dmesh_extract_create(2, comm);

  PDM_dmesh_extract_selected_gnum_set(dme,
                                      PDM_MESH_ENTITY_FACE,
                                      dn_extract,
                                      selected);

  PDM_dmesh_extract_dmesh_nodal_set(dme, dmn);
  PDM_dmesh_extract_compute(dme);

  PDM_dmesh_nodal_t* extract_dmn = NULL;
  PDM_dmesh_extract_dmesh_nodal_get(dme, &extract_dmn, PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(extract_dmn, PDM_GEOMETRY_KIND_SURFACIC , "out_extract_surfacic");
  }

  PDM_dmesh_extract_free(dme);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
