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
  PDM_Mesh_nodal_elt_t t_elt   = PDM_MESH_NODAL_HEXA8;
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

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dedge_vtx_idx;
  PDM_g_num_t *dedge_vtx;
  int dn_edge = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                           &dedge_vtx,
                                           &dedge_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx ::");
  }

  int         *dbound_edge_idx = NULL;
  PDM_g_num_t *dbound_edge     = NULL;
  int n_group = PDM_dmesh_bound_get(dmesh,
                                    PDM_BOUND_TYPE_EDGE,
                                    &dbound_edge,
                                    &dbound_edge_idx,
                                    PDM_OWNERSHIP_KEEP);

  /* Extract */
  PDM_dmesh_extract_t *dme = PDM_dmesh_extract_create(1, comm);


  PDM_dmesh_extract_dn_entity_set(dme, PDM_MESH_ENTITY_EDGE, dn_edge);
  PDM_dmesh_extract_dn_entity_set(dme, PDM_MESH_ENTITY_VTX,  dn_vtx);
  PDM_dmesh_extract_dconnectivity_set(dme,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      dedge_vtx,
                                      dedge_vtx_idx);
  PDM_dmesh_extract_vtx_coord_set(dme, dvtx_coord);


  PDM_dmesh_extract_dmesh_bound_set(dme,
                                    PDM_BOUND_TYPE_EDGE,
                                    n_group,
                                    dbound_edge,
                                    dbound_edge_idx);

  PDM_dmesh_extract_selected_gnum_set(dme,
                                      PDM_MESH_ENTITY_EDGE,
                                      dbound_edge_idx[n_group],
                                      dbound_edge);

  PDM_dmesh_extract_compute(dme);

  /*
   * Post-traitement
   */
  PDM_dmesh_t* dmesh_extract = NULL;
  PDM_dmesh_extract_dmesh_get(dme,
                              &dmesh_extract,
                              PDM_OWNERSHIP_KEEP);

  int         *dextract_edge_vtx_idx = NULL;
  PDM_g_num_t *dextract_edge_vtx     = NULL;
  int dn_extract_edge = PDM_dmesh_connectivity_get(dmesh_extract, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                   &dextract_edge_vtx,
                                                   &dextract_edge_vtx_idx,
                                                   PDM_OWNERSHIP_KEEP);

  assert(dextract_edge_vtx_idx == NULL);
  dextract_edge_vtx_idx = malloc( (dn_extract_edge+1) * sizeof(int));
  for(int i_edge = 0; i_edge < dn_extract_edge+1; ++i_edge) {
    dextract_edge_vtx_idx[i_edge] = 2*i_edge;
  }


  double *dextract_vtx_coord = NULL;
  PDM_dmesh_vtx_coord_get(dmesh_extract,
                          &dextract_vtx_coord,
                          PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *distrib_edge_extract = PDM_compute_entity_distribution(comm, dn_extract_edge);

  PDM_g_num_t* extract_edge_ln_to_gn = malloc(dn_extract_edge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_extract_edge; ++i) {
    extract_edge_ln_to_gn[i] = distrib_edge_extract[i_rank] + i + 1;
  }

  int dn_extract_vtx  = PDM_dmesh_dn_entity_get(dmesh_extract, PDM_MESH_ENTITY_VTX);
  PDM_g_num_t *extract_vtx_distribution = PDM_compute_entity_distribution(comm, dn_extract_vtx);

  int pn_extract_vtx = -1;
  PDM_g_num_t *pextract_vtx_ln_to_gn = NULL;
  int         *pextract_edge_vtx_idx = NULL;
  int         *pextract_edge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge_extract,
                                                           dextract_edge_vtx_idx,
                                                           dextract_edge_vtx,
                                                           dn_extract_edge,
                                                           extract_edge_ln_to_gn,
                                                           &pn_extract_vtx,
                                                           &pextract_vtx_ln_to_gn,
                                                           &pextract_edge_vtx_idx,
                                                           &pextract_edge_vtx);

  double** tmp_pextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        extract_vtx_distribution,
                                        dextract_vtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pextract_vtx_ln_to_gn,
                                        &tmp_pextract_vtx_coord);
  double* pextract_vtx_coord = tmp_pextract_vtx_coord[0];
  free(tmp_pextract_vtx_coord);
  free(extract_vtx_distribution);

  if (post) {
    char filename[999];
    sprintf(filename, "export_edge_%i.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx,
                           pextract_vtx_coord,
                           pextract_vtx_ln_to_gn,
                           dn_extract_edge,
                           pextract_edge_vtx_idx,
                           pextract_edge_vtx,
                           extract_edge_ln_to_gn,
                           NULL);
  }

  free(pextract_vtx_ln_to_gn);
  free(pextract_edge_vtx_idx);
  free(pextract_edge_vtx    );
  free(pextract_vtx_coord   );
  free(dextract_edge_vtx_idx);

  free(distrib_edge_extract);
  free(extract_edge_ln_to_gn);

  PDM_dmesh_extract_free(dme);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
