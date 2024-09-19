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
#include "pdm_dmesh_to_dmesh_nodal.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_logging.h"

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
           double                *length)
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
  // PDM_Mesh_nodal_elt_t t_elt   = PDM_MESH_NODAL_PRISM6;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &t_elt,
             &length);

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

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
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

  int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);
  PDM_dmesh_nodal_to_dmesh_transform_t       transform_kind       = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE;
  PDM_dmesh_nodal_to_dmesh_translate_group_t transform_group_kind = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE;
  if(dim == 2) {
    transform_kind       = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE;
    transform_group_kind = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE;
  }
  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   transform_kind,
                                   transform_group_kind);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  int dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                           &dface_vtx,
                                           &dface_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  int         *dface_cell_idx = NULL;
  PDM_g_num_t *tmp_dface_cell = NULL;
  dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                       &tmp_dface_cell,
                                       &dface_cell_idx,
                                       PDM_OWNERSHIP_KEEP);
  assert(dface_cell_idx == NULL);

  int         *dcell_face_idx;
  PDM_g_num_t *dcell_face;
  int dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                           &dcell_face,
                                           &dcell_face_idx,
                                           PDM_OWNERSHIP_KEEP);

  int         *dedge_face_idx;
  PDM_g_num_t *dedge_face;
  int dn_edge = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                           &dedge_face,
                                           &dedge_face_idx,
                                           PDM_OWNERSHIP_KEEP);

  int         *dface_edge_idx;
  PDM_g_num_t *dface_edge;
  PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_KEEP);

  int         *dedge_vtx_idx;
  PDM_g_num_t *dedge_vtx;
  int dn_edge2 = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                           &dedge_vtx,
                                           &dedge_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);
  assert(dn_edge == dn_edge2);

  if(0 == 1) {
    PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx ::");
    // PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");
    PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dn_face, "dface_edge ::");
  }

  PDM_UNUSED(dn_face);
  PDM_UNUSED(dn_cell);
  PDM_UNUSED(dn_vtx);
  PDM_UNUSED(dvtx_coord);

  /* Recreate dmesh_nodal */
  PDM_dmesh_to_dmesh_nodal_t *dm_to_dmn = PDM_dmesh_to_dmesh_nodal_create(1, comm);

  PDM_g_num_t *distrib_cell = NULL;
  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  // PDM_g_num_t *distrib_vtx  = NULL;

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_CELL, &distrib_cell);
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &distrib_edge);
  // PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_VTX, &distrib_vtx );

  double *dvtx_coords = PDM_DMesh_nodal_vtx_get(dmn);

  PDM_dmesh_to_dmesh_nodal_distribution_set(dm_to_dmn,
                                            0,
                                            distrib_cell,
                                            distrib_face,
                                            distrib_edge,
                                            vtx_distrib );
  if(0 == 1) {
    PDM_dmesh_to_dmesh_nodal_connectivity_set(dm_to_dmn,
                                              0,
                                              dcell_face_idx,
                                              dcell_face,
                                              dface_edge_idx,
                                              dface_edge,
                                              dedge_vtx,
                                              dface_vtx_idx,
                                              dface_vtx,
                                              dvtx_coords);
  } else {

    if(dim == 2) {
      PDM_dmesh_to_dmesh_nodal_connectivity_set(dm_to_dmn,
                                                0,
                                                dcell_face_idx,
                                                dcell_face,
                                                dface_edge_idx,
                                                dface_edge,
                                                dedge_vtx,
                                                dface_vtx_idx,
                                                dface_vtx,
                                                dvtx_coords);
    } else {
      PDM_dmesh_to_dmesh_nodal_connectivity_set(dm_to_dmn,
                                                0,
                                                dcell_face_idx,
                                                dcell_face,
                                                NULL,
                                                NULL,
                                                NULL,
                                                dface_vtx_idx,
                                                dface_vtx,
                                                dvtx_coords);
    }
  }

  int         *dbound_face_idx = NULL;
  PDM_g_num_t *dbound_face     = NULL;
  int         *dbound_edge_idx = NULL;
  PDM_g_num_t *dbound_edge     = NULL;
  int n_group_face = PDM_dmesh_bound_get(dmesh,
                                         PDM_BOUND_TYPE_FACE,
                                         &dbound_face,
                                         &dbound_face_idx,
                                         PDM_OWNERSHIP_KEEP);

  int n_group_edge = PDM_dmesh_bound_get(dmesh,
                                         PDM_BOUND_TYPE_EDGE,
                                         &dbound_edge,
                                         &dbound_edge_idx,
                                         PDM_OWNERSHIP_KEEP);

  PDM_dmesh_to_dmesh_nodal_group_set(dm_to_dmn,
                                     0,
                                     PDM_BOUND_TYPE_FACE,
                                     n_group_face,
                                     dbound_face_idx,
                                     dbound_face);

  PDM_dmesh_to_dmesh_nodal_group_set(dm_to_dmn,
                                     0,
                                     PDM_BOUND_TYPE_EDGE,
                                     n_group_edge,
                                     dbound_edge_idx,
                                     dbound_edge);

  PDM_dmesh_to_dmesh_nodal_compute(dm_to_dmn);

  PDM_dmesh_nodal_t* dmn_out = NULL;
  PDM_dmesh_to_dmesh_nodal_dmesh_nodal_get(dm_to_dmn,
                                           0,
                                           &dmn_out,
                                           PDM_OWNERSHIP_KEEP);

  int dn_vtx2 = PDM_DMesh_nodal_n_vtx_get(dmn);

  PDM_DMesh_nodal_coord_set(dmn_out,
                            dn_vtx2,
                            dvtx_coords,
                            PDM_OWNERSHIP_USER);


  if(0 == 1) {
    if(dim == 3) {
      PDM_dmesh_nodal_dump_vtk(dmn_out, PDM_GEOMETRY_KIND_VOLUMIC , "out_post_volumic");
    }
    PDM_dmesh_nodal_dump_vtk(dmn_out, PDM_GEOMETRY_KIND_SURFACIC, "out_post_surfacic");
    PDM_dmesh_nodal_dump_vtk(dmn_out, PDM_GEOMETRY_KIND_RIDGE   , "out_post_ridge");
  }


  PDM_dmesh_to_dmesh_nodal_free(dm_to_dmn);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
