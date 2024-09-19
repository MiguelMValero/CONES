/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_poly_surf_gen.h"
#include "pdm_part.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_part_extension.h"
#include "pdm_vtk.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"


#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MIN(a,b) ((a) < (b) ?  (a) : (b))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

typedef enum {
    CWP_VERSION_OLD,
    CWP_VERSION_NEW
} CWP_Version_t;


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n           <> Number of vertices in band length.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                   argc,
  char                **argv,
  CWP_Version_t        *version,
  int                  *n_vtx_seg1,
  int                  *n_vtx_seg2,
  double               *length,
  double               *separation_x,
  double               *separation_y,
  double               *separation_z,
  int                  *deform,
  double               *tolerance,
  int                  *randomize,
  int                  *nProcData,
  PDM_split_dual_t     *part_method,
  CWP_Spatial_interp_t *loc_method,
  char                **output_filename,
  int                  *verbose,
  int                  *extension_depth_tgt,
  int                  *extension_depth_src,
  PDM_Mesh_nodal_elt_t *elt_type,
  int                  *use_gnum,
  int                  *interlaced_field
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-new") == 0) {
      *version = CWP_VERSION_NEW;
    }
    else if (strcmp(argv[i], "-old") == 0) {
      *version = CWP_VERSION_OLD;
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-length") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-output") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *output_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nProcData = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-ext_depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_tgt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_tgt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ext_depth_src") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *extension_depth_src = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-no_gnum") == 0) {
      *use_gnum = 0;
    }
    else if (strcmp(argv[i], "-interlaced") == 0) {
      *interlaced_field = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static double R[3][3] =
{
  {-0.14547275709949994,  0.8415293589391187 , -0.5202557207618055 },
  { 0.9893622576902102 ,  0.12373586628506748, -0.07649678720582984},
  { 0.                 , -0.5258495730132333 , -0.8505775840931856 }
};

static void
_rotate(const int n_pts, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[j][0] * x + R[j][1] * y + R[j][2] * z;
    }
  }
}

static void
_unrotate(const int n_pts, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[0][j] * x + R[1][j] * y + R[2][j] * z;
    }
  }
}



static void
_cube_mesh
(
 const int                     activeRank,
 const PDM_MPI_Comm            comm,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const PDM_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     deform,
 const int                     part_extension_depth,
 const PDM_Mesh_nodal_elt_t    elt_type,
 int                         **pn_cell,
 int                         **pn_face,
 int                         **pn_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 int                        ***pcell_vtx,
 int                        ***pface_vtx_idx,
 int                        ***pface_vtx,
 double                     ***pvtx_coord,
 PDM_g_num_t                ***pcell_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pvtx_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  if (activeRank) {

    // Utiliser pdm_t_part_dcube_nodal_cache_blocking !!!!!
    //          pdm_partitioning_nodal_algorith


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
                                                           PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build (dcube);


    PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];


    PDM_dmesh_nodal_to_dmesh_t *dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

    PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

    PDM_dmesh_t *dmesh2 = NULL;
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh2);

    int dn_cell = 0;
    int dn_face = 0;
    int dn_edge = -1;
    int n_face_group = 0;

    int         *dcell_face_idx  = NULL;
    PDM_g_num_t *dcell_face      = NULL;
    int         *dface_vtx_idx   = NULL;
    PDM_g_num_t *dface_vtx       = NULL;
    int         *dface_cell_idx  = NULL;
    PDM_g_num_t *dface_cell      = NULL;
    PDM_g_num_t *dface_group     = NULL;
    int         *dface_group_idx = NULL;

    dn_face = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         &dface_vtx,
                                         &dface_vtx_idx,
                                         PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               &dface_cell,
                               &dface_cell_idx,
                               PDM_OWNERSHIP_KEEP);
    assert(dface_cell_idx == NULL);

    dn_cell = PDM_dmesh_connectivity_get(dmesh2, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                         &dcell_face,
                                         &dcell_face_idx,
                                         PDM_OWNERSHIP_KEEP);

    n_face_group = PDM_dmesh_bound_get(dmesh2,
                                       PDM_BOUND_TYPE_FACE,
                                       &dface_group,
                                       &dface_group_idx,
                                       PDM_OWNERSHIP_KEEP);



    if (deform) {
      /*for (int i = 0; i < dn_vtx; i++) {
        double x = dvtx_coord[3*i];
        double z = dvtx_coord[3*i + 2];

        dvtx_coord[3*i]     += 0.1 * z * z;
        dvtx_coord[3*i + 2] += 0.2 * cos(PDM_PI * x);
        }*/
      _rotate (dn_vtx,
               dvtx_coord);
    }

    /*
     *  Create mesh partitiions
     */

    /* Initialize multipart */
    PDM_multipart_t *mpart = PDM_multipart_create(1,
                                                  &n_part,
                                                  PDM_FALSE,
                                                  part_method,
                                                  PDM_PART_SIZE_HOMOGENEOUS,
                                                  NULL,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                           dn_cell,
                                           dn_face,
                                           dn_edge,
                                           dn_vtx,
                                           comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_USER);

    PDM_multipart_dmesh_set (mpart, 0, dmesh);

    /* Run */
    PDM_multipart_compute (mpart);


    // PDM_dcube_gen_free(dcube);


    PDM_part_extension_t *part_ext = NULL;


    if (part_extension_depth > 0) {

      // TODO: Part extension et nodal sont-ils compatibles ?
      // assert(0);


      part_ext = PDM_part_extension_create(1,
                                           &n_part,
                                           PDM_EXTEND_FROM_FACE,
                                           part_extension_depth,
                                           comm,
                                           PDM_OWNERSHIP_KEEP);


      for (int i_part = 0; i_part < n_part; i_part++) {

        PDM_g_num_t* cell_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        &cell_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        int *cell_face     = NULL;
        int *cell_face_idx = NULL;
        int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                         &cell_face_idx,
                                                         &cell_face,
                                                         PDM_OWNERSHIP_KEEP);


        int *face_cell     = NULL;
        int *face_cell_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                            &face_cell_idx,
                                            &face_cell,
                                            PDM_OWNERSHIP_KEEP);
        assert(face_cell_idx == NULL);

        int *face_vtx     = NULL;
        int *face_vtx_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                            &face_vtx_idx,
                                            &face_vtx,
                                            PDM_OWNERSHIP_KEEP);

        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                         &face_edge_idx,
                                                         &face_edge,
                                                         PDM_OWNERSHIP_KEEP);

        PDM_g_num_t* face_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                         &edge_vtx_idx,
                                                         &edge_vtx,
                                                         PDM_OWNERSHIP_KEEP);
        assert(edge_vtx_idx == NULL);
        PDM_g_num_t* edge_ln_to_gn = NULL;
        int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      0,
                                                      i_part,
                                                      PDM_MESH_ENTITY_EDGE,
                                                      &edge_ln_to_gn,
                                                      PDM_OWNERSHIP_KEEP);
        assert(n_edge2 == n_edge);

        PDM_g_num_t* vtx_ln_to_gn = NULL;
        int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    0,
                                                    i_part,
                                                    PDM_MESH_ENTITY_VTX,
                                                    &vtx_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);
        int          pn_face_group        = 0;
        int         *face_group_idx      = NULL;
        int         *face_group          = NULL;
        PDM_g_num_t *face_group_ln_to_gn = NULL;
        PDM_multipart_group_get(mpart,
                                0,
                                i_part,
                                PDM_MESH_ENTITY_FACE,
                                &pn_face_group,
                                &face_group_idx,
                                &face_group,
                                &face_group_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);

        int *vtx_part_bound_proc_idx = NULL;
        int *vtx_part_bound_part_idx = NULL;
        int *vtx_part_bound          = NULL;
        PDM_multipart_part_graph_comm_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_VTX,
                                          &vtx_part_bound_proc_idx,
                                          &vtx_part_bound_part_idx,
                                          &vtx_part_bound,
                                          PDM_OWNERSHIP_KEEP);

        int n_face_part_bound = 0;
        int *face_part_bound_proc_idx = NULL;
        int *face_part_bound_part_idx = NULL;
        int *face_part_bound          = NULL;
        PDM_multipart_part_graph_comm_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_FACE,
                                          &face_part_bound_proc_idx,
                                          &face_part_bound_part_idx,
                                          &face_part_bound,
                                          PDM_OWNERSHIP_KEEP);
        double *vtx = NULL;
        PDM_multipart_part_vtx_coord_get(mpart,
                                         0,
                                         i_part,
                                         &vtx,
                                         PDM_OWNERSHIP_KEEP);

        PDM_part_extension_set_part(part_ext,
                                    0,
                                    i_part,
                                    n_cell,
                                    n_face,
                                    n_face_part_bound,
                                    n_face_group,
                                    0,   // n_edge
                                    n_vtx,
                                    cell_face_idx,
                                    cell_face,
                                    face_cell,
                                    NULL, // face_edge_idx
                                    NULL, // face_edge
                                    face_vtx_idx,
                                    face_vtx,
                                    NULL, //edge_vtx
                                    face_group_idx,
                                    face_group,
                                    NULL, // face_join_idx
                                    NULL, // face_join
                                    face_part_bound_proc_idx,
                                    face_part_bound_part_idx,
                                    face_part_bound,
                                    NULL, // vtx_part_bound_proc_idx
                                    NULL, // vtx_part_bound_part_idx
                                    NULL, // vtx_part_bound
                                    cell_ln_to_gn,
                                    face_ln_to_gn,
                                    NULL, // edge_ln_to_gn
                                    vtx_ln_to_gn,
                                    face_group_ln_to_gn,
                                    vtx);
      }

      PDM_part_extension_compute(part_ext);
    }


    *pn_cell        = (int *)          malloc(sizeof(int *)          * n_part);
    *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
    *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
    *pcell_face_idx = (int **)         malloc(sizeof(int **)         * n_part);
    *pcell_face     = (int **)         malloc(sizeof(int **)         * n_part);
    *pcell_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
    *pface_vtx_idx  = (int **)         malloc(sizeof(int **)         * n_part);
    *pface_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
    *pcell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
    *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
    *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
    *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face_idx,
                                                       &cell_face,
                                                       PDM_OWNERSHIP_USER);

      PDM_g_num_t *cell_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_USER);

      int *face_vtx_idx  = NULL;
      int *face_vtx      = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX ,
                                                       &face_vtx_idx,
                                                       &face_vtx,
                                                       PDM_OWNERSHIP_USER);

      PDM_g_num_t *face_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_USER);


      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_USER);

      PDM_g_num_t *vtx_ln_to_gn  = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_USER);

      int n_ext_vtx  = 0;
      int n_ext_face = 0;
      int n_ext_cell = 0;
      double      *ext_vtx_coord     = NULL;
      PDM_g_num_t *ext_vtx_ln_to_gn  = NULL;
      int         *ext_cell_face     = NULL;
      int         *ext_cell_face_idx = NULL;
      PDM_g_num_t *ext_cell_ln_to_gn = NULL;
      int         *ext_face_vtx      = NULL;
      int         *ext_face_vtx_idx  = NULL;
      PDM_g_num_t *ext_face_ln_to_gn = NULL;

      if (part_extension_depth > 0) {
        /* Vertices */
        n_ext_vtx = PDM_part_extension_vtx_coord_get(part_ext,
                                                     0,
                                                     i_part,
                                                     &ext_vtx_coord);

        PDM_part_extension_ln_to_gn_get(part_ext,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &ext_vtx_ln_to_gn);


        /* Cells */
        n_ext_cell = PDM_part_extension_connectivity_get(part_ext,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                         &ext_cell_face,
                                                         &ext_cell_face_idx);

        PDM_part_extension_ln_to_gn_get(part_ext,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        &ext_cell_ln_to_gn);


        /* Faces */
        n_ext_face = PDM_part_extension_connectivity_get(part_ext,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                         &ext_face_vtx,
                                                         &ext_face_vtx_idx);
        PDM_part_extension_ln_to_gn_get(part_ext,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &ext_face_ln_to_gn);
      }


      *(pn_cell)[i_part] = n_cell + n_ext_cell;
      *(pn_face)[i_part] = n_face + n_ext_face;
      *(pn_vtx)[i_part]  = n_vtx  + n_ext_vtx;

      /* Vertices */
      (*pvtx_coord)[i_part] = vtx_coord;
      (*pvtx_ln_to_gn)[i_part] = vtx_ln_to_gn;


      /* Cells */
      (*pcell_face_idx)[i_part] = cell_face_idx;

      // int s_cell_face = cell_face_idx[n_cell];
      // if (part_extension_depth > 0) {
      //   s_cell_face += ext_cell_face_idx[n_ext_cell];
      // }
      (*pcell_face)[i_part] = cell_face;
      (*pcell_ln_to_gn)[i_part] = cell_ln_to_gn;


      /* Faces */
      (*pface_vtx_idx)[i_part] = face_vtx_idx;

      // int s_face_vtx = face_vtx_idx[n_face];
      // if (part_extension_depth > 0) {
      //   s_face_vtx += ext_face_vtx_idx[n_ext_face];
      // }
      (*pface_vtx)[i_part] = face_vtx;

      (*pface_ln_to_gn)[i_part] = face_ln_to_gn;


      if (part_extension_depth > 0) {
        /* Vertices */
        memcpy((*pvtx_coord)[i_part] + 3*n_vtx, ext_vtx_coord, sizeof(double) * 3 * n_ext_vtx);
        memcpy((*pvtx_ln_to_gn)[i_part] + n_vtx, ext_vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_vtx);

        /* Cells */
        for (int i = 1; i <= n_ext_cell; i++) {
          (*pcell_face_idx)[i_part][n_cell + i] = cell_face_idx[n_cell] + ext_cell_face_idx[i];
        }

        memcpy((*pcell_face)[i_part] + cell_face_idx[n_cell],
               ext_cell_face,
               sizeof(int) * ext_cell_face_idx[n_ext_cell]);

        memcpy((*pcell_ln_to_gn)[i_part] + n_cell, ext_cell_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_cell);

        /* Faces */
        for (int i = 1; i <= n_ext_face; i++) {
          (*pface_vtx_idx)[i_part][n_face + i] = face_vtx_idx[n_face] + ext_face_vtx_idx[i];
        }

        memcpy((*pface_vtx)[i_part] + face_vtx_idx[n_face],
               ext_face_vtx,
               sizeof(int) * ext_face_vtx_idx[n_ext_face]);

        memcpy((*pface_ln_to_gn)[i_part] + n_face, ext_face_ln_to_gn, sizeof(PDM_g_num_t) * n_ext_face);
      }

    }

    PDM_Mesh_nodal_t *nodal = PDM_Mesh_nodal_create (n_part, comm);
    PDM_l_num_t      **face_vtx_nb  = malloc (sizeof(int*) * n_part);
    PDM_l_num_t      **cell_face_nb = malloc (sizeof(int*) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_Mesh_nodal_coord_set (nodal,
                                i_part,
                                (*pn_vtx)[i_part],
                                (*pvtx_coord)[i_part],
                                (*pvtx_ln_to_gn)[i_part],
                                PDM_OWNERSHIP_USER);

      face_vtx_nb[i_part]  = malloc (sizeof(int) * (*pn_face)[i_part]);
      cell_face_nb[i_part] = malloc (sizeof(int) * (*pn_cell)[i_part]);

      for (int i = 0; i < (*pn_cell)[i_part]; i++) {
        cell_face_nb[i_part][i] = (*pcell_face_idx)[i_part][i + 1] - (*pcell_face_idx)[i_part][i];
      }

      for (int i = 0; i < (*pn_face)[i_part]; i++) {
        face_vtx_nb[i_part][i] = (*pface_vtx_idx)[i_part][i + 1] - (*pface_vtx_idx)[i_part][i];
      }

      PDM_Mesh_nodal_cell3d_cellface_add (nodal,
                                          i_part,
                                          (*pn_cell)[i_part],
                                          (*pn_face)[i_part],
                                          (*pface_vtx_idx)[i_part],
                                          face_vtx_nb[i_part],
                                          (*pface_vtx)[i_part],
                                          NULL,
                                          (*pcell_face_idx)[i_part],
                                          cell_face_nb[i_part],
                                          (*pcell_face)[i_part],
                                          (*pcell_ln_to_gn)[i_part],
                                          PDM_OWNERSHIP_USER);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      free (cell_face_nb[i_part]);
      free (face_vtx_nb[i_part]);
    }

    free (cell_face_nb);
    free (face_vtx_nb);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_Mesh_nodal_block_std_get (nodal,
                                    0,
                                    i_part,
                                    &((*pcell_vtx)[i_part]));

      int *parent_num = PDM_Mesh_nodal_block_parent_num_get(nodal,
                                                            0,
                                                            i_part);
      free(parent_num);

      PDM_g_num_t *numabs = PDM_Mesh_nodal_g_num_get(nodal,
                                                     0,
                                                     i_part);
      free(numabs);
    }

    PDM_Mesh_nodal_free(nodal);

    PDM_multipart_free (mpart);
    PDM_dmesh_free (dmesh);
    PDM_part_extension_free (part_ext);

    PDM_dmesh_nodal_to_dmesh_free(dmntodm);
    PDM_dcube_nodal_gen_free(dcube);
  }
  else {
    *pn_cell = (int *) malloc (sizeof(int *) * n_part);
    *pn_face = (int *) malloc (sizeof(int *) * n_part);
    *pn_vtx = (int *) malloc (sizeof(int *) * n_part);
    *pcell_face_idx = (int **) malloc (sizeof(int *) * n_part);
    *pcell_face = (int **) malloc (sizeof(int *) * n_part);
    *pcell_vtx = (int **) malloc (sizeof(int *) * n_part);
    *pface_vtx_idx = (int **) malloc (sizeof(int *) * n_part);
    *pface_vtx = (int **) malloc (sizeof(int *) * n_part);
    *pvtx_coord = (double **) malloc (sizeof(int *) * n_part);
    *pcell_ln_to_gn = (PDM_g_num_t **) malloc (sizeof(int *) * n_part);
    *pface_ln_to_gn = (PDM_g_num_t **) malloc (sizeof(int *) * n_part);
    *pvtx_ln_to_gn = (PDM_g_num_t **) malloc (sizeof(int *) * n_part);
    for (int ipart = 0 ; ipart < n_part ; ipart++) {
      int _nFace = 0;
      int _nCell = 0;
      int _nVtx = 0;
      (*pn_cell)[ipart] = _nCell;
      (*pn_face)[ipart] = _nFace;
      (*pn_vtx)[ipart] = _nVtx;
      (*pcell_face_idx)[ipart] = (int *) malloc(sizeof(int) * (_nCell + 1));
      (*pcell_face_idx)[ipart][0] = 0;
      (*pcell_face)[ipart] = (int *) malloc(sizeof(int) * 0);
      (*pcell_vtx)[ipart] = (int *) malloc(sizeof(int) * 0);
      (*pface_vtx_idx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*pface_vtx_idx)[ipart][0] = 0;
      (*pface_vtx)[ipart] = (int *) malloc(sizeof(int) * 0);
      (*pvtx_coord)[ipart] = (double *) malloc(sizeof(double) * 0);
      (*pcell_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nCell);
      (*pface_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);
      (*pvtx_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
    }
  }
//  printf("pncel, pnface, pnvtx : %d %d %d \n", (*pn_cell)[0], (*pn_face)[0], (*pn_vtx)[0]);
}


static int
_set_rank_has_mesh(const MPI_Comm comm, const int nProcData, PDM_MPI_Comm *meshComm) {
  int current_rank_has_mesh = 1;
  int rank;
  int commSize;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commSize);

  PDM_MPI_Comm _comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &comm);

  if (nProcData > 0 && nProcData < commSize) {
    int rankInNode = PDM_io_mpi_node_rank(_comm);

    int nNode = 0;
    int iNode = -1;
    int masterRank = (rankInNode == 0);

    int *rankInNodes = malloc(sizeof(int) * commSize);

    MPI_Allreduce(&masterRank, &nNode, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allgather(&rankInNode, 1, MPI_INT, rankInNodes, 1, MPI_INT, comm);

    current_rank_has_mesh = 0;

    for (int i = 0 ; i < rank ; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && masterRank) {
        current_rank_has_mesh = 1;
      }
    }
    else {
      if (rankInNode < (nProcData / nNode)) {
        current_rank_has_mesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        current_rank_has_mesh = 1;
      }
    }

    PDM_MPI_Comm_split(_comm, current_rank_has_mesh, rank, meshComm);
    free(rankInNodes);
  }

  return current_rank_has_mesh;
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  // Read args from command line
  CWP_Version_t version           = CWP_VERSION_NEW;
  int n_vtx_seg1                  = 4;
  int n_vtx_seg2                  = 4;
  int randomize                   = 1;
  int n_proc_data                 = -1;

  CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif
  int verbose                     = 0;

  double length                   = 20.;
  int deform                      = 0;

  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;

  double      tolerance           = 1e-2;

  char* output_filename           = NULL;
  int filedump                    = 0;

  int         extension_depth_tgt = 0;
  int         extension_depth_src = 0;
  int         use_gnum            = 1;
  int         interlaced_field    = 0;

  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_TETRA4;

  _read_args (argc,
              argv,
             &version,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &deform,
             &tolerance,
             &randomize,
             &n_proc_data,
             &part_method,
             &loc_method,
             &output_filename,
             &verbose,
             &extension_depth_tgt,
             &extension_depth_src,
             &elt_type,
             &use_gnum,
             &interlaced_field);

  if (output_filename !=NULL) {
    filedump = 1;
  }

  CWP_Field_storage_t field_storage = CWP_FIELD_STORAGE_INTERLEAVED;
  if (interlaced_field) {
    field_storage = CWP_FIELD_STORAGE_INTERLACED;
  }


  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }

  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  if (version == CWP_VERSION_OLD) {
    cwipi_init(MPI_COMM_WORLD, code_name[0], intra_comm);
  }

  else {
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);
  }

  if (verbose && rank == 0) {
    printf("CWIPI Init OK\n");
  }

  // Create coupling
  const char *coupling_name = "c_new_vs_old_vol";

  if (version == CWP_VERSION_OLD) {
    cwipi_create_coupling(coupling_name,                             // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          coupled_code_name[0],                      // Coupled application id
                          3,                                         // Geometric entities dimension
                          tolerance,                                 // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          CWIPI_SOLVER_CELL_VERTEX,                  // Solver type
                          -1,                                        // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");
  }
  else {
    CWP_Cpl_create(code_name[0],
                   coupling_name,
                   coupled_code_name[0],
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    // CWP_Visu_set (code_name[0],
    //               coupling_name,
    //               1,
    //               CWP_VISU_FORMAT_ENSIGHT,
    //               "text");

  }

  if (verbose) {
    printf("Create coupling OK\n");
  }

  // Define mesh
  PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

  int _n_proc_data = n_proc_data;
  if (n_proc_data > 0) {
    if (code_id == 1) {
      _n_proc_data /= 2;
    }
    else {
      _n_proc_data -= n_proc_data / 2;
    }
  }
  int current_rank_has_mesh = _set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  int true_n_proc_data;
  MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("nb procs with mesh data = %d\n", true_n_proc_data);
  }

  double xmin = -0.5 * length;
  double ymin = -0.5 * length;
  double zmin = -0.5 * length;

  double xyz_min[3] = {xmin, ymin, zmin};
  double xyz_max[3] = {xmin + length, ymin + length, zmin + length};

//  int init_random = (int) time(NULL);
  int init_random = 5;

  PDM_MPI_Comm code_mesh_comm;
  PDM_MPI_Comm_split(mesh_comm, code_id, rank, &code_mesh_comm);

  if (code_id == 2) {
    init_random++;
    xmin += separation_x;
    ymin += separation_y;
    zmin += separation_z;
  }

  int                         *pn_cell;
  int                         *pn_face;
  int                         *pn_vtx;
  int                        **pcell_face_idx;
  int                        **pcell_face;
  int                        **pcell_vtx;
  int                        **pface_vtx_idx;
  int                        **pface_vtx;
  double                     **pvtx_coord;
  PDM_g_num_t                **pcell_ln_to_gn;
  PDM_g_num_t                **pface_ln_to_gn;
  PDM_g_num_t                **pvtx_ln_to_gn;

  _cube_mesh (current_rank_has_mesh,
              code_mesh_comm,
              1,
              part_method,
              n_vtx_seg,
              xmin,
              ymin,
              zmin,
              length,
              deform,
              extension_depth_src,
              elt_type,
             &pn_cell,
             &pn_face,
             &pn_vtx,
             &pcell_face_idx,
             &pcell_face,
             &pcell_vtx,
             &pface_vtx_idx,
             &pface_vtx,
             &pvtx_coord,
             &pcell_ln_to_gn,
             &pface_ln_to_gn,
             &pvtx_ln_to_gn);

  int *cellVtxIdx = malloc(sizeof(int) * (pn_cell[0] + 1));
  cellVtxIdx[0] = 0;

  CWP_Block_t block_type;

  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    block_type = CWP_BLOCK_CELL_TETRA4;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 4 + cellVtxIdx[i];
    }
  }

  else if (elt_type == PDM_MESH_NODAL_HEXA8) {
    block_type = CWP_BLOCK_CELL_HEXA8;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 8 + cellVtxIdx[i];
    }
  }

  else if (elt_type == PDM_MESH_NODAL_PYRAMID5) {
    block_type = CWP_BLOCK_CELL_PYRAM5;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 5 + cellVtxIdx[i];
    }
  }

  else if (elt_type == PDM_MESH_NODAL_PRISM6) {
    block_type = CWP_BLOCK_CELL_PRISM6;
    for (int i = 0; i < pn_cell[0]; i++) {
      cellVtxIdx[i+1] = 6 + cellVtxIdx[i];
    }
  }

  // Set interface mesh
  if (version == CWP_VERSION_OLD) {
    cwipi_define_mesh(coupling_name, pn_vtx[0], pn_cell[0], pvtx_coord[0], cellVtxIdx, pcell_vtx[0]);
  }
  else {
    PDM_g_num_t *_pvtx_ln_to_gn  = NULL;
    PDM_g_num_t *_pcell_ln_to_gn = NULL;
    if (use_gnum) {
      _pvtx_ln_to_gn  = pvtx_ln_to_gn[0];
      _pcell_ln_to_gn = pcell_ln_to_gn[0];
    }
    CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, pn_vtx[0], pvtx_coord[0], _pvtx_ln_to_gn);

    int block_id = CWP_Mesh_interf_block_add (code_name[0],
                                              coupling_name,
                                              block_type);


    CWP_Mesh_interf_block_std_set (code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   pn_cell[0],
                                   pcell_vtx[0],
                                   _pcell_ln_to_gn);


    CWP_Mesh_interf_finalize(code_name[0], coupling_name);
  }

  if (verbose && rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  double *send_val  = NULL;
  double *send_val2 = NULL;
  double *recv_val  = NULL;
  double *recv_val2 = NULL;

  const char *field_name  = "coo";
  const char *field_name2 = "coocooY";

  if (code_id == 1) {
    send_val  = (double *) malloc(sizeof(double) * pn_vtx[0] * 3);
    send_val2 = (double *) malloc(sizeof(double) * pn_vtx[0]);

    for (int i = 0 ; i < pn_vtx[0]; i++) {
      send_val2[i] = pvtx_coord[0][3*i+1];
    }

    if (interlaced_field || version == CWP_VERSION_OLD) {
      for (int i = 0 ; i < 3*pn_vtx[0]; i++) {
        send_val[i] = pvtx_coord[0][i];
      }

    }
    else {
      for (int i = 0 ; i < pn_vtx[0]; i++) {
        for (int j = 0; j < 3; j++) {
          send_val[pn_vtx[0]*j + i] = pvtx_coord[0][3*i + j];
        }
      }
    }

    if (0) {
      log_trace("send_val = \n");
      if (interlaced_field || version == CWP_VERSION_OLD) {
        for (int i = 0 ; i < pn_vtx[0]; i++) {
          log_trace("%f %f %f\n", send_val[3*i], send_val[3*i+1], send_val[3*i+2]);
        }
      } else {
        for (int j = 0; j < 3; j++) {
          for (int i = 0 ; i < pn_vtx[0]; i++) {
            log_trace("%f ", send_val[pn_vtx[0]*j + i]);
          }
          log_trace("\n");
        }
      }
    }
  }
  else {
    recv_val  = (double *) malloc(sizeof(double) * pn_vtx[0] * 3);
    recv_val2 = (double *) malloc(sizeof(double) * pn_vtx[0]);
  }

  if (version == CWP_VERSION_NEW) {
    CWP_Status_t visu_status = CWP_STATUS_ON;
    MPI_Barrier(MPI_COMM_WORLD);

    if (code_id == 1) {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       field_storage,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       field_storage,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val2);
    }
    else {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       field_storage,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       field_storage,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val2);
    }

    CWP_Time_step_beg(code_name[0],
                      0.0);
  }

  if (verbose && rank == 0) {
    printf("Fields OK\n");
  }

  // Perform geometric algorithm
  PDM_timer_t *timer = PDM_timer_create();
  double t_start, t_end;

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_init(timer);

  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  int n_unlocated = 0;
  int n_located = 0;
  const int *located = NULL;
  if (version == CWP_VERSION_OLD) {
    cwipi_locate(coupling_name);

    if (code_id != 1) {
      n_unlocated = cwipi_get_n_not_located_points(coupling_name);
      n_located = cwipi_get_n_located_points(coupling_name);
      located = cwipi_get_located_points(coupling_name);
    }

  }
  else {
    PDM_part_to_block_global_statistic_reset();
    PDM_block_to_part_global_statistic_reset();

    char char_tol[99];
    sprintf(char_tol, "%e", tolerance);
    CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", CWP_DOUBLE, char_tol);
    CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

    if (code_id != 1) {
      n_unlocated = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      n_located = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      located =  CWP_Computed_tgts_get(code_name[0], coupling_name, field_name2, 0);
    }

    double min_elaps_create_ptb;
    double max_elaps_create_ptb;
    double min_cpu_create_ptb;
    double max_cpu_create_ptb;
    double min_elaps_create2_ptb;
    double max_elaps_create2_ptb;
    double min_cpu_create2_ptb;
    double max_cpu_create2_ptb;
    double min_elaps_exch_ptb;
    double max_elaps_exch_ptb;
    double min_cpu_exch_ptb;
    double max_cpu_exch_ptb;

    PDM_part_to_block_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_ptb,
                                       &max_elaps_create_ptb,
                                       &min_cpu_create_ptb,
                                       &max_cpu_create_ptb,
                                       &min_elaps_create2_ptb,
                                       &max_elaps_create2_ptb,
                                       &min_cpu_create2_ptb,
                                       &max_cpu_create2_ptb,
                                       &min_elaps_exch_ptb,
                                       &max_elaps_exch_ptb,
                                       &min_cpu_exch_ptb,
                                       &max_cpu_exch_ptb);

    double min_elaps_create_btp;
    double max_elaps_create_btp;
    double min_cpu_create_btp;
    double max_cpu_create_btp;
    double min_elaps_exch_btp;
    double max_elaps_exch_btp;
    double min_cpu_exch_btp;
    double max_cpu_exch_btp;

    PDM_block_to_part_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_btp,
                                       &max_elaps_create_btp,
                                       &min_cpu_create_btp,
                                       &max_cpu_create_btp,
                                       &min_elaps_exch_btp,
                                       &max_elaps_exch_btp,
                                       &min_cpu_exch_btp,
                                       &max_cpu_exch_btp);

    int min_exch_rank_send_ptb;
    int min_exch_rank_recv_ptb;
    int max_exch_rank_send_ptb;
    int max_exch_rank_recv_ptb;
    unsigned long long min_exch_data_send_ptb;
    unsigned long long min_exch_data_recv_ptb;
    unsigned long long max_exch_data_send_ptb;
    unsigned long long max_exch_data_recv_ptb;

    PDM_part_to_block_global_statistic_get(PDM_MPI_COMM_WORLD,
                                           &min_exch_rank_send_ptb,
                                           &min_exch_rank_recv_ptb,
                                           &max_exch_rank_send_ptb,
                                           &max_exch_rank_recv_ptb,
                                           &min_exch_data_send_ptb,
                                           &min_exch_data_recv_ptb,
                                           &max_exch_data_send_ptb,
                                           &max_exch_data_recv_ptb);

    int min_exch_rank_send_btp;
    int min_exch_rank_recv_btp;
    int max_exch_rank_send_btp;
    int max_exch_rank_recv_btp;
    unsigned long long min_exch_data_send_btp;
    unsigned long long min_exch_data_recv_btp;
    unsigned long long max_exch_data_send_btp;
    unsigned long long max_exch_data_recv_btp;

    PDM_block_to_part_global_statistic_get(PDM_MPI_COMM_WORLD,
                                           &min_exch_rank_send_btp,
                                           &min_exch_rank_recv_btp,
                                           &max_exch_rank_send_btp,
                                           &max_exch_rank_recv_btp,
                                           &min_exch_data_send_btp,
                                           &min_exch_data_recv_btp,
                                           &max_exch_data_send_btp,
                                           &max_exch_data_recv_btp);

    if (rank == 0) {
      printf("Global time in PDM_part_to_block : \n");
      printf("   - ptb min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_ptb,
             max_elaps_create_ptb);
      printf("   - ptb min max elaps create2 : %12.5e %12.5e\n",
             min_elaps_create2_ptb,
             max_elaps_create2_ptb);
      printf("   - ptb min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_ptb,
             max_elaps_exch_ptb);

      printf("   - ptb min max send rank     : %d %d\n",
             min_exch_rank_send_ptb,
             max_exch_rank_send_ptb);
      printf("   - ptb min max recv rank     : %d %d\n",
             min_exch_rank_recv_ptb,
             max_exch_rank_recv_ptb);
      printf("   - ptb min max send data     : %llu %llu\n",
             min_exch_data_send_ptb,
             max_exch_data_send_ptb);
      printf("   - ptb min max recv data     : %llu %llu\n",
             min_exch_data_recv_ptb,
             max_exch_data_recv_ptb);
      fflush(stdout);

      printf("Global time in PDM_block_to_part : \n");
      printf("   - btp min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_btp,
             max_elaps_create_btp);
      printf("   - btp min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_btp,
             max_elaps_exch_btp);

      printf("   - btp min max send rank     : %d %d\n",
             min_exch_rank_send_btp,
             max_exch_rank_send_btp);
      printf("   - btp min max recv rank     : %d %d\n",
             min_exch_rank_recv_btp,
             max_exch_rank_recv_btp);
      printf("   - btp min max send data     : %llu %llu\n",
             min_exch_data_send_btp,
             max_exch_data_send_btp);
      printf("   - btp min max recv data     : %llu %llu\n",
             min_exch_data_recv_btp,
             max_exch_data_recv_btp);
      fflush(stdout);

      char output_filename_ptbbtp[50];
      FILE* output_ptbbtp = stdout;
      if (filedump) {
        sprintf(output_filename_ptbbtp, "ptbbtp_%s", output_filename);
        output_ptbbtp = fopen(output_filename_ptbbtp, "w");
      }

      // fprintf(output_ptbbtp,
      //         "%d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
      //         comm_world_size,
      //         min_elaps_create_ptb,
      //         max_elaps_create_ptb,
      //         min_elaps_create2_ptb,
      //         max_elaps_create2_ptb,
      //         min_elaps_exch_ptb,
      //         max_elaps_exch_ptb);


//      fprintf(output_ptbbtp,
//              "%d %d %d %d %llu %llu %llu %llu ",
//              min_exch_rank_send_ptb,
//              max_exch_rank_send_ptb,
//              min_exch_rank_recv_ptb,
//              max_exch_rank_recv_ptb,
//              min_exch_data_send_ptb,
//              max_exch_data_send_ptb,
//              min_exch_data_recv_ptb,
//              max_exch_data_recv_ptb);

      // fprintf(output_ptbbtp,
      //       "%12.5e %12.5e %12.5e %12.5e\n",
      //       min_elaps_create_btp,
      //       max_elaps_create_btp,
      //       min_elaps_exch_btp,
      //       max_elaps_exch_btp);

//      fprintf(output_ptbbtp,
//              "%d %d %d %d %llu %llu %llu %llu ",
//              min_exch_rank_send_btp,
//              max_exch_rank_send_btp,
//              min_exch_rank_recv_btp,
//              max_exch_rank_recv_btp,
//              min_exch_data_send_btp,
//              max_exch_data_send_btp,
//              min_exch_data_recv_btp,
//              max_exch_data_recv_btp);

      if (filedump) {
        fclose(output_ptbbtp);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  double geom_time = t_end - t_start;
  double max_geom_time;
  MPI_Reduce(&geom_time, &max_geom_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  FILE *output = stdout;
  if (filedump) {
    output = fopen(output_filename, "w");
  }

  if (rank == 0) {
    printf("\n\nGeometric algorithm :%12.5es\n", max_geom_time);

    if (filedump) {
      fprintf(output, "%d %12.5e ", comm_world_size, max_geom_time);
    }
  }

  if (verbose && rank == 0) {
    printf("Geometric algorithm OK\n");

  }

  printf("n_unlocated : %d\n", n_unlocated);

  //  Exchange interpolated fields 1
  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  int request;
  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 3, 1, 0.1, field_name, send_val, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 3, 1, 0.1, field_name, recv_val, &request);
    }
  }

  else {
    if (code_id == 1) {
      CWP_Field_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_irecv(code_name[0], coupling_name, field_name);
    }
  }

  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* PDM_timer_hang_on(timer); */
  /* t_end = PDM_timer_elapsed(timer); */
  /* double exch_time1 = t_end - t_start; */
  /* double max_exch_time1; */
  /* t_start = t_end; */
  /* MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); */
  /* PDM_timer_resume(timer); */

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
      CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  double exch_time = t_end - t_start;
  double max_exch_time;
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    /* printf("Exchange 1 fields issend/irecv: %12.5es\n", max_exch_time1); */
    /* printf("Exchange 1 fields wait        : %12.5es\n", max_exch_time); */
    /* printf("Total Exchange 1              : %12.5es\n", max_exch_time1 + max_exch_time); */
    printf("Total Exchange 1              : %12.5es\n", max_exch_time);

    if (filedump) {
      fprintf(output, "%12.5e ", max_exch_time);
    }
  }

  /* double exch_time1 = 0; */
  /* double redondance_geom = max_exch_time1; */
  /* max_geom_time += max_exch_time1; */

  //  Exchange interpolated fields 2
  PDM_timer_resume(timer);

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, send_val2, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, recv_val2, &request);
    }
  }

  else {
    if (code_id == 1) {
     CWP_Field_issend(code_name[0], coupling_name, field_name2);
    }
    else {
     CWP_Field_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* PDM_timer_hang_on(timer); */
  /* t_end = PDM_timer_elapsed(timer); */
  /* exch_time1 = t_end - t_start; */

  /* CWP_UNUSED (max_exch_time1); */

  /* t_start = t_end; */
  /* MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); */

  /* PDM_timer_resume(timer); */

  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* PDM_timer_hang_on(timer); */
  /* t_start = PDM_timer_elapsed(timer); */
  /* PDM_timer_resume(timer); */

  /* redondance_geom += -max_exch_time1; */

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
      CWP_Field_wait_issend(code_name[0], coupling_name, field_name2);
    }
    else {
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  exch_time = t_end - t_start;
  /* CWP_UNUSED (max_exch_time1); */
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    /* printf("Exchange 2 fields issend/irecv :%12.5es\n", max_exch_time1); */
    /* printf("Exchange 2  fields wait        :%12.5es\n", max_exch_time); */
    printf("Total exchange 2               :%12.5es\n", max_exch_time);

    printf("\n\nTemps geometrie                            : %12.5es\n", max_geom_time);
    /* printf("Temps geometrie escompte (sans redondance) : %12.5es\n", */
    /*        max_geom_time - redondance_geom); */
    /* printf("Temps un Echange aux noeuds                : %12.5es\n", */
    /*        max_exch_time1 + max_exch_time); */

    /* if (filedump) { */
    /*   fprintf(output, "%12.5e\n", max_exch_time1 + max_exch_time); */
    /*   fclose(output); */
    /* } */
  }

  //  Check

  if (loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE ||
      loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) {
    double max_err = 0.;
    PDM_g_num_t n_wrong = 0;
    if (code_id == 2) {

      double *coord = malloc(sizeof(double) * pn_vtx[0] * 3);
      memcpy(coord, pvtx_coord[0], sizeof(double) * pn_vtx[0] * 3);
      if (loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE ||
          loc_method == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) {
        if (deform) {
          _unrotate(pn_vtx[0], coord);
        }
        for (int i = 0 ; i < pn_vtx[0]; i++) {
          for (int j = 0; j < 3; j++) {
            coord[3*i+j] = MIN(MAX(coord[3*i+j], xyz_min[j]), xyz_max[j]);
          }
        }
        if (deform) {
          _rotate(pn_vtx[0], coord);
        }
      }

      if (0) {
        log_trace("recv_val : %p\n", (void *) recv_val);
        // PDM_log_trace_array_double(recv_val, 3 * n_located, "recv_val     : ");
        log_trace("recv_val / coord = \n");
        if (interlaced_field || version == CWP_VERSION_OLD) {
          for (int i = 0 ; i < n_located; i++) {
            int ivtx = located[i] - 1;
            log_trace("%d ("PDM_FMT_G_NUM"): %3.3f %3.3f %3.3f / %3.3f %3.3f %3.3f\n",
                      located[i],
                      pvtx_ln_to_gn[0][located[i]-1],
                      recv_val[3*i], recv_val[3*i+1], recv_val[3*i+2],
                      coord[3*ivtx], coord[3*ivtx+1], coord[3*ivtx+2]);
          }
        } else {
          for (int j = 0; j < 3; j++) {
            for (int i = 0 ; i < n_located; i++) {
              log_trace("%3.3f ", recv_val[n_located*j + i]);
            }
            log_trace("\n");
          }

          log_trace("/\n");
          for (int j = 0; j < 3; j++) {
            for (int i = 0 ; i < n_located; i++) {
              log_trace("%3.3f ", coord[3*(located[i] - 1) + j]);
            }
            log_trace("\n");
          }
        }
      }


      for (int i = 0 ; i < n_located ; i++) {
        int is_wrong = 0;

        for (int j = 0; j < 3; j++) {
          double _recv_val;
          if (interlaced_field || version == CWP_VERSION_OLD) {
            _recv_val = recv_val[3*i + j];
          } else {
            _recv_val = recv_val[n_located*j + i];
          }

          double true_val = coord[3 * (located[i] -1) + j];

          double err = ABS (_recv_val - true_val);
          if (err > 1.e-4) {
            is_wrong = 1;
            log_trace("[%d] !! vtx "PDM_FMT_G_NUM" %d err = %g (coord#%d = %f (%f), recv = %f)\n",
                   rank, pvtx_ln_to_gn[0][(located[i] - 1)], located[i], err, j, coord[3*(located[i]-1) + j], true_val, _recv_val);
            int ivtx = located[i] - 1;
            log_trace(" coord = %f %f %f\n", pvtx_coord[0][3*ivtx], pvtx_coord[0][3*ivtx+1], pvtx_coord[0][3*ivtx+2]);

          }
          if (err > max_err) {
            max_err = err;
          }
        }

        n_wrong += is_wrong;
      }

      if (coord != pvtx_coord[0]) {
        free(coord);
      }
    }

    double global_max_err = 0.;
    PDM_g_num_t global_n_wrong = 0;
    MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    PDM_MPI_Reduce(&n_wrong, &global_n_wrong, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, 0, PDM_MPI_COMM_WORLD);
    if (rank == 0) {
      printf("Max error = %g ("PDM_FMT_G_NUM" wrong points)\n", global_max_err, global_n_wrong);
    }
  }

  //  Delete interface mesh
  if (version == CWP_VERSION_NEW) {
    CWP_Mesh_interf_del(code_name[0], coupling_name);
  }
  //  Delete coupling
  if (version == CWP_VERSION_OLD) {
    cwipi_delete_coupling(coupling_name);
  }
  else {
    CWP_Time_step_end(code_name[0]);
    CWP_Cpl_del(code_name[0], coupling_name);
  }
  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(intra_comm);

  if (current_rank_has_mesh) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free (pcell_face_idx[ipart]);
      free (pcell_face[ipart]);
      free (pcell_vtx[ipart]);
      free (pface_vtx_idx[ipart]);
      free (pface_vtx[ipart]);
      free (pvtx_coord[ipart]);
      free (pcell_ln_to_gn[ipart]);
      free (pface_ln_to_gn[ipart]);
      free (pvtx_ln_to_gn[ipart]);
    }
  }

  free(pn_vtx);
  free(pn_cell);
  free(pn_face);
  free(pvtx_coord);
  free(pvtx_ln_to_gn);
  free(pcell_face_idx);
  free(pcell_face);
  free(pcell_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(cellVtxIdx);

  if (code_id == 1) {
    free(send_val);
    free(send_val2);
  }
  else {
    free(recv_val);
    free(recv_val2);
  }
  PDM_timer_free(timer);

  //  Finalize CWIPI
  if (version == CWP_VERSION_OLD) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }

  if (version == CWP_VERSION_NEW) {
    double min_elaps_create_ptb;
    double max_elaps_create_ptb;
    double min_cpu_create_ptb;
    double max_cpu_create_ptb;
    double min_elaps_create2_ptb;
    double max_elaps_create2_ptb;
    double min_cpu_create2_ptb;
    double max_cpu_create2_ptb;
    double min_elaps_exch_ptb;
    double max_elaps_exch_ptb;
    double min_cpu_exch_ptb;
    double max_cpu_exch_ptb;

    PDM_part_to_block_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_ptb,
                                       &max_elaps_create_ptb,
                                       &min_cpu_create_ptb,
                                       &max_cpu_create_ptb,
                                       &min_elaps_create2_ptb,
                                       &max_elaps_create2_ptb,
                                       &min_cpu_create2_ptb,
                                       &max_cpu_create2_ptb,
                                       &min_elaps_exch_ptb,
                                       &max_elaps_exch_ptb,
                                       &min_cpu_exch_ptb,
                                       &max_cpu_exch_ptb);

    double min_elaps_create_btp;
    double max_elaps_create_btp;
    double min_cpu_create_btp;
    double max_cpu_create_btp;
    double min_elaps_exch_btp;
    double max_elaps_exch_btp;
    double min_cpu_exch_btp;
    double max_cpu_exch_btp;

    PDM_block_to_part_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_btp,
                                       &max_elaps_create_btp,
                                       &min_cpu_create_btp,
                                       &max_cpu_create_btp,
                                       &min_elaps_exch_btp,
                                       &max_elaps_exch_btp,
                                       &min_cpu_exch_btp,
                                       &max_cpu_exch_btp);

    if (rank == 0) {
      printf("Global time in PDM_part_to_block : \n");
      printf("   - ptb min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_ptb,
             max_elaps_create_ptb);
      printf("   - ptb min max elaps create2 : %12.5e %12.5e\n",
             min_elaps_create2_ptb,
             max_elaps_create2_ptb);
      printf("   - ptb min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_ptb,
             max_elaps_exch_ptb);
      fflush(stdout);

      printf("Global time in PDM_block_to_part : \n");
      printf("   - btp min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_btp,
             max_elaps_create_btp);
      printf("   - btp min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_btp,
             max_elaps_exch_btp);
      fflush(stdout);
    }
  }

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
