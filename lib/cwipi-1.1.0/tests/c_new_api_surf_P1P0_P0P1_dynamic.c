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

#include "grid_mesh.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"

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
  int                    argc,
  char                 **argv,
  int                   *n_vtx_seg1,
  int                   *n_vtx_seg2,
  PDM_split_dual_t      *part_method,
  double                *tolerance,
  int                   *randomize,
  int                   *variable_mesh
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
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
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
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
    else if (strcmp(argv[i], "-variable_mesh") == 0) {
      *variable_mesh = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}












static void
_gen_mesh
(
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
 const PDM_g_num_t         n_vtx_seg,
 const int                 randomize,
 const int                 random_seed,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 PDM_g_num_t            ***pface_ln_to_gn,
 PDM_g_num_t            ***pvtx_ln_to_gn
 )
{
  /* Generate a distributed polygonal mesh */
  PDM_g_num_t  ng_face         = 0;
  PDM_g_num_t  ng_vtx          = 0;
  PDM_g_num_t  ng_edge         = 0;
  int          dn_vtx          = 0;
  double      *dvtx_coord      = NULL;
  int          dn_face         = 0;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  PDM_g_num_t *dface_edge      = NULL;
  int          dn_edge         = 0;
  PDM_g_num_t *dedge_vtx       = NULL;
  PDM_g_num_t *dedge_face      = NULL;
  int          n_edge_group    = 0;
  int         *dedge_group_idx = NULL;
  PDM_g_num_t *dedge_group     = NULL;

  PDM_poly_surf_gen(comm,
                    0.,
                    0.1,
                    0.,
                    0.1,
                    randomize,
                    random_seed,
                    n_vtx_seg,
                    n_vtx_seg,
                    &ng_face,
                    &ng_vtx,
                    &ng_edge,
                    &dn_vtx,
                    &dvtx_coord,
                    &dn_face,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dface_edge,
                    &dn_edge,
                    &dedge_vtx,
                    &dedge_face,
                    &n_edge_group,
                    &dedge_group_idx,
                    &dedge_group);

  /* Spit the mesh */
  int n_zone = 1;
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);
  free(n_part_zones);

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                        0,
                                        dn_face,
                                        dn_edge,
                                        dn_vtx,
                                        comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord,
                          PDM_OWNERSHIP_USER);

  int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             dedge_vtx,
                             dedge_vtx_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                             dedge_face,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_bound_set(dmesh,
                      PDM_BOUND_TYPE_EDGE,
                      n_edge_group,
                      dedge_group,
                      dedge_group_idx,
                      PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(mpart, 0, dmesh);

  /* Run */
  PDM_multipart_compute(mpart);



  /* Get partitioned mesh */
  *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
  *pface_vtx_idx  = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int *face_edge     = NULL;
    int *face_edge_idx = NULL;
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge_idx,
                                                     &face_edge,
                                                     PDM_OWNERSHIP_USER);

    PDM_g_num_t* face_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_USER);


    int *edge_vtx     = NULL;
    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &edge_vtx,
                                        PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 0,
                                                 i_part,
                                                 &vtx_coord,
                                                 PDM_OWNERSHIP_USER);

    PDM_g_num_t* vtx_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &vtx_ln_to_gn,
                                    PDM_OWNERSHIP_USER);

    (*pn_face)[i_part] = n_face;
    (*pn_vtx)[i_part]  = n_vtx;

    /* Vertices */
    (*pvtx_coord)[i_part] = vtx_coord;
    (*pvtx_ln_to_gn)[i_part] = vtx_ln_to_gn;


    /* Faces */
    (*pface_vtx_idx)[i_part] = face_edge_idx;

    PDM_compute_face_vtx_from_face_and_edge_unsigned(n_face,
                                                     face_edge_idx,
                                                     face_edge,
                                                     edge_vtx,
                                                     *pface_vtx + i_part);
    free(face_edge);

    (*pface_ln_to_gn)[i_part] = face_ln_to_gn;
  }
  PDM_multipart_free(mpart);
  PDM_dmesh_free(dmesh);

  free(dvtx_coord);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dface_edge);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face);
  free(dedge_group_idx);
  free(dedge_group);
}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // Read args from command line
  int    n_vtx_seg1            = 4;
  int    n_vtx_seg2            = 4;
  int    randomize             = 1;
  double tolerance             = 1e-2;
  int    variable_mesh         = 0;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif


  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &part_method,
             &tolerance,
             &randomize,
             &variable_mesh);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);


  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id[2];
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
    code_id[0] = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id[0] = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }


  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);


  if (rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  CWP_Dynamic_mesh_t displacement = CWP_DYNAMIC_MESH_DEFORMABLE;
  if (variable_mesh) {
    displacement = CWP_DYNAMIC_MESH_VARIABLE;
  }
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1_dynamic";
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_create(code_name[i_code],              // Code name
                   cpl_name,                       // Coupling id
                   coupled_code_name[i_code],      // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,         // Coupling type
                   spatial_interp,
                   n_part,                         // Partition number
                   displacement,                   // Mesh displacement type
                   CWP_TIME_EXCH_USER_CONTROLLED); // Postprocessing frequency
  }


  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 2,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  if (rank == 0) {
    printf("Create coupling OK\n");
  }



  // Mesh definition
  int          **pn_face        = (int          **) malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = (int          **) malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = (int         ***) malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = (int         ***) malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = (double      ***) malloc(sizeof(double      **) * n_code);
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc(sizeof(PDM_g_num_t **) * n_code);


  for (int i_code = 0 ; i_code < n_code ; i_code++) {


    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

    _gen_mesh(mesh_comm,
              n_part,
              part_method,
              n_vtx_seg,
              randomize,
              code_id[i_code],
              &pn_face[i_code],
              &pn_vtx[i_code],
              &pface_vtx_idx[i_code],
              &pface_vtx[i_code],
              &pvtx_coord[i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn[i_code]);


    if (!variable_mesh) {
      CWP_Mesh_interf_del(code_name[i_code],
                          cpl_name);
      CWP_Mesh_interf_vtx_set(code_name[i_code],
                              cpl_name,
                              0,
                              pn_vtx[i_code][0],
                              pvtx_coord[i_code][0],
                              pvtx_ln_to_gn[i_code][0]);

      int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                               cpl_name,
                                               CWP_BLOCK_FACE_POLY);

      CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                       cpl_name,
                                       0,
                                       block_id,
                                       pn_face[i_code][0],
                                       pface_vtx_idx[i_code][0],
                                       pface_vtx[i_code][0],
                                       pface_ln_to_gn[i_code][0]);

      CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
    }
  }

  if (rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX_t0";
  const char *field_name2 = "code2_elt_gnum";
  double **send_val = (double **) malloc(sizeof(double *) * n_code);
  double **recv_val = (double **) malloc(sizeof(double *) * n_code);
  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (code_id[i_code] == 1) {
      send_val[i_code] = (double *) malloc(sizeof(double) * 3 * pn_vtx[i_code][0]);
      recv_val[i_code] = (double *) malloc(sizeof(double) * pn_face[i_code][0]);
      for (int i = 0 ; i < pn_vtx[i_code][0] ; i++) {
        send_val[i_code][3*i  ] = pvtx_coord[i_code][0][3*i  ];
        send_val[i_code][3*i+1] = pvtx_coord[i_code][0][3*i+1];
        send_val[i_code][3*i+2] = pvtx_coord[i_code][0][3*i+2];
      }
    }
    else {
      send_val[i_code] = (double *) malloc(sizeof(double) * pn_face[i_code][0]);
      recv_val[i_code] = (double *) malloc(sizeof(double) * 3 * pn_vtx[i_code][0]);
      for (int i = 0 ; i < pn_face[i_code][0] ; i++) {
        send_val[i_code][i] = (double) pface_ln_to_gn[i_code][0][i];
      }
    }
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    if (code_id[i_code] == 1) {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val[i_code]);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val[i_code]);
    }
    else {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name1,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val[i_code]);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[i_code],
                         cpl_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val[i_code]);
    }
  }





  double recv_time = 0.;
  for (int step = 0; step < 10; step++) {

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      printf("  Step %d\n", step);
    }

    // Mesh rotation and new localisation
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      if (code_id[i_code] == 1) {
        mesh_rotate(pvtx_coord[i_code][0], pn_vtx[i_code][0], recv_time);
      } else {
        mesh_rotate(pvtx_coord[i_code][0], pn_vtx[i_code][0], 3 * recv_time);
      }

      // Start time step
      CWP_Time_step_beg(code_name[i_code],
                        recv_time);

      if (variable_mesh && step%2 == 0) {
        CWP_Mesh_interf_del(code_name[i_code],
                            cpl_name);

        CWP_Mesh_interf_vtx_set(code_name[i_code],
                                cpl_name,
                                0,
                                pn_vtx[i_code][0],
                                pvtx_coord[i_code][0],
                                pvtx_ln_to_gn[i_code][0]);

        int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                                 cpl_name,
                                                 CWP_BLOCK_FACE_POLY);

        CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                         cpl_name,
                                         0,
                                         block_id,
                                         pn_face[i_code][0],
                                         pface_vtx_idx[i_code][0],
                                         pface_vtx[i_code][0],
                                         pface_ln_to_gn[i_code][0]);

        CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
      }
    }

    // Separate loops to avoid deadlock if multiple codes on same MPI rank
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      CWP_Spatial_interp_property_set(code_name[i_code],
                                      cpl_name,
                                      "n_neighbors",
                                      CWP_INT,
                                      "1");

      CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (step%3 == 0) {
      for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 1) {
          CWP_Field_issend(code_name[i_code], cpl_name, field_name1);
          CWP_Field_irecv (code_name[i_code], cpl_name, field_name2);
        }
        else {
          CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
          CWP_Field_issend(code_name[i_code], cpl_name, field_name2);
        }
      }


      for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 1) {
          CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
          CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name2);
        }
        else {
          CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
          CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name2);
        }
      }
    }

    // End time step
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
      CWP_Time_step_end(code_name[i_code]);
    }

    // Increase
    recv_time += 1.;

  }


  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Mesh_interf_del(code_name[i_code], cpl_name);
  }

  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    CWP_Cpl_del(code_name[i_code], cpl_name);
  }



  for (int i_code = 0 ; i_code < n_code ; i_code++) {
    for (int i_part = 0 ; i_part < n_part ; i_part++) {
      free(pface_vtx_idx[i_code][i_part]);
      free(pface_vtx[i_code][i_part]);
      free(pvtx_coord[i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn[i_code][i_part]);
    }
    free(pn_face[i_code]);
    free(pn_vtx[i_code]);
    free(pface_vtx_idx[i_code]);
    free(pface_vtx[i_code]);
    free(pvtx_coord[i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn[i_code]);

    free(send_val[i_code]);
    free(recv_val[i_code]);
  }
  free(pn_face);
  free(pn_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);

  free(send_val);
  free(recv_val);


  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
