/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "grid_mesh.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_io.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_generate_mesh.h"
#include "client_server/client.h"

#include "cwp_priv.h"

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int                    argc,
 char                 **argv,
 int                    code_n_rank[],
 int                    code_n_part[],
 PDM_g_num_t            code_n_vtx_seg[],
 PDM_Mesh_nodal_elt_t   code_elt_type[],
 char                 **config  // filename for server ip adresses + ports
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);
    else if (strcmp(argv[i], "-n_rank1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_part[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_part[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_vtx_seg[0] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_vtx_seg[1] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-e1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_elt_type[0] = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-e2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_elt_type[1] = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
)
{
  // default
  char                 *config = NULL;
  int                   code_n_rank[2]    = {-1, -1};
  int                   code_n_part[2]    = {1, 1};
  PDM_g_num_t           code_n_vtx_seg[2] = {10, 5};
  PDM_Mesh_nodal_elt_t  code_elt_type[2]  = {PDM_MESH_NODAL_TRIA3, PDM_MESH_NODAL_QUAD4};
  _read_args(argc,
             argv,
             code_n_rank,
             code_n_part,
             code_n_vtx_seg,
             code_elt_type,
             &config);

  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

  if (code_n_rank[0] < 0) {
    if (code_n_rank[1] < 0) {
      code_n_rank[0] = n_rank / 2;
    }
    else {
      code_n_rank[0] = n_rank - code_n_rank[1];
    }
  }
  if (code_n_rank[1] < 0) {
    code_n_rank[1] = n_rank - code_n_rank[0];
  }

  int is_code1 = (i_rank < code_n_rank[0]);

  if (config == NULL) {
    if (is_code1) {
      config = (char *) "client_new_api_surf_P1P0_P0P1_o/code1/cwp_config_srv.txt";
    }
    else {
      config = (char *) "client_new_api_surf_P1P0_P0P1_o/code2/cwp_config_srv.txt";
    }
  }

  // assert(n_rank % 2 == 0);

  // launch server

  if (i_rank == 0) {
    char cmd[999];
    system("mkdir -p client_new_api_surf_P1P0_P0P1_o/code1");
    system("mkdir -p client_new_api_surf_P1P0_P0P1_o/code2");
    system("rm -f ./client_new_api_surf_P1P0_P0P1_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_surf_P1P0_P0P1_o/code2/cwp_config_srv.txt");

    sprintf(cmd, "mpirun -n %d cwp_server -cn code1 -p %d %d -c \"client_new_api_surf_P1P0_P0P1_o/code1/cwp_config_srv.txt\" \
                  : -n %d  cwp_server -cn code2 -p %d %d -c \"client_new_api_surf_P1P0_P0P1_o/code2/cwp_config_srv.txt\" &",
                  code_n_rank[0], 52100, 52100 + code_n_rank[0] - 1,
                  code_n_rank[1], 52100 + code_n_rank[0], 52100 + code_n_rank[0] + code_n_rank[1] - 1);
    system(cmd);
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(10);

  MPI_Barrier(comm);

  // Initialization
  const char   *code_name         = NULL;
  const char   *coupled_code_name = NULL;
  CWP_Status_t  is_active_rank    = CWP_STATUS_ON;

  int           n_part;
  PDM_g_num_t   n_vtx_seg;
  PDM_Mesh_nodal_elt_t elt_type;
  if (is_code1) {
    code_name         = "code1";
    coupled_code_name = "code2";
    n_part            = code_n_part[0];
    n_vtx_seg         = code_n_vtx_seg[0];
    elt_type          = code_elt_type[0];
  }
  else {
    code_name         = "code2";
    coupled_code_name = "code1";
    n_part            = code_n_part[1];
    n_vtx_seg         = code_n_vtx_seg[1];
    elt_type          = code_elt_type[1];
  }
  assert(n_part > 0);

  MPI_Comm intra_comm;
  MPI_Comm_split(comm, is_code1, i_rank, &intra_comm);

  int intra_rank;
  int n_intra_rank;
  MPI_Comm_rank(comm, &intra_rank);
  MPI_Comm_size(comm, &n_intra_rank);

  CWP_client_Init(intra_comm,
                  config,
                  code_name,
                  is_active_rank);

  // EXIT_SUCCESS ?
  int exit_check = 0;

  const char *cpl_name = "client_new_api_surf_P1P0_P0P1";
  CWP_client_Cpl_create(code_name,                                             // Code name
                        cpl_name,                                              // Coupling id
                        coupled_code_name,                                     // Coupled application id
                        CWP_INTERFACE_SURFACE,
                        CWP_COMM_PAR_WITH_PART,                                // Coupling type
                        CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                        n_part,                                                // Number of partitions
                        CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                        CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

  CWP_client_Visu_set(code_name,            // Code name
                      cpl_name,                // Coupling id
                      1,                       // Postprocessing frequency
                      CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                      "text");                 // Postprocessing option

  // Barrier on coupling communicator
  CWP_client_Cpl_barrier(code_name,
                         cpl_name);

  // Mesh definition
  int          *pn_vtx         = NULL;
  int          *pn_edge        = NULL;
  int          *pn_face        = NULL;
  double      **pvtx_coord     = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  CWP_g_num_t **pvtx_ln_to_gn  = NULL;
  CWP_g_num_t **pedge_ln_to_gn = NULL;
  CWP_g_num_t **pface_ln_to_gn = NULL;
  int         **pface_vtx      = NULL;

  PDM_generate_mesh_rectangle_ngon(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm),
                                   elt_type,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   n_vtx_seg,
                                   n_vtx_seg,
                                   n_part,
                                   PDM_SPLIT_DUAL_WITH_HILBERT,
                                   0,
                                   &pn_vtx,
                                   &pn_edge,
                                   &pn_face,
                                   &pvtx_coord,
                                   &pedge_vtx,
                                   &pface_edge_idx,
                                   &pface_edge,
                                   &pface_vtx,
                                   &pvtx_ln_to_gn,
                                   &pedge_ln_to_gn,
                                   &pface_ln_to_gn);

  int block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                  cpl_name,
                                                  CWP_BLOCK_FACE_POLY);

  for (int ipart = 0; ipart < n_part; ipart++) {

    CWP_client_Mesh_interf_vtx_set(code_name,
                                   cpl_name,
                                   ipart,
                                   pn_vtx       [ipart],
                                   pvtx_coord   [ipart],
                                   pvtx_ln_to_gn[ipart]);

    CWP_client_Mesh_interf_f_poly_block_set(code_name,
                                            cpl_name,
                                            ipart,
                                            block_id,
                                            pn_face       [ipart],
                                            pface_edge_idx[ipart],
                                            pface_vtx     [ipart],
                                            pface_ln_to_gn[ipart]);

    // --> check
    int GETnElts = -1;
    int *GETeltsConnecPointer  = NULL;
    int *GETeltsConnec         = NULL;
    CWP_g_num_t *GETglobal_num = NULL;
    CWP_client_Mesh_interf_f_poly_block_get(code_name,
                                            cpl_name,
                                            ipart,
                                            block_id,
                                            &GETnElts,
                                            &GETeltsConnecPointer,
                                            &GETeltsConnec,
                                            &GETglobal_num);

    if (GETnElts != pn_face[ipart]) {
      exit_check = 1;
    }

    for (int i = 0; i <= pn_face[ipart]; i++) {
      if (GETeltsConnecPointer[i] != pface_edge_idx[ipart][i]) {
        exit_check = 1;
        break;
      }
    }

    for (int i = 0; i <= pn_face[ipart]; i++) {
      if (GETeltsConnec[i] != pface_vtx[ipart][i]) {
        exit_check = 1;
        break;
      }
    }

    for (int i = 0; i < pn_face[ipart]; i++) {
      if (GETglobal_num[i] != pface_ln_to_gn[ipart][i]) {
        exit_check = 1;
        break;
      }
    }

    if (exit_check) {
      break;
    }
  }


  CWP_client_Mesh_interf_finalize(code_name, cpl_name);

  double **send_values = malloc(sizeof(double *) * n_part);
  double **recv_values = malloc(sizeof(double *) * n_part);

  if (is_code1) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      send_values[ipart] = malloc(sizeof(double) * pn_vtx [ipart]);
      recv_values[ipart] = malloc(sizeof(double) * pn_face[ipart]);

      for (int i = 0 ; i < pn_vtx[ipart] ; i++) {
        send_values[ipart][i] = pvtx_coord[ipart][3 * i];
      }
    }
  }
  else {
    for (int ipart = 0; ipart < n_part; ipart++) {
      send_values[ipart] = malloc(sizeof(double) * pn_face[ipart]);
      recv_values[ipart] = malloc(sizeof(double) * pn_vtx [ipart]);

      for (int i = 0 ; i < pn_face[ipart] ; i++) {
        send_values[ipart][i] = i_rank;
      }
    }
  }

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX";
  const char *field_name2 = "i_rank";

  CWP_Status_t visu_status = CWP_STATUS_ON;
  if (is_code1) {
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name1,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_NODE,
                            CWP_FIELD_EXCH_SEND,
                            visu_status);
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name2,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_CELL_CENTER,
                            CWP_FIELD_EXCH_RECV,
                            visu_status);

    for (int ipart = 0; ipart < n_part; ipart++) {
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name1,
                                0,
                                CWP_FIELD_MAP_SOURCE,
                                send_values[ipart]);
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name2,
                                0,
                                CWP_FIELD_MAP_TARGET,
                                recv_values[ipart]);
    }
  }
  else {
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name1,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_NODE,
                            CWP_FIELD_EXCH_RECV,
                            visu_status);
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name2,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_CELL_CENTER,
                            CWP_FIELD_EXCH_SEND,
                            visu_status);

    for (int ipart = 0; ipart < n_part; ipart++) {
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name2,
                                ipart,
                                CWP_FIELD_MAP_SOURCE,
                                send_values[ipart]);
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name1,
                                ipart,
                                CWP_FIELD_MAP_TARGET,
                                recv_values[ipart]);
    }
  }

  CWP_client_Spatial_interp_weights_compute(code_name, cpl_name);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_uncomputed = 0;

    if (is_code1) {
      n_uncomputed = CWP_client_N_uncomputed_tgts_get(code_name, cpl_name, field_name2, ipart);
    }
    else {
      n_uncomputed = CWP_client_N_uncomputed_tgts_get(code_name, cpl_name, field_name1, ipart);
    }

    // --> check
    if (!(n_uncomputed == 0)) {
      exit_check = 1;
    }
  }

  if (is_code1) {
    CWP_client_Field_issend(code_name, cpl_name, field_name1);
    CWP_client_Field_irecv (code_name, cpl_name, field_name2);
  }
  else {
    CWP_client_Field_irecv (code_name, cpl_name, field_name1);
    CWP_client_Field_issend(code_name, cpl_name, field_name2);
  }

  if (is_code1) {
    CWP_client_Field_wait_issend(code_name, cpl_name, field_name1);
    CWP_client_Field_wait_irecv (code_name, cpl_name, field_name2);
  }
  else {
    CWP_client_Field_wait_irecv (code_name, cpl_name, field_name1);
    CWP_client_Field_wait_issend(code_name, cpl_name, field_name2);
  }

  // Mesh_interf_from_faceedge_set
  for (int ipart = 0; ipart < n_part; ipart++) {
    CWP_client_Mesh_interf_from_faceedge_set(code_name,
                                             cpl_name,
                                             ipart,
                                             pn_face       [ipart],
                                             pface_edge_idx[ipart],
                                             pface_edge    [ipart],
                                             pn_edge       [ipart],
                                             pedge_vtx     [ipart],
                                             pface_ln_to_gn[ipart]);
  }


  CWP_client_Mesh_interf_del(code_name, cpl_name);

  CWP_client_Cpl_del(code_name, cpl_name);

  // Freeing memory
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pvtx_coord    [ipart]);
    free(pvtx_ln_to_gn [ipart]);
    free(pedge_vtx     [ipart]);
    free(pedge_ln_to_gn[ipart]);
    free(pface_edge_idx[ipart]);
    free(pface_edge    [ipart]);
    free(pface_vtx     [ipart]);
    free(pface_ln_to_gn[ipart]);
    free(send_values   [ipart]);
    free(recv_values   [ipart]);
  }
  free(pn_vtx        );
  free(pn_edge       );
  free(pn_face       );
  free(pvtx_coord    );
  free(pvtx_ln_to_gn );
  free(pedge_vtx     );
  free(pedge_ln_to_gn);
  free(pface_edge_idx);
  free(pface_edge    );
  free(pface_vtx     );
  free(pface_ln_to_gn);
  free(send_values   );
  free(recv_values   );

  // Finalize
  CWP_client_Finalize();

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  return exit_check;
}
