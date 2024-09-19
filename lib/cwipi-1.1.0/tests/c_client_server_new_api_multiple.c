/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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
     "  -n_rank1 \n\n"
     "  -n_rank2 \n\n"
     "  -h       This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int                    argc,
 char                 **argv,
 int                    code_n_rank[]
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
  int                   code_n_rank[2]    = {-1, -1};
  _read_args(argc,
             argv,
             code_n_rank);

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

  char *config = NULL;
  if (is_code1) {
    config = (char *) "client_new_api_multiple_o/code1/cwp_config_srv.txt";
  }
  else {
    config = (char *) "client_new_api_multiple_o/code2/cwp_config_srv.txt";
  }


  // Launch server
  if (i_rank == 0) {
    char cmd[999];
    system("mkdir -p client_new_api_multiple_o/code1");
    system("mkdir -p client_new_api_multiple_o/code2");
    system("rm -f ./client_new_api_multiple_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_multiple_o/code2/cwp_config_srv.txt");

    sprintf(cmd, "mpirun -n %d cwp_server -cn code1 -p %d %d -c \"client_new_api_multiple_o/code1/cwp_config_srv.txt\" \
                  : -n %d  cwp_server -cn code2 -p %d %d -c \"client_new_api_multiple_o/code2/cwp_config_srv.txt\" &",
                  code_n_rank[0], 51100, 51100 + code_n_rank[0] - 1,
                  code_n_rank[1], 51100 + code_n_rank[0], 51100 + code_n_rank[0] + code_n_rank[1] - 1);
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
  int           n_part            = 1;

  if (is_code1) {
    code_name         = "code1";
    coupled_code_name = "code2";
  }
  else {
    code_name         = "code2";
    coupled_code_name = "code1";
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

  // Create coupling
  const char *cpl_name = "client_new_api_multiple";
  CWP_client_Cpl_create(code_name,
                        cpl_name,
                        coupled_code_name,
                        CWP_INTERFACE_SURFACE,
                        CWP_COMM_PAR_WITH_PART,
                        CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                        n_part,
                        CWP_DYNAMIC_MESH_VARIABLE,
                        CWP_TIME_EXCH_USER_CONTROLLED);

  // Begin Time step
  CWP_client_Time_step_beg(code_name, 0.0);

  // Create mesh with several blocks
  int n_vtx = 9;
  double coords[27] = {0., 0., 0.,   1., 0., 0.,   2., 0., 0.,   0., 1., 0.,
                       1., 1., 0.,   2., 1., 0.,   0., 2., 0.,   1., 2., 0.,   2., 2., 0.};
  CWP_client_Mesh_interf_vtx_set(code_name,
                                 cpl_name,
                                 0,
                                 n_vtx,
                                 coords,
                                 NULL);

  int first_block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                        cpl_name,
                                                        CWP_BLOCK_FACE_TRIA3);
  int n_first_elts = 4;
  int first_connec[12] = {1, 2, 5,   1, 5, 4,   5, 6, 9,   5, 9, 8};
  CWP_client_Mesh_interf_block_std_set(code_name,
                                       cpl_name,
                                       0,
                                       first_block_id,
                                       n_first_elts,
                                       first_connec,
                                       NULL);

  int second_block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                         cpl_name,
                                                         CWP_BLOCK_FACE_QUAD4);

  int n_second_elts = 2;
  int second_connec[8] = {2, 3, 6, 5,   4, 5, 8, 7};
  CWP_client_Mesh_interf_block_std_set(code_name,
                                       cpl_name,
                                       0,
                                       second_block_id,
                                       n_second_elts,
                                       second_connec,
                                       NULL);

  CWP_client_Mesh_interf_finalize(code_name, cpl_name);

  // Delete mesh
  CWP_client_Mesh_interf_del(code_name, cpl_name);

  // End time step
  CWP_client_Time_step_end(code_name);

  // Begin Time step
  CWP_client_Time_step_beg(code_name, 1.0);

  // Create new mesh
  int n_second_vtx = 6;
  double second_coodrs[18] = {0., 0., 0.,   1., 0., 0.,   2., 0., 0.,   0., 1., 0.,
                              1., 1., 0.,   2., 1., 0.};
  CWP_client_Mesh_interf_vtx_set(code_name,
                                 cpl_name,
                                 0,
                                 n_second_vtx,
                                 second_coodrs,
                                 NULL);

  int third_block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                        cpl_name,
                                                        CWP_BLOCK_FACE_QUAD4);

  int n_third_elts = 2;
  int third_connec[8] = {1, 2, 5, 4,  2, 3, 6, 5};
  CWP_client_Mesh_interf_block_std_set(code_name,
                                       cpl_name,
                                       0,
                                       third_block_id,
                                       n_third_elts,
                                       third_connec,
                                       NULL);

  // End time step
  CWP_client_Time_step_end(code_name);

  // Delete mesh
  CWP_client_Mesh_interf_del(code_name, cpl_name);

  // Delete coupling
  CWP_client_Cpl_del(code_name, cpl_name);

  // Finalize
  CWP_client_Finalize();

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  return exit_check;
}
