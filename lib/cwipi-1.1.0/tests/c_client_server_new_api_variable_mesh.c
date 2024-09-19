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
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
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
    config = (char *) "client_new_api_variable_mesh_o/code1/cwp_config_srv.txt";
  }
  else {
    config = (char *) "client_new_api_variable_mesh_o/code2/cwp_config_srv.txt";
  }


  // Launch server
  if (i_rank == 0) {
    char cmd[999];
    system("mkdir -p client_new_api_variable_mesh_o/code1");
    system("mkdir -p client_new_api_variable_mesh_o/code2");
    system("rm -f ./client_new_api_variable_mesh_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_variable_mesh_o/code2/cwp_config_srv.txt");

    sprintf(cmd, "mpirun -n %d cwp_server -cn code1 -p %d %d -c \"client_new_api_variable_mesh_o/code1/cwp_config_srv.txt\" \
                  : -n %d  cwp_server -cn code2 -p %d %d -c \"client_new_api_variable_mesh_o/code2/cwp_config_srv.txt\" &",
                  code_n_rank[0], 53100, 53100 + code_n_rank[0] - 1,
                  code_n_rank[1], 53100 + code_n_rank[0], 53100 + code_n_rank[0] + code_n_rank[1] - 1);
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

  // Create coupling
  const char *cpl_name = "client_new_api_variable_mesh";
  CWP_client_Cpl_create(code_name,
                        cpl_name,
                        coupled_code_name,
                        CWP_INTERFACE_SURFACE,
                        CWP_COMM_PAR_WITH_PART,
                        CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                        n_part,
                        CWP_DYNAMIC_MESH_VARIABLE,
                        CWP_TIME_EXCH_USER_CONTROLLED);

  // Set coupling visualisation:
  CWP_client_Visu_set(code_name,
                      cpl_name,
                      1,
                      CWP_VISU_FORMAT_ENSIGHT,
                      "text");



  const char *field_name = "field";
  double *send_field = NULL;
  double *recv_field = NULL;
  if (is_code1) {
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_NODE,
                            CWP_FIELD_EXCH_SEND,
                            CWP_STATUS_ON);

  } else {
    CWP_client_Field_create(code_name,
                            cpl_name,
                            field_name,
                            CWP_DOUBLE,
                            CWP_FIELD_STORAGE_INTERLEAVED,
                            1,
                            CWP_DOF_LOCATION_NODE,
                            CWP_FIELD_EXCH_RECV,
                            CWP_STATUS_ON);
  }

  // Iteration
  const int    itdeb = 1;
  const int    itend = 5;
  const double freq  = 0.20;
  const double ampl  = 0.1;
  const double phi   = 0.1;
  double       time  = 0.0;
  double       dt    = 0.1;

  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;

  double omega = 2.0*acos(-1.0)*freq;

  for (int it = itdeb; it <= itend; it ++) {

    time = (it-itdeb)*dt;

    CWP_client_Time_step_beg(code_name, time);


    // Create mesh
    if (coords != NULL) {
      free(coords);
    }
    if (elt_vtx_idx != NULL) {
      free(elt_vtx_idx);
    }
    if (elt_vtx != NULL) {
      free(elt_vtx);
    }
    PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm),
                                           it+2,
                                           &n_vtx,
                                           &n_elt,
                                           &coords,
                                           &elt_vtx_idx,
                                           &elt_vtx);

    send_field = realloc(send_field, sizeof(double) * n_vtx);
    recv_field = realloc(recv_field, sizeof(double) * n_vtx);

    for (int i = 0; i < n_vtx; i++) {
      coords[3 * i + 2]  = ampl * (coords[3 * i]*coords[3 * i]+coords[1 + 3 * i]*coords[1 + 3 * i])*cos(omega*time+phi);
      send_field[i] = coords[3 * i + 2];
    }

    // Create new mesh
    CWP_client_Mesh_interf_del(code_name, cpl_name);

    CWP_client_Mesh_interf_vtx_set(code_name,
                                   cpl_name,
                                   0,
                                   n_vtx,
                                   coords,
                                   NULL);



    if (is_code1) {
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name,
                                0,
                                CWP_FIELD_MAP_SOURCE,
                                send_field);
    } else {
      CWP_client_Field_data_set(code_name,
                                cpl_name,
                                field_name,
                                0,
                                CWP_FIELD_MAP_TARGET,
                                recv_field);
    }



    int idx = n_elt / 2;

    int block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                    cpl_name,
                                                    CWP_BLOCK_FACE_POLY);
    CWP_client_Mesh_interf_f_poly_block_set(code_name,
                                            cpl_name,
                                            0,
                                            block_id,
                                            idx,
                                            elt_vtx_idx,
                                            elt_vtx,
                                            NULL);


    block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                cpl_name,
                                                CWP_BLOCK_FACE_TRIA3);
    CWP_client_Mesh_interf_block_std_set(code_name,
                                         cpl_name,
                                         0,
                                         block_id,
                                         n_elt - idx,
                                         elt_vtx + elt_vtx_idx[idx],
                                         NULL);

    CWP_client_Mesh_interf_finalize(code_name, cpl_name);




    MPI_Barrier(comm);

    CWP_client_Spatial_interp_weights_compute(code_name, cpl_name);

    // Exchange field
    if (is_code1) {
      CWP_client_Field_issend(code_name, cpl_name, field_name);
    }
    else {
      CWP_client_Field_irecv (code_name, cpl_name, field_name);
    }

    if (is_code1) {
      CWP_client_Field_wait_issend(code_name, cpl_name, field_name);
    }
    else {
      CWP_client_Field_wait_irecv (code_name, cpl_name, field_name);
    }

    CWP_client_Time_step_end(code_name);

    MPI_Barrier(comm);

  }

  // Delete mesh
  CWP_client_Mesh_interf_del(code_name, cpl_name);

  // Delete coupling
  CWP_client_Cpl_del(code_name, cpl_name);

  // free
  free(coords     );
  free(elt_vtx_idx);
  free(elt_vtx    );
  free(send_field );
  free(recv_field );

  // Finalize
  CWP_client_Finalize();

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
