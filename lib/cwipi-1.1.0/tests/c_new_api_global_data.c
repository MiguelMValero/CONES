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
#include <assert.h>
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

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
  int                  *verbose,
  int                  *is_mixed,
  int                  *is_joint
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-mixed") == 0) {
      *is_mixed = 1;
    }
    else if (strcmp(argv[i], "-joint") == 0) {
      *is_joint = 1;
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
main(int argc, char *argv[]) {

  int is_mixed = 1;
  int is_joint = 0;
  int verbose  = 0;

  _read_args(argc,
             argv,
             &verbose,
             &is_mixed,
             &is_joint);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  // Choice

  if (is_mixed) {

    // Initialize CWIPI
    int n_part = 1;
    int n_code = 1;
    int code_id;
    const char **code_name = malloc(sizeof(char *) * n_code);
    const char **coupled_code_name = malloc(sizeof(char *) * n_code);
    CWP_Status_t is_active_rank = CWP_STATUS_ON;

    if (rank < comm_world_size / 2) {
      code_id = 1;
      code_name[0] = "code1";
      coupled_code_name[0] = "code2";
    }
    else {
      code_id = 2;
      code_name[0] = "code2";
      coupled_code_name[0] = "code1";
    }

    MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);

    if (verbose && rank == 0) {
      printf("CWIPI Init OK\n");
    }

    // Create coupling
    const char *coupling_name = "c_surf_cpl_P1P1";
    CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
    CWP_Cpl_create(code_name[0],
                   coupling_name,
                   coupled_code_name[0],
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    if (verbose) {
      printf("Create coupling OK between %s and %s \n", code_name[0], coupled_code_name[0]);
    }

    // Get communicator
    // OLD CWIPI :
    // MPI_Comm   cpl_comm;
    // int       *cpl_ranks = NULL;
    // int size = cwipi_coupling_comm_get(coupling_name,
    //                                    &cpl_comm,
    //                                    &cpl_ranks);

    // if (verbose) {
    //   printf("cpl_comm = %d\n", cpl_comm);
    //   printf("cpl_ranks : \n");
    //   printf("size = %d\n", size);
    //   for (int i =0; i < size; i++) {
    //     printf("cpl_ranks[%d] = %d\n", i, cpl_ranks[i]);
    //   }

    //   int recv = -1;
    //   MPI_Allreduce(&rank, &recv, 1, MPI_INT, MPI_SUM, cpl_comm);
    //   printf("recv : %d\n", recv);
    //   fflush(stdout);
    // }

    // NEW CWIPI :
    // MPI_Comm   cpl_comm;
    // int       *cpl_ranks = NULL;
    // int size = CWP_Cpl_comm_get(code_name[0],
    //                             coupling_name,
    //                             &cpl_comm,
    //                             &cpl_ranks);

    // if (verbose) {
    //   printf("cpl_comm = %d\n", cpl_comm);
    //   printf("cpl_ranks : \n");
    //   printf("size = %d\n", size);
    //   for (int i =0; i < size; i++) {
    //     printf("cpl_ranks[%d] = %d\n", i, cpl_ranks[i]);
    //   }

    //   int recv = -1;
    //   MPI_Allreduce(&rank, &recv, 1, MPI_INT, MPI_SUM, cpl_comm);
    //   printf("recv : %d\n", recv);
    //   fflush(stdout);
    // }

    // if (verbose && rank == 0) {
    //   printf("Get communicator OK\n");
    //   fflush(stdout);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    // Exchange vectors
    const char *global_data_name1  = "lapin";
    size_t  s_entity1 = sizeof(double);
    int     stride1   = 2;
    int     n_entity1 = 2;

    size_t  s_entity2 = sizeof(double);
    int     stride2   = 2;
    int     n_entity2 = 2;
    // TO DO : Exchange these 3 values via GlobalDatas...

    double *send_data1 = malloc(s_entity1 * stride1 * n_entity1);
    send_data1[0] = 42.42;
    send_data1[1] = 13.10;
    send_data1[2] = 1959.07;
    send_data1[3] = 1954.02;


    // size_t  s_recv_entity1 = 0;
    // int     recv_stride1   = -1;
    // int     n_recv_entity1 = -1;
    double *recv_data1 = malloc(s_entity1 * stride1 * n_entity1);//NULL;

    const char *global_data_name2  = "capybara";
    double *send_data2 = malloc(s_entity2 * stride2 * n_entity2);
    send_data2[0] = 0.1;
    send_data2[1] = 0.2;
    send_data2[2] = 0.3;
    send_data2[3] = 0.4;

    // size_t  s_recv_entity2 = 0;
    // int     recv_stride2   = -1;
    // int     n_recv_entity2 = -1;
    double *recv_data2 = malloc(s_entity2 * stride2 * n_entity2);//NULL;

    if (code_id == 1) {
      CWP_Global_data_issend(code_name[0],
                             coupling_name,
                             global_data_name2,
                             s_entity2,
                             stride2,
                             n_entity2,
                    (void *) send_data2);
    }

    if (code_id == 2) {
      CWP_Global_data_irecv(code_name[0],
                            coupling_name,
                            global_data_name1,
                            s_entity1,
                            stride1,
                            n_entity1,
                   (void *) recv_data1);
    }

    if (code_id == 1) {
      CWP_Global_data_issend(code_name[0],
                             coupling_name,
                             global_data_name1,
                             s_entity1,
                             stride1,
                             n_entity1,
                    (void *) send_data1);
    }

    if (code_id == 2) {
      CWP_Global_data_irecv(code_name[0],
                            coupling_name,
                            global_data_name2,
                            s_entity2,
                            stride2,
                            n_entity2,
                   (void *) recv_data2);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (code_id == 2) {
      CWP_Global_data_wait_irecv(code_name[0],
                                 coupling_name,
                                 global_data_name2);
    }

    if (code_id == 1) {
      CWP_Global_data_wait_issend(code_name[0],
                                  coupling_name,
                                  global_data_name1);
    }

    if (code_id == 1) {
      CWP_Global_data_wait_issend(code_name[0],
                                  coupling_name,
                                  global_data_name2);
    }

    if (code_id == 2) {
      CWP_Global_data_wait_irecv(code_name[0],
                                 coupling_name,
                                 global_data_name1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose) {
      for (int i = 0; i < 4; i++) {
        if (code_id == 1) {
          printf("rank %d -- send[%d] : %f\n", rank, i, send_data1[i]);
          fflush(stdout);
          printf("rank %d -- send[%d] : %f\n", rank, i, send_data2[i]);
          fflush(stdout);
        } else {
          printf("rank %d -- recv[%d] : %f\n", rank, i, recv_data1[i]);
          fflush(stdout);
          printf("rank %d -- recv[%d] : %f\n", rank, i, recv_data2[i]);
          fflush(stdout);
        }
      }
    }

    // free
    free(send_data2);
    free(recv_data2);

    // Delete coupling
    CWP_Cpl_del(code_name[0], coupling_name);

    // Free memory
    free(code_name);
    free(coupled_code_name);
    free(intra_comm);

    // Finalize cwipi
    CWP_Finalize();

    // free
    free(send_data1);
    free(recv_data1);

  } // end is mixed

  if (is_joint) {

    // Initialize CWIPI
    int n_part = 1;
    int n_code = 2;
    const char **code_name = malloc(sizeof(char *) * n_code);
    const char **coupled_code_name = malloc(sizeof(char *) * n_code);
    CWP_Status_t is_active_rank = CWP_STATUS_ON;

    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    code_name[1] = "code2";
    coupled_code_name[1] = "code1";

    MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);

    if (verbose && rank == 0) {
      printf("CWIPI Init OK\n");
    }

    // Create coupling
    const char *coupling_name = "c_surf_cpl_P1P1";

    CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
    CWP_Cpl_create(code_name[0],
                   coupling_name,
                   coupled_code_name[0],
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    CWP_Cpl_create(code_name[1],
                   coupling_name,
                   coupled_code_name[1],
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    // Exchange vector
    const char *global_data_name  = "lapin";
    size_t  s_entity = sizeof(double);
    int     stride   = 2;
    int     n_entity = 2;

    double *send_data = malloc(s_entity * stride * n_entity);
    send_data[0] = 42.42;
    send_data[1] = 13.10;
    send_data[2] = 1959.07;
    send_data[3] = 1954.02;

    // size_t  s_recv_entity = 0;
    // int     recv_stride   = -1;
    // int     n_recv_entity = -1;
    double *recv_data = malloc(s_entity * stride * n_entity);//NULL;

    CWP_Global_data_issend(code_name[0],
                           coupling_name,
                           global_data_name,
                           s_entity,
                           stride,
                           n_entity,
                  (void *) send_data);

    CWP_Global_data_irecv(code_name[1],
                          coupling_name,
                          global_data_name,
                          s_entity,
                          stride,
                          n_entity,
                 (void *) recv_data);

    MPI_Barrier(MPI_COMM_WORLD);

    CWP_Global_data_wait_issend(code_name[0],
                                coupling_name,
                                global_data_name);

    CWP_Global_data_wait_irecv(code_name[1],
                               coupling_name,
                               global_data_name);

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose) {
      for (int i = 0; i < 4; i++) {
        printf("rank %d -- send[%d] : %f\n", rank, i, send_data[i]);
        fflush(stdout);
        printf("rank %d -- recv[%d] : %f\n", rank, i, recv_data[i]);
        fflush(stdout);
      }
    }

    // free
    free(send_data);
    free(recv_data);

    // Delete coupling
    CWP_Cpl_del(code_name[0], coupling_name);
    CWP_Cpl_del(code_name[1], coupling_name);

    // Free memory
    free(code_name);
    free(coupled_code_name);
    free(intra_comm);

    // Finalize cwipi
    CWP_Finalize();

  } // end joint

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

