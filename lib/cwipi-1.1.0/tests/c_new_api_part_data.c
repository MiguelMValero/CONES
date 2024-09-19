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
         "  -v              verbose\n\n"
         "  -o              Test choice\n\n"
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
  int                  *option
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

    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *option = atoi(argv[i]);
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
main(int argc, char *argv[]) {

  int joint    = 0;
  int verbose  = 0;

  _read_args (argc,
              argv,
             &verbose,
             &joint);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  if (joint) {

    assert (comm_world_size == 2);

    // Initialize CWIPI
    int n_part = 1;
    int n_code;
    if (rank == 0) {
      n_code = 1;
    }

    if (rank == 1) {
      n_code = 2;
    }
    const char **code_name = malloc(sizeof(char *) * n_code);
    const char **coupled_code_name = malloc(sizeof(char *) * n_code);
    CWP_Status_t is_active_rank = CWP_STATUS_ON;

    if (rank == 0) {
      code_name[0] = "code1";
      coupled_code_name[0] = "code2";
    }

    if (rank == 1) {
      code_name[0] = "code1";
      coupled_code_name[0] = "code2";
      code_name[1] = "code2";
      coupled_code_name[1] = "code1";
    }

    MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);

    // Create coupling
    const char *coupling_name = "couplage";
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

    if (rank == 1) {
      CWP_Cpl_create(code_name[1],
                     coupling_name,
                     coupled_code_name[1],
                     CWP_INTERFACE_VOLUME,
                     CWP_COMM_PAR_WITH_PART,
                     loc_method,
                     n_part,
                     CWP_DYNAMIC_MESH_STATIC,
                     CWP_TIME_EXCH_USER_CONTROLLED);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // create
    int n_elt_send = 6;
    int n_elt_recv = 6;

    int *send_n_elts = malloc(sizeof(int) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      send_n_elts[i_part] = n_elt_send;
    }


    int *recv_n_elts = NULL;
    if (rank == 1) {
      recv_n_elts = malloc(sizeof(int) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        recv_n_elts[i_part] = n_elt_recv;
      }
    }

    // gnum
    CWP_g_num_t **gnum_elt_send = malloc(sizeof(CWP_g_num_t *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      gnum_elt_send[i_part] = malloc(sizeof(CWP_g_num_t) * send_n_elts[i_part]);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < send_n_elts[i_part]; i++) {
        gnum_elt_send[i_part][i] = i + 1; // send_n_elts[i_part] * rank +
      }
    }
    if (verbose) {
      PDM_log_trace_array_long(gnum_elt_send[0], send_n_elts[0], "gnum_elt_send: ");
    }

    CWP_g_num_t **gnum_elt_recv = NULL;
    if (rank == 1) {
      gnum_elt_recv = malloc(sizeof(CWP_g_num_t *) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        gnum_elt_recv[i_part] = malloc(sizeof(CWP_g_num_t) * recv_n_elts[i_part]);
      }

      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < recv_n_elts[i_part]; i++) {
          gnum_elt_recv[i_part][i] = i + 1;
        }
      }
      if (verbose) {
        PDM_log_trace_array_long(gnum_elt_recv[0], recv_n_elts[0], "gnum_elt_recv: ");
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Partitionned data exchange
    const char *part_data_name = "schtroumpf";

    // pd create
    if (rank == 0 || rank == 1) {
      CWP_Part_data_create(code_name[0],
                           coupling_name,
                           part_data_name,
                           CWP_PARTDATA_SEND,
                           gnum_elt_send,
                           send_n_elts,
                           n_part);
    }

    if (rank == 1) {
      CWP_Part_data_create(code_name[1],
                           coupling_name,
                           part_data_name,
                           CWP_PARTDATA_RECV,
                           gnum_elt_recv,
                           recv_n_elts,
                           n_part);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // --> exchange
    double **send_data = NULL;
    double **recv_data = NULL;
    int n_comp = 3;

    if (rank == 0 || rank == 1) {

      send_data = malloc(sizeof(double *) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        send_data[i_part] = malloc(sizeof(double) * send_n_elts[i_part] * n_comp);
      }

      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < send_n_elts[i_part]; i++) {
          for (int i_comp = 0; i_comp < n_comp; i_comp++) {
            send_data[i_part][3*i + i_comp] = rank + 0.1*i + i_comp;
          }
        }
      }

      CWP_Part_data_issend(code_name[0],
                           coupling_name,
                           part_data_name,
                           0,
                           sizeof(double),
                           n_comp,
                 (void **) send_data);
    }

    if (rank == 1) {

      recv_data = malloc(sizeof(double *) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        recv_data[i_part] = malloc(sizeof(double) * recv_n_elts[i_part] * n_comp);
      }

      CWP_Part_data_irecv(code_name[1],
                          coupling_name,
                          part_data_name,
                          0,
                          sizeof(double),
                          n_comp,
                (void **) recv_data);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    // --> wait
    if (rank == 0 || rank == 1) {

      CWP_Part_data_wait_issend(code_name[0],
                                coupling_name,
                                part_data_name,
                                0);
    }

    if (rank == 1) {

      CWP_Part_data_wait_irecv(code_name[1],
                               coupling_name,
                               part_data_name,
                               0);
    }

    // --> check
    if (verbose) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        if (rank == 0 || rank == 1) {
          for (int i = 0; i < send_n_elts[i_part]; i++) {
            for (int i_comp = 0; i_comp < n_comp; i_comp++) {
              log_trace("%d - %ld -> s[%d][%d][%d] : %f\n", rank, gnum_elt_send[i_part][i], i_part, i, i_comp, send_data[i_part][3*i + i_comp]);
            }
          }
        }

          if (rank == 1) {
          for (int i = 0; i < n_elt_send; i++) { // recv_n_elts[i_part]
            for (int i_comp = 0; i_comp < n_comp; i_comp++) {
              log_trace("%d - %ld -> r[%d][%d][%d] : %f\n", rank, gnum_elt_recv[i_part][i], i_part, i, i_comp, recv_data[i_part][3*i + i_comp]);
            }
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // second send

    if (rank == 0 || rank == 1) {

      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < send_n_elts[i_part]; i++) {
          for (int i_comp = 0; i_comp < n_comp; i_comp++) {
            send_data[i_part][3*i + i_comp] = rank * 100 + 3*i + i_comp;
          }
        }
      }

      CWP_Part_data_issend(code_name[0],
                           coupling_name,
                           part_data_name,
                           1,
                           sizeof(double),
                           n_comp,
                 (void **) send_data);
    }

    if (rank == 1) {

      // if (recv_data != NULL) {
      //   for (int i_part = 0; i_part < n_part; i_part++) {
      //     if (recv_data[i_part] != NULL) free(recv_data[i_part]);
      //   }
      //   free(recv_data);
      // }
      // recv_data = NULL;

      CWP_Part_data_irecv(code_name[1],
                          coupling_name,
                          part_data_name,
                          1,
                          sizeof(double),
                          n_comp,
                (void **) recv_data);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    // --> wait
    if (rank == 0 || rank == 1) {

      CWP_Part_data_wait_issend(code_name[0],
                                coupling_name,
                                part_data_name,
                                1);
    }

    if (rank == 1) {

      CWP_Part_data_wait_irecv(code_name[1],
                               coupling_name,
                               part_data_name,
                               1);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    // --> check
    if (verbose) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        if (rank == 0 || rank == 1) {
          for (int i = 0; i < send_n_elts[i_part]; i++) {
            for (int i_comp = 0; i_comp < n_comp; i_comp++) {
              log_trace("%d - %ld -> s2[%d][%d][%d] : %f\n", rank, gnum_elt_send[i_part][i], i_part, i, i_comp, send_data[i_part][3*i + i_comp]);
            }
          }
        }

        if (rank == 1) {
          for (int i = 0; i < n_elt_send; i++) { // recv_n_elts[i_part]
            for (int i_comp = 0; i_comp < n_comp; i_comp++) {
              log_trace("%d - %ld -> r2[%d][%d][%d] : %f\n", rank, gnum_elt_recv[i_part][i], i_part, i, i_comp, recv_data[i_part][3*i + i_comp]);
            }
          }
        }
      }
    }

    // del
    if (rank == 0) {
      CWP_Part_data_del(code_name[0],
                        coupling_name,
                        part_data_name);
    }
    CWP_Part_data_del(code_name[1],
                      coupling_name,
                      part_data_name);


    MPI_Barrier(MPI_COMM_WORLD);

    // Delete coupling
    if (rank == 0 || rank == 1) {
      CWP_Cpl_del(code_name[0], coupling_name);
    }
    if (rank == 1) {
      CWP_Cpl_del(code_name[1], coupling_name);
    }

    // Free memory
    free(code_name);
    free(coupled_code_name);
    free(intra_comm);

    if (rank == 0 || rank == 1) {
      if (send_data != NULL) {
        for (int i_part = 0; i_part < n_part; i_part++) {
          free(send_data[i_part]);
        }
        free(send_data);
      }
      if (gnum_elt_send != NULL) {
        for (int i_part = 0; i_part < n_part; i_part++) {
          if (gnum_elt_send[i_part] != NULL) free(gnum_elt_send[i_part]);
        }
        free(gnum_elt_send);
      }
      if (send_n_elts != NULL) free(send_n_elts);
    }
    if (rank == 1) {
      if (recv_data != NULL) {
        for (int i_part = 0; i_part < n_part; i_part++) {
          if (recv_data[i_part] != NULL) free(recv_data[i_part]);
        }
        free(recv_data);
      }
      if (gnum_elt_recv != NULL) {
        for (int i_part = 0; i_part < n_part; i_part++) {
          if (gnum_elt_recv[i_part] != NULL) free(gnum_elt_recv[i_part]);
        }
        free(gnum_elt_recv);
      }
      if (recv_n_elts != NULL) free(recv_n_elts);
    }

    // Finalize cwipi
    CWP_Finalize();

  } else {

    assert (comm_world_size == 4);

    // Initialize CWIPI
    int code_id;
    if (rank == 0) {
      code_id = 1;
    } else {
      code_id = 2;
    }
    int n_part;
     if (code_id == 1) {
       n_part = 1;
    } else {
       n_part = 2;
    }
    int n_code = 1;
    const char **code_name = malloc(sizeof(char *) * n_code);
    const char **coupled_code_name = malloc(sizeof(char *) * n_code);
    CWP_Status_t is_active_rank = CWP_STATUS_ON;

    if (code_id == 1) {
      code_name[0] = "code1";
      coupled_code_name[0] = "code2";
    }
    else {
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

    // Partitionned data exchange
    const char *part_data_name = "schtroumpf";

    // --> create
    CWP_PartData_exch_t side;
    if (code_id == 1) {
      side = CWP_PARTDATA_SEND;
    }
    else {
      side = CWP_PARTDATA_RECV;
    }

    int n_elt;
    if (code_id == 1) {
       n_elt = 12;
    } else {
       n_elt = 2;
    }
    int *n_elts = malloc(sizeof(int *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      n_elts[i_part] = n_elt;
    }

    CWP_g_num_t **gnum_elt = malloc(sizeof(CWP_g_num_t *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      gnum_elt[i_part] = malloc(sizeof(CWP_g_num_t) * n_elt);
    }

    if (verbose) {
      log_trace("code_id %d\n", code_id);
    }
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_elt; i++) {
        if (code_id == 1) {
          gnum_elt[i_part][i] = i + 1;
        } else {
          gnum_elt[i_part][i] = n_part * n_elt * (rank - 1) + i_part * n_elt + i + 1;
        }
      }
      if (verbose) {
        log_trace("part %d : ", i_part);
        PDM_log_trace_array_long(gnum_elt[i_part], n_elt, "gnum_elt[i_part] : ");
      }
    }

    CWP_Part_data_create(code_name[0],
                         coupling_name,
                         part_data_name,
                         side,
                         gnum_elt,
                         n_elts,
                         n_part);

    // --> exchange
    double **send_data = NULL;
    double **recv_data = NULL;
    int n_comp = 3;

    if (code_id == 1) {

      send_data = malloc(sizeof(double *) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        send_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
      }

      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < n_elt; i++) {
          for (int i_comp = 0; i_comp < n_comp; i_comp++) {
            send_data[i_part][3*i + i_comp] = 0.1*i + i_comp;
          }
        }
      }

      CWP_Part_data_issend(code_name[0],
                           coupling_name,
                           part_data_name,
                           0,
                           sizeof(double),
                           n_comp,
                 (void **) send_data);
    }
    else {

      recv_data = malloc(sizeof(double *) * n_part);
      for (int i_part = 0; i_part < n_part; i_part++) {
        recv_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
      }

      CWP_Part_data_irecv(code_name[0],
                          coupling_name,
                          part_data_name,
                          0,
                          sizeof(double),
                          n_comp,
                (void **) recv_data);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    // --> wait
    if (code_id == 1) {

      CWP_Part_data_wait_issend(code_name[0],
                                coupling_name,
                                part_data_name,
                                0);
    }
    else {

      CWP_Part_data_wait_irecv(code_name[0],
                               coupling_name,
                               part_data_name,
                               0);
    }

    // --> check
    if (verbose) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < n_elt; i++) {
          for (int i_comp = 0; i_comp < n_comp; i_comp++) {
            if (code_id == 1) {
              log_trace("%d - %ld -> s[%d][%d][%d] : %f\n", rank, gnum_elt[i_part][i], i_part, i, i_comp, send_data[i_part][3*i + i_comp]);
            }
            else {
              log_trace("%d - %ld -> r[%d][%d][%d] : %f\n", rank, gnum_elt[i_part][i], i_part, i, i_comp, recv_data[i_part][3*i + i_comp]);
            }
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Delete part_data object
    CWP_Part_data_del(code_name[0],
                      coupling_name,
                      part_data_name);

    // Delete coupling
    CWP_Cpl_del(code_name[0], coupling_name);

    // Free memory
    free(code_name);
    free(coupled_code_name);
    free(intra_comm);

    if (send_data != NULL) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        free(send_data[i_part]);
      }
      free(send_data);
    }
    if (recv_data != NULL) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        if (recv_data[i_part] != NULL) free(recv_data[i_part]);
      }
      free(recv_data);
    }
    if (gnum_elt != NULL) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        if (gnum_elt[i_part] != NULL) free(gnum_elt[i_part]);
      }
      free(gnum_elt);
    }
    if (n_elts != NULL) free(n_elts);

    // Finalize cwipi
    CWP_Finalize();

  }

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

