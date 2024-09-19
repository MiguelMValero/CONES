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

#include "pdm_sphere_surf_gen.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#define ABS(a) ((a) < 0 ? -(a) : (a))

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
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  PDM_g_num_t           *n,
  int                   *n_part1,
  int                   *n_part2,
  PDM_split_dual_t      *part_method,
  int                   *disjoint_comm,
  int                   *verbose,
  int                   *swap_codes,
  int                   *exchange_fields
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
        long _n = atol(argv[i]);
        *n = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part2 = atoi(argv[i]);
      }
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
    else if (strcmp(argv[i], "-disjoint") == 0) {
      *disjoint_comm = 1;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-swap_codes") == 0) {
      *swap_codes = 1;
    }
    else if (strcmp(argv[i], "-no_exchange") == 0) {
      *exchange_fields = 0;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : part_data
 *
 *---------------------------------------------------------------------*/
int main
(
 int   argc,
 char *argv[]
 )
{
  /* Set default values */
  PDM_g_num_t      n               = 3;
  int              n_part1         = 1;
  int              n_part2         = 1;
  PDM_split_dual_t part_method     = PDM_SPLIT_DUAL_WITH_HILBERT;
  int              disjoint_comm   = 0;
  int              verbose         = 0;
  int              swap_codes      = 0;
  int              exchange_fields = 1;
  _read_args(argc,
             argv,
             &n,
             &n_part1,
             &n_part2,
             &part_method,
             &disjoint_comm,
             &verbose,
             &swap_codes,
             &exchange_fields);

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

  /* Initialize CWIPI */
  const char *all_code_names[2] = {"code1", "code2"};
  int all_n_part[2] = {n_part1, n_part2};

  int has_code[2] = {0, 0};
  if (disjoint_comm) {
    has_code[0] = i_rank < n_rank/2;
    has_code[1] = !has_code[0];
  }
  else {
    has_code[0] = (i_rank < (2*n_rank) / 3);
    has_code[1] = (i_rank >= n_rank / 3);
  }

  if (swap_codes) {
    int tmp = has_code[0];
    has_code[0] = has_code[1];
    has_code[1] = tmp;
  }

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  int           *n_part            = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_part           [n_code] = all_n_part    [icode];

      if (verbose) {
        log_trace("Running %s, coupled with %s, n_part = %d\n",
                  code_name[n_code], coupled_code_name[n_code], n_part[n_code]);
      }
      n_code++;
    }
  }

  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }

  /* Create coupling */
  const char *cpl_name = "c_new_api_part_data2";

  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   spatial_interp,
                   n_part[icode],
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

  }

  // !! this must be performed in 2 separate loops if the intra comms do overlap
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Visu_set(code_name[icode],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }


  /* Define mesh */
  int          **pn_face        = malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = malloc(sizeof(double      **) * n_code);
  PDM_g_num_t ***pface_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]);

    PDM_sphere_surf_icosphere_gen_part(mesh_comm,
                                       n,
                                       0.,
                                       0.,
                                       0.,
                                       1.,
                                       n_part[icode],
                                       part_method,
                                       &pn_vtx        [icode],
                                       &pvtx_coord    [icode],
                                       &pvtx_ln_to_gn [icode],
                                       &pn_face       [icode],
                                       &pface_vtx_idx [icode],
                                       &pface_vtx     [icode],
                                       &pface_ln_to_gn[icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx       [icode][ipart],
                              pvtx_coord   [icode][ipart],
                              pvtx_ln_to_gn[icode][ipart]);

      CWP_Mesh_interf_f_poly_block_set(code_name[icode],
                                       cpl_name,
                                       ipart,
                                       block_id,
                                       pn_face       [icode][ipart],
                                       pface_vtx_idx [icode][ipart],
                                       pface_vtx     [icode][ipart],
                                       pface_ln_to_gn[icode][ipart]);
    }

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  /* Create and set fields */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "field";

  int stride = 1;
  double **send_val  = NULL;
  double **recv_val  = NULL;
  double **field_ptr = NULL;

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Field_exch_t exch_type;
    CWP_Field_map_t  map_type;
    if (code_id[icode] == 1) {
      send_val = malloc(sizeof(double *) * n_part[icode]);
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        send_val[ipart] = malloc(sizeof(double) * pn_face[icode][ipart] * stride);
        for (int i = 0; i < pn_face[icode][ipart]; i++) {
          for (int j = 0; j < stride; j++) {
            send_val[ipart][stride*i+j] = (double) (j+1)*pface_ln_to_gn[icode][ipart][i];
          }
        }
      }
      exch_type = CWP_FIELD_EXCH_SEND;
      map_type  = CWP_FIELD_MAP_SOURCE;
      field_ptr = send_val;
    }
    else {
      recv_val = malloc(sizeof(double *) * n_part[icode]);
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        recv_val[ipart] = malloc(sizeof(double) * pn_face[icode][ipart] * stride);
      }
      exch_type = CWP_FIELD_EXCH_RECV;
      map_type  = CWP_FIELD_MAP_TARGET;
      field_ptr = recv_val;
    }

    CWP_Field_create(code_name[icode],
                     cpl_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     stride,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     exch_type,
                     visu_status);

    CWP_Time_step_beg(code_name[icode],
                      0.0);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Field_data_set(code_name[icode],
                         cpl_name,
                         field_name,
                         ipart,
                         map_type,
                         field_ptr[ipart]);
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create fields OK\n");
    fflush(stdout);
  }

  /* Map source to target */
  for (int icode = 0; icode < n_code; icode++) {
    if (spatial_interp == CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES) {
      CWP_Spatial_interp_property_set(code_name[icode],
                                      cpl_name,
                                      "n_neighbors",
                                      CWP_INT,
                                      "1");
    }
    CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
  }


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Interpolation weights computation OK\n");
    fflush(stdout);
  }


  int error = 0;
  if (exchange_fields) {
    for (int icode = 0; icode < n_code; icode++) {
      if (code_id[icode] == 1) {
        CWP_Field_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_irecv (code_name[icode], cpl_name, field_name);
      }

      if (code_id[icode] == 1) {
        CWP_Field_wait_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name);
      }
    }

    MPI_Barrier(comm);
    if (i_rank == 0) {
      printf("Exchange fields OK\n");
      fflush(stdout);
    }


    /* Check recv field */
    for (int icode = 0; icode < n_code; icode++) {
      if (code_id[icode] == 1) {
        if (verbose) {
          for (int ipart = 0; ipart < n_part[icode]; ipart++) {
            log_trace("\n-- src part %d --\n", ipart);
            for (int i = 0; i < pn_face[icode][ipart]; i++) {
              log_trace(PDM_FMT_G_NUM" sends:\n", pface_ln_to_gn[icode][ipart][i]);
              for (int j = 0; j < stride; j++) {
                log_trace("  %f\n",
                          send_val[ipart][stride*i+j]);
              }
            }
          }
        }
      }
      else {
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {
          if (verbose) {
            log_trace("\n-- tgt part %d --\n", ipart);
          }
          for (int i = 0; i < pn_face[icode][ipart]; i++) {
            if (verbose) {
              log_trace(PDM_FMT_G_NUM" received:\n", pface_ln_to_gn[icode][ipart][i]);
            }
            for (int j = 0; j < stride; j++) {
              double expected = (double) (j+1)*pface_ln_to_gn[icode][ipart][i];
              if (verbose) {
                log_trace("  %f (expected %f)\n",
                          recv_val[ipart][stride*i+j],
                          expected);
              }

              if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
                if (ABS(recv_val[ipart][stride*i+j] - expected) > 1e-12) {
                  error = 1;
                  printf("[%d] error for "PDM_FMT_G_NUM" : received %e, expected %e\n",
                         i_rank, pface_ln_to_gn[icode][ipart][i],
                         recv_val[ipart][stride*i+j], expected);
                  fflush(stdout);
                }
              }
            }
          }
        }
      }
    }

    MPI_Barrier(comm);
    if (i_rank == 0) {
      printf("Check fields OK\n");
      fflush(stdout);
    }
  }



  /* Part data */
  const char *part_data_name = "part_data";

  for (int icode = 0; icode < n_code; icode++) {
    CWP_PartData_exch_t exch_type;
    if (code_id[icode] == 1) {
      exch_type = CWP_PARTDATA_SEND;
    }
    else {
      // reset recv data
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        for (int i = 0; i < pn_face[icode][ipart] * stride; i++) {
          recv_val[ipart][i] = -1234;
        }
      }
      exch_type = CWP_PARTDATA_RECV;
    }

    CWP_Part_data_create(code_name[icode],
                         cpl_name,
                         part_data_name,
                         exch_type,
                         pface_ln_to_gn[icode],
                         pn_face       [icode],
                         n_part        [icode]);
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create part data OK\n");
    fflush(stdout);
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Part_data_issend(code_name[icode],
                           cpl_name,
                           part_data_name,
                           0,
                           sizeof(double),
                           stride,
                 (void **) send_val);
    }
    else {
      CWP_Part_data_irecv(code_name[icode],
                          cpl_name,
                          part_data_name,
                          0,
                          sizeof(double),
                          stride,
                (void **) recv_val);
    }
  }



  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Part_data_wait_issend(code_name[icode],
                                cpl_name,
                                part_data_name,
                                0);
    }
    else {
      CWP_Part_data_wait_irecv(code_name[icode],
                               cpl_name,
                               part_data_name,
                               0);
    }
  }

  /* Check recv part data */
  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 2) {
      if (verbose) {
        log_trace("\n--- PartData ---\n");
      }
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        if (verbose) {
          log_trace("\n-- tgt part %d --\n", ipart);
        }
        for (int i = 0; i < pn_face[icode][ipart]; i++) {
          if (verbose) {
            log_trace(PDM_FMT_G_NUM" received:\n", pface_ln_to_gn[icode][ipart][i]);
          }
          for (int j = 0; j < stride; j++) {
            double expected = (double) (j+1)*pface_ln_to_gn[icode][ipart][i];
            if (verbose) {
              log_trace("  %f (expected %f)\n",
                        recv_val[ipart][stride*i+j],
                        expected);
            }

            if (ABS(recv_val[ipart][stride*i+j] - expected) > 0) {
              error = 1;
              printf("[%d] error for "PDM_FMT_G_NUM" : received %e, expected %e\n",
                     i_rank, pface_ln_to_gn[icode][ipart][i],
                     recv_val[ipart][stride*i+j], expected);
              fflush(stdout);
            }
          }
        }
      }
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Check received part data OK\n");
    fflush(stdout);
  }


  for (int icode = 0; icode < n_code; icode++) {
    CWP_Part_data_del(code_name[icode],
                      cpl_name,
                      part_data_name);
  }


  /* Global data */
  const char *global_data_name = "global_data";
  int global_stride   = 3;
  int global_n_entity = 4;
  int *global_data = malloc(sizeof(int) * global_stride * global_n_entity);

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Global_data_irecv(code_name[icode],
                            cpl_name,
                            global_data_name,
                            sizeof(int),
                            global_stride,
                            global_n_entity,
                            global_data);
    }
    else {
      for (int i = 0; i < global_n_entity; i++) {
        for (int j = 0; j < global_stride; j++) {
          global_data[global_stride*i + j] = (i+1) * (j+1);
        }
      }

      CWP_Global_data_issend(code_name[icode],
                             cpl_name,
                             global_data_name,
                             sizeof(int),
                             global_stride,
                             global_n_entity,
                             global_data);
    }
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      if (verbose) {
        log_trace("\n--- GlobalData ---\n");
      }
      CWP_Global_data_wait_irecv(code_name[icode],
                                 cpl_name,
                                 global_data_name);
      for (int i = 0; i < global_n_entity; i++) {
        if (verbose) {
          log_trace("global entity %d received ", i);
          PDM_log_trace_array_int(global_data + global_stride*i,
                                  global_stride,
                                  "");
        }
        for (int j = 0; j < global_stride; j++) {
          int expected = (i+1) * (j+1);
          if (global_data[global_stride*i + j] != expected) {
            error = 1;
            printf("[%d] error global entity %d comp %d : received %d (expected %d)\n",
                   i_rank, i, j, global_data[global_stride*i + j], expected);
            fflush(stdout);
          }
        }
      }
    }
    else {
      CWP_Global_data_wait_issend(code_name[icode],
                                  cpl_name,
                                  global_data_name);
    }
  }


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Check received global data OK\n");
    fflush(stdout);
  }



  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("End\n");
    fflush(stdout);
  }

  /* Free memory */
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Time_step_end(code_name[icode]);

    CWP_Mesh_interf_del(code_name[icode], cpl_name);

    CWP_Cpl_del(code_name[icode], cpl_name);
  }

  free(global_data);

  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      free(pface_vtx_idx [icode][ipart]);
      free(pface_vtx     [icode][ipart]);
      free(pvtx_coord    [icode][ipart]);
      free(pface_ln_to_gn[icode][ipart]);
      free(pvtx_ln_to_gn [icode][ipart]);
    }
    free(pn_face       [icode]);
    free(pn_vtx        [icode]);
    free(pface_vtx_idx [icode]);
    free(pface_vtx     [icode]);
    free(pvtx_coord    [icode]);
    free(pface_ln_to_gn[icode]);
    free(pvtx_ln_to_gn [icode]);

    if (code_id[icode] == 1) {
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        free(send_val[ipart]);
      }
      free(send_val);
    }
    else {
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        free(recv_val[ipart]);
      }
      free(recv_val);
    }
  }
  free(pn_face       );
  free(pn_vtx        );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pvtx_coord    );
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );

  free(code_id);
  free(n_part);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  /* Finalize CWIPI */
  CWP_Finalize();

  /* Finalize MPI */
  MPI_Finalize();

  // return error;
  return error;
}
