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

#include "pdm_sphere_surf_gen.h"
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
  int                   *swap_codes
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : part_data3
 *
 *---------------------------------------------------------------------*/
int main
(
 int   argc,
 char *argv[]
 )
{
  /* Set default values */
  PDM_g_num_t      n             = 0;
  int              n_part1       = 1;
  int              n_part2       = 1;
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_HILBERT;
  int              disjoint_comm = 0;
  int              verbose       = 0;
  int              swap_codes    = 0;
  _read_args(argc,
             argv,
             &n,
             &n_part1,
             &n_part2,
             &part_method,
             &disjoint_comm,
             &verbose,
             &swap_codes);

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
  const char *cpl_name = "c_new_api_part_data3";

  for (int icode = 0; icode < n_code; icode++) {
    if (verbose) {
      log_trace("Cpl_create %s\n", code_name[icode]);
    }
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   CWP_SPATIAL_INTERP_FROM_IDENTITY, // unused
                   n_part[icode],
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

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
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }



  /* Part data */
  const char *part_data_name1 = "part_data_vtx";
  const char *part_data_name2 = "part_data_face";

  double **send_val1[2][2]; // #icode, #exch
  double **recv_val1[2][2];
  int    **send_val2[2][2]; // #icode, #exch
  int    **recv_val2[2][2];

  for (int icode = 0; icode < n_code; icode++) {

    for (int i = 0; i < 2; i++) {
      send_val1[icode][i] = malloc(sizeof(double *) * n_part[icode]);
      recv_val1[icode][i] = malloc(sizeof(double *) * n_part[icode]);
      send_val2[icode][i] = malloc(sizeof(int    *) * n_part[icode]);
      recv_val2[icode][i] = malloc(sizeof(int    *) * n_part[icode]);
    }

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      send_val1[icode][0][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart] * 3);
      for (int i = 0; i < pn_vtx[icode][ipart] * 3; i++) {
        send_val1[icode][0][ipart][i] = pvtx_coord[icode][ipart][i];
      }

      send_val1[icode][1][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
      for (int i = 0; i < pn_vtx[icode][ipart]; i++) {
        send_val1[icode][1][ipart][i] = pvtx_coord[icode][ipart][3*i];
      }


      send_val2[icode][0][ipart] = malloc(sizeof(int) * pn_face[icode][ipart]);
      for (int i = 0; i < pn_face[icode][ipart]; i++) {
        send_val2[icode][0][ipart][i] = (int) pface_ln_to_gn[icode][ipart][i];
      }

      send_val2[icode][1][ipart] = malloc(sizeof(int) * pn_face[icode][ipart] * 2);
      for (int i = 0; i < pn_face[icode][ipart]; i++) {
        send_val2[icode][1][ipart][2*i  ] = (int)  pface_ln_to_gn[icode][ipart][i];
        send_val2[icode][1][ipart][2*i+1] = (int) -pface_ln_to_gn[icode][ipart][i];
      }

      recv_val1[icode][0][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart] * 3);
      recv_val1[icode][1][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
      recv_val2[icode][0][ipart] = malloc(sizeof(int   ) * pn_face[icode][ipart]);
      recv_val2[icode][1][ipart] = malloc(sizeof(int   ) * pn_face[icode][ipart] * 2);
    }

    if (code_id[icode] == 1) {
      CWP_Part_data_create(code_name[icode],
                           cpl_name,
                           part_data_name1,
                           CWP_PARTDATA_SEND,
                           pvtx_ln_to_gn[icode],
                           pn_vtx       [icode],
                           n_part       [icode]);

      CWP_Part_data_create(code_name[icode],
                           cpl_name,
                           part_data_name2,
                           CWP_PARTDATA_RECV,
                           pface_ln_to_gn[icode],
                           pn_face       [icode],
                           n_part        [icode]);
    }
    else {
      CWP_Part_data_create(code_name[icode],
                           cpl_name,
                           part_data_name1,
                           CWP_PARTDATA_RECV,
                           pvtx_ln_to_gn[icode],
                           pn_vtx       [icode],
                           n_part       [icode]);

      CWP_Part_data_create(code_name[icode],
                           cpl_name,
                           part_data_name2,
                           CWP_PARTDATA_SEND,
                           pface_ln_to_gn[icode],
                           pn_face       [icode],
                           n_part        [icode]);
    }

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
                           part_data_name1,
                           0,
                           sizeof(double),
                           3,
                 (void **) send_val1[icode][0]);

      CWP_Part_data_issend(code_name[icode],
                           cpl_name,
                           part_data_name1,
                           1,
                           sizeof(double),
                           1,
                 (void **) send_val1[icode][1]);

      CWP_Part_data_irecv (code_name[icode],
                           cpl_name,
                           part_data_name2,
                           0,
                           sizeof(int),
                           1,
                 (void **) recv_val2[icode][0]);

      CWP_Part_data_irecv (code_name[icode],
                           cpl_name,
                           part_data_name2,
                           1,
                           sizeof(int),
                           2,
                 (void **) recv_val2[icode][1]);
    }
    else {
      CWP_Part_data_irecv (code_name[icode],
                           cpl_name,
                           part_data_name1,
                           0,
                           sizeof(double),
                           3,
                 (void **) recv_val1[icode][0]);

      CWP_Part_data_irecv (code_name[icode],
                           cpl_name,
                           part_data_name1,
                           1,
                           sizeof(double),
                           1,
                 (void **) recv_val1[icode][1]);

      CWP_Part_data_issend(code_name[icode],
                           cpl_name,
                           part_data_name2,
                           0,
                           sizeof(int),
                           1,
                 (void **) send_val2[icode][0]);

      CWP_Part_data_issend(code_name[icode],
                           cpl_name,
                           part_data_name2,
                           1,
                           sizeof(int),
                           2,
                 (void **) send_val2[icode][1]);
    }
  }


  for (int icode = 0; icode < n_code; icode++) {

    if (code_id[icode] == 1) {
      CWP_Part_data_wait_issend(code_name[icode],
                                cpl_name,
                                part_data_name1,
                                0);

      CWP_Part_data_wait_issend(code_name[icode],
                                cpl_name,
                                part_data_name1,
                                1);

      CWP_Part_data_wait_irecv (code_name[icode],
                                cpl_name,
                                part_data_name2,
                                0);

      CWP_Part_data_wait_irecv (code_name[icode],
                                cpl_name,
                                part_data_name2,
                                1);
    }
    else {
      CWP_Part_data_wait_irecv (code_name[icode],
                                cpl_name,
                                part_data_name1,
                                0);

      CWP_Part_data_wait_irecv (code_name[icode],
                                cpl_name,
                                part_data_name1,
                                1);

      CWP_Part_data_wait_issend(code_name[icode],
                                cpl_name,
                                part_data_name2,
                                0);

      CWP_Part_data_wait_issend(code_name[icode],
                                cpl_name,
                                part_data_name2,
                                1);
    }

  }




  /* Check recv part data */
  int error = 0;

  for (int icode = 0; icode < n_code; icode++) {

    if (code_id[icode] == 1) {

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        for (int i = 0; i < pn_face[icode][ipart]; i++) {
          int expected = (int) pface_ln_to_gn[icode][ipart][i];

          if (recv_val2[icode][0][ipart][i    ] !=  expected ||
              recv_val2[icode][1][ipart][2*i  ] !=  expected ||
              recv_val2[icode][1][ipart][2*i+1] != -expected) {
            error = 1;
            printf("[%d] error for "PDM_FMT_G_NUM" : received %d %d %d, expected %d %d %d\n",
                     i_rank,
                     pvtx_ln_to_gn[icode][ipart][i],
                     recv_val2[icode][0][ipart][i    ],
                     recv_val2[icode][1][ipart][2*i  ],
                     recv_val2[icode][1][ipart][2*i+1],
                     expected,
                     expected,
                     -expected);
              fflush(stdout);
          }
        }
      }
    }

    else {

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        for (int i = 0; i < pn_vtx[icode][ipart]; i++) {
          double expected_x = pvtx_coord[icode][ipart][3*i  ];
          double expected_y = pvtx_coord[icode][ipart][3*i+1];
          double expected_z = pvtx_coord[icode][ipart][3*i+2];

CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
          if (recv_val1[icode][0][ipart][3*i  ] != expected_x ||
              recv_val1[icode][0][ipart][3*i+1] != expected_y ||
              recv_val1[icode][0][ipart][3*i+2] != expected_z ||
              recv_val1[icode][1][ipart][i    ] != expected_x) {
CWP_GCC_SUPPRESS_WARNING_POP
            error = 1;
            printf("[%d] error for "PDM_FMT_G_NUM" : received %f %f %f %f, expected %f %f %f %f\n",
                   i_rank,
                   pvtx_ln_to_gn[icode][ipart][i],
                   recv_val1[icode][0][ipart][3*i  ],
                   recv_val1[icode][0][ipart][3*i+1],
                   recv_val1[icode][0][ipart][3*i+2],
                   recv_val1[icode][1][ipart][i    ],
                   expected_x,
                   expected_x,
                   expected_y,
                   expected_z);
            fflush(stdout);
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
                      part_data_name1);

    CWP_Part_data_del(code_name[icode],
                      cpl_name,
                      part_data_name2);
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

    for (int i = 0; i < 2; i++) {
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        free(send_val1[icode][i][ipart]);
        free(recv_val1[icode][i][ipart]);
        free(send_val2[icode][i][ipart]);
        free(recv_val2[icode][i][ipart]);
      }
      free(send_val1[icode][i]);
      free(recv_val1[icode][i]);
      free(send_val2[icode][i]);
      free(recv_val2[icode][i]);
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
