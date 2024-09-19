/*
  This file is part of the CWIPI library.

  Copyright (C) 2022  ONERA

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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "client_server/client.h"
#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "cwp.h"

#include "cwp_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
 int            argc,
 char         **argv,
 char         **config  // filename for server ip adresses + ports
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

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

/*=============================================================================
 * Main
 *============================================================================*/

int main ( int argc, char *argv[] )
{
  // default
  char *config = NULL;

  _read_args(argc,
             argv,
             &config);

  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  if (config == NULL) {
    if (i_rank == 0) {
      config = (char *) "client_main_o/code1/cwp_config_srv.txt";
    }
    else {
      config = (char *) "client_main_o/code2/cwp_config_srv.txt";
    }
  }

  assert (n_rank == 2);

  // launch server
  if (i_rank == 0) {
    system("mkdir -p client_main_o/code1");
    system("mkdir -p client_main_o/code2");
    system("rm -f ./client_main_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_main_o/code2/cwp_config_srv.txt");
    system("mpirun -n 1 cwp_server -cn code0 -p 60100 60100 -c \"client_main_o/code1/cwp_config_srv.txt\" : -n 1  cwp_server -cn code1 -p 60101 60101 -c \"client_main_o/code2/cwp_config_srv.txt\" &");
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(5);

  // CWP_Init
  int id_code = 0;

  const char  *code_name       = NULL;
  CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

  if (i_rank == 0) {
    id_code = 0;
    code_name = "code1";
  }

  if (i_rank == 1) {
    id_code = 1;
    code_name = "code2";
  }

  // Outputfile
  if (id_code == 0) {
    FILE *f = fopen("output_file_code1.txt", "w");
    CWP_client_Output_file_set(f);
  }

  if (id_code == 1) {
    FILE *f = fopen("output_file_code2.txt", "w");
    CWP_client_Output_file_set(f);
  }

  // EXIT_SUCCESS ?
  int exit_check = 0;

  // Intra-communicator
  MPI_Comm intra_comm;
  MPI_Comm_split(comm, id_code, i_rank, &intra_comm); // smallest i_rank becomes rank 0 in intra_comm

  CWP_client_Init(intra_comm,
                  config,
                  code_name,
                  is_coupled_rank);

  // CWP_User_structure_*
  CWP_client_User_structure_set("code1", NULL);
  CWP_client_User_structure_get("code1");

  // CWP_Codes_*
  int    n_codes   = CWP_client_Codes_nb_get();
  char **codeNames = (char **) CWP_client_Codes_list_get();
  // --> check
  if (!(n_codes == 2 && strcmp(codeNames[0], "code1") == 0 && strcmp(codeNames[1], "code2") == 0)) {
    exit_check = 1;
  }


  // CWP_Loc_codes_*
  int    n_Loc_codes  = CWP_client_Loc_codes_nb_get();
  char **LoccodeNames = (char **) CWP_client_Loc_codes_list_get();
  // --> check
  if (!(n_Loc_codes == 1 && strcmp(LoccodeNames[0], code_name) == 0)) {
    exit_check = 1;
  }

  // Properties_dump
  CWP_client_Properties_dump();

  // State_update
  if (id_code == 0) {
    CWP_client_State_update("code1", CWP_STATE_IN_PROGRESS);
    CWP_State_t state = CWP_client_State_get("code1");
    // --> check
    if (!(state == CWP_STATE_IN_PROGRESS)) {
      exit_check = 1;
    }
  }

  // CWP_Param_*
  int    toto = 42;
  double tata = 0.99;
  double tota = 0.55;
  if (id_code == 0) {
    CWP_client_Param_lock("code1");
    CWP_client_Param_add("code1", "toto", CWP_INT, &toto);
    CWP_client_Param_add("code1", "tata", CWP_DOUBLE, &tata);
    CWP_client_Param_set("code1", "tata", CWP_DOUBLE, &tota);
    CWP_client_Param_unlock("code1");
  }

  const char *A = "Bonjour code 1 !";
  if (id_code == 1) {
    CWP_client_Param_lock("code2");
    CWP_client_Param_add("code2", "toto2", CWP_CHAR, &A);
    CWP_client_Param_unlock("code2");
  }

  MPI_Barrier(comm);

  double titi1;
  CWP_client_Param_get("code1", "tata", CWP_DOUBLE, &titi1);
  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (!(titi1 == tota)) {
    exit_check = 1;
  }
  CWP_GCC_SUPPRESS_WARNING_POP
  const char *titi2 = NULL;
  CWP_client_Param_get("code2", "toto2", CWP_CHAR, &titi2);
  // --> check
  if (!(strcmp(titi2, A) == 0)) {
    exit_check = 1;
  }
  int titi;
  CWP_client_Param_get("code1", "toto", CWP_INT, &titi);
  // --> check
  if (!(titi == toto)) {
    exit_check = 1;
  }

  int code1_n_double = -1;

  if (id_code == 0) {
    int code1_n_int = CWP_client_Param_n_get("code1", CWP_INT);
    // --> check
    if (!(code1_n_int == 2)) {
      exit_check = 1;
    }
    code1_n_double = CWP_client_Param_n_get("code1", CWP_DOUBLE);
    // --> check
    if (!(code1_n_double == 2)) {
      exit_check = 1;
    }
  }

  double tatic = 107.52;
  if (id_code == 0) {
    CWP_client_Param_lock("code1");
    CWP_client_Param_add("code1", "tatic", CWP_DOUBLE, &tatic);
    CWP_client_Param_unlock("code1");
  }

  double totic = 33.50;
  if (id_code == 1) {
    CWP_client_Param_lock("code2");
    CWP_client_Param_add("code2", "tatic", CWP_DOUBLE, &totic);
    CWP_client_Param_unlock("code2");
  }

  if (id_code == 0) {
    CWP_client_Param_lock("code1");
    CWP_client_Param_del("code1", "toto", CWP_INT);
    CWP_client_Param_unlock("code1");
  }

  MPI_Barrier(comm);

  double tita;
  CWP_client_Param_get("code1", "tatic", CWP_DOUBLE, &tita);
  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (!(tita == tatic)) {
    exit_check = 1;
  }
  CWP_GCC_SUPPRESS_WARNING_POP
  double tito;
  CWP_client_Param_get("code2", "tatic", CWP_DOUBLE, &tito);
  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (!(tito == totic)) {
    exit_check = 1;
  }
  CWP_GCC_SUPPRESS_WARNING_POP

  // reduce
  double res = 0;
  char **g_code_names = malloc(sizeof(char *) * 2);
  g_code_names[0] = malloc(sizeof(char) * 99);
  g_code_names[1] = malloc(sizeof(char) * 99);
  char test1[6] = "code1";
  strcpy(g_code_names[0], test1);
  char test2[6] = "code2";
  strcpy(g_code_names[1], test2);
  g_code_names[0][6] = '\0';
  g_code_names[1][6] = '\0';

  CWP_client_Param_reduce(CWP_OP_MAX, "tatic", CWP_DOUBLE, &res, 2, (const char **) g_code_names);

  free(g_code_names[1]);
  free(g_code_names[0]);
  free(g_code_names);

  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (!(res == tatic)) {
    exit_check = 1;
  }
  CWP_GCC_SUPPRESS_WARNING_POP

  if (id_code == 0) {
    char **param_names = NULL;
    CWP_client_Param_list_get("code1", CWP_DOUBLE, &code1_n_double, &param_names);

    // --> check
    if (!(strcmp(param_names[0], "time") == 0 && strcmp(param_names[1], "tata") == 0)
        && strcmp(param_names[1], "tatic") == 0) {
      exit_check = 1;
    }

    for (int i = 0; i < code1_n_double; i++) {
      free(param_names[i]);
    }
    free(param_names);

  }

  if (id_code == 1) {
    int code2_n_char = -1;
    char **param_names = NULL;
    CWP_client_Param_list_get("code2", CWP_CHAR, &code2_n_char, &param_names);

    // --> check
    if (!(strcmp(param_names[0], "toto2") == 0) ){
      exit_check = 1;
    }

    for (int i = 0; i < code2_n_char; i++) {
      free(param_names[i]);
    }
    free(param_names);

  }

  if (id_code == 1) {
    int is_param = CWP_client_Param_is("code2", "toto2", CWP_CHAR);
    // --> check
    if (!(is_param == 1)) {
      exit_check = 1;
    }

    is_param = CWP_client_Param_is("code2", "tambouille", CWP_CHAR);
    // --> check
    if (!(is_param == 0)) {
      exit_check = 1;
    }
  }

  MPI_Barrier(comm);

  char cpl_id1[] = "cpl1_code1_code2";
  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  if (id_code == 0) {
    CWP_client_Cpl_create("code1",
                          cpl_id1,
                          "code2",
                          CWP_INTERFACE_SURFACE,
                          CWP_COMM_PAR_WITH_PART,
                          interp_method,
                          1,
                          CWP_DYNAMIC_MESH_STATIC,
                          CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (id_code == 1) {
    CWP_client_Cpl_create("code2",
                          cpl_id1,
                          "code1",
                          CWP_INTERFACE_SURFACE,
                          CWP_COMM_PAR_WITH_PART,
                          interp_method,
                          1,
                          CWP_DYNAMIC_MESH_STATIC,
                          CWP_TIME_EXCH_USER_CONTROLLED);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  // GLOBAL DATA exchange
  const char *global_data_name  = "andouillette";
  size_t      s_entity          = sizeof(double);
  int         stride            = 2;
  int         n_entity          = 2;

  double *send_data = malloc(s_entity * stride * n_entity);
  send_data[0] = 0.1;
  send_data[1] = 0.2;
  send_data[2] = 0.3;
  send_data[3] = 0.4;

  if (id_code == 0) {
    CWP_client_Global_data_issend("code1",
                                  cpl_id1,
                                  global_data_name,
                                  s_entity,
                                  stride,
                                  n_entity,
                                  (void *) send_data);
  }

  double *recv_data = malloc(s_entity * stride * n_entity);
  if (id_code == 1) {
    CWP_client_Global_data_irecv("code2",
                                 cpl_id1,
                                 global_data_name,
                                 s_entity,
                                 stride,
                                 n_entity,
                       (void **) recv_data);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  if (id_code == 0) {
    CWP_client_Global_data_wait_issend("code1",
                                       cpl_id1,
                                       global_data_name);
  }

  if (id_code == 1) {
    CWP_client_Global_data_wait_irecv("code2",
                                      cpl_id1,
                                      global_data_name);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (id_code == 1) {
    for (int i = 0; i < 4; i++) {
      if (send_data[i] != recv_data[i]) {
        exit_check = 1;
      }
    }
  }
  CWP_GCC_SUPPRESS_WARNING_POP

  // PART DATA exchange
  const char *part_data_name = "Fifi Brindacier";

  int n_part = 1;
  int n_elt = 12;
  int *n_elts = malloc(sizeof(int *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    n_elts[i_part] = n_elt;
  }
  CWP_g_num_t **gnum_elt = malloc(sizeof(CWP_g_num_t *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum_elt[i_part] = malloc(sizeof(CWP_g_num_t) * n_elt);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_elt; i++) {
      gnum_elt[i_part][i] = i + 1;
    }
  }

  if (id_code == 0) {
    CWP_client_Part_data_create("code1",
                                cpl_id1,
                                part_data_name,
                                CWP_PARTDATA_SEND,
                                gnum_elt,
                                n_elts,
                                n_part);
  }

  if (id_code == 1) {
    CWP_client_Part_data_create("code2",
                                cpl_id1,
                                part_data_name,
                                CWP_PARTDATA_RECV,
                                gnum_elt,
                                n_elts,
                                n_part);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  double **send_part_data = NULL;
  double **recv_part_data = NULL;
  int n_comp = 3;

  send_part_data = malloc(sizeof(double *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    send_part_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_elt; i++) {
      for (int i_comp = 0; i_comp < n_comp; i_comp++) {
        send_part_data[i_part][3*i + i_comp] = 0.1*i + i_comp;
      }
    }
  }

  if (id_code == 0) {

    CWP_client_Part_data_issend("code1",
                                cpl_id1,
                                part_data_name,
                                0,
                                sizeof(double),
                                n_comp,
                      (void **) send_part_data);
  }

  recv_part_data = malloc(sizeof(double *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    recv_part_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
  }

  if (id_code == 1) {

    CWP_client_Part_data_irecv("code2",
                               cpl_id1,
                               part_data_name,
                               0,
                               sizeof(double),
                               n_comp,
                     (void **) recv_part_data);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  if (id_code == 0) {

    CWP_client_Part_data_wait_issend("code1",
                                     cpl_id1,
                                     part_data_name,
                                     0);
  }

  if (id_code == 1) {

    CWP_client_Part_data_wait_irecv("code2",
                                    cpl_id1,
                                    part_data_name,
                                    0);
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  // --> check
  CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (id_code == 1) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_elt; i++) {
        for (int i_comp = 0; i_comp < n_comp; i_comp++) {
          if (send_part_data[i_part][3*i + i_comp] != recv_part_data[i_part][3*i + i_comp]) {
            exit_check = 1;
          }
        }
      }
    }
  }
  CWP_GCC_SUPPRESS_WARNING_POP

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  CWP_client_Part_data_del(code_name,
                           cpl_id1,
                           part_data_name);

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  // HO mesh
  int block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                  cpl_id1,
                                                  CWP_BLOCK_FACE_TRIAHO);

  double *node_coord = malloc(sizeof(double) * 3 * 6);
  node_coord[0] = 0.;
  node_coord[1] = 0.;
  node_coord[2] = 0.;

  node_coord[3] = 0.;
  node_coord[4] = 1.;
  node_coord[5] = 0.;

  node_coord[6] = 0.;
  node_coord[7] = 0.;
  node_coord[8] = 1.;

  node_coord[9]  = 0.;
  node_coord[10] = 0.5;
  node_coord[11] = 0.;

  node_coord[12] = 0.;
  node_coord[13] = 0.;
  node_coord[14] = 0.5;

  node_coord[15] = 0.5;
  node_coord[16] = 0.5;
  node_coord[17] = 0.;
  CWP_client_Mesh_interf_vtx_set(code_name,
                                 cpl_id1,
                                 0,
                                 6,
                                 node_coord,
                                 NULL);

  int *elt_node = malloc(sizeof(int) * 6);
  for (int i = 0; i < 6; i++) {
    elt_node[i] = i+1;
  }

  CWP_client_Mesh_interf_block_ho_set(code_name,
                                      cpl_id1,
                                      0,
                                      block_id,
                                      1,
                                      2,
                                      elt_node,
                                      NULL);

  int *ijk = malloc(sizeof(int) * 6 * 2);
  ijk[0] = 0;
  ijk[1] = 0;

  ijk[2] = 2;
  ijk[3] = 0;

  ijk[4] = 0;
  ijk[5] = 2;

  ijk[6] = 1;
  ijk[7] = 0;

  ijk[8] = 0;
  ijk[9] = 1;

  ijk[10] = 1;
  ijk[11] = 1;
  CWP_client_Mesh_interf_ho_ordering_from_IJK_set(code_name,
                                                  cpl_id1,
                                                  CWP_BLOCK_FACE_TRIAHO,
                                                  2,
                                                  6,
                                                  ijk);

  CWP_client_Mesh_interf_finalize(code_name, cpl_id1);

  MPI_Barrier(comm);

  int            Nelts     = 0;
  int            order      = 0;
  int           *connec     = NULL;
  CWP_g_num_t   *global_num = NULL;
  CWP_client_Mesh_interf_block_ho_get(code_name,
                                      cpl_id1,
                                      0,
                                      block_id,
                                      &Nelts,
                                      &order,
                                      &connec,
                                      &global_num);
  // --> check
  if (id_code == 1) {
    if (Nelts != 1) {
      exit_check = 1;
    }
    for (int i = 0; i < 6; i++) {
      if (connec[i] != i+1) {
        exit_check = 1;
      }
    }
  }

  // Barrier on coupling communicator >>>>
  if (id_code == 0) {
    CWP_client_Cpl_barrier("code1",
                           cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Cpl_barrier("code2",
                           cpl_id1);
  }
  // <<<<

  if (id_code == 0) {
    CWP_client_Mesh_interf_del("code1", cpl_id1);
    CWP_client_Cpl_del("code1", cpl_id1);
  }

  if (id_code == 1) {
    CWP_client_Mesh_interf_del("code2", cpl_id1);
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  // CWP_Finalize
  CWP_client_Finalize();

  // free
  free(node_coord);
  free(elt_node);
  free(ijk);
  free(n_elts);
  for (int i_part = 0; i_part < n_part; i_part++) {
    if (gnum_elt[i_part] != NULL) free(gnum_elt[i_part]);
    if (send_part_data[i_part] != NULL) free(send_part_data[i_part]);
    if (recv_part_data[i_part] != NULL) free(recv_part_data[i_part]);
  }
  if (gnum_elt  != NULL) free(gnum_elt);
  if (send_data != NULL) free(send_data);
  if (recv_data != NULL) free(recv_data);
  if (send_part_data != NULL) free(send_part_data);
  if (recv_part_data != NULL) free(recv_part_data);

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  return exit_check;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
