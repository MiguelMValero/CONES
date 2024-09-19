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
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "cwp.h"


/*----------------------------------------------------------------------
 *
 * Main : linear coupling test
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank, comm_world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  bool cond_code1 = rank == 0 || rank == 1 || rank == 2 || rank == 5 || rank == 7;
  bool cond_code2 = rank == 0 || rank == 2 || rank == 6 || rank == 7 || rank == 9;
  bool cond_code3 = rank == 2 || rank == 3 || rank == 4 || rank == 5 || rank == 7 || rank == 9;
  bool cond_code4 = rank == 2 || rank == 4 || rank == 8;

  // Initialization
  int n_code = 0;
  const char **code_names = NULL;
//  double *times_init = NULL;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  if (rank == 0) {
    n_code = 2;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    code_names[1] = "code2";
  }
  else if (rank == 1) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
  }
  else if (rank == 2) {
    n_code = 4;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    code_names[1] = "code2";
    code_names[2] = "code3";
    code_names[3] = "code4";
  }
  else if (rank == 3) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code3";
  }
  else if (rank == 4) {
    n_code = 2;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code3";
    code_names[1] = "code4";
  }
  else if (rank == 5) {
    n_code = 2;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    code_names[1] = "code3";
  }
  else if (rank == 6) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
  }
  else if (rank == 7) {
    n_code = 3;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    code_names[1] = "code2";
    code_names[2] = "code3";
  }
  else if (rank == 8) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code4";
  }
  else if (rank == 9) {
    n_code = 2;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
    code_names[1] = "code3";
  }

  char fileName[19];
  sprintf(fileName, "c_new_api_000%d.txt", rank);
  outputFile = fopen(fileName, "w");

  MPI_Comm *localComm = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_names,
           is_active_rank,
           localComm);

  // Output redirection
  int currentRank;
  int localCommSize;

  for (int i = 0 ; i < n_code ; i++) {
    MPI_Comm_rank(localComm[i], &currentRank);
    MPI_Comm_size(localComm[i], &localCommSize);
  }

  // Finalize
  if (cond_code1) {
    int toto = 111;
    CWP_Param_lock("code1");
    CWP_Param_add("code1", "toto", CWP_INT, &toto);
    const char *A = "Bonjour code 1 !";
    CWP_Param_add("code1", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock("code1");
  }

  if (cond_code2) {
    int toto = 222;
    CWP_Param_lock("code2");
    CWP_Param_add("code2", "toto", CWP_INT, &toto);
    const char *A = "Bonjour code 2!";
    CWP_Param_add("code2", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock("code2");
  }

  if (cond_code3) {
    int toto = 333;
    CWP_Param_lock("code3");
    CWP_Param_add("code3", "toto", CWP_INT, &toto);
    const char *A = "Bonjour code 3!";
    CWP_Param_add("code3", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock("code3");
  }

  if (cond_code4) {
    int toto = 444;
    CWP_Param_lock("code4");
    CWP_Param_add("code4", "toto", CWP_INT, &toto);
    const char *A = "Bonjour code 4!";
    CWP_Param_add("code4", "toto2", CWP_CHAR, &A);
    CWP_Param_unlock("code4");
  }

  int titi;
  char *titi2;

  while (1) {
    if (CWP_Param_is("code4", "toto", CWP_INT) == 1) {
      CWP_Param_get("code4", "toto", CWP_INT, &titi);
      printf("%d : code 4 : toto : %d\n", rank, titi);
      fflush(stdout);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code4", "toto2", CWP_CHAR) == 1) {
      CWP_Param_get("code4", "toto2", CWP_CHAR, &titi2);
      printf("code 4 : toto2 : %s\n", titi2);
      free (titi2);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code3", "toto", CWP_INT) == 1) {
      CWP_Param_get("code3", "toto", CWP_INT, &titi);
      printf("code 3 : toto : %d\n", titi);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code3", "toto2", CWP_CHAR) == 1) {
      CWP_Param_get("code3", "toto2", CWP_CHAR, &titi2);
      printf("code 3 : toto2 : %s\n", titi2);
      free (titi2);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code2", "toto", CWP_INT) == 1) {
      CWP_Param_get("code2", "toto", CWP_INT, &titi);
      printf("code 2 : toto : %d\n", titi);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code2", "toto2", CWP_CHAR) == 1) {
      CWP_Param_get("code2", "toto2", CWP_CHAR, &titi2);
      printf("code 2 : toto2 : %s\n", titi2);
      free (titi2);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code1", "toto", CWP_INT) == 1) {
      CWP_Param_get("code1", "toto", CWP_INT, &titi);
      printf("code 1 : toto : %d\n", titi);
      break;
    }
  }

  while (1) {
    if (CWP_Param_is("code1", "toto2", CWP_CHAR) == 1) {
      CWP_Param_get("code1", "toto2", CWP_CHAR, &titi2);
      printf("code 1 : toto2 : %s\n", titi2);
      free (titi2);
      break;
    }
  }

  assert(titi == 111);

  char cpl_id1[] = "cpl1_code1_code2";
  char cpl_id2[] = "cpl2_code1_code3";
  char cpl_id3[] = "cpl3_code2_code3";
  char cpl_id4[] = "cpl4_code4_code3";
  char cpl_id5[] = "cpl5_code1_code4";
  char cpl_id6[] = "cpl6_code2_code4";

  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  // cpl1: code1 (0, 1, 2, 5, 7) <-> code2 (0, 2, 6, 7, 9)
  if (cond_code1) {
    int v = -1;
    if (rank == 0) {
      v = 11;
    }
    MPI_Bcast(&v, 1, MPI_INT, 0, localComm[0]);
    printf("code 1 v : %d\n", v);
    fflush(stdout);
    CWP_Cpl_create("code1",
                   cpl_id1,
                   "code2",
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }
  if (cond_code2) {
    CWP_Cpl_create("code2",
                   cpl_id1,
                   "code1",
                   CWP_INTERFACE_VOLUME,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
    int v = -2;
    if (rank == 0) {
      v = 21;
    }
    MPI_Comm intraComm;
    if (rank == 0 || rank == 2 || rank == 7) {
      intraComm = localComm[1];
    }
    else if (rank == 6 || rank == 9) {
      intraComm = localComm[0];
    }
    else {
      intraComm = MPI_COMM_NULL;
    }

    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    printf("code 2 v : %d\n", v);
    fflush(stdout);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // cpl2: code1 (0, 1, 2, 5, 7) <-> code3 (2, 3, 4, 5, 7, 9)
  if (cond_code1) {
    CWP_Cpl_create("code1",
                   cpl_id2,
                   "code3",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }
  if (cond_code3) {
    CWP_Cpl_create("code3",
                   cpl_id2,
                   "code1",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
    int v = -2;
    if (rank == 2) {
      v = 31;
    }
    MPI_Comm intraComm;
    if (rank == 2 || rank == 7) {
      intraComm = localComm[2];
    }
    else if (rank == 3 || rank == 4) {
      intraComm = localComm[0];
    }
    else if (rank == 5 || rank == 9) {
      intraComm = localComm[1];
    }
    else {
      intraComm = MPI_COMM_NULL;
    }

    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    printf("code 3 v : %d\n", v);
    fflush(stdout);
  }

  // cpl3: code2 (0, 2, 6, 7, 9) <-> code3 (2, 3, 4, 5, 7, 9)
  if (cond_code2) {
    CWP_Cpl_create("code2",
                   cpl_id3,
                   "code3",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }
  if (cond_code3) {
    CWP_Cpl_create("code3",
                   cpl_id3,
                   "code2",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }

  // cpl4: code4 (2, 4, 8) <-> code3 (2, 3, 4, 5, 7, 9)
  if (cond_code4) {
    CWP_Cpl_create("code4",
                   cpl_id4,
                   "code3",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
    int v = -2;
    if (rank == 2) {
      v = 41;
    }
    MPI_Comm intraComm;
    if (rank == 2) {
      intraComm = localComm[3];
    }
    else if (rank == 4) {
      intraComm = localComm[1];
    }
    else if (rank == 8) {
      intraComm = localComm[0];
    }
    else {
      intraComm = MPI_COMM_NULL;
    }

    MPI_Bcast(&v, 1, MPI_INT, 0, intraComm);
    //    printf("code 4 v : %d\n", v);
    //    fflush(stdout);
  }

  if (cond_code3) {
    CWP_Cpl_create("code3",
                   cpl_id4,
                   "code4",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }

  // cpl5: code1 (0, 1, 2, 5, 7) <-> code4 (2, 4, 8)
  if (cond_code1) {
    CWP_Cpl_create("code1",
                   cpl_id5,
                   "code4",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }
  if (cond_code4) {
    CWP_Cpl_create("code4",
                   cpl_id5,
                   "code1",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }

  // cpl6: code2 (0, 2, 6, 7, 9) <-> code4 (2, 4, 8)
  if (cond_code2) {
    CWP_Cpl_create("code2",
                   cpl_id6,
                   "code4",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }
  if (cond_code4) {
    CWP_Cpl_create("code4",
                   cpl_id6,
                   "code2",
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   interp_method,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }

  // Create Visu
  if (cond_code1) {
    CWP_Visu_set("code1", cpl_id1, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code2) {
    CWP_Visu_set("code2", cpl_id1, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  if (cond_code1) {
    CWP_Visu_set("code1", cpl_id2, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code3) {
    CWP_Visu_set("code3", cpl_id2, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  if (cond_code2) {
    CWP_Visu_set("code2", cpl_id3, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code3) {
    CWP_Visu_set("code3", cpl_id3, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  if (cond_code3) {
    CWP_Visu_set("code3", cpl_id4, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code4) {
    CWP_Visu_set("code4", cpl_id4, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  if (cond_code1) {
    CWP_Visu_set("code1", cpl_id5, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code4) {
    CWP_Visu_set("code4", cpl_id5, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  if (cond_code2) {
    CWP_Visu_set("code2", cpl_id6, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }
  if (cond_code4) {
    CWP_Visu_set("code4", cpl_id6, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
  }

  // Begin time step
  if (cond_code1) {
    CWP_Time_step_beg("code1", 0.0);
  }
  if (cond_code2) {
    CWP_Time_step_beg("code2", 0.0);
  }
  if (cond_code3) {
    CWP_Time_step_beg("code3", 0.0);
  }
  if (cond_code4) {
    CWP_Time_step_beg("code4", 0.0);
  }

  // End time step
  if (cond_code1) {
    CWP_Time_step_end("code1");
  }
  if (cond_code2) {
    CWP_Time_step_end("code2");
  }
  if (cond_code3) {
    CWP_Time_step_end("code3");
  }
  if (cond_code4) {
    CWP_Time_step_end("code4");
  }

  // Delete coupling
  if (cond_code1) {
    CWP_Cpl_del("code1", cpl_id1);
  }
  if (cond_code2) {
    CWP_Cpl_del("code2", cpl_id1);
  }

  if (cond_code1) {
    CWP_Cpl_del("code1", cpl_id2);
  }
  if (cond_code3) {
    CWP_Cpl_del("code3", cpl_id2);
  }

  if (cond_code2) {
    CWP_Cpl_del("code2", cpl_id3);
  }
  if (cond_code3) {
    CWP_Cpl_del("code3", cpl_id3);
  }

  if (cond_code3) {
    CWP_Cpl_del("code3", cpl_id4);
  }
  if (cond_code4) {
    CWP_Cpl_del("code4", cpl_id4);
  }

  if (cond_code1) {
    CWP_Cpl_del("code1", cpl_id5);
  }
  if (cond_code4) {
    CWP_Cpl_del("code4", cpl_id5);
  }

  if (cond_code2) {
    CWP_Cpl_del("code2", cpl_id6);
  }
  if (cond_code4) {
    CWP_Cpl_del("code4", cpl_id6);
  }

  printf("All done for rank %d\n", rank);

  CWP_Finalize();

  MPI_Finalize();

  free(localComm);
  free(code_names);
  fclose(outputFile);

  return 0;
}
