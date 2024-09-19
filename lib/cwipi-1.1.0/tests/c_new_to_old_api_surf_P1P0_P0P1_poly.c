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
#include <math.h>
#include <time.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_timer.h"
#include "pdm.h"

#include "cwipi.h"


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_display_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -nx     <val>    Global number of vertices in the side of domaine.\n\n"
         "  -part   <val>    Part of active ranks in the coupling.\n\n"
         "  -s      <val>    Size of domain.\n\n"
         "  -new             Use new CWIPI API.\n\n"
         "  -old             Use old CWIPI API.\n\n"
         "  -h               this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nx             --> Global number of vertices in the side of domaine
 *   part           --> Part of active ranks in the coupling
 *   s              --> Size of domain
 *   new            --> New CWIPI API is used
 *   old            --> Old CWIPI API is used
 *
 *---------------------------------------------------------------------*/

static void
_read_args(int argc, char **argv, PDM_g_num_t *nx, double *part, double *s, int *new, int *old,
           int *n_compute) {
  int i = 1;
  int isNew = 0;
  int isOld = 0;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _display_usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *nx = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-part") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *part = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *s = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-old") == 0) {
      if (isNew != 0) {
        printf("Error : new CWIPI is already selected\n");
        exit(1);
      }
      isOld = 1;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *old = 1;
        *new = 0;
      }
    }
    else if (strcmp(argv[i], "-nc") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *n_compute = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-new") == 0) {
      if (isOld != 0) {
        printf("Error : old CWIPI is already selected\n");
        exit(1);
      }
      isNew = 1;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *new = 1;
        *old = 0;
      }
    }
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
  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  srand(rank + time(0));

  // Read args from command line
  PDM_g_num_t nx = 10;
  double part = 1.;
  double s = 10.;
  int new = 0;
  int old = 1;
  const int nPart = 1;
  int n_compute = 10;
  const double dev_limit = 0.05;

  _read_args(argc, argv, &nx, &part, &s, &new, &old, &n_compute);

  // Init + create coupling
  const char *codeName;
  int codeId;
  const char *codeCoupledName;

  assert (commWorldSize >= 2);

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  assert ((old == 1) || (new == 1));

  MPI_Comm localComm = MPI_COMM_NULL;

  const int nb_part = 1;

  const char *cpl_name;

  if (old) {
    cwipi_init(MPI_COMM_WORLD, codeName, &localComm);

    cwipi_solver_type_t solver_type;

    if (codeId == 1) {
      solver_type = CWIPI_SOLVER_CELL_VERTEX;
    }
    else {
      solver_type = CWIPI_SOLVER_CELL_CENTER;
    }

    cpl_name = "old_cpl";
    cwipi_create_coupling(cpl_name,                                  // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          codeCoupledName,                           // Coupled application id
                          2,                                         // Geometric entities dimension
                          0.01,                                      // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          solver_type,                               // Solver type
                          1,                                         // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");                                   // Postprocessing option
  }
  else {
    const int n_code = 1;
    const CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) &(codeName),
             is_coupled_rank,
             &localComm);

    cpl_name = "new_cpl";
    CWP_Cpl_create(codeName,                                              // Code name
                   cpl_name,                                              // Coupling id
                   codeCoupledName,                                       // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   nb_part,                                               // Partition number
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                   CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

    CWP_Visu_set(codeName,                // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  char gen_name[] = "generator";
  CWP_surf_gen_init(gen_name, (int) nx, (int) nx, nPart, &localComm, part, s, (double) codeId);
  CWP_surf_gen_compute(gen_name);

  int nVtx = 0;
  double *coords = NULL;
  CWP_g_num_t *vtxGnum = NULL;

  int nElts = 0;

  int localRank;
  MPI_Comm_rank(localComm, &localRank);

  PDM_timer_t *timer = PDM_timer_create();
  PDM_timer_t *timer2 = PDM_timer_create();
  PDM_timer_init(timer);

  int n_int = 1;
  double compute_time[n_compute];
  double compute_exch_time[n_int];

  if (old) {
    int *eltsConnecIndex = NULL;
    int *eltsConnec = NULL;
    CWP_g_num_t *eltsGnum = NULL;

    CWP_surf_gen_one_connectivity_get(gen_name,
                                      0,
                                      &nVtx,
                                      &coords,
                                      &vtxGnum,
                                      &nElts,
                                      &eltsConnecIndex,
                                      &eltsConnec,
                                      &eltsGnum);

    cwipi_define_mesh(cpl_name, nVtx, nElts, coords, eltsConnecIndex, eltsConnec);

    MPI_Barrier(MPI_COMM_WORLD);

    PDM_timer_init(timer);
    PDM_timer_init(timer2);

    double mean = 0.0;
    double std_dev = 0.0;
    double mean2;

    int n_it = 0;
    PDM_timer_resume(timer);
    for (int i = 0 ; i < n_compute ; i++) {
      cwipi_update_location(cpl_name);
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);
      cwipi_locate(cpl_name);
      PDM_timer_hang_on(timer2);
      compute_time[i] = PDM_timer_elapsed(timer2);

      mean += compute_time[i];
      MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

      std_dev = 0.0;
      //printf("compute_time[i] %10.5e mmean %10.5e\n",compute_time[i],mean2);
      for (int h = 0 ; h <= i ; h++) {
        std_dev += pow((compute_time[h] - mean2) / mean2, 2);
      }
      std_dev = sqrt(std_dev) / (double) (i + 1);

      double std2;
      MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std_dev = std2 / (double) commWorldSize;
      n_it = i;
      if (i > 3 && std_dev < dev_limit) {
        i = n_compute + 1;
      }

      if (mean2 * n_compute > 600) {
        n_compute = n_compute / 2;
      }

      if (localRank == 0) {
        printf("Survey localization %i %5.4e\n", i, std_dev);
      }
    }
    PDM_timer_hang_on(timer);
    if (localRank == 0) {
      printf("Old localization time %5.4e codeName %s deviation %5.4e nb_it %i\n",
             mean2,
             codeName,
             std_dev,
             n_it);
    }
  }
  else {
    int n_tri = 0;
    int n_quad = 0;
    int n_poly2d = 0;

    int *eltsConnecQuad = NULL;
    CWP_g_num_t *eltsGnumQuad = NULL;

    int *eltsConnecTri = NULL;
    CWP_g_num_t *eltsGnumTri = NULL;

    int *eltsConnecPolyIndex = NULL;
    int *eltsConnecPoly = NULL;
    CWP_g_num_t *eltsGnumPoly = NULL;

    CWP_surf_gen_by_block_get(gen_name,
                              0,
                              &nVtx,
                              &coords,
                              &vtxGnum,
                              &nElts,
                              &n_tri,
                              &eltsConnecTri,
                              &eltsGnumTri,
                              &n_quad,
                              &eltsConnecQuad,
                              &eltsGnumQuad,
                              &n_poly2d,
                              &eltsConnecPolyIndex,
                              &eltsConnecPoly,
                              &eltsGnumPoly);

    CWP_Mesh_interf_vtx_set(codeName, cpl_name, 0, nVtx, coords, vtxGnum);

    int block_id = CWP_Mesh_interf_block_add(codeName, cpl_name, CWP_BLOCK_FACE_TRIA3);

    CWP_Mesh_interf_block_std_set(codeName,
                                  cpl_name,
                                  0,
                                  block_id,
                                  n_tri,
                                  eltsConnecTri,
                                  eltsGnumTri);

    block_id = CWP_Mesh_interf_block_add(codeName, cpl_name, CWP_BLOCK_FACE_QUAD4);

    CWP_Mesh_interf_block_std_set(codeName,
                                  cpl_name,
                                  0,
                                  block_id,
                                  n_quad,
                                  eltsConnecQuad,
                                  eltsGnumQuad);

    block_id = CWP_Mesh_interf_block_add(codeName, cpl_name, CWP_BLOCK_FACE_POLY);

    CWP_Mesh_interf_f_poly_block_set(codeName,
                                     cpl_name,
                                     0,
                                     block_id,
                                     n_poly2d,
                                     eltsConnecPolyIndex,
                                     eltsConnecPoly,
                                     eltsGnumPoly);


    CWP_Mesh_interf_finalize(codeName, cpl_name);
  }

  // Fields exchange
  printf("        Exchange Code1 <-> Code2 %i\n", rank);

  double *sendValues;
  double *recvValues;

  if (codeId == 1) {
    sendValues = (double *) malloc(sizeof(double) * nVtx);
    recvValues = (double *) malloc(sizeof(double) * nElts);
    for (int i = 0 ; i < nVtx ; i++) {
      sendValues[i] = coords[3 * i];
    }
  }

  else {
    sendValues = (double *) malloc(sizeof(double) * nElts);
    recvValues = (double *) malloc(sizeof(double) * nVtx);
    for (int i = 0 ; i < nElts ; i++) {
      sendValues[i] = rank;
    }
  }

  // Exchange
  const char *fieldName1 = "cooX";
  //  const char *fieldName2 = "rank";

  if (old) {
    int sRequest, rRequest;
    int tag = 1;

    double mean = 0.0;
    double std_dev = 0.0;
    double mean2 = 0.0;

    PDM_timer_init(timer);
    PDM_timer_resume(timer);

    for (int i = 0 ; i < n_int ; i++) {
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);

      if (codeId == 1) {
        cwipi_issend(cpl_name, "ech1", tag, 1, 1, 0.1, fieldName1, sendValues, &sRequest);
        cwipi_wait_issend(cpl_name, sRequest);

        //        cwipi_irecv(cpl_name, "ech2", tag, 1, 1, 0.1, fieldName2, recvValues, &rRequest);
        //        cwipi_wait_irecv(cpl_name, rRequest);
      }

      else {
        cwipi_irecv(cpl_name, "ech1", tag, 1, 1, 0.1, fieldName1, recvValues, &rRequest);
        cwipi_wait_irecv(cpl_name, rRequest);

        //        cwipi_issend(cpl_name, "ech2", tag, 1, 1, 0.1, fieldName2, sendValues, &sRequest);
        //        cwipi_wait_issend(cpl_name, sRequest);
      }

      PDM_timer_hang_on(timer2);
      compute_exch_time[i] = PDM_timer_elapsed(timer2);

      mean += compute_exch_time[i];
      MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

      std_dev = 0.0;
      for (int h = 0 ; h <= i ; h++) {
        std_dev += pow((compute_exch_time[h] - mean2) / mean2, 2);
      }
      std_dev = sqrt(std_dev) / (double) (i + 1);
      /*printf("compute_exch_time[i] %10.5e mmean %10.5e std_dev  %10.5e compute_time[i] - mean2 %10.5e\n",
              compute_exch_time[i],mean2,std_dev,compute_time[i] - mean2);
        */
      double std2;
      MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std_dev = std2 / (double) commWorldSize;

      if (i > 3 && std_dev < dev_limit) {
        i = n_int + 1;
      }
      if (mean2 * n_compute > 600) {
        n_compute = n_compute / 2;
      }
      if (localRank == 0) {
        printf("Survey localization %i %5.4e\n", i, std_dev);
      }
    }

    PDM_timer_hang_on(timer);

    if (localRank == 0) {
      printf("Old exchange time for %i iterations %5.4e s codeName %s deviation %5.4e\n",
             0,
             mean2,
             codeName,
             std_dev);
    }
  }
  else {
    CWP_Status_t visu_status = CWP_STATUS_ON;

    if (strcmp(codeName, "code1") == 0) {
      CWP_Field_create(codeName,
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(codeName,
                        0.0);

      CWP_Field_data_set(codeName, cpl_name, fieldName1, 0, CWP_FIELD_MAP_SOURCE, sendValues);
      //      CWP_Field_create(codeName,
      //                       cpl_name,
      //                       fieldName2,
      //                       CWP_DOUBLE,
      //                       CWP_FIELD_STORAGE_BLOCK,
      //                       1,
      //                       CWP_DOF_LOCATION_CELL_CENTER,
      //                       CWP_FIELD_EXCH_RECV,
      //                       visu_status);
      //      CWP_Field_data_set(codeName, cpl_name, fieldName2, 0, recvValues);
    }
    else {
      CWP_Field_create(codeName,
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(codeName,
                        0.0);

      CWP_Field_data_set(codeName, cpl_name, fieldName1, 0, CWP_FIELD_MAP_TARGET, recvValues);
      //      CWP_Field_create(codeName,
      //                       cpl_name,
      //                       fieldName2,
      //                       CWP_DOUBLE,
      //                       CWP_FIELD_STORAGE_BLOCK,
      //                       1,
      //                       CWP_DOF_LOCATION_CELL_CENTER,
      //                       CWP_FIELD_EXCH_SEND,
      //                       visu_status);
      //      CWP_Field_data_set(codeName, cpl_name, fieldName2, 0, CWP_FIELD_MAP_SOURCE, sendValues);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    PDM_timer_init(timer);
    PDM_timer_resume(timer);
    PDM_timer_init(timer);
    PDM_timer_init(timer2);

    PDM_timer_resume(timer);
    double mean = 0.0;
    double mean2 = 0.0;
    double std_dev = 0.0;
    int n_it = 0;
    for (int i = 0 ; i < n_compute ; i++) {
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);
      CWP_Spatial_interp_weights_compute(codeName, cpl_name);;
      PDM_timer_hang_on(timer2);
      compute_time[i] = PDM_timer_elapsed(timer2);

      mean += compute_time[i];
      MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

      std_dev = 0.0;
      for (int h = 0 ; h <= i ; h++) {
        std_dev += pow((compute_time[h] - mean2) / mean2, 2);
      }
      std_dev = sqrt(std_dev) / (double) (i + 1);

      double std2;
      MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std_dev = std2 / (double) commWorldSize;

      n_it = i;
      if (i > 3 && std_dev < dev_limit) {
        i = n_compute + 1;
      }

      if (localRank == 0) {
        printf("Survey exchange %i %5.4e\n", i, std_dev);
      }
    }
    PDM_timer_hang_on(timer);

    if (localRank == 0 && new == 1) {
      printf("New localization time %5.4e codeName %s deviation %5.4e nb_it %i \n",
             mean2,
             codeName,
             std_dev,
             n_it);
    }

    PDM_timer_init(timer);
    PDM_timer_resume(timer);

    mean = 0.0;
    std_dev = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0 ; i < n_int ; i++) {
      PDM_timer_init(timer2);
      PDM_timer_resume(timer2);

      if (strcmp(codeName, "code1") == 0) {
        CWP_Field_issend(codeName, cpl_name, fieldName1);
        CWP_Field_wait_issend(codeName, cpl_name, fieldName1);
        //        CWP_Field_irecv(codeName, cpl_name, fieldName2);
        //        CWP_Field_wait_irecv(codeName, cpl_name, fieldName2);
      }
      else {
        CWP_Field_irecv(codeName, cpl_name, fieldName1);
        CWP_Field_wait_irecv(codeName, cpl_name, fieldName1);
        //        CWP_Field_issend(codeName, cpl_name, fieldName2);
        //        CWP_Field_wait_issend(codeName, cpl_name, fieldName2);
      }

      PDM_timer_hang_on(timer2);
      compute_exch_time[i] = PDM_timer_elapsed(timer2);

      mean += compute_exch_time[i];
      MPI_Allreduce(&mean, &mean2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mean2 = mean2 / ((double) (i + 1) * (double) commWorldSize);

      std_dev = 0.0;
      for (int h = 0 ; h <= i ; h++) {
        std_dev += pow((compute_exch_time[h] - mean2) / mean2, 2);
      }
      std_dev = sqrt(std_dev) / (double) (i + 1);

      double std2;
      MPI_Allreduce(&std_dev, &std2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std_dev = std2 / (double) commWorldSize;

      if (i > 3 && std_dev < dev_limit) {
        i = n_int + 1;
      }
      if (localRank == 0) {
        printf("Survey exchange %i %5.4e\n", i, std_dev);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    PDM_timer_hang_on(timer);

    if (localRank == 0) {
      printf("New exchange time for %i iterations %5.4e s codeName %s deviation %5.4e\n",
             n_int,
             mean2,
             codeName,
             std_dev);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (old) {
    cwipi_delete_coupling(cpl_name);
  }
  else {
    CWP_Time_step_end(codeName);
    CWP_Mesh_interf_del(codeName, cpl_name);
    CWP_Cpl_del(codeName, cpl_name);
  }

  free(sendValues);
  free(recvValues);

  // Finalize
  if (old) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
}
