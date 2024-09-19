/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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
#include <stddef.h>
#include <string.h>

#include <mpi.h>

#include "cwipi.h"

/*----------------------------------------------------------------------
 *
 * Dump status exchange
 *
 * parameters:
 *   status              <-- Exchange status
 *---------------------------------------------------------------------*/

static void _dumpStatus(FILE* outputFile, cwipi_exchange_status_t status)
{
  switch(status) {
  case CWIPI_EXCHANGE_OK :
    fprintf(outputFile, "Exchange Ok\n");
    break;
  case CWIPI_EXCHANGE_BAD_RECEIVING :
    fprintf(outputFile, "Bad receiving\n");
    break;
  default :
    printf("Error : bad exchange status\n");
    exit(1);
  }
}

/*----------------------------------------------------------------------
 *
 * Dump not located points
 *
 * parameters:
 *   coupling_id         <-- Coupling id
 *   nNotLocatedPoints   <-- Number of not located points
 *---------------------------------------------------------------------*/

static void _dumpNotLocatedPoints(FILE* outputFile,
                                  const char *coupling_id,
                                  const int nNotLocatedPoints)
{
  if ( nNotLocatedPoints > 0) {
    fprintf(outputFile, "Not located points :\n");
    const int* notLocatedPoints = cwipi_get_not_located_points(coupling_id);
    for(int i = 0; i < nNotLocatedPoints; i++)
     fprintf(outputFile, "%i ", notLocatedPoints[i]);
    fprintf(outputFile, "\n");
  }
}

/*----------------------------------------------------------------------
 *
 * Main : linear coupling test
 *
 *---------------------------------------------------------------------*/

int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  char *srcName = (char *) malloc (sizeof(char) * (strlen(__FILE__) + 1));
  strcpy(srcName, __FILE__);
  char *srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '.');
  if (srcBaseName != NULL)
    *srcBaseName = '\0';
  srcBaseName = NULL;
  srcBaseName = strrchr(srcName, '/');
  if (srcBaseName != NULL)
    srcBaseName += 1;
  else
    srcBaseName = srcName;

  if (rank == 0)
    printf("\nSTART: %s\n", srcBaseName);

  if (comm_world_size != 2) {
    if (rank == 0)
      printf("      Not executed : only available for 2 processus\n");
    MPI_Finalize();
    return 1;
  }

  /* Initialization
   * -------------- */

  const char *codeName;
  const char *codeCoupledName;

  if (rank == 0) {
    codeName="code1";
    codeCoupledName="code2";
  }
  else {
    codeName="code2";
    codeCoupledName="code1";
  }

  const char* fileName = NULL;
  if (rank == 0)
    fileName="c_linear_coupling_0000.txt";
  else
    fileName="c_linear_coupling_0001.txt";

  outputFile = fopen(fileName,"w");

  cwipi_set_output_listing(outputFile);

  MPI_Comm localComm;
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localComm);

  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);

  fprintf(outputFile, "Linear coupling P0 <-> P1\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");
  cwipi_dump_application_properties();


  if (rank == 0)
      printf("        Create coupling\n");

  cwipi_solver_type_t solver_type;

  if (rank == 0)
    solver_type = CWIPI_SOLVER_CELL_CENTER;
  else
    solver_type = CWIPI_SOLVER_CELL_VERTEX;

  /* Coupling creation
   * ----------------- */

  cwipi_create_coupling("test2D_0",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        1,                                         // Geometric entities dimension
                        100,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        1,                                         // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option

  /* Mesh definition
   * --------------- */

  if (rank == 0)
    printf("        Create mesh\n");

  int nVertex = 3;
  int nElts = 2;

  double *coords = (double *) malloc(3 * nVertex * sizeof(double));
  coords[0] = 0.;
  coords[1] = 0.;
  coords[2] = 0.;

  if (rank == 0)
    coords[3] = 1.;
  else
    coords[3] = 0.5;
  coords[4] = 0.;
  coords[5] = 0.;

  coords[6] = 2.;
  coords[7] = 0.;
  coords[8] = 0.;

  int *connecIdx = (int *) malloc((nElts + 1) * sizeof(int));
  connecIdx[0] = 0;
  connecIdx[1] = 2;
  connecIdx[2] = 4;

  int *connec = (int *) malloc(connecIdx[nElts] * sizeof(int));
  connec[0] = 1;
  connec[1] = 2;

  connec[2] = 2;
  connec[3] = 3;


  cwipi_define_mesh("test2D_0",
                    nVertex,
                    nElts,
                    coords,
                    connecIdx,
                    connec);


  /* Fields exchange
   * --------------- */

  double *values = NULL;
  double *values1 = NULL;

  if (rank == 0) {
    values = (double *) malloc(nElts * sizeof(double));
    values1 = (double *) malloc(nElts * sizeof(double));

    values[0] = (coords[0]+coords[3])/2.;
    values[1] = (coords[3]+coords[6])/2.;

    values1[0] = 0.;
    values1[1] = 0.;
  }
  else {
    values = (double *) malloc(nVertex * sizeof(double));
    values1 = (double *) malloc(nVertex * sizeof(double));

    values[0] = coords[0];
    values[1] = coords[3];
    values[2] = coords[6];

    values1[0] = 0.;
    values1[1] = 0.;
    values1[2] = 0.;
  }

  if (rank == 0)
    printf("        Exchange\n");

  int nNotLocatedPoints = 0;
  cwipi_exchange_status_t status = cwipi_exchange("test2D_0",
                                                  "echange1",
                                                  1,
                                                  1,     // n_step
                                                  0.1,   // physical_time
                                                  "cooX",
                                                  values,
                                                  "cooX",
                                                  values1,
                                                  &nNotLocatedPoints);

  _dumpStatus(outputFile, status);
  _dumpNotLocatedPoints(outputFile, "test2D_0", nNotLocatedPoints);

  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("test2D_0");

  /* Free
   * ---- */

  free(coords);
  free(connecIdx);
  free(connec);
  free(values);
  free(values1);
  free(srcName);

  /* Finalize
   * -------- */

  cwipi_finalize();

  MPI_Finalize();

  fclose(outputFile);

  return 0;
}
