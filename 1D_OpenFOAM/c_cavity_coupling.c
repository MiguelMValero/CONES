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
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwipi.h"

char cl_coupling_name[20], cl_exchange_name[20];
char output_format[20], output_format_option[20];
int il_error = 0;
/*----------------------------------------------------------------------
 *                                                                     
 * Main : surface coupling test :  
 *
 *---------------------------------------------------------------------*/
 
int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

char cl_sending_field_name[20], cl_receiving_field_name[20];

int nvertex = 11;
int nelts = 5;
int grank;
int rank;
int commWorldSize;
int coord_id;
int nNotLocatedPoints = 0;

int stride;
float geom_tol = 0.1;

double *coords = NULL;
coords = (double*) malloc(sizeof(double) * 3 * nvertex);
int *connecindex = NULL;
connecindex = (int*) malloc(sizeof(int) * (nelts + 1));
int *connec = NULL;
connec = (int*) malloc(sizeof(int) * 21);

double *values = NULL;
values = (double*) malloc(sizeof(double) * nvertex);
double *localvalues = NULL;
localvalues = (double*) malloc(sizeof(double) * nvertex);

// MPI Initilization
MPI_Comm localcomm;
char *codeName;
int codeId;
char *codeCoupledName;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &grank);
MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

if (grank < commWorldSize / 2) {
  codeName = "code1";
  coord_id = 0; // coordinate to send
  codeCoupledName = "code2";
}
else {
  codeName = "code2";
  coord_id = 1; // coordinate to send
  codeCoupledName = "code1";
}

//******************************************************************** To fill
// Initialization of the coupling
cwipi_init(MPI_COMM_WORLD,
          codeName,
          &localcomm);

//******************************************************************** End To fill

MPI_Comm_rank(localcomm, &rank);

//******************************************************************** To fill
// Create coupling
if (rank == 0) printf("        Create coupling\n");
cwipi_solver_type_t solver_type;
solver_type = CWIPI_SOLVER_CELL_VERTEX;
sprintf(cl_coupling_name,"cpl1");
sprintf(output_format,"EnSight Gold");
sprintf(output_format_option,"text");
cwipi_create_coupling(cl_coupling_name,
                      CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                      codeCoupledName,                           // Coupled application id
                      2,                                         // Geometric entities dimension
                      geom_tol,                                  // Geometric tolerance
                      CWIPI_STATIC_MESH,                         // Mesh type
                      solver_type,                               // Solver type
                      1,                                         // Postprocessing frequency
                      output_format,                             // Postprocessing format
                      output_format_option);

//******************************************************************** End To fill
  
// Mesh definition
if (rank == 0) printf("        Create mesh\n");

coords[0] = 0, coords[1] = 0, coords[2] = 0;
coords[3] = 1, coords[4] = 0, coords[5] = 0;
coords[6] = 2, coords[7] = 0, coords[8] = 0;
coords[9] = 3, coords[10] = 0, coords[11] = 0;
coords[12] = 0, coords[13] = 1, coords[14] = 0;
coords[15] = 2, coords[16] = 1, coords[17] = 0;
coords[18] = 3, coords[19] = 1, coords[20] = 0;
coords[21] = 1, coords[22] = 2, coords[23] = 0;
coords[24] = 0, coords[25] = 3, coords[26] = 0;
coords[27] = 2, coords[28] = 3, coords[29] = 0;
coords[30] = 3, coords[31] = 3, coords[32] = 0;

connecindex[0] = 0;
connecindex[1] = 3;
connecindex[2] = 7;
connecindex[3] = 11;
connecindex[4] = 16;
connecindex[5] = 21;

connec[0] = 1, connec[1] = 2, connec[2] = 5;
connec[3] = 3, connec[4] = 4, connec[5] = 7, connec[6] = 6;
connec[7] = 5, connec[8] = 8, connec[9] = 10, connec[10] = 9;
connec[11] = 5, connec[12] = 2, connec[13] = 3, connec[14] = 6, connec[15] = 8;
connec[16] = 6, connec[17] = 7, connec[18] = 11, connec[19] = 10, connec[20] = 8;

//******************************************************************** To fill
// Define mesh
cwipi_define_mesh(cl_coupling_name,
                  nvertex,
                  nelts,
                  coords,
                  connecindex,
                  connec);

//******************************************************************** End To fill

//******************************************************************** To fill
// Send receive
if (rank == 0) printf("        Exchange\n");

printf("I arrived here 1\n");
for (int i = 0; i < nvertex; ++i){
  values[i] = coords[3 * i + coord_id];
}
  
stride = 1;

sprintf(cl_exchange_name, "exch1");
if (coord_id == 0){
  sprintf(cl_sending_field_name, "coox");
}
else{
  sprintf(cl_sending_field_name, "cooy");
}
sprintf(cl_receiving_field_name, "recv");

printf("I arrived here 2\n");
cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                cl_exchange_name,
                                stride,
                                1,     // n_step
                                0.1,   // physical_time
                                cl_sending_field_name,
                                values,
                                cl_receiving_field_name,
                                localvalues,
                                &nNotLocatedPoints);

//******************************************************************** End To fill

if (rank == 0) printf("        Delete coupling\n");

//******************************************************************** To fill
// Delete coupling
cwipi_delete_coupling(cl_coupling_name);

//******************************************************************** End To fill

// Freeing memory
free(coords);
free(connecindex);
free(connec);
free(values);
free(localvalues);

//******************************************************************** To fill
// Finalize
cwipi_finalize();
MPI_Finalize();

//******************************************************************** End To fill
  return EXIT_SUCCESS;
}
