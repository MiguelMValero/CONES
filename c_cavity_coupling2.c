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
char cl_coupling_name[20], cl_exchange_name[20];
char output_format[20], output_format_option[20];
int il_error = 0;

int nvertex = 9;
int nelts = 4;
int grank;
int rank;
int commWorldSize;
int coord_id;
int nNotLocatedPoints = 0;

int stride;
float geom_tol = 0.1;

double *coords = NULL;
coords = (double*) malloc(sizeof(double) * 27);
int *connecindex = NULL;
connecindex = (int*) malloc(sizeof(int) * (nelts + 1));
int *connec = NULL;
connec = (int*) malloc(sizeof(int) * 16);

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
                      output_format,                            // Postprocessing format
                      output_format_option);

//******************************************************************** End To fill
  
// Mesh definition

if (rank == 0) printf("        Create mesh\n");

// Coordinates table
// int l = 0;
//   for (int j = 0; j < ny; ++j){
//     for (int i = 0; i < nx; ++i){
//       coords[l] = xmin + (xmax - xmin) / (nx - 1) * i; // X
//       coords[l + 1] = ymin + (ymax - ymin) / (ny - 1) * j; // Y
//       coords[l + 2] = 0; // Z
//       l = l + 3;
//     }
//   }

coords[0] = 0, coords[1] = 0, coords[2] = 0;
coords[3] = 1.5, coords[4] = 0, coords[5] = 0;
coords[6] = 3, coords[7] = 0, coords[8] = 0;
coords[9] = 0, coords[10] = 1.5, coords[11] = 0;
coords[12] = 1.5, coords[13] = 1.5, coords[14] = 0;
coords[15] = 3, coords[16] = 1.5, coords[17] = 0;
coords[18] = 0, coords[19] = 3, coords[20] = 0;
coords[21] = 1.5, coords[22] = 3, coords[23] = 0;
coords[24] = 3, coords[25] = 3, coords[26] = 0;

// Connectivity
// connecindex[0] = 0;
//   for (int i = 0; i < nelts; ++i){
//     connecindex[i + 1] = i * 4 + 4;
//   }
connecindex[0] = 0;
connecindex[1] = 4;
connecindex[2] = 8;
connecindex[3] = 12;
connecindex[4] = 16;

// int m = 0;
//   for (int j = 0; j < (ny - 1); ++j){
//     for (int i = 0; i < (nx - 1); ++i){
//       connec[m] = i + j * nx;
//       connec[m + 1] = i + 1 + j * nx;
//       connec[m + 2] = i + 1 + (j + 1) * nx;
//       connec[m + 3] = i + (j + 1) * nx;
//     }
//   }

connec[0] = 0;
connec[1] = 1;
connec[2] = 4;
connec[3] = 3;
connec[4] = 1;
connec[5] = 2;
connec[6] = 5;
connec[7] = 4;
connec[8] = 3;
connec[9] = 4;
connec[10] = 7;
connec[11] = 6;
connec[12] = 4;
connec[13] = 5;
connec[14] = 8;
connec[15] = 7;

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
if (rank == 0) printf("        Exchange Code 1 <-> Code2\n");

printf("I arrived here 1\n");
for (int i = 0; i < nvertex; ++i){
  values[i] = coords[3 * i + coord_id];
  localvalues[i] = 0;
}

stride = 1;
printf("I arrived here 2\n");
sprintf(cl_exchange_name,"exch1");
if (coord_id == 0){
  sprintf(cl_sending_field_name, "coox");
}
else{
  sprintf(cl_sending_field_name, "cooy");
}
sprintf(cl_receiving_field_name,"VarRcv");

printf("I arrived here 3\n");
cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                cl_exchange_name,
                                stride,
                                1,     
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
