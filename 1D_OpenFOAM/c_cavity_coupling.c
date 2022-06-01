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
int itdeb;
int itend;
float ttime;
float dt;

int stride;
float geom_tol = 0.001;

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
  coord_id = 0; // coordinate to send (X)
  codeCoupledName = "code2";
}
else {
  codeName = "code2";
  coord_id = 1; // coordinate to send (Y)
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


itdeb = 1; // Initial iteration
itend = 84; // Final iteration
ttime = 0.2; // Physical time
dt = 0.2; // time step

for (int it = itdeb; it <= itend; ++it) {


//******************************************************************** To fill
// Define mesh
if (it == itdeb) {
  cwipi_define_mesh(cl_coupling_name,
                  nvertex,
                  nelts,
                  coords,
                  connecindex,
                  connec);
}
else {
  cwipi_update_location(cl_coupling_name);
    }

//******************************************************************** End To fill

//******************************************************************** To fill
// Send receive
if (rank == 0) printf("        Exchange\n");

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

cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                cl_exchange_name,
                                stride,
                                it,     // n_step
                                ttime,   // physical_time
                                cl_sending_field_name,
                                values,
                                cl_receiving_field_name,
                                localvalues,
                                &nNotLocatedPoints);

//******************************************************************** End To fill
  ttime = ttime + dt;
}

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
