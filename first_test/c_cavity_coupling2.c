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

#include <cwipi.h>

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
float geom_tol = 1;

double *coords = NULL;
coords = (double*) malloc(sizeof(double) * 3 * nvertex);
int *connecindex = NULL;
connecindex = (int*) malloc(sizeof(int) * (nelts + 1));
int *connec = NULL;
connec = (int*) malloc(sizeof(int) * 16);

double *values = NULL;
values = (double*) malloc(sizeof(double) * 3 * 400);

// MPI Initilization
printf("cav here 1 \n");

MPI_Comm localcomm;
char *codeName;
int codeId;
char *codeCoupledName;
MPI_Init(&argc, &argv);

printf("cav here 2 \n");

MPI_Comm_rank(MPI_COMM_WORLD, &grank);
MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

printf("cav here 3 \n");

codeName = "FOAM_APE";
codeCoupledName = "cwipiFoam";

printf("cav here 4 \n");


//******************************************************************** To fill
// Initialization of the coupling
cwipi_init(MPI_COMM_WORLD,
          codeName,
          &localcomm);

printf("cav here 4.5 \n");

//******************************************************************** End To fill

MPI_Comm_rank(localcomm, &rank);
 
//******************************************************************** To fill
// Create coupling
printf("Create coupling\n");
cwipi_solver_type_t solver_type;
solver_type = CWIPI_SOLVER_CELL_CENTER;
sprintf(cl_coupling_name,"cwipiFoamCoupling");
sprintf(output_format,"EnSight Gold");
sprintf(output_format_option,"text");

printf("cav here 5 \n");

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
                      
printf("cav here 6 \n");

//******************************************************************** End To fill
  
// Mesh definition

if (rank == 0) printf("        Create mesh\n");

coords[0] = 0, coords[1] = 0, coords[2] = 0;
coords[3] = 0.05, coords[4] = 0, coords[5] = 0;
coords[6] = 0.1, coords[7] = 0, coords[8] = 0;
coords[9] = 0, coords[10] = 0.05, coords[11] = 0;
coords[12] = 0.05, coords[13] = 0.05, coords[14] = 0;
coords[15] = 0.1, coords[16] = 0.05, coords[17] = 0;
coords[18] = 0, coords[19] = 0.1, coords[20] = 0;
coords[21] = 0.05, coords[22] = 0.1, coords[23] = 0;
coords[24] = 0.1, coords[25] = 0.1, coords[26] = 0;

connecindex[0] = 0;
connecindex[1] = 4;
connecindex[2] = 8;
connecindex[3] = 12;
connecindex[4] = 16;

connec[0] = 1;
connec[1] = 2;
connec[2] = 5;
connec[3] = 4;
connec[4] = 2;
connec[5] = 3;
connec[6] = 6;
connec[7] = 5;
connec[8] = 4;
connec[9] = 5;
connec[10] = 8;
connec[11] = 7;
connec[12] = 5;
connec[13] = 6;
connec[14] = 9;
connec[15] = 8;

//******************************************************************** To fill
// Define mesh
cwipi_define_mesh(cl_coupling_name,
                  nvertex,
                  nelts,
                  coords,
                  connecindex,
                  connec);

cwipi_locate("cwipiFoamCoupling");
//******************************************************************** End To fill

printf("cav here 7 \n");
 //******************************************************************** To fill

static int recvTag;
char recv_field_name[20];
sprintf(recv_field_name,"recv_field");
static int status;

printf("cav here 7b \n");
cwipi_irecv(cl_coupling_name, "ex1", recvTag, 3, 1, 0, recv_field_name, values, &status);

printf("cav here 8 \n");

switch(status)
{
	case CWIPI_EXCHANGE_OK :
	printf("Receive Ok\n");
	break;
	case CWIPI_EXCHANGE_BAD_RECEIVING :
	printf("Bad receiving\n");
	break;
	default :
	printf("Error : bad exchange status\n");
}

cwipi_wait_irecv(cl_coupling_name,status);

static int sendTag2;
static int status2;

double *fieldsToSend2 = NULL;
fieldsToSend2 = (double*) malloc(sizeof(double) * nvertex);

int j = 0;
for (int i = 0; i < nvertex; ++i){
  fieldsToSend2[i] = coords[j];
  j = j + 3;
}

cwipi_issend("cwipiFoamCoupling","ex2",sendTag2,1,1,0,"u0,v0,w02",fieldsToSend2,&status2);

switch(status2)
  {
    case CWIPI_EXCHANGE_OK :
    printf("Send2 Ok\n");
    break;
    case CWIPI_EXCHANGE_BAD_RECEIVING :
    printf("Bad receiving\n");
    break;
    default :
    printf("Error : bad exchange status\n");
  }

printf("Before wait 2\n");
free(fieldsToSend2);
cwipi_wait_issend("cwipiFoamCoupling",status2);
printf("After wait 2\n");

//******************************************************************** End To fill   

if (rank == 0) printf("        Delete coupling\n");

//******************************************************************** To fill
// Delete coupling
cwipi_delete_coupling(cl_coupling_name);

printf("cav here 9 \n");

//******************************************************************** End To fill

// Freeing memory
free(coords);
printf("cav here 9a \n");
free(connecindex);
printf("cav here 9b \n");
free(connec);
printf("cav here 9c \n");
free(values);

printf("cav here 10 \n");
//******************************************************************** To fill
// Finalize
cwipi_finalize();
MPI_Finalize();

//******************************************************************** End To fill
  return EXIT_SUCCESS;
}
