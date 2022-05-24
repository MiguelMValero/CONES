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

int nvertex;
int nelts;
int grank;
int rank;
int commWorldSize;
int coord_id;
int nNotLocatedPoints = 0;

int stride;
float geom_tol = 0.1;

// MPI Initilization
MPI_Comm localcomm;
char *codeName;
int codeId;
char *codeCoupledName;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &grank);
MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

codeName = "code2";
coord_id = 1; // coordinate to send (Y coordinate, since it is a 1D problem)
codeCoupledName = "code1";

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
                      3,                                         // Geometric entities dimension
                      geom_tol,                                  // Geometric tolerance
                      CWIPI_STATIC_MESH,                         // Mesh type
                      solver_type,                               // Solver type
                      1,                                         // Postprocessing frequency
                      output_format,                             // Postprocessing format
                      output_format_option);

//******************************************************************** End To fill
