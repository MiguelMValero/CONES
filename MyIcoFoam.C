/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwipi.h"

int il_err = 0;
char cla_obj[20], cla_space[20], cl_coupling_name[20];
char output_format[20], output_format_option[20], cl_exchange_name[20];
char cl_varsend[20], cl_varrec[20];

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

  int grank;
  int rank;
  int commWorldSize;
  int stride;
  int itdeb;
  int itend;
  float geom_tol;
  float ttime;
  float dt;
  int it;

  MPI_Comm localcomm;
  char *codeName;
  int codeId;
  char *codeCoupledName;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &grank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  codeName = "code2";
  codeCoupledName = "code1";

  // Initialization of the coupling
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localcomm);

  MPI_Comm_rank(localcomm, &rank);

const int nElemSeg = 2;
const int nVertexSeg = nElemSeg + 1;

// Create coupling
if (rank == 0) printf("        Create coupling\n");
cwipi_solver_type_t solver_type;
solver_type = CWIPI_SOLVER_CELL_CENTER;
sprintf(cl_coupling_name,"surf_cpl");
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
                      output_format,                            // Postprocessing format
                      output_format_option);

cwipi_dump_application_properties();

  /* Mesh definition */

  if (rank == 0) printf("        Create mesh\n");

  int nVertex = 0;               // Number of vertex
  double *coords = NULL;         // Vertex coordinates
  int nElts = 0;                 // Number of elements
  int *eltsConnecPointer = NULL; // Connectivity index
  int *eltsConnec = NULL;        // Connectivity
  
  /* Domain bounds */

  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 1;
  const double zmin = 0;
  const double zmax = 0.1;

  nVertex = (nElemSeg + 1) * (nElemSeg + 1) * 2;   // 2 vertices (1 element) in the Z direction
  nElts = nElemSeg * nElemSeg;

  coords = (double *) malloc(sizeof(double) * 3 * nElts ); // 3 coordinates (X, Y, Z) for each element
  eltsConnecPointer = (int *) malloc(sizeof(int) * (nElts + 1)); // always equal to the number of elements + 1
  eltsConnec = (int *) malloc(sizeof(int) * 8 * nElts); // 8 because as we are in 3D the elements are hexahedra

// coordinates table
int l = 0;
  for (int j = 0; j < nElemSeg; ++j){
      for (int i = 0; i < nElemSeg; ++i){
        coords[l] = (xmax - xmin) / (nElemSeg * 2) + i * (xmax - xmin) / nElemSeg;
        coords[l + 1] = (ymax - ymin) / (nElemSeg * 2) + j * (ymax - ymin) / nElemSeg;
        coords[l + 2] = (zmax - zmin) / 2;
        l = l + 3;
    }
  }

// connectivity
eltsConnecPointer[0] = 0;
int m = 0;
  for (int i = 0; i < nElts; ++i){
    eltsConnecPointer[i + 1] = m + 8;
    m = m + 8;
  }

int n = 0;
  for (int j = 0; j < nElemSeg; ++j){
    for (int i = 0; i < nElemSeg; ++i){
      eltsConnec[n] = j * nVertexSeg + i;
      eltsConnec[n + 1] = j * nVertexSeg + i + 1;
      eltsConnec[n + 2] = j * nVertexSeg + nVertexSeg * nVertexSeg + i + 1;
      eltsConnec[n + 3] = j * nVertexSeg + nVertexSeg * nVertexSeg + i;
      eltsConnec[n + 4] = (j + 1) * nVertexSeg + i;
      eltsConnec[n + 5] = (j + 1) * nVertexSeg + i + 1;
      eltsConnec[n + 6] = (j + 1) * nVertexSeg + nVertexSeg * nVertexSeg + i + 1;
      eltsConnec[n + 7] = (j + 1) * nVertexSeg + nVertexSeg * nVertexSeg + i;
      n = n + 8;
    }
  }

    if (rank == 0) printf("        Exchange Code2 <-> Code1\n");

  double *sendValues = NULL;
  double *recvValues = NULL;
  
  sendValues = (double *) malloc(sizeof(double) * nElts);
  recvValues = (double *) malloc(sizeof(double) * nElts);

  int nNotLocatedPoints = 0;

  for (int i = 0; i < nElts; i++) {
    sendValues[i] = coords[3 * i + 1];    // For the moment, we only send the value of the "y" coordinate
  }

// Define mesh
  cwipi_define_mesh(cl_coupling_name,
                    nVertex,
                    nElts,
                    coords,
                    eltsConnecPointer,
                    eltsConnec);

// Send receive
  stride = 1;
  it = itdeb;
  ttime = itdeb;
  sprintf(cl_exchange_name, "echange1");
  sprintf(cl_varsend, "coord y");
  sprintf(cl_varrec, "coord x");

  cwipi_exchange_status_t status = cwipi_exchange(cl_coupling_name,
                                                  cl_exchange_name,
                                                  stride,
                                                  it,     // n_step
                                                  ttime,   // physical_time
                                                  cl_varsend,
                                                  sendValues,
                                                  cl_varrec,
                                                  recvValues,
                                                  &nNotLocatedPoints);

     /* Coupling deletion
   * ----------------- */

  if (rank == 0) printf("        Delete coupling\n");

  // Delete coupling
  cwipi_delete_coupling(cl_coupling_name);

    /* Freeing memory
   * -------------- */
  
  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec);
  free(sendValues);
  free(recvValues);

  // Finalize
  cwipi_finalize();

  MPI_Finalize();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Info<< "\nStarting time loop\n" << endl;

    // while (runTime.loop())
    // {
    //     Info<< "Time = " << runTime.timeName() << nl << endl;

    //     #include "CourantNo.H"

    //     // Momentum predictor

    //     fvVectorMatrix UEqn
    //     (
    //         fvm::ddt(U)
    //       + fvm::div(phi, U)
    //       - fvm::laplacian(nu, U)
    //     );

    //     if (piso.momentumPredictor())
    //     {
    //         solve(UEqn == -fvc::grad(p));
    //     }

    //     // --- PISO loop
    //     while (piso.correct())
    //     {
    //         volScalarField rAU(1.0/UEqn.A());
    //         volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    //         surfaceScalarField phiHbyA
    //         (
    //             "phiHbyA",
    //             fvc::flux(HbyA)
    //           + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
    //         );

    //         adjustPhi(phiHbyA, U, p);

    //         // Update the pressure BCs to ensure flux consistency
    //         constrainPressure(p, U, phiHbyA, rAU);

    //         // Non-orthogonal pressure corrector loop
    //         while (piso.correctNonOrthogonal())
    //         {
    //             // Pressure corrector

    //             fvScalarMatrix pEqn
    //             (
    //                 fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    //             );

    //             pEqn.setReference(pRefCell, pRefValue);

    //             pEqn.solve();

    //             if (piso.finalNonOrthogonalIter())
    //             {
    //                 phi = phiHbyA - pEqn.flux();
    //             }
    //         }

    //         #include "continuityErrs.H"

    //         U = HbyA - rAU*fvc::grad(p);
    //         U.correctBoundaryConditions();
    //     }

    //     runTime.write();

    //     Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    //         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    //         << nl << endl;
    // }

    // Info<< "End\n" << endl;

    return EXIT_SUCCESS;
}


// ************************************************************************* //
