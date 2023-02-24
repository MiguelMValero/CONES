/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "cwipiPstream.H"
#include "fvCFD.H"
#include "pisoControl.H"
#include "interpolationCellPointWallModified.H"
#include <iostream>
#include <fstream>
#include <algorithm>

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "cwipiCreateTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //========== Declaration of cwipi variables ==========
    #include "cwipiVariables.H"

    //========== Create cwipi coupling and control parameters ==========
    if (cwipiSwitch)
    {
        if (cwipiVerbose == 1) Info << "Here we are" << nl << endl;
        addControlParams(numberCwipiPhase, runTime.deltaTValue(), runTime.value(), nbParts, partsRepart[1]);
        cwipiCoupling(mesh, pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, c2fconnec_size, fconnec_size, cwipiVerbose, geom_tol);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

        //========= Sending Velocity Field and Parameters and creating sampled velocities ==========
        if (cwipiSwitch && cwipiTimestep == cwipiStep)
        {
            if (cwipiParamsObs == 0) UInterpolation(U, mesh, cwipiObsU, cwipiVerbose, globalRootPath);
            else if (cwipiParamsObs == 1) pInterpolation(p, mesh, cwipiVerbose, globalRootPath);
            else if (cwipiParamsObs == 2) UpInterpolation(U, p, mesh, cwipiObsU, cwipiObsp, cwipiVerbose, globalRootPath);
            
            cwipiSend(mesh, U, runTime, cwipiIteration, cwipiVerbose);
            cwipiSendParams(mesh, U, runTime, cwipiIteration, cwipiParams, nbParts, partsRepart[1], cwipiVerbose);

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecv(mesh, U, runTime, cwipiIteration, cwipiVerbose);
            cwipiRecvParams(mesh, U, cwipiParams, nbParts, partsRepart[1], cwipiVerbose);

            // ========== We correct the pressure after the DA cycle 
            //(solve a Poisson equation for the approximate pressure taking into account the
            //updated source term)==========
            
            Info<< "Out of the receive parameters function" << endl;
            volScalarField magSqrU_DA(magSqr(U));
            volSymmTensorField FF(sqr(U)/(magSqrU_DA + small*average(magSqrU_DA)));
            volScalarField divDivUU_DA
            (
                fvc::div
                (
                    FF & fvc::div(phi, U),
                    "div(div(phi,U))"
                )
            );
            fvScalarMatrix pEqn_DA
            (
                fvm::laplacian(p) + divDivUU_DA
            );
            pEqn_DA.setReference(pRefCell, pRefValue);
            pEqn_DA.solve();

            if (cwipiVerbose == 1) Info<< "Pressure updated after DA analysis" << endl;

            cwipiPhaseCheck = 0;
            cwipiIteration = cwipiIteration + 1;
        }
        //=========================================================

        runTime.write();

    }

    //========== Delete Cwipi Coupling and allocated arrays ===========
    if (cwipiSwitch)
    {
        cwipideleteCoupling(pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, cwipiVerbose);
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
