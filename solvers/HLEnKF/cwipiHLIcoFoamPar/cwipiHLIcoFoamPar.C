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
/**
 * @dir ./solvers/HLEnKF
 * @brief **HLEnKF solvers.**
 */
/**
 * @dir ./solvers/HLEnKF/cwipiHLIcoFoamPar
 * @brief Transient solver for incompressible and laminar flows.
 * 
 */
/**
 * @file cwipiHLIcoFoamPar.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Main solver file.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "cwipiPstreamPar.H"
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
    #include "createTime.H"
    //#include "cwipiCreateTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    //========== Declaration of cwipi variables ==========
    #include "cwipiVariables.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //========== Create cwipi coupling and control parameters ==========
    if (cwipiSwitch)
    {
        addControlParams(numberCwipiPhase, runTime.deltaTValue(), runTime.value());
        cwipiCoupling(mesh, pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, 
        c2fconnec_size, fconnec_size, subdomains, cwipiVerbose, geometricToleranceCoarse);
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
            if (Pstream::master()) tic();
            if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The remainder between the rank and the number of partitions by simulation " << myGlobalRank << " is " << myGlobalRank % subdomains << nl << endl;
            if (obsType == "U"){
                UInterpolation(U, mesh, runTime, cwipiObsU, subdomains, cwipiVerbose, globalRootPath, globalCasePath, obsCoordFile);
            }
            else if (obsType == "p"){
                pInterpolation(p, mesh, runTime, cwipiObsp, subdomains, cwipiVerbose, globalRootPath, globalCasePath, obsCoordFile);
            }
            else if (obsType == "Up"){
                UpInterpolation(U, p, mesh, runTime, cwipiObsU, cwipiObsp, subdomains, cwipiVerbose, globalRootPath, globalCasePath, obsCoordFile);
            }
            cwipiSendVolVectorField(mesh, U, runTime, cwipiIteration, subdomains, cwipiVerbose);

            cwipiSendParams(mesh, U, runTime, cwipiIteration, cwipiParams, subdomains, cwipiVerbose); 
            // cwipiSendParams_sin(mesh, U, cwipiParams, subdomains, cwipiVerbose);   //For sinusoidal inlets

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecvVolVectorField(mesh, U, runTime, cwipiIteration, subdomains, cwipiVerbose);
            
            cwipiRecvParams(mesh, U, cwipiParams, subdomains, cwipiVerbose);
            // cwipiRecvParams_sin(mesh, U, cwipiParams, subdomains, cwipiVerbose);       //For sinusoidal inlets

            // ========== We correct the pressure after the DA cycle 
            //(solve a Poisson equation for the approximate pressure taking into account the
            //updated source term)==========
            
            if (cwipiVerbose) Pout << "Out of the receive parameters function" << endl;
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

            if (cwipiVerbose) Pout << "Pressure updated after DA analysis" << endl;

            cwipiPhaseCheck = 0;
            cwipiIteration = cwipiIteration + 1;

            if (Pstream::master()) toc();
        }
        //=========================================================

        runTime.write();
        if (cwipiVerbose) Pout << "After write in the solver" << endl;
    }

    //========== Delete Cwipi Coupling and allocated arrays ===========
    if (cwipiSwitch)
    {
        cwipideleteCoupling(pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, subdomains, cwipiVerbose);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
