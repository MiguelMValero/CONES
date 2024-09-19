/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/
/**
 * @dir ./solvers/HLEnKF/cwipiHLSimpleFoamkEpsPar
 * @brief Steady-state solver for incompressible, turbulent flows. Used for k-epsilon turbulence model in CONES.
 * 
 */
/**
 * @file cwipiHLSimpleFoamkEpsPar.C
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
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "simpleControl.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include <iostream>
#include <fstream>
#include <cstring>

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"

    #include "createTime.H"

    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //========== Declaration of cwipi variables ==========
    #include "cwipiVariables.H"

    //========== Create cwipi coupling and control parameters ==========
    if (cwipiSwitch)
    {
        addControlParams(numberCwipiPhase, runTime.deltaTValue(), runTime.value());
        cwipiCoupling(mesh, pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, c2fconnec_size, fconnec_size, subdomains, cwipiVerbose, geometricTolerance);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        turbulence->read(); // Re-read the model coefficients in case modified by DA phase (read only performed when the dictionary is modified)

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

        //========= Sending Velocity Field and Parameters and creating sampled velocities ==========
        if (cwipiSwitch && cwipiTimestep == cwipiStep)
        {
            if (cwipiVerbose) Foam::Pout<< "The remainder between the rank and the number of partitions by simulation " << myGlobalRank << " is " << myGlobalRank % subdomains << nl << endl;
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
            cwipiSendParamsKEps(mesh, turbulence(), runTime, cwipiIteration, cwipiParams, subdomains, cwipiVerbose);

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;

        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecvVolVectorField(mesh, U, runTime, cwipiIteration, subdomains, cwipiVerbose);
            cwipiRecvParamsKEps(mesh, turbulence(), cwipiParams, subdomains, cwipiVerbose, globalRootPath);
            
            // ========== We correct the pressure after the DA cycle 
            //(solve a Poisson equation for the approximate pressure taking into account the
            //updated source term)==========
            
            if (cwipiVerbose) Info<< "Out of the receive parameters function" << endl;
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

            if (cwipiVerbose) Info<< "Pressure updated after DA analysis" << endl;

            cwipiPhaseCheck = 0;
            cwipiIteration = cwipiIteration + 1;
        }
        //=========================================================

        runTime.write();
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
