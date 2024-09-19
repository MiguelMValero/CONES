/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    rhoPimpleFoam

Description
    Transient solver for turbulent flow of compressible fluids for HVAC and
    similar applications, with optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/
/**
 * @dir ./solvers/HLEnKF/cwipiHLRhoPimpleFoam
 * @brief Compressible pressure-based solver for transient simulations 
 * 
 */
/**
 * @file cwipiHLRhoPimpleFoam.C
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
#include "dynamicFvMesh.H"
#include "fluidThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "interpolationCellPointWallModified.H"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //========== Declaration of cwipi variables ==========
    #include "cwipiVariables.H"

    //========== Create cwipi coupling and control parameters ==========
    if (cwipiSwitch)
    {
        addControlParams(numberCwipiPhase, runTime.deltaTValue(), runTime.value());
        cwipiCoupling(mesh, pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, 
        c2fconnec_size, fconnec_size, subdomains, cwipiVerbose, geometricToleranceCoarse);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho*U);
                }

                fvModels.preUpdateMesh();

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if
            (
                !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
            {
                #include "rhoEqn.H"
            }

            fvModels.correct();

            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

        //========= Sending Velocity Field and Parameters and creating sampled velocities ==========
        if (cwipiSwitch && cwipiTimestep == cwipiStep)
        {
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
            // cwipiSendParams_OFRsin(mesh, U, cwipiParams, subdomains, cwipiVerbose);   //For sinusoidal inlets

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecvVolVectorField(mesh, U, runTime, cwipiIteration, subdomains, cwipiVerbose);
            
            cwipiRecvParams(mesh, U, cwipiParams, subdomains, cwipiVerbose);
            // cwipiRecvParams_OFRsin(mesh, U, cwipiParams, subdomains, cwipiVerbose);       //For sinusoidal inlets

            //== We solve the energy equation after the velocity update and update the pressure ==
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
    
            //== We update Ck in the files and turbulence model ==
            
            // Ck.write();
            // turbulence->read();

            cwipiPhaseCheck = 0;
            cwipiIteration = cwipiIteration + 1;
        }
        //=========================================================

        // #include "RSTcalculation.H"

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
