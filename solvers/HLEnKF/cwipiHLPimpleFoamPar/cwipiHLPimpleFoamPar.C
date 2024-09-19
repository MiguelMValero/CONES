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
    cwipiPimpleFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/
/**
 * @dir ./solvers/HLEnKF/cwipiHLPimpleFoamPar
 * @brief Turbulent and transient simulations.
 * 
 */
/**
 * @file cwipiHLPimpleFoamPar.C
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
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
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
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
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
    // Info << "After if" << nl << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                fvModels.preUpdateMesh();

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

            fvModels.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

        //========= Sending Velocity Field and Parameters and creating sampled velocities ==========
        if (cwipiSwitch && cwipiTimestep == cwipiStep)
        {
            tic();
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
            cwipiSendParamsChannel(mesh, Ck, runTime, cwipiIteration, cwipiParams, subdomains, cwipiVerbose);

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecvVolVectorField(mesh, U, runTime, cwipiIteration, subdomains, cwipiVerbose);
            cwipiRecvParamsChannel(mesh, Ck, cwipiParams, subdomains, cwipiVerbose);

            cwipiPhaseCheck = 0;
            cwipiIteration = cwipiIteration + 1;

            //== We correct the pressure after the DA cycle ==
            // (solve a Poisson equation for the approximate pressure taking into account the updated source term)

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

            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn_DA
                (
                    fvm::laplacian(p) + divDivUU_DA
                );
                pEqn_DA.setReference(pressureReference.refCell(),pressureReference.refValue());
                pEqn_DA.solve();
            }
    
            //== We update Ck in the files and turbulence model ==
            
            Ck.write();
            turbulence->read();
            toc();
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
