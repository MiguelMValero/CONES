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
    cwipiPimplePenFoamCyl

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

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

#include "fvMesh.H"
#include <mpi.h>
// C++ headers
#include <new>


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
        c2fconnec_size, fconnec_size, nbParts, cwipiVerbose, geom_tol);
    }
    // Info << "After if" << nl << endl;

    //========== Initial force F calculation ==========//
    scalar penaltyCoeff = readScalar(runTime.controlDict().lookup("coeffD"));
    scalar diameter = readScalar(runTime.controlDict().lookup("Diameter"));
    scalar zCoord = readScalar(runTime.controlDict().lookup("zCoord"));

    //======== Some parameters for the lift and drag coefficients =========//
    scalar Uref = readScalar(runTime.controlDict().lookup("magUInf"));
    scalar rhoref = readScalar(runTime.controlDict().lookup("rhoInf"));

    #include "cylinderCells.H"
    #include "force_DarcyForchheimer.H"

    scalar drag = 0;
    scalar lift = 0;
    Info<< "The number of cells affected by the IBM force is " << forces.size() << "\n";
    for (size_t i = 0; i < forces.size(); ++i)
    {
        F[forces[i]][0] = -D[forces[i]][0]*U[forces[i]][0];
        F[forces[i]][1] = -D[forces[i]][1]*U[forces[i]][1];
        F[forces[i]][2] = 0;

        drag = drag + F[forces[i]][0]*rhoref*V[forces[i]];
        lift = lift + F[forces[i]][1]*rhoref*V[forces[i]];
        if (i == 0) Info<< "The normalized volume of the cell is " << V[forces[i]]/zCoord << endl;
    }
    #include "writeOutputs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
            if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The remainder between the rank and the number of partitions by simulation " << myGlobalRank << " is " << myGlobalRank % nbParts << nl << endl;
            if (cwipiParamsObs == 0){
                UInterpolation(U, mesh, runTime, cwipiObsU, nbParts, cwipiVerbose, globalRootPath, globalCasePath);
            }
            else if (cwipiParamsObs == 1){
                pInterpolation(p, mesh, runTime, cwipiObsp, nbParts, cwipiVerbose, globalRootPath, globalCasePath);
            }
            else if (cwipiParamsObs == 2){
                UpInterpolation(U, p, mesh, runTime, cwipiObsU, cwipiObsp, nbParts, cwipiVerbose, globalRootPath, globalCasePath);
            }
            
            cwipiSend(mesh, U, runTime, cwipiIteration, nbParts, cwipiVerbose);
            //cwipiSendParamsChannel(mesh, dummy_fld, runTime, cwipiIteration, cwipiParams, nbParts, cwipiVerbose);

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecv(mesh, U, runTime, cwipiIteration, nbParts, cwipiVerbose);
            //cwipiRecvParamsChannel(mesh, dummy_fld, cwipiParams, nbParts, cwipiVerbose);

        #include "force_DarcyForchheimer.H"
        drag = 0, lift = 0;
        for (size_t i = 0; i < forces.size(); ++i)
        {
            F[forces[i]][0] = -D[forces[i]][0]*U[forces[i]][0];
            F[forces[i]][1] = -D[forces[i]][1]*U[forces[i]][1];
            F[forces[i]][2] = 0;

            drag = drag + F[forces[i]][0]*rhoref*V[forces[i]];
            lift = lift + F[forces[i]][1]*rhoref*V[forces[i]];
	}
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

//  Commented 19 Feb. 2024 to update pressure after EnKF (Sarp)
//            while (pimple.correctNonOrthogonal()) 
//            {
                fvScalarMatrix pEqn_DA
                (
                    fvm::laplacian(p) + divDivUU_DA - fvc::div(F)  // added div(F) forcing due to IBM (IBM)
                );
                pEqn_DA.setReference(pressureReference.refCell(),pressureReference.refValue());
                pEqn_DA.solve();
//            }

        }

        runTime.write();
        #include "writeOutputs.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    //========== Delete Cwipi Coupling and allocated arrays ===========
    if (cwipiSwitch)
    {
        cwipideleteCoupling(pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, nbParts, cwipiVerbose);
    }

    delete[] cylinderCellID;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
