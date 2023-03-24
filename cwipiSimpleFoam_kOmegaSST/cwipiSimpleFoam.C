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

#include "cwipiPstream.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include <iostream>
#include <fstream>
#include <cstring>

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"

    #include "cwipiCreateTime.H"
    //#include "createTime.H"

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
        if (cwipiVerbose == 1) Info << "Here we are" << nl << endl;
        addControlParams(numberCwipiPhase, runTime.deltaTValue(), runTime.value(), nbParts, partsRepart[1]);
        cwipiCoupling(mesh, pointCoords, face_index, face_connectivity_index, cell_to_face_connectivity, face_connectivity, c2fconnec_size, fconnec_size, cwipiVerbose, geom_tol);
    }
    // Info << "After if" << nl << endl;

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
            if (cwipiParamsObs == 0) UInterpolation(U, mesh, cwipiObsU, cwipiVerbose, globalRootPath);
            else if (cwipiParamsObs == 1) pInterpolation(p, mesh, cwipiVerbose, globalRootPath);
            else if (cwipiParamsObs == 2) UpInterpolation(U, p, mesh, cwipiObsU, cwipiObsp, cwipiVerbose, globalRootPath);

            cwipiSend(mesh, U, runTime, cwipiIteration, cwipiVerbose);
            cwipiSendParamsKEps(mesh, turbulence(), runTime, cwipiIteration, cwipiParams, nbParts, partsRepart[1], cwipiVerbose);
            //cwipiSendParamsKOmegaSST(mesh, turbulence(), runTime, cwipiIteration, cwipiParams, nbParts, partsRepart[1], cwipiVerbose);

            cwipiTimestep = 0;
            cwipiPhaseCheck = 1;
        }

        cwipiTimestep = cwipiTimestep + 1;


        //========= Receiving back updated Velocity Field and parameters ==========
        if (cwipiSwitch && cwipiPhaseCheck == 1)
        {
            cwipiRecv(mesh, U, runTime, cwipiIteration, cwipiVerbose);
            cwipiRecvParamsKEps(mesh, turbulence(), cwipiParams, nbParts, partsRepart[1], cwipiVerbose, globalRootPath);
            //cwipiRecvParamsKOmegaSST(mesh, turbulence(), cwipiParams, nbParts, partsRepart[1], cwipiVerbose, globalRootPath);

            // We correct the pressure after the DA cycle

            // MRF.correctBoundaryVelocity(U);
            // tmp<fvVectorMatrix> tU_DAEqn
            // (
            //     fvm::div(phi, U)
            //   + MRF.DDt(U)
            //   + turbulence->divDevSigma(U)
            //     ==
            //     fvOptions(U)
            // );
            // fvVectorMatrix& U_DAEqn = tU_DAEqn.ref();
            // U_DAEqn.relax();
            // fvOptions.constrain(U_DAEqn);

            // solve(U_DAEqn == -fvc::grad(p));

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
