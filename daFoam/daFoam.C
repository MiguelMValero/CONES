/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
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
    daFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "initContinuityErrs.H"
    #include "createDA.H"
    #include "error.H"
    


    pisoControl piso(mesh);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        while (piso.correct())
        {
            while (piso.correctNonOrthogonal())
            {
                fvVectorMatrix Ueqn
                (
                    fvm::ddt(U)
                    + fvm::div(phi, U)
                    - fvm::laplacian(nu, U)
                );
                Ueqn.solve();
            }
        }
        #include "continuityErrs.H"
        runTime.write();
        if(rankObs >= 0){
            List<vector> Usim(obsCount);
            int n = 0;
            forAll(cellObs, x){
                //Pout << "CELL " << cellObs[x] << endl;
                Usim[n] = U[cellObs[x]];
                Pout << "CELL " << cellObs[x] << " U sim = " << Usim[n] << " U obs = " << Uobs[x] << endl;
                n++;
            }
       }
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
