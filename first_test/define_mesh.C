#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <cwipi.h>
// #include "KF_coupling_functions.h"

void define_mesh(double* pointCoords,
                int* connecIdx,
                int* connec,
				const fvMesh& mesh)
{
    // //char indexChar[50];
    // //sprintf(indexChar, "%i", index+1);

    // //char folder[50] = {"runCase_"};
    // char casePath[250] = {"/home/miguel/CWIPI_OpenFOAM"};
    // //char RunCase[250] = "";

    // char RunCase[250] = {"first_test"};

    // //strcat(RunCase,folder);
    // //strcat(RunCase,indexChar);

    // //-- Declaration des variables openFoam --
    // Foam::Time runTime(Foam::Time::controlDictName, casePath, RunCase);
    
    // Foam::fvMesh mesh
    // (
    //     Foam::IOobject
    //     (
    //         Foam::fvMesh::defaultRegion,
    //         runTime.timeName(),
    //         runTime,
    //         Foam::IOobject::MUST_READ
    //     )
    // );

    scalar nCells = mesh.nCells();
    scalar nPoints = mesh.nPoints();
    
    Info << "CC 1" << nl << endl;

    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }
    
    Info << "CC 2" << nl << endl;

    connecIdx[0]=0;
    forAll(mesh.cells(),i)
    {
        connecIdx[i+1]=connecIdx[i]+8;
    }

    Info << "CC 3" << nl << endl;

    forAll(mesh.cells(),i)
    {
        forAll(mesh.cellShapes()[i],j)
        {
            connec[8*i+j]=mesh.cellShapes()[i][j]+1;
        }
    }

    Info << "CC 4" << nl << endl;
    
    cwipi_dump_application_properties();
    cwipi_define_mesh("cwipiFoamCoupling", nPoints, nCells, pointCoords,connecIdx,connec);
    cwipi_locate("cwipiFoamCoupling");
}