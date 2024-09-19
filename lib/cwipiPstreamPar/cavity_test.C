/*--------------------------------*- C++ -*----------------------------------*\
|                         __________  _   _____________                       |
|                        / ____/ __ \/ | / / ____/ ___/                       |
|                       / /   / / / /  |/ / __/  \__ \                        |
|                      / /___/ /_/ / /|  / /___ ___/ /                        |
|                      \____/\____/_/ |_/_____//____/                         |
|                                                                             |
| C.  : Coupling               |    C.O.N.Es. version : 2405                  |
| O.  : OpenFOAM (with)        |    OpenFOAM version  : OpenFOAM-org v9       |
| N.  : Numerical              |    Lucas Villanueva   /   Miguel M.Valero    |
| Es. : EnvironmentS           |    Tom Moussie        /   Sarp Er            |
|      ===============         |    Paolo Errante      /   Marcello Meldi     |
\*---------------------------------------------------------------------------*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**
 * @file cavity_test.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Cavity test case file.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "cwipiPstreamPar.H"
#include <cwipi.h>

namespace Foam
{

/**
 * @brief Send parameters for cavity test case.
 * 
 * @param mesh 
 * @param vf 
 * @param runTime 
 * @param cwipiIteration 
 * @param cwipiParams 
 * @param nbParts 
 * @param cwipiVerbose 
 */
void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a velocity boundary condition ===
    
    static int sendTag_params = 11; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "Before sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl; 

    label top = mesh.boundaryMesh().findPatchID("movingWall");

    int cellsCount = 0;
    forAll(vf.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    Pout << "What is inside movingWall boundary " << cellsCount << endl;

    if(cellsCount > 0)
    {
        double movingWallU = vf.boundaryField()[top][0].component(0);
        ParamsToSend[0] = movingWallU;
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After affection if find ok " << endl;

    label n=Pstream::nProcs();
    int tempFindMovingWall = 0;
    double tempParam = 0;
    if (Pstream::master())
    {
        for(int i=1; i<n; i++)
        {
            IPstream streamFromSlaves(Pstream::commsTypes::blocking, i);
            streamFromSlaves >> tempFindMovingWall;
            streamFromSlaves >> tempParam;
            if(tempFindMovingWall > 0)
            {
                ParamsToSend[0] = tempParam;
            }
        }

        //* We use a basic MPI primitive because we don't want any interpolation performed by
        //the cwipi primitive, the destination rank is always 0 because KF_coupling is supposed
        // to always be 0 when we launch a calculation (first in the command line)*
        MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    }
    else
    {
        // Create the input stream from processor i
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        stream2Main << cellsCount;
        stream2Main << ParamsToSend[0];
    }


    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
    

    delete[] ParamsToSend;
}

/**
 * @brief 
 * 
 * @param mesh 
 * @param vf 
 * @param runTime 
 * @param cwipiIteration 
 * @param cwipiParams 
 * @param nbParts 
 * @param cwipiVerbose 
 */
void cwipiSendParamsMGEnKF(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a velocity boundary condition ===
    
    static int sendTag_params_MGEnKF = 21; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "Before sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl; 

    label top = mesh.boundaryMesh().findPatchID("movingWall");

    int cellsCount = 0;
    forAll(vf.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    Pout << "What is inside movingWall boundary " << cellsCount << endl;

    if(cellsCount > 0)
    {
        double movingWallU = vf.boundaryField()[top][0].component(0);
        ParamsToSend[0] = movingWallU;
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After affection if find ok " << endl;

    label n=Pstream::nProcs();
    int tempFindMovingWall = 0;
    double tempParam = 0;
    if (Pstream::master())
    {
        for(int i=1; i<n; i++)
        {
            IPstream streamFromSlaves(Pstream::commsTypes::blocking, i);
            streamFromSlaves >> tempFindMovingWall;
            streamFromSlaves >> tempParam;
            if(tempFindMovingWall > 0)
            {
                ParamsToSend[0] = tempParam;
            }
        }

        //* We use a basic MPI primitive because we don't want any interpolation performed by
        //the cwipi primitive, the destination rank is always 0 because KF_coupling is supposed
        // to always be 0 when we launch a calculation (first in the command line)*
        MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params_MGEnKF, MPI_COMM_WORLD);
    }
    else
    {
        // Create the input stream from processor i
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        stream2Main << cellsCount;
        stream2Main << ParamsToSend[0];
    }


    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
    

    delete[] ParamsToSend;
}

/**
 * @brief Receive parameters for cavity test case.
 * 
 * @param mesh 
 * @param U 
 * @param cwipiParams 
 * @param nbParts 
 * @param cwipiVerbose 
 */
void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the boundary velocity send function ===
    
    static int recvTag_params = 13; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];
    
    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1; // Normally works for any decomposition
    
    label top = mesh.boundaryMesh().findPatchID("movingWall");

    int cellsCount = 0;
    forAll(U.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    if (cwipiVerbose) if (Pstream::master()) Pout << "What is inside movingWall boundary " << cellsCount << endl;

    label n=Pstream::nProcs();
    if (Pstream::master()){
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

        for(int i=1; i<n; i++)
        {
            // Create the input stream from processor i
            OPstream stream2Slaves(Pstream::commsTypes::blocking, i);
            stream2Slaves << paramsToRecv[0];
        }
        
        if(cellsCount > 0)
        {
            fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(movingWallU,faceI)
            {
                movingWallU[faceI]=vector(paramsToRecv[0],0,0);
            }
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }
    else{
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        // create the stream to send to the main proc
        IPstream streamFromMain(Pstream::commsTypes::blocking, 0);
        streamFromMain >> paramsToRecv[0];

        if(cellsCount > 0)
        {
            fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(movingWallU,faceI)
            {
                movingWallU[faceI]=vector(paramsToRecv[0],0,0);
            }
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}



void cwipiRecvParamsMGEnKF(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the boundary velocity send function ===
    
    static int recvTag_params_MGEnKF = 23; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status8;
    double* paramsToRecv = new double[cwipiParams];
    
    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1; // Normally works for any decomposition
    
    label top = mesh.boundaryMesh().findPatchID("movingWall");

    int cellsCount = 0;
    forAll(U.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    if (cwipiVerbose) if (Pstream::master()) Pout << "What is inside movingWall boundary " << cellsCount << endl;

    label n=Pstream::nProcs();
    if (Pstream::master()){
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params_MGEnKF, MPI_COMM_WORLD, &status8);

        for(int i=1; i<n; i++)
        {
            // Create the input stream from processor i
            OPstream stream2Slaves(Pstream::commsTypes::blocking, i);
            stream2Slaves << paramsToRecv[0];
        }
        
        if(cellsCount > 0)
        {
            fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(movingWallU,faceI)
            {
                movingWallU[faceI]=vector(paramsToRecv[0],0,0);
            }
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }
    else{
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        // create the stream to send to the main proc
        IPstream streamFromMain(Pstream::commsTypes::blocking, 0);
        streamFromMain >> paramsToRecv[0];

        if(cellsCount > 0)
        {
            fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];
            if (cwipiVerbose) Pout << movingWallU[0].component(0) << endl;
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(movingWallU,faceI)
            {
                movingWallU[faceI]=vector(paramsToRecv[0],0,0);
            }
            if (cwipiVerbose) Pout << "HH Les paramètres reçus du EnKF pour la simulation fine est = " << movingWallU[0].component(0) << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}


}