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
 * @file channel_test.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Plane channel test case file.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "cwipiPstreamPar.H"
#include <random>
#include <chrono>
#include <cwipi.h>

namespace Foam
{
/**
 * @brief Send parameters for plane channel test case.
 * 
 * @param mesh 
 * @param Ck 
 * @param runTime 
 * @param cwipiIteration 
 * @param cwipiParams 
 * @param nbParts 
 * @param cwipiVerbose 
 */
void cwipiSendParamsChannel(const fvMesh& mesh, const volScalarField& Ck, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a field of an uniform constant, we can imagine modifying 
    // this with a non uniform field ===

    static int sendTag_params = 11; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (Pstream::master())
    {
        if (cwipiVerbose == 1) Pout << "Before sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl; 
        // if (cwipiVerbose) Pout << "CK[0] = " << Ck[0] << endl;

        ParamsToSend[0]=Ck[0];

        //* We use a basic MPI primitive because we don't want any interpolation performed by
        //the cwipi primitive, the destination rank is always 0 because KF_coupling is supposed
        // to always be 0 when we launch a calculation (first in the command line)*
        MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;

    delete[] ParamsToSend;
}

/**
 * @brief Receive parameters for plane channel test case.
 * 
 * @param mesh 
 * @param Ck 
 * @param cwipiParams 
 * @param nbParts 
 * @param cwipiVerbose 
 */
void cwipiRecvParamsChannel(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the constant uniform field  send function ===
    
    static int recvTag_params = 13; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];
    
    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition

    label n=Pstream::nProcs();
    if (Pstream::master()){
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        if (cwipiVerbose) Info << Ck[0] << endl;

        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

        //* Here we can easaly imagine not having a uniform field *
        const double mean = 0;
        const double stddev = 0.11;
        std::normal_distribution<double> dist(mean, stddev);
        float gaussample;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);

        if (paramsToRecv[0] < 0.005) 
        {
            do
            {
                gaussample=dist(generator);
            }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));

            paramsToRecv[0] = 0.01 + 0.01*gaussample;
        }

        for(int i=1; i<n; i++)
        {
            // Create the input stream from processor i
            OPstream stream2Slaves(Pstream::commsTypes::blocking, i);
            stream2Slaves << paramsToRecv[0];
        }

        forAll(Ck,cellI)
        {
            Ck[cellI]=paramsToRecv[0];
            // Ck[cellI]=5;
        }
        

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }
    else{
        // if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        // create the stream to send to the main proc
        IPstream streamFromMain(Pstream::commsTypes::blocking, 0);
        streamFromMain >> paramsToRecv[0];

        forAll(Ck,cellI)
        {
            Ck[cellI]=paramsToRecv[0];
            // Ck[cellI]=5;
        }

        // if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}
}