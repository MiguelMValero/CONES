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
 * @file pitzDaily_test.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Pitz daily test case file.
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

void cwipiSendParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameters the optimize, in this case model coefficients ===

    static int sendTag_params = 11; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (Pstream::master())
    {

        if (cwipiVerbose) std::cout << "Before sending params to KF" << std::endl; 

        //* We can read the coefficient of a turbulence model from the object turbulence. Be careful,
        // here "turbulence" is the object "turbulence()" and not the pointer "turbulence->" , because we 
        // gave "turbulence()" to the function in the solver code *
        ParamsToSend[0] = readScalar(turbulence.coeffDict().lookup("Cmu"));
        ParamsToSend[1] = readScalar(turbulence.coeffDict().lookup("C1"));
        ParamsToSend[2] = readScalar(turbulence.coeffDict().lookup("C2"));
        ParamsToSend[3] = readScalar(turbulence.coeffDict().lookup("sigmak"));
        ParamsToSend[4] = readScalar(turbulence.coeffDict().lookup("sigmaEps"));

        MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;

    delete[] ParamsToSend;
}

void cwipiSendParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameters the optimize, in this case model coefficients ===

    static int sendTag_params = 11; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (Pstream::master())
    {

        if (cwipiVerbose) std::cout << "Before sending params to KF" << std::endl; 

        //* We can read the coefficient of a turbulence model from the object turbulence. Be careful,
        // here "turbulence" is the object "turbulence()" and not the pointer "turbulence->" , because we 
        // gave "turbulence()" to the function in the solver code *
        ParamsToSend[0] = readScalar(turbulence.coeffDict().lookup("alphaK1"));
        ParamsToSend[1] = readScalar(turbulence.coeffDict().lookup("alphaK2"));
        ParamsToSend[2] = readScalar(turbulence.coeffDict().lookup("alphaOmega1"));
        ParamsToSend[3] = readScalar(turbulence.coeffDict().lookup("alphaOmega2"));
        ParamsToSend[4] = readScalar(turbulence.coeffDict().lookup("gamma1"));
        ParamsToSend[5] = readScalar(turbulence.coeffDict().lookup("gamma2"));
        ParamsToSend[6] = readScalar(turbulence.coeffDict().lookup("beta1"));
        ParamsToSend[7] = readScalar(turbulence.coeffDict().lookup("beta2"));
        ParamsToSend[8] = readScalar(turbulence.coeffDict().lookup("betaStar"));

        MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;

    delete[] ParamsToSend;
}

void cwipiRecvParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath)
{
    //=== Receive equivalent of the model coefficients send function, here we need to overwrite the momentumTransport dictionary 
    // And then read again the coefficient in the solver (turbulence->read()). We need to declare a OFstream containing the stream
    // of the file to overwrite. The stream is then closed at the end of this function. ===
    
    static int recvTag_params = 13; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    if (cwipiVerbose) if (Pstream::master()) Pout << "nMyProc is: " << nMyProc << endl;
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    if (cwipiVerbose) if (Pstream::master()) Pout << "myGlobalRank is: " << myGlobalRank << endl;
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1; // Normally works for any decomposition

    if (Pstream::master()){
        std::string dictPath = globalRootPath+"/member"+std::to_string(appSuffix)+"/constant/momentumTransport";

        if (cwipiVerbose) Pout<< "dictPath = " << dictPath << endl; 

        Foam::OFstream FileStream(dictPath);
    
        if (cwipiVerbose) Pout << "Before Re-receive params" << endl;

        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);
        const double mean = 0.0;
        const double stddev = 0.01;
        std::normal_distribution<double> dist(mean, stddev);
        float gaussample;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);

        for (int i=0; i<cwipiParams; i++)
        {
            if (paramsToRecv[i] < 0.05)
            {
            	scalar oldRecv = paramsToRecv[i];
                do
                {
                    gaussample=dist(generator);
                }   while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));

                paramsToRecv[i] = 0.05 + 0.1*gaussample;
                if (cwipiVerbose) Pout << "Parameter " << i << " from ensemble " << appSuffix << " replaced from " << oldRecv << " to " << paramsToRecv[i] << endl;
            }
        }

        if (cwipiVerbose) Pout << turbulence.subDict("RAS") << endl;

        turbulence.subDict("RAS").lookup("Cmu")[0] = paramsToRecv[0];
        turbulence.subDict("RAS").lookup("C1")[0] = paramsToRecv[1];
        turbulence.subDict("RAS").lookup("C2")[0] = paramsToRecv[2];
        turbulence.subDict("RAS").lookup("sigmak")[0] = paramsToRecv[3];
        turbulence.subDict("RAS").lookup("sigmaEps")[0] = paramsToRecv[4];

        turbulence.writeHeader(FileStream);  //First overwrite with the OF header
        turbulence.dictionary::write(FileStream, false); //Then overwrite the coefficients. "false" allows to write everything and not just the coeffs

        if (cwipiVerbose) Pout << turbulence.subDict("RAS") << endl;
    
        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}

void cwipiRecvParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath)
{
    //=== Receive equivalent of the model coefficients send function, here we need to overwrite the momentumTransport dictionary 
    // And then read again the coefficient in the solver (turbulence->read()). We need to declare a OFstream containing the stream
    // of the file to overwrite. The stream is then closed at the end of this function. ===
    
    static int recvTag_params = 13; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    if (cwipiVerbose) if (Pstream::master()) Pout << "nMyProc is: " << nMyProc << endl;
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    if (cwipiVerbose) if (Pstream::master()) Pout << "myGlobalRank is: " << myGlobalRank << endl;
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1; // Normally works for any decomposition

    if (Pstream::master()){
        std::string dictPath = globalRootPath+"/member"+std::to_string(appSuffix)+"/constant/momentumTransport";

        if (cwipiVerbose) Pout<< "dictPath = " << dictPath << endl; 

        Foam::OFstream FileStream(dictPath);
    
        if (cwipiVerbose) Pout << "Before Re-receive params" << endl;

        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);
        const double mean = 0.0;
        const double stddev = 0.11;
        std::normal_distribution<double> dist(mean, stddev);
        float gaussample;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
    
        for (int i=1; i<cwipiParams; i++)
        {
            if (paramsToRecv[i] < 0.05)
            {
                do
                {
                    gaussample=dist(generator);
                }   while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));

                paramsToRecv[i] = 0.1 + 0.1*gaussample;
            }
        }

        if (cwipiVerbose) Pout << turbulence.subDict("RAS") << endl;

        turbulence.subDict("RAS").lookup("alphaK1")[0] = paramsToRecv[0];
        turbulence.subDict("RAS").lookup("alphaK2")[0] = paramsToRecv[1];
        turbulence.subDict("RAS").lookup("alphaOmega1")[0] = paramsToRecv[2];
        turbulence.subDict("RAS").lookup("alphaOmega2")[0] = paramsToRecv[3];
        turbulence.subDict("RAS").lookup("gamma1")[0] = paramsToRecv[4];
        turbulence.subDict("RAS").lookup("gamma2")[0] = paramsToRecv[5];
        turbulence.subDict("RAS").lookup("beta1")[0] = paramsToRecv[6];
        turbulence.subDict("RAS").lookup("beta2")[0] = paramsToRecv[7];
        turbulence.subDict("RAS").lookup("betaStar")[0] = paramsToRecv[8];

        turbulence.writeHeader(FileStream);  //First overwrite with the OF header
        turbulence.dictionary::write(FileStream, false); //Then overwrite the coefficients. "false" allows to write everything and not just the coeffs

        if (cwipiVerbose) Pout << turbulence.subDict("RAS") << endl;
    
        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}
}
