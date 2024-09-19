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
 * @file cavity_sin_test.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Sinusoidal cavity test case file.
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
 * @brief Send parameters for sinusoidal cavity test case.
 * 
 * @param mesh OpenFOAM mesh.
 * @param U Velocity field.
 * @param cwipiParams Number of parameters.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiSendParams_sin(const fvMesh& mesh, const volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a velocity boundary condition ===
    static int sendTag_params = 3; /*!< Tag for send parameters from OpenFOAM to EnKF */
    double* ParamsToSend = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "Before sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl; 
    label top = mesh.boundaryMesh().findPatchID("movingWall");


    int cellsCount = 0;
    forAll(U.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    if (cwipiVerbose) Pout << "What is inside BC " << cellsCount << endl;

    if(cellsCount > 0)
    {
        const fvPatchVectorField& movingWallU = U.boundaryField()[top];
        OStringStream BCdict;
        BCdict << movingWallU;
        string stringDict = BCdict.str();
        IStringStream ISSdict(stringDict);
        dictionary boundaryDict(ISSdict);
        ParamsToSend[0] = boundaryDict.lookup<vector>("amplitude1").component(0);
        ParamsToSend[1] = boundaryDict.lookup<scalar>("frequency1");
        ParamsToSend[2] = boundaryDict.lookup<scalar>("phase1");
        ParamsToSend[3] = boundaryDict.lookup<vector>("amplitude2").component(0);
        ParamsToSend[4] = boundaryDict.lookup<scalar>("frequency2");
        ParamsToSend[5] = boundaryDict.lookup<scalar>("phase2");
        ParamsToSend[6] = boundaryDict.lookup<vector>("amplitude3").component(0);
        ParamsToSend[7] = boundaryDict.lookup<scalar>("frequency3");
        ParamsToSend[8] = boundaryDict.lookup<scalar>("phase3");
        ParamsToSend[9] = boundaryDict.lookup<vector>("amplitude4").component(0);
        ParamsToSend[10] = boundaryDict.lookup<scalar>("frequency4");
        ParamsToSend[11] = boundaryDict.lookup<scalar>("phase4");
        ParamsToSend[12] = boundaryDict.lookup<vector>("offset").component(0);
    }

    if (cwipiVerbose) if (Pstream::master()) Pout << "After affection if find ok " << endl;

    label n=Pstream::nProcs();
    int tempFindBC = 0;
    double* tempParams = new double[cwipiParams];
    if (Pstream::master())
    {
        for(int i=1; i<n; i++)
        {
            IPstream streamFromSlaves(Pstream::commsTypes::blocking, i);
            streamFromSlaves >> tempFindBC;
            for (int j=0; j<cwipiParams; j++)
            {
                streamFromSlaves >> tempParams[j];
            }
            if(tempFindBC > 0)
            {
                for (int j=0; j<cwipiParams; j++)
                {
                    ParamsToSend[j] = tempParams[j];
                }
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
        for (int j=0; j<cwipiParams; j++)
        {
            stream2Main << ParamsToSend[j];
        }
    }


    if (cwipiVerbose) if (Pstream::master()) Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
    

    delete[] ParamsToSend;
    delete[] tempParams;
}

/**
 * @brief Receive parameters for sinusoidal cavity test case.
 * 
 * @param mesh OpenFOAM mesh.
 * @param U Velocity field.
 * @param cwipiParams Number of parameters.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiRecvParams_sin(fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the boundary velocity send function ===
    static int recvTag_params = 4; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];

    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    label top = mesh.boundaryMesh().findPatchID("movingWall");

    int cellsCount = 0;
    forAll(U.boundaryField()[top],cells)
    {
        cellsCount = cellsCount + 1;
    }
    if (cwipiVerbose) Pout << "What is inside BC " << cellsCount << endl;

    label n=Pstream::nProcs();
    if (Pstream::master()){
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

        for(int i=1; i<n; i++)
        {
            // Create the input stream from processor i
            OPstream stream2Slaves(Pstream::commsTypes::blocking, i);
            for (int j=0; j<cwipiParams; j++)
            {
                stream2Slaves << paramsToRecv[j];
            }
        }
        
        if(cellsCount > 0)
        {
            double t = U.db().time().timeOutputValue();
            double sin1 = paramsToRecv[0]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[1]*t + paramsToRecv[2]);
            double sin2 = paramsToRecv[3]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[4]*t + paramsToRecv[5]);
            double sin3 = paramsToRecv[6]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[7]*t + paramsToRecv[8]);
            double sin4 = paramsToRecv[9]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[10]*t + paramsToRecv[11]);
            double inletUt = sin1 + sin2 + sin3 + sin4 + paramsToRecv[12];
            string dictRecvString = "{type multiSinInletVelocity;\
            value uniform ("+ std::to_string(inletUt) + " 0 0);\
            amplitude1   ("+ std::to_string(paramsToRecv[0]) + " 0 0);\
            frequency1   "+ std::to_string(paramsToRecv[1]) + ";\
            phase1       "+ std::to_string(paramsToRecv[2]) + ";\
            amplitude2   ("+ std::to_string(paramsToRecv[3]) + " 0 0);\
            frequency2   "+ std::to_string(paramsToRecv[4]) + ";\
            phase2       "+ std::to_string(paramsToRecv[5]) + ";\
            amplitude3   ("+ std::to_string(paramsToRecv[6]) + " 0 0);\
            frequency3   "+ std::to_string(paramsToRecv[7]) + ";\
            phase3       "+ std::to_string(paramsToRecv[8]) + ";\
            amplitude4   ("+ std::to_string(paramsToRecv[9]) + " 0 0);\
            frequency4   "+ std::to_string(paramsToRecv[10]) + ";\
            phase4       "+ std::to_string(paramsToRecv[11]) + ";\
            offset       ("+ std::to_string(paramsToRecv[12]) + " 0 0);}";
            IStringStream ISSdictRecv(dictRecvString);
            // Info << "Here we are " << endl;
            dictionary boundaryDictRecv(ISSdictRecv);
            // Info << boundaryDictRecv << endl;
            U.boundaryFieldRef().set(top,fvPatchField<vector>::New(mesh.boundary()[top], U, boundaryDictRecv));
            // if (cwipiVerbose) if (Pstream::master()) Pout << "U.boundaryFieldRef()[top] : " << U.boundaryFieldRef()[top] << endl;
        }
        else{
            string dictRecvString = "{type multiSinInletVelocity;\
            value nonuniform List<vector> 0();\
            amplitude1   ("+ std::to_string(paramsToRecv[0]) + " 0 0);\
            frequency1   "+ std::to_string(paramsToRecv[1]) + ";\
            phase1       "+ std::to_string(paramsToRecv[2]) + ";\
            amplitude2   ("+ std::to_string(paramsToRecv[3]) + " 0 0);\
            frequency2   "+ std::to_string(paramsToRecv[4]) + ";\
            phase2       "+ std::to_string(paramsToRecv[5]) + ";\
            amplitude3   ("+ std::to_string(paramsToRecv[6]) + " 0 0);\
            frequency3   "+ std::to_string(paramsToRecv[7]) + ";\
            phase3       "+ std::to_string(paramsToRecv[8]) + ";\
            amplitude4   ("+ std::to_string(paramsToRecv[9]) + " 0 0);\
            frequency4   "+ std::to_string(paramsToRecv[10]) + ";\
            phase4       "+ std::to_string(paramsToRecv[11]) + ";\
            offset       ("+ std::to_string(paramsToRecv[12]) + " 0 0);}";
            IStringStream ISSdictRecv(dictRecvString);
            // Info << "Here we are " << endl;
            dictionary boundaryDictRecv(ISSdictRecv);
            // Info << boundaryDictRecv << endl;
            U.boundaryFieldRef().set(top,fvPatchField<vector>::New(mesh.boundary()[top], U, boundaryDictRecv));
            // if (cwipiVerbose) if (Pstream::master()) Pout << "U.boundaryFieldRef()[top] : " << U.boundaryFieldRef()[top] << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }
    else{
        if (cwipiVerbose) Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        // create the stream to send to the main proc
        IPstream streamFromMain(Pstream::commsTypes::blocking, 0);
        for (int j=0; j<cwipiParams; j++)
        {
            streamFromMain >> paramsToRecv[j];
        }

        if(cellsCount > 0)
        {
            double t = U.db().time().timeOutputValue();
            double sin1 = paramsToRecv[0]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[1]*t + paramsToRecv[2]);
            double sin2 = paramsToRecv[3]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[4]*t + paramsToRecv[5]);
            double sin3 = paramsToRecv[6]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[7]*t + paramsToRecv[8]);
            double sin4 = paramsToRecv[9]*Foam::sin(2*Foam::constant::mathematical::pi*paramsToRecv[10]*t + paramsToRecv[11]);
            double inletUt = sin1 + sin2 + sin3 + sin4 + paramsToRecv[12];
            string dictRecvString = "{type multiSinInletVelocity;\
            value uniform ("+ std::to_string(inletUt) + " 0 0);\
            amplitude1   ("+ std::to_string(paramsToRecv[0]) + " 0 0);\
            frequency1   "+ std::to_string(paramsToRecv[1]) + ";\
            phase1       "+ std::to_string(paramsToRecv[2]) + ";\
            amplitude2   ("+ std::to_string(paramsToRecv[3]) + " 0 0);\
            frequency2   "+ std::to_string(paramsToRecv[4]) + ";\
            phase2       "+ std::to_string(paramsToRecv[5]) + ";\
            amplitude3   ("+ std::to_string(paramsToRecv[6]) + " 0 0);\
            frequency3   "+ std::to_string(paramsToRecv[7]) + ";\
            phase3       "+ std::to_string(paramsToRecv[8]) + ";\
            amplitude4   ("+ std::to_string(paramsToRecv[9]) + " 0 0);\
            frequency4   "+ std::to_string(paramsToRecv[10]) + ";\
            phase4       "+ std::to_string(paramsToRecv[11]) + ";\
            offset       ("+ std::to_string(paramsToRecv[12]) + " 0 0);}";
            IStringStream ISSdictRecv(dictRecvString);
            // Info << "Here we are " << endl;
            dictionary boundaryDictRecv(ISSdictRecv);
            // Info << boundaryDictRecv << endl;
            U.boundaryFieldRef().set(top,fvPatchField<vector>::New(mesh.boundary()[top], U, boundaryDictRecv));
            // if (cwipiVerbose) if (Pstream::master()) Pout << "U.boundaryFieldRef()[top] : " << U.boundaryFieldRef()[top] << endl;
        }
        else{
            string dictRecvString = "{type multiSinInletVelocity;\
            value nonuniform List<vector> 0();\
            amplitude1   ("+ std::to_string(paramsToRecv[0]) + " 0 0);\
            frequency1   "+ std::to_string(paramsToRecv[1]) + ";\
            phase1       "+ std::to_string(paramsToRecv[2]) + ";\
            amplitude2   ("+ std::to_string(paramsToRecv[3]) + " 0 0);\
            frequency2   "+ std::to_string(paramsToRecv[4]) + ";\
            phase2       "+ std::to_string(paramsToRecv[5]) + ";\
            amplitude3   ("+ std::to_string(paramsToRecv[6]) + " 0 0);\
            frequency3   "+ std::to_string(paramsToRecv[7]) + ";\
            phase3       "+ std::to_string(paramsToRecv[8]) + ";\
            amplitude4   ("+ std::to_string(paramsToRecv[9]) + " 0 0);\
            frequency4   "+ std::to_string(paramsToRecv[10]) + ";\
            phase4       "+ std::to_string(paramsToRecv[11]) + ";\
            offset       ("+ std::to_string(paramsToRecv[12]) + " 0 0);}";
            IStringStream ISSdictRecv(dictRecvString);
            // Info << "Here we are " << endl;
            dictionary boundaryDictRecv(ISSdictRecv);
            // Info << boundaryDictRecv << endl;
            U.boundaryFieldRef().set(top,fvPatchField<vector>::New(mesh.boundary()[top], U, boundaryDictRecv));
            // if (cwipiVerbose) if (Pstream::master()) Pout << "U.boundaryFieldRef()[top] : " << U.boundaryFieldRef()[top] << endl;
        }

        if (cwipiVerbose) Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}

}