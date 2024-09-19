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
 * @file airfoil_naca0012.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief NACA0012 airfoil test case file.
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
 * @brief Send parameters for NACA0012 airfoil test case.
 * 
 * @param mesh OpenFOAM mesh.
 * @param vf OpenFOAM volume vector field.
 * @param runTime OpenFOAM runTime.
 * @param cwipiIteration 
 * @param cwipiParams Number of parameters.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiSendParamsAirfoil(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    static int sendTag_params = 11;
    double* ParamsToSend = new double[cwipiParams];

    int nMyProc = Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "Before sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl; 

    label inlet = mesh.boundaryMesh().findPatchID("inlet");

    const fvPatchVectorField& inletPatchField = vf.boundaryField()[inlet];
    // const dictionary& bcDict = inletPatchField.dict();

    OStringStream BCdict;
    BCdict << inletPatchField;
    string stringDict = BCdict.str();
    IStringStream ISSdict(stringDict);
    dictionary boundaryDict(ISSdict);

    double angleOfAttack = boundaryDict.lookup<scalar>("angleOfAttack");

    int cellsCount = 0;
    forAll(vf.boundaryField()[inlet],cells)
    {
        cellsCount = cellsCount + 1;
    }
    Foam::Pout << "What is inside inlet boundary " << cellsCount << endl;

    if(cellsCount > 0)
    {
        // double uInlet = vf.boundaryField()[inlet][0][].component(0);
        ParamsToSend[0] = angleOfAttack;
        Foam::Info << "THE SENT PARAMETER IN cwipiPstream is " << ParamsToSend[0] << Foam::endl;
    }

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "After affection if find ok " << endl;

    label n=Pstream::nProcs();
    int tempFindAngleOfAttack = 0;
    double tempParam = 0;
    if (Pstream::master())
    {
        for(int i=1; i<n; i++)
        {
            IPstream streamFromSlaves(Pstream::commsTypes::blocking, i);
            streamFromSlaves >> tempFindAngleOfAttack;
            streamFromSlaves >> tempParam;
            if(tempFindAngleOfAttack > 0)
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

    Foam::Info << "THE SENT PARAMETER IN cwipiPstream is " << ParamsToSend[0] << Foam::endl;
    
    if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "After sending params to KF from member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
    
    delete[] ParamsToSend;
}

/**
 * @brief Receive parameters for NACA0012 airfoil.
 * 
 * @param mesh OpenFOAM mesh.
 * @param U Velocity field.
 * @param cwipiParams Number of parameters.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiRecvParamsAirfoil(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the boundary velocity send function ===
    
    static int recvTag_params = 13; /*!< Tag for receive parameters from EnKF to OpenFOAM */
    MPI_Status status4;
    double* paramsToRecv = new double[cwipiParams];
    double pi = Foam::constant::mathematical::pi;
    
    label nMyProc=Pstream::myProcNo();
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1; // Normally works for any decomposition
    
    label inlet = mesh.boundaryMesh().findPatchID("inlet");
    const fvPatchVectorField& inletPatchField = U.boundaryField()[inlet];

    OStringStream BCdict;
    BCdict << inletPatchField;
    string stringDict = BCdict.str();
    IStringStream ISSdict(stringDict);
    dictionary boundaryDict(ISSdict);

    float magU = boundaryDict.lookup<scalar>("magU");

    int cellsCount = 0;
    forAll(U.boundaryField()[inlet],cells)
    {
        cellsCount = cellsCount + 1;
    }
    if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "What is inside inlet boundary " << cellsCount << Foam::endl;

    label n=Pstream::nProcs();
    if (Pstream::master()){
        if (cwipiVerbose) Foam::Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

        for(int i=1; i<n; i++)
        {
            // Create the input stream from processor i
            OPstream stream2Slaves(Pstream::commsTypes::blocking, i);
            stream2Slaves << paramsToRecv[0];
        }
        
        if(cellsCount > 0)
        {
            
            // fvPatchScalarField& angleOfAttack = U.boundaryFieldRef()[inlet];
            fvPatchVectorField& inletPatchFieldRef = U.boundaryFieldRef()[inlet]; //[0];
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(inletPatchFieldRef,faceI)
            {
                float newUx = magU * cos(paramsToRecv[0] * ((pi)/180));
                float newUy = magU * sin(paramsToRecv[0] * ((pi)/180));

                inletPatchFieldRef[faceI]=vector(newUx,newUy,0);
            }
            if (cwipiVerbose) Foam::Pout << paramsToRecv[0] << Foam::endl;
            
        }

        if (cwipiVerbose) Foam::Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }
    else{
        if (cwipiVerbose) Foam::Pout << "Before receive back params in member " << appSuffix  << " from processor " << Foam::Pstream::myProcNo() << Foam::endl;
        // create the stream to send to the main proc
        IPstream streamFromMain(Pstream::commsTypes::blocking, 0);
        streamFromMain >> paramsToRecv[0];

        if(cellsCount > 0)
        {
            fvPatchVectorField& inletPatchFieldRef = U.boundaryFieldRef()[inlet];
            if (cwipiVerbose) Foam::Pout << inletPatchFieldRef << Foam::endl;
            //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
            // folders indicate a uniform boundary condition automatically *
            forAll(inletPatchFieldRef,faceI)
            {
                float newUx = magU * cos(paramsToRecv[0] * ((pi)/180));
                float newUy = magU * sin(paramsToRecv[0] * ((pi)/180));

                inletPatchFieldRef[faceI]=vector(newUx,newUy,0);
            }
            // if (cwipiVerbose) Foam::Pout << angleOfAttack[1]<< Foam::endl;
        }

        if (cwipiVerbose) Foam::Pout << "After receive back params in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << Foam::endl << "\n";
    }

    delete[] paramsToRecv;
}

}