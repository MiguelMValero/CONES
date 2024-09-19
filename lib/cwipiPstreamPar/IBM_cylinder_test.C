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

#include "cwipiPstreamPar.H"
#include <random>
#include <chrono>
#include <cwipi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace Foam
{

void cwipiSendParamsCylinder_IBM(const fvMesh& mesh, int cwipiParams, int nbParts, int cwipiVerbose, std::string parPath)
{
    //========== Send the expansion coefficients to be optimized =============

    double* ParamsToSend = new double[cwipiParams];
	static int sendTag_params = 11;

	//========== Path of the correct OF instance ==========
    char parGen[500];
    strcpy(parGen, parPath.c_str());
    strcat(parGen, "/UpdatedVariables.txt");   

    if (cwipiVerbose) Info<< "The path with the model's parameters is " << parGen << endl;

    std::string fileParamsSTRING;
    fileParamsSTRING += parGen; 
    std::ifstream file(fileParamsSTRING);

	int count = 0;
	if (file.is_open())
    {
        std::string line;
        while( std::getline(file,line) )
        {
            std::stringstream ss(line);
            std::string valueOpt;
            std::getline(ss, valueOpt, ' ');
            ParamsToSend[count] = std::stod(valueOpt); 
            ++count;
        }
        file.close();
    }
    else{
        Info<< "Couldn't open the file with the parameters to optimize" << endl;
    }

    MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    if (cwipiVerbose) Info<< "After sending the parameters to KF" << endl;

    delete[] ParamsToSend;
}

void cwipiRecvParamsCylinder_IBM(int cwipiParams, double* paramsToRecv, int cwipiVerbose, std::string parPath)
{
    //========= Receive coefficients from expansion =========

    static int recvTag_params = 13;
    MPI_Status status4;
    MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    if (cwipiVerbose) Info<< "After receiving the " << cwipiParams << " parameters from KF" << endl;

    //========= We rewrite the new received values in a "UpdatedVariables" file =========
	//========== Path of the correct OF instance ==========

    char parGen[500];
    strcpy(parGen, parPath.c_str());
    strcat(parGen, "/UpdatedVariables.txt"); 
    
    std::string fileParamsSTRING;
    fileParamsSTRING += parGen;
    std::ofstream file;
    file.open(fileParamsSTRING, std::ios::out);

    for (int i = 0; i < cwipiParams; ++i){
        file << paramsToRecv[i] << "\n"; 
    }
    file.close();
}

}