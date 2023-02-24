#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "cwipiPstream.H"
#include <cwipi.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstdio>

#include <random>
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

void cwipiSendParamsChannel_IBM(const fvMesh& mesh, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath)
{
    //========== Send the Fourier coefficients to be optimized =============

    double* ParamsToSend = new double[cwipiParams];
	static int sendTag_params = 3;

	//========== Path of the correct OF instance ==========

	int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
	std::string simOF = globalRootPath + "/";
	std::ifstream file;
	char ensemble[50];
	std::string varOpt;
	sprintf(ensemble, "processor%i/UpdatedVariables", myGlobalRank-1);
	varOpt += simOF;    // Append the first character
    varOpt += ensemble; // Append the second character
    file.open(varOpt);

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
    }
    else{
        Info<< "Couldn't open the file with the parameters to optimize" << endl;
    }

    MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    if (cwipiVerbose == 1) Info<< "After sending the parameters to KF" << endl;

    delete[] ParamsToSend;
}

void cwipiRecvParamsChannel_IBM(const fvMesh& mesh, volVectorField& U, volScalarField& viscosity, volVectorField& F, volVectorField& D, int cwipiParams, int Fourier_m, int Fourier_n, int Fourier_o, 
double lambdaX, double lambdaY, double lambdaZ, double deltaChannel, int DarcyForchheimer, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath)
{
    //========== Receive coefficients of Fourier expansion ==========
    
    double* paramsToRecv = new double[cwipiParams];
	static int recvTag_params = 4;
    const volVectorField& C = mesh.C();
    MPI_Status status4;
    MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    if (cwipiVerbose == 1) Info<< "After receiving the parameters from KF" << endl;

    //========== Here we should read the already optimized coefficients 
    //(alpha, beta, gamma, delta) from cwipiParams ==========

    double *alpha, *beta, *gamma, *delta = NULL;

    alpha = (double*) malloc(sizeof(double) * (Fourier_m + 1) * (Fourier_n + 1) * (Fourier_o + 1) * 3); // 3 for the 3 components of the force
    beta = (double*) malloc(sizeof(double) * (Fourier_m + 1) * (Fourier_n + 1) * Fourier_o * 3);
    gamma = (double*) malloc(sizeof(double) * Fourier_m * (Fourier_n + 1) * (Fourier_o + 1) * 3);
    delta = (double*) malloc(sizeof(double) * Fourier_m * (Fourier_n + 1) * Fourier_o * 3);
  
    if (DarcyForchheimer == 2){
        alpha = (double*) realloc(alpha, sizeof(double) * (Fourier_m + 1) * (Fourier_n + 1) * (Fourier_o + 1) * 6);
        beta = (double*) realloc(beta, sizeof(double) * (Fourier_m + 1) * (Fourier_n + 1) * Fourier_o * 6);
        gamma = (double*) realloc(gamma, sizeof(double) * Fourier_m * (Fourier_n + 1) * (Fourier_o + 1) * 6);
        delta = (double*) realloc(delta, sizeof(double) * Fourier_m * (Fourier_n + 1) * Fourier_o * 6);
    }

    int alpha_sep = (Fourier_m + 1)*(Fourier_n + 1)*(Fourier_o + 1);
    int beta_sep = (Fourier_m + 1)*(Fourier_n + 1)*Fourier_o;
    int gamma_sep = Fourier_m*(Fourier_n + 1)*(Fourier_o + 1);
    int delta_sep = Fourier_m*(Fourier_n + 1)* Fourier_o;

    //========== Path of the correct OF instance:
    // We rewrite the new received values in a "UpdatedVariables" file ==========

    char simOFGen[500];
    char rankChar[500];
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    sprintf(rankChar, "%i", myGlobalRank-1);
    std::string simOF = globalRootPath + "/";
    char ensemble[50];
	std::string varOpt;
	sprintf(ensemble, "processor%i/UpdatedVariables", myGlobalRank-1);
	varOpt += simOF;    // Append the first character
    varOpt += ensemble; // Append the second character

    int lChar = varOpt.length();
    char UpdatedVariables[lChar];
    strcpy(UpdatedVariables, varOpt.c_str());

    if (remove(UpdatedVariables) != 0)
         if (cwipiVerbose == 1) perror( "Error deleting file with updated model's parameters" );
    else
         if (cwipiVerbose == 1) Info<< "File successfully deleted with my processor " << rankChar << nl << endl;

    if (DarcyForchheimer == 0 || DarcyForchheimer == 1){
	    for (int i = 0; i < alpha_sep; ++i){
		    alpha[i] = paramsToRecv[i];
		    alpha[i + alpha_sep] = paramsToRecv[i + cwipiParams/3];
		    alpha[i + 2*alpha_sep] = paramsToRecv[i + 2*cwipiParams/3];
	    }
	    for (int i = 0; i < beta_sep; ++i){
		    beta[i] = paramsToRecv[i + alpha_sep];
		    beta[i + beta_sep] = paramsToRecv[i + alpha_sep + cwipiParams/3];
		    beta[i + 2*beta_sep] = paramsToRecv[i + alpha_sep + 2*cwipiParams/3];
	    }
	    for (int i = 0; i < gamma_sep; ++i){
		    gamma[i] = paramsToRecv[i + alpha_sep + beta_sep];
		    gamma[i + gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + cwipiParams/3];
		    gamma[i + 2*gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + 2*cwipiParams/3];
	    }
	    for (int i = 0; i < delta_sep; ++i){
		    delta[i] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep];
		    delta[i + delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + cwipiParams/3];
		    delta[i + 2*delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + 2*cwipiParams/3];
	    }
    }
    else if (DarcyForchheimer == 2){
	    for (int i = 0; i < alpha_sep; ++i){
		    alpha[i] = paramsToRecv[i];
		    alpha[i + alpha_sep] = paramsToRecv[i + cwipiParams/6];
		    alpha[i + 2*alpha_sep] = paramsToRecv[i + cwipiParams/3];
            alpha[i + 3*alpha_sep] = paramsToRecv[i + cwipiParams/2];
            alpha[i + 4*alpha_sep] = paramsToRecv[i + 2*cwipiParams/3];
            alpha[i + 5*alpha_sep] = paramsToRecv[i + 5*cwipiParams/6];
	    }
	    for (int i = 0; i < beta_sep; ++i){
		    beta[i] = paramsToRecv[i + alpha_sep];
		    beta[i + beta_sep] = paramsToRecv[i + alpha_sep + cwipiParams/6];
		    beta[i + 2*beta_sep] = paramsToRecv[i + alpha_sep + cwipiParams/3];
            beta[i + 3*beta_sep] = paramsToRecv[i + alpha_sep + cwipiParams/2];
            beta[i + 4*beta_sep] = paramsToRecv[i + alpha_sep + 2*cwipiParams/3];
            beta[i + 5*beta_sep] = paramsToRecv[i + alpha_sep + 5*cwipiParams/6];
	    }
	    for (int i = 0; i < gamma_sep; ++i){
		    gamma[i] = paramsToRecv[i + alpha_sep + beta_sep];
		    gamma[i + gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + cwipiParams/6];
		    gamma[i + 2*gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + cwipiParams/3];
            gamma[i + 3*gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + cwipiParams/2];
            gamma[i + 4*gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + 2*cwipiParams/3];
            gamma[i + 5*gamma_sep] = paramsToRecv[i + alpha_sep + beta_sep + 5*cwipiParams/6];
	    }
	    for (int i = 0; i < delta_sep; ++i){
		    delta[i] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep];
		    delta[i + delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + cwipiParams/6];
		    delta[i + 2*delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + cwipiParams/3];
            delta[i + 3*delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + cwipiParams/2];
            delta[i + 4*delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + 2*cwipiParams/3];
            delta[i + 5*delta_sep] = paramsToRecv[i + alpha_sep + beta_sep + gamma_sep + 5*cwipiParams/6];
	    }
    }

    std::ofstream myfile;
    myfile.open(UpdatedVariables, std::ios::out | std::ios::app);

	for (int i = 0; i < cwipiParams; ++i){
		myfile << paramsToRecv[i] << "\n";
	}
	myfile.close();

    //========== Our volVectorField F is equal to 0, with the exception of those cells in which y=0 and y=2 
    //(and their adjacent cells). We read those cells affected by it ==========

    //========== Bottom wall ==========//
    std::string bottomCells_file = globalRootPath + "/bottomCells.txt";
    std::ifstream bottomCells;
    bottomCells.open(bottomCells_file);

	// Read the size of the bottomCells.txt from the file created by topoSet
    std::string bottomCells_topoSet = globalRootPath + "/constant/polyMesh/sets/bottomCells";
    std::fstream file1(bottomCells_topoSet);
    GotoLine(file1, 19);

    int bottomCellsnumber;
    file1 >> bottomCellsnumber;

    if (cwipiVerbose == 1) Info<< "The number of cells at the bottom wall is " << bottomCellsnumber << nl << endl;

    int bottomCellID[bottomCellsnumber];
    double Xbot[bottomCellsnumber], Ybot[bottomCellsnumber], Zbot[bottomCellsnumber];

    if (bottomCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(bottomCells, line) )
        {
            std::stringstream ss(line);
            std::string bottomCellID_str;
            std::getline(ss, bottomCellID_str, ',');

            bottomCellID[count] = std::stoi(bottomCellID_str);
            Xbot[count] = C[bottomCellID[count]][0];
            Ybot[count] = C[bottomCellID[count]][1];
            Zbot[count] = C[bottomCellID[count]][2];

            ++count;
	    }
    }
    else
    {
        Info<< "Couldn't find bottomCells.txt" << nl << endl;
    }
    bottomCells.close();

    //========== TOP WALL ==========//

    std::string topCells_file = globalRootPath + "/topCells.txt";
    std::ifstream topCells;
    topCells.open(topCells_file);

    if (cwipiVerbose == 1) Info<< "The path with my file containing the cell IDs of the top wall is: " << topCells_file << nl << endl;

    // Read the size of the topCells.txt from the file created by topoSet
    std::string topCells_topoSet = globalRootPath + "/constant/polyMesh/sets/topCells";
    std::fstream file2(topCells_topoSet);
    GotoLine(file2, 19);

    int topCellsnumber;
    file2 >> topCellsnumber;

    if (cwipiVerbose == 1) Info<< "The number of cells at the top wall is " << topCellsnumber << nl << endl;

    int topCellID[topCellsnumber];
    double Xtop[topCellsnumber], Ytop[topCellsnumber], Ztop[topCellsnumber];

    if (topCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(topCells, line) )
        {
            std::stringstream ss(line);
            std::string topCellID_str;
            std::getline(ss, topCellID_str, ',');

            topCellID[count] = std::stoi(topCellID_str);
            Xtop[count] = C[topCellID[count]][0];
            Ytop[count] = C[topCellID[count]][1];
            Ztop[count] = C[topCellID[count]][2];

            ++count;
	    }
    }
    else
    {
        Info<< "Couldn't find topCells.txt" << nl << endl;
    }
    topCells.close();

    //========== We initialise the force F with a value equal to 0 everywhere ==========//
    forAll(mesh.cells(), cellID)
    { 
        F[cellID][0] = 0;
        F[cellID][1] = 0;
        F[cellID][2] = 0;
        D[cellID][0] = 0;
        D[cellID][1] = 0;
        D[cellID][2] = 0;
    }

    //========= Fourier expansion at the bottom wall ==========//
    int count;
    for (int i = 0; i < bottomCellsnumber; ++i){
    
        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double cx = std::cos(2*M_PI*j*Xbot[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*Ybot[i] / lambdaY);
                    double cz = std::cos(2*M_PI*l*Zbot[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[bottomCellID[i]][0] += alpha[count] * cx * cy * cz;
                        F[bottomCellID[i]][1] += alpha[count+alpha_sep] * cx * cy * cz;
                        F[bottomCellID[i]][2] += alpha[count+2*alpha_sep] * cx * cy * cz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[bottomCellID[i]][0] += alpha[count] * cx * cy * cz;
                        D[bottomCellID[i]][1] += alpha[count+alpha_sep] * cx * cy * cz;
                        D[bottomCellID[i]][2] += alpha[count+2*alpha_sep] * cx * cy * cz;
                        F[bottomCellID[i]][0] += -viscosity[0] * alpha[count] * cx * cy * cz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -viscosity[0] * alpha[count+alpha_sep] * cx * cy * cz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -viscosity[0] * alpha[count+2*alpha_sep] * cx * cy * cz * U[bottomCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[bottomCellID[i]][0] + U[bottomCellID[i]][1] + U[bottomCellID[i]][2]);
                        F[bottomCellID[i]][0] += -(viscosity[0] * alpha[count] + 0.5*UMag*alpha[count+3*alpha_sep]) * cx * cy * cz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -(viscosity[0] * alpha[count+alpha_sep] + 0.5*UMag*alpha[count+4*alpha_sep]) * cx * cy * cz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -(viscosity[0] * alpha[count+2*alpha_sep] + 0.5*UMag*alpha[count+5*alpha_sep]) * cx * cy * cz * U[bottomCellID[i]][2];
                    }

                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double cx = std::cos(2*M_PI*j*Xbot[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*Ybot[i] / lambdaY);
                    double sz = std::sin(2*M_PI*l*Zbot[i] / lambdaZ);

                    if (DarcyForchheimer == 0) {
                        F[bottomCellID[i]][0] += beta[count] * cx * cy * sz;
                        F[bottomCellID[i]][1] += beta[count+beta_sep] * cx * cy * sz;
                        F[bottomCellID[i]][2] += beta[count+2*beta_sep] * cx * cy * sz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[bottomCellID[i]][0] += beta[count] * cx * cy * sz;
                        D[bottomCellID[i]][1] += beta[count+beta_sep] * cx * cy * sz;
                        D[bottomCellID[i]][2] += beta[count+2*beta_sep] * cx * cy * sz;
                        F[bottomCellID[i]][0] += -viscosity[0] * beta[count] * cx * cy * sz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -viscosity[0] * beta[count+beta_sep] * cx * cy * sz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -viscosity[0] * beta[count+2*beta_sep] * cx * cy * sz * U[bottomCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[bottomCellID[i]][0] + U[bottomCellID[i]][1] + U[bottomCellID[i]][2]);
                        F[bottomCellID[i]][0] += -(viscosity[0] * beta[count] + 0.5*UMag*beta[count+3*beta_sep]) * cx * cy * sz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -(viscosity[0] * beta[count+beta_sep] + 0.5*UMag*beta[count+4*beta_sep]) * cx * cy * sz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -(viscosity[0] * beta[count+2*beta_sep] + 0.5*UMag*beta[count+5*beta_sep]) * cx * cy * sz * U[bottomCellID[i]][2];
                    }

                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double sx = std::sin(2*M_PI*j*Xbot[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*Ybot[i] / lambdaY);
                    double cz = std::cos(2*M_PI*l*Zbot[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[bottomCellID[i]][0] += gamma[count] * sx * cy * cz;
                        F[bottomCellID[i]][1] += gamma[count+gamma_sep] * sx * cy * cz;
                        F[bottomCellID[i]][2] += gamma[count+2*gamma_sep] * sx * cy * cz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[bottomCellID[i]][0] += gamma[count] * sx * cy * cz;
                        D[bottomCellID[i]][1] += gamma[count+gamma_sep] * sx * cy * cz;
                        D[bottomCellID[i]][2] += gamma[count+2*gamma_sep] * sx * cy * cz;
                        F[bottomCellID[i]][0] += -viscosity[0] * gamma[count] * sx * cy * cz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -viscosity[0] * gamma[count+gamma_sep] * sx * cy * cz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -viscosity[0] * gamma[count+2*gamma_sep] * sx * cy * cz * U[bottomCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[bottomCellID[i]][0] + U[bottomCellID[i]][1] + U[bottomCellID[i]][2]);
                        F[bottomCellID[i]][0] += -(viscosity[0] * gamma[count] + 0.5*UMag*gamma[count+3*gamma_sep]) * sx * cy * cz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -(viscosity[0] * gamma[count+gamma_sep] + 0.5*UMag*gamma[count+4*gamma_sep]) * sx * cy * cz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -(viscosity[0] * gamma[count+2*gamma_sep] + 0.5*UMag*gamma[count+5*gamma_sep]) * sx * cy * cz * U[bottomCellID[i]][2];
                    }
                                                            
                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;
            
                    double sx = std::sin(2*M_PI*j*Xbot[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*Ybot[i] / lambdaY);
                    double sz = std::sin(2*M_PI*l*Zbot[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[bottomCellID[i]][0] += delta[count] * sx * cy * sz;
                        F[bottomCellID[i]][1] += delta[count+delta_sep] * sx * cy * sz;
                        F[bottomCellID[i]][2] += delta[count+2*delta_sep] * sx * cy * sz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[bottomCellID[i]][0] += delta[count] * sx * cy * sz;
                        D[bottomCellID[i]][1] += delta[count+delta_sep] * sx * cy * sz;
                        D[bottomCellID[i]][2] += delta[count+2*delta_sep] * sx * cy * sz;
                        F[bottomCellID[i]][0] += -viscosity[0] * delta[count] * sx * cy * sz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -viscosity[0] * delta[count+delta_sep] * sx * cy * sz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -viscosity[0] * delta[count+2*delta_sep] * sx * cy * sz * U[bottomCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[bottomCellID[i]][0] + U[bottomCellID[i]][1] + U[bottomCellID[i]][2]);
                        F[bottomCellID[i]][0] += -(viscosity[0] * delta[count] + 0.5*UMag*delta[count+3*delta_sep]) * sx * cy * sz * U[bottomCellID[i]][0];
                        F[bottomCellID[i]][1] += -(viscosity[0] * delta[count+delta_sep] + 0.5*UMag*delta[count+4*delta_sep]) * sx * cy * sz * U[bottomCellID[i]][1];
                        F[bottomCellID[i]][2] += -(viscosity[0] * delta[count+2*delta_sep] + 0.5*UMag*delta[count+5*delta_sep]) * sx * cy * sz * U[bottomCellID[i]][2];
                    }

                    ++count;
                }
            }
        }
    }

    //========= Fourier expansion at the top wall ==========//                
    for (int i = 0; i < topCellsnumber; ++i){
    
        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double cx = std::cos(2*M_PI*j*Xtop[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*(2*deltaChannel-Ytop[i]) / lambdaY);
                    double cz = std::cos(2*M_PI*l*Ztop[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[topCellID[i]][0] += alpha[count] * cx * cy * cz;
                        F[topCellID[i]][1] += -alpha[count+alpha_sep] * cx * cy * cz;
                        F[topCellID[i]][2] += alpha[count+2*alpha_sep] * cx * cy * cz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[topCellID[i]][0] += alpha[count] * cx * cy * cz;
                        D[topCellID[i]][1] += alpha[count+alpha_sep] * cx * cy * cz;
                        D[topCellID[i]][2] += alpha[count+2*alpha_sep] * cx * cy * cz;
                        F[topCellID[i]][0] += -viscosity[0] * alpha[count] * cx * cy * cz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -viscosity[0] * alpha[count+alpha_sep] * cx * cy * cz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -viscosity[0] * alpha[count+2*alpha_sep] * cx * cy * cz * U[topCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[topCellID[i]][0] + U[topCellID[i]][1] + U[topCellID[i]][2]);
                        F[topCellID[i]][0] += -(viscosity[0] * alpha[count] + 0.5*UMag*alpha[count+3*alpha_sep]) * cx * cy * cz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -(viscosity[0] * alpha[count+alpha_sep] + 0.5*UMag*alpha[count+4*alpha_sep]) * cx * cy * cz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -(viscosity[0] * alpha[count+2*alpha_sep] + 0.5*UMag*alpha[count+5*alpha_sep]) * cx * cy * cz * U[topCellID[i]][2];
                    }

                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double cx = std::cos(2*M_PI*j*Xtop[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*(2*deltaChannel-Ytop[i]) / lambdaY);
                    double sz = std::sin(2*M_PI*l*Ztop[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[topCellID[i]][0] += beta[count] * cx * cy * sz;
                        F[topCellID[i]][1] += -beta[count+beta_sep] * cx * cy * sz;
                        F[topCellID[i]][2] += beta[count+2*beta_sep] * cx * cy * sz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[topCellID[i]][0] += beta[count] * cx * cy * sz;
                        D[topCellID[i]][1] += beta[count+beta_sep] * cx * cy * sz;
                        D[topCellID[i]][2] += beta[count+2*beta_sep] * cx * cy * sz;
                        F[topCellID[i]][0] += -viscosity[0] * beta[count] * cx * cy * sz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -viscosity[0] * beta[count+beta_sep] * cx * cy * sz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -viscosity[0] * beta[count+2*beta_sep] * cx * cy * sz * U[topCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[topCellID[i]][0] + U[topCellID[i]][1] + U[topCellID[i]][2]);
                        F[topCellID[i]][0] += -(viscosity[0] * beta[count] + 0.5*UMag*beta[count+3*beta_sep]) * cx * cy * sz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -(viscosity[0] * beta[count+beta_sep] + 0.5*UMag*beta[count+4*beta_sep]) * cx * cy * sz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -(viscosity[0] * beta[count+2*beta_sep] + 0.5*UMag*beta[count+5*beta_sep]) * cx * cy * sz * U[topCellID[i]][2];
                    }

                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;

                    double sx = std::sin(2*M_PI*j*Xtop[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*(2*deltaChannel-Ytop[i]) / lambdaY);
                    double cz = std::cos(2*M_PI*l*Ztop[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[topCellID[i]][0] += gamma[count] * sx * cy * cz;
                        F[topCellID[i]][1] += -gamma[count+gamma_sep] * sx * cy * cz;
                        F[topCellID[i]][2] += gamma[count+2*gamma_sep] * sx * cy * cz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[topCellID[i]][0] += gamma[count] * sx * cy * cz;
                        D[topCellID[i]][1] += gamma[count+gamma_sep] * sx * cy * cz;
                        D[topCellID[i]][2] += gamma[count+2*gamma_sep] * sx * cy * cz;
                        F[topCellID[i]][0] += -viscosity[0] * gamma[count] * sx * cy * cz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -viscosity[0] * gamma[count+gamma_sep] * sx * cy * cz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -viscosity[0] * gamma[count+2*gamma_sep] * sx * cy * cz * U[topCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[topCellID[i]][0] + U[topCellID[i]][1] + U[topCellID[i]][2]);
                        F[topCellID[i]][0] += -(viscosity[0] * gamma[count] + 0.5*UMag*gamma[count+3*gamma_sep]) * sx * cy * cz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -(viscosity[0] * gamma[count+gamma_sep] + 0.5*UMag*gamma[count+4*gamma_sep]) * sx * cy * cz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -(viscosity[0] * gamma[count+2*gamma_sep] + 0.5*UMag*gamma[count+5*gamma_sep]) * sx * cy * cz * U[topCellID[i]][2];
                    }
                                                            
                    ++count;
                }
            }
        }

        count = 0;
        for (int j = 0; j < (Fourier_m + 1); ++j){
            for (int k = 0; k < (Fourier_n + 1); ++k){
                for (int l = 0; l < (Fourier_o + 1); ++l){
                    //Info<< "count " << count << endl;
            
                    double sx = std::sin(2*M_PI*j*Xtop[i] / lambdaX);
                    double cy = std::cos(2*M_PI*k*(2*deltaChannel-Ytop[i]) / lambdaY);
                    double sz = std::sin(2*M_PI*l*Ztop[i] / lambdaZ);

                    if (DarcyForchheimer == 0){
                        F[topCellID[i]][0] += delta[count] * sx * cy * sz;
                        F[topCellID[i]][1] += -delta[count+delta_sep] * sx * cy * sz;
                        F[topCellID[i]][2] += delta[count+2*delta_sep] * sx * cy * sz;
                    }
                    else if (DarcyForchheimer == 1){
                        D[topCellID[i]][0] += delta[count] * sx * cy * sz;
                        D[topCellID[i]][1] += delta[count+delta_sep] * sx * cy * sz;
                        D[topCellID[i]][2] += delta[count+2*delta_sep] * sx * cy * sz;
                        F[topCellID[i]][0] += -viscosity[0] * delta[count] * sx * cy * sz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -viscosity[0] * delta[count+delta_sep] * sx * cy * sz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -viscosity[0] * delta[count+2*delta_sep] * sx * cy * sz * U[topCellID[i]][2];
                    }
                    else if (DarcyForchheimer == 2){
                        scalar UMag = std::sqrt(U[topCellID[i]][0] + U[topCellID[i]][1] + U[topCellID[i]][2]);
                        F[topCellID[i]][0] += -(viscosity[0] * delta[count] + 0.5*UMag*delta[count+3*delta_sep]) * sx * cy * sz * U[topCellID[i]][0];
                        F[topCellID[i]][1] += -(viscosity[0] * delta[count+delta_sep] + 0.5*UMag*delta[count+4*delta_sep]) * sx * cy * sz * U[topCellID[i]][1];
                        F[topCellID[i]][2] += -(viscosity[0] * delta[count+2*delta_sep] + 0.5*UMag*delta[count+5*delta_sep]) * sx * cy * sz * U[topCellID[i]][2];
                    }

                    ++count;
                }
            }
        }
    };

    delete[] paramsToRecv;
    free(alpha);
    free(beta);
    free(gamma);
    free(delta);
    if (cwipiVerbose == 1) Info<< "Updated Force term in Fourier expansion " << endl;
}

void cwipiRecvParamsChannelNOFourier(const fvMesh& mesh, volVectorField& U, volScalarField& viscosity, volScalarField& nu, volVectorField& F, volVectorField& D, int cwipiParams, double deltaChannel, int DarcyForchheimer, 
int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath)
{
    //========== Receive coefficients of Fourier expansion ==========
    
    double* paramsToRecv = new double[cwipiParams];
	static int recvTag_params = 4;
    const volVectorField& C = mesh.C();
    MPI_Status status4;
    MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    if (cwipiVerbose == 1) Info<< "After receiving the parameters from KF" << endl;

    //========== We are imposing the coefficients of D to be positive ==========
    const double mean = 0;
    const double stddev = 10;
    std::normal_distribution<double> dist(mean, stddev);
    float gaussample;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    for (int i = 0; i<cwipiParams; i++){
        if (paramsToRecv[i] < 0){
            do
            {
                gaussample = dist(generator);
            }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));
            paramsToRecv[i] = 10 + 10*gaussample;
        }
    }


    //========== Path of the correct OF instance:
    // We rewrite the new received values in a "UpdatedVariables" file ==========

    char simOFGen[500];
    char rankChar[500];
    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    sprintf(rankChar, "%i", myGlobalRank-1);
    std::string simOF = globalRootPath + "/";
    char ensemble[50];
	std::string varOpt;
	sprintf(ensemble, "processor%i/UpdatedVariables", myGlobalRank-1);
	varOpt += simOF;    // Append the first character
    varOpt += ensemble; // Append the second character

    int lChar = varOpt.length();
    char UpdatedVariables[lChar];
    strcpy(UpdatedVariables, varOpt.c_str());

    if (remove(UpdatedVariables) != 0)
         if (cwipiVerbose == 1) perror( "Error deleting file with updated model's parameters" );
    else
         if (cwipiVerbose == 1) Info<< "File successfully deleted with my processor " << rankChar << nl << endl;

    double *DD = NULL;
    DD = (double*) malloc(sizeof(double) * 9);
    std::ofstream myfile;
    myfile.open(UpdatedVariables, std::ios::out | std::ios::app);

	for (int i = 0; i < cwipiParams; ++i){
		myfile << paramsToRecv[i] << "\n";
        DD[i] = paramsToRecv[i];
	}
	myfile.close();

    //========== We read those cells affected by the force ==========
    // The volumeVectorField F is 0 in all cells with the exception of the cells where
    // y = 0, y = 2 and their adjacent cells

    //========== BOTTOM WALL ==========//

    std::string bottomCells_file = globalRootPath + "/bottomCells.txt";
    std::ifstream bottomCells;
    bottomCells.open(bottomCells_file);

    if (cwipiVerbose == 1) Info<< "The path with my file containing the cell IDs of the bottom wall is: " << bottomCells_file << nl << endl;

    // Read the size of the bottomCells.txt from the file created by topoSet
    std::string bottomCells_topoSet = globalRootPath + "/constant/polyMesh/sets/bottomCells";
    std::fstream file1(bottomCells_topoSet);
    GotoLine(file1, 19);
    int bottomCellsnumber;
    file1 >> bottomCellsnumber;

    if (cwipiVerbose == 1) Info<< "The number of cells at the bottom wall is " << bottomCellsnumber << nl << endl;

    int bottomCellID[bottomCellsnumber];
    double Xbot[bottomCellsnumber], Ybot[bottomCellsnumber], Zbot[bottomCellsnumber];

    if (bottomCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(bottomCells, line) )
        {
            std::stringstream ss(line);
            std::string bottomCellID_str;
            std::getline(ss, bottomCellID_str, ',');

            bottomCellID[count] = std::stoi(bottomCellID_str);
            Xbot[count] = C[bottomCellID[count]][0];
            Ybot[count] = C[bottomCellID[count]][1];
            Zbot[count] = C[bottomCellID[count]][2];

            ++count;
	    }
    }
    else
    {
        Info<< "Couldn't find bottomCells.txt" << nl << endl;
    }
    bottomCells.close();

    //========== BELOW BOTTOM WALL ==========//

    std::string belowbottomCells_file = globalRootPath + "/belowbottomCells.txt";
    std::ifstream belowbottomCells;
    belowbottomCells.open(belowbottomCells_file);

    if (cwipiVerbose == 1) Info<< "The path with my file containing the cell IDs below the bottom wall is: " << belowbottomCells_file << nl << endl;

    // Read the size of the belowbottomCells.txt from the file created by topoSet
    std::string belowbottomCells_topoSet = globalRootPath + "/constant/polyMesh/sets/belowbottomCells";
    std::fstream file2(belowbottomCells_topoSet);
    GotoLine(file2, 19);

    int belowbottomCellsnumber;
    file2 >> belowbottomCellsnumber;

    if (cwipiVerbose == 1) Info<< "The number of cells below my bottom wall is " << belowbottomCellsnumber << nl << endl;

    int belowbottomCellID[belowbottomCellsnumber];
    double Xbelowbottom[belowbottomCellsnumber], Ybelowbottom[belowbottomCellsnumber], Zbelowbottom[belowbottomCellsnumber];

    if (belowbottomCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(belowbottomCells, line) )
        {
            std::stringstream ss(line);
            std::string belowbottomCellID_str;
            std::getline(ss, belowbottomCellID_str, ',');

            belowbottomCellID[count] = std::stoi(belowbottomCellID_str);
            Xbelowbottom[count] = C[belowbottomCellID[count]][0];
            Ybelowbottom[count] = C[belowbottomCellID[count]][1];
            Zbelowbottom[count] = C[belowbottomCellID[count]][2];

            ++count;
	    }
    }
    else
    {
        Info<< "Couldn't find belowbottomCells.txt" << nl << endl;
    }   
    belowbottomCells.close();

    //========== ABOVE BOTTOM WALL ==========//

    std::string abovebottomCells_file = globalRootPath + "/abovebottomCells.txt";
    std::ifstream abovebottomCells;
    abovebottomCells.open(abovebottomCells_file);

    if (cwipiVerbose == 1) Info<< "The path with my file containing the cell IDs above the bottom wall is: " << abovebottomCells_file << nl << endl;

    // Read the size of the abovebottomCells.txt from the file created by topoSet
    std::string abovebottomCells_topoSet = globalRootPath + "/constant/polyMesh/sets/abovebottomCells";
    std::fstream file3(abovebottomCells_topoSet);
    GotoLine(file3, 19);

    int abovebottomCellsnumber;
    file3 >> abovebottomCellsnumber;

    if (cwipiVerbose == 1) Info<< "The number of cells above my bottom wall is " << abovebottomCellsnumber << nl << endl;

    int abovebottomCellID[abovebottomCellsnumber];
    double Xabovebottom[abovebottomCellsnumber], Yabovebottom[abovebottomCellsnumber], Zabovebottom[abovebottomCellsnumber];

    if (abovebottomCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(abovebottomCells, line) )
        {
            std::stringstream ss(line);
            std::string abovebottomCellID_str;
            std::getline(ss, abovebottomCellID_str, ',');

            abovebottomCellID[count] = std::stoi(abovebottomCellID_str);
            Xabovebottom[count] = C[abovebottomCellID[count]][0];
            Yabovebottom[count] = C[abovebottomCellID[count]][1];
            Zabovebottom[count] = C[abovebottomCellID[count]][2];

            ++count;
	    }
    }
    else
    {
        Info<< "Couldn't find abovebottomCells.txt" << nl << endl;
    }
    abovebottomCells.close();

    if (cwipiVerbose == 1)Info<< "Before initialization of the force" << endl;

    //========== We initialise the force F with a value equal to 0 everywhere ==========//
    forAll(mesh.cells(), cellID)
    { 
        F[cellID][0] = 0;
        F[cellID][1] = 0;
        F[cellID][2] = 0;
        D[cellID][0] = 0;
        D[cellID][1] = 0;
        D[cellID][2] = 0;
        nu[cellID] = viscosity[0];
    }

    //========= Calculation of F and D near the walls ==========//
    point PSymbottom, PSymbelowbottom, PSymabovebottom;
    for (int i = 0; i < bottomCellsnumber; ++i){
    
        D[bottomCellID[i]][0] = DD[1];
        D[bottomCellID[i]][1] = DD[4];
        D[bottomCellID[i]][2] = DD[7];

        F[bottomCellID[i]][0] = viscosity[0]*D[bottomCellID[i]][0]*U[bottomCellID[i]][0];
        F[bottomCellID[i]][1] = viscosity[0]*D[bottomCellID[i]][1]*U[bottomCellID[i]][1];
        F[bottomCellID[i]][2] = viscosity[0]*D[bottomCellID[i]][2]*U[bottomCellID[i]][2];

        D[belowbottomCellID[i]][0] = DD[0];
        D[belowbottomCellID[i]][1] = DD[3];
        D[belowbottomCellID[i]][2] = DD[6];

        F[belowbottomCellID[i]][0] = viscosity[0]*D[belowbottomCellID[i]][0]*U[belowbottomCellID[i]][0];
        F[belowbottomCellID[i]][1] = viscosity[0]*D[belowbottomCellID[i]][1]*U[belowbottomCellID[i]][1];
        F[belowbottomCellID[i]][2] = viscosity[0]*D[belowbottomCellID[i]][2]*U[belowbottomCellID[i]][2];

        D[abovebottomCellID[i]][0] = DD[2];
        D[abovebottomCellID[i]][1] = DD[5];
        D[abovebottomCellID[i]][2] = DD[8];

        F[abovebottomCellID[i]][0] = viscosity[0]*D[abovebottomCellID[i]][0]*U[abovebottomCellID[i]][0];
        F[abovebottomCellID[i]][1] = viscosity[0]*D[abovebottomCellID[i]][1]*U[abovebottomCellID[i]][1];
        F[abovebottomCellID[i]][2] = viscosity[0]*D[abovebottomCellID[i]][2]*U[abovebottomCellID[i]][2];

        //========== Find the symmetric points ==========//
        PSymbottom.x() = Xbot[i];
        PSymbottom.y() = 2*(deltaChannel - Ybot[i]) + Ybot[i];
        PSymbottom.z() = Zbot[i];
        label Symbottom = mesh.findCell(PSymbottom);

        D[Symbottom][0] = D[bottomCellID[i]][0];
        D[Symbottom][1] = D[bottomCellID[i]][1];
        D[Symbottom][2] = D[bottomCellID[i]][2];

        F[Symbottom][0] = viscosity[0]*D[Symbottom][0]*U[Symbottom][0];
        F[Symbottom][1] = viscosity[0]*D[Symbottom][1]*U[Symbottom][1];
        F[Symbottom][2] = viscosity[0]*D[Symbottom][2]*U[Symbottom][2];

        PSymbelowbottom.x() = Xbelowbottom[i];
        PSymbelowbottom.y() = 2*(deltaChannel - Ybelowbottom[i]) + Ybelowbottom[i];
        PSymbelowbottom.z() = Zbelowbottom[i];
        label Symbelowbottom = mesh.findCell(PSymbelowbottom);

        D[Symbelowbottom][0] = D[belowbottomCellID[i]][0];
        D[Symbelowbottom][1] = D[belowbottomCellID[i]][1];
        D[Symbelowbottom][2] = D[belowbottomCellID[i]][2];

        F[Symbelowbottom][0] = viscosity[0]*D[Symbelowbottom][0]*U[Symbelowbottom][0];
        F[Symbelowbottom][1] = viscosity[0]*D[Symbelowbottom][1]*U[Symbelowbottom][1];
        F[Symbelowbottom][2] = viscosity[0]*D[Symbelowbottom][2]*U[Symbelowbottom][2];

        PSymabovebottom.x() = Xabovebottom[i];
        PSymabovebottom.y() = 2*(deltaChannel - Yabovebottom[i]) + Yabovebottom[i];
        PSymabovebottom.z() = Zabovebottom[i];
        label Symabovebottom = mesh.findCell(PSymabovebottom);

        D[Symabovebottom][0] = D[abovebottomCellID[i]][0];
        D[Symabovebottom][1] = D[abovebottomCellID[i]][1];
        D[Symabovebottom][2] = D[abovebottomCellID[i]][2];

        F[Symabovebottom][0] = viscosity[0]*D[Symabovebottom][0]*U[Symabovebottom][0];
        F[Symabovebottom][1] = viscosity[0]*D[Symabovebottom][1]*U[Symabovebottom][1];
        F[Symabovebottom][2] = viscosity[0]*D[Symabovebottom][2]*U[Symabovebottom][2];          
}

    delete[] paramsToRecv;
    free(DD);

    if (cwipiVerbose == 1) Info<< "Coefficients of Darcy-Forchheimer defined" << nl << endl;
}
}