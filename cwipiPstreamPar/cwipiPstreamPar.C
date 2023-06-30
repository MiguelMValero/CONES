/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

    Functions for cwipi coupling for OpenFOAM. Some functions are specific
    to some solvers and type of parameters to optimize. Use with caution 

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "interpolationCellPointWallModified.H"
#include "cwipiPstreamPar.H"
#include <cwipi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Declaration of variables for tags and status of the exchanges
static int sendTag = 1;
static int sendTag_params = 3;
static int recvTag = 2;
static int recvTag_params = 4;
static int status;
static int status2;
MPI_Status status4;

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime)
{
    // Add control paramaters that must be sent to KF_coupling code

    cwipi_add_local_int_control_parameter("numberCwipiPhase", numberCwipiPhase); // Number of DA phases
    cwipi_add_local_double_control_parameter("deltaT", deltaT); // DeltaT of CFD loop
    cwipi_add_local_double_control_parameter("currentTime", currentTime); // Start time of CFD simulation
}

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, float cwipiVerbose, double geom_tol)
{
    //=== Function that deals with the coupling of each OF instance ===
    
    //* Calculations for declaring the coupling mesh *
    
    int nCells = mesh.nCells();
    int nFaces = mesh.nFaces();
    int nPoints = mesh.nPoints();

    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }

    face_index[0]=0;
    face_connectivity_index[0]=0;

    int *face_index_temp = new int[nCells];
    int *face_index_temp_IDs = new int[nCells];
    int *face_connectivity_index_temp = new int[nFaces];
    int *face_connectivity_index_temp_IDs = new int[nFaces];

    forAll(mesh.cells(), i)
    {
        face_index_temp[i] = mesh.cells()[i].size();
        face_index_temp_IDs[i] = i;
    }

	forAll(mesh.cells(), i)
    {
        face_index[i+1] = face_index[i] + face_index_temp[i];
    }

    int faceOwner = 0;
    int faceNeighbour = 0;
    int faceID = 0;
    int faces_count = 0;
    for(int i = 0; i < nCells; i++){
        for (int j = 0; j < face_index_temp[i] ; j++) {
            faceID = mesh.cells()[face_index_temp_IDs[i]][j];
            if (mesh.isInternalFace(faceID)){
                faceOwner = mesh.faceOwner()[faceID];
                faceNeighbour = mesh.faceNeighbour()[faceID];
                if (faceOwner == face_index_temp_IDs[i]){
                    if (faceOwner > faceNeighbour)
                    {
                        cell_to_face_connectivity[faces_count+j] = -(faceID + 1);
                    }
                    else
                    {
                        cell_to_face_connectivity[faces_count+j] = faceID + 1;
                    }
                }
                else if (faceNeighbour == face_index_temp_IDs[i]){
                    if (faceOwner > faceNeighbour)
                    {
                        cell_to_face_connectivity[faces_count+j] = faceID + 1;
                    }
                    else
                    {
                        cell_to_face_connectivity[faces_count+j] = -(faceID + 1);
                    }
                }
                else
                    break;
            }
            else
            {
                cell_to_face_connectivity[faces_count+j] = faceID + 1;
            }
        }
        faces_count = faces_count + face_index_temp[i];
    }

	//=========== Faces calculations ===========

	forAll(mesh.faces(), i)
    {
        face_connectivity_index_temp[i] = mesh.faces()[i].size();
        face_connectivity_index_temp_IDs[i] = i;
    }

    forAll(mesh.faces(), i)
    {    
        face_connectivity_index[i+1] = face_connectivity_index[i] + face_connectivity_index_temp[i];
    }

	int vertices_cout = 0;
    for(int i = 0; i < nFaces; i++)
    {
        for (int j = 0; j < face_connectivity_index_temp[i] ; j++)
        {
        face_connectivity[vertices_cout+j] = mesh.faces()[face_connectivity_index_temp_IDs[i]][j] + 1;
        }
		vertices_cout = vertices_cout + face_connectivity_index_temp[i];
    }

    //* Change the name of the coupling for each OF instance taking the rank of the proc *

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    /* Creation of the coupling and coupling mesh
    
    Options :
    CWIPI_COUPLING_SEQUENTIAL
    CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
    CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING
    */

   if (cwipiVerbose) Foam::Pout<< "The name of my coupling is " << couplingName <<endl;   

    cwipi_create_coupling(couplingName,
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          "KF",
                          3,
                          geom_tol,
                          CWIPI_STATIC_MESH,
                          CWIPI_SOLVER_CELL_CENTER,
                          0,
                          "Ensight Gold",
                          "text");

    cwipi_synchronize_control_parameter("KF");
    
    //if (cwipiVerbose == 2) cwipi_dump_application_properties();

    int t[2];
    t[0] = 0;
    t[1] = 0;

    cwipi_define_mesh(couplingName, nPoints, 0, pointCoords,t,t);
    cwipi_add_polyhedra(couplingName, 
                        nCells,
                        face_index,
                        cell_to_face_connectivity,
                        nFaces,
                        face_connectivity_index,
                        face_connectivity);
    
    printf("%s: Create mesh in OF \n", __func__);
    cwipi_locate(couplingName);

    if (cwipiVerbose) Foam::Pout << "Cwipi locate passed" << nl << endl;
    delete[] face_index_temp;
	delete[] face_index_temp_IDs;
	delete[] face_connectivity_index_temp;
	delete[] face_connectivity_index_temp_IDs;
}

void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int mainsubDomain, interpolationCellPointWallModified<vector> triangulateCellsU,
float cwipiVerbose, std::string globalPath, std::string UIntPath)
{
    //========== Produce a file for each OF instance containing the sampled velocities 
    //(H.x term in the kalman gain calculation) ========== 
    
    //========== Path of the correct OF instance ==========
    std::string proj_file = "/home/miguel/parallelizationCavity2/obs_coordinates.txt";
    if (cwipiVerbose) Pout << "The path where I have my observation coordinates is " << proj_file << endl;
    char UIntGen[500];
    strcpy(UIntGen, UIntPath.c_str());

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    strcat(UIntGen, "/UInt");

    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    char* UInt = strcat(UIntGen, appSuffixChar);

    if (cwipiVerbose) Foam::Pout<< "The file where UInt is created is " << UInt << endl;

    //* Open file to write inside the result of the velocities interpolation on the observation coordinates *
    std::ifstream file;
    file.open(proj_file);
    string data = "";
    int columns = 0;
    int parameters = 3; // 3 coordinates of the probes (only velocity field)

    if (cwipiVerbose) Foam::Pout<< "Creating the sampling files..." << endl;

    point pointCoord;
    double Ux[cwipiObsU], Uy[cwipiObsU], Uz[cwipiObsU];
    int cellID[cwipiObsU];

    //remove(UInt);

    if (file.is_open())
    {
        // printing the success message
        std::cout << "File exists" << std::endl;
    }
    else
    {
        // printing the error message
        fprintf(stderr, "File does not exist\n");
    }

    if (cwipiVerbose) Foam::Pout<< "Creating the sampling files 2..." << endl;

    if (file.is_open())
    {
        if (cwipiVerbose) std::cout<< "Inside the observation coordinates file" << std::endl;
        std::string line;
        int countProbes = 0;
        while( std::getline(file,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
            pointCoord.x() = std::stod(coordX);
            std::getline(ss, coordY,',');
            pointCoord.y() = std::stod(coordY);
            std::getline(ss, coordZ,','); 
            pointCoord.z() = std::stod(coordZ);

            label ownCell = mesh.findCell(pointCoord);
            if (ownCell == -1){
                std::cerr<< "Some sensors are not located at the main subdomain" << endl;
            }
            vector UatSP(vector::zero);
            UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

            Ux[countProbes] = UatSP.x();
            //std::cout << "The value of Ux for observation " << countProbes << " is " << Ux[countProbes] << std::endl;
            Uy[countProbes] = UatSP.y();
            //std::cout << "The value of Uy for observation " << countProbes << " is " << Uy[countProbes] << std::endl;
            Uz[countProbes] = UatSP.z();
            //std::cout << "The value of Uz for observation " << countProbes << " is " << Uz[countProbes] << std::endl;
            cellID[countProbes] = ownCell;
            //std::cout << "The value of the cellID for observation " << countProbes << " is " << cellID[countProbes] << std::endl;

            ++columns;
            if (columns == parameters){
                columns = 0;
            }
        ++countProbes;
        }
    }
    else{
        Pout<< " Not able to open the obs_coordinates.txt file " << endl;
    }
    file.close();

                    /* WE WRITE THE NEW VALUES IN A TXT FILE */
    
    std::ofstream myfile;
    myfile.open(UInt, std::ios::out | std::ios::app);

    Pout << "I am about to write the sampling file" << endl;
    if (myfile.is_open()){
        for (int i = 0; i < cwipiObsU; ++i){
            myfile << Ux[i] << "\n";
        }
        for (int i = 0; i < cwipiObsU; ++i){
            myfile << Uy[i] << "\n";
        }
        for (int i = 0; i < cwipiObsU; ++i){
            myfile << Uz[i] << "\n";
        }
    }
    else{
        std::cerr << "Not able to open the sampling file" << std::endl;
    }
    myfile.close();

    //========== Retrieve the cellIDs of my probes ==========//
    // Change "0" by initTime
    if (myGlobalRank == mainsubDomain && (runTime.value() - 0) < 1e-9){
        std::ofstream myfile_2;
        std::string cellSensors = globalPath + "/results/cellIDs";

        //remove(cellSensors.c_str());
        myfile_2.open(cellSensors, std::ios::out | std::ios::app);
        for (int i = 0; i < cwipiObsU; ++i){
            myfile_2 << cellID[i] << "\n";
        }
        myfile_2.close();
    }
}

void pInterpolation(volScalarField& p, fvMesh& mesh, int cwipiObsp, float cwipiVerbose, std::string pIntPath)
{
    //=== Produce a file for each OF instance containing the sampled pressures (H.x term in the kalman gain calculation) === 
    
    //* Path of the correct OF instance *
    std::string proj_file = pIntPath + "/obs_coordinates.txt";
    char pIntGen[500];
    strcpy(pIntGen, pIntPath.c_str());
    strcat(pIntGen, "/processor");

    if (cwipiVerbose) std::cout << "pIntPath = " << pIntPath << std::endl;

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank-1);
    strcat(pIntGen, rankChar);
    strcat(pIntGen, "/pInt");
    char* pInt = strcat(pIntGen, rankChar);

    //* Open file to write inside the result of the velocities interpolation on the observation coordinates *
    std::ifstream file;
    file.open(proj_file);
    string data = "";
    int columns = 0;
    int parameters = 1; // The pressure is a scalar

    interpolationCellPoint<double> triangulateCellsp(p);
    point pointCoord;
    int cellID[cwipiObsp];

    if (Pstream::master()){
        if (remove(pInt) != 0)
            perror( "Error deleting file with interpolated pressures" );
        else
            if (cwipiVerbose) std::cout << "Sampling file " << myGlobalRank << " successfully deleted " << std::endl;
    }
    
    std::ofstream myfile;
    myfile.open(pInt, std::ios::out | std::ios::app);

    std::ifstream inFile(proj_file);
        if (inFile.is_open())
        {
            int countProbes = 0;
            std::string line;
            while( std::getline(inFile,line) )
            {
                std::stringstream ss(line);

                std::string coordX, coordY, coordZ;
                std::getline(ss, coordX,',');
                pointCoord.x() = stod(coordX);
                std::getline(ss, coordY,',');
                pointCoord.y() = stod(coordY);
                std::getline(ss, coordZ,','); 
                pointCoord.z() = stod(coordZ);

                label ownCell = mesh.findCell(pointCoord);
                cellID[countProbes] = ownCell;
                // if (ownCell == -1){
                //     // Foam::Pout<< "Cell ID not found in my processor" << Foam::nl << Foam::endl;
                //     continue;
                // }
                scalar patSP = triangulateCellsp.interpolate(pointCoord, ownCell, -1);

                /* WE WRITE THE NEW VALUES IN A TXT FILE */

                myfile << patSP << "\n";
                ++columns;
                if (columns == parameters){
                    columns = 0;
                }

            ++countProbes;
            }
        }

    myfile.close();
    file.close();

        //========== Retrieve the cellIDs of my probes ==========//
    if (myGlobalRank == 1){
        std::ofstream myfile_2;
        std::string cellSensors = pIntPath + "/results/cellIDs";

        if (remove(cellSensors.c_str()) != 0)
            perror( "Error deleting file with cell IDs of the sensors" );
        myfile_2.open(cellSensors, std::ios::out | std::ios::app);
            for (int i = 0; i < cwipiObsp; ++i){
                myfile_2 << cellID[i] << "\n";
            }
    }
}

void UpInterpolation(volVectorField& U, volScalarField& p, fvMesh& mesh, int cwipiObsU, int cwipiObsp, float cwipiVerbose, std::string UpIntPath)
{
    //=== Produce a file for each OF instance containing the sampled velocities and pressures (H.x term in the kalman gain calculation) === 
    
    //* Path of the correct OF instance *
    std::string proj_file = UpIntPath + "/obs_coordinates.txt";
    char UpIntGen[500];
    strcpy(UpIntGen, UpIntPath.c_str());
    strcat(UpIntGen, "/processor");

    if (cwipiVerbose) std::cout << "UpintPath = " << UpIntPath << std::endl;

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank-1);
    strcat(UpIntGen, rankChar);
    strcat(UpIntGen, "/UpInt");
    char* UpInt = strcat(UpIntGen, rankChar);

    //* Open file to write inside the result of the velocities interpolation on the observation coordinates *
    std::ifstream file;
    file.open(proj_file);
    string data = "";
    int columns = 0;
    int parameters = 3; // Always coordinates here 

    interpolationCellPointWallModified<vector> triangulateCellsU(U);
    interpolationCellPoint<double> triangulateCellsp(p);
    point pointCoord;
    int cellID[cwipiObsU + cwipiObsp];
    double Ux[cwipiObsU], Uy[cwipiObsU], Uz[cwipiObsU], pInt[cwipiObsp];

    if (Pstream::master()){
        if (remove(UpInt) != 0)
            perror( "Error deleting file with interpolated values" );
        else
            if (cwipiVerbose) std::cout << "Sampling file " << myGlobalRank << " successfully deleted " << std::endl;
    }

    std::ifstream inFile(proj_file);
        if (inFile.is_open())
        {
            std::string line;
            int countProbes = 0;
            while( std::getline(inFile,line) )
            {
                std::stringstream ss(line);

                std::string coordX, coordY, coordZ;
                std::getline(ss, coordX,',');
                pointCoord.x() = stod(coordX);
                std::getline(ss, coordY,',');
                pointCoord.y() = stod(coordY);
                std::getline(ss, coordZ,','); 
                pointCoord.z() = stod(coordZ);

                label ownCell = mesh.findCell(pointCoord);
                cellID[countProbes] = ownCell;
                // if (ownCell == -1){
                //     // Foam::Pout<< "Cell ID not found in my processor" << Foam::nl << Foam::endl;
                //     continue;
                // }
                if (countProbes < cwipiObsU) {
                    vector UatSP(vector::zero);
                    UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

                    Ux[countProbes] = UatSP.x();
                    Uy[countProbes] = UatSP.y();
                    Uz[countProbes] = UatSP.z();
                }
                else {
                    scalar patSP = triangulateCellsp.interpolate(pointCoord, ownCell, -1);
                    pInt[countProbes - cwipiObsU] = patSP;
                }

                // myfile << patSP << "\n";
                ++columns;
                if (columns == parameters){
                    columns = 0;
                }
            ++countProbes;
            }
        }
    file.close();

                    /* WE WRITE THE NEW VALUES IN A TXT FILE */
    std::ofstream myfile;
    myfile.open(UpInt, std::ios::out | std::ios::app);

    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Ux[i] << "\n";
    }
    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Uy[i] << "\n";
    }
    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Uz[i] << "\n";
    }
    for (int i = 0; i < cwipiObsp; ++i){
        myfile << pInt[i] << "\n";
    }

    myfile.close();

    //========== Retrieve the cellIDs of my probes ==========//
    if (myGlobalRank == 1){
        std::ofstream myfile_2;
        std::string cellSensors = UpIntPath + "/results/cellIDs";

        if (remove(cellSensors.c_str()) != 0)
            perror( "Error deleting file with cell IDs of the sensors" );
        myfile_2.open(cellSensors, std::ios::out | std::ios::app);
            for (int i = 0; i < (cwipiObsU + cwipiObsp); ++i){
                myfile_2 << cellID[i] << "\n";
            }
    }
}

void UCfInterpolation(volVectorField& U, double deltaChannel, volScalarField& nu, fvMesh& mesh, int cwipiObsU, float cwipiVerbose, std::string UIntPath)
{
    //========== Produce a file for each OF instance containing the sampled velocities 
    //(H.x term in the kalman gain calculation) ========== 
    
    //========== Path of the correct OF instance ==========
    std::string proj_file = UIntPath + "/obs_coordinates.txt";
    char UIntGen[500];
    strcpy(UIntGen, UIntPath.c_str());
    strcat(UIntGen, "/processor");

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank-1);
    strcat(UIntGen, rankChar);
    strcat(UIntGen, "/UCfInt");
    char* UCfInt = strcat(UIntGen, rankChar);

    //* Open file to write inside the result of the velocities interpolation on the observation coordinates *
    std::ifstream file;
    file.open(proj_file);
    string data = "";
    int columns = 0;
    int parameters = 3; // 3 coordinates of the probes (only velocity field)

    if (cwipiVerbose) Info<< "Creating the sampling files..." << endl;

    interpolationCellPointWallModified<vector> triangulateCellsU(U);
    point pointCoord;
    double Ux[cwipiObsU], Uy[cwipiObsU], Uz[cwipiObsU];
    int cellID[cwipiObsU];

    if (Pstream::master()){
        if (remove(UCfInt) != 0)
            perror( "Error deleting file with interpolated velocities" );
        else
            if (cwipiVerbose) std::cout << "Sampling file " << myGlobalRank << " successfully deleted " << std::endl;
    }

    std::ifstream inFile(proj_file);
        if (inFile.is_open())
        {
            std::string line;
            int countProbes = 0;
            while( std::getline(inFile,line) )
            {
                std::stringstream ss(line);

                std::string coordX, coordY, coordZ;
                std::getline(ss, coordX,',');
                pointCoord.x() = std::stod(coordX);
                std::getline(ss, coordY,',');
                pointCoord.y() = std::stod(coordY);
                std::getline(ss, coordZ,','); 
                pointCoord.z() = std::stod(coordZ);

                label ownCell = mesh.findCell(pointCoord);
                vector UatSP(vector::zero);
                UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

                Ux[countProbes] = UatSP.x();
                Uy[countProbes] = UatSP.y();
                Uz[countProbes] = UatSP.z();
                cellID[countProbes] = ownCell;

                ++columns;
                if (columns == parameters){
                    columns = 0;
                }
            ++countProbes;
            }
        }
    file.close();

    // Read the size of the abovebottomCells.txt from the file created by topoSet
    std::string abovebottomCells_topoSet = UIntPath + "/constant/polyMesh/sets/abovebottomCells";
    std::fstream file1(abovebottomCells_topoSet);
    GotoLine(file1, 19);

    int abovebottomCellsnumber;
    file1 >> abovebottomCellsnumber;
    //double uTau[abovebottomCellsnumber];
    double shear[abovebottomCellsnumber];

    //========== We calculate the friction coefficient Cf from the bottom wall ==========//
    std::string Cf_cellsIDs = UIntPath + "/abovebottomCells.txt"; 
    std::ifstream fileCellIDs;
    fileCellIDs.open(Cf_cellsIDs);
    point pointCoord1;
    point pointCoord2;
    const volVectorField& C = mesh.C();

    if (fileCellIDs.is_open())
    {
        std::string line;
        int countCells = 0;
        while( std::getline(fileCellIDs,line) )
        {
            std::stringstream ss(line);
            std::string abovebottomcellID;
            std::getline(ss, abovebottomcellID,',');
            label cell2 = std::stoi(abovebottomcellID);
            pointCoord2.y() = C[cell2][1];

            pointCoord1.x() = C[cell2][0];
            //pointCoord2.y() = pointCoord1.y() + deltaChannel/ReTau;
            pointCoord1.y() = 0;
            pointCoord1.z() = C[cell2][2];
            label cell1 = mesh.findCell(pointCoord1);
            shear[countCells] = (U[cell2][0] - U[cell1][0])/(pointCoord2.y()-pointCoord1.y());
            ++countCells;
        }
    }
    fileCellIDs.close();
    double shearGlobal = 0;
    for (int i = 0; i < abovebottomCellsnumber; ++i){
        shearGlobal = shearGlobal + shear[i];
    }

    shearGlobal =  shearGlobal / abovebottomCellsnumber;

    double uTauGlobal = std::sqrt(nu[0] * shearGlobal);

    //========= Calculate the friction coefficient Cf from uTau ==========//
    double Cf = 2*pow(uTauGlobal,2);

    if (cwipiVerbose == 1) Info<< "The Cf is equal to " << Cf << nl << endl;

    /* WE WRITE THE NEW VALUES IN A TXT FILE */
    
    std::ofstream myfile;
    myfile.open(UCfInt, std::ios::out | std::ios::app);

    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Ux[i] << "\n";
    }
    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Uy[i] << "\n";
    }
    for (int i = 0; i < cwipiObsU; ++i){
        myfile << Uz[i] << "\n";
    }

    myfile << Cf << "\n";
    myfile.close();

    //========== Retrieve the cellIDs of my probes ==========//
    if (myGlobalRank == 1){
        std::ofstream myfile_2;
        std::string cellSensors = UIntPath + "/results/cellIDs";

        if (remove(cellSensors.c_str()) != 0)
            perror( "Error deleting file with cell IDs of the sensors" );
        myfile_2.open(cellSensors, std::ios::out | std::ios::app);
            for (int i = 0; i < cwipiObsU; ++i){
                myfile_2 << cellID[i] << "\n";
            }
    }
}

void cwipiSend(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, float cwipiVerbose)
{
    //=== Basic sent with the velocity field, can be automatically interpolated by cwipi if the coupling meshes are different
    //between the OF instance and the KF_coupling code === 

    //* Retrieve the velocity field *
    double* fieldsToSend = new double[3*mesh.nCells()];
    if (cwipiVerbose) Pout << "number of cells: " << mesh.nCells() << endl;
    forAll(mesh.cellCentres(),i)
    {
        fieldsToSend[3*i+0]=vf[i].component(0);
        fieldsToSend[3*i+1]=vf[i].component(1);
        fieldsToSend[3*i+2]=vf[i].component(2);
    }

    //* Use the correct coupling name depending on the OF instance *
    scalar t = runTime.value();

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName,appSuffixChar);

    //* Send and wait for the receive *
    if (cwipiVerbose) Pout<< "Before sending to KF in ensemble "<< appSuffix << endl;  
    cwipi_issend(couplingName,"ex1",sendTag,3,cwipiIteration,t,"u0,v0,w0",fieldsToSend,&status);
    cwipi_wait_issend(couplingName,status);
    //if (cwipiVerbose) Pout<< "After sending to KF in ensemble "<< appSuffix << endl;

    switch(status)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) Info << "Sent Ok in ensemble " << appSuffix << endl << "\n";
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) Info << "Bad receiving in ensemble " << appSuffix << endl << "\n";
        break;
        default :
        if (cwipiVerbose) Info << "Error: bad exchange status in ensemble " << appSuffix << endl << "\n";
    }
    
    delete[] fieldsToSend;
}

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a velocity boundary condition ===
    
    if (cwipiVerbose) Info << "Here we are 11" << endl;
    // double* ParamsToSend = new double[mesh.nCells()];
    double* ParamsToSend = new double[cwipiParams];

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains
    if (cwipiVerbose == 1) Pout << "Before sending params to KF from ensemble " << appSuffix << endl; 
    label top = mesh.boundaryMesh().findPatchID("movingWall");
    double movingWallU = vf.boundaryField()[top][0].component(0);
    
    ParamsToSend[0] = movingWallU;

    //* We use a basic MPI primitive because we don't want any interpolation performed by
    //the cwipi primitive, the destination rank is always 0 because KF_coupling is supposed
    // to always be 0 when we launch a calculation (first in the command line)*
    MPI_Send(ParamsToSend, cwipiParams, MPI_DOUBLE, 0, sendTag_params, MPI_COMM_WORLD);
    if (cwipiVerbose == 1) Pout << "After sending params to KF from ensemble " << appSuffix << endl;

    delete[] ParamsToSend;
}

void cwipiSendParamsChannel(const fvMesh& mesh, const volScalarField& Ck, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose)
{
    //=== Send the parameter to optimize if it is a field of an uniform constant, we can imagine modifying 
    // this with a non uniform field ===

    // double* ParamsToSend = new double[mesh.nCells()];
    double* ParamsToSend = new double[cwipiParams];
    // Info << "number of cells: " << mesh.nCells() << endl;

    if (cwipiVerbose) std::cout << "Before sending params to KF" << std::endl; 
    
    if (cwipiVerbose) std::cout << "CK[0] = " << Ck[0] << std::endl;

    ParamsToSend[0]=Ck[0];

    MPI_Send(ParamsToSend,1,MPI_DOUBLE,0,sendTag_params,MPI_COMM_WORLD);
    if (cwipiVerbose) std::cout << "After sending params to KF" << std::endl;

    delete[] ParamsToSend;
}

void cwipiSendParamsChannelFourrier(int cwipiParams, float cwipiVerbose, std::string stringRootPath)
{
    //=== Send the parameter to optimize for a non uniform field ===

    double* ParamsToSend = new double[cwipiParams];

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    std::string file_path= stringRootPath + "/processor" + std::to_string(myGlobalRank-1) + "/fCoeffs.txt";
    if (cwipiVerbose) std::cout << "fCoeffs path = " << file_path << std::endl;
 
	std::vector<std::vector<string>> file_content;
	std::vector<string> file_row;
	std::string file_line, file_word;
 
	std::fstream file_stream (file_path, std::ios::in);
	if(file_stream.is_open())
	{
		while(getline(file_stream, file_line))
		{
			file_row.clear();
 
			std::stringstream file_str(file_line);
 
			while(getline(file_str, file_word, ','))
				file_row.push_back(file_word);
			file_content.push_back(file_row);
		}
	}
	else
		std::cout << "Could not open the file" << std::endl << "\n";

    if (cwipiVerbose) std::cout << "Before sending params to KF" << std::endl; 

    for (int i=0; i < cwipiParams; i++)
    {
        ParamsToSend[i]=std::stod(file_content[i][0]);
        if (cwipiVerbose) std::cout << "ParamsToSend[i]" << ParamsToSend[i] << std::endl; 
    }

    MPI_Send(ParamsToSend,cwipiParams,MPI_DOUBLE,0,sendTag_params,MPI_COMM_WORLD);

    if (cwipiVerbose) std::cout << "After sending params to KF" << std::endl;

    delete[] ParamsToSend;
}

void cwipiSendParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose)
{
    //=== Send the parameters the optimize, in this case model coefficients ===

    // double* ParamsToSend = new double[mesh.nCells()];
    double* ParamsToSend = new double[cwipiParams];
    // Info << "number of cells: " << mesh.nCells() << endl;

    if (cwipiVerbose) std::cout << "Before sending params to KF" << std::endl; 

    //* We can read the coefficient of a turbulence model from the object turbulence. Be careful,
    // here "turbulence" is the object "turbulence()" and not the pointer "turbulence->" , because we 
    // gave "turbulence()" to the function in the solver code *
    ParamsToSend[0] = readScalar(turbulence.coeffDict().lookup("Cmu"));
    ParamsToSend[1] = readScalar(turbulence.coeffDict().lookup("C1"));
    ParamsToSend[2] = readScalar(turbulence.coeffDict().lookup("C2"));
    ParamsToSend[3] = readScalar(turbulence.coeffDict().lookup("sigmak"));
    ParamsToSend[4] = readScalar(turbulence.coeffDict().lookup("sigmaEps"));

    MPI_Send(ParamsToSend,cwipiParams,MPI_DOUBLE,0,sendTag_params,MPI_COMM_WORLD);
    if (cwipiVerbose) std::cout << "After sending params to KF" << std::endl;

    delete[] ParamsToSend;
}

void cwipiSendParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose)
{
    //=== Send the parameters the optimize, in this case model coefficients ===

    // double* ParamsToSend = new double[mesh.nCells()];
    double* ParamsToSend = new double[cwipiParams];
    // Info << "number of cells: " << mesh.nCells() << endl;

    if (cwipiVerbose == 1) std::cout << "Before sending params to KF" << std::endl; 

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

    // ParamsToSend[0]=Ck[0];

    MPI_Send(ParamsToSend,cwipiParams,MPI_DOUBLE,0,sendTag_params,MPI_COMM_WORLD);
    if (cwipiVerbose) std::cout << "After sending params to KF" << std::endl;

    delete[] ParamsToSend;
}

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, float cwipiVerbose)
{
    //=== Just as the send equivalent, that's the basic receive of the velocity field, can be automatically 
    // interpolated by cwipi if the coupling meshes are different between the OF instance and the KF_coupling code === 
    
    double* fieldsToRecv = new double[3 * mesh.nCells()];
    char cl_receiving_field_name[20];
    sprintf(cl_receiving_field_name, "U_a, V_a, W_a");

    if (cwipiVerbose) Info << "Before Re-receive" << endl;

    scalar t = runTime.value();

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    strcat(couplingName,rankChar);

    cwipi_irecv(couplingName,"ex2",recvTag,3,cwipiIteration,t, cl_receiving_field_name,fieldsToRecv,&status2);

    cwipi_wait_irecv(couplingName, status2);

    switch(status2)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) Info << "Re-receive Ok" << endl;
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) printf("Bad re-receiving\n");
        break;
        default :
        if (cwipiVerbose) printf("Error : bad exchange status\n");
    }

    //* Here the we change the velocity field of the calculation, be careful not to a const volVectorField 
    // in the arguments of the function or it doesn't work *
    forAll(mesh.cells(), i)
    {
        U[i].component(0) = fieldsToRecv[3*i+0];
        U[i].component(1) = fieldsToRecv[3*i+1];
        U[i].component(2) = fieldsToRecv[3*i+2];
    }
    
    if (cwipiVerbose) Info << "After Re-receive" << endl << "\n";

    delete[] fieldsToRecv;
}

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose)
{
    //=== Receive equivalent of the boundary velocity send function ===
    
    double* paramsToRecv = new double[cwipiParams];
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Recv(paramsToRecv,1, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    label top = mesh.boundaryMesh().findPatchID("movingWall");
    fvPatchVectorField& movingWallU = U.boundaryFieldRef()[top];

    if (cwipiVerbose) Info << movingWallU[0].component(0) << endl;
    
    //* If all the velocities are the same througout the entire boundary, the U file in the timeStep resulting
    // folders indicate a uniform boundary condition automatically *
    forAll(movingWallU,faceI)
    {
        movingWallU[faceI]=vector(paramsToRecv[0],0,0);
    }

    if (cwipiVerbose) Info << movingWallU[0].component(0) << endl;
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipiRecvParamsChannel(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose)
{
    //=== Receive equivalent of the constant uniform field  send function ===
    
    double* paramsToRecv = new double[cwipiParams];
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    MPI_Recv(paramsToRecv,1, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    if (cwipiVerbose) Info << Ck[0] << endl;

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

    forAll(Ck,cellI)
    {
        Ck[cellI]=paramsToRecv[0];
        // Ck[cellI]=5;
    }

    if (cwipiVerbose) Info << Ck[0] << endl;
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipiRecvParamsChannelFourrier(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, float cwipiVerbose, std::string stringRootPath, int m, int n, int o)
{
    //=== Receive equivalent of the constant uniform field  send function ===
    
    double* paramsToRecv = new double[cwipiParams];
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    MPI_Recv(paramsToRecv,cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    //Write the new results in a file

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    std::string file_path= stringRootPath + "/processor" + std::to_string(myGlobalRank-1) + "/fCoeffs.txt";
    if (cwipiVerbose) std::cout << "fCoeffs path = " << file_path << std::endl;

    const char* file_path_char = file_path.c_str();

    if (Pstream::master()){
        if (remove(file_path_char) != 0)
            perror( "Error deleting file with coeffs of fourrier expension" );
        else
            if (cwipiVerbose) std::cout << "Fourrier coeffs file " << myGlobalRank << " successfully deleted " << std::endl;
    }

    std::ofstream myfile;
    myfile.open(file_path.c_str(), std::ios::out | std::ios::app);
    if (cwipiVerbose) std::cout << "fCoeffs = "; 
    for (int i = 0; i < cwipiParams; ++i){
        if (i==1 || i==3 || i==5 || i==7 || i==9){
            paramsToRecv[i] = 0;
        }
        myfile << paramsToRecv[i] << "\n";
        if (cwipiVerbose) std::cout << paramsToRecv[i];
    }
    if (cwipiVerbose) std::cout << std::endl;
    myfile.close();

    const double mean = 0;
    const double stddev = 0.11;
    std::normal_distribution<double> dist(mean, stddev);
    float gaussample;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    do
    {
        gaussample=dist(generator);
    }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));

    scalar lambdax = 3*Foam::constant::mathematical::pi;
    scalar lambday = 4;
    scalar lambdaz = Foam::constant::mathematical::pi;

    double alphaTot = 0;
    double betaTot = 0;
    double gammaTot = 0;
    double deltaTot = 0;

    int fourrierindex = 0;
    
    scalar x = 0;
    scalar y = 0;
    scalar z = 0;

    //Calculate the new Ck field
    forAll(mesh.cells(), cellI)
    {
        x = mesh.C()[cellI][0];
        y = mesh.C()[cellI][1];
        z = mesh.C()[cellI][2];

        // if (cwipiVerbose == 1) Info << "x = " << x << " and " << "y = " << y << " and " << "z = " << z << endl;
        
        for (int i=0; i<=m; i++){
            for (int j=1; j<=n; j++){
                for (int k=0; k<=o; k++){
                    alphaTot = alphaTot + paramsToRecv[fourrierindex]*Foam::cos(2*Foam::constant::mathematical::pi*i*x/lambdax)*Foam::sin(2*Foam::constant::mathematical::pi*j*y/lambday)*Foam::cos(2*Foam::constant::mathematical::pi*k*z/lambdaz);
                    // alphaTot = alphaTot + paramsToRecv[fourrierindex]*Foam::cos(2*Foam::constant::mathematical::pi*i*x/lambdax)*Foam::sin(2*Foam::constant::mathematical::pi*j*1/lambday)*Foam::cos(2*Foam::constant::mathematical::pi*k*z/lambdaz);
                    // alphaTot = alphaTot + paramsToRecv[fourrierindex];
                    // if (cwipiVerbose == 1) Info << "paramsToRecv " << paramsToRecv[fourrierindex] << endl;
                    fourrierindex = fourrierindex +1;
                }
            }
        }

        for (int i=0; i<=m; i++){
            for (int j=1; j<=n; j++){
                for (int k=1; k<=o; k++){
                    betaTot = betaTot + paramsToRecv[fourrierindex]*Foam::cos(2*Foam::constant::mathematical::pi*i*x/lambdax)*Foam::sin(2*Foam::constant::mathematical::pi*j*y/lambday)*Foam::sin(2*Foam::constant::mathematical::pi*k*z/lambdaz);
                    fourrierindex = fourrierindex +1;
                }
            }
        }

        for (int i=1; i<=m; i++){
            for (int j=1; j<=n; j++){
                for (int k=0; k<=o; k++){
                    gammaTot = gammaTot + paramsToRecv[fourrierindex]*Foam::sin(2*Foam::constant::mathematical::pi*i*x/lambdax)*Foam::sin(2*Foam::constant::mathematical::pi*j*y/lambday)*Foam::cos(2*Foam::constant::mathematical::pi*k*z/lambdaz);
                    fourrierindex = fourrierindex +1;
                }
            }
        }

        for (int i=1; i<=m; i++){
            for (int j=1; j<=n; j++){
                for (int k=1; k<=o; k++){
                    deltaTot = deltaTot + paramsToRecv[fourrierindex]*Foam::sin(2*Foam::constant::mathematical::pi*i*x/lambdax)*Foam::sin(2*Foam::constant::mathematical::pi*j*y/lambday)*Foam::sin(2*Foam::constant::mathematical::pi*k*z/lambdaz);
                    fourrierindex = fourrierindex +1;
                }
            }
        }

        // if (cwipiVerbose == 1) Info << "alphaTot = " << alphaTot << " and " << "betaTot = " << betaTot << " and " << "gammaTot = " << gammaTot << " and " << "deltaTot = " << deltaTot << endl;

        Ck[cellI] = alphaTot + betaTot + gammaTot + deltaTot;
        // Ck[cellI]=5;

        alphaTot = 0;
        betaTot = 0;
        gammaTot = 0;
        deltaTot = 0;
        fourrierindex = 0;

        if (Ck[cellI] < 0.005) 
        {
            Ck[cellI] = 0.01 + 0.01*gaussample;
        }
    }

    // std::ofstream myfile;
    // myfile.open(file_path.c_str(), std::ios::out | std::ios::app);
    // if (cwipiVerbose == 1) std::cout << "fCoeffs = "; 
    // for (int i = 0; i < cwipiParams; ++i){
    //     myfile << Ck[10] << "\n";
    //     if (cwipiVerbose == 1) std::cout << Ck[10];
    // }
    // if (cwipiVerbose == 1) std::cout << std::endl;
    // myfile.close();
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipiRecvParamsChannelGaussian(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, float cwipiVerbose, std::string stringRootPath)
{
    //=== Receive equivalent of the constant uniform field send function ===
    
    double* paramsToRecv = new double[cwipiParams]; //Number of parameters is 2*number of gaussian functions
    // double* yDistrib = new double[cwipiParams]; //Number of points in the yDistrib is 2*number of functions since we apply symetry on the functions
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    MPI_Recv(paramsToRecv,cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    //*** Write the new params in a file ***

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    std::string file_path= stringRootPath + "/processor" + std::to_string(myGlobalRank-1) + "/gCoeffs.txt";
    if (cwipiVerbose) std::cout << "gCoeffs path = " << file_path << std::endl;

    const char* file_path_char = file_path.c_str();

    if (Pstream::master()){
        if (remove(file_path_char) != 0)
            perror( "Error deleting file with coeffs of fourrier expension" );
        else
            if (cwipiVerbose) std::cout << "Gaussian coeffs file " << myGlobalRank << " successfully deleted " << std::endl;
    }

    std::ofstream myfile;
    myfile.open(file_path.c_str(), std::ios::out | std::ios::app);
    if (cwipiVerbose) std::cout << "gCoeffs = "; 
    for (int i = 0; i < cwipiParams; ++i){
        myfile << paramsToRecv[i] << "\n";
        if (cwipiVerbose) std::cout << paramsToRecv[i];
    }
    if (cwipiVerbose) std::cout << std::endl;
    myfile.close();

    //*** Calculation of the new Ck field ***

    //== We fisrt open the y distribution file for the gaussian function ==
    // if (cwipiVerbose == 1) std::cout << "Opening y distribution file " << std::endl;
    
    // std::string yDistrib_file= stringRootPath + "/yDistrib.txt";
    // if (cwipiVerbose == 1) std::cout << "yDistrib path = " << stringRootPath << std::endl;
 
	// std::vector<std::vector<string>> yDistrib_content;
	// std::vector<string> yDistrib_row;
	// std::string yDistrib_line, yDistrib_word;
 
	// std::fstream yDistrib_stream (yDistrib_file, std::ios::in);
	// if(yDistrib_stream.is_open())
	// {
	// 	while(getline(yDistrib_stream, yDistrib_line))
	// 	{
	// 		yDistrib_row.clear();
 
	// 		std::stringstream yDistrib_str(yDistrib_line);
 
	// 		while(getline(yDistrib_str, yDistrib_word, ','))
	// 			yDistrib_row.push_back(yDistrib_word);
	// 		yDistrib_content.push_back(yDistrib_row);
	// 	}
	// }
	// else
	// 	std::cout << "Could not open the file" << std::endl << "\n";

    // if (cwipiVerbose == 1) std::cout << "yDistrib file is open, begin to fill array" << std::endl;
 
    // for(int i = 0; i < cwipiParams/2; i++)
    // {
    //     yDistrib[i]=std::stod(yDistrib_content[i][0]);
    // }
    // // if (cwipiVerbose == 1) std::cout << "yDistrib[i] = " << yDistrib[0] << std::endl;

    // for(int i = 0; i < cwipiParams/2; i++)
    // {
    //     yDistrib[cwipiParams/2+i]=2-std::stod(yDistrib_content[cwipiParams/2-i-1][0]);
    //     // if (cwipiVerbose == 1) std::cout << "yDistrib[i] = " << yDistrib[cwipiParams/2+i] << std::endl;
    // }

    // //== Then we create a randomazier in case Ck values are negatives ==
    const double mean = 0;
    const double stddev = 0.11;
    std::normal_distribution<double> dist(mean, stddev);
    float gaussample;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    do
    {
        gaussample=dist(generator);
    }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));


    // //== We calculate the new Ck field ==
    double CkTemp = 0;
    scalar y = 0;

    forAll(mesh.cells(), cellI)
    {
        y = mesh.C()[cellI][1];
        // if (cwipiVerbose == 1) Info << "y = " << y << endl;
        
        for (int i=0; i<cwipiParams/3; i++){
            // CkTemp = CkTemp + Foam::exp(paramsToRecv[2*i]-Foam::pow((y-yDistrib[i]),2)/pow(paramsToRecv[2*i+1],2));
            // CkTemp = CkTemp + Foam::exp(paramsToRecv[2*i]-Foam::pow((y-yDistrib[cwipiParams-i-1]),2)/pow(paramsToRecv[2*i+1],2));
            CkTemp = CkTemp + Foam::exp(paramsToRecv[3*i]-Foam::pow((y-paramsToRecv[3*i+2]),2)/pow(paramsToRecv[3*i+1],2));
            CkTemp = CkTemp + Foam::exp(paramsToRecv[3*i]-Foam::pow((y-(2-paramsToRecv[3*i+2])),2)/pow(paramsToRecv[3*i+1],2));
        }

        Ck[cellI] = CkTemp;
        // Ck[cellI]=5;
        if (Ck[cellI] < 0.005) 
        {
            Ck[cellI] = 0.01 + 0.01*gaussample;
        }
        CkTemp = 0;
    }
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
    // delete[] yDistrib;
}

void cwipiRecvParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose, std::string globalRootPath)
{
    //=== Receive equivalent of the model coefficients send function, here we need to overwrite the momentumTransport dictionary 
    // And then read again the coefficient in the solver (turbulence->read()). We need to declare a OFstream containing the stream
    // of the file to overwrite. The stream is then closed at the end of this function. ===
    
    double* paramsToRecv = new double[cwipiParams];

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    std::string dictPath = globalRootPath+"/processor"+std::to_string(myGlobalRank-1)+"/constant/momentumTransport";

    if (cwipiVerbose) Info << "dictPath = " << dictPath << endl; 

    Foam::OFstream FileStream(dictPath);
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    MPI_Recv(paramsToRecv,cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);
    const double mean = 0.0;
    const double stddev = 0.11;
    std::normal_distribution<double> dist(mean, stddev);
    float gaussample;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    
    if (paramsToRecv[0] < 0)
    {
        do
        {
            gaussample=dist(generator);
        }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));
        paramsToRecv[0] = 0.01 + 0.01*gaussample;
    }
    for (int i=1; i<cwipiParams; i++)
    {
        if (paramsToRecv[i] < 0.05)
        {
            do
            {
                gaussample=dist(generator);
            }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));
            paramsToRecv[i] = 0.1 + 0.1*gaussample;
        }
    }

    if (cwipiVerbose) Info << turbulence.subDict("RAS") << endl;

    turbulence.subDict("RAS").lookup("Cmu")[0] = paramsToRecv[0];
    turbulence.subDict("RAS").lookup("C1")[0] = paramsToRecv[1];
    turbulence.subDict("RAS").lookup("C2")[0] = paramsToRecv[2];
    turbulence.subDict("RAS").lookup("sigmak")[0] = paramsToRecv[3];
    turbulence.subDict("RAS").lookup("sigmaEps")[0] = paramsToRecv[4];

    turbulence.writeHeader(FileStream);  //First overwrite with the OF header
    turbulence.dictionary::write(FileStream, false); //Then overwrite the coefficients. "false" allows to write everything and not just the coeffs

    if (cwipiVerbose) Info << turbulence.subDict("RAS") << endl;
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipiRecvParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, int partsReparty, float cwipiVerbose, std::string globalRootPath)
{
    //=== Receive equivalent of the model coefficients send function, here we need to overwrite the momentumTransport dictionary 
    // And then read again the coefficient in the solver (turbulence->read()). We need to declare a OFstream containing the stream
    // of the file to overwrite. The stream is then closed at the end of this function. ===
    
    double* paramsToRecv = new double[cwipiParams];

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    std::string dictPath = globalRootPath+"/processor"+std::to_string(myGlobalRank-1)+"/constant/momentumTransport";

    if (cwipiVerbose) Info << "dictPath = " << dictPath << endl; 

    Foam::OFstream FileStream(dictPath);
    
    if (cwipiVerbose) Info << "Before Re-receive params" << endl;

    MPI_Status status4;
    MPI_Recv(paramsToRecv,cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

    const double mean = 0.0;
    const double stddev = 0.11;
    std::normal_distribution<double> dist(mean, stddev);
    float gaussample;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    for (int i=0; i<cwipiParams; i++)
    {
        if (paramsToRecv[i] < 0)
        {
            do
            {
                gaussample=dist(generator);
            }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));
            paramsToRecv[i] = 0.01 + 0.01*gaussample;
        }
    }
    if (cwipiVerbose) Info << turbulence.subDict("RAS") << endl;

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

    if (cwipiVerbose) Info << turbulence.subDict("RAS") << endl;
    
    if (cwipiVerbose) Info << "After Re-receive params" << endl << "\n";

    delete[] paramsToRecv;
}

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, float cwipiVerbose)
{
    //=== Delete the coupling for each OF instance and delete the mesh arrays ===
    
    if (cwipiVerbose) Info << "Delete Cwipi coupling from OF" << endl << "\n";

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    strcat(couplingName,rankChar);

    cwipi_delete_coupling(couplingName);
    delete[] pointCoords;
    // delete[] connecIdx;
    // delete[] connec;
    delete[] face_index;
    delete[] cell_to_face_connectivity;
    delete[] face_connectivity_index;
    delete[] face_connectivity;
}
}

// ************************************************************************* //
