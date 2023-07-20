/**
 * \file cwipiPstreamPar.C
 * \brief Defines all functions needed to employ CWIPI routines
 */

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

/**
 * Declaration of variables for tags and status of the exchanges
 */ 

static int sendTag = 1; /*!< Tag for send field from OpenFOAM to EnKF */
static int sendTag_params = 3; /*!< Tag for send parameters from OpenFOAM to EnKF */
static int recvTag = 2; /*!< Tag for receive field from EnKF to OpenFOAM */
static int recvTag_params = 4; /*!< Tag for receive parameters from EnKF to OpenFOAM */
static int status; /*!< Status for send field from OpenFOAM to EnKF */
static int status2;
MPI_Status status4;

/*! \fn void addControlParams(int numberCwipiPhase, double deltaT, double currentTime)
    \brief Add control parameters that must be sent to KF_coupling code
    \param numberCwipiPhase Number of DA phases
    \param deltaT Time step of CFD loop
    \param currentTime Start time of CFD simulation
*/

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime)
{

    cwipi_add_local_int_control_parameter("numberCwipiPhase", numberCwipiPhase);
    cwipi_add_local_double_control_parameter("deltaT", deltaT);
    cwipi_add_local_double_control_parameter("currentTime", currentTime);
}

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, 
int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol)
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
    int appSuffix = floor((myGlobalRank + 1)/nbParts);

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

void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int mainsubDomain, int nbParts, 
interpolationCellPointWallModified<vector> triangulateCellsU, float cwipiVerbose, std::string globalPath, std::string UIntPath)
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

    int appSuffix = floor((myGlobalRank + 1)/nbParts); // Change 2 by the number of subdomains

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

    remove(UInt);

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
    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    cwipi_irecv(couplingName,"ex2", recvTag, 3, cwipiIteration, t, cl_receiving_field_name, fieldsToRecv, &status2);

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

    MPI_Recv(paramsToRecv, cwipiParams, MPI_DOUBLE, 0, recvTag_params, MPI_COMM_WORLD, &status4);

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

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, float cwipiVerbose)
{
    //=== Delete the coupling for each OF instance and delete the mesh arrays ===
    
    if (cwipiVerbose) Info << "Delete Cwipi coupling from OF" << endl << "\n";

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank+1)/2); // Change 2 by the number of subdomains

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

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
