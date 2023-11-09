/**
 * \file cwipiPstreamPar.C
 * \brief Defines all functions needed to employ CWIPI routines
 */

#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "interpolationCellPointWallModified.H"
// #include "interpolationCellPointFace.H"
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
static int recvTag = 2; /*!< Tag for receive field from EnKF to OpenFOAM */
static int status; /*!< Status for send field from OpenFOAM to EnKF */
static int status2;

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
    label nMyProc=Pstream::myProcNo();
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    /* Creation of the coupling and coupling mesh
    
    Options :
    CWIPI_COUPLING_SEQUENTIAL
    CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
    CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING

    */

   if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The name of my coupling is " << couplingName <<endl;   

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
    
    if (cwipiVerbose) if (Pstream::master()) printf("%s: Create mesh in OF \n", __func__);
    cwipi_locate(couplingName);

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "Cwipi locate passed in member " << appSuffix << Foam::nl << Foam::endl;
    delete[] face_index_temp;
	delete[] face_index_temp_IDs;
	delete[] face_connectivity_index_temp;
	delete[] face_connectivity_index_temp_IDs;
}

void cwipiSend(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Basic sent with the velocity field, can be automatically interpolated by cwipi if the coupling meshes are different
    //between the OF instance and the KF_coupling code === 

    //* Retrieve the velocity field *
    double* fieldsToSend = new double[3*mesh.nCells()];
    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiSend function: " << mesh.nCells() << endl;
    forAll(mesh.cellCentres(),i)
    {
        fieldsToSend[3*i+0]=vf[i].component(0);
        fieldsToSend[3*i+1]=vf[i].component(1);
        fieldsToSend[3*i+2]=vf[i].component(2);
    }

    //* Use the correct coupling name depending on the OF instance *
    scalar t = runTime.value();

    int myGlobalRank;
    label nMyProc=Pstream::myProcNo();
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition

    char couplingName[250] = {"cwipiFoamCoupling"};


    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    //* Send and wait for the receive *
    if (cwipiVerbose) if (Pstream::master()) Pout<< "Before sending to KF in member "<< appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;  
    cwipi_issend(couplingName, "ex1", sendTag, 3, cwipiIteration, t, "u0,v0,w0", fieldsToSend, &status);
    cwipi_wait_issend(couplingName, status);
    //if (cwipiVerbose) Pout<< "After sending to KF in member "<< appSuffix << endl;

    switch(status)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) if (Pstream::master()) Pout << "Sent Ok in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl << "\n";
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) Pout << "Bad receiving in member " << appSuffix << endl << "\n";
        break;
        default :
        if (cwipiVerbose) Pout << "Error: bad exchange status in member " << appSuffix << endl << "\n";
    }
    
    delete[] fieldsToSend;
}

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Just as the send equivalent, that's the basic receive of the velocity field, can be automatically 
    // interpolated by cwipi if the coupling meshes are different between the OF instance and the KF_coupling code === 
    
    double* fieldsToRecv = new double[3 * mesh.nCells()];
    char cl_receiving_field_name[20];
    sprintf(cl_receiving_field_name, "U_a, V_a, W_a");

    scalar t = runTime.value();

    int myGlobalRank;
    label nMyProc=Pstream::myProcNo();

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    if (cwipiVerbose) if (Pstream::master()) Pout << "Before receive back state in member " << appSuffix << " the coupling name for the receive is " << couplingName << endl;


    cwipi_irecv(couplingName,"ex2", recvTag, 3, cwipiIteration, t, cl_receiving_field_name, fieldsToRecv, &status2);

    cwipi_wait_irecv(couplingName, status2);

    switch(status2)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) Pout << "Back-receive Ok in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) Pout << "Bad Back-receiving in member " << appSuffix << "\n";
        break;
        default :
        if (cwipiVerbose) Pout << "Error : bad exchange status in member " << appSuffix << "\n";
    }

    //* Here the we change the velocity field of the calculation, be careful not to a const volVectorField 
    // in the arguments of the function or it doesn't work *
    forAll(mesh.cells(), i)
    {
        U[i].component(0) = fieldsToRecv[3*i+0];
        U[i].component(1) = fieldsToRecv[3*i+1];
        U[i].component(2) = fieldsToRecv[3*i+2];
    }
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "After receive back state in member " << appSuffix << endl << "\n";

    delete[] fieldsToRecv;
}

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int nbParts, float cwipiVerbose)
{
    //=== Delete the coupling for each OF instance and delete the mesh arrays ===
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "Delete Cwipi coupling from OF" << endl << "\n";

    int myGlobalRank;
    label nMyProc=Pstream::myProcNo();
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;
    

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
