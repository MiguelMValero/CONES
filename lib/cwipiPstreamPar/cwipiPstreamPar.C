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
 * @dir ./lib/cwipiPstreamPar
 * @brief **Receiving and sending data, parallel computing and test case dependencies.**
 */

/**
 * @file cwipiPstreamPar.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel Martinez Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Defines all functions needed to employ CWIPI routines
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
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

static int sendTag = 10;     /*!< Tag for send field from OpenFOAM to HLEnKF */
static int recvTag = 12;     /*!< Tag for receive field from HLEnKF to OpenFOAM */

// static int sendTag_MGEnKF = 20; /*!< Tag for send field from OpenFOAM to MGEnKF */
// static int recvTag_MGEnKF = 22; /*!< Tag for receive field from MGEnKF to OpenFOAM */

// static int sendTag_MFEnKF = 30; /*!< Tag for send field from OpenFOAM to MFEnKF */
// static int recvTag_MFEnKF = 32; /*!< Tag for receive field from MFEnKF to OpenFOAM */

static int status;   // Status for send field from OpenFOAM to EnKF */
static int status2;
static int status5; // cwipi status for send
static int status7; // cwipi status for recv


// === MFEnKF status for send for principal projection to coarse grid === //
static int status100; // status for send
static int status200; // status for recv

// static int sendTag_proj = 101; 
// static int recvTag_proj = 201; 



//===========================================================================
void tic(int mode) {
        static std::chrono::_V2::system_clock::time_point t_start;

        if (mode==0)
                t_start = std::chrono::high_resolution_clock::now();
        else {
              	auto t_end = std::chrono::high_resolution_clock::now();
                std::cout << "Elapsed time for total analysis phase is " << (t_end - t_start).count()*1E-9 << " seconds\n";
        }
}

void toc() { tic(1); }

//===========================================================================
/*!
    @brief Add control parameters that must be sent to KF_coupling code.
    @param numberCwipiPhase Number of data assimilation phases.
    @param deltaT Time step of CFD loop.
    @param currentTime Start time of CFD simulation.
*/
void addControlParams(int numberCwipiPhase, double deltaT, double currentTime)
{

    cwipi_add_local_int_control_parameter("numberCwipiPhase", numberCwipiPhase);
    cwipi_add_local_double_control_parameter("deltaT", deltaT);
    cwipi_add_local_double_control_parameter("currentTime", currentTime);
}


void deleteControlParams()
{

    cwipi_delete_local_int_control_parameter("numberCwipiPhase");
    cwipi_delete_local_double_control_parameter("deltaT");
    cwipi_delete_local_double_control_parameter("currentTime");
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


/* **************************************************** */
/* ** ** ** ** MGEnKF specific functions ** ** ** ** ** */
/* **************************************************** */

void cwipiCouplingMGEnKF(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol_fine, double geom_tol_coarse, int subdomains, int cwipiMembers)
{
    //=== Function that deals with the coupling of each OF instance ===
    
    //* Calculations for declaring the coupling mesh *
    
    int nCells = mesh.nCells();
    Foam::Info << "DD Number of cells from cwipiCouplingMGEnKF func = " << nCells << Foam::endl;
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

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "PP The name of my coupling is " << couplingName <<endl;   

    Foam::Info << "About to create EnKF coupling in cwipiPstreamPar.C..." << Foam::endl;

    cwipi_create_coupling(couplingName,
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                        "MGEnKF",
                        3,
                        geom_tol_coarse,
                        CWIPI_STATIC_MESH,
                        CWIPI_SOLVER_CELL_CENTER,
                        1,
                        "Ensight Gold",
                        "text");

    cwipi_synchronize_control_parameter("MGEnKF");

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

void cwipiSendVolVectorFieldMGEnKF(const fvMesh& mesh, const volVectorField& vField, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Basic sent with the velocity field, can be automatically interpolated by cwipi if the coupling meshes are different
    //between the OF instance and the KF_coupling code === 

    // ===================================================================== //
    // Send fields from OpenFOAM (cwipiMultigridIcoFoamPar) to Kalman Filter //
    // ===================================================================== //

    //* Retrieve the velocity field *
    double* fieldsToSendFine = new double[3*mesh.nCells()];
    // double* fieldsToSend_2 = new double[3*mesh.nCells()];
    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiSendVolVectorFieldMGEnKF function: " << mesh.nCells() << endl;

    forAll(mesh.cellCentres(),i)
    {
        fieldsToSendFine[3*i+0]=vField[i].component(0);
        fieldsToSendFine[3*i+1]=vField[i].component(1);
        fieldsToSendFine[3*i+2]=vField[i].component(2);

        // Foam::Pout << "Indice i pour fieldsToSend dans cwipiSendVolVectorFieldMGEnKF = " << i << Foam::endl;
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

    Info << "THE COUPLING NAME FROM THE OF SIDE IS " << couplingName << endl;

    //* Send and wait for the receive *
    if (cwipiVerbose) if (Pstream::master()) Pout<< "Before sending to MGEnKF in member "<< appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;  
    cwipi_issend(couplingName, "ex1", sendTag, 3, cwipiIteration, t, "uf0,vf0,wf0", fieldsToSendFine, &status5);

    // std::cout << "fieldsToSendFine =" << fieldsToSendFine << std::endl;

    if (cwipiVerbose) if (Pstream::master()) Pout<< "cwipi_issend passed for fine member" << endl; 

    cwipi_wait_issend(couplingName, status5);

    if (cwipiVerbose) if (Pstream::master()) Pout<< "cwipi_wait_issend passed for fine member" << endl; 

    if (cwipiVerbose) Pout<< "After sending to KF in member "<< appSuffix << endl;

    switch(status5)
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
    
    delete[] fieldsToSendFine;
}


void cwipiRecvVolVectorFieldMGEnKF(const fvMesh& mesh, volVectorField& vField, volVectorField& stateMatricesDiff, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Just as the send equivalent, that's the basic receive of the velocity field, can be automatically 
    // interpolated by cwipi if the coupling meshes are different between the OF instance and the KF_coupling code === 

    // ================================================================== //
    // Receives fields from MGEnKF to OpenFOAM (cwipiMultigridIcoFoamPar) //
    // ================================================================== //
    
    double* fieldsToRecvFine = new double[3 * mesh.nCells()];
    // double* fieldsToRecv_2 = new double[3 * mesh.nCells()];

    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiRecvVolVectorFieldMGEnKF function: " << mesh.nCells() << endl;

    char cl_receiving_field_name[20];
    // char cl_receiving_field_name_2[60];
    sprintf(cl_receiving_field_name, "Uf_a, Vf_a, Wf_a");
    // sprintf(cl_receiving_field_name_2, "U_fineForecast, V_fineForecast, W_fineForecast");

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

    if (cwipiVerbose) if (Pstream::master()) Pout << "[MGEnKF] Before receive back state in member " << appSuffix << " the coupling name for the receive is " << couplingName << endl;

    cwipi_irecv(couplingName,"ex3", recvTag, 3, cwipiIteration, t, cl_receiving_field_name, fieldsToRecvFine, &status7);

    Foam::Info << "cwipi_irecv func from cwipiRecvVolVectorFieldMGEnKF function from cwipiPstreamPar.C passed" << Foam::endl;

    cwipi_wait_irecv(couplingName, status7);

    Foam::Info << "cwipi_wait_irecv func from cwipiRecvVolVectorFieldMGEnKF function from cwipiPstreamPar.C passed" << Foam::endl;

    switch(status7)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) if (Pstream::master()) Pout << "Back-receive Ok in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
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
        stateMatricesDiff[i].component(0) = fieldsToRecvFine[3*i+0];
        stateMatricesDiff[i].component(1) = fieldsToRecvFine[3*i+1];
        stateMatricesDiff[i].component(2) = fieldsToRecvFine[3*i+2];

        vField[i].component(0) = vField[i].component(0) - stateMatricesDiff[i].component(0);
        vField[i].component(1) = vField[i].component(1) - stateMatricesDiff[i].component(1);
        vField[i].component(2) = vField[i].component(2) - stateMatricesDiff[i].component(2);
    }
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "After receive back state in member " << appSuffix << endl << "\n";

    delete[] fieldsToRecvFine;
}



/* **************************************************** */
/* ** ** ** ** MFEnKF specific functions ** ** ** ** ** */
/* **************************************************** */

/**
 * 
 * @brief Coupling function for principal members projection to control members.
 * 
 * @param mesh OpenFOAM mesh.
 * @param pointCoords Point coordiantes.
 * @param face_index Face indices.
 * @param face_connectivity_index Face connectivity indices.
 * @param cell_to_face_connectivity Cell to face connectivity.
 * @param face_connectivity Face connectivity.
 * @param c2fconnec_size Cell to face connectivity size.
 * @param fconnec_size face connectivity size.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 * @param geom_tol Geometric tolerance.
 */
void cwipiCouplingProjection(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, 
int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol)
{
    //=== Function that deals with the coupling of each OF instance ===
    
    // Calculations for declaring the coupling mesh *
    
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
    int myGlobalSize;
    label nMyProc=Pstream::myProcNo();
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &myGlobalSize);
    char couplingName[250] = {"cwipiFoamCoupling_forPrincipalVariateProjection"};
    int nMembers = (myGlobalSize / nbParts);
    std::cout << "nMembers from cwipiCouplingProjection = " << nMembers << std::endl;
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    // Creation of the coupling and coupling mesh
    
    /*Options :
    CWIPI_COUPLING_SEQUENTIAL
    CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
    CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING*/

    

   if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The name of my coupling is " << couplingName <<endl;   

    cwipi_create_coupling(couplingName,
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          "MFEnKF",
                          3,
                          geom_tol,
                          CWIPI_STATIC_MESH,
                          CWIPI_SOLVER_CELL_CENTER,
                          1,
                          "Ensight Gold",
                          "text");

    cwipi_synchronize_control_parameter("MFEnKF");
    cwipi_dump_application_properties();
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

/**
 * @brief Coupling function for MFEnKF.
 * 
 * @param mesh OpenFOAM mesh.
 * @param pointCoords Point coordinates.
 * @param face_index Face indices.
 * @param face_connectivity_index face connectivity indices.
 * @param cell_to_face_connectivity Cell to face connectivity.
 * @param face_connectivity Face connectivity.
 * @param c2fconnec_size Cell to face conenctivity size.
 * @param fconnec_size Face connectivity size.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 * @param geom_tol_fine Fine geometric tolerance.
 * @param geom_tol_coarse Coarse geometric tolerance.
 * @param allMembers All ensemble members (principal + control + ancillary).
 */
void cwipiCouplingMFEnKF(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol_fine, double geom_tol_coarse, int allMembers)
{
    //=== Function that deals with the coupling of each OF instance ===
    
    //* Calculations for declaring the coupling mesh *
    
    int nCells = mesh.nCells();
    Foam::Info << "Number of cells from cwipiCouplingMGEnKF func = " << nCells << Foam::endl;
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
    std::cout << "myGlobalRank = " << myGlobalRank << std::endl;
    label nMyProc=Pstream::myProcNo();
    std::cout << "nMyProc = " << nMyProc << std::endl;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling"};
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;
    std::cout << "appSuffix = " << appSuffix << std::endl;

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    /* Creation of the coupling and coupling mesh
    
    Options :
    CWIPI_COUPLING_SEQUENTIAL
    CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING
    CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING

    */

   if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "PP The name of my coupling is " << couplingName <<endl;   

    if (myGlobalRank > (1+(nbParts*(allMembers+1)) - nbParts)) 
    {
        
        Foam::Info << "About to create EnKF coupling in cwipiPstreamPar.C..." << Foam::endl;

        cwipi_create_coupling(couplingName,
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "MFEnKF",
                            3,
                            geom_tol_coarse,
                            CWIPI_STATIC_MESH,
                            CWIPI_SOLVER_CELL_CENTER,
                            1,
                            "Ensight Gold",
                            "text");

        cwipi_synchronize_control_parameter("MFEnKF");
        cwipi_dump_application_properties();
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
    
    else {

        Foam::Info << "About to create MFEnKF coupling in cwipiPstreamPar.C..." << Foam::endl;
        
        cwipi_create_coupling(couplingName,
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            "MFEnKF",
                            3,
                            geom_tol_fine,
                            CWIPI_STATIC_MESH,
                            CWIPI_SOLVER_CELL_CENTER,
                            1,
                            "Ensight Gold",
                            "text");

        cwipi_synchronize_control_parameter("MFEnKF");
        cwipi_dump_application_properties();
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
    
    //if (cwipiVerbose == 2) cwipi_dump_application_properties();

}

/**
 * @brief Send volume vector fields.
 * 
 * @param mesh OpenFOAM mesh.
 * @param vField Volume vector field.
 * @param runTime Current runTime.
 * @param cwipiIteration 
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiSendVolVectorField(const fvMesh& mesh, const volVectorField& vField, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Basic sent with the velocity field, can be automatically interpolated by cwipi if the coupling meshes are different
    //between the OF instance and the KF_coupling code === 


    // ============================================================ //
    // Send fields from OpenFOAM (cwipiIcoFoamPar) to Kalman Filter //
    // ============================================================ //

    //* Retrieve the velocity field *
    double* fieldsToSend = new double[3*mesh.nCells()];

    std::cout << "Current number of cells = " << mesh.nCells() << std::endl;

    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiSendVolVectorField function: " << mesh.nCells() << endl;

    forAll(mesh.cellCentres(),i)
    {
        fieldsToSend[3*i+0]=vField[i].component(0);
        fieldsToSend[3*i+1]=vField[i].component(1);
        fieldsToSend[3*i+2]=vField[i].component(2);
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
    if (cwipiVerbose) if (Pstream::master()) Pout<< "Before sending to MFEnKF in member "<< appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;


    cwipi_issend(couplingName, "ex1", sendTag, 3, cwipiIteration, t, "u0,v0,w0", fieldsToSend, &status);

    Foam::Info << "cwipi_issend function from cwipiPstreamPar.C passed" << Foam::endl;

    cwipi_wait_issend(couplingName, status);

    Foam::Info << "cwipi_wait_issend function passed" << Foam::endl;
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

/**
 * @brief Send volume vector fields for intermediate projection from principal to control members.
 * 
 * @param mesh OpenFOAM mesh.
 * @param vField Volume vector field.
 * @param runTime Current runTime.
 * @param cwipiIteration 
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiSendVolVectorFieldMFEnKF(const fvMesh& mesh, const volVectorField& vField, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Basic sent with the velocity field, can be automatically interpolated by cwipi if the coupling meshes are different
    //between the OF instance and the KF_coupling code === 

    // ===================================================================== //
    // Send fields from OpenFOAM (cwipiMultigridIcoFoamPar) to Kalman Filter //
    // ===================================================================== //

    //* Retrieve the velocity field *
    double* fieldsToSendFine = new double[3*mesh.nCells()];
    // double* fieldsToSend_2 = new double[3*mesh.nCells()];
    if (cwipiVerbose) Pout << "Number of cells in the cwipiSendVolVectorFieldMFEnKF function: " << mesh.nCells() << endl;

    forAll(mesh.cellCentres(),i)
    {
        fieldsToSendFine[3*i+0]=vField[i].component(0);
        fieldsToSendFine[3*i+1]=vField[i].component(1);
        fieldsToSendFine[3*i+2]=vField[i].component(2);

        // Foam::Pout << "Indice i pour fieldsToSend dans cwipiSendVolVectorFieldMGEnKF = " << i << Foam::endl;
    }

    //* Use the correct coupling name depending on the OF instance *
    scalar t = runTime.value();

    int myGlobalRank;
    label nMyProc=Pstream::myProcNo();
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition

    char couplingName[250] = {"cwipiFoamCoupling_forPrincipalVariateProjection"};


    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    Info << "THE COUPLING NAME FROM THE OF SIDE IS " << couplingName << endl;

    //* Send and wait for the receive *
    if (cwipiVerbose) if (Pstream::master()) Pout<< "Before sending to MFEnKF in member "<< appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;  
    cwipi_issend(couplingName, "ex1", sendTag, 3, cwipiIteration, t, "u_proj,v_proj,w_proj", fieldsToSendFine, &status100);

    // std::cout << "fieldsToSendFine =" << fieldsToSendFine << std::endl;

    if (cwipiVerbose) if (Pstream::master()) Pout<< "cwipi_issend passed for projection" << endl; 

    cwipi_wait_issend(couplingName, status100);

    if (cwipiVerbose) if (Pstream::master()) Pout<< "cwipi_wait_issend passed for projection" << endl; 

    if (cwipiVerbose) Pout<< "After sending to KF in member "<< appSuffix << endl;

    switch(status100)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) if (Pstream::master()) Pout << "Sent Ok projection in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl << "\n";
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) Pout << "Bad receiving projection in member " << appSuffix << endl << "\n";
        break;
        default :
        if (cwipiVerbose) Pout << "Error: bad exchange status projection in member " << appSuffix << endl << "\n";
    }
    
    delete[] fieldsToSendFine;
}

/**
 * @brief Receive volume vector fields.
 * 
 * @param mesh OpenFOAM mesh.
 * @param vField Volume vector field.
 * @param runTime Current runTime.
 * @param cwipiIteration 
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiRecvVolVectorField(const fvMesh& mesh, volVectorField& vField, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Just as the send equivalent, that's the basic receive of the velocity field, can be automatically 
    // interpolated by cwipi if the coupling meshes are different between the OF instance and the KF_coupling code === 

    // ======================================================= //
    // Receives fields from EnKF to OpenFOAM (cwipiIcoFoamPar) //
    // ======================================================= //
    
    double* fieldsToRecv = new double[3 * mesh.nCells()];
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiRecvVolVectorField function: " << mesh.nCells() << endl;

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

    if (cwipiVerbose) if (Pstream::master()) Pout << "[EnKF] Before receive back state in member " << appSuffix << " the coupling name for the receive is " << couplingName << endl;


    cwipi_irecv(couplingName,"ex2", recvTag, 3, cwipiIteration, t, cl_receiving_field_name, fieldsToRecv, &status2);

    Foam::Info << "cwipi_irecv func from cwipiRecvVolVectorField function from cwipiPstreamPar.C passed" << Foam::endl;

    cwipi_wait_irecv(couplingName, status2);

    Foam::Info << "cwipi_wait_irecv func from cwipiRecvVolVectorField function from cwipiPstreamPar.C passed" << Foam::endl;


    switch(status2)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) if (Pstream::master()) Pout << "Back-receive Ok in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
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
        vField[i].component(0) = fieldsToRecv[3*i+0];
        vField[i].component(1) = fieldsToRecv[3*i+1];
        vField[i].component(2) = fieldsToRecv[3*i+2];
    }
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "After receive back state in member " << appSuffix << endl << "\n";

    delete[] fieldsToRecv;
}

/**
 * @brief Receive field for intermediate projection from principal to control members.
 * 
 * @param mesh OpenFOAM mesh.
 * @param vField Velocity field.
 * @param runTime Current runTime.
 * @param cwipiIteration 
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
void cwipiRecvVolVectorFieldMFEnKF(const fvMesh& mesh, volVectorField& vField, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose)
{
    //=== Just as the send equivalent, that's the basic receive of the velocity field, can be automatically 
    // interpolated by cwipi if the coupling meshes are different between the OF instance and the KF_coupling code === 

    // ================================================================== //
    // Receives fields from MFEnKF to OpenFOAM (cwipiMultigridIcoFoamPar) //
    // ================================================================== //
    
    double* fieldsToRecvFine = new double[3 * mesh.nCells()];
    std::cout << "nCells for cwipiRecvVolVectorFieldMFEnKF = " << mesh.nCells() << std::endl;
    // double* fieldsToRecv_2 = new double[3 * mesh.nCells()];

    if (cwipiVerbose) if (Pstream::master()) Pout << "Number of cells in the cwipiRecvVolVectorFieldMFEnKF function: " << mesh.nCells() << endl;

    char cl_receiving_field_name[60];
    // char cl_receiving_field_name_2[60];
    sprintf(cl_receiving_field_name, "U_projected, V_projected, W_projected");
    // sprintf(cl_receiving_field_name_2, "U_fineForecast, V_fineForecast, W_fineForecast");

    scalar t = runTime.value();

    int myGlobalRank;
    label nMyProc=Pstream::myProcNo();

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    int appSuffix = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;  // Normally works for any decomposition
    
    char rankChar[50];
    sprintf(rankChar, "%i", myGlobalRank);
    char couplingName[250] = {"cwipiFoamCoupling_forPrincipalVariateProjection"};

    char appSuffixChar[50];
    sprintf(appSuffixChar, "%i", appSuffix);
    strcat(couplingName, appSuffixChar);

    if (cwipiVerbose) if (Pstream::master()) Pout << "[MFEnKF] Before receive back state projection " << appSuffix << " the coupling name for the receive is " << couplingName << endl;

    cwipi_irecv(couplingName,"ex_proj", recvTag, 3, cwipiIteration, t, cl_receiving_field_name, fieldsToRecvFine, &status200);

    Foam::Info << "cwipi_irecv func from cwipiRecFEnKF function from cwipiPstreamPar.C passed" << Foam::endl;

    Foam::Info << "About to pass through cwipi_wait_irecv func from cwipiRecvVolVectorFieldMFEnKF function from cwipiPstreamPar.C..." << Foam::endl;

    cwipi_wait_irecv(couplingName, status200);

    Foam::Info << "cwipi_wait_irecv func from cwipiRecvVolVectorFieldMFEnKF function from cwipiPstreamPar.C passed" << Foam::endl;

    switch(status200)
    {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose) if (Pstream::master()) Pout << "Back-receive Ok projection in member " << appSuffix << " from processor " << Foam::Pstream::myProcNo() << endl;
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose) Pout << "Bad Back-receiving projection in member " << appSuffix << "\n";
        break;
        default :
        if (cwipiVerbose) Pout << "Error : bad exchange status projection in member " << appSuffix << "\n";
    }

    //* Here the we change the velocity field of the calculation, be careful not to a const volVectorField 
    // in the arguments of the function or it doesn't work *
    
    forAll(mesh.cells(), i)
    {
        vField[i].component(0) = fieldsToRecvFine[3*i+0];
        vField[i].component(1) = fieldsToRecvFine[3*i+1];
        vField[i].component(2) = fieldsToRecvFine[3*i+2];

    }
    
    if (cwipiVerbose) if (Pstream::master()) Pout << "After receive back state projection in member " << appSuffix << endl << "\n";

    // std::cout << "Dimensions of statematricesDiff = (" << stateMatricesDiff.rows() << ", " << stateMatricesDiff.cols() << ")" << std::endl;

    delete[] fieldsToRecvFine;
}

/**
 * @brief Delete coupling.
 * 
 * @param pointCoords Point coordinates.
 * @param face_index Face indices.
 * @param face_connectivity_index Face connectivity indices.
 * @param cell_to_face_connectivity Cell to face connectivity.
 * @param face_connectivity Face connectivity.
 * @param nbParts Number of subdomains.
 * @param cwipiVerbose Verbose.
 */
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
