#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include <time.h>
#include <math.h>
#include <libgen.h>         // dirname
#include <unistd.h>         // readlink
#include <linux/limits.h>   // PATH_MAX
#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <cwipi.h>
#include "KF_coupling_functions.h"

using namespace Eigen;

/*----------------------------------------------------------------------
 *                                                                     
 * Main : Couling code for EnKF
 *
 *---------------------------------------------------------------------*/
 
int main(int argc, char *argv[])
{
    //========== Parameters from config file ==========

    std::cout << "Beginning of the file" << std::endl << "\n";

    std::ifstream cFile ("cwipiConfig");

    double configValues[12]={0};
    int k = 0;

    if (cFile.is_open())
    {
        std::string line;
        while(std::getline(cFile, line))
        {
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());
            if( line.empty() || line[0] == '#' )
            {
                continue;
            }
            auto delimiterPos = line.find("=");
            // auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);

            configValues[k] = stod(value);
            // std::cout << values[k] << '\n';
            k = k+1;
        }
    }
    else 
    {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    int cwipiStep = configValues[0];       // Number of time step between each DA phase
    int cwipiMembers = configValues[1];    // Number of members in the ensemble
    int cwipiObsU = configValues[2];        // Number of observation probes for velocity
    int cwipiObsp = configValues[3];        // Number of observation probes for pressure
    int cwipiParams = configValues[4];     // Number of parameters to optimize with the DA
    double geom_tol = configValues[5];     // Geometric tolerance of the cwipi coupling meshes
    int cwipiOutputNb = configValues[6];   // Number of txt file written beginning from the last one
    double sigmaUserU = configValues[7];    // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for velocity)
    double sigmaUserp = configValues[8];    // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for pressure)
    int cwipiVerbose = configValues[9];    // Print all the debuging messages or not, 1 printed 0 nothing
    int cwipiTimedObs = configValues[10];   // Switch to for the obs : 1 = obs depends on time
    int cwipiParamsObs = configValues[11];   // Definition of observation parameter : 0 = vel, 1 = pres, 2 = both

    int cwipiObs = 0;
    if (cwipiParamsObs == 0) cwipiObs = 3*cwipiObsU;
    else if (cwipiParamsObs == 1) cwipiObs = cwipiObsp;
    else if (cwipiParamsObs == 2) cwipiObs = 3*cwipiObsU + cwipiObsp;

    if (cwipiVerbose == 1) std::cout << "End of config" << cwipiTimedObs << std::endl << "\n";

    //========== OpenFOAM environment for mesh definition ==========

    char result[PATH_MAX]={};
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    const char *path;
    if (count != -1) {
        path = dirname(result);
    }
    
    std::string stringRootPath = path;

    if (cwipiVerbose == 1) std::cout << "stringRootPath = " << stringRootPath << std::endl << "\n";

    // char casePath[250] = {path};
    //char RunCase[250] = "";

    char RunCase[250] = {"/processor0"};
      //-- Declaration des variables openFoam --
    Foam::Time runTime(Foam::Time::controlDictName, path, RunCase);
    
    if (cwipiVerbose == 1) std::cout << "Runtime created" << std::endl << "\n";
    
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    int nb_cells = mesh.nCells();
    if (cwipiVerbose == 1) std::cout << "nombre de cellules côté KF " << mesh.nCells() << std::endl << "\n";


    //========== Declaration of variables ==========

    char cl_sending_field_name[20], cl_receiving_field_name[20];
    char cl_coupling_name[20], cl_exchange_name[20];
    char output_format[20], output_format_option[20];
    char codeName[20]={"KF"}; 
    char codeCoupledName[20];
    char recv_field_name[20]={"U,V,W"};
    char recv_field_name_params[20]={"params"};
    char indexChar[20];
    int il_error = 0;

    static int recvTag = 1;
    static int recvTag_params = 3;
    static int stride;

    static int status;
    static int sendTag = 2;
    static int sendTag_params = 4;
    static int status2;
    MPI_Status status3;

    int grank;
    int rank;
    int commWorldSize;
    int coord_id;
    int nNotLocatedPoints = 0;

    double *values = NULL;
    values = (double*) malloc(sizeof(double) * 3*nb_cells);

    double *sendValues = NULL;
    sendValues = (double*) malloc(sizeof(double) * 3*nb_cells);

    double *paramsValues = NULL;
    paramsValues = (double*) malloc(sizeof(double) * cwipiParams);

    double *paramsSendValues = NULL;
    paramsSendValues = (double*) malloc(sizeof(double) * cwipiParams);

    double* pointCoords = new double[3*mesh.nPoints()];
    // int* connecIdx = new int[mesh.nCells()+1];
    // int* connec = new int[mesh.nCells()*8];
    int* face_index = new int[mesh.nCells()+1];
    int* face_connectivity_index = new int[mesh.nFaces()+1];

    int c2fconnec_size = 0;
    forAll(mesh.cells(),i){
      c2fconnec_size = c2fconnec_size + mesh.cells()[i].size();
    }
    int* cell_to_face_connectivity = new int[c2fconnec_size];

    int fconnec_size = 0;
    forAll(mesh.faces(),i){
      fconnec_size = fconnec_size + mesh.faces()[i].size();
    }
    int* face_connectivity = new int[fconnec_size];

    MatrixXf stateVector = MatrixXf::Zero(nb_cells*3+cwipiParams,cwipiMembers);

    //========== MPI Initilization ==========

    MPI_Comm localcomm;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);


    //========== CWIPI Initialization ==========

    cwipi_init(MPI_COMM_WORLD, codeName, &localcomm);

    MPI_Comm_rank(localcomm, &rank);


    //========== Create coupling ==========

    printf("Create coupling\n");
    cwipi_solver_type_t solver_type;
    solver_type = CWIPI_SOLVER_CELL_CENTER;
    sprintf(output_format,"EnSight Gold");
    sprintf(output_format_option,"text");


    //* Coupling declaration in a loop in order to have a coupling with each OF instance *
    for (int j = 1; j < cwipiMembers+1; j++)
    {
      sprintf(cl_coupling_name,"cwipiFoamCoupling");
      sprintf(codeCoupledName,"cwipiFoam");
      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name,indexChar);
      strcat(codeCoupledName,indexChar);

      cwipi_create_coupling(cl_coupling_name,
                            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                            codeCoupledName,                           // Coupled application id
                            3,                                         // Geometric entities dimension
                            geom_tol,                                  // Geometric tolerance
                            CWIPI_STATIC_MESH,                         // Mesh type
                            solver_type,                               // Solver type
                            0,                                         // Postprocessing frequency
                            output_format,                            // Postprocessing format
                            output_format_option);

      cwipi_synchronize_control_parameter(codeCoupledName);   


      //========== Mesh definition ============

      if (rank == 0) printf("Create mesh\n");

      define_mesh(pointCoords, face_index, cell_to_face_connectivity, face_connectivity_index, face_connectivity, c2fconnec_size, fconnec_size, mesh, cl_coupling_name, cwipiVerbose);
    }


    //========== Receive values ==========  

    //* Receive the values of the control parameters from the first OF instance *
    double time = cwipi_get_distant_double_control_parameter("cwipiFoam1","currentTime");
    int numberCwipiPhase = cwipi_get_distant_int_control_parameter("cwipiFoam1","numberCwipiPhase");
    double deltaT = cwipi_get_distant_double_control_parameter("cwipiFoam1","deltaT");
    int nbParts = cwipi_get_distant_int_control_parameter("cwipiFoam1", "nbParts");
    int partsReparty = cwipi_get_distant_int_control_parameter("cwipiFoam1","partsReparty");

    if (cwipiVerbose == 1) std::cout << "numberCwipiPhase: " << numberCwipiPhase << "\n";

    //* We receive and send back on a loop for the number of exchanges througout calculation, and another loop 
    // for all the members (Of instances) in the ensemble *

    for(int i = 0 ; i < numberCwipiPhase; i++)
    {
      time = time + cwipiStep*deltaT;
      for (int j = 1; j < cwipiMembers+1; j++)
      {
        sprintf(cl_coupling_name,"cwipiFoamCoupling");
        sprintf(indexChar, "%i", j);
        strcat(cl_coupling_name,indexChar);

        //**** Receive velocity field ****
        if (cwipiVerbose == 1) printf("Before receive \n");
        cwipi_irecv(cl_coupling_name, "ex1", recvTag, 3, i+1, time, recv_field_name, values, &status);
        cwipi_wait_irecv(cl_coupling_name,status);
        if (cwipiVerbose == 1) printf("After receive \n");

        switch(status)
        {
          case CWIPI_EXCHANGE_OK :
          if (cwipiVerbose == 1) printf("Receive Ok\n");
          break;
          case CWIPI_EXCHANGE_BAD_RECEIVING :
          if (cwipiVerbose == 1) printf("Bad receiving\n");
          break;
          default :
          if (cwipiVerbose == 1) printf("Error : bad exchange status\n");
        }

        //**** Receive parameters ****
        if (cwipiVerbose == 1) printf("Before receive params \n");
        MPI_Recv(paramsValues, cwipiParams, MPI_DOUBLE, j, recvTag_params, MPI_COMM_WORLD, &status3);
        if (cwipiVerbose == 1) printf("After receive params \n");

        if (cwipiVerbose == 1) std::cout << "==== Les parametres sont bien reçus ==== parametre 1 = " << paramsValues[0] << std::endl << "\n";

        for (int k = 0; k < nb_cells; k++)
        {
          stateVector(k,j-1) = values[3*k + 0];
          stateVector(k+nb_cells,j-1) = values[3*k + 1];
          stateVector(k+2*nb_cells,j-1) = values[3*k + 2];
        }
        
        //* The parameters are added at the end of the state vector *
        for (int k = 0; k < cwipiParams; k++)
        {
          stateVector(3*nb_cells+k,j-1) = paramsValues[k];
        }
      }


      //====== Writing values to a txt file =========
      print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateVector, "UMat");


      //======== Kalman filter code =======

      //* The observation Matrix will be different depending on the high-fidelity data parameters
      ArrayXf obsMatrix(cwipiObs); // By default our parameter is the velocity

      if (cwipiTimedObs == 1) obsMatrix = obs_Data_timed(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, i, cwipiParamsObs, cwipiVerbose, stringRootPath);
      else if (cwipiTimedObs == 0) obsMatrix = obs_Data(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, cwipiParamsObs, cwipiVerbose, stringRootPath);
     
      MatrixXf sampMatrix = samp_Data(cwipiMembers, cwipiObs, cwipiObsU, cwipiParamsObs, cwipiVerbose, stringRootPath);

      MatrixXf UptMatrix = KF(stateVector, obsMatrix, sampMatrix, cwipiMembers, nb_cells, cwipiObs, cwipiObsU, sigmaUserU, sigmaUserp, cwipiParams, cwipiParamsObs, cwipiVerbose);


      //====== Writing updated values to a txt file =========
      print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, UptMatrix, "UMat_upt");


      //================== Send back ======================

      for (int j = 1; j < cwipiMembers+1; j++)
      {
        KF_output(sendValues, paramsSendValues, UptMatrix, cwipiMembers, nb_cells, time, cwipiParams, j, cwipiVerbose);

        sprintf(cl_coupling_name,"cwipiFoamCoupling");
        sprintf(indexChar, "%i", j);
        strcat(cl_coupling_name,indexChar);

        if (cwipiVerbose == 1) printf("Before re-send \n");
        cwipi_issend(cl_coupling_name, "ex2", sendTag, 3, i+1, time, recv_field_name, sendValues, &status2);
        cwipi_wait_issend(cl_coupling_name,status2);

        switch(status2)
        {
        case CWIPI_EXCHANGE_OK :
        if (cwipiVerbose == 1) printf("Re-sent Ok\n");
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        if (cwipiVerbose == 1) printf("Bad receiving\n");
        break;
        default :
        if (cwipiVerbose == 1) printf("Error : bad exchange status\n");
        }
        if (cwipiVerbose == 1) printf("After re-send \n");

        MPI_Send(paramsSendValues, cwipiParams, MPI_DOUBLE, j, sendTag_params, MPI_COMM_WORLD);
      }
    }

    //========= Delete coupling ==========

    if (rank == 0) printf("Delete Cwipi coupling from c++ file\n");

    for (int j = 1; j < cwipiMembers+1; j++)
    {
      sprintf(cl_coupling_name,"cwipiFoamCoupling");
      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name,indexChar);

      cwipi_delete_coupling(cl_coupling_name);
    }

    if (cwipiVerbose == 1) printf("After cwipi delete coupling from c++ file");

    //========== Freeing memory ==========

    delete[] pointCoords;
    // delete[] connecIdx;
    // delete[] connec;
    delete[] face_index;
    delete[] cell_to_face_connectivity;
    delete[] face_connectivity_index;
    delete[] face_connectivity;

    free(sendValues);
    free(values);
    free(paramsValues);
    free(paramsSendValues);

    //========== Finalize ==========

    cwipi_finalize();
    MPI_Finalize();

    return EXIT_SUCCESS;
}
