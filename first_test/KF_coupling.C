#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <cwipi.h>
#include "KF_coupling_functions.h"

using namespace Eigen;

/*----------------------------------------------------------------------
 *                                                                     
 * Main : surface coupling test :
 *
 *---------------------------------------------------------------------*/
 
int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{
    char casePath[250] = {"/home/villanul/OpenFOAM/villanul-8/run/cwipi_tests"};
    //char RunCase[250] = "";

    char RunCase[250] = {"first_test"};
      //-- Declaration des variables openFoam --
    Foam::Time runTime(Foam::Time::controlDictName, casePath, RunCase);
    
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
    std::cout << "nombre de cellules côté KF " << mesh.nCells() << std::endl << "\n";

    char cl_sending_field_name[20], cl_receiving_field_name[20];
    char cl_coupling_name[20], cl_exchange_name[20];
    char output_format[20], output_format_option[20];
    int il_error = 0;

    int grank;
    int rank;
    int commWorldSize;
    int coord_id;
    int nNotLocatedPoints = 0;

    float geom_tol = 1.0;

    double *values = NULL;
    values = (double*) malloc(sizeof(double) * 3*nb_cells);

    double *sendValues = NULL;
    sendValues = (double*) malloc(sizeof(double) * 3*nb_cells);

    double* pointCoords = new double[3*mesh.nPoints()];
    int* connecIdx = new int[mesh.nCells()+1];
    int* connec = new int[mesh.nCells()*8];

    //========== MPI Initilization ==========
    printf("cav here 1 \n");

    MPI_Comm localcomm;
    char *codeName;
    int codeId;
    char *codeCoupledName;
    MPI_Init(&argc, &argv);

    printf("cav here 2 \n");

    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

    printf("cav here 3 \n");

    codeName = "FOAM_APE";
    codeCoupledName = "cwipiFoam";

    printf("cav here 4 \n");

    //========== Initialization of the coupling ==========
    cwipi_init(MPI_COMM_WORLD,
              codeName,
              &localcomm);

    printf("cav here 4.5 \n");
    MPI_Comm_rank(localcomm, &rank);

    //========== Create coupling ==========
    printf("Create coupling\n");
    cwipi_solver_type_t solver_type;
    solver_type = CWIPI_SOLVER_CELL_CENTER;
    sprintf(cl_coupling_name,"cwipiFoamCoupling");
    sprintf(output_format,"EnSight Gold");
    sprintf(output_format_option,"text");

    printf("cav here 5 \n");

    // cwipi_add_local_int_control_parameter("receiveTag",1);

    cwipi_create_coupling(cl_coupling_name,
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          codeCoupledName,                           // Coupled application id
                          3,                                         // Geometric entities dimension
                          geom_tol,                                  // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          solver_type,                               // Solver type
                          1,                                         // Postprocessing frequency
                          output_format,                            // Postprocessing format
                          output_format_option);

    cwipi_synchronize_control_parameter("cwipiFoam");
                          
    printf("cav here 6 \n");
      
    //========== Mesh definition ============
    if (rank == 0) printf("Create mesh\n");

    define_mesh(pointCoords,
                connecIdx,
                connec,
                mesh);

    printf("cav here 7 \n");

    //========== Receive values ==========  
    static int recvTag = 1;
    static int recvTag_params = 3;
    static int stride;

    char recv_field_name[20];
    char recv_field_name_params[20];
    // const char *recv_field_name;
    // recv_field_name = cwipi_get_distant_string_control_parameter("cwipiFoam", "sendFieldNames");
    sprintf(recv_field_name,"U,V,W");
    sprintf(recv_field_name_params,"params");

    static int status;
    static int sendTag = 2;
    static int sendTag_params = 4;
    static int status2;
    MPI_Status status3;

    double time = cwipi_get_distant_double_control_parameter("cwipiFoam","currentTime");
    int numberCwipiPhase = cwipi_get_distant_int_control_parameter("cwipiFoam","numberCwipiPhase");
    int cwipiStep = cwipi_get_distant_int_control_parameter("cwipiFoam","cwipiStep");
    double deltaT = cwipi_get_distant_double_control_parameter("cwipiFoam","deltaT");
    int cwipiMembers = cwipi_get_distant_int_control_parameter("cwipiFoam", "cwipiMembers");
    int cwipiObs = cwipi_get_distant_int_control_parameter("cwipiFoam", "cwipiObs");
    int cwipiParams = cwipi_get_distant_int_control_parameter("cwipiFoam","cwipiParams");
    int nbParts = cwipi_get_distant_int_control_parameter("cwipiFoam", "nbParts");
    int partsReparty = cwipi_get_distant_int_control_parameter("cwipiFoam","partsReparty");

    double *paramsValues = NULL;
    paramsValues = (double*) malloc(sizeof(double) * cwipiParams);

    double *paramsSendValues = NULL;
    paramsSendValues = (double*) malloc(sizeof(double) * cwipiParams);

    std::cout << "numberCwipiPhase: " << numberCwipiPhase << "\n";
    for(int i = 0 ; i < numberCwipiPhase; i++)
    {
      //**** Receive velocity field ****
      printf("Before receive \n");
      time = time + cwipiStep*deltaT;
      cwipi_irecv(cl_coupling_name, "ex1", recvTag, 3, i+1, time, recv_field_name, values, &status);
      cwipi_wait_irecv(cl_coupling_name,status);
      printf("After receive \n");

      switch(status)
      {
        case CWIPI_EXCHANGE_OK :
        printf("Receive Ok\n");
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        printf("Bad receiving\n");
        break;
        default :
        printf("Error : bad exchange status\n");
      }

      // //**** Receive parameters ****
      printf("Before receive params \n");
      MPI_Recv(paramsValues,1,MPI_DOUBLE,nbParts,recvTag_params,MPI_COMM_WORLD, &status3);
      printf("After receive params \n");

      std::cout << "==== Le paramètre est bien reçu ==== " << paramsValues[0] << std::endl << "\n";

      //====== Writing values to a txt file =========
      std::string stepNumber = std::to_string(time);
      std::string filename = "U"+stepNumber;
      std::fstream file_out;
      
      file_out.open(filename, std::ios_base::out);
      
      if (!file_out.is_open())
      {
          std::cout << "failed to open " << filename << '\n';
      } 
      else 
      {
        for(int i=0;i<nb_cells;i++)
        {
          file_out << values[3*i+0] << ' ' << values[3*i+1] << ' ' << values[3*i+2] << ' ' << std::endl;
        }
        file_out << "\n";
      }
      std::cout << "Done Writing Updated Matrix on txt file!" << std::endl << "\n";

      //======== Kalman filter code =======

      MatrixXf stateVector = repro_Members(cwipiMembers, nb_cells, time, values, paramsValues, cwipiParams);

      ArrayXf obsMatrix = obs_Data(cwipiMembers, nb_cells, cwipiObs, time);
     
      MatrixXf sampMatrix = samp_Data(cwipiMembers , cwipiObs);

      MatrixXf UptMatrix = KF(stateVector, obsMatrix, sampMatrix, cwipiMembers, nb_cells, cwipiObs, 0.1, cwipiParams);

      KF_output(sendValues, paramsSendValues, stateVector, cwipiMembers, nb_cells, cwipiObs, time, cwipiParams);

      //================== Send back ======================
      printf("Before re-send \n");
      cwipi_issend(cl_coupling_name, "ex2", sendTag, 3, i+1, time, recv_field_name, sendValues, &status2);
      cwipi_wait_issend(cl_coupling_name,status2);

      switch(status2)
      {
      case CWIPI_EXCHANGE_OK :
      printf("Re-sent Ok\n");
      break;
      case CWIPI_EXCHANGE_BAD_RECEIVING :
      printf("Bad receiving\n");
      break;
      default :
      printf("Error : bad exchange status\n");
      }
      printf("After re-send \n");

      int testRank = nbParts/partsReparty + nbParts/partsReparty*(partsReparty-2);

      for(int i = testRank+1; i < nbParts+1; i++)
      {
        MPI_Send(paramsSendValues,1,MPI_DOUBLE,i,sendTag_params,MPI_COMM_WORLD);
      }
    }

    //========= Delete coupling ==========
    if (rank == 0) printf("Delete Cwipi coupling from c++ file\n");

    cwipi_delete_coupling(cl_coupling_name);

    printf("After cwipi delete coupling from c++ file");

    //========== Freeing memory ==========
    delete[] pointCoords;
    delete[] connecIdx;
    delete[] connec;
    free(sendValues);
    free(values);
    free(paramsValues);
    free(paramsSendValues);

    //========== Finalize ==========
    cwipi_finalize();
    MPI_Finalize();

    return EXIT_SUCCESS;
}
