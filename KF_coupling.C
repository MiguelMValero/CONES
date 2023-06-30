#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include <time.h>
#include <math.h>
#include <libgen.h>       // dirname
#include <unistd.h>       // readlink
#include <linux/limits.h> // PATH_MAX
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

  double *configValues = NULL;
  configValues = (double *)malloc(sizeof(double) * 34);
  configuration(configValues);

  int cwipiStep = configValues[0];         // Number of time step between each DA phase
  int cwipiMembers = configValues[1];      // Number of members in the ensemble
  int subdomains = configValues[2];        // Number of subdomains
  int mainsubdom = configValues[3];        // Main subdomain (where the probes are located)
  int cwipiObsU = configValues[4];         // Number of observation probes for velocity
  int cwipiObsp = configValues[5];         // Number of observation probes for pressure
  int cwipiObsCf = configValues[6];        // Number of observation probes for friction coefficient (default = 1)
  int cwipiParams = configValues[7];       // Number of parameters to optimize with the DA
  double geom_tol = configValues[8];       // Geometric tolerance of the cwipi coupling meshes
  int cwipiOutputNb = configValues[9];     // Number of txt file written beginning from the last one
  double sigmaUserU = configValues[10];    // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for velocity)
  double sigmaUserp = configValues[11];    // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for pressure)
  double sigmaUserCf = configValues[12];   // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for friction coefficient)
  float cwipiVerbose = configValues[13];   // Print all the debuging messages or not: 1 = printed, 0 = nothing
  int cwipiTimedObs = configValues[14];    // Switch to for the obs: 1 = obs depends on time, 0 = obs does not depend on time
  double obsTimeStep = configValues[15];   // The time step of the observations if the case is unsteady
  int cwipiParamsObs = configValues[16];   // Definition of observation parameter : 0 = vel, 1 = pres, 2 = both, 3 = vel+cf
  double stateInfl = configValues[17];     // State inflation (default = 1)
  double paramsInfl = configValues[18];    // Parameters inflation (default = 1)
  float typeInfl = configValues[19];       // Definition of inflation (0 = stochastic, 1 = deterministic)
  float clippingSwitch = configValues[20]; // Switch for the clipping
  float localSwitch = configValues[21];    // Switch for the localisation (clippingSwitch == 1)
  float paramEstSwitch = configValues[22]; // Switch for parameter estimation
  float Ux = configValues[23];             // Specification if Ux is read or not (cwipiParamsObs needs to be either 0, 2 or 3)
  float Uy = configValues[24];             // Specification if Uy is read or not (cwipiParamsObs needs to be either 0, 2 or 3)
  float Uz = configValues[25];             // Specification if Uz is read or not (cwipiParamsObs needs to be either 0, 2 or 3)
  int typeInputs = configValues[26];       // Inputs for R are given in absolute values (0), percentage (1) or potential function (2) (default = 0)
  double sigmaLocX = configValues[27];     // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in X direction)
  double sigmaLocY = configValues[28];     // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in Y direction)
  double sigmaLocZ = configValues[29];     // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in Z direction)
  double epsilon = configValues[30];       // Value used to calculate the NSRMD
  double sigmaUserUa = configValues[31];   // sigma of the EnKF given by a+by^c
  double sigmaUserUb = configValues[32];   // sigma of the EnKF given by a+by^c
  double sigmaUserUc = configValues[33];   // sigma of the EnKF given by a+by^c

  if (cwipiVerbose) std::cout << "Beginning of the configuration file" << std::endl << "\n";
  
  //========== We define the 7 possibilities in case we are dealing with some components of the velocity =========//

  int velocityCase = 0;
  if (cwipiParamsObs != 1)
  {
    if (Ux == 1 && Uy == 0 && Uz == 0)
      velocityCase = 1;
    else if (Ux == 0 && Uy == 1 && Uz == 0)
      velocityCase = 2;
    else if (Ux == 0 && Uy == 0 && Uz == 1)
      velocityCase = 3;
    else if (Ux == 1 && Uy == 1 && Uz == 0)
      velocityCase = 4;
    else if (Ux == 1 && Uy == 0 && Uz == 1)
      velocityCase = 5;
    else if (Ux == 0 && Uy == 1 && Uz == 1)
      velocityCase = 6;
    else if (Ux == 1 && Uy == 1 && Uz == 1)
      velocityCase = 7;
    else
    {
      std::cerr << "Check your inputs. Either cwipiParamsObs or the components of the velocity are not properly defined.\n";
    }
  }

  if (cwipiVerbose)
    std::cout << "Your case is velocityCase " << velocityCase << ".\n";
  int cwipiObs = 0;
  if (cwipiParamsObs == 0)
  {
    switch (velocityCase)
    {
    case 1:
    case 2:
    case 3:
      cwipiObs = cwipiObsU;
      break;
    case 4:
    case 5:
    case 6:
      cwipiObs = 2 * cwipiObsU;
      break;
    case 7:
      cwipiObs = 3 * cwipiObsU;
      break;
    }
  }
  else if (cwipiParamsObs == 1)
    cwipiObs = cwipiObsp;
  else if (cwipiParamsObs == 2)
  {
    switch (velocityCase)
    {
    case 1:
    case 2:
    case 3:
      cwipiObs = cwipiObsU + cwipiObsp;
      break;
    case 4:
    case 5:
    case 6:
      cwipiObs = 2 * cwipiObsU + cwipiObsp;
      break;
    case 7:
      cwipiObs = 3 * cwipiObsU + cwipiObsp;
      break;
    }
  }
  else if (cwipiParamsObs == 3)
  {
    switch (velocityCase)
    {
    case 1:
    case 2:
    case 3:
      cwipiObs = cwipiObsU + cwipiObsCf;
      break;
    case 4:
    case 5:
    case 6:
      cwipiObs = 2 * cwipiObsU + cwipiObsCf;
      break;
    case 7:
      cwipiObs = 3 * cwipiObsU + cwipiObsCf;
      break;
    }
  }

  //========== Declaration of variables ==========

  char cl_coupling_name[20];
  char output_format[20], output_format_option[20];
  char codeName[20] = {"KF"};
  char codeCoupledName[20];
  char recv_field_name[20] = {"U,V,W"};
  char indexChar[20];

  static int recvTag = 1;
  static int recvTag_params = 3;

  static int status;
  static int sendTag = 2;
  static int sendTag_params = 4;
  static int status2;
  MPI_Status status3;

  int grank;
  int rank;
  int commWorldSize;
  // int nNotLocatedPoints = 0;

  //========== MPI Initilization ==========

  MPI_Comm localcomm;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &grank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  //========== CWIPI Initialization ==========

  cwipi_init(MPI_COMM_WORLD, codeName, &localcomm);
  MPI_Comm_rank(localcomm, &rank);

  //========== OpenFOAM environment for mesh definition ==========

  char result[PATH_MAX] = {};
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  const char *path;
  if (count != -1)
  {
    path = dirname(result);
  }

  std::string stringRootPath = path;
  // char RunCase[250]= {"/processor0"};
  char RunCase[250] = "/cavity_testPar1";


  // Foam::Time runTime(Foam::Time::controlDictName, path, RunCase);
  Foam::Time runTime(Foam::Time::controlDictName, path, {RunCase});

  if (cwipiVerbose)
    std::cout << "Runtime created" << std::endl
              << "\n";

  Foam::fvMesh mesh(
      Foam::IOobject(
          Foam::fvMesh::defaultRegion,
          runTime.timeName(),
          runTime,
          Foam::IOobject::MUST_READ));

  int nb_cells = mesh.nCells();
  if (cwipiVerbose)
    std::cout << "The number of cells from EnKF is " << mesh.nCells() << std::endl << "\n";

  double *values = (double *)calloc(sizeof(double), 3 * nb_cells);
  double *sendValues = (double *)calloc(sizeof(double), 3 * nb_cells);
  double *paramsValues = (double *)calloc(sizeof(double), cwipiParams);
  double *paramsSendValues = (double *)calloc(sizeof(double), cwipiParams);

  double *pointCoords = new double[3 * mesh.nPoints()];
  int *face_index = new int[mesh.nCells() + 1];
  int *face_connectivity_index = new int[mesh.nFaces() + 1];

  int c2fconnec_size = 0;
  forAll(mesh.cells(), i)
  {
    c2fconnec_size = c2fconnec_size + mesh.cells()[i].size();
  }
  int *cell_to_face_connectivity = new int[c2fconnec_size];

  int fconnec_size = 0;
  forAll(mesh.faces(), i)
  {
    fconnec_size = fconnec_size + mesh.faces()[i].size();
  }
  int *face_connectivity = new int[fconnec_size];

  //========== Create coupling ==========

  if (cwipiVerbose)
    printf("Create coupling\n");
  cwipi_solver_type_t solver_type;
  solver_type = CWIPI_SOLVER_CELL_CENTER;
  sprintf(output_format, "EnSight Gold");
  sprintf(output_format_option, "text");

  //* Coupling declaration in a loop in order to have a coupling with each OF instance *
  for (int j = 1; j < (cwipiMembers + 1); j++)
  {
    sprintf(cl_coupling_name, "cwipiFoamCoupling");
    sprintf(codeCoupledName, "cwipiFoam");
    sprintf(indexChar, "%i", j);
    strcat(cl_coupling_name, indexChar);
    strcat(codeCoupledName, indexChar);

    if (cwipiVerbose) std::cout<< "From the KF side the name of my coupling is " << 
    cl_coupling_name << std::endl;

    cwipi_create_coupling(cl_coupling_name,
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          codeCoupledName,   // Coupled application id
                          3,                 // Geometric entities dimension
                          geom_tol,          // Geometric tolerance
                          CWIPI_STATIC_MESH, // Mesh type
                          solver_type,       // Solver type
                          0,                 // Postprocessing frequency
                          output_format,     // Postprocessing format
                          output_format_option);

    cwipi_synchronize_control_parameter(codeCoupledName);

    if (cwipiVerbose) std::cout<< "Coupling called " << cl_coupling_name << " correctly created " << std::endl;

    //========== Mesh definition ============

    if (cwipiVerbose)
      printf("%s: Create mesh \n", __func__);

    define_mesh(pointCoords, face_index, cell_to_face_connectivity, face_connectivity_index, face_connectivity, c2fconnec_size, fconnec_size, mesh, cl_coupling_name, cwipiVerbose);
  }

  //========== Receive values ==========

  //* Receive the values of the control parameters from the first OF instance *
  double time = cwipi_get_distant_double_control_parameter("cwipiFoam1", "currentTime");
  int numberCwipiPhase = cwipi_get_distant_int_control_parameter("cwipiFoam1", "numberCwipiPhase");
  double deltaT = cwipi_get_distant_double_control_parameter("cwipiFoam1", "deltaT");
  // int nbParts = cwipi_get_distant_int_control_parameter("cwipiFoam1", "nbParts");
  // int partsReparty = cwipi_get_distant_int_control_parameter("cwipiFoam1","partsReparty");

  int firstCwipiPhase = (time / deltaT) / cwipiStep;
  int inv_nb_cells = nb_cells;

  if (cwipiVerbose)
    std::cout << "numberCwipiPhase: " << numberCwipiPhase << "\n";

  //* We receive and send back on a loop for the number of exchanges througout calculation, and another loop
  // for all the members (Of instances) in the ensemble *

  for (int i = firstCwipiPhase; i < numberCwipiPhase; i++)
  {
    configuration(configValues);

    cwipiVerbose = configValues[13]; // Print all the debuging messages or not, 1 printed 0 nothing
    stateInfl = configValues[17];    // Standard deviation for state inflation
    paramsInfl = configValues[18];   // Standard deviation for parameters inflation
    typeInfl = configValues[19];     // Definition of inflation (0 = stochastic, 1 = deterministic)

    time = time + cwipiStep * deltaT;

    MatrixXf stateVector = MatrixXf::Zero(nb_cells * 3 + cwipiParams, cwipiMembers);
    MatrixXf invStateVector;
    if (clippingSwitch)
    {
      invStateVector = MatrixXf::Zero(nb_cells * 3 + cwipiParams, cwipiMembers);
    }

    if (cwipiVerbose)
      std::cout << "Phase " << i << " going from " << firstCwipiPhase << " to " << numberCwipiPhase - 1 << std::endl
                << "\n";

    for (int j = 1; j < (cwipiMembers + 1); j++)
    {
      sprintf(cl_coupling_name, "cwipiFoamCoupling");
      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name, indexChar);

      //**** Receive velocity field ****
      if (cwipiVerbose == 1)
        std::cout << "Before receive in ensemble " << j << "\n";
      cwipi_irecv(cl_coupling_name, "ex1", recvTag, 3, i + 1, time, recv_field_name, values, &status);
      cwipi_wait_irecv(cl_coupling_name, status);
      if (cwipiVerbose == 1)
        std::cout << "After receive in ensemble " << j << "\n";

      switch (status)
      {
      case CWIPI_EXCHANGE_OK:
        if (cwipiVerbose)
          std::cout << "Receive OK in ensemble " << j << "\n";
        break;
      case CWIPI_EXCHANGE_BAD_RECEIVING:
        if (cwipiVerbose)
          std::cout << "Bad receiving in ensemble " << j << "\n";
        break;
      default:
        if (cwipiVerbose)
          std::cout << "Error: bad exchange status in ensemble " << j << "\n";
      }

      //**** Receive parameters ****
      if (cwipiVerbose)
        std::cout << "Before receive the parameters in ensemble " << j << "\n";
      
      //========= For the cavity testcase, the rank from where the parameters come is j*subdomains  + (mainsubDomain - 1) (we need
      // to find a generic formula for this) =========//
      int coming_rank = j*subdomains + (mainsubdom - 1);
      MPI_Recv(paramsValues, cwipiParams, MPI_DOUBLE, coming_rank, recvTag_params, MPI_COMM_WORLD, &status3);

      if (cwipiVerbose){
        std::cout << "After receive the parameters in ensemble " << j << "\n";
        std::cout << "==== Parameters well received ==== Parameter 1 = " << paramsValues[0] << std::endl << "\n";
        std::cout << "The number of cells in the side of EnKF is " << nb_cells << std::endl;
      }

      for (int k = 0; k < nb_cells; k++)
      {
        stateVector(k, j-1) = values[3*k + 0];
        stateVector(k + nb_cells, j-1) = values[3*k + 1];
        stateVector(k + 2*nb_cells, j-1) = values[3*k + 2];
      }

      //* The parameters are added at the end of the state vector *
      for (int k = 0; k < cwipiParams; k++)
      {
        stateVector(3*nb_cells + k, j-1) = paramsValues[k];
      }
    }

    //====== Writing values to a txt file =========
    print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateVector, "UMat");

    if (clippingSwitch == 1)
    {
      stateVector = doClipping(stateVector, invStateVector, nb_cells, cwipiParams, cwipiMembers, cwipiVerbose, stringRootPath);

      //====== Writing values to a txt file =========
      print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateVector, "UMat_Clip");
    }
    //======== Kalman filter code =======

    //========= The observation Matrix will be different depending on the high-fidelity observations
    ArrayXf obsMatrix(cwipiObs); // By default our parameter is the velocity

    if (cwipiTimedObs)
    {
      int obsIndex = time / obsTimeStep - 1;
      obsMatrix = obs_Data_timed(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, obsIndex, cwipiParamsObs, cwipiVerbose, stringRootPath);
    }
    else
      obsMatrix = obs_Data(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, cwipiParamsObs, cwipiVerbose, stringRootPath);

    MatrixXf sampMatrix = samp_Data(cwipiMembers, cwipiObs, cwipiObsU, velocityCase, cwipiParamsObs, cwipiVerbose, stringRootPath);

    MatrixXf UptMatrix = KF(stateVector, obsMatrix, sampMatrix, cwipiMembers, nb_cells, cwipiObs, cwipiObsU, sigmaUserU, sigmaUserp, sigmaUserCf,
                            sigmaLocX, sigmaLocY, sigmaLocZ, localSwitch, clippingSwitch, cwipiParams, cwipiParamsObs, stateInfl, paramsInfl, typeInfl,
                            typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, paramEstSwitch, mesh, cwipiVerbose, stringRootPath);

    if (clippingSwitch)
    {

      //====== Writing updated values to a txt file =========
      print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, UptMatrix, "UMat_Clip_upt");
      UptMatrix = undoClipping(UptMatrix, invStateVector, inv_nb_cells, nb_cells, cwipiParams, cwipiMembers, cwipiVerbose, stringRootPath);
    }

    print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, UptMatrix, "UMat_upt");

    //================== Send back ======================

    for (int j = 1; j < cwipiMembers + 1; j++)
    {
      KF_output(sendValues, paramsSendValues, values, UptMatrix, sampMatrix, obsMatrix, cwipiMembers, nb_cells, time, cwipiParams, cwipiObsU, cwipiObs, cwipiParamsObs, velocityCase,
                j, subdomains, mainsubdom, epsilon, cwipiVerbose);
      sprintf(cl_coupling_name, "cwipiFoamCoupling");
      sprintf(indexChar, "%i", j);
      strcat(cl_coupling_name, indexChar);

      if (cwipiVerbose)
        std::cout << "Before re-receive the parameters in ensemble " << j << "\n";

      cwipi_issend(cl_coupling_name, "ex2", sendTag, 3, i + 1, time, recv_field_name, sendValues, &status2);
      cwipi_wait_issend(cl_coupling_name, status2);

      switch (status2)
      {
        case CWIPI_EXCHANGE_OK:
          if (cwipiVerbose)
            std::cout << "Re-sent Ok in ensemble " << j << "\n";
          break;
        case CWIPI_EXCHANGE_BAD_RECEIVING:
          if (cwipiVerbose)
            std::cout << "Bad re-sent in ensemble " << j << "\n";
          break;
        default:
          if (cwipiVerbose)
            std::cout << "Error: bad re-sent status in ensemble " << j << "\n";
      }
      // if (cwipiVerbose)
      //   printf("After re-send \n");

      int going_rank = j*subdomains + (mainsubdom - 1);
      MPI_Send(paramsSendValues, cwipiParams, MPI_DOUBLE, going_rank, sendTag_params, MPI_COMM_WORLD);
    }
  }

  //========= Delete coupling ==========

  if (rank == 0)
    printf("Delete Cwipi coupling from c++ file\n");

  for (int j = 1; j < (cwipiMembers*subdomains) + 1; j++)
  {
    sprintf(cl_coupling_name, "cwipiFoamCoupling");
    sprintf(indexChar, "%i", j);
    strcat(cl_coupling_name, indexChar);

    cwipi_delete_coupling(cl_coupling_name);
  }

  if (cwipiVerbose)
    std::cout << "After cwipi delete coupling from c++ file" << std::endl;

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
  free(configValues);

  if (cwipiVerbose)
    std::cout << "Arrays deallocated in KF" << std::endl;

  //========== Finalize ==========

  cwipi_finalize();

  if (cwipiVerbose)
    std::cout << "Cwipi environment finalized from KF " << std::endl;

  MPI_Finalize();

  if (cwipiVerbose)
    std::cout << "MPI environment finalized from KF" << std::endl;

  return EXIT_SUCCESS;
}
