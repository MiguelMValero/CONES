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
 * @dir ./lib/CONES_interface/HLEnKF
 * @brief **Hyper-localisation.** 
 * 
 */
/**
 * @file HLEnKF_main.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Main HLEnKF file.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "fvCFD.H"
#include <mpi.h>

#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <cwipi.h>
#include "src/DAFunctions/HLEnKF_functions.h"
#include "src/DAClasses/region/region.h"
#include "src/DAClasses/observations/observations.h"

using namespace Eigen;

/*----------------------------------------------------------------------
 *
 * Main : Coupling code for EnKF
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  //========== Initialization of the path, Foam environment and configuration of the hyperparameters ==========
  #include "init_Path_Foam.h"
  #include "init_configuration.h"  

  //========== Declaration of the cwipi coupling =======
  #include "coupling_declaration.h"

  //========== Receive values ==========

  //* Receive the values of the control parameters from the first OF instance *
  double time = cwipi_get_distant_double_control_parameter("cwipiFoam1", "currentTime");
  int numberCwipiPhase = cwipi_get_distant_int_control_parameter("cwipiFoam1", "numberCwipiPhase");
  double deltaT = cwipi_get_distant_double_control_parameter("cwipiFoam1", "deltaT");

  int firstCwipiPhase = round(time/deltaT)/cwipiStep;

  if (cwipiVerbose) std::cout << "numberCwipiPhase: " << numberCwipiPhase << "\n";

  //* We receive and send back on a loop for the number of exchanges througout calculation, and another loop
  // for all the members (Of instances) in the ensemble *

  // Testing implementation of observations class
  observations obsData(cwipiObsU, cwipiObsp, cwipiObsCf, obsVelocityComponents, sigmaUserU, sigmaUserp, sigmaUserCf, typeInputs, cwipiTimedObs, obsTimeStep, obsStartTime, obsType, stringRootPath, obsCoordFile, obsDataFile);
  // obsData.info_obsCoord();
  // obsData.info_obsData();
  // std::cout << "First probe degrees of freedom : " << obsData.get_probeVar(0) << "\n";
  // std::cout << "First probe components : " << obsData.get_probeVarLabels(0) << "\n";
  // std::cout << "First probe coordinates: " << obsData.get_probeCoord(0)[0] << " " << obsData.get_probeCoord(0)[1] << " " << obsData.get_probeCoord(0)[2] << "\n";
  // std::cout << "First probe data first element : " << obsData.get_probeData(0)[0][0] << " and second element " << obsData.get_probeData(0)[1][0] << std::endl << "\n";

  // Testing implementation of region class 
  std::vector<region> regionList = initialize_DA_regions(stringRootPath, clippingFile, cwipiVerbose, clippingSwitch, nb_cells, obsData);


  // end of observation class test
  for (int i = firstCwipiPhase; i < numberCwipiPhase; i++)
  {
    // ========== Reload configurations if needed ===========
    // configuration(configValues);
    // cwipiVerbose = configValues[12]; // Print all the debuging messages or not, 1 printed 0 nothing
    // stateInfl = configValues[16];    // Standard deviation for state inflation
    // paramsInfl = configValues[17];   // Standard deviation for parameters inflation
    // typeInfl = configValues[18];     // Definition of inflation (0 = stochastic, 1 = deterministic)

    time = time + cwipiStep * deltaT;

    if (cwipiVerbose)
      std::cout << "Phase " << i << " going from " << firstCwipiPhase << " to " << numberCwipiPhase - 1 << std::endl
                << "\n";

    MatrixXd stateMatrix = MatrixXd::Zero(nb_cells * 3 + cwipiParams, cwipiMembers);
    MatrixXd stateMatrixUpt = MatrixXd::Zero(nb_cells * 3 + cwipiParams, cwipiMembers);
    
    // =================== Receive procedure ======================
    #include "receive_procedure.h"

    //** Writing values to a txt file if needed **
    print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateMatrix, "UMat");

    //==================== Kalman filter code =====================

    // tic();
    stateMatrixUpt = HLEnKF(stateMatrix, mesh, cwipiMembers, nb_cells, regionList, obsData, sigmaLocX, sigmaLocY, sigmaLocZ, localSwitch, cwipiParams, stateInfl, preStateInfl, paramsInfl, typeInfl, paramEstSwitch, stateEstSwitch, cwipiVerbose, stringRootPath, time, epsilon);
    // toc("HLEnKF analysis (from reading to outputs)");

    //** Writing values to a txt file if needed **
    print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateMatrixUpt, "UMat_upt");

    //===================== Send back procedure ========================
    #include "send_procedure.h"
  }

  //===================== Finalization procedure =======================
  #include "finalize_calculation.h"

  return EXIT_SUCCESS;
}
