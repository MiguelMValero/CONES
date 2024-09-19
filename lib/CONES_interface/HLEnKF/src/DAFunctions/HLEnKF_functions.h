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
 * @file HLEnKF_functions.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Data Assimilation function declarations.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <vector>
#include <numeric>  
#include <algorithm>

#include <cwipi.h>

#include "../DAClasses/region/region.h"


//===========================================================================
void tic(int mode=0, std::string procedure="default");

void toc(std::string procedure="default");


//===========================================================================
void define_mesh(double* pointsCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], float cwipiVerbose);

//===========================================================================
Eigen::MatrixXd retrieve_regionState(const Eigen::Ref<const Eigen::MatrixXd>& invStateMatrix, region currentRegion, int nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose);

//===========================================================================
void update_globalState(Eigen::Ref<Eigen::MatrixXd> stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXd>& regionState, region currentRegion, int nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose);

//===========================================================================
Eigen::ArrayXd retrieve_obsVector(region currentRegion, observations obsData, double time, float cwipiVerbose);

//========================== Randomize observations =========================
Eigen::MatrixXd randomize_Obs(const Eigen::Ref<const Eigen::ArrayXd>& obsVector, int cwipiMembers, const observations obsData, const region currentRegion, float cwipiVerbose);

//===========================================================================
Eigen::MatrixXd read_SampData(int cwipiMembers, const observations obsData, float cwipiVerbose, std::string stringRootPath);

//===========================================================================
Eigen::MatrixXd retrieve_regionSampData(const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, int cwipiMembers, observations obsData, region currentRegion, float cwipiVerbose, std::string stringRootPath);

//============= Calculate measurments error covariance matrix ===============
Eigen::MatrixXd calculate_R(const Eigen::Ref<const Eigen::ArrayXd>& obsVector, const region currentRegion, const observations obsData, float cwipiVerbose);

//========================== Calculate Kalman gain ==========================
Eigen::MatrixXd calculate_K(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, const Eigen::Ref<const Eigen::MatrixXd>& obsMatrix, const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, const Eigen::Ref<const Eigen::MatrixXd>& R, int cwipiMembers, region currentRegion, int cwipiParams, float cwipiVerbose);

//===========================================================================
void EnKF_outputs(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, const observations obsData,int cwipiMembers, int nb_cells, double time, int cwipiParams, float cwipiVerbose, double epsilon);

//========================== Simple localization ============================
Eigen::MatrixXd localisation(const Eigen::Ref<const Eigen::MatrixXd>& KnoCorr, int cwipiParams, region currentRegion, double sigmaLocX, double sigmaLocY, double sigmaLocZ, observations obsData, const fvMesh& mesh, float cwipiVerbose);

//======================== State Inflation function =========================
Eigen::MatrixXd stateInflation(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, int cwipiMembers, const region currentRegion, int cwipiParams, double stateInfl, Foam::word typeInfl, float cwipiVerbose);

//====================== Parameters Inflation function ======================
Eigen::MatrixXd paramInflation(const Eigen::Ref<const Eigen::MatrixXd>& parameters, int cwipiMembers, int cwipiParams, double paramsInfl, Foam::word typeInfl, float cwipiVerbose);

//=========================== HLEnKF algorithm ==============================
Eigen::MatrixXd HLEnKF(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, const fvMesh& mesh, int cwipiMembers, int nb_cells, std::vector<region> regionList, observations obsData,double sigmaLocX, double sigmaLocY, double sigmaLocZ, bool localSwitch, int cwipiParams, double stateInfl, double preStateInfl,double paramsInfl, Foam::word typeInfl, bool paramEstSwitch, bool stateEstSwitch, float cwipiVerbose, std::string stringRootPath, double time, double epsilon);

//===========================================================================
void prepare_sendBack(double *sendValues, double *paramsSendValues, double *values, const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, int nb_cells, int cwipiParams, int index, float cwipiVerbose);

//===========================================================================
void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, std::string name);