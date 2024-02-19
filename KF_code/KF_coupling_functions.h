#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <numeric>  
#include <algorithm>
#include <Eigen/Dense>

#include <cwipi.h>

//===========================================================================
void configuration(double* configValues);

//===========================================================================
Eigen::MatrixXf doClipping(const Eigen::Ref<const Eigen::MatrixXf>& invStateVector, int& nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose, std::string stringRootPath);

//===========================================================================
Eigen::MatrixXf undoClipping(Eigen::MatrixXf stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXf>& invStateVector, int inv_nb_cells, int& nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose, std::string stringRootPath); 

//===========================================================================
void define_mesh(double* pointsCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], float cwipiVerbose);

//===========================================================================
Eigen::ArrayXf obs_Data(int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

//===========================================================================
Eigen::ArrayXf obs_Data_timed(int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, int obsIndex, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

//===========================================================================
Eigen::MatrixXf samp_Data(int cwipiMembers, int cwipiObs, int cwipiObsU, int velocityCase, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

//===========================================================================
void EnKF_outputs(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiMembers, int nb_cells, double time, int cwipiParams, int cwipiObsU, int cwipiObs, int cwipiParamsObs, int velocityCase, int index, float cwipiVerbose, double epsilon);

//========================== Simple localization ============================
Eigen::MatrixXf localisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, int nb_cells, int cwipiParams, int cwipiObs, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float clippingSwitch, float hyperlocSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath);

//========================== Hyper-localization =============================
Eigen::MatrixXf hyperlocalisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, const Eigen::Ref<const Eigen::ArrayXf>& clippingCells, int nb_clipCells, int clipCellIndex, int cwipiParams, int count_obs, int obs_hyperloc, double sigmaLocX, double sigmaLocY, double sigmaLocZ, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath);

//========================== Inflation function =============================
Eigen::MatrixXf inflation(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, int cwipiMembers, int nb_cells, int cwipiParams, double stateInfl, double paramsInfl, float typeInfl, float paramEstSwitch, float cwipiVerbose, std::string stringRootPath);

//========================== Randomize observations =============================
Eigen::MatrixXf randomize_Obs(const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiMembers, int cwipiObs, int cwipiObsU, int cwipiParamsObs, double sigmaUserU, double sigmaUserp, double sigmaUserCf, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float cwipiVerbose);

//========================== Randomize observations =============================
Eigen::MatrixXf calculate_R(const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiObs, int cwipiObsU, int cwipiParamsObs, double sigmaUserU, double sigmaUserp, double sigmaUserCf, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float cwipiVerbose);

//**********================= Base EnKF =============================*********
Eigen::MatrixXf calculate_K(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const Eigen::Ref<const Eigen::MatrixXf>& obsMatrix, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::MatrixXf>& R, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiParams, float cwipiVerbose);

//***********=============== Hyper EnKF ============================**********
Eigen::MatrixXf EnKF_hyperloc(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsArray, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, double sigmaUserU, double sigmaUserp, double sigmaUserCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, float clippingSwitch, float hyperlocSwitch, int cwipiParams, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float paramEstSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath);

// *************************** MAIN EnKF Function *******************************
Eigen::MatrixXf mainEnKF(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const fvMesh& mesh, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, double sigmaUserU, double sigmaUserp, double sigmaUserCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, float clippingSwitch, float hyperlocSwitch, int cwipiParams, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float paramEstSwitch, float stateEstSwitch, float cwipiVerbose, std::string stringRootPath, int cwipiTimedObs, double obsTimeStep, double time, const double init_time, double epsilon);

//===========================================================================
void prepare_sendBack(double *sendValues, double *paramsSendValues, double *values, const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, int nb_cells, int cwipiParams, int index, float cwipiVerbose);

//===========================================================================
void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, std::string name);
