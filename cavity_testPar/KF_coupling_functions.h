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

void configuration(double* configValues);

Eigen::MatrixXf doClipping(Eigen::MatrixXf stateVector, Eigen::Ref<Eigen::MatrixXf> invStateVector, int& nb_cells, int nb_p, int nb_e, float cwipiVerbose, std::string stringRootPath);

Eigen::MatrixXf undoClipping(Eigen::MatrixXf UptMatrix, const Eigen::Ref<const Eigen::MatrixXf>& invStateVector, int inv_nb_cells, int& nb_cells, int nb_p, int nb_e, float cwipiVerbose, std::string stringRootPath); 

void define_mesh(double* pointsCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], float cwipiVerbose);

Eigen::ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, int nb_oU, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

Eigen::ArrayXf obs_Data_timed(int nb_e, int nb_cells, int nb_o, int nb_oU, int obsIndex, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

Eigen::MatrixXf samp_Data(int nb_e, int nb_o, int nb_oU, int velocityCase, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath);

Eigen::MatrixXf KF(const Eigen::Ref<const Eigen::MatrixXf>& velo_field_mat, const Eigen::Ref<const Eigen::ArrayXf>& obs_field, const Eigen::Ref<const Eigen::MatrixXf>& proj_field_mat, int nb_e, int nb_cells, int nb_o, int nb_oU, double sigma_userU, double sigma_userp, double sigma_userCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, 
float clippingSwitch, int nb_p, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float paramEstSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath);

void KF_output(double *sendValues, double *paramsSendValues, const Eigen::Ref<const Eigen::MatrixXf>& UptMatrix, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsMatrix, int nb_e, int nb_cells, double time, int nb_p, int nb_oU, int nb_o, int cwipiParamsObs, int velocityCase, int index, double epsilon, float cwipiVerbose);

void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXf>& stateVector, std::string name);

Eigen::MatrixXf localisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, int nb_cells, int nb_p, int nb_o, int nb_oU, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float clippingSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath);

