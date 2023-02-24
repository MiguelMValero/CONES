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

void heapify(int arr[], int arr_index[], int n, int i);
void heapSort(int arr[], int arr_index[], int n);
void printArray(int arr[], int arr_index[], int n);
void define_mesh(double* pointsCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], int cwipiVerbose);

Eigen::MatrixXf repro_Members(const int nb_e, int nb_cells, double time, double* values, double* paramsValues, int nb_p, int cwipiVerbose);

Eigen::ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, int nb_oU, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath);

Eigen::ArrayXf obs_Data_timed(int nb_e, int nb_cells, int nb_o, int nb_oU, int i, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath);

Eigen::MatrixXf samp_Data(int nb_e, int nb_o, int nb_oU, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath);

Eigen::MatrixXf KF(Eigen::MatrixXf velo_field_mat, Eigen::ArrayXf obs_field, Eigen::MatrixXf proj_field_mat, int nb_e, int nb_cells, int nb_o, int nb_oU, double sigma_userU, double sigma_userp, int nb_p, int cwipiParamsObs, int cwipiVerbose);

void  KF_output(double *sendValues, double *paramsSendValues, Eigen::MatrixXf stateVector, int nb_e, int nb_cells, double time, int nb_p, int index, int cwipiVerbose);

void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, Eigen::MatrixXf stateVector, std::string name);