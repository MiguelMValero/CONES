#include "fvCFD.H"
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include <cwipi.h>

void define_mesh(double* pointsCoords, int* connecIdx, int* connec, const fvMesh& mesh);

Eigen::MatrixXf repro_Members(const int nb_e, int nb_cells, double time, double* values, double* paramsValues, int nb_p);

Eigen::ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, double time);

Eigen::MatrixXf samp_Data(int nb_e, int nb_o);

Eigen::MatrixXf KF(Eigen::MatrixXf velo_field_mat, Eigen::ArrayXf obs_field, Eigen::MatrixXf proj_field_mat, int nb_e, int nb_cells, int nb_o, double sigma_user, int nb_p);

void  KF_output(double *sendValues, double *paramsSendValues, Eigen::MatrixXf stateVector, int nb_e, int nb_cells, int nb_o, double time, int nb_p);