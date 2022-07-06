//#include "fvCFD.H"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <Eigen/Dense>
//#include "interfaceC_daicofoam.h"

using namespace Eigen;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C" int kf_unit_()
{
    std::cout << "Opening kf_unit" << std::endl << "\n";
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // --- Initialisation des variables ---
    // int nb_cells=2;
    int nb_cells_cmpnt=nb_cells*3;
    int nb_o_cmpnt=nb_o*3; //inférieur à nb_cells_cmpnt/3
    float gaussample;

    //Create velocity matrix to receive fields
    MatrixXf velo_field_mat(nb_cells_cmpnt+1,nb_e);
    double velo_field_mat_inp[nb_cells_cmpnt+1][nb_e];

    ArrayXf obs_field(nb_o_cmpnt);
    double obs_field_inp[nb_o_cmpnt];

    MatrixXf proj_field_mat(nb_o_cmpnt,nb_e);
    double proj_field_mat_inp[nb_o*3][nb_e];


    //==== Receiving velocity field ====
    //Add variables to receive velo_field
    int il_tag_matinp = PL_NO_TAG;
    int il_time_matinp = PL_NO_TIME;
    int il_err_matinp = 0 ;
    char cla_obj_matinp[PL_LNAME], cla_space_matinp[PL_LNAME];
    
    // //Receive "velo_field"
    sprintf (cla_space_matinp,"mat2d");
    sprintf (cla_obj_matinp,"velo_field_mat_inp");
    il_err_matinp=PALM_Get(cla_space_matinp,cla_obj_matinp,&il_time_matinp,&il_tag_matinp,&velo_field_mat_inp);

    //==== Receiving observation field ====
    //Add variables to receive obs_field
    int il_tag_arr = PL_NO_TAG;
    int il_time_arr = PL_NO_TIME;
    int il_err_arr = 0 ;
    char cla_obj_arr[PL_LNAME], cla_space_arr[PL_LNAME];
    
    //Create velocity matrix to receive OF field 
    sprintf (cla_space_arr,"vect1d");
    sprintf (cla_obj_arr,"obs_field_inp");
    il_err_arr=PALM_Get(cla_space_arr,cla_obj_arr,&il_time_arr,&il_tag_arr,&obs_field_inp);


    //==== Receive projected fields matrix ====
    int il_tag_pfm = PL_NO_TAG;
    int il_time_pfm = PL_NO_TIME;
    int il_err_pfm = 0 ;
    char cla_obj_pfm[PL_LNAME], cla_space_pfm[PL_LNAME];
    
    sprintf (cla_space_pfm,"mat2d_proj");
    sprintf (cla_obj_pfm,"proj_field_mat");
    il_err_pfm=PALM_Get(cla_space_pfm,cla_obj_pfm,&il_time_pfm,&il_tag_pfm,&proj_field_mat_inp);

    //==== Copying fields to eigen matrix ====
    for (int i=0; i<nb_cells_cmpnt+1; i++)
    {   
        for (int j=0; j<nb_e; j++)
        {
            velo_field_mat(i,j)=velo_field_mat_inp[i][j];
        } 
    }

    for (int i=0; i<nb_o_cmpnt; i++)
    {
        obs_field(i)=obs_field_inp[i];
        for (int j=0; j<nb_e; j++)
        {
            proj_field_mat(i,j)=proj_field_mat_inp[i][j];
        }
    } 


    MatrixXf velo_field_mat_upt(nb_cells_cmpnt+1,nb_e);
    ArrayXf velo_field_mean(nb_cells_cmpnt+1);
    MatrixXf velo_field_mat_norm(nb_cells_cmpnt+1,nb_e);

    MatrixXf obs_field_mat(nb_o_cmpnt,nb_e);
    ArrayXf obs_field_mean(nb_o_cmpnt);
    MatrixXf obs_field_mat_norm(nb_o_cmpnt,nb_e);

    ArrayXf proj_field_mean(nb_o_cmpnt);
    MatrixXf proj_field_mat_norm(nb_o_cmpnt,nb_e);

    MatrixXf err_mat(nb_o_cmpnt,nb_e);
    ArrayXf err_mat_mean(nb_o_cmpnt);
    MatrixXf err_mat_norm(nb_o_cmpnt,nb_e);

    MatrixXf K;
    MatrixXf K1;
    MatrixXf K2;
    MatrixXf K21;

    // Define random generator with Gaussian distribution
    const double mean_obs = 0;
    const double stddev_obs = sigma_user;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> dist_obs(mean_obs, stddev_obs);
    
    MatrixXf R(nb_o_cmpnt,nb_o_cmpnt);
    R.setIdentity();
    R=R*sigma_user;

    // Add Gaussian noise to create random velocity field and observation matrix
    for (int j=0;j<nb_e;j++)
    {
        for (int i=0;i<nb_o_cmpnt;i++)
        {
            do
            {
                gaussample=dist_obs(generator);
            }while (gaussample<(mean_obs-3*stddev_obs) || gaussample>(mean_obs+3*stddev_obs));
            err_mat(i,j)=gaussample;
        }

        obs_field_mat.col(j)=obs_field;
    }
    
    obs_field_mat=obs_field_mat+err_mat;

    //proj_field_mat=H_mod*velo_field_mat;

    // std::cout << "velo_field_mat :" << std::endl;
    // std::cout << velo_field_mat << std::endl << "\n";

    // std::cout << "err_mat :" << std::endl;
    // std::cout << err_mat << std::endl << "\n";

    // std::cout << "obs_field_mat :" << std::endl;
    // std::cout << obs_field_mat << std::endl << "\n";

    // std::cout << "proj_field_mat :" << std::endl;
    // std::cout << proj_field_mat << std::endl << "\n";

    // --- Ensemble means ---/
    for (int i=0;i<nb_cells_cmpnt+1;i++)
    {
        velo_field_mean(i)=velo_field_mat.row(i).mean();
    }

    for (int i=0;i<nb_o_cmpnt;i++)
    {
        err_mat_mean(i)=err_mat.row(i).mean();
        obs_field_mean(i)=obs_field_mat.row(i).mean();
        proj_field_mean(i)=proj_field_mat.row(i).mean();
    }
    
    // std::cout << "velo_field_mean :" << std::endl;
    // std::cout << velo_field_mean << std::endl << "\n";

    // std::cout << "err_mat_mean :" << std::endl;
    // std::cout << err_mat_mean << std::endl << "\n";

    // std::cout << "proj_field_mean :" << std::endl;
    // std::cout << proj_field_mean << std::endl << "\n";



    // --- Normalized anomalies (Ue_mat-Ume)/sqrt(Ne-1) ---
    for (int i=0;i<nb_cells_cmpnt+1;i++)
    {
        for (int j=0;j<nb_e;j++)
        {
            velo_field_mat_norm(i,j)=(velo_field_mat(i,j)-velo_field_mean(i))/std::sqrt(nb_e-1);
        }
    }

    for (int i=0;i<nb_o_cmpnt;i++)
    {
        for (int j=0;j<nb_e;j++)
        {
            proj_field_mat_norm(i,j)=(proj_field_mat(i,j)-proj_field_mean(i))/std::sqrt(nb_e-1);
            err_mat_norm(i,j)=(err_mat(i,j)-err_mat_mean(i))/std::sqrt(nb_e-1);
            obs_field_mat_norm(i,j)=(obs_field_mat(i,j)-obs_field_mean(i))/std::sqrt(nb_e-1);
        }
    }

    // std::cout << "velo_field_mat_norm :" << std::endl;
    // std::cout << velo_field_mat_norm << std::endl << "\n";

    // std::cout << "proj_field_mat_norm :" << std::endl;
    // std::cout << proj_field_mat_norm << std::endl << "\n";

    // --- Kalman gain NAe_mat*trans(NAo_mat)*inv(NAo_mat*trans(NAo_mat))--- 
   
    // K21=proj_field_mat_norm*(proj_field_mat_norm.transpose())+obs_field_mat_norm*(obs_field_mat_norm.transpose());
    K21=proj_field_mat_norm*(proj_field_mat_norm.transpose())+R;
    K2=K21.inverse();
    //K22=K21.determinant();

    K1=velo_field_mat_norm*(proj_field_mat_norm.transpose());

    K=K1*K2;
    
    // std::cout << "K2 :" << std::endl;
    // std::cout << K2 << std::endl << "\n";

    // std::cout << "K1 :" << std::endl;
    // std::cout << K1 << std::endl << "\n";

    // std::cout << "K :" << std::endl;
    // std::cout << K << std::endl << "\n";
    
    // std::cout << "ve1 :" << std::endl;
    // std::cout << velo_field_mat_norm.col(j) << std::endl << "\n";

    // std::cout << "proj1 :" << std::endl;
    // std::cout << (proj_field_mat_norm.transpose()).row(j) << std::endl << "\n";

    // --- Update of velocities Ue_mat + K_mat*(Uo_mat-H_mat*(Ue_mat))---
    std::cout << "Debut de la phase d'update" << std::endl;
    for (int j=0;j<nb_e;j++)
    {
        // std::cout << "======== Boucle numero " << j+1 << " ========" << std::endl << "\n";
        velo_field_mat_upt.col(j)=velo_field_mat.col(j)+K*(obs_field_mat.col(j)-proj_field_mat.col(j));

        // std::cout << "velo_field_mat :" << std::endl;
        // std::cout << velo_field_mat << std::endl << "\n";

        // std::cout << "velo_field_mat_upt :" << std::endl;
        // std::cout << velo_field_mat_upt << std::endl << "\n";
    }

    
    //==== OpenPalm Output ====
    double velo_field_mat_upt_out[nb_cells_cmpnt+1][nb_e];

    for (int i=0; i<nb_cells_cmpnt+1; i++)
    {   
        for (int j=0; j<nb_e; j++)
        {
            velo_field_mat_upt_out[i][j]=velo_field_mat_upt(i,j);
        } 
    }
    
    //Add variables to receive velo_field
    int il_tag_matout = PL_NO_TAG;
    int il_time_matout = PL_NO_TIME;
    int il_err_matout = 0 ;
    char cla_obj_matout[PL_LNAME], cla_space_matout[PL_LNAME];
    
    // //Receive "velo_field"
    sprintf (cla_space_matout,"mat2d");
    sprintf (cla_obj_matout,"velo_field_mat_upt_out");
    il_err_matout=PALM_Put(cla_space_matout,cla_obj_matout,&il_time_matout,&il_tag_matout,&velo_field_mat_upt_out);

    std::cout << "Filtre de Kalman termine" << std::endl << "\n";
    
    return 0;
}


// ************************************************************************* //
