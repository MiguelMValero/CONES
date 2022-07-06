//#include "fvCFD.H"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <cmath>

#include "fvCFD.H"
#include <mpi.h>

#include <cwipi.h>

// using namespace std;
using namespace Eigen;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void define_mesh(double* pointCoords, int* connecIdx, int* connec, const fvMesh& mesh)
{
    // //char indexChar[50];
    // //sprintf(indexChar, "%i", index+1);

    // //char folder[50] = {"runCase_"};
    // char casePath[250] = {"/home/miguel/CWIPI_OpenFOAM"};
    // //char RunCase[250] = "";

    // char RunCase[250] = {"first_test"};

    // //strcat(RunCase,folder);
    // //strcat(RunCase,indexChar);

    // //-- Declaration des variables openFoam --
    // Foam::Time runTime(Foam::Time::controlDictName, casePath, RunCase);
    
    // Foam::fvMesh mesh
    // (
    //     Foam::IOobject
    //     (
    //         Foam::fvMesh::defaultRegion,
    //         runTime.timeName(),
    //         runTime,
    //         Foam::IOobject::MUST_READ
    //     )
    // );

    Info << "Definition du maillage KF" << nl << endl;

    scalar nCells = mesh.nCells();
    scalar nPoints = mesh.nPoints();

    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }
    
    // Info << "CC 2" << nl << endl;

    connecIdx[0]=0;
    forAll(mesh.cells(),i)
    {
        connecIdx[i+1]=connecIdx[i]+8;
    }

    // Info << "CC 3" << nl << endl;

    forAll(mesh.cells(),i)
    {
        forAll(mesh.cellShapes()[i],j)
        {
            connec[8*i+j]=mesh.cellShapes()[i][j]+1;
        }
    }

    // Info << "CC 4" << nl << endl;
    
    cwipi_dump_application_properties();
    cwipi_define_mesh("cwipiFoamCoupling", nPoints, nCells, pointCoords,connecIdx,connec);
    cwipi_locate("cwipiFoamCoupling");
}

//===========================================================================

MatrixXf repro_Members(const int nb_e, int nb_cells, double time, double* values, double* paramsValues, int nb_p)
{    
    MatrixXf velo_field_mat = MatrixXf::Zero(nb_cells*3 + nb_p, nb_e);

    // std::cout << "values[0] =" << values[0] << std::endl;
    std::cout << "Reproduction phase"<< std::endl;

    float gaussample;

    //==== Define random generator with Gaussian distribution ====
    const double mean = 0;
    const double stddev = 0.5;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> dist(mean, stddev);

    // std::cout << "seed =" << seed << std::endl;

    for (int j=0;j<nb_e;j++)
    {
        if (j>0)
        {
            do
            {
                gaussample=dist(generator);
            }while (gaussample<(mean-2*stddev) || gaussample>(mean+2*stddev));
        }

        // std::cout << "gaussample =" << gaussample << std::endl;

        for (int i = 0; i < nb_cells; i++)
        {
            velo_field_mat(i , j) = values[3*i+0] + values[3*i+0] * (gaussample);
            velo_field_mat(i + nb_cells , j) = values[3*i+1] + values[3*i+1] * (gaussample);
            velo_field_mat(i + 2*nb_cells , j) = values[3*i+2] + values[3*i+2] * (gaussample);
        }
        velo_field_mat(3*nb_cells, j) = paramsValues[0]+paramsValues[0]*gaussample;
    }

    // std::cout << "velo_field_mat[0] =" << velo_field_mat(0,0) << std::endl;

    //==== Writing txt file with matrix values ====
    std::string stepNumber = std::to_string(time);
    std::string filename = "UMat"+stepNumber;
    std::fstream file_out;
    
    file_out.open(filename, std::ios_base::out);
    
    if (!file_out.is_open()) 
    {
        std::cout << "failed to open " << filename << '\n';
    } 
    else 
    {
        for (int i=0; i<nb_cells*3+nb_p; i++)
        {
            for(int j=0;j<nb_e;j++)
            {
                file_out << velo_field_mat(i , j) << ' ';
	        }
            file_out << "\n";
        }
    }
    std::cout << "Done Writing Matrix on txt file!" << std::endl << "\n";

    return velo_field_mat;
}

//===========================================================================

ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, double time)
{        
    //==== Extract observation data from file ====
    std::cout << "Opening observation file" << std::endl;
    std::string obs_file="/home/miguel/OpenFOAM/miguel-8/run/CWIPI_OpenFOAM/first_test/observation_velocities.txt";

 
	std::vector<std::vector<string>> obs_content;
	std::vector<string> obs_row;
	std::string obs_line, obs_word;
 
	std::fstream obs_stream (obs_file, std::ios::in);
	if(obs_stream.is_open())
	{
		while(getline(obs_stream, obs_line))
		{
			obs_row.clear();
 
			std::stringstream obs_str(obs_line);
 
			while(getline(obs_str, obs_word, ','))
				obs_row.push_back(obs_word);
			obs_content.push_back(obs_row);
		}
	}
	else
		std::cout << "Could not open the file" << std::endl << "\n";


    ArrayXf obs_field(nb_o*3);
    // ArrayXf coord_obs_field(nb_o*3);

    std::cout << "Obs file is open, begin to fill arrays" << std::endl;
 
	for(int i=0; i<nb_o; i++)
	{
        obs_field(i)=stod(obs_content[i][0]);
        obs_field(i+nb_o)=stod(obs_content[i][1]);
        obs_field(i+2*nb_o)=stod(obs_content[i][2]);
        
        // coord_obs_field(i)=stod(obs_content[i][3]);
        // coord_obs_field(i+nb_o)=stod(obs_content[i][4]);
        // coord_obs_field(i+2*nb_o)=stod(obs_content[i][5]);       
	}

    std::cout << "obs_field created" << std::endl << "\n";
    // for (int i=0; i<nb_o; i++)
    // {
    //     cout << obs_field(i) << " " << obs_field(i+nb_o)<< " " << obs_field(i+2*nb_o) << "   for loop " << i << endl ; 
    // }
    // cout << endl << "\n";

    return obs_field;
}

//===========================================================================

Eigen::MatrixXf samp_Data(int nb_e, int nb_o)
{
    ArrayXf sampArray = ArrayXf::Zero(nb_o*3);
    // sampArray << 1.1364, 1.20305, 1.24147, 1.25221, 1.23345, 0.0322568, 0.0189137, 0.00712773, -0.00441124, -0.0168068, 0, 0, 0, 0, 0;
    std::string UInt2 = "/home/miguel/OpenFOAM/miguel-8/run/CWIPI_OpenFOAM/first_test/UInt.txt";
    std::ifstream file;
    file.open(UInt2);
    std::string data = "";
    int count = 0;

    if (file.is_open())
    {
        std::string line;
        while( std::getline(file,line) )
        {
            std::stringstream ss(line);

            std::string compX, compY, compZ;
            std::getline(ss, compX,',');
			sampArray(3*count) = stod(compX);
            std::getline(ss, compY,',');
			sampArray(1 + 3*count) = stod(compY);
            std::getline(ss, compZ,','); 
			sampArray(2 + 3*count) = stod(compZ);

        }
        ++count;
    }

    file.close();
    
    MatrixXf sampMatrix = MatrixXf::Zero(nb_o*3, nb_e);

    float gaussample;

    //==== Define random generator with Gaussian distribution ====
    const double mean = 0;
    const double stddev = 0.05;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> dist(mean, stddev);

    for (int j=0;j<nb_e;j++)
    {
        if (j>0)
        {
            do
            {
                gaussample=dist(generator);
            }while (gaussample<(mean-3*stddev) || gaussample>(mean+3*stddev));
        }

        for (int i = 0; i < nb_o; i++)
        {
            sampMatrix(i , j) = sampArray(i) + sampArray(i) * (gaussample);
            sampMatrix(i + nb_o , j) = sampArray(i+nb_o) + sampArray(i+nb_o) * (gaussample);
            sampMatrix(i + 2*nb_o , j) = sampArray(i+2*nb_o) + sampArray(i+2*nb_o) * (gaussample);
        }
    }
    
    return sampMatrix;
}

//===========================================================================

MatrixXf KF(MatrixXf velo_field_mat, ArrayXf obs_field, MatrixXf proj_field_mat, int nb_e, int nb_cells, int nb_o, double sigma_user, int nb_p)
{
    std::cout << "Opening KF" << std::endl;
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // --- Initialisation des variables ---
    // int nb_cells=2;
    int nb_cells_cmpnt=nb_cells*3;
    int nb_o_cmpnt=nb_o*3; //inférieur à nb_cells_cmpnt/3
    float gaussample;

    //========== Initialisation des matrices ==========
    MatrixXf velo_field_mat_upt(nb_cells_cmpnt + nb_p,nb_e);
    ArrayXf velo_field_mean(nb_cells_cmpnt + nb_p);
    MatrixXf velo_field_mat_norm(nb_cells_cmpnt + nb_p,nb_e);

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

    std::cout << "Initialisation done" << std::endl;

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

    std::cout << "Perturbations added" << std::endl;

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
    std::cout << "nb_cells_cmpnt = " << nb_cells_cmpnt << std::endl;


    for (int i=0;i<nb_cells_cmpnt+nb_p;i++)
    {
        velo_field_mean(i)=velo_field_mat.row(i).mean();
    }

    std::cout << "Velo means calculated" << std::endl;

    std::cout << "nb_o_cmpnt = " << nb_o_cmpnt << std::endl;

    for (int i=0;i<nb_o_cmpnt;i++)
    {
        err_mat_mean(i)=err_mat.row(i).mean();
        obs_field_mean(i)=obs_field_mat.row(i).mean();
        proj_field_mean(i)=proj_field_mat.row(i).mean();
    }

    std::cout << "Means calculated" << std::endl;
    
    // std::cout << "velo_field_mean :" << std::endl;
    // std::cout << velo_field_mean << std::endl << "\n";

    // std::cout << "err_mat_mean :" << std::endl;
    // std::cout << err_mat_mean << std::endl << "\n";

    // std::cout << "proj_field_mean :" << std::endl;
    // std::cout << proj_field_mean << std::endl << "\n";



    // --- Normalized anomalies (Ue_mat-Ume)/sqrt(Ne-1) ---
    for (int i=0;i<nb_cells_cmpnt+nb_p;i++)
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

    std::cout << "Anomaly matrices calculated" << std::endl;

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

    std::cout << "Kalman gain calculated" << std::endl;
    
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
    
    std::cout << "Filtre de Kalman termine" << std::endl << "\n";
    
    return velo_field_mat_upt;
}

//===========================================================================

void KF_output(double *sendValues, double *paramsSendValues, MatrixXf stateVector, int nb_e, int nb_cells, int nb_o, double time, int nb_p)
{
    std::cout << "Opening output" << std::endl;
    
    int nb_cells_cmpnt=nb_cells*3;

    ArrayXf stateVector_mean(nb_cells_cmpnt + nb_p);

    //==== Calculating mean ====
    for (int i=0 ; i<nb_cells_cmpnt + nb_p ; i++)
    {
        stateVector_mean(i)=stateVector.row(i).mean();
    }

    //==== Writing txt file with matrix values ====
    std::string stepNumber = std::to_string(time);
    std::string filename = "UMat_upt"+stepNumber;
    std::fstream file_out;
    
    file_out.open(filename, std::ios_base::out);
    
    if (!file_out.is_open()) 
    {
        std::cout << "failed to open " << filename << '\n';
    } 
    else 
    {
        for (int i=0; i<nb_cells_cmpnt + nb_p; i++)
        {
            for(int j=0;j<nb_e;j++)
            {
                file_out << stateVector(i , j) << ' ';
	        }
            file_out << "\n";
        }
    }
    std::cout << "Done Writing Updated Matrix on txt file!" << std::endl << "\n";
    

    //==== Writing txt file with mean array values ====
    std::string filename_mean = "UMean_upt"+stepNumber;
    std::fstream file_out_mean;
    
    file_out_mean.open(filename_mean, std::ios_base::out);
    
    if (!file_out_mean.is_open()) 
    {
        std::cout << "failed to open " << filename_mean << '\n';
    } 
    else 
    {
        for(int i=0 ; i<nb_cells ; i++)
        {
            file_out_mean << stateVector_mean(i) << ' ' << stateVector_mean(i+nb_cells) << ' ' << stateVector_mean(i+2*nb_cells) << '\n';
	    }
        file_out_mean << stateVector_mean[nb_cells_cmpnt];
        std::cout << "Done Writing output mean on txt file!" << std::endl << "\n";
    }

    for (int i=0; i<nb_cells; i++)
    {   
        sendValues[3*i + 0] = stateVector_mean(i);
        sendValues[3*i + 1] = stateVector_mean(i + nb_cells);    
        sendValues[3*i + 2] = stateVector_mean(i + 2*nb_cells);    
    }
    paramsSendValues[0] = stateVector_mean(nb_cells_cmpnt);
}

// ************************************************************************* //