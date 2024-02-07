//#include "fvCFD.H"
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

#include "fvCFD.H"
#include <mpi.h>

#include <cwipi.h>

// using namespace std;
using namespace Eigen;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void configuration(double* configValues){

    std::ifstream cFile ("cwipiConfig");
    int k = 0;

    if (cFile.is_open())
    {
        std::string line;
        while(std::getline(cFile, line))
        {
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());
            if( line.empty() || line[0] == '#' )
            {
                continue;
            }
            auto delimiterPos = line.find("=");
            auto value = line.substr(delimiterPos + 1);

            configValues[k] = stod(value);
            k = k+1;
        }
    }
    else 
    {
        std::cerr << "Couldn't open configuration file for reading.\n";
    }
}

//===========================================================================
void tic(int mode=0) {
	static std::chrono::_V2::system_clock::time_point t_start;

	if (mode==0)
		t_start = std::chrono::high_resolution_clock::now();
	else {
		auto t_end = std::chrono::high_resolution_clock::now();
		std::cout << "Elapsed time for the EnKF selected procedure is " << (t_end - t_start).count()*1E-9 << " seconds\n";
	}
}

void toc() { tic(1); }

//===========================================================================
MatrixXf doClipping(const Eigen::Ref<const Eigen::MatrixXf>& invStateMatrix, int& nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose, std::string stringRootPath)
{
    if (cwipiVerbose) std::cout << "Entering the doClipping function" << std::endl;
    std::string clipping_file = stringRootPath + "/clippingCells.txt";

    std::ifstream clippingCells;
    clippingCells.open(clipping_file);

    ArrayXf clippingCellID;
    if (clippingCells.is_open())
    {
        std::string line;
        int count = 0;

        while( std::getline(clippingCells, line) )
        {
            std::stringstream ss(line);
            std::string clippingCellID_str;
            std::getline(ss, clippingCellID_str, ',');

            if (count == 0){
                clippingCellID.resize(std::stoi(clippingCellID_str));
            }
            else{
                clippingCellID(count-1) = std::stoi(clippingCellID_str);
            }
            ++count;
		}
	}
    else
    {
        Info<< "Couldn't open clippingCells.txt" << endl;
    }
    clippingCells.close();

    int nb_cellsClipping = clippingCellID.size();
    MatrixXf stateMatrixClipped = MatrixXf::Zero(3*nb_cellsClipping + cwipiParams, cwipiMembers);

    int clipIndex = 0;

    for (int i = 0; i < cwipiMembers; ++i){
        for (int j = 0; j < nb_cellsClipping; ++j){
            clipIndex = clippingCellID(j);
            stateMatrixClipped(j, i) = invStateMatrix(clipIndex, i);
            stateMatrixClipped(nb_cellsClipping + j, i) = invStateMatrix(clipIndex + nb_cells, i);
            stateMatrixClipped(2*nb_cellsClipping + j, i) = invStateMatrix(clipIndex + 2*nb_cells, i);
        }

        for (int j = 0; j < cwipiParams; ++j){
            stateMatrixClipped(3*nb_cellsClipping+j, i) = invStateMatrix(3*nb_cells+j, i);
        }
    }

    nb_cells = nb_cellsClipping;
    if (cwipiVerbose) std::cout << "Leaving the doClipping function" << std::endl;

    return stateMatrixClipped;
}

//===========================================================================
MatrixXf undoClipping(Eigen::MatrixXf stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXf>& invStateMatrix, int inv_nb_cells, int& nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose, std::string stringRootPath)
{
    std::string clipping_file = stringRootPath + "/clippingCells.txt";
    std::ifstream clippingCells;
    clippingCells.open(clipping_file);
    if (cwipiVerbose) std::cout << "Entering the undoClipping function" << std::endl;

    ArrayXf clippingCellID;

    if (clippingCells.is_open())
    {
        std::string line;
        int count = 0;
        while( std::getline(clippingCells, line) )
            {
            std::stringstream ss(line);
            std::string clippingCellID_str;
            std::getline(ss, clippingCellID_str, ',');

            if (count == 0){
                clippingCellID.resize(std::stoi(clippingCellID_str));
            }
            else{
                clippingCellID(count-1) = std::stoi(clippingCellID_str);
            }
            ++count;
		}
	}
    else
    {
        Info<< "Couldn't open clippingCells.txt" << endl;
    }
    clippingCells.close();

    MatrixXf tempStateMatrixUpt = stateMatrixUpt;
    int nb_cellsClipping = clippingCellID.size();
    stateMatrixUpt.resize(inv_nb_cells + cwipiParams, cwipiMembers);
    stateMatrixUpt = invStateMatrix;

    int clipIndex = 0;

    for (int i = 0; i < cwipiMembers; ++i){
        for (int j = 0; j < nb_cellsClipping; ++j){
            clipIndex = clippingCellID(j);
            stateMatrixUpt(clipIndex, i) = tempStateMatrixUpt(j, i);
            stateMatrixUpt(inv_nb_cells + clipIndex, i) = tempStateMatrixUpt(j + nb_cellsClipping, i);
            stateMatrixUpt(2*inv_nb_cells + clipIndex, i) = tempStateMatrixUpt(j + 2*nb_cellsClipping, i);
        }

        for (int j = 0; j < cwipiParams; ++j){
            stateMatrixUpt(3*inv_nb_cells+j, i) = tempStateMatrixUpt(3*nb_cellsClipping+j, i);
        }
    }

    nb_cells = inv_nb_cells;

    if (cwipiVerbose) std::cout << "Leaving the undoClipping function" << std::endl;
    return stateMatrixUpt;
}

//===========================================================================
void define_mesh(double* pointCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], float cwipiVerbose)
{
    //========== Coupling mesh definition ==========
    
    if (cwipiVerbose) Info << "Implementation of mesh in EnKF side" << nl << endl;

    int nCells = mesh.nCells();
    int nPoints = mesh.nPoints();
    int nFaces = mesh.nFaces();

    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }

    // connecIdx[0]=0;
    // forAll(mesh.cells(),i)
    // {
    //     connecIdx[i+1]=connecIdx[i]+8;
    // }

    // forAll(mesh.cells(),i)
    // {
    //     forAll(mesh.cellShapes()[i],j)
    //     {
    //         connec[8*i+j]=mesh.cellShapes()[i][j]+1;
    //     }
    // }

    face_index[0]=0;
    face_connectivity_index[0]=0;

    int *face_index_temp = new int[nCells];
    int *face_index_temp_IDs = new int[nCells];
    int *face_connectivity_index_temp = new int[nFaces];
    int *face_connectivity_index_temp_IDs = new int[nFaces];

	forAll(mesh.cells(), i)
    {
        face_index_temp[i] = mesh.cells()[i].size();
        face_index_temp_IDs[i] = i;
    }

	forAll(mesh.cells(), i)
    {
        face_index[i+1] = face_index[i] + face_index_temp[i];
    }

    int faceOwner = 0;
    int faceNeighbour = 0;
    int faceID = 0;
	int faces_count = 0;
	for(int i = 0; i < nCells; i++){
        for (int j = 0; j < face_index_temp[i] ; j++) {
            faceID = mesh.cells()[face_index_temp_IDs[i]][j];
            if (mesh.isInternalFace(faceID)){
                faceOwner = mesh.faceOwner()[faceID];
                faceNeighbour = mesh.faceNeighbour()[faceID];
                if (faceOwner == face_index_temp_IDs[i]){
                    if (faceOwner > faceNeighbour)
                    {
                        cell_to_face_connectivity[faces_count+j] = -(faceID + 1);
                    }
                    else
                    {
                        cell_to_face_connectivity[faces_count+j] = faceID + 1;
                    }
                }
                else if (faceNeighbour == face_index_temp_IDs[i]){
                    if (faceOwner > faceNeighbour)
                    {
                        cell_to_face_connectivity[faces_count+j] = faceID + 1;
                    }
                    else
                    {
                        cell_to_face_connectivity[faces_count+j] = -(faceID + 1);
                    }
                }
                else
                    break;
            }
            else
            {
                cell_to_face_connectivity[faces_count+j] = faceID + 1;
            }
        }
        faces_count = faces_count + face_index_temp[i];
    }

	forAll(mesh.faces(), i)
    {
        face_connectivity_index_temp[i] = mesh.faces()[i].size();
        face_connectivity_index_temp_IDs[i] = i;
    }

    forAll(mesh.faces(), i)
    {    
        face_connectivity_index[i+1] = face_connectivity_index[i] + face_connectivity_index_temp[i];
    }

	int vertices_cout = 0;
    for(int i = 0; i < nFaces; i++)
    {
        for (int j = 0; j < face_connectivity_index_temp[i] ; j++)
        {
            face_connectivity[vertices_cout+j] = mesh.faces()[face_connectivity_index_temp_IDs[i]][j] + 1;
        }
		vertices_cout = vertices_cout + face_connectivity_index_temp[i];
    }

    int t[2];
    t[0] = 0;
    t[1] = 0;

    cwipi_define_mesh(cl_coupling_name, nPoints, 0, pointCoords,t,t);

    cwipi_add_polyhedra(cl_coupling_name, 
                        nCells,
                        face_index,
                        cell_to_face_connectivity,
                        nFaces,
                        face_connectivity_index,
                        face_connectivity);
 
    if (cwipiVerbose) Info << "Mesh defined in EnKF side" << nl << endl;

    cwipi_locate(cl_coupling_name);    

    delete[] face_index_temp;
    delete[] face_index_temp_IDs;
    delete[] face_connectivity_index_temp;
    delete[] face_connectivity_index_temp_IDs;
}

//===========================================================================
ArrayXf obs_Data(int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{        
    //========== Not depending on time : Read the observation data from a file named "obs_field.txt" only and return the values
    // in a eigen array. The txt file has to contain 
    // 1) velocity values (first all values for X, secondly all values for Y and finally all values for Z) 
    // 2) pressure values 
    // 3) velocity and pressure values ==========
    
    //=========== Extract observation data from file ===========
    std::string obs_file = stringRootPath + "/obs_field.txt";
 
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
	else {
		std::cerr << "Couldn't open observation file.\n";
    }

    ArrayXf obsArray(cwipiObs);
	for(int i=0; i<cwipiObs; i++)
	{
        obsArray(i) = std::stod(obs_content[i][0]);         
	}

    if (cwipiVerbose) std::cout << "observation field created" << std::endl << "\n";
    // if (cwipiVerbose == 1){
    // for (int i=0; i<cwipiObs; i++)
    // {
    //     std::cout << obsArray(i) << "   for loop " << i << std::endl << "\n" ;
    // }
    // std::cout << std::endl << "\n";
    // }

    return obsArray;
}

//===========================================================================
ArrayXf obs_Data_timed(int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, int obsIndex, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{        
    //========== Depending on time : Read the observations data from a file named "obs_field.txt" only and return the values 
    // of the specific time step in a eigen array. The txt file has to contain velocity values with X Y and Z in the same column
    // (first all the X, then Y then Z). Number of column = number of DA phases ==========
    
    //=========== Extract observation data from file ===========
    std::string obs_file= stringRootPath + "/obs_field.txt";
 
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
	else {
		std::cerr << "Couldn't open observation file.\n";
    }

    ArrayXf obsArray(cwipiObs); 
	for(int i=0; i<cwipiObs; i++)
	{
        obsArray(i) = std::stod(obs_content[i][obsIndex]);
    }
    
    if (cwipiVerbose) std::cout << "obsArray created" << std::endl << "\n";
    //if (cwipiVerbose == 1){
    // for (int i=0; i<cwipiObs; i++)
    // {
    //     std::cout << obsArray(i) << "   for loop " << i << std::endl << "\n" ;
    // }
    // std::cout << std::endl << "\n";
    // }

    return obsArray;
}

//===========================================================================
MatrixXf samp_Data(int cwipiMembers, int cwipiObs, int cwipiObsU, int velocityCase, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{
    //========== Read the sampled velocities from the files created by each OF instance. Return of matrix 
    // with all the sampled velocities organised in a X Y Z manner, each colunm corresponding to 
    // a member ==========     
    MatrixXf sampMatrix = MatrixXf::Zero(cwipiObs, cwipiMembers); // By default our parameter is the velocity
    std::string IntGen = stringRootPath + "/member";
    std::ifstream file;
    char* IntGenChar = new char[1000];
    strcpy(IntGenChar, IntGen.c_str());

    for (int i = 0; i < cwipiMembers; ++i)
    {
        char ensemble[50];
        char member[50];
        if (cwipiParamsObs == 0){
            std::string UInt;
            sprintf(IntGenChar, "%i", i+1);
            sprintf(ensemble, "/UInt%i", i+1);
            sprintf(member, "%i", i+1);
            UInt += IntGen; // Append the path
            UInt += member; // Append the member index
            UInt += ensemble; // Append file name
            if (cwipiVerbose) std::cout << "Opening " << UInt << std::endl;
            file.open(UInt);
        }
        else if (cwipiParamsObs == 1){
            std::string pInt;
            sprintf(ensemble, "/pInt%i", i+1);
            sprintf(member, "%i", i+1);
            pInt += IntGen; // Append the path
            pInt += member; // Append the member index
            pInt += ensemble; // Append file name
            if (cwipiVerbose) std::cout << "Opening " << pInt << std::endl;
            file.open(pInt);
        }
        else if (cwipiParamsObs == 2){
            std::string UpInt;
            sprintf(ensemble, "/UpInt%i", i+1);
            sprintf(member, "%i", i+1);
            UpInt += IntGen; // Append the path
            UpInt += member; // Append the member index
            UpInt += ensemble; // Append file name
            if (cwipiVerbose) std::cout << "Opening " << UpInt << std::endl;
            file.open(UpInt);
        }
        else if (cwipiParamsObs == 3){
            std::string UCfInt;
            sprintf(ensemble, "/UCfInt%i", i);
            sprintf(member, "%i", i);
            UCfInt += IntGen; // Append the path
            UCfInt += member; // Append the member index
            UCfInt += ensemble; // Append file name
	    if (cwipiVerbose) std::cout << "Opening " << UCfInt << std::endl;
            file.open(UCfInt);
        }

        int count = 0;
        if (file.is_open())
        {
            std::string line;
            while( std::getline(file,line) )
            {
                std::stringstream ss(line);
                std::string IntValue;
                std::getline(ss, IntValue, ' ');
                
                if (cwipiParamsObs == 0){
                    switch (velocityCase){
                        case 1 :
                            if (count < cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            break;
                        case 2 :
                            if (count >= cwipiObsU && count < 2*cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 3 :
                            if (count >= 2*cwipiObsU) sampMatrix(count-2*cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 4 :
                            if (count < 2*cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            break;
                        case 5 :
                            if (count < cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 2*cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 6 :
                            if (count >= cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 7 :
                            sampMatrix(count, i) = std::stod(IntValue);
                            break;
                    }
                }
                else if (cwipiParamsObs == 1){
                    sampMatrix(count, i) = std::stod(IntValue);
                }
                else if (cwipiParamsObs == 2 || cwipiParamsObs == 3){
                    switch (velocityCase){
                        case 1 :
                            if (count < cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 3*cwipiObsU) sampMatrix(count-2*cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 2 :
                            if (count >= cwipiObsU && count < 2*cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            else if (count >= 3*cwipiObsU) sampMatrix(count-2*cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 3 :
                            if (count >= 2*cwipiObsU) sampMatrix(count-2*cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 4 :
                            if (count < 2*cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 3*cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 5 :
                            if (count < cwipiObsU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 2*cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 6 :
                            if (count >= cwipiObsU) sampMatrix(count-cwipiObsU, i) = std::stod(IntValue);
                            break;
                        case 7 :
                            sampMatrix(count, i) = std::stod(IntValue);
                            break;
                    }
                }
                ++count;
            }

        }
	    else {
		    std::cerr << "Couldn't open any of the sampling files.\n";
        }
        file.close();
    }
    
    if (cwipiVerbose) std::cout << "Sampling matrix created" << std::endl << "\n";
    // if (cwipiVerbose == 1){
    // for (int i=0; i< cwipiObs; i++)
    // { 
    //     for (int j=0; j<cwipiMembers; j++)
    //     {
    //         std::cout << sampMatrix(i, j) << " ";
    //     }
    //     std::cout << std::endl << "\n";
    // }
    // std::cout << std::endl << "\n";
    // }

    return sampMatrix;
}

//===========================================================================
void EnKF_outputs(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiMembers, int nb_cells, double time, int cwipiParams, int cwipiObsU, int cwipiObs, int cwipiParamsObs, int velocityCase, float cwipiVerbose, double epsilon)
{
    //========== We do an output of all optimized coefficients to evaluate convergence ==========//
    std::fstream file_coeffs_out;
    std::string filename_coeffs = "results/UpdatedCoefficients";
    file_coeffs_out.open(filename_coeffs, std::ios_base::out | std::ios_base::app);
    file_coeffs_out << time << ' ';
    for (int i = 0; i < cwipiParams; ++i){
        double average = stateMatrixUpt.row(nb_cells*3 + i).mean();
        double std_dev_temp = 0;
        for (int j = 0; j < cwipiMembers; ++j){
            std_dev_temp = std_dev_temp + std::pow((stateMatrixUpt(nb_cells*3 + i, j) - average), 2);
        }
        double std_dev = std::sqrt(std_dev_temp/cwipiMembers);
        file_coeffs_out << average << ' ' << std_dev << ' ';
    }
    file_coeffs_out << "\n";
    file_coeffs_out.close();

    std::fstream file_RMS_out;
    std::string filename_RMS = "results/NRMSD";
    file_RMS_out.open(filename_RMS, std::ios_base::out | std::ios_base::app);
    file_RMS_out << time << ' ';
    double RMSD1, RMSD2, RMSD3, RMSD4, RMSD5;
    double NRMSD1, NRMSD2, NRMSD3, NRMSD4, NRMSD5;
    if (cwipiParamsObs == 0){
        ArrayXf MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU);
        ArrayXf MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU);
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                // NRMSD1 = std::sqrt(MSvals1.sum()/nb_oU)/obsArray.mean(); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << "\n";
                break;
            case 4 :
            case 5 :
            case 6 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << "\n";
                break;
            case 7 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                    MSvals3(i) = std::pow(sampMatrix.row(i+2*cwipiObsU).mean() - obsArray(i+2*cwipiObsU), 2);
                    MSvobs3(i) = std::pow(obsArray(i+2*cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD3 = std::sqrt(MSvals3.sum())/cwipiObsU; //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                file_RMS_out << "\n";
                break;
        }
    }
    else if (cwipiParamsObs == 1){
        ArrayXf MSvals4(cwipiObs-3*cwipiObsU);
        ArrayXf MSvobs4(cwipiObs-3*cwipiObsU);
        for (int i = 0; i < (cwipiObs-cwipiObsU); ++i){
            MSvals4(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
            MSvobs4(i) = std::pow(obsArray(i), 2) + epsilon;
        }
        RMSD4 = std::sqrt(MSvals4.sum())/cwipiObsU; //Root Mean Square Deviation
        NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(cwipiObs-cwipiObsU); //Normalized Root Mean Square Deviation
        file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
        file_RMS_out << "\n";
    }
    else if (cwipiParamsObs == 2){
        ArrayXf MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU), MSvals4(cwipiObs-3*cwipiObsU);
        ArrayXf MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU), MSvobs4(cwipiObs-3*cwipiObsU);
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                }
                for (int i = 0; i < (cwipiObs-cwipiObsU); ++i){
                    MSvals4(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs4(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD4 = std::sqrt(MSvals4.sum())/cwipiObsU; //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(cwipiObs-cwipiObsU); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                file_RMS_out << "\n";
                break;
            case 4 :
            case 5 :
            case 6 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                }
                for (int i = 0; i < (cwipiObs-2*cwipiObsU); ++i){
                    MSvals4(i) = std::pow(sampMatrix.row(i+2*cwipiObsU).mean() - obsArray(i+2*cwipiObsU), 2);
                    MSvobs4(i) = std::pow(obsArray(i+2*cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD4 = std::sqrt(MSvals4.sum())/(cwipiObs-2*cwipiObsU); //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(cwipiObs-2*cwipiObsU); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                file_RMS_out << "\n";
                break;
            case 7 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                    MSvals3(i) = std::pow(sampMatrix.row(i+2*cwipiObsU).mean() - obsArray(i+2*cwipiObsU), 2);
                    MSvobs3(i) = std::pow(obsArray(i+2*cwipiObsU), 2) + epsilon;
                }
                for (int i = 0; i < (cwipiObs-3*cwipiObsU); ++i){
                    MSvals4(i) = std::pow(sampMatrix.row(i+3*cwipiObsU).mean() - obsArray(i+3*cwipiObsU), 2);
                    MSvobs4(i) = std::pow(obsArray(i+3*cwipiObsU), 2) + epsilon; 
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD3 = std::sqrt(MSvals3.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD4 = std::sqrt(MSvals4.sum())/(cwipiObs-3*cwipiObsU); //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(cwipiObs-3*cwipiObsU); //Normalized Root Mean Square Deviation                 
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                file_RMS_out << "\n";
                break;
        }
    }
    else if (cwipiParamsObs == 3){
        ArrayXf MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU);
        ArrayXf MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU);
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD5 = -(sampMatrix.row(cwipiObsU).mean() - obsArray(cwipiObsU)); //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD5 = -(sampMatrix.row(cwipiObsU).mean() - obsArray(cwipiObsU))/ std::sqrt(std::pow(obsArray(cwipiObsU), 2) + epsilon); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                file_RMS_out << "\n";
                break;
            case 4 :
            case 5 :
            case 6 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD5 = -(sampMatrix.row(cwipiObsU).mean() - obsArray(cwipiObsU)); //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD5 = -(sampMatrix.row(2*cwipiObsU).mean() - obsArray(2*cwipiObsU))/ std::sqrt(std::pow(obsArray(2*cwipiObsU), 2) + epsilon); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                file_RMS_out << "\n";
                break;
            case 7 :
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                    MSvals2(i) = std::pow(sampMatrix.row(i+cwipiObsU).mean() - obsArray(i+cwipiObsU), 2);
                    MSvobs2(i) = std::pow(obsArray(i+cwipiObsU), 2) + epsilon;
                    MSvals3(i) = std::pow(sampMatrix.row(i+2*cwipiObsU).mean() - obsArray(i+2*cwipiObsU), 2);
                    MSvobs3(i) = std::pow(obsArray(i+2*cwipiObsU), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD2 = std::sqrt(MSvals2.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD3 = std::sqrt(MSvals3.sum())/cwipiObsU; //Root Mean Square Deviation
                RMSD5 = -(sampMatrix.row(cwipiObsU).mean() - obsArray(cwipiObsU)); //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                NRMSD5 = -(sampMatrix.row(3*cwipiObsU).mean() - obsArray(3*cwipiObsU))/ std::sqrt(std::pow(obsArray(3*cwipiObsU), 2) + epsilon); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                file_RMS_out << "\n";
                break;
        }
    file_RMS_out.close();
    }
}

//========================== Simple localization ============================
MatrixXf localisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, int nb_cells, int cwipiParams, int cwipiObs, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float clippingSwitch, float hyperlocSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath)
{
    //========= Definition of the cells affected by localisation
    std::string obs_file = stringRootPath + "/obs_coordinates.txt";
    std::ifstream observation;
    observation.open(obs_file);
    MatrixXf K(3*nb_cells+cwipiParams, cwipiObs);
    MatrixXf loc = MatrixXf::Zero(3*nb_cells+cwipiParams, cwipiObs);

    //========== We impose no localisation for the parameters to update ==========
    for (int i = 0; i < cwipiObs; ++i){
        for (int j = 0; j < cwipiParams; ++j){
            loc(3*nb_cells+j, i) = 1;
        }
    }

    double coordXdom, coordYdom, coordZdom;
	std::vector<double> coordXobs; std::vector<double> coordYobs; std::vector<double> coordZobs;

    const Foam::volVectorField& C = mesh.C();
    if (observation.is_open()){
        std::string line;
        // if (cwipiVerbose) std::cout << "Reading the probes for localisation..." << std::endl;
        int count = 0;
        while( std::getline(observation,line) )
        {
            std::stringstream ss(line);
            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX, ',');
            coordXobs.push_back(std::stod(coordX));

            std::getline(ss, coordY, ',');
            coordYobs.push_back(std::stod(coordY));

            std::getline(ss, coordZ, ',');
            coordZobs.push_back(std::stod(coordZ));
            ++count;
        }
    }
    else
    {
        std::cerr << "Couldn't open observation file for localisation.\n";
    }

    if (clippingSwitch == 1 && hyperlocSwitch!=1){
        std::string clipping_file = stringRootPath + "/clippingCells.txt";
        std::ifstream clippingCells;
        clippingCells.open(clipping_file);
        ArrayXf clippingCellID;
        if (clippingCells.is_open()){
            std::string line;
            int count = 0;
            while( std::getline(clippingCells, line) )
            {
                std::stringstream ss(line);
                std::string clippingCellID_str;
                std::getline(ss, clippingCellID_str, ',');
                if (count == 0){
                    clippingCellID.resize(std::stoi(clippingCellID_str));
                }
                else{
                    clippingCellID(count-1) = std::stoi(clippingCellID_str);
                }
            ++count;
            }
        }
        else
        {
            std::cerr<< "Couldn't open clippingCells.txt" << std::endl;
        }

        int clipIndex = 0;
        for(int i = 0; i < nb_cells; ++i){
            clipIndex = clippingCellID[i],
            coordXdom = C[clipIndex][0];
            coordYdom = C[clipIndex][1];
            coordZdom = C[clipIndex][2];
            for (int j = 0; j < cwipiObs; ++j){
                double deltaX = coordXobs[j] - coordXdom;
                double deltaY = coordYobs[j] - coordYdom;
                double deltaZ = coordZobs[j] - coordZdom;
                loc(i, j) = std::exp(-std::pow(deltaX, 2)/sigmaLocX);
                loc(i + nb_cells, j) = std::exp(-std::pow(deltaY, 2)/sigmaLocY);
                loc(i + 2*nb_cells, j) = std::exp(-std::pow(deltaZ, 2)/sigmaLocZ);
            }
        }
    }
    else{
        for(int i = 0; i < nb_cells; ++i){
            coordXdom = C[i][0];
            coordYdom = C[i][1];
            coordZdom = C[i][2];
            for (int j = 0; j < cwipiObs; ++j){
                double deltaX = coordXobs[j] - coordXdom;
                double deltaY = coordYobs[j] - coordYdom;
                double deltaZ = coordZobs[j] - coordZdom;
                loc(i, j) = std::exp(-std::pow(deltaX, 2)/sigmaLocX);
                loc(i + nb_cells, j) = std::exp(-std::pow(deltaY, 2)/sigmaLocY);
                loc(i + 2*nb_cells, j) = std::exp(-std::pow(deltaZ, 2)/sigmaLocZ);
            }
        }
    }
    // std::cout << "The localisation matrix is " << std::endl;
    // std::cout << loc << std::endl << "\n";
    //========== Multiplication of both matrices =========
    K = loc.cwiseProduct(KnoCorr);
    observation.close();
    return K;
}

//========================== Hyper-localization =============================
MatrixXf hyperlocalisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, const Eigen::Ref<const Eigen::ArrayXf>& clippingCells, int nb_clipCells, int clipCellIndex, int cwipiParams, int count_obs, int obs_hyperloc, double sigmaLocX, double sigmaLocY, double sigmaLocZ, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath)
{
    //========= Definition of the cells affected by localisation
    std::string obs_file = stringRootPath + "/obs_coordinates.txt";
    std::ifstream observation;
    observation.open(obs_file);
    MatrixXf K(3*nb_clipCells+cwipiParams, obs_hyperloc);
    MatrixXf loc = MatrixXf::Zero(3*nb_clipCells+cwipiParams, obs_hyperloc);

    //========== We impose no localisation for the parameters to update ==========
    for (int i = 0; i < obs_hyperloc; ++i){
        for (int j = 0; j < cwipiParams; ++j){
            loc(3*nb_clipCells+j, i) = 1;
        }
    }

    double coordXdom, coordYdom, coordZdom;
	std::vector<double> coordXobs; std::vector<double> coordYobs; std::vector<double> coordZobs;

    const Foam::volVectorField& C = mesh.C();
    if (observation.is_open()){
        std::string line;
        // if (cwipiVerbose) std::cout << "Reading the probes for localisation..." << std::endl;
        int count = 0;
        while( std::getline(observation,line) )
        {
            std::stringstream ss(line);
            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX, ',');
            coordXobs.push_back(std::stod(coordX));

            std::getline(ss, coordY, ',');
            coordYobs.push_back(std::stod(coordY));

            std::getline(ss, coordZ, ',');
            coordZobs.push_back(std::stod(coordZ));
            ++count;
        }
    }
    else
    {
        std::cerr << "Couldn't open observation file for localisation.\n";
    }

    int clipIndex = 0;
    for(int i = 0; i < nb_clipCells; ++i){
        clipIndex = clippingCells[clipCellIndex - nb_clipCells + i],
        coordXdom = C[clipIndex][0];
        coordYdom = C[clipIndex][1];
        coordZdom = C[clipIndex][2];
        double deltaX = coordXobs[count_obs] - coordXdom;
        double deltaY = coordYobs[count_obs] - coordYdom;
        double deltaZ = coordZobs[count_obs] - coordZdom;
        for (int j = 0; j < obs_hyperloc; ++j){
            loc(i, j) = std::exp(-std::pow(deltaX, 2)/sigmaLocX);
            loc(i + nb_clipCells, j) = std::exp(-std::pow(deltaY, 2)/sigmaLocY);
            loc(i + 2*nb_clipCells, j) = std::exp(-std::pow(deltaZ, 2)/sigmaLocZ);
        }
    }

    // std::cout << "The localisation matrix is " << std::endl;
    // std::cout << loc << std::endl << "\n";
    //========== Multiplication of both matrices =========
    K = loc.cwiseProduct(KnoCorr);
    observation.close();
    return K;
}

//========================== Inflation function =============================
MatrixXf inflation(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, int cwipiMembers, int nb_cells, int cwipiParams, double stateInfl, double paramsInfl, float typeInfl, float paramEstSwitch, float cwipiVerbose, std::string stringRootPath)
{
    //========== Initialisation of the variables ==========
    int nb_cells_cmpnt=nb_cells*3;
    float gaussample;

    MatrixXf stateMatrixInflated = stateMatrixUpt;

    // ========== Define random generator with Gaussian distribution ==========
    const double mean_state = 0;
    const double mean_params = 0;
    const double stddev_state = stateInfl;
    const double stddev_params = paramsInfl;
    std::normal_distribution<double> dist_state(mean_state, stddev_state);
    std::normal_distribution<double> dist_params(mean_params, stddev_params);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);


    if (stateInfl>0.00001){
    if(cwipiVerbose) std::cout << "Adding state inflation" << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {

            for (int i = 0; i < nb_cells_cmpnt; i++)
            {
                do
                {
                    gaussample=dist_state(generator);
                }while (gaussample<(mean_state-2*stddev_state) || gaussample>(mean_state+2*stddev_state));

                if (typeInfl){
                    stateMatrixInflated(i,j) = stateMatrixUpt.row(i).mean() + (1+stddev_state+(0.1*gaussample))*(stateMatrixUpt(i,j) - stateMatrixUpt.row(i).mean());  //Deterministic
                }
                else stateMatrixInflated(i,j) = stateMatrixUpt(i,j) + gaussample*stateMatrixUpt(i,j);     //Stochastic
            }
        }    
    }

    if (paramsInfl>0.00001){
    if(cwipiVerbose) std::cout << "Adding parameters inflation" << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {
            if (paramEstSwitch){

                for (int i = 0; i < cwipiParams; i++)
                {
                    do
                    {
                        gaussample=dist_params(generator);
                    }while (gaussample<(mean_params-2*stddev_params) || gaussample>(mean_params+2*stddev_params));
                    
                    if (typeInfl){
                        stateMatrixInflated(i + nb_cells_cmpnt, j) = stateMatrixUpt.row(i + nb_cells_cmpnt).mean() + (1+stddev_params+(0.1*gaussample))*(stateMatrixUpt(i + nb_cells_cmpnt, j) - stateMatrixUpt.row(i + nb_cells_cmpnt).mean()); //Deterministic
                    }
                    else stateMatrixInflated(i + nb_cells_cmpnt, j) = stateMatrixUpt(i + nb_cells_cmpnt, j) + gaussample*stateMatrixUpt(i + nb_cells_cmpnt, j);  //Stochastic
                }
            }
            else {
                std::string IntGen = stringRootPath + "/processor";
                std::string ensemble = "/UpdatedVariables";
                std::ifstream file;
                char member[50];
                std::string UpdatedVar;

                sprintf(member, "%i", j);
                UpdatedVar += IntGen; // Append the path
                UpdatedVar += member; // Append the member index
                UpdatedVar += ensemble; // Append file name
                file.open(UpdatedVar);

                MatrixXf parameters(cwipiParams, cwipiMembers);
                int count = 0;

                if (file.is_open())
                {
                    std::string line;
                    while( std::getline(file,line) )
                    {
                        std::stringstream ss(line);
                        std::string parNOUpt;
                        std::getline(ss, parNOUpt, ' ');
                        parameters(count, j) = std::stod(parNOUpt);

                        ++count;
                    }

                }
                else {
                    std::cerr << "Couldn't open any of the UpdatedVariables files.\n";
                }
                file.close();

                for (int i = 0; i < cwipiParams; i++)
                {
                    stateMatrixInflated(i + nb_cells_cmpnt, j) = parameters(i, j);    
                }
            }
        }
    }
   
    return stateMatrixInflated;
}

//========================== Randomize observations =========================
MatrixXf randomize_Obs(const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiMembers, int cwipiObs, int cwipiObsU, int cwipiParamsObs, double sigmaUserU, double sigmaUserp, double sigmaUserCf, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    float gaussample;

    //========== Initialisation of the matrices ==========
    MatrixXf obsMatrix(cwipiObs,cwipiMembers);
    MatrixXf err_mat(cwipiObs,cwipiMembers);


    //========== Define random generator with Gaussian distribution ==========
    const double mean_obs = 0;
    const double stddev_obsU = sigmaUserU;
    const double stddev_obsp = sigmaUserp;
    const double stddev_obsCf = sigmaUserCf;
    std::normal_distribution<double> dist_obsU(mean_obs, stddev_obsU);
    std::normal_distribution<double> dist_obsp(mean_obs, stddev_obsp);
    std::normal_distribution<double> dist_obsCf(mean_obs, stddev_obsCf);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    //========== Add Gaussian noise to create random velocity field and observation matrix ===========
    for (int j=0;j<cwipiMembers;j++)
    {
        if (cwipiParamsObs == 0){
            for (int i=0;i<cwipiObs;i++)
            {
                do
                {
                    gaussample=dist_obsU(generator);
                }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                err_mat(i,j)=gaussample;
                obsMatrix(i,j)=obsArray(i);
            }
        }
        else if (cwipiParamsObs == 1){
            for (int i=0;i<cwipiObs;i++)
            {
                do
                {
                    gaussample=dist_obsp(generator);
                }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                err_mat(i,j)=gaussample;
                obsMatrix(i,j)=obsArray(i);
            }
        }
        else if (cwipiParamsObs == 2){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i=0;i<cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample = dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j) = gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(cwipiObsU+i,j)=gaussample;
                        obsMatrix(cwipiObsU+i,j)=obsArray(cwipiObsU+i);
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i=0;i<2*cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - 2*cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(2*cwipiObsU+i,j)=gaussample;
                        obsMatrix(2*cwipiObsU+i,j)=obsArray(2*cwipiObsU+i);
                    }
                    break;
                case 7 :
                    for (int i=0;i<3*cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - 3*cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(3*cwipiObsU+i,j)=gaussample;
                        obsMatrix(3*cwipiObsU+i,j)=obsArray(3*cwipiObsU+i);
                    }
                    break;
            }
        }
        else if (cwipiParamsObs == 3){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i=0;i<cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample = dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j) = gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(cwipiObsU+i,j)=gaussample;
                        obsMatrix(cwipiObsU+i,j)=obsArray(cwipiObsU+i);
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i=0;i<2*cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - 2*cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(2*cwipiObsU+i,j)=gaussample;
                        obsMatrix(2*cwipiObsU+i,j)=obsArray(2*cwipiObsU+i);
                    }
                    break;
                case 7 :
                    for (int i=0;i<3*cwipiObsU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obsMatrix(i,j)=obsArray(i);
                    }
                    for (int i=0;i<(cwipiObs - 3*cwipiObsU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(3*cwipiObsU+i,j)=gaussample;
                        obsMatrix(3*cwipiObsU+i,j)=obsArray(3*cwipiObsU+i);
                    }
                    break;
            }
            //obsMatrix.col(j) = obsArray;
        }
    }
    if (typeInputs == 0) obsMatrix = obsMatrix + err_mat;
    else if (typeInputs == 1) obsMatrix = obsMatrix + obsMatrix.cwiseProduct(err_mat);
    // if (cwipiVerbose) std::cout << "obsMatrix = "<< obsMatrix << "\n" << std::endl;

    return obsMatrix;
}

//============= Calculate measurments error covariance matrix ===============
MatrixXf calculate_R(const Eigen::Ref<const Eigen::ArrayXf>& obsArray, int cwipiObs, int cwipiObsU, int cwipiParamsObs, double sigmaUserU, double sigmaUserp, double sigmaUserCf, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    MatrixXf R(cwipiObs,cwipiObs);
    R.setIdentity();

    // === Not Used ===
    // // **** Read the file with observation coordinates to apply to formula a+by^c ****
	// std::string obs_file = stringRootPath + "/obs_coordinates.txt";
    // std::ifstream observation;
    // observation.open(obs_file);
    // double coordXdom, coordYdom, coordZdom;
    // double coordXobs[cwipiObsU], coordYobs[cwipiObsU], coordZobs[cwipiObsU];

    // if (observation.is_open()){
    //     std::string line;
    //     if (cwipiVerbose == 1) std::cout << "Reading the probes for localisation..." << std::endl;
    //     int count = 0;
    //     while( std::getline(observation,line) )
    //     {
    //         std::stringstream ss(line);
    //         std::string coordX, coordY, coordZ;
    //         std::getline(ss, coordX, ',');
    //         coordXobs[count] = std::stod(coordX);
    //         std::getline(ss, coordY, ',');
    //         coordYobs[count] = std::stod(coordY);
    //         std::getline(ss, coordZ, ',');
    //         coordZobs[count] = std::stod(coordZ);

    //         ++count;
    //     }
    // }
    // else
    // {
    //     std::cerr << "Couldn't open observation file for localisation.\n";
    // }


    // if (cwipiParamsObs == 0){
    //     for(int i=0; i<cwipiObsU; i++){
    //         if(coordYobs[i]<1){
    //             R(i,i) = std::pow(sigmaUserUa+sigmaUserUb*std::pow(coordYobs[i],sigmaUserUc),2);
    //         }
    //         else{
    //             R(i,i) = std::pow(sigmaUserUa+sigmaUserUb*std::pow(2-coordYobs[i],sigmaUserUc),2);
    //         }
    //         R(i+cwipiObsU, i+cwipiObsU) = R(i,i);
    //         R(i+2*cwipiObsU, i+2*cwipiObsU) = R(i,i);
    //     }
    // } 

    // === Calculation of R ===
    if (cwipiParamsObs == 0){ 
        if (typeInputs == 0) R = R*sigmaUserU*sigmaUserU;
        else if (typeInputs == 1) R.diagonal() = obsArray.cwiseProduct(obsArray)*sigmaUserU*sigmaUserU;
    }
    else if (cwipiParamsObs == 1){
        if (typeInputs == 0) R = R*sigmaUserp*sigmaUserp;
        else if (typeInputs == 1) R.diagonal() = obsArray.cwiseProduct(obsArray)*sigmaUserp*sigmaUserp;
    }
    else if (cwipiParamsObs == 2){
        if (typeInputs == 0) R = R*sigmaUserp*sigmaUserp;
        else if (typeInputs == 1) R.diagonal() = obsArray.cwiseProduct(obsArray)*sigmaUserp*sigmaUserp;
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                R.block(0,0, cwipiObsU, cwipiObsU) = R.block(0,0, cwipiObsU, cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserp*sigmaUserp);
                break;
            case 4 :
            case 5 :
            case 6 :
                R.block(0,0, 2*cwipiObsU, 2*cwipiObsU) = R.block(0,0, 2*cwipiObsU, 2*cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserp*sigmaUserp);
                break;
            case 7 :
                R.block(0,0, 3*cwipiObsU, 3*cwipiObsU) = R.block(0,0, 3*cwipiObsU, 3*cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserp*sigmaUserp);
                break;
        }
    } 
    else if (cwipiParamsObs == 3){
        if (typeInputs == 0) R = R*sigmaUserCf*sigmaUserCf;
        else if (typeInputs == 1) R.diagonal() = obsArray.cwiseProduct(obsArray)*sigmaUserCf*sigmaUserCf;
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                R.block(0,0, cwipiObsU, cwipiObsU) = R.block(0,0, cwipiObsU, cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserCf*sigmaUserCf);
                break;
            case 4 :
            case 5 :
            case 6 :
                R.block(0,0, 2*cwipiObsU, 2*cwipiObsU) = R.block(0,0, 2*cwipiObsU, 2*cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserCf*sigmaUserCf);
                break;
            case 7 :
                R.block(0,0, 3*cwipiObsU, 3*cwipiObsU) = R.block(0,0, 3*cwipiObsU, 3*cwipiObsU)*sigmaUserU*sigmaUserU/(sigmaUserCf*sigmaUserCf);
                break;
        }
    }

    // if (cwipiVerbose) std::cout << "R :" << std::endl;
    // if (cwipiVerbose) std::cout << R << std::endl << "\n";

    return R;
}

//************============== Calculate Kalman gain =============*************
MatrixXf calculate_K(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const Eigen::Ref<const Eigen::MatrixXf>& obsMatrix, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::MatrixXf>& R, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiParams, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    int nb_cells_cmpnt=nb_cells*3;
    float gaussample;

    //========== Initialisation of the matrices ==========
    MatrixXf stateMatrix_upt(nb_cells_cmpnt + cwipiParams,cwipiMembers);
    ArrayXf stateMatrix_mean(nb_cells_cmpnt + cwipiParams);
    MatrixXf stateMatrix_norm(nb_cells_cmpnt + cwipiParams,cwipiMembers);

    ArrayXf obsMatrix_mean(cwipiObs);
    MatrixXf obsMatrix_norm(cwipiObs,cwipiMembers);

    ArrayXf sampMatrix_mean(cwipiObs);
    MatrixXf sampMatrix_norm(cwipiObs,cwipiMembers);

    MatrixXf K(3*nb_cells+cwipiParams, cwipiObs);
    MatrixXf K1;
    MatrixXf K2;
    MatrixXf K21;

    // ========== Ensemble means ===========
    for (int i=0;i<nb_cells_cmpnt+cwipiParams;i++)
    {
        stateMatrix_mean(i)=stateMatrix.row(i).mean();
    }

    for (int i=0;i<cwipiObs;i++)
    {
        obsMatrix_mean(i)=obsMatrix.row(i).mean();
        sampMatrix_mean(i)=sampMatrix.row(i).mean();
    }
    // if (cwipiVerbose) std::cout << "All means calculated" << std::endl;


    // ========== Normalized anomalies (Ue_mat-Ume)/sqrt(Ne-1) ==========
    for (int i=0;i<nb_cells_cmpnt+cwipiParams;i++)
    {
        for (int j=0;j<cwipiMembers;j++)
        {
            stateMatrix_norm(i,j)=(stateMatrix(i,j)-stateMatrix_mean(i))/std::sqrt(cwipiMembers-1);
        }
    }

    // if (cwipiVerbose) std::cout << "Anomaly state matrix calculated" << std::endl;

    // if (cwipiVerbose) std::cout << "sampMatrix = "<< sampMatrix << "\n" << std::endl;
    // if (cwipiVerbose) std::cout << "obsMatrix = "<< obsMatrix << "\n" << std::endl;

    for (int i=0;i<cwipiObs;i++)
    {
        // if (cwipiVerbose) std::cout << "i = " << i << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {
            // if (cwipiVerbose) std::cout << "j = " << j << std::endl;
            sampMatrix_norm(i,j)=(sampMatrix(i,j)-sampMatrix_mean(i))/std::sqrt(cwipiMembers-1);
            obsMatrix_norm(i,j)=(obsMatrix(i,j)-obsMatrix_mean(i))/std::sqrt(cwipiMembers-1);
        }
    }
    // if (cwipiVerbose) std::cout << "Anomaly matrices calculated" << std::endl;


    // ========== Kalman gain NAe_mat*trans(NAo_mat)*inv(NAo_mat*trans(NAo_mat)) ==========

    K21=sampMatrix_norm*(sampMatrix_norm.transpose())+R;

    // if (cwipiVerbose) std::cout << "Inverting the matrix for Kalman Gain" << std::endl;
    K2=K21.inverse();
    K1=stateMatrix_norm*(sampMatrix_norm.transpose());

    K=K1*K2;

    // if (cwipiVerbose) std::cout << "Kalman gain calculated (without localization) \n" << std::endl;

    
    return K;
}

//***********=============== Hyper EnKF ============================**********
MatrixXf EnKF_hyperloc(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsArray, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, double sigmaUserU, double sigmaUserp, double sigmaUserCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, float clippingSwitch, float hyperlocSwitch, int cwipiParams, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float paramEstSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath)
{
    if (cwipiVerbose) std::cout << "Entering the Hyperlocalized EnKF function" << std::endl;
    
    //========== Read clippingCells.txt file ===========
    std::string clipping_file = stringRootPath + "/clippingCells.txt";
 
    int count_clippingCells = 0;

    std::vector<std::vector<string>> clips_content;
	std::vector<string> clips_row;
	std::string clips_line, clips_word;

	std::fstream clips_stream (clipping_file, std::ios::in);
	if(clips_stream.is_open())
	{
		while(getline(clips_stream, clips_line))
		{
			clips_row.clear();
 
			std::stringstream clips_str(clips_line);
 
			while(getline(clips_str, clips_word, ','))
				clips_row.push_back(clips_word);
			clips_content.push_back(clips_row);
            count_clippingCells = count_clippingCells + 1; 
		}
	}
	else {
		std::cerr << "Couldn't open observation file.\n";
    }

    if (cwipiVerbose) std::cout << "clippingCells.txt file succesfully read" << std::endl;
    if (cwipiVerbose) std::cout << "count_clippingCells = " << count_clippingCells << std::endl;


    //======== Create an eigen array containing the ID of the clippings ========
    ArrayXf clippingCells(count_clippingCells); 
	for(int i=0; i<count_clippingCells; i++)
	{
        clippingCells(i) = std::stod(clips_content[i][0]);
    }

    if (cwipiVerbose) std::cout << "Array of the clippingCells created" << std::endl;
    // if (cwipiVerbose) std::cout << "clippingCells = " << clippingCells << std::endl;


    //======== Declare the good amount of observations =========
    MatrixXf temp_obsArray;
    MatrixXf temp_sampMatrix;
    int obs_hyperloc = 0;

    switch (velocityCase)
    {
        case 1:
        case 2:
        case 3:
            obs_hyperloc = 1;
            temp_obsArray = ArrayXf::Zero(obs_hyperloc); 
            temp_sampMatrix = MatrixXf::Zero(obs_hyperloc,cwipiMembers);
            
        break;
        case 4:
        case 5:
        case 6:
            obs_hyperloc = 2;
            temp_obsArray = ArrayXf::Zero(obs_hyperloc);
            temp_sampMatrix = MatrixXf::Zero(obs_hyperloc,cwipiMembers);
        break;
        case 7:
            obs_hyperloc = 3;
            temp_obsArray = ArrayXf::Zero(obs_hyperloc);
            temp_sampMatrix = MatrixXf::Zero(obs_hyperloc,cwipiMembers);
        break;
    }

    if (cwipiVerbose) std::cout << "Observations components declared" << std::endl;


    //======== Declaration of variables and arrays before the hyperlocalization loop =========
    int count_obs = 0;              // Count the observation use in the hyperlocalization
    int count_cellsClipsDone = 0;   // Count the number of cells of all the clipping already done in the hyperlocalization
    int count_nbClips = 0;          // Count the number of clippings 
    int cellIndex = 0;              // Retrieve cell index in the clipping file for the loop
    int IDclippingCell = 0;         // Current cell index of the clipping file


    MatrixXf stateMatrixUpt = stateMatrix;  // Save orginial state matrix to reconstruct it after hyperlocalization process

    int nb_clipCells = clippingCells(0);
    MatrixXf temp_state(nb_clipCells*3 + cwipiParams,cwipiMembers);
    MatrixXf temp_state_upt;
    MatrixXf temp_params = MatrixXf::Zero(cwipiParams,cwipiMembers);

    if (cwipiVerbose) std::cout << "Variables for the hyperlocalization loop declared " << std::endl;
    if (cwipiVerbose) std::cout << "Performing hyperlocalization loop ..." << std::endl;

    //========= Hyperlocalization loop goes over all the ID and perform EnKF once temp_state is filled with an entire clipping ===========
    for (int i=0; i<count_clippingCells; i++)
    {
        // if (cwipiVerbose) std::cout << " inside count_clippingCells loop i = " << i << std::endl;
        if (i==(nb_clipCells + count_cellsClipsDone)){
            count_nbClips = count_nbClips +1;
            // if (cwipiVerbose) std::cout << "Entered clipping if " << std::endl;
            // if (cwipiVerbose) std::cout << "count_cellsClipsDone value is = " << count_cellsClipsDone << std::endl;

            //** Fill the parameters now that we have all the values**
            for (int j = 0; j < cwipiMembers; j++){
                for (int k = 0; k < cwipiParams; k++){
                    temp_state(3*nb_clipCells + k, j) = stateMatrix(3*nb_cells + k, j);
                }
            }
            // if (cwipiVerbose) std::cout << "Params retrieved " << std::endl;
            // if (cwipiVerbose) std::cout << "temp_state = \n " << temp_state << std::endl;

            switch (velocityCase)
            {
                case 1:
                case 2:
                case 3:
                    temp_obsArray(0) = obsArray(count_obs);
                    for (int j = 0; j < cwipiMembers; j++){
                        temp_sampMatrix(0,j) = sampMatrix(count_obs,j);
                    }
                    break;
                case 4:
                case 5:
                case 6:
                    // if (cwipiVerbose) std::cout << "Inside case 6 " << temp_obsArray << std::endl;
                    temp_obsArray(0) = obsArray(count_obs);
                    temp_obsArray(1) = obsArray(cwipiObsU + count_obs);
                    // if (cwipiVerbose) std::cout << "temp_obsMatrix filled = \n " << temp_obsArray << std::endl;

                    for (int j = 0; j < cwipiMembers; j++){
                        temp_sampMatrix(0,j) = sampMatrix(count_obs,j);
                        temp_sampMatrix(1,j) = sampMatrix(cwipiObsU + count_obs,j);
                    }
                    // if (cwipiVerbose) std::cout << "temp_sampMatrix filled = \n " << temp_sampMatrix << std::endl;
                    break;
                case 7:
                    temp_obsArray(0) = obsArray(count_obs);
                    temp_obsArray(1) = obsArray(cwipiObsU + count_obs);
                    temp_obsArray(2) = obsArray(2*cwipiObsU + count_obs);
                    for (int j = 0; j < cwipiMembers; j++){
                        temp_sampMatrix(0,j) = sampMatrix(count_obs,j);
                        temp_sampMatrix(1,j) = sampMatrix(cwipiObsU + count_obs,j);
                        temp_sampMatrix(2,j) = sampMatrix(2*cwipiObsU + count_obs,j);
                    }
                    break;
            }

            // if (cwipiVerbose) std::cout << "Observations and samples selected according to velocityCase " << std::endl;
            // ** Extend observations to all the members via perturbations addition **
            MatrixXf temp_obsMatrix = randomize_Obs(temp_obsArray, cwipiMembers, obs_hyperloc, cwipiObsU, cwipiParamsObs, sigmaUserU, sigmaUserp, sigmaUserCf, typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, cwipiVerbose);
            // if (cwipiVerbose) std::cout << "temp_obsMatrix filled = \n " << temp_obsMatrix << std::endl;

            //** Calculate R **
            MatrixXf temp_R = calculate_R(temp_obsArray, obs_hyperloc, cwipiObsU, cwipiParamsObs, sigmaUserU, sigmaUserp, sigmaUserCf, typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, cwipiVerbose);
            // if (cwipiVerbose) std::cout << "temp_R filled = \n " << temp_R << std::endl;

            // //** Calculate the Kalman gain for temp_state **
            MatrixXf temp_K = calculate_K(temp_state, temp_obsMatrix, temp_sampMatrix, temp_R, cwipiMembers, nb_clipCells, obs_hyperloc, cwipiParams, cwipiVerbose);
            // if (cwipiVerbose) std::cout << "temp_K filled = \n " << temp_K << std::endl;
            
            // //** Perform covariance localization for temporary Kalman gain **
            temp_K = hyperlocalisation(temp_K, clippingCells, nb_clipCells, i, cwipiParams, count_obs, obs_hyperloc, sigmaLocX, sigmaLocY, sigmaLocZ, mesh, cwipiVerbose, stringRootPath);
            // if (cwipiVerbose) std::cout << "temp_K_loc filled = \n " << temp_K << std::endl;

            //** Update temporary state
            temp_state_upt = temp_state + temp_K*(temp_obsMatrix - temp_sampMatrix);
            // if (cwipiVerbose) std::cout << "temp_state_upt  = \n " << temp_state_upt << std::endl;

            //** Perform inflation on temp_sate withtout parameter inflation **
            temp_state_upt = inflation(temp_state_upt, cwipiMembers, nb_clipCells, cwipiParams, stateInfl, 0, typeInfl, paramEstSwitch, cwipiVerbose, stringRootPath);
            // if (cwipiVerbose) std::cout << "temp_state_upt_infl  = \n " << temp_state_upt << std::endl;

            //** Save updated params to perform an average on the values updated by each clipping **
            for (int j = 0; j < cwipiMembers; ++j){
                for (int k = 0; k < cwipiParams; ++k){
                    temp_params(k, j) = temp_params(k, j) + temp_state_upt(3*nb_clipCells + k, j);
                }
            }
            // if (cwipiVerbose) std::cout << "Params saved = " << temp_params << std::endl;

            //** Update the values giving the updated temp state **
            for (int j = 0; j < cwipiMembers; j++){
                for (int k = 0; k < nb_clipCells; k++){
                    // if (cwipiVerbose) std::cout << "k when updating the values = " << k << std::endl;
                    cellIndex = clippingCells(k + count_cellsClipsDone + 1);
                    stateMatrixUpt(cellIndex, j) = temp_state_upt(k, j);
                    stateMatrixUpt(nb_cells + cellIndex, j) = temp_state_upt(nb_clipCells + k, j);
                    stateMatrixUpt(2*nb_cells + cellIndex, j) = temp_state_upt(2*nb_clipCells + k, j);
                }
            }
            if (i<(count_clippingCells-1)){                                           // If clippings are not done
                count_cellsClipsDone = count_cellsClipsDone + nb_clipCells +1;        // Update the count of the clipping already performed
                nb_clipCells = clippingCells(i+1);                                    // Update the number of cells in next clipping
                count_obs = count_obs + 1;                                            // New observation index
                temp_state.resize(nb_clipCells*3 + cwipiParams,cwipiMembers);         // Resize temp_state over the new number of cells in the clipping
            }

        }
        else{
            IDclippingCell = clippingCells(i+1);
            //** Retrieve the correct values until reaching the next clipping **
            for (int j = 0; j < cwipiMembers; j++){
                temp_state(i - count_cellsClipsDone, j) = stateMatrixUpt(IDclippingCell, j);
                temp_state(nb_clipCells + i - count_cellsClipsDone, j) = stateMatrixUpt(nb_cells + IDclippingCell, j);
                temp_state(2*nb_clipCells + i - count_cellsClipsDone, j) = stateMatrixUpt(2*nb_cells + IDclippingCell, j);
            }
        }
    }

    if (cwipiVerbose) std::cout << "Hyperlocalization loop done, beginning to finalize the parameters  " << std::endl;
    temp_params = temp_params/count_nbClips;
    // if (cwipiVerbose) std::cout << "Params saved avg = " << temp_params << std::endl;
    for (int j = 0; j < cwipiMembers; j++){
        for (int k = 0; k < cwipiParams; k++){
            stateMatrixUpt(3*nb_cells + k, j) = temp_params(k, j);
        }
    }

    //** Perform inflation on stateMatrixUpt withtout state inflation **
    stateMatrixUpt = inflation(stateMatrixUpt, cwipiMembers, nb_cells, cwipiParams, 0, paramsInfl, typeInfl, paramEstSwitch, cwipiVerbose, stringRootPath);

    return stateMatrixUpt;
}

// ====************************ MAIN EnKF Function ************************===
MatrixXf mainEnKF(const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, const fvMesh& mesh, int cwipiMembers, int nb_cells, int cwipiObs, int cwipiObsU, double sigmaUserU, double sigmaUserp, double sigmaUserCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, float clippingSwitch, float hyperlocSwitch, int cwipiParams, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, double sigmaUserUc, float paramEstSwitch, float stateEstSwitch, float cwipiVerbose, std::string stringRootPath, int cwipiTimedObs, double obsTimeStep, double time, double epsilon)
{
    if(cwipiVerbose) std::cout << "** ENTERING MAIN EnKF FUNCTION **" << std::endl;
    //** The observation Matrix will be different depending on the high-fidelity observations **
    int inv_nb_cells = nb_cells;
    ArrayXf obsArray(cwipiObs); // By default our parameter is the velocity
    MatrixXf stateMatrixUpt = stateMatrix;

    // ================= Retrieve observations informations ================
    if (cwipiVerbose) std::cout << "Opening observation file function " << std::endl;
    if (cwipiTimedObs){
      int obsIndex = time / obsTimeStep - 1;
      obsArray = obs_Data_timed(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, obsIndex, cwipiParamsObs, cwipiVerbose, stringRootPath);
    }
    else obsArray = obs_Data(cwipiMembers, nb_cells, cwipiObs, cwipiObsU, cwipiParamsObs, cwipiVerbose, stringRootPath);

    // ================= Retrieve sampled informations ===================
    if (cwipiVerbose) std::cout << "Opening sampled data function " << std::endl;
    MatrixXf sampMatrix = samp_Data(cwipiMembers, cwipiObs, cwipiObsU, velocityCase, cwipiParamsObs, cwipiVerbose, stringRootPath);

    // ================= Begin calculation of the EnKF ===================
    tic();
    // **** Case using hyperlocalization 
    if (hyperlocSwitch == 1){
        stateMatrixUpt = EnKF_hyperloc(stateMatrix, obsArray, sampMatrix, cwipiMembers, nb_cells, cwipiObs, cwipiObsU, sigmaUserU, sigmaUserp, sigmaUserCf, sigmaLocX, sigmaLocY, sigmaLocZ, localSwitch, clippingSwitch, hyperlocSwitch, cwipiParams, cwipiParamsObs, stateInfl, paramsInfl, typeInfl, typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, paramEstSwitch, mesh, cwipiVerbose, stringRootPath);
    }
    // **** Normal cases
    else{
        // ** Extend observations to all the members via perturbations addition
        if (cwipiVerbose) std::cout << "Randomizing observations ... " << std::endl;
        MatrixXf obsMatrix = randomize_Obs(obsArray, cwipiMembers, cwipiObs, cwipiObsU, cwipiParamsObs, sigmaUserU, sigmaUserp, sigmaUserCf, typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, cwipiVerbose);
        
        // ** Calculate the measurment errors covariance matrix
        if (cwipiVerbose) std::cout << "Calculating measurements error covariance matrix ... " << std::endl;
        MatrixXf R = calculate_R(obsArray, cwipiObs, cwipiObsU, cwipiParamsObs, sigmaUserU, sigmaUserp, sigmaUserCf, typeInputs, velocityCase, sigmaUserUa, sigmaUserUb, sigmaUserUc, cwipiVerbose);

        // ** For basic localization clippings
        if (clippingSwitch == 1 && hyperlocSwitch!=1){
            if (cwipiVerbose) std::cout << "Creating the clipping ... " << std::endl;
            stateMatrixUpt = doClipping(stateMatrix, nb_cells, cwipiParams, cwipiMembers, cwipiVerbose, stringRootPath);
            if (cwipiVerbose) std::cout << "The number of values in each member is " << nb_cells*3 << std::endl;
            //** Writing values to a txt file **
            // print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateVector, "UMat_Clip");
        }

        // ** Calculation of the Kalman Gain
        if (cwipiVerbose) std::cout << "Calculating the Kalman gain ... " << std::endl;
        MatrixXf K = calculate_K(stateMatrixUpt, obsMatrix, sampMatrix, R, cwipiMembers, nb_cells, cwipiObs, cwipiParams, cwipiVerbose);
        // if (cwipiVerbose) std::cout << "Knoloc = " <<  K << std::endl;
        
        // ** If covariance localization needed
        if (localSwitch == 1 && hyperlocSwitch !=1){
            if (cwipiVerbose) std::cout << "Applying localization on the Kalman gain ... " << std::endl;
            K = localisation(K, nb_cells, cwipiParams, cwipiObs, sigmaLocX, sigmaLocY, sigmaLocZ, clippingSwitch, hyperlocSwitch, mesh, cwipiVerbose, stringRootPath);
            // if (cwipiVerbose) std::cout << "New K = " <<  K << std::endl;
        }
        else{
            if (cwipiVerbose) std::cout << "No localization implemented" << std::endl;
        }


        // ** Update of velocities Ue_mat + K_mat*(Uo_mat-H_mat*(Ue_mat))
        if (cwipiVerbose) std::cout << "Updating the values ... " << std::endl;
        stateMatrixUpt = stateMatrixUpt + K*(obsMatrix-sampMatrix);

        // ** If inflation needed
        if (cwipiVerbose) std::cout << "Applying inflation if needed ... " << std::endl;
        stateMatrixUpt = inflation(stateMatrixUpt, cwipiMembers, nb_cells, cwipiParams, stateInfl, paramsInfl, typeInfl, paramEstSwitch, cwipiVerbose, stringRootPath);

        // ** Reconstructing state if clipping
        if (clippingSwitch == 1 && hyperlocSwitch!=1){
            if (cwipiVerbose) std::cout << "Unclipping the domain ... " << std::endl;
            //** Writing updated values to a txt file **
            // print_matrix_on_txt(i, numberCwipiPhase, cwipiOutputNb, cwipiMembers, nb_cells, cwipiParams, time, stateMatrixUpt, "UMat_Clip_upt");
            stateMatrixUpt = undoClipping(stateMatrixUpt, stateMatrix, inv_nb_cells, nb_cells, cwipiParams, cwipiMembers, cwipiVerbose, stringRootPath);
        }
    }

    // ** If State estimation desactivated return to original values
    if(stateEstSwitch == 0){
        if (cwipiVerbose) std::cout << "State estimation desactivated, returning to original values ..." << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {
            for (int i = 0; i < inv_nb_cells*3; i++)
            {
                stateMatrixUpt(i,j) = stateMatrix(i,j);
            }
        }
    }

    toc();

    if (cwipiVerbose) std::cout << "Starting writing outputs from the EnKF..." << std::endl;
    EnKF_outputs(stateMatrixUpt, sampMatrix, obsArray, cwipiMembers, nb_cells, time, cwipiParams, cwipiObsU, cwipiObs, cwipiParamsObs, velocityCase, cwipiVerbose, epsilon);

    if (cwipiVerbose) std::cout << "** End of Kalman Filter **" << std::endl << "\n";

    return stateMatrixUpt;
}

//===========================================================================
void prepare_sendBack(double *sendValues, double *paramsSendValues, double *values, const Eigen::Ref<const Eigen::MatrixXf>& stateMatrixUpt, int nb_cells, int cwipiParams, int index, float cwipiVerbose)
{
    //========== Take the right column with the good values to send to the corresponding OF instance ==========

    if (cwipiVerbose) std::cout << "Sending back informations from EnKF to OF simulations ..." << std::endl;

    for (int i=0; i<nb_cells; i++)
    {   
        sendValues[3*i + 0] = stateMatrixUpt(i, index-1);
        sendValues[3*i + 1] = stateMatrixUpt(i + nb_cells, index-1);    
        sendValues[3*i + 2] = stateMatrixUpt(i + 2*nb_cells, index-1);    
    }

    for (int i = 0; i < cwipiParams; i++)
    {
        paramsSendValues[i] = stateMatrixUpt(i + 3*nb_cells, index-1);
    }
}

//===========================================================================
void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXf>& stateMatrix, std::string name)
{
    //========== Write the state vector in a txt file ==========
    
    if (phaseIndex >= numberCwipiPhase - cwipiOutputNb)
    {
    std::string stepNumber = std::to_string(time);
    std::string filename = "results/"+name+stepNumber;
    std::fstream file_out;
    
    file_out.open(filename, std::ios_base::out);
    
    if (!file_out.is_open()) 
    {
        std::cerr << "Failed to open " << filename << '\n';
    } 
    else 
    {
        for (int i=0; i<nb_cells*3+cwipiParams; i++)
        {
            for(int j=0;j<cwipiMembers;j++)
            {
                file_out << stateMatrix(i , j) << ' ';
            }
            file_out << "\n";
        }
    }
    // std::cout << "Done Writing Matrix on txt file!" << std::endl << "\n";
    }
}

// ************************************************************************* //
