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
		std::cout << "Elapsed time for the multiplication is " << (t_end - t_start).count()*1E-9 << " seconds\n";
	}
}

void toc() { tic(1); }

//===========================================================================

MatrixXf doClipping(Eigen::MatrixXf stateVector, Eigen::Ref<Eigen::MatrixXf> invStateVector, int&
nb_cells, int nb_p, int nb_e, float cwipiVerbose, std::string stringRootPath)
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

    invStateVector = stateVector;
    int nb_cellsClipping = clippingCellID.size();
    stateVector.resize(3*nb_cellsClipping + nb_p, nb_e);

    int clipIndex = 0;

    for (int i = 0; i < nb_e; ++i){
        for (int j = 0; j < nb_cellsClipping; ++j){
            clipIndex = clippingCellID(j);
            stateVector(j, i) = invStateVector(clipIndex, i);
            stateVector(nb_cellsClipping + j, i) = invStateVector(clipIndex + nb_cells, i);
            stateVector(2*nb_cellsClipping + j, i) = invStateVector(clipIndex + 2*nb_cells, i);
        }

        for (int j = 0; j < nb_p; ++j){
            stateVector(3*nb_cellsClipping+j, i) = invStateVector(3*nb_cells+j, i);
        }
    }

    nb_cells = nb_cellsClipping;
    if (cwipiVerbose) std::cout << "Leaving the doClipping function" << std::endl;

    return stateVector;
}

//===========================================================================
MatrixXf undoClipping(Eigen::MatrixXf UptMatrix, const Eigen::Ref<const Eigen::MatrixXf>& invStateVector, int inv_nb_cells, int& nb_cells, int nb_p, int nb_e, float cwipiVerbose, std::string stringRootPath)
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

    MatrixXf tempUptMatrix = UptMatrix;
    int nb_cellsClipping = clippingCellID.size();
    UptMatrix.resize(inv_nb_cells + nb_p, nb_e);
    UptMatrix = invStateVector;

    int clipIndex = 0;

    for (int i = 0; i < nb_e; ++i){
        for (int j = 0; j < nb_cellsClipping; ++j){
            clipIndex = clippingCellID(j);
            UptMatrix(clipIndex, i) = tempUptMatrix(j, i);
            UptMatrix(inv_nb_cells + clipIndex, i) = tempUptMatrix(j + nb_cellsClipping, i);
            UptMatrix(2*inv_nb_cells + clipIndex, i) = tempUptMatrix(j + 2*nb_cellsClipping, i);
        }

        for (int j = 0; j < nb_p; ++j){
            UptMatrix(3*inv_nb_cells+j, i) = tempUptMatrix(3*nb_cellsClipping+j, i);
        }
    }

    nb_cells = inv_nb_cells;

    if (cwipiVerbose) std::cout << "Leaving the undoClipping function" << std::endl;
    return UptMatrix;
}


//===========================================================================
MatrixXf localisation(const Eigen::Ref<const Eigen::MatrixXf>& KnoCorr, int nb_cells, int nb_p, int nb_o, int nb_oU, double sigmaLocX, double sigmaLocY, double sigmaLocZ, 
float clippingSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath)
{
    //========= Definition of the cells affected by localisation
    // (in the case of the plane channel ALL CELLS) =========//
    std::string obs_file = stringRootPath + "/obs_coordinates.txt";
    std::ifstream observation;
    observation.open(obs_file);
    MatrixXf K(3*nb_cells+nb_p, nb_o);
    MatrixXf loc = MatrixXf::Zero(3*nb_cells+nb_p, nb_o);

    //========== We impose no localisation for the parameters to update ==========
    for (int i = 0; i < nb_o; ++i){
        for (int j = 0; j < nb_p; ++j){
            loc(3*nb_cells+j, i) = 1;
        }
    }

    double coordXdom, coordYdom, coordZdom;
    double coordXobs[nb_oU], coordYobs[nb_oU], coordZobs[nb_oU];
    const Foam::volVectorField& C = mesh.C();
    if (observation.is_open()){
        std::string line;
        if (cwipiVerbose) std::cout << "Reading the probes for localisation..." << std::endl;
        int count = 0;
        while( std::getline(observation,line) )
        {
            std::stringstream ss(line);
            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX, ',');
            coordXobs[count] = std::stod(coordX);
            std::getline(ss, coordY, ',');
            coordYobs[count] = std::stod(coordY);
            std::getline(ss, coordZ, ',');
            coordZobs[count] = std::stod(coordZ);
            ++count;
        }
    }
    else
    {
        std::cerr << "Couldn't open observation file for localisation.\n";
    }

    if (clippingSwitch){
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
            for (int j = 0; j < nb_oU; ++j){
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
            for (int j = 0; j < nb_oU; ++j){
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
ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, int nb_oU, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{        
    //========== Not depending on time : Read the observation data from a file named "obs_field.txt" only and return the values
    // in a eigen array. The txt file has to contain 
    // 1) velocity values (first all values for X, secondly all values for Y and finally all values for Z) 
    // 2) pressure values 
    // 3) velocity and pressure values ==========
    
    //=========== Extract observation data from file ===========
    if (cwipiVerbose) std::cout << "Opening observation file function " << std::endl;

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

    ArrayXf obs_field(nb_o);
	for(int i=0; i<nb_o; i++)
	{
        obs_field(i) = std::stod(obs_content[i][0]);         
	}

    if (cwipiVerbose) std::cout << "observation field created" << std::endl << "\n";
    // if (cwipiVerbose == 1){
    // for (int i=0; i<nb_o; i++)
    // {
    //     std::cout << obs_field(i) << "   for loop " << i << std::endl << "\n" ;
    // }
    // std::cout << std::endl << "\n";
    // }

    return obs_field;
}

//===========================================================================
ArrayXf obs_Data_timed(int nb_e, int nb_cells, int nb_o, int nb_oU, int obsIndex, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{        
    //========== Depending on time : Read the observations data from a file named "obs_field.txt" only and return the values 
    // of the specific time step in a eigen array. The txt file has to contain velocity values with X Y and Z in the same column
    // (first all the X, then Y then Z). Number of column = number of DA phases ==========
    
    //=========== Extract observation data from file ===========
    if (cwipiVerbose) std::cout << "Opening observation file function" << std::endl;

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

    ArrayXf obs_field(nb_o); 
	for(int i=0; i<nb_o; i++)
	{
        obs_field(i) = std::stod(obs_content[i][obsIndex]);
    }
    
    if (cwipiVerbose) std::cout << "obs_field created" << std::endl << "\n";
    //if (cwipiVerbose == 1){
    // for (int i=0; i<nb_o; i++)
    // {
    //     std::cout << obs_field(i) << "   for loop " << i << std::endl << "\n" ;
    // }
    // std::cout << std::endl << "\n";
    // }

    return obs_field;
}

//===========================================================================

MatrixXf samp_Data(int nb_e, int nb_o, int nb_oU, int velocityCase, int cwipiParamsObs, float cwipiVerbose, std::string stringRootPath)
{
    //========== Read the sampled velocities from the files created by each OF instance. Return of matrix 
    // with all the sampled velocities organised in a X Y Z manner, each colunm corresponding to 
    // a member ========== 
    if (cwipiVerbose) std::cout << "Opening sampled data function " << std::endl;
    
    MatrixXf sampMatrix = MatrixXf::Zero(nb_o, nb_e); // By default our parameter is the velocity
    std::string IntGen = stringRootPath + "/cavity_testPar";
    std::ifstream file;
    char* IntGenChar = new char[1000];
    strcpy(IntGenChar, IntGen.c_str());

    for (int i = 0; i < nb_e; ++i)
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
            sprintf(ensemble, "/pInt%i", i);
            sprintf(member, "%i", i);
            pInt += IntGen; // Append the path
            pInt += member; // Append the member index
            pInt += ensemble; // Append file name
            if (cwipiVerbose) std::cout << "Opening " << pInt << std::endl;
            file.open(pInt);
        }
        else if (cwipiParamsObs == 2){
            std::string UpInt;
            sprintf(ensemble, "/UpInt%i", i);
            sprintf(member, "%i", i);
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
                            if (count < nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            break;
                        case 2 :
                            if (count >= nb_oU && count < 2*nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
                            break;
                        case 3 :
                            if (count >= 2*nb_oU) sampMatrix(count-2*nb_oU, i) = std::stod(IntValue);
                            break;
                        case 4 :
                            if (count < 2*nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            break;
                        case 5 :
                            if (count < nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 2*nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
                            break;
                        case 6 :
                            if (count >= nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
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
                            if (count < nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 3*nb_oU) sampMatrix(count-2*nb_oU, i) = std::stod(IntValue);
                            break;
                        case 2 :
                            if (count >= nb_oU && count < 2*nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
                            else if (count >= 3*nb_oU) sampMatrix(count-2*nb_oU, i) = std::stod(IntValue);
                            break;
                        case 3 :
                            if (count >= 2*nb_oU) sampMatrix(count-2*nb_oU, i) = std::stod(IntValue);
                            break;
                        case 4 :
                            if (count < 2*nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 3*nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
                            break;
                        case 5 :
                            if (count < nb_oU) sampMatrix(count, i) = std::stod(IntValue);
                            else if (count >= 2*nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
                            break;
                        case 6 :
                            if (count >= nb_oU) sampMatrix(count-nb_oU, i) = std::stod(IntValue);
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
    // for (int i=0; i< nb_o; i++)
    // { 
    //     for (int j=0; j<nb_e; j++)
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

MatrixXf KF(const Eigen::Ref<const Eigen::MatrixXf>& velo_field_mat, const Eigen::Ref<const Eigen::ArrayXf>& obs_field, const Eigen::Ref<const Eigen::MatrixXf>& proj_field_mat, 
int nb_e, int nb_cells, int nb_o, int nb_oU, double sigma_userU, double sigma_userp, double sigma_userCf, double sigmaLocX, double sigmaLocY, double sigmaLocZ, float localSwitch, 
float clippingSwitch, int nb_p, int cwipiParamsObs, double stateInfl, double paramsInfl, float typeInfl, int typeInputs, int velocityCase, double sigmaUserUa, double sigmaUserUb, 
double sigmaUserUc, float paramEstSwitch, const fvMesh& mesh, float cwipiVerbose, std::string stringRootPath)
{
    
    //========== Initialisation of the variables ==========
    int nb_cells_cmpnt=nb_cells*3;
    if (cwipiVerbose) std::cout << "The size of my state vector is " << nb_cells_cmpnt << std::endl;
    int nb_o_cmpnt = nb_o;
    float gaussample;

    //========== Initialisation of the matrices ==========
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

    MatrixXf K(3*nb_cells+nb_p, nb_o);
    MatrixXf K1;
    MatrixXf K2;
    MatrixXf K21;

    //========== Define random generator with Gaussian distribution ==========
    const double mean_obs = 0;
    const double stddev_obsU = sigma_userU;
    const double stddev_obsp = sigma_userp;
    const double stddev_obsCf = sigma_userCf;
    std::normal_distribution<double> dist_obsU(mean_obs, stddev_obsU);
    std::normal_distribution<double> dist_obsp(mean_obs, stddev_obsp);
    std::normal_distribution<double> dist_obsCf(mean_obs, stddev_obsCf);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    // // **** Read the file with observation coordinates to apply to formula a+by^c ****
	// std::string obs_file = stringRootPath + "/obs_coordinates.txt";
    // std::ifstream observation;
    // observation.open(obs_file);
    // double coordXdom, coordYdom, coordZdom;
    // double coordXobs[nb_oU], coordYobs[nb_oU], coordZobs[nb_oU];

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

    MatrixXf R(nb_o_cmpnt,nb_o_cmpnt);
    R.setIdentity();
    // if (cwipiParamsObs == 0){
    //     for(int i=0; i<nb_oU; i++){
    //         if(coordYobs[i]<1){
    //             R(i,i) = std::pow(sigmaUserUa+sigmaUserUb*std::pow(coordYobs[i],sigmaUserUc),2);
    //         }
    //         else{
    //             R(i,i) = std::pow(sigmaUserUa+sigmaUserUb*std::pow(2-coordYobs[i],sigmaUserUc),2);
    //         }
    //         R(i+nb_oU, i+nb_oU) = R(i,i);
    //         R(i+2*nb_oU, i+2*nb_oU) = R(i,i);
    //     }
    // } 
    if (cwipiParamsObs == 0){ 
        if (typeInputs == 0) R = R*stddev_obsU*stddev_obsU;
        else if (typeInputs == 1) R.diagonal() = obs_field.cwiseProduct(obs_field)*stddev_obsU*stddev_obsU;
    }
    else if (cwipiParamsObs == 1){
        if (typeInputs == 0) R = R*stddev_obsp*stddev_obsp;
        else if (typeInputs == 1) R.diagonal() = obs_field.cwiseProduct(obs_field)*stddev_obsp*stddev_obsp;
    }
    else if (cwipiParamsObs == 2){
        if (typeInputs == 0) R = R*stddev_obsp*stddev_obsp;
        else if (typeInputs == 1) R.diagonal() = obs_field.cwiseProduct(obs_field)*stddev_obsp*stddev_obsp;
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                R.block(0,0, nb_oU, nb_oU) = R.block(0,0, nb_oU, nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsp*stddev_obsp);
                break;
            case 4 :
            case 5 :
            case 6 :
                R.block(0,0, 2*nb_oU, 2*nb_oU) = R.block(0,0, 2*nb_oU, 2*nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsp*stddev_obsp);
                break;
            case 7 :
                R.block(0,0, 3*nb_oU, 3*nb_oU) = R.block(0,0, 3*nb_oU, 3*nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsp*stddev_obsp);
                break;
        }
    } 
    else if (cwipiParamsObs == 3){
        if (typeInputs == 0) R = R*stddev_obsCf*stddev_obsCf;
        else if (typeInputs == 1) R.diagonal() = obs_field.cwiseProduct(obs_field)*stddev_obsCf*stddev_obsCf;
        switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
                R.block(0,0, nb_oU, nb_oU) = R.block(0,0, nb_oU, nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsCf*stddev_obsCf);
                break;
            case 4 :
            case 5 :
            case 6 :
                R.block(0,0, 2*nb_oU, 2*nb_oU) = R.block(0,0, 2*nb_oU, 2*nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsCf*stddev_obsCf);
                break;
            case 7 :
                R.block(0,0, 3*nb_oU, 3*nb_oU) = R.block(0,0, 3*nb_oU, 3*nb_oU)*stddev_obsU*stddev_obsU/(stddev_obsCf*stddev_obsCf);
                break;
        }
    }

    // if (cwipiVerbose) std::cout << "R :" << std::endl;
    // if (cwipiVerbose) std::cout << R << std::endl << "\n";

    //========== Add Gaussian noise to create random velocity field and observation matrix ===========
    for (int j=0;j<nb_e;j++)
    {
        if (cwipiParamsObs == 0){
            for (int i=0;i<nb_o_cmpnt;i++)
            {
                do
                {
                    gaussample=dist_obsU(generator);
                }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                err_mat(i,j)=gaussample;
                obs_field_mat(i,j)=obs_field(i);
            }
        }
        else if (cwipiParamsObs == 1){
            for (int i=0;i<nb_o_cmpnt;i++)
            {
                do
                {
                    gaussample=dist_obsp(generator);
                }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                err_mat(i,j)=gaussample;
                obs_field_mat(i,j)=obs_field(i);
            }
        }
        else if (cwipiParamsObs == 2){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i=0;i<nb_oU;i++)
                    {
                        do
                        {
                            gaussample = dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j) = gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(nb_oU+i,j)=gaussample;
                        obs_field_mat(nb_oU+i,j)=obs_field(nb_oU+i);
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i=0;i<2*nb_oU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - 2*nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(2*nb_oU+i,j)=gaussample;
                        obs_field_mat(2*nb_oU+i,j)=obs_field(2*nb_oU+i);
                    }
                    break;
                case 7 :
                    for (int i=0;i<3*nb_oU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - 3*nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsp(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                        err_mat(3*nb_oU+i,j)=gaussample;
                        obs_field_mat(3*nb_oU+i,j)=obs_field(3*nb_oU+i);
                    }
                    break;
            }
        }
        else if (cwipiParamsObs == 3){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i=0;i<nb_oU;i++)
                    {
                        do
                        {
                            gaussample = dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j) = gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(nb_oU+i,j)=gaussample;
                        obs_field_mat(nb_oU+i,j)=obs_field(nb_oU+i);
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i=0;i<2*nb_oU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - 2*nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(2*nb_oU+i,j)=gaussample;
                        obs_field_mat(2*nb_oU+i,j)=obs_field(2*nb_oU+i);
                    }
                    break;
                case 7 :
                    for (int i=0;i<3*nb_oU;i++)
                    {
                        do
                        {
                            gaussample=dist_obsU(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                        err_mat(i,j)=gaussample;
                        obs_field_mat(i,j)=obs_field(i);
                    }
                    for (int i=0;i<(nb_o - 3*nb_oU);i++)
                    {
                        do
                        {
                            gaussample=dist_obsCf(generator);
                        }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                        err_mat(3*nb_oU+i,j)=gaussample;
                        obs_field_mat(3*nb_oU+i,j)=obs_field(3*nb_oU+i);
                    }
                    break;
            }
            //obs_field_mat.col(j) = obs_field;
        }
    }
    
    // std::cout << "obs_field_mat :" << std::endl;
    // std::cout << obs_field_mat << std::endl << "\n";
    obs_field_mat = obs_field_mat + err_mat;

    if (cwipiVerbose) std::cout << "Perturbations added" << std::endl;

    //proj_field_mat=H_mod*velo_field_mat;

    // std::cout << "velo_field_mat :" << std::endl;
    // std::cout << velo_field_mat << std::endl << "\n";

    // std::cout << "err_mat :" << std::endl;
    // std::cout << err_mat << std::endl << "\n";

    // ========== Ensemble means ===========

    for (int i=0;i<nb_cells_cmpnt+nb_p;i++)
    {
        velo_field_mean(i)=velo_field_mat.row(i).mean();
    }

    for (int i=0;i<nb_o_cmpnt;i++)
    {
        err_mat_mean(i)=err_mat.row(i).mean();
        obs_field_mean(i)=obs_field_mat.row(i).mean();
        proj_field_mean(i)=proj_field_mat.row(i).mean();
    }

    if (cwipiVerbose) std::cout << "All means calculated" << std::endl;
    
    // std::cout << "velo_field_mean :" << std::endl;
    // std::cout << velo_field_mean << std::endl << "\n";

    // std::cout << "err_mat_mean :" << std::endl;
    // std::cout << err_mat_mean << std::endl << "\n";

    // std::cout << "proj_field_mean :" << std::endl;
    // std::cout << proj_field_mean << std::endl << "\n";

    // ========== Normalized anomalies (Ue_mat-Ume)/sqrt(Ne-1) ==========
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

    if (cwipiVerbose) std::cout << "Anomaly matrices calculated" << std::endl;

    // std::cout << "velo_field_mat_norm :" << std::endl;
    // std::cout << velo_field_mat_norm << std::endl << "\n";

    // std::cout << "proj_field_mat_norm :" << std::endl;
    // std::cout << proj_field_mat_norm << std::endl << "\n";

    // ========== Kalman gain NAe_mat*trans(NAo_mat)*inv(NAo_mat*trans(NAo_mat)) ==========
    
    tic();
    K21=proj_field_mat_norm*(proj_field_mat_norm.transpose())+R;

    if (cwipiVerbose) std::cout << "Inverting the matrix for Kalman Gain" << std::endl;
    K2=K21.inverse();
    K1=velo_field_mat_norm*(proj_field_mat_norm.transpose());

    K=K1*K2;
    toc();

    if (cwipiVerbose) std::cout << "Kalman gain calculated (without localization)" << std::endl;

    // ========== We apply localization ==========

    if (localSwitch){
    K = localisation(K, nb_cells, nb_p, nb_o, nb_oU, sigmaLocX, sigmaLocY, sigmaLocZ, clippingSwitch, mesh, cwipiVerbose, stringRootPath);
    }
    else {
        if (cwipiVerbose) std::cout << "No localization implemented" << std::endl;
    }

    if (cwipiVerbose) std::cout << "Kalman gain calculated" << std::endl;
    
    // std::cout << "K2 :" << std::endl;
    // std::cout << K2 << std::endl << "\n";

    // std::cout << "K1 :" << std::endl;
    // std::cout << K1 << std::endl << "\n";

    // if (cwipiVerbose) std::cout << "K :" << std::endl;
    // if (cwipiVerbose) std::cout << K << std::endl << "\n";
    
    // std::cout << "ve1 :" << std::endl;
    // std::cout << velo_field_mat_norm.col(j) << std::endl << "\n";

    // std::cout << "proj1 :" << std::endl;
    // std::cout << (proj_field_mat_norm.transpose()).row(j) << std::endl << "\n";

    // ========== Update of velocities Ue_mat + K_mat*(Uo_mat-H_mat*(Ue_mat)) ==========

    // ========== Define random generator with Gaussian distribution ==========
    const double mean_state = stateInfl;
    const double mean_params = paramsInfl;
    const double stddev_state = (stateInfl-1)*0.1;
    const double stddev_params = (paramsInfl-1)*0.1;
    std::normal_distribution<double> dist_state(mean_state, stddev_state);
    std::normal_distribution<double> dist_params(mean_params, stddev_params);

    velo_field_mat_upt = velo_field_mat + K*(obs_field_mat-proj_field_mat);

    if (cwipiVerbose) std::cout << "The updated vector calculated. Adding inflation..." << std::endl;

    for (int j=0;j<nb_e;j++)
    {
        do
        {
            gaussample=dist_state(generator);
        }while (gaussample<(mean_state-2*stddev_state) || gaussample>(mean_state+2*stddev_state));

        for (int i = 0; i < nb_cells_cmpnt; i++)
        {
            if (typeInfl){
                velo_field_mat_upt(i,j) = velo_field_mat_upt.row(i).mean() + gaussample*(velo_field_mat_upt(i,j) - velo_field_mat_upt.row(i).mean());  //Deterministic
            }
            else velo_field_mat_upt(i,j) = velo_field_mat_upt(i,j) + gaussample*velo_field_mat_upt(i,j);     //Stochastic
        }

        if (paramEstSwitch){
            do
            {
                gaussample=dist_params(generator);
            }while (gaussample<(mean_params-2*stddev_params) || gaussample>(mean_params+2*stddev_params));

            for (int i = 0; i < nb_p; i++)
            {
                if (typeInfl){
                    velo_field_mat_upt(i + nb_cells_cmpnt, j) = velo_field_mat_upt.row(i + nb_cells_cmpnt).mean() + gaussample*(velo_field_mat_upt(i + nb_cells_cmpnt, j) - velo_field_mat_upt.row(i + nb_cells_cmpnt).mean());
                }
                else velo_field_mat_upt(i + nb_cells_cmpnt, j) = velo_field_mat_upt(i + nb_cells_cmpnt, j) + gaussample*velo_field_mat_upt(i + nb_cells_cmpnt, j);  //Stochastic
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

            MatrixXf parameters(nb_p, nb_e);
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

            for (int i = 0; i < nb_p; i++)
            {
                velo_field_mat_upt(i + nb_cells_cmpnt, j) = parameters(i, j);    
            }
        }
    }
    
    if (cwipiVerbose) std::cout << "End of Kalman Filter" << std::endl << "\n";
    
    return velo_field_mat_upt;
}

//===========================================================================

void KF_output(double *sendValues, double *paramsSendValues, double *values, const Eigen::Ref<const Eigen::MatrixXf>& UptMatrix, const Eigen::Ref<const Eigen::MatrixXf>& sampMatrix, const Eigen::Ref<const Eigen::ArrayXf>& obsMatrix, 
int nb_e, int nb_cells, double time, int nb_p, int nb_oU, int nb_o, int cwipiParamsObs, int velocityCase, int index, int subdomains, int mainsubdom, double epsilon, float cwipiVerbose)
{
    //========== Take the right column with the good values to send to the corresponding OF instance ==========

    if (cwipiVerbose) std::cout << "Starting writing outputs from the EnKF..." << std::endl;

    for (int i=0; i<nb_cells; i++)
    {   
        sendValues[3*i + 0] = UptMatrix(i, index-1);
        sendValues[3*i + 1] = UptMatrix(i + nb_cells, index-1);    
        sendValues[3*i + 2] = UptMatrix(i + 2*nb_cells, index-1);    
    }

    for (int i = 0; i < nb_p; i++)
    {
        paramsSendValues[i] = UptMatrix(i + 3*nb_cells, index-1);
    }
    
    //========== We do an output of all optimized coefficients to evaluate convergence ==========//

    if (index == 1){
        std::fstream file_coeffs_out;
        std::string filename_coeffs = "results/UpdatedCoefficients";
        file_coeffs_out.open(filename_coeffs, std::ios_base::out | std::ios_base::app);
        file_coeffs_out << time << ' ';
        for (int i = 0; i < nb_p; ++i){
            double average = UptMatrix.row(nb_cells*3 + i).mean();
            double std_dev_temp = 0;
            for (int j = 0; j < nb_e; ++j){
                std_dev_temp = std_dev_temp + std::pow((UptMatrix(nb_cells*3 + i, j) - average), 2);
            }
            double std_dev = std::sqrt(std_dev_temp/nb_e);
            file_coeffs_out << average << ' ' << std_dev << ' ';
        }
        file_coeffs_out << "\n";
        file_coeffs_out.close();

        std::fstream file_RMS_out;
        std::string filename_RMS = "results/NRMSD";
        file_RMS_out.open(filename_RMS, std::ios_base::out | std::ios_base::app);
        file_RMS_out << time << ' ';
        ArrayXf MSvals1(nb_oU), MSvals2(nb_oU), MSvals3(nb_oU), MSvals4(nb_o-nb_oU);
        ArrayXf MSvobs1(nb_oU), MSvobs2(nb_oU), MSvobs3(nb_oU), MSvobs4(nb_o-nb_oU);
        double RMSD1, RMSD2, RMSD3, RMSD4, RMSD5;
        double NRMSD1, NRMSD2, NRMSD3, NRMSD4, NRMSD5;
        if (cwipiParamsObs == 0){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 7 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                        MSvals3(i) = std::pow(sampMatrix.row(i+2*nb_oU).mean() - obsMatrix(i+2*nb_oU), 2);
                        MSvobs3(i) = std::pow(obsMatrix(i+2*nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD3 = std::sqrt(MSvals3.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                    file_RMS_out << "\n";
                    break;
            }
        }
        else if (cwipiParamsObs == 1){
            for (int i = 0; i < (nb_o-nb_oU); ++i){
                MSvals4(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                MSvobs4(i) = std::pow(obsMatrix(i), 2) + epsilon;
            }
            RMSD4 = std::sqrt(MSvals4.sum())/nb_oU; //Root Mean Square Deviation
            NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(nb_o-nb_oU); //Normalized Root Mean Square Deviation
            file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
            file_RMS_out << "\n";
        }
        else if (cwipiParamsObs == 2){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                    }
                    for (int i = 0; i < (nb_o-nb_oU); ++i){
                        MSvals4(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs4(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD4 = std::sqrt(MSvals4.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(nb_o-nb_oU); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                    }
                    for (int i = 0; i < (nb_o-nb_oU); ++i){
                        MSvals4(i) = std::pow(sampMatrix.row(i+2*nb_oU).mean() - obsMatrix(i+2*nb_oU), 2);
                        MSvobs4(i) = std::pow(obsMatrix(i+2*nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD4 = std::sqrt(MSvals4.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(nb_o-nb_oU); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 7 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                        MSvals3(i) = std::pow(sampMatrix.row(i+2*nb_oU).mean() - obsMatrix(i+2*nb_oU), 2);
                        MSvobs3(i) = std::pow(obsMatrix(i+2*nb_oU), 2) + epsilon;
                    }
                    for (int i = 0; i < (nb_o-nb_oU); ++i){
                        MSvals4(i) = std::pow(sampMatrix.row(i+3*nb_oU).mean() - obsMatrix(i+3*nb_oU), 2);
                        MSvobs4(i) = std::pow(obsMatrix(i+3*nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD3 = std::sqrt(MSvals3.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD4 = std::sqrt(MSvals4.sum())/nb_oU; //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/(nb_o-nb_oU); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                    file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
                    file_RMS_out << "\n";
                    break;
            }
        }
        else if (cwipiParamsObs == 3){
            switch (velocityCase){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD5 = -(sampMatrix.row(nb_oU).mean() - obsMatrix(nb_oU)); //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD5 = -(sampMatrix.row(nb_oU).mean() - obsMatrix(nb_oU))/ std::sqrt(std::pow(obsMatrix(nb_oU), 2) + epsilon); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD5 = -(sampMatrix.row(nb_oU).mean() - obsMatrix(nb_oU)); //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD5 = -(sampMatrix.row(2*nb_oU).mean() - obsMatrix(2*nb_oU))/ std::sqrt(std::pow(obsMatrix(2*nb_oU), 2) + epsilon); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                    file_RMS_out << "\n";
                    break;
                case 7 :
                    for (int i = 0; i < nb_oU; ++i){
                        MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsMatrix(i), 2);
                        MSvobs1(i) = std::pow(obsMatrix(i), 2) + epsilon;
                        MSvals2(i) = std::pow(sampMatrix.row(i+nb_oU).mean() - obsMatrix(i+nb_oU), 2);
                        MSvobs2(i) = std::pow(obsMatrix(i+nb_oU), 2) + epsilon;
                        MSvals3(i) = std::pow(sampMatrix.row(i+2*nb_oU).mean() - obsMatrix(i+2*nb_oU), 2);
                        MSvobs3(i) = std::pow(obsMatrix(i+2*nb_oU), 2) + epsilon;
                    }
                    RMSD1 = std::sqrt(MSvals1.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD2 = std::sqrt(MSvals2.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD3 = std::sqrt(MSvals3.sum())/nb_oU; //Root Mean Square Deviation
                    RMSD5 = -(sampMatrix.row(nb_oU).mean() - obsMatrix(nb_oU)); //Root Mean Square Deviation
                    NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD2 = std::sqrt(MSvals2.sum()/MSvobs2.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD3 = std::sqrt(MSvals3.sum()/MSvobs3.sum())/nb_oU; //Normalized Root Mean Square Deviation
                    NRMSD5 = -(sampMatrix.row(3*nb_oU).mean() - obsMatrix(3*nb_oU))/ std::sqrt(std::pow(obsMatrix(3*nb_oU), 2) + epsilon); //Normalized Root Mean Square Deviation
                    file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                    file_RMS_out << RMSD2 << ' ' << NRMSD2 << ' ';
                    file_RMS_out << RMSD3 << ' ' << NRMSD3 << ' ';
                    file_RMS_out << RMSD5 << ' ' << NRMSD5 << ' ';
                    file_RMS_out << "\n";
                    break;
            }
        }
    file_RMS_out.close();
    }
}

//===========================================================================

void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXf>& stateVector, std::string name)
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
                file_out << stateVector(i , j) << ' ';
            }
            file_out << "\n";
        }
    }
    // std::cout << "Done Writing Matrix on txt file!" << std::endl << "\n";
    }
}

// ************************************************************************* //
