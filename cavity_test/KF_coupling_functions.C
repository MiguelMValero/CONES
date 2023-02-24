//#include "fvCFD.H"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <Eigen/Dense>
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

// To heapify a subtree rooted with node i which is an index in arr[]. n is size of heap
void heapify(int arr[], int arr_index[], int n, int i)
{
	int largest = i; // Initialize largest as root
	int l = 2 * i + 1; // left = 2*i + 1
	int r = 2 * i + 2; // right = 2*i + 2

	// If left child is larger than root
	if (l < n && arr[l] > arr[largest])
		largest = l;

	// If right child is larger than largest so far
	if (r < n && arr[r] > arr[largest])
		largest = r;

	// If largest is not root
	if (largest != i) {
		std::swap(arr[i], arr[largest]);
		std::swap(arr_index[i], arr_index[largest]);

		// Recursively heapify the affected sub-tree
		heapify(arr, arr_index, n, largest);
	}
}

// main function to do heap sort
void heapSort(int arr[], int arr_index[], int n)
{
	// Build heap (rearrange array)
	for (int i = n / 2 - 1; i >= 0; i--)
		heapify(arr, arr_index, n, i);

	// One by one extract an element from heap
	for (int i = n - 1; i >= 0; i--) {
		// Move current root to end
		std::swap(arr[0], arr[i]);
		std::swap(arr_index[0], arr_index[i]);

		// call max heapify on the reduced heap
		heapify(arr, arr_index, i, 0);
	}
}

/* A utility function to print array of size n */
void printArray(int arr[], int arr_index[],  int n)
{
	for (int i = 0; i < n; ++i)
		std::cout << arr[i] << " ";
	std::cout << "\n";

	for (int i = 0; i < n; ++i)
		std::cout << arr_index[i] << " ";
	std::cout << "\n";
}

void define_mesh(double* pointCoords, int* face_index, int* cell_to_face_connectivity, int* face_connectivity_index, int* face_connectivity, int c2fconnec_size, int fconnec_size, const fvMesh& mesh, char cl_coupling_name[20], int cwipiVerbose)
{
    //=== Coupling mesh definition ===
    
    if (cwipiVerbose == 1) Info << "Definition du maillage KF" << nl << endl;

    int nCells = mesh.nCells();
    int nPoints = mesh.nPoints();
    int nFaces = mesh.nFaces();

    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }
    
    Info << "CC 2" << nl << endl;

    // connecIdx[0]=0;
    // forAll(mesh.cells(),i)
    // {
    //     connecIdx[i+1]=connecIdx[i]+8;
    // }

    Info << "CC 3" << nl << endl;

    // forAll(mesh.cells(),i)
    // {
    //     forAll(mesh.cellShapes()[i],j)
    //     {
    //         connec[8*i+j]=mesh.cellShapes()[i][j]+1;
    //     }
    // }

    Info << "After declaration of numbers" << endl;

	Info << "nCells = " << nCells << endl;
	Info << "nPoints = " << nPoints << endl;
	Info << "nFaces = " << nFaces << endl;
	Info << "c2fconnec_size = " << c2fconnec_size << endl;
	Info << "fconnec_size = " << fconnec_size << endl;

    face_index[0]=0;
    face_connectivity_index[0]=0;

	Info << "After init of arrays" << endl;

    int *face_index_temp = new int[nCells];
	// Info << "Size of face_index_temp = " << sizeof(face_index_temp)/sizeof(face_index_temp[0]) << endl;

	int *face_index_temp_IDs = new int[nCells];
	// Info << "Size of face_index_temp_IDs = " << sizeof(face_index_temp_IDs)/sizeof(face_index_temp_IDs[0]) << endl;

	int *face_connectivity_index_temp = new int[nFaces];
	// Info << "Size of face_connectivity_index_temp = " << sizeof(face_connectivity_index_temp)/sizeof(face_connectivity_index_temp[0]) << endl;

	int *face_connectivity_index_temp_IDs = new int[nFaces];
	// Info << "Size of face_connectivity_index_temp_IDs = " << sizeof(face_connectivity_index_temp_IDs)/sizeof(face_connectivity_index_temp_IDs[0]) << endl;

    //int* face_connectivity_index_NoCORR = new int[mesh.nFaces()+1];
	Info << "We are here" << endl;

	forAll(mesh.cells(), i)
    {
        face_index_temp[i] = mesh.cells()[i].size();
        face_index_temp_IDs[i] = i;
    }

	Info << "We are here 1" << endl;

	//heapSort(face_index_temp, face_index_temp_IDs, nCells);

	Info << "We are here 2" << endl;

	forAll(mesh.cells(), i)
    {
        face_index[i+1] = face_index[i] + face_index_temp[i];
    }

	Info << "We are here 3" << endl;

    int faceOwner = 0;
    int faceNeighbour = 0;
    int faceID = 0;
	int faces_count = 0;
	for(int i = 0; i < nCells; i++){
        for (int j = 0; j < face_index_temp[i] ; j++) {
            // Info << "The size of the cell " << i << " is : " <<  face_index_temp[i] << endl;
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
            //if (cell_to_face_connectivity[faces_count+j] < 0) std::cout << cell_to_face_connectivity[faces_count+j] << " ";
        }
        faces_count = faces_count + face_index_temp[i];
    }

	Info << "We are here 4" << endl;

	//Faces calculations

	forAll(mesh.faces(), i)
    {
        face_connectivity_index_temp[i] = mesh.faces()[i].size();
        face_connectivity_index_temp_IDs[i] = i;
    }

	Info << "We are here 5" << endl;

	//heapSort(face_connectivity_index_temp, face_connectivity_index_temp_IDs, nFaces);

	Info << "We are here 6" << endl;

    //face_connectivity_index_NoCORR[0]=0;
    forAll(mesh.faces(), i)
    {    
        face_connectivity_index[i+1] = face_connectivity_index[i] + face_connectivity_index_temp[i];
    }

	Info << "We are here 7" << endl;

	int vertices_cout = 0;
    for(int i = 0; i < nFaces; i++)
    {
		// Info << "i = " << i << endl;
        for (int j = 0; j < face_connectivity_index_temp[i] ; j++)
        {
        face_connectivity[vertices_cout+j] = mesh.faces()[face_connectivity_index_temp_IDs[i]][j] + 1;
        }
		vertices_cout = vertices_cout + face_connectivity_index_temp[i];
    }

	Info << "We are here 8" << endl;

	Info << "The number of vertices of the first face_connectivity_index_temp_IDs is : " << mesh.faces()[face_connectivity_index_temp_IDs[0]] << endl;
	Info << "We can find the IDs of thoses vertices in the 6 firsts values of face_connectivity : " << face_connectivity[0] << " " << face_connectivity[1] << " " << face_connectivity[2] << " " << face_connectivity[3] << " " << endl;
	Info << "The number of faces of the lasts face_connectivity_index_temp_IDs is : " << mesh.faces()[face_connectivity_index_temp_IDs[nFaces-1]] << endl;
	Info << "We can find the IDs of thoses vertices in the 18 lasts values of face_connectivity : " << face_connectivity[fconnec_size-7] << " " << face_connectivity[fconnec_size-6] << " " << face_connectivity[fconnec_size-5] << " " << face_connectivity[fconnec_size-4] << " " << face_connectivity[fconnec_size-3]<< " " << face_connectivity[fconnec_size-2] << " " << face_connectivity[fconnec_size-1]<< endl;
    
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
 
    if (cwipiVerbose == 1) Info << "Fin de definition du maillage KF" << nl << endl;

    cwipi_locate(cl_coupling_name);
    

    delete[] face_index_temp;
	delete[] face_index_temp_IDs;
	delete[] face_connectivity_index_temp;
	delete[] face_connectivity_index_temp_IDs;
}

//===========================================================================

MatrixXf repro_Members(const int nb_e, int nb_cells, double time, double* values, double* paramsValues, int nb_p, int cwipiVerbose)
{    
    //=== NOT USED ANYMORE : Function to reproduce velocity fields in the fake EnKF ===
    
    MatrixXf velo_field_mat = MatrixXf::Zero(nb_cells*3 + nb_p, nb_e);

    // std::cout << "values[0] =" << values[0] << std::endl;
    std::cout << "Reproduction phase"<< std::endl;

    float gaussample;

    //==== Define random generator with Gaussian distribution ====
    const double mean = 0;
    const double stddev = 1;
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

ArrayXf obs_Data(int nb_e, int nb_cells, int nb_o, int nb_oU, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath)
{        
    //=== Not depending on time : Read the observation data from a file named "obs_field.txt" only and return the values
    // in a eigen array. The txt file has to contain 
    // 1) velocity values (first all values for X, secondly all values for Y and finally all values for Z) 
    // 2) pressure values 
    // 3) velocity and pressure values ===
    
    //==== Extract observation data from file ====
    if (cwipiVerbose == 1) std::cout << "Opening observation file" << std::endl;
    
    std::string obs_file= stringRootPath + "/obs_field.txt";
    if (cwipiVerbose == 1) std::cout << "obs_file path = " << stringRootPath << std::endl;
 
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

    ArrayXf obs_field(nb_o); // By default our parameter is the velocity
    // ArrayXf coord_obs_field(nb_o*3);

    if (cwipiVerbose == 1) std::cout << "Obs file is open, begin to fill arrays" << std::endl;
 
    if (cwipiParamsObs == 0){
	    for(int i=0; i<nb_oU; i++)
	    {
        obs_field(i)=stod(obs_content[i][0]);
        obs_field(i+nb_oU)=stod(obs_content[i + nb_oU][0]);

        }
    }
    else if (cwipiParamsObs == 1){
    	for(int i=0; i<nb_o; i++)
	    {
        obs_field(i)=stod(obs_content[i][0]);

        }
    }
    else if (cwipiParamsObs == 2){
        for(int i=0; i<nb_oU; i++)
        {
        obs_field(i)=stod(obs_content[i][0]);
        obs_field(i+nb_oU)=stod(obs_content[i + nb_oU][0]);
        obs_field(i+2*nb_oU)=stod(obs_content[i + 2*nb_oU][0]);
        }
        for(int i=0; i<(nb_o - 3*nb_oU); i++)
        {
        obs_field(i+3*nb_oU)=stod(obs_content[i + 3*nb_oU][0]);
        }
       
	}

    if (cwipiVerbose == 1) std::cout << "obs_field created" << std::endl << "\n";
    if (cwipiVerbose == 1){
    for (int i=0; i<nb_o; i++)
    {
        std::cout << obs_field(i) << "   for loop " << i << std::endl << "\n" ;
    }
    std::cout << std::endl << "\n";
    }

    return obs_field;
}

//===========================================================================

ArrayXf obs_Data_timed(int nb_e, int nb_cells, int nb_o, int nb_oU, int i, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath)
{        
    //=== Depending on time : Read the observations data from a file named "obs_field.txt" only and return the values 
    // of the specific time step in a eigen array. The txt file has to contain velocity values with X Y and Z in the same column
    // (first all the X, then Y then Z). Number of column = number of DA phases ===
    
    //==== Extract observation data from file ====
    if (cwipiVerbose == 1) std::cout << "Opening observation file" << std::endl;

    std::string obs_file= stringRootPath + "/obs_field.txt";
    if (cwipiVerbose == 1) std::cout << "obs_file path = " << stringRootPath << std::endl;

    int timeIndex = i;

    if (cwipiVerbose == 1) std::cout << "timeIndex = " << timeIndex << std::endl << "\n";
    //cout << obs_file;
 
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
		std::cout<<"Could not open the file\n";

    ArrayXf obs_field(nb_o); // By default our parameter is the velocity

    // ArrayXf coord_obs_field(nb_o*3);

    if (cwipiVerbose == 1) std::cout << "Obs file is open, begin to fill arrays" << endl;
    // for (int i=0; i<nb_o*3; i++)
    // {
    //     for (int j=0; j<250; j++)
    //     {
    //         cout << stod(obs_content[i][j]) << " ";
    //     }
    //     cout << endl;
    // }
        
    if (cwipiParamsObs == 0){
	    for(int i=0; i<nb_oU; i++)
	    {
            obs_field(i)=stod(obs_content[i][timeIndex]);
            obs_field(i+nb_oU)=stod(obs_content[i+nb_oU][timeIndex]);

        }
    }
    else if (cwipiParamsObs == 1){
    	for(int i=0; i<nb_o; i++)
	    {
            obs_field(i)=stod(obs_content[i][timeIndex]);
        }
    }
    else if (cwipiParamsObs == 2){
        for(int i=0; i<nb_oU; i++)
        {
            obs_field(i)=stod(obs_content[i][timeIndex]);
            obs_field(i+nb_oU)=stod(obs_content[i + nb_oU][timeIndex]);

        }
        for(int i=0; i<(nb_o - 2*nb_oU); i++)
        {
            obs_field(i+2*nb_oU)=stod(obs_content[i + 2*nb_oU][timeIndex]);    
        } 
    }

    //std::cout << ""obs_field = " << obs_field  << std::endl << "\n" ;

    if (cwipiVerbose == 1) std::cout << "obs_field created" << std::endl << "\n";

    return obs_field;
}

//===========================================================================

MatrixXf samp_Data(int nb_e, int nb_o, int nb_oU, int cwipiParamsObs, int cwipiVerbose, std::string stringRootPath)
{
    //=== Read the sampled velocities from the files created by each OF instance. Return of matrix 
    // with all the sampled velocities organised in a X Y Z manner, each colunm corresponding to 
    // a member. 
    
    MatrixXf sampMatrix = MatrixXf::Zero(nb_o, nb_e); // By default our parameter is the velocity

    // std::cout << "sampMatrix size is " << sampMatrix.size() << "\n";

    // char UIntGen[500] = "/home/villanul/OpenFOAM/villanul-8/run/cwipi_tests/pitzDaily_cwipi/";
    std::string IntGen = stringRootPath + "/processor";
    if (cwipiVerbose == 1) std::cout << "Dans KF UIntGen = " << IntGen << std::endl;

    std::ifstream file;

    for (int i = 0; i < nb_e; ++i)
    {
        char ensemble[50];
        char member[50];
        if (cwipiParamsObs == 0){
            std::string UInt;
            sprintf(ensemble, "/UInt%i", i);
            sprintf(member, "%i", i);
            UInt += IntGen; // Append the first character
            UInt += member;
            UInt += ensemble; // Append the second character
            std::cout << "Opening " << UInt << std::endl;
            file.open(UInt);
        }
        else if (cwipiParamsObs == 1){
            std::string pInt;
            sprintf(ensemble, "/pInt%i", i);
            sprintf(member, "%i", i);
            pInt += IntGen; // Append the first character
            pInt += member;
            pInt += ensemble; // Append the second character
            file.open(pInt);
        }
        else if (cwipiParamsObs == 2){
            std::string UpInt;
            sprintf(ensemble, "/UpInt%i", i);
            sprintf(member, "%i", i);
            UpInt += IntGen; // Append the first character
            UpInt += member;
            UpInt += ensemble; // Append the second character
            // std::cout << "Opening " << UpInt << std::endl;
            file.open(UpInt);
        }

        if (cwipiVerbose == 1) std::cout << ensemble << std::endl;
        int count = 0;

        if (file.is_open())
        {
            std::string line;
            while( std::getline(file,line) )
            {
                std::stringstream ss(line);

                std::string IntValue;
                std::getline(ss, IntValue, ' ');
                if (count < nb_oU){
                    sampMatrix(count, i) = stod(IntValue);
                    // std::cout << "count inside if is: " << count << "\n";
                }

                if (count >= 2*nb_oU){
                    sampMatrix(count - nb_oU, i) = stod(IntValue);
                    // std::cout << "count inside if is: " << count << "\n";
                }

                // std::cout << "count is: " << count << "\n";
                ++count;
            }

        }
        else
        {
            std::cout << "Couldn't open the file" << std::endl << "\n";
        }

        file.close();
    }
    
    // std::cout << "We are here after file" << std::endl << "\n";
    if (cwipiVerbose == 1){
    for (int i=0; i< nb_o; i++)
    { 
        for (int j=0; j<nb_e; j++)
        {
            std::cout << sampMatrix(i, j) << " ";
        }
        std::cout << std::endl << "\n";
    }
    std::cout << std::endl << "\n";
    }

    return sampMatrix;
}

//===========================================================================

MatrixXf KF(MatrixXf velo_field_mat, ArrayXf obs_field, MatrixXf proj_field_mat, int nb_e, int nb_cells, int nb_o, int nb_oU, double sigma_userU, double sigma_userp, int nb_p, int cwipiParamsObs, int cwipiVerbose)
{
    //=== Perform the actual Ensemble Kalman Filter and return the updated state vector ===
    
    if (cwipiVerbose == 1) std::cout << "Opening KF" << std::endl;
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // --- Initialisation des variables ---
    // int nb_cells=2;
    int nb_cells_cmpnt=nb_cells*3;
    int nb_o_cmpnt = nb_o;
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

    if (cwipiVerbose == 1) std::cout << "Initialisation done" << std::endl;

    // Define random generator with Gaussian distribution
    const double mean_obs = 0;
    const double stddev_obsU = sigma_userU;
    const double stddev_obsp = sigma_userp;
    std::normal_distribution<double> dist_obsU(mean_obs, stddev_obsU);
    std::normal_distribution<double> dist_obsp(mean_obs, stddev_obsp);
    
    // const double stddev_obs = sigma_userU;
    // const double stddev_obs = sigma_userp;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    MatrixXf R(nb_o_cmpnt,nb_o_cmpnt);
    R.setIdentity();
    if (cwipiParamsObs == 0) R=R*stddev_obsU;
    else if (cwipiParamsObs == 1) R=R*stddev_obsp;
    else if (cwipiParamsObs == 2){
        R = R*stddev_obsp;
        R.block(0,0, 3*nb_oU, 3*nb_oU) = R.block(0,0, 3*nb_oU, 3*nb_oU)*stddev_obsU; 

        std::cout << "R :" << std::endl;
        std::cout << R << std::endl << "\n"; 
    }

    // Add Gaussian noise to create random velocity field and observation matrix
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
        }
        }
     else if (cwipiParamsObs == 2){
        for (int i=0;i<3*nb_oU;i++)
        {
            do
            {
                gaussample=dist_obsU(generator);
            }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
            err_mat(i,j)=gaussample;
        }
        for (int i=0;i<(nb_o - 3*nb_oU);i++)
        {
            do
            {
                gaussample=dist_obsp(generator);
            }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
            err_mat(3*nb_oU+i,j)=gaussample;
        }

        }


        obs_field_mat.col(j)=obs_field;
    }
    
    obs_field_mat=obs_field_mat+err_mat;

    if (cwipiVerbose == 1) std::cout << "Perturbations added" << std::endl;

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
    if (cwipiVerbose == 1) std::cout << "nb_cells_cmpnt = " << nb_cells_cmpnt << std::endl;
    if (cwipiVerbose == 1) std::cout << "nb_p = " << nb_p << std::endl;
    if (cwipiVerbose == 1) std::cout << "StateVector = " << velo_field_mat.size() << std::endl;

    for (int i=0;i<nb_cells_cmpnt+nb_p;i++)
    {
        velo_field_mean(i)=velo_field_mat.row(i).mean();
    }

    if (cwipiVerbose == 1) std::cout << "Means of the chosen parameter calculated" << std::endl;

    if (cwipiVerbose == 1) std::cout << "nb_o_cmpnt = " << nb_o_cmpnt << std::endl;

    for (int i=0;i<nb_o_cmpnt;i++)
    {
        err_mat_mean(i)=err_mat.row(i).mean();
        obs_field_mean(i)=obs_field_mat.row(i).mean();
        proj_field_mean(i)=proj_field_mat.row(i).mean();
    }

    if (cwipiVerbose == 1) std::cout << "Means calculated" << std::endl;
    
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

    if (cwipiVerbose == 1) std::cout << "Anomaly matrices calculated" << std::endl;

    // std::cout << "velo_field_mat_norm :" << std::endl;
    // std::cout << velo_field_mat_norm << std::endl << "\n";

    // std::cout << "proj_field_mat_norm :" << std::endl;
    // std::cout << proj_field_mat_norm << std::endl << "\n";

    // --- Kalman gain NAe_mat*trans(NAo_mat)*inv(NAo_mat*trans(NAo_mat))--- 
   
    // K21=proj_field_mat_norm*(proj_field_mat_norm.transpose())+obs_field_mat_norm*(obs_field_mat_norm.transpose());
    
    K21=proj_field_mat_norm*(proj_field_mat_norm.transpose())+R;
    
    std::cout << "Inverting the matrix for Kalman Gain" << std::endl;
    K2=K21.inverse();
    //K22=K21.determinant();

    K1=velo_field_mat_norm*(proj_field_mat_norm.transpose());

    K=K1*K2;

    if (cwipiVerbose == 1) std::cout << "Kalman gain calculated" << std::endl;
    
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
    if (cwipiVerbose == 1) std::cout << "Debut de la phase d'update" << std::endl;
    for (int j=0;j<nb_e;j++)
    {
        // std::cout << "======== Boucle numero " << j+1 << " ========" << std::endl << "\n";
        velo_field_mat_upt.col(j)=velo_field_mat.col(j)+K*(obs_field_mat.col(j)-proj_field_mat.col(j));

        // std::cout << "velo_field_mat :" << std::endl;
        // std::cout << velo_field_mat << std::endl << "\n";

        // std::cout << "velo_field_mat_upt :" << std::endl;
        // std::cout << velo_field_mat_upt << std::endl << "\n";
    }
    
    if (cwipiVerbose == 1) std::cout << "Filtre de Kalman termine" << std::endl << "\n";
    
    return velo_field_mat_upt;
}

//===========================================================================

void KF_output(double *sendValues, double *paramsSendValues, MatrixXf UptMatrix, int nb_e, int nb_cells, double time, int nb_p, int index, int cwipiVerbose)
{
    //=== Take the right column with the good values to send to the corresponding OF instance ===

    if (cwipiVerbose == 1) std::cout << "Opening output" << std::endl;

    for (int i=0; i<nb_cells; i++)
    {   
        sendValues[3*i + 0] = UptMatrix(i, index-1);
        sendValues[3*i + 1] = UptMatrix(i + nb_cells, index-1);    
        sendValues[3*i + 2] = UptMatrix(i + 2*nb_cells, index-1);    
    }

    for (int i = 0; i < nb_p; i++)
    {
        paramsSendValues[i] = UptMatrix(nb_cells*3+i, index-1);
    }

    //=== We do an output of all optimized coefficients to evaluate convergence ===//

    if (index == 1){
        std::fstream file_out;
        std::string filename = "UpdatedCoefficients";
        file_out.open(filename, std::ios_base::out | std::ios_base::app);
        file_out << time << ' ';
        
        for (int i = 0; i < nb_p; ++i){
            double sum = 0;
            for (int j = 0; j < nb_e; ++j){
                sum = sum + UptMatrix(nb_cells*3 + i, j);
            }
            double average = sum / nb_e;
            file_out << average << ' ';
        }
        file_out << "\n";
        file_out.close();
    }

}

//===========================================================================

void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, Eigen::MatrixXf stateVector, std::string name)
{
    //=== Write the state vector in a txt file ===
    
    if (phaseIndex >= numberCwipiPhase - cwipiOutputNb)
    {
    std::string stepNumber = std::to_string(time);
    std::string filename = name+stepNumber;
    std::fstream file_out;
    
    file_out.open(filename, std::ios_base::out);
    
    if (!file_out.is_open()) 
    {
        std::cout << "failed to open " << filename << '\n';
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
