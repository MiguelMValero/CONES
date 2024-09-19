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
 * @dir ./lib/CONES_interface/HLEnKF/src
 * @brief **Data Assimilation HLEnKF source files.**
 * 
 */

/**
 * @dir ./lib/CONES_interface/HLEnKF/src/DAFunctions
 * @brief **Data Assimilation HLEnKF functions declarations and defintions.**
 */

/**
 * @file HLEnKF_functions.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Data Assimilation function definitions.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "HLEnKF_functions.h"

// using namespace std;
using namespace Eigen;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//===========================================================================
void tic(int mode, std::string procedure) {
        static std::chrono::_V2::system_clock::time_point t_start;

        if (mode==0)
                t_start = std::chrono::high_resolution_clock::now();
        else {
                  auto t_end = std::chrono::high_resolution_clock::now();
                std::cout << "Elapsed time for the " << procedure << " is " << (t_end - t_start).count()*1E-9 << " seconds" << std::endl;
        }
}

void toc(std::string procedure) { tic(1, procedure); }

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
MatrixXd retrieve_regionState(const Eigen::Ref<const Eigen::MatrixXd>& invStateMatrix, region currentRegion, int nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose)
{
    /* Retrieve the state information of the current region */
    // if (cwipiVerbose) std::cout << "Entering the retrieve_regionState function" << std::endl;
    
    int nCellsRegion = currentRegion.get_nCells();
    ArrayXi cellIDs = currentRegion.get_cellIDs();
    MatrixXd regionState = MatrixXd::Zero(3 * nCellsRegion + cwipiParams, cwipiMembers);

    for (int i = 0; i < cwipiMembers; ++i){
        for (int j = 0; j < nCellsRegion; ++j){
            regionState(j, i) = invStateMatrix(cellIDs(j), i);
            regionState(nCellsRegion + j, i) = invStateMatrix(cellIDs(j) + nb_cells, i);
            regionState(2*nCellsRegion + j, i) = invStateMatrix(cellIDs(j) + 2*nb_cells, i);
        }

        for (int j = 0; j < cwipiParams; ++j){
            regionState(3*nCellsRegion + j, i) = invStateMatrix(3*nb_cells + j, i);
        }
    }

    return regionState;
}

//===========================================================================
void update_globalState(Eigen::Ref<Eigen::MatrixXd> stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXd>& regionState, region currentRegion, int nb_cells, int cwipiParams, int cwipiMembers, float cwipiVerbose)
{
    /* Update the global state  with information of the current region */
    // if (cwipiVerbose) std::cout << "Entering the update_globalState function" << std::endl;
    
    int nCellsRegion = currentRegion.get_nCells();
    ArrayXi cellIDs = currentRegion.get_cellIDs();

    for (int i = 0; i < cwipiMembers; ++i){
        for (int j = 0; j < nCellsRegion; ++j){
            stateMatrixUpt(cellIDs(j), i) = regionState(j, i);
            stateMatrixUpt(nb_cells + cellIDs(j), i) = regionState(j + nCellsRegion, i);
            stateMatrixUpt(2*nb_cells + cellIDs(j), i) = regionState(j + 2*nCellsRegion, i);
        }
    }
}

//===========================================================================
ArrayXd retrieve_obsVector(region currentRegion, observations obsData, double time, float cwipiVerbose)
{        
    /* ========= Read the observation data object return the values of the correct time in a eigen array ========== */
    // if (cwipiVerbose) std::cout << "Entering the retrieve_obsVector function" << std::endl;

    int timeObsIndex = 0;
    if (obsData.get_timeDependency())
    {
        timeObsIndex = (time-obsData.get_obsStartTime()) / obsData.get_obsTimeStep() -1;
    }

    ArrayXi region_probesIDs = currentRegion.get_probesIDs();
    int region_nProbes = currentRegion.get_nProbes();
    // if (cwipiVerbose) std::cout << "Number of region_nProbes (probes) in the region " << region_nProbes << std::endl;
    ArrayXd obsVector(currentRegion.get_nObs());
    // if (cwipiVerbose) std::cout << "Number of observations in the region " << currentRegion.get_nObs() << std::endl;
    std::string obsType = obsData.get_obsType();
    int velocityCaseSwitch = obsData.get_velocityCaseSwitch();
    int nObsProbesVel = currentRegion.get_nObsProbesVel();

    if (obsType == "U"){
		switch (velocityCaseSwitch){
			case 1 :
			case 2 :
			case 3 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
					std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
                    obsVector(i) = probeData[0][timeObsIndex];
				}
        		break;
			case 4 :
			case 5 :
			case 6 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
                    obsVector(i) = probeData[0][timeObsIndex];
                    obsVector(i + nObsProbesVel) = probeData[1][timeObsIndex];
				}
        		break;
			case 7 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					obsVector(i) = probeData[0][timeObsIndex];
                    obsVector(i + nObsProbesVel) = probeData[1][timeObsIndex];
                    obsVector(i + 2 * nObsProbesVel) = probeData[2][timeObsIndex];
				}
        		break;
		}
	}
	else if (obsType == "p"){
        for (int i = 0 ; i < region_nProbes; i++)
        {
            std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
            obsVector(i) = probeData[0][timeObsIndex];
        }
	}
	else if (obsType == "Up"){
        int velCounter = 0;
        int PresCounter = 0;
		switch (velocityCaseSwitch){
			case 1 :
			case 2 :
			case 3 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(PresCounter + nObsProbesVel) = probeData[0][timeObsIndex];
                        PresCounter = PresCounter + 1;
                    }
				}
				break;
			case 4 :
			case 5 :
			case 6 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        obsVector(velCounter + nObsProbesVel) = probeData[1][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(PresCounter + 2*nObsProbesVel) = probeData[0][timeObsIndex];
                        PresCounter = PresCounter + 1;
                    }
				}
        		break;
			case 7 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        obsVector(velCounter + nObsProbesVel) = probeData[1][timeObsIndex];
                        obsVector(velCounter + 2 * nObsProbesVel) = probeData[2][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(PresCounter + 3*nObsProbesVel) = probeData[0][timeObsIndex];
                        PresCounter = PresCounter + 1;
                    }
				}
        		break;
		}
	}
	else if (obsType == "UCf"){
		int velCounter = 0;
        int CfCounter = 0;
		switch (velocityCaseSwitch){
			case 1 :
			case 2 :
			case 3 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="Cf")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(CfCounter + nObsProbesVel) = probeData[0][timeObsIndex];
                        CfCounter = CfCounter + 1;
                    }
				}
				break;
			case 4 :
			case 5 :
			case 6 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="Cf")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        obsVector(velCounter + nObsProbesVel) = probeData[1][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(CfCounter + 2*nObsProbesVel) = probeData[0][timeObsIndex];
                        CfCounter = CfCounter + 1;
                    }
				}
        		break;
			case 7 :
      	        for (int i = 0 ; i < region_nProbes; i++)
				{
                    std::vector<std::vector<double>> probeData = obsData.get_probeData(region_probesIDs(i));
					if (obsData.get_probeVarLabels(region_probesIDs(i))!="Cf")
                    {
                        obsVector(velCounter) = probeData[0][timeObsIndex];
                        obsVector(velCounter + nObsProbesVel) = probeData[1][timeObsIndex];
                        obsVector(velCounter + 2 * nObsProbesVel) = probeData[2][timeObsIndex];
                        velCounter = velCounter + 1;
                    }
                    else
                    {
                        obsVector(CfCounter + 3*nObsProbesVel) = probeData[0][timeObsIndex];
                        CfCounter = CfCounter + 1;
                    }
				}
        		break;
		}
	}

    return obsVector;
}

//========================== Randomize observations =========================
MatrixXd randomize_Obs(const Eigen::Ref<const Eigen::ArrayXd>& obsVector, int cwipiMembers, const observations obsData, const region currentRegion, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    float gaussample;
    // int nProbesRegionVel = currentRegion.get_nObsProbesVel();
    int nProbesRegionPres = currentRegion.get_nObsProbesPres();
    int nProbesRegionCf = currentRegion.get_nObsProbesCf();
    int nObsRegion = currentRegion.get_nObs();
    std::string obsType = obsData.get_obsType();

    //========== Initialisation of the matrices ==========
    MatrixXd obsMatrix(nObsRegion,cwipiMembers);
    MatrixXd err_mat(nObsRegion,cwipiMembers);

    //========== Define random generator with Gaussian distribution ==========
    const double mean_obs = 0;
    const double stddev_obsU = obsData.get_std_u();
    const double stddev_obsp = obsData.get_std_p();
    const double stddev_obsCf = obsData.get_std_Cf();
    std::normal_distribution<double> dist_obsU(mean_obs, stddev_obsU);
    std::normal_distribution<double> dist_obsp(mean_obs, stddev_obsp);
    std::normal_distribution<double> dist_obsCf(mean_obs, stddev_obsCf);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    //========== Add Gaussian noise to create random velocity field and observation matrix ===========
    for (int j=0;j<cwipiMembers;j++)
    {
        if (obsType == "U"){
            for (int i=0;i<nObsRegion;i++)
            {
                do
                {
                    gaussample=dist_obsU(generator);
                }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                err_mat(i,j)=gaussample;
                obsMatrix(i,j)=obsVector(i);
            }
        }
        else if (obsType == "p"){
            for (int i=0;i<nObsRegion;i++)
            {
                do
                {
                    gaussample=dist_obsp(generator);
                }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                err_mat(i,j)=gaussample;
                obsMatrix(i,j)=obsVector(i);
            }
        }
        else if (obsType == "Up"){
            int nObsRegionU = nObsRegion - nProbesRegionPres;
            for (int i=0; i < (nObsRegionU); i++)
            {
                do
                {
                    gaussample = dist_obsU(generator);
                }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                err_mat(i,j) = gaussample;
                obsMatrix(i,j)=obsVector(i);
            }
            for (int i=0; i < (nProbesRegionPres); i++)
            {
                do
                {
                    gaussample=dist_obsp(generator);
                }while (gaussample<(mean_obs-3*stddev_obsp) || gaussample>(mean_obs+3*stddev_obsp));
                err_mat(nObsRegionU + i,j)=gaussample;
                obsMatrix(nObsRegionU + i,j)=obsVector(nObsRegionU + i);
            }
        }
        else if (obsType == "UCf"){
            int nObsRegionU = nObsRegion - nProbesRegionCf;
            for (int i=0; i < nObsRegionU; i++)
            {
                do
                {
                    gaussample = dist_obsU(generator);
                }while (gaussample<(mean_obs-3*stddev_obsU) || gaussample>(mean_obs+3*stddev_obsU));
                err_mat(i,j) = gaussample;
                obsMatrix(i,j)=obsVector(i);
            }
            for (int i=0; i < (nProbesRegionCf); i++)
            {
                do
                {
                    gaussample=dist_obsCf(generator);
                }while (gaussample<(mean_obs-3*stddev_obsCf) || gaussample>(mean_obs+3*stddev_obsCf));
                err_mat(nObsRegionU + i,j)=gaussample;
                obsMatrix(nObsRegionU + i,j)=obsVector(nObsRegionU + i);
            }
        }
    }
    if (obsData.get_obsNoiseType() == "absVal") obsMatrix = obsMatrix + err_mat;
    else if (obsData.get_obsNoiseType() == "relVal") obsMatrix = obsMatrix + obsMatrix.cwiseProduct(err_mat);
    // if (cwipiVerbose) std::cout << "obsMatrix = "<< obsMatrix << "\n" << std::endl;

    return obsMatrix;
}

//===========================================================================
MatrixXd read_SampData(int cwipiMembers, const observations obsData, float cwipiVerbose, std::string stringRootPath)
{
    //========== Read the sampled velocities (or else) from the files created by each OF instance. Return of matrix 
    // with all the sampled velocities organised in a X Y Z manner, each colunm corresponding to 
    // a member ==========     
    std::string obsType = obsData.get_obsType();
    std::vector<std::vector<double>> temp_sampData(cwipiMembers);
    std::string IntGen = stringRootPath + "/member";
    std::ifstream file;
    char* IntGenChar = new char[1000];
    strcpy(IntGenChar, IntGen.c_str());

    for (int i = 0; i < cwipiMembers; ++i)
    {
        char ensemble[50];
        char member[50];
        std::string sampFilePath;
        if (obsType == "U"){
            sprintf(IntGenChar, "%i", i+1);
            sprintf(ensemble, "/UInt%i", i+1);
            sprintf(member, "%i", i+1);
            sampFilePath += IntGen; // Append the path
            sampFilePath += member; // Append the member index
            sampFilePath += ensemble; // Append file name
        }
        else if (obsType == "p"){
            sprintf(ensemble, "/pInt%i", i+1);
            sprintf(member, "%i", i+1);
            sampFilePath += IntGen; // Append the path
            sampFilePath += member; // Append the member index
            sampFilePath += ensemble; // Append file name
        }
        else if (obsType == "Up"){
            sprintf(ensemble, "/UpInt%i", i+1);
            sprintf(member, "%i", i+1);
            sampFilePath += IntGen; // Append the path
            sampFilePath += member; // Append the member index
            sampFilePath += ensemble; // Append file name
        }
        else if (obsType == "UCf"){
            sprintf(ensemble, "/UCfInt%i", i);
            sprintf(member, "%i", i);
            sampFilePath += IntGen; // Append the path
            sampFilePath += member; // Append the member index
            sampFilePath += ensemble; // Append file name
        }
        // if (cwipiVerbose) std::cout << "Opening " << sampFilePath << std::endl;

        //=========== Extract sampled data from file ===========
        std::vector<std::vector<double>> samp_content;
        std::vector<double> samp_row;
        std::string samp_line, samp_word;
    
        std::fstream samp_stream (sampFilePath, std::ios::in);
        if(samp_stream.is_open())
        {
            while(getline(samp_stream, samp_line))
            {
                samp_row.clear();
    
                std::stringstream samp_str(samp_line);
    
                while(getline(samp_str, samp_word, ','))
                    samp_row.push_back(std::stod(samp_word));
                samp_content.push_back(samp_row);
            }
        }
        else {
            std::cerr << "Couldn't open observation file.\n";
        }
        
        // if (cwipiVerbose) std::cout << "Before affecting to Eigen array in read sampling samp_content.size() = " << samp_content.size() << std::endl;
        for (unsigned int j=0; j<samp_content.size(); j++)
        {
            temp_sampData[i].push_back(samp_content[j][0]);
        }
    }
    
    MatrixXd sampMatrix = MatrixXd::Zero(temp_sampData[0].size(), cwipiMembers);
    if (cwipiVerbose) std::cout << "Sampling matrix created" << std::endl;
    if (cwipiVerbose == 1){
    for (unsigned int i=0; i< temp_sampData[0].size(); i++)
    { 
        for (int j=0; j<cwipiMembers; j++)
        {
            sampMatrix(i, j) = temp_sampData[j][i];
            // std::cout << sampMatrix(i, j) << " ";
        }
        // std::cout << std::endl;
    }
    // std::cout << std::endl << "\n";
    }

    return sampMatrix;
}

//===========================================================================
MatrixXd retrieve_regionSampData(const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, int cwipiMembers, observations obsData, region currentRegion, float cwipiVerbose, std::string stringRootPath)
{
    /* ========= Read the samp matrix and return the values of the region probes only ========== */
    std::string obsType = obsData.get_obsType();
    int velocityCaseSwitch = obsData.get_velocityCaseSwitch();
    int totalObsProbesVel = obsData.get_nObsProbesVel();
    int regionObsProbesVel = currentRegion.get_nObsProbesVel();
    ArrayXi region_probesIDs = currentRegion.get_probesIDs();
    int region_nProbes = currentRegion.get_nProbes();
    MatrixXd regionSampMatrix = MatrixXd::Zero(currentRegion.get_nObs(), cwipiMembers);

    for (int j = 0; j < cwipiMembers; j++)
    {
        if (obsType == "U"){
            switch (velocityCaseSwitch){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        regionSampMatrix(i,j) = sampMatrix(region_probesIDs(i),j);
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        regionSampMatrix(i,j) = sampMatrix(region_probesIDs(i),j);
                        regionSampMatrix(i + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                    }
                    break;
                case 7 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        regionSampMatrix(i,j) = sampMatrix(region_probesIDs(i),j);
                        regionSampMatrix(i + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                        regionSampMatrix(i + 2*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + 2*totalObsProbesVel,j);
                    }
                    break;
            }
        }
        else if (obsType == "p"){
            for (int i = 0 ; i < region_nProbes; i++)
            {
                regionSampMatrix(i,j) = sampMatrix(region_probesIDs(i),j);
            }
        }
        else if (obsType == "Up"){
            int velCounter = 0;
            int PresCounter = 0;
            switch (velocityCaseSwitch){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i) ,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(PresCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            PresCounter = PresCounter + 1;
                        }
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i),j);
                            regionSampMatrix(velCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(PresCounter + 2*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            PresCounter = PresCounter + 1;
                        }
                    }
                    break;
                case 7 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i),j);
                            regionSampMatrix(velCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                            regionSampMatrix(velCounter + 2*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + 2*totalObsProbesVel,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(PresCounter + 3*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            PresCounter = PresCounter + 1;
                        }
                    }
                    break;
            }
        }
        else if (obsType == "UCf"){
            int velCounter = 0;
            int CfCounter = 0;
            switch (velocityCaseSwitch){
                case 1 :
                case 2 :
                case 3 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i) ,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(CfCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            CfCounter = CfCounter + 1;
                        }
                    }
                    break;
                case 4 :
                case 5 :
                case 6 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i),j);
                            regionSampMatrix(velCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(CfCounter + 2*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            CfCounter = CfCounter + 1;
                        }
                    }
                    break;
                case 7 :
                    for (int i = 0 ; i < region_nProbes; i++)
                    {
                        if (obsData.get_probeVarLabels(region_probesIDs(i))!="p")
                        {
                            regionSampMatrix(velCounter,j) = sampMatrix(region_probesIDs(i),j);
                            regionSampMatrix(velCounter + regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + totalObsProbesVel,j);
                            regionSampMatrix(velCounter + 2*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i) + 2*totalObsProbesVel,j);
                            velCounter = velCounter + 1;
                        }
                        else
                        {
                            regionSampMatrix(CfCounter + 3*regionObsProbesVel,j) = sampMatrix(region_probesIDs(i),j);
                            CfCounter = CfCounter + 1;
                        }
                    }
                    break;
            }
        }
    }

    return regionSampMatrix;
}

//============= Calculate measurments error covariance matrix ===============
MatrixXd calculate_R(const Eigen::Ref<const Eigen::ArrayXd>& obsVector, const region currentRegion, const observations obsData, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    MatrixXd R(currentRegion.get_nObs(),currentRegion.get_nObs());
    R.setIdentity();
    double sigmaUserU = obsData.get_std_u();
    double sigmaUserp = obsData.get_std_p();
    double sigmaUserCf = obsData.get_std_Cf();
    int cwipiObsU = currentRegion.get_nObsProbesVel();
    std::string obsType = obsData.get_obsType();

    // === Calculation of R ===
    if (obsType == "U"){ 
        if (obsData.get_obsNoiseType() == "absVal") R = R*sigmaUserU*sigmaUserU;
        else if (obsData.get_obsNoiseType() == "relVal") R.diagonal() = obsVector.cwiseProduct(obsVector)*sigmaUserU*sigmaUserU;
    }
    else if (obsType == "p"){
        if (obsData.get_obsNoiseType() == "absVal") R = R*sigmaUserp*sigmaUserp;
        else if (obsData.get_obsNoiseType() == "relVal") R.diagonal() = obsVector.cwiseProduct(obsVector)*sigmaUserp*sigmaUserp;
    }
    else if (obsType == "Up"){
        if (obsData.get_obsNoiseType() == "absVal") R = R*sigmaUserp*sigmaUserp;
        else if (obsData.get_obsNoiseType() == "relVal") R.diagonal() = obsVector.cwiseProduct(obsVector)*sigmaUserp*sigmaUserp;
        switch (obsData.get_velocityCaseSwitch()){
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
    else if (obsType == "UCf"){
        if (obsData.get_obsNoiseType() == "absVal") R = R*sigmaUserCf*sigmaUserCf;
        else if (obsData.get_obsNoiseType() == "relVal") R.diagonal() = obsVector.cwiseProduct(obsVector)*sigmaUserCf*sigmaUserCf;
        switch (obsData.get_velocityCaseSwitch()){
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

//========================== Calculate Kalman gain ==========================
MatrixXd calculate_K(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, 
                    const Eigen::Ref<const Eigen::MatrixXd>& obsMatrix, 
                    const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, 
                    const Eigen::Ref<const Eigen::MatrixXd>& R, 
                    int cwipiMembers, 
                    region currentRegion, 
                    int cwipiParams, 
                    float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    int nb_cells_cmpnt = currentRegion.get_nCells()*3;
    int cwipiObs = currentRegion.get_nObs();

    //========== Initialisation of the matrices ==========
    MatrixXd stateMatrix_upt(nb_cells_cmpnt + cwipiParams,cwipiMembers);
    ArrayXd stateMatrix_mean(nb_cells_cmpnt + cwipiParams);
    MatrixXd stateMatrix_norm(nb_cells_cmpnt + cwipiParams,cwipiMembers);

    ArrayXd obsMatrix_mean(cwipiObs);
    MatrixXd obsMatrix_norm(cwipiObs,cwipiMembers);

    ArrayXd sampMatrix_mean(cwipiObs);
    MatrixXd sampMatrix_norm(cwipiObs,cwipiMembers);

    MatrixXd K(nb_cells_cmpnt + cwipiParams, cwipiObs);
    MatrixXd K1;
    MatrixXd K2;
    MatrixXd K21;

    // ========== Ensemble means ===========
    for (int i=0; i < (nb_cells_cmpnt + cwipiParams); i++)
    {
        stateMatrix_mean(i)=stateMatrix.row(i).mean();
    }

    for (int i=0; i < cwipiObs; i++)
    {
        obsMatrix_mean(i)=obsMatrix.row(i).mean();
        sampMatrix_mean(i)=sampMatrix.row(i).mean();
    }
    // if (cwipiVerbose) std::cout << "All means calculated" << std::endl;


    // ========== Normalized anomalies (Ue_mat-Ume)/sqrt(Ne-1) ==========
    for (int i=0; i < (nb_cells_cmpnt + cwipiParams); i++)
    {
        for (int j=0; j < cwipiMembers; j++)
        {
            stateMatrix_norm(i,j)=(stateMatrix(i,j)-stateMatrix_mean(i))/std::sqrt(cwipiMembers-1);
        }
    }

    // if (cwipiVerbose) std::cout << "Anomaly state matrix calculated" << std::endl;

    // if (cwipiVerbose) std::cout << "sampMatrix = "<< sampMatrix << "\n" << std::endl;
    // if (cwipiVerbose) std::cout << "obsMatrix = "<< obsMatrix << "\n" << std::endl;

    for (int i=0; i < cwipiObs; i++)
    {
        // if (cwipiVerbose) std::cout << "i = " << i << std::endl;
        for (int j=0; j < cwipiMembers; j++)
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

//===========================================================================
void EnKF_outputs(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, const Eigen::Ref<const Eigen::MatrixXd>& sampMatrix, const observations obsData,int cwipiMembers, int nb_cells, double time, int cwipiParams, float cwipiVerbose, double epsilon)
{
    int cwipiObsU = obsData.get_nObsProbesVel();
    int cwipiObs = obsData.get_nObs();
    Foam::word obsType = obsData.get_obsType();
    int velocityCase = obsData.get_velocityCaseSwitch();
    ArrayXd obsArray(cwipiObs);

    int timeObsIndex = 0;
    if (obsData.get_timeDependency())
    {
        timeObsIndex = (time-obsData.get_obsStartTime()) / obsData.get_obsTimeStep() -1;
    }

    std::vector<std::vector<double>> obsRawData = obsData.get_obs_raw_data();
    for (int i = 0; i < cwipiObs; i++)
    {
        obsArray(i) = obsRawData[i][timeObsIndex];
    }

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
    if (obsType == "U"){
                switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
{
                ArrayXd MSvals1(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU); 
                for (int i = 0; i < cwipiObsU; ++i){
                    MSvals1(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
                    MSvobs1(i) = std::pow(obsArray(i), 2) + epsilon;
                }
                RMSD1 = std::sqrt(MSvals1.sum())/cwipiObsU; //Root Mean Square Deviation
                NRMSD1 = std::sqrt(MSvals1.sum()/MSvobs1.sum())/cwipiObsU; //Normalized Root Mean Square Deviation
                // NRMSD1 = std::sqrt(MSvals1.sum()/nb_oU)/obsArray.mean(); //Normalized Root Mean Square Deviation
                file_RMS_out << RMSD1 << ' ' << NRMSD1 << ' ';
                file_RMS_out << "\n";
}
                break;
            case 4 :
            case 5 :
            case 6 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU);
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
}
                break;
            case 7 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU);
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
}
                break;
        }
    }
    else if (obsType == "p"){
        ArrayXd MSvals4(cwipiObs);
        ArrayXd MSvobs4(cwipiObs);
        for (int i = 0; i < (cwipiObs-cwipiObsU); ++i){
            MSvals4(i) = std::pow(sampMatrix.row(i).mean() - obsArray(i), 2);
            MSvobs4(i) = std::pow(obsArray(i), 2) + epsilon;
        }
        RMSD4 = std::sqrt(MSvals4.sum())/cwipiObs; //Root Mean Square Deviation
        NRMSD4 = std::sqrt(MSvals4.sum()/MSvobs4.sum())/cwipiObs; //Normalized Root Mean Square Deviation
        file_RMS_out << RMSD4 << ' ' << NRMSD4 << ' ';
        file_RMS_out << "\n";
    }
    else if (obsType == "Up"){
                switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals4(cwipiObs-cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs4(cwipiObs-cwipiObsU);
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
}
                break;
            case 4 :
            case 5 :
            case 6 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals4(cwipiObs-2*cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs4(cwipiObs-2*cwipiObsU);
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
}
                break;
            case 7 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU), MSvals4(cwipiObs-3*cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU), MSvobs4(cwipiObs-3*cwipiObsU);
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
}
                break;
        }
    }
    else if (obsType == "UCf"){
                switch (velocityCase){
            case 1 :
            case 2 :
            case 3 :
{
                ArrayXd MSvals1(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU);
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
}
                break;
            case 4 :
            case 5 :
            case 6 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU);
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
}
                break;
            case 7 :
{
                ArrayXd MSvals1(cwipiObsU), MSvals2(cwipiObsU), MSvals3(cwipiObsU);
                ArrayXd MSvobs1(cwipiObsU), MSvobs2(cwipiObsU), MSvobs3(cwipiObsU);
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
}
                break;
        }
    file_RMS_out.close();
    }
}

//========================== Covariance localization ============================
MatrixXd localisation(const Eigen::Ref<const Eigen::MatrixXd>& KnoCorr, int cwipiParams, region currentRegion, double sigmaLocX, double sigmaLocY, double sigmaLocZ, observations obsData, const fvMesh& mesh, float cwipiVerbose)
{
    /* Weighting the gain depending on the distance between the observations and the cells */
    // if(cwipiVerbose) std::cout << "Entering localizaztion function" << std::endl;

    int cwipiObs = currentRegion.get_nObs();
    int nb_cells = currentRegion.get_nCells();
    std::string obsType = obsData.get_obsType();

    //========= Definition of the cells affected by localisation
    MatrixXd K(3*nb_cells+cwipiParams, cwipiObs);
    MatrixXd loc = MatrixXd::Zero(3*nb_cells+cwipiParams, cwipiObs);

    //========== We impose no localisation for the parameters to update ==========
    for (int i = 0; i < cwipiObs; ++i){
        for (int j = 0; j < cwipiParams; ++j){
            loc(3*nb_cells+j, i) = 1;
        }
    }

    double coordXdom, coordYdom, coordZdom;
  
    const Foam::volVectorField& C = mesh.C();

    // tic();
    ArrayXi cellIDs = currentRegion.get_cellIDs();
    MatrixXd obsCoords = currentRegion.get_obsCoords();
    // std::cout << "obsCoords = " << obsCoords << std::endl;
    for(int i = 0; i < nb_cells; ++i){
        // if(cwipiVerbose) std::cout << "cellID treated = " << cellID << std::endl;
        coordXdom = C[cellIDs(i)][0];
        coordYdom = C[cellIDs(i)][1];
        coordZdom = C[cellIDs(i)][2];
        for (int j = 0; j < cwipiObs; ++j){
            double deltaX = obsCoords(j,0) - coordXdom;
            double deltaY = obsCoords(j,1) - coordYdom;
            double deltaZ = obsCoords(j,2) - coordZdom;
            double deltaRNorm = std::pow(deltaX/sigmaLocX, 2) + std::pow(deltaY/sigmaLocY, 2) + std::pow(deltaZ/sigmaLocZ, 2);
            double lw = std::exp(-0.5*deltaRNorm);
            loc(i, j)              = lw;
            loc(i + nb_cells, j)   = lw;
            loc(i + 2*nb_cells, j) = lw;
        }
    }
    // toc("localization loop of the current region");
    
    // std::cout << "The localisation matrix is " << std::endl;
    // std::cout << loc << std::endl << "\n";
    //========== Multiplication of both matrices =========
    K = loc.cwiseProduct(KnoCorr);
    return K;
}

//======================== State Inflation function =========================
MatrixXd stateInflation(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, int cwipiMembers, const region currentRegion, int cwipiParams, double stateInfl, Foam::word typeInfl, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    int nb_cells_cmpnt=currentRegion.get_nCells()*3;
    float gaussample;

    MatrixXd stateMatrixInflated = stateMatrixUpt;

    // ========== Define random generator with Gaussian distribution ==========
    const double mean_state = 0;
    const double stddev_state = stateInfl;
    std::normal_distribution<double> dist_state(mean_state, stddev_state);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    if (stateInfl>0.000001){
    // if(cwipiVerbose) std::cout << "Applying " << typeInfl << " inflation on current region state" << std::endl << "\n"; 
        for (int j=0;j<cwipiMembers;j++)
        {

            for (int i = 0; i < nb_cells_cmpnt; i++)
            {
                do
                {
                    gaussample=dist_state(generator);
                }while (gaussample<(mean_state-2*stddev_state) || gaussample>(mean_state+2*stddev_state));

                if (typeInfl == "deterministic"){
                    stateMatrixInflated(i,j) = stateMatrixUpt.row(i).mean() + (1+stddev_state+(0.1*gaussample))*(stateMatrixUpt(i,j) - stateMatrixUpt.row(i).mean());  //Deterministic
                }
                else stateMatrixInflated(i,j) = stateMatrixUpt(i,j) + gaussample*stateMatrixUpt(i,j);     //Stochastic
            }
        }    
    }
   
    return stateMatrixInflated;
}

//====================== Parameters Inflation function ======================
MatrixXd paramInflation(const Eigen::Ref<const Eigen::MatrixXd>& parameters, int cwipiMembers, int cwipiParams, double paramsInfl, Foam::word typeInfl, float cwipiVerbose)
{
    //========== Initialisation of the variables ==========
    float gaussample;
    MatrixXd parametersInflated = parameters;

    // ========== Define random generator with Gaussian distribution ==========
    const double mean_params = 0;
    const double stddev_params = paramsInfl;
    std::normal_distribution<double> dist_params(mean_params, stddev_params);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    if (paramsInfl>0.000001){
    if(cwipiVerbose) std::cout << "Applying " << typeInfl << " inflation on parameters" << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {
            for (int i = 0; i < cwipiParams; i++)
            {
                do
                {
                    gaussample=dist_params(generator);
                }while (gaussample<(mean_params-2*stddev_params) || gaussample>(mean_params+2*stddev_params));
                
                if (typeInfl == "deterministic"){
                    // if(cwipiVerbose) std::cout << "typeInfl is in if condition " << typeInfl << std::endl; 
                    parametersInflated(i, j) = parameters.row(i).mean() + (1+stddev_params+(0.1*gaussample))*(parameters(i, j) - parameters.row(i).mean()); //Deterministic
                }
                else parametersInflated(i, j) = parameters(i, j) + gaussample*parameters(i, j);  //Stochastic
            }
        }
    }
   
    return parametersInflated;
}

//=========================== HLEnKF algorithm ==============================
MatrixXd HLEnKF(const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, const fvMesh& mesh, int cwipiMembers, int nb_cells, std::vector<region> regionList, observations obsData,double sigmaLocX, double sigmaLocY, double sigmaLocZ, bool localSwitch, int cwipiParams, double stateInfl, double preStateInfl, double paramsInfl, Foam::word typeInfl, bool paramEstSwitch, bool stateEstSwitch, float cwipiVerbose, std::string stringRootPath, double time, double epsilon)
{
    /* ================= Begin calculation of the EnKF =================== */ 
    if(cwipiVerbose) std::cout << "** ENTERING MAIN EnKF FUNCTION **" << std::endl;
    
    MatrixXd stateMatrixUpt = stateMatrix;
    MatrixXd parameters = MatrixXd::Zero(cwipiParams, cwipiMembers);

    //** Create the global sampMatrix containing all the forecasted observations **
    MatrixXd sampMatrix = read_SampData(cwipiMembers, obsData, cwipiVerbose, stringRootPath);

    tic();
    if (cwipiVerbose) std::cout <<  "Beginning regions loop" << std::endl;
    for (unsigned int i = 0; i < regionList.size(); i++)
    {
        // if (cwipiVerbose) std::cout <<  "Treatment of region " << i << std::endl;
        //** retrieve the state matrix of the region **
        // tic();
        MatrixXd regionState = retrieve_regionState(stateMatrix, regionList[i], nb_cells, cwipiParams, cwipiMembers, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "regionState created " << std::endl;
        // toc("regionState");
        // if (cwipiVerbose) std::cout << "\t regionState = \n " << regionState << std::endl << "\n";

        //** Perform a pre-inflation on temp_sate withtout parameter inflation (considered as model error) **
        regionState = stateInflation(regionState, cwipiMembers, regionList[i], cwipiParams, preStateInfl, typeInfl, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "regionState pre inflated " << std::endl;
        // if (cwipiVerbose) std::cout << "\t regionState = \n " << regionState << std::endl << "\n";
        
        //** retrieve observations of the current region **
        // tic();
        ArrayXd obsVector = retrieve_obsVector(regionList[i], obsData, time, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "obsVector created " << std::endl;
        // toc("obsVector");
        // if (cwipiVerbose) std::cout << "\t obsVector = \n " << obsVector << std::endl << "\n";

        //** Extend observations to all the members via perturbations addition **
        // tic();
        MatrixXd regionObsMatrix = randomize_Obs(obsVector, cwipiMembers, obsData, regionList[i], cwipiVerbose);
        // if(cwipiVerbose) std::cout << "regionObsMatrix created " << std::endl;
        // toc("regionObsMatrix");
        // if (cwipiVerbose) std::cout << "\t regionObsMatrix filled = \n " << regionObsMatrix << std::endl << "\n";

        //** Retrieve forecasted obersvations of the current region **
        // tic();
        MatrixXd regionSampMatrix = retrieve_regionSampData(sampMatrix, cwipiMembers, obsData, regionList[i], cwipiVerbose, stringRootPath);
        // if(cwipiVerbose) std::cout << "regionSampMatrix created " << std::endl;
        // toc("regionSampMatrix");
        // if (cwipiVerbose) std::cout << "\t regionSampMatrix filled = \n " << regionSampMatrix << std::endl << "\n";

        //** Calculate R **
        // tic();
        MatrixXd R = calculate_R(obsVector, regionList[i], obsData, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "R created " << std::endl;
        // toc("R");
        // if (cwipiVerbose) std::cout << "\t R filled = \n " << R << std::endl << "\n";

        //** Calculate the Kalman gain for temp_state **
        // tic();
        MatrixXd K = calculate_K(regionState, regionObsMatrix, regionSampMatrix, R, cwipiMembers, regionList[i], cwipiParams, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "K created " << std::endl;
        // toc("K");
        // if (cwipiVerbose) std::cout << "\t K filled = \n " << K << std::endl << "\n";

        //** Perform covariance localization for temporary Kalman gain **
        // tic();
        if(localSwitch) K = localisation(K, cwipiParams, regionList[i], sigmaLocX, sigmaLocY, sigmaLocZ, obsData, mesh, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "K localized " << std::endl;
        // toc("localization of K");
        // if (cwipiVerbose) std::cout << "\t K localizezd = \n " << K << std::endl << "\n";

        //** Update temporary state
        // tic();
        // MatrixXd correction = K*(obsMatrix - sampMatrix);
        // regionState = regionState + correction
        regionState = regionState + K*(regionObsMatrix - regionSampMatrix);
        // if(cwipiVerbose) std::cout << "regionState updated " << std::endl;
        // toc("regionState");
        // if (cwipiVerbose) std::cout << "\t regionState updated  = \n " << regionState << std::endl << "\n";

        //** Perform inflation on temp_sate withtout parameter inflation **
        // tic();
        regionState = stateInflation(regionState, cwipiMembers, regionList[i], cwipiParams, stateInfl, typeInfl, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "regionState inflated " << std::endl;
        // toc("regionState");
        // if (cwipiVerbose) std::cout << "\t regionState inflated  = \n " << regionState << std::endl << "\n";

        //** Update the global state matrix **
        // tic();
        update_globalState(stateMatrixUpt, regionState, regionList[i], nb_cells, cwipiParams, cwipiMembers, cwipiVerbose);
        // if(cwipiVerbose) std::cout << "stateMatrixUpt updated with updated region " << std::endl;
        // toc("stateMatrixUpt");

        //** Add up updated params to perform an average on the values updated by each region **
        // tic();
        int region_nCells = regionList[i].get_nCells();
        int region_nProbes = regionList[i].get_nProbes();
        for (int j = 0; j < cwipiMembers; ++j){
            for (int k = 0; k < cwipiParams; ++k){
                parameters(k, j) = parameters(k, j) + regionState(3*region_nCells + k, j)*region_nProbes;
            }
        }
        // toc("update params");
        
    }
    if (cwipiVerbose) std::cout <<  "Regions loop done" << std::endl;
    toc("regions loop");

    //** If parameter estimation is activated average, inflate and affect to global state **
    if (paramEstSwitch)
    {
        parameters = parameters/obsData.get_nObsProbes();
        // if (cwipiVerbose) std::cout << "\t Averaged parameters  = \n " << parameters << std::endl << "\n";
        parameters = paramInflation(parameters, cwipiMembers, cwipiParams, paramsInfl, typeInfl, cwipiVerbose);
        // if (cwipiVerbose) std::cout << "\t Averaged parameters inflated  = \n " << parameters << std::endl << "\n";
        for (int j = 0; j < cwipiMembers; j++){
            for (int k = 0; k < cwipiParams; k++){
                stateMatrixUpt(3*nb_cells + k, j) = parameters(k, j);
            }
        }
    }

    //** If State estimation desactivated return to original values **
    if(!stateEstSwitch){
        if (cwipiVerbose) std::cout << "State estimation desactivated, returning to original values ..." << std::endl;
        for (int j=0;j<cwipiMembers;j++)
        {
            for (int i = 0; i < nb_cells*3; i++)
            {
                stateMatrixUpt(i,j) = stateMatrix(i,j);
            }
        }
    }


    if (cwipiVerbose) std::cout << "Starting writing outputs from the EnKF..." << std::endl;
    EnKF_outputs(stateMatrixUpt, sampMatrix, obsData, cwipiMembers, nb_cells, time, cwipiParams, cwipiVerbose, epsilon);
    if (cwipiVerbose) std::cout << "** End of Kalman Filter **" << std::endl << "\n";

    return stateMatrixUpt;
}

//===========================================================================
void prepare_sendBack(double *sendValues, double *paramsSendValues, double *values, const Eigen::Ref<const Eigen::MatrixXd>& stateMatrixUpt, int nb_cells, int cwipiParams, int index, float cwipiVerbose)
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
void print_matrix_on_txt(int phaseIndex, int numberCwipiPhase, int cwipiOutputNb, int cwipiMembers, int nb_cells, int cwipiParams, double time, const Eigen::Ref<const Eigen::MatrixXd>& stateMatrix, std::string name)
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