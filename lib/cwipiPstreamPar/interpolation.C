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
 * @file interpolation.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Defines interpolations.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "volPointInterpolation.H"
#include "interpolationCellPointWallModified.H"
#include "cwipiPstreamPar.H"
#include <cwipi.h>

namespace Foam
{

/**
 * @brief 
 * 
 * @param U 
 * @param mesh 
 * @param runTime 
 * @param cwipiObsU 
 * @param nbParts 
 * @param cwipiVerbose 
 * @param globalPath 
 * @param UIntPath 
 * @param filename 
 */
void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int nbParts, float cwipiVerbose, std::string globalPath, std::string UIntPath, const std::string filename)
{
    //========== Produce a file for each OF instance containing the sampled velocities 
    //(H.x term in the kalman gain calculation) ==========
    interpolationCellPoint<vector> triangulateCellsU(U); 
    
    //========== Path of the correct OF instance ==========
    std::string obs_file_path = globalPath + "/" + filename;
    if (cwipiVerbose) if (Pstream::master()) Pout << "The path where I have my observation coordinates is " << obs_file_path << endl;
    char UIntGen[500];
    strcpy(UIntGen, UIntPath.c_str());

    int myGlobalRank;
    char procChar[50];
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    label nMyProc=Pstream::myProcNo();
    std::string procString = std::to_string(nMyProc);
    Info << "nMyProc = " << nMyProc << endl; 

    sprintf(procChar, "%i", nMyProc);
    // strcat(UIntGen, "/processor");
    // strcat(UIntGen, procChar);
    strcat(UIntGen, "/UInt");

    label n=Pstream::nProcs();
    int ens = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;
    
    std::string ensString = std::to_string(ens);
    // std::string sampString = ensString + procString;
    // char sampChar[500];
    // strcpy(sampChar, sampString.c_str());
    char ensChar[500];
    strcpy(ensChar, ensString.c_str());

    // char* UIntPart_path = strcat(UIntGen, sampChar);
    char* UIntEns_path = strcat(UIntGen, ensChar);

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The file where UInt is created is " << UIntEns_path << endl;

    //* Configuration of streams to open the files : observation coordinates, sampling partioned file, member sampling file *
    std::ifstream obs_file;
    obs_file.open(obs_file_path);

    // std::ofstream UintPart_file;
    // UintPart_file.open(UIntPart_path, std::ios::out | std::ios::app);

    // string data = "";
    // int columns = 0;
    // int parameters = 3;

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "Creating the sampling files..." << endl;

    point pointCoord;
    double Ux[cwipiObsU] = {0}; 
    double Uy[cwipiObsU] = {0};
    double Uz[cwipiObsU] = {0};
    // int cellID[cwipiObsU];
    int cellIndex[cwipiObsU] = {0};

    // remove(UIntPart_path);

    if (obs_file.is_open())
    {
        // if (cwipiVerbose) std::cout<< "Inside the observation coordinates file" << std::endl;
        std::string line;
        int countProbes = 0;
        while( std::getline(obs_file,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
            pointCoord.x() = std::stod(coordX);
            std::getline(ss, coordY,',');
            pointCoord.y() = std::stod(coordY);
            std::getline(ss, coordZ,','); 
            pointCoord.z() = std::stod(coordZ);

            label ownCell = mesh.findCell(pointCoord);
            if (ownCell != -1){
                vector UatSP(vector::zero);
                UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

                Ux[countProbes] = UatSP.x();
                //std::cout << "The value of Ux for observation " << countProbes << " is " << Ux[countProbes] << std::endl;
                Uy[countProbes] = UatSP.y();
                //std::cout << "The value of Uy for observation " << countProbes << " is " << Uy[countProbes] << std::endl;
                Uz[countProbes] = UatSP.z();
                //std::cout << "The value of Uz for observation " << countProbes << " is " << Uz[countProbes] << std::endl;
                // cellID[countProbes] = ownCell;
                // //std::cout << "The value of the cellID for observation " << countProbes << " is " << cellID[countProbes] << std::endl;
                cellIndex[countProbes] = countProbes+1; //index in the list of observations, we do not use 0 because it corresponds to no cell in the domain

                // ++columns;
                // if (columns == parameters){
                //     columns = 0;
                // }
                // Pout<< "probe number " << countProbes + 1 << " found here" << endl;
            }
            ++countProbes;
        }
    }
    else{
        Foam::Pout << " Not able to open the obs_coordinates.txt file " << endl;
        // printing the error message
        fprintf(stderr, "File does not exist\n");
    }
    obs_file.close();


    /* Gather all the values from the procs and then write the new values in a txt file */
    if (Pstream::master()){
        /* Gather the values */
        int tempCellIndex = 0;
        double tempUx = 0; 
        double tempUy = 0;
        double tempUz = 0;

        for(int i=1; i<n; i++)
        {
            // Pout << "Receiving sampled data from proc " << i << endl;  
                // create the input stream from processor i
            IPstream streamFromMain(Pstream::commsTypes::blocking, i);
            for (int i = 0; i < cwipiObsU; ++i){
                streamFromMain >> tempCellIndex;
                streamFromMain >> tempUx;
                streamFromMain >> tempUy;
                streamFromMain >> tempUz;
                cellIndex[i] = cellIndex[i] + tempCellIndex;
                Ux[i] = Ux[i] + tempUx;
                Uy[i] = Uy[i] + tempUy;
                Uz[i] = Uz[i] + tempUz;
            }
        }

        /* Write them */
        remove(UIntEns_path);
        std::ofstream UintEns_file;
        UintEns_file.open(UIntEns_path, std::ios::out | std::ios::app);
        if (cwipiVerbose) Foam::Pout << "I am about to write the sampling file" << endl;
        if (UintEns_file.is_open()){
            // for (int i = 0; i < cwipiObsU; ++i){
            //     UintEns_file << cellIndex[i] << "\n";
            // }
            for (int i = 0; i < cwipiObsU; ++i){
                
                UintEns_file << Ux[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                
                UintEns_file << Uy[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                
                UintEns_file << Uz[i] << "\n";
            }
        }
        else{
            std::cerr << "Not able to open the sampling file" << std::endl;
        }
        UintEns_file.close();
    }
    else{
        // create the stream to send to the main proc
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        for (int i = 0; i < cwipiObsU; ++i){
            stream2Main << cellIndex[i];
            stream2Main << Ux[i];
            stream2Main << Uy[i];
            stream2Main << Uz[i];
        }
    }

    //========== Retrieve the cellIDs of my probes ==========//
    // if (myGlobalRank == 1 && (runTime.value() - runTime.deltaTValue()) == 0){
    //     std::ofstream myfile_2;
    //     std::string cellSensors = globalPath + "/results/cellIDs";

    //     remove(cellSensors.c_str());
    //     myfile_2.open(cellSensors, std::ios::out | std::ios::app);
    //     for (int i = 0; i < cwipiObsU; ++i){
    //         myfile_2 << cellID[i] << "\n";
    //     }
    //     myfile_2.close();
    // }

    if (cwipiVerbose) if (Pstream::master()) Pout << "Sampling files created " << endl;
}

/**
 * @brief 
 * 
 * @param p 
 * @param mesh 
 * @param runTime 
 * @param cwipiObsp 
 * @param nbParts 
 * @param cwipiVerbose 
 * @param globalPath 
 * @param pIntPath 
 * @param filename 
 */
void pInterpolation(volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string pIntPath, const std::string filename)
{
    //========== Produce a file for each OF instance containing the sampled pressures 
    //(H.x term in the kalman gain calculation) ========== 
    interpolationCellPoint<double> triangulateCellsp(p);

    //========== Path of the correct OF instance ==========
    std::string obs_file_path = globalPath + "/" + filename;
    if (cwipiVerbose) if (Pstream::master()) Pout << "The path where I have my observation coordinates is " << obs_file_path << endl;
    char pIntGen[500];
    strcpy(pIntGen, pIntPath.c_str());

    int myGlobalRank;
    char procChar[50];
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    label nMyProc=Pstream::myProcNo();
    std::string procString = std::to_string(nMyProc);
    Info << "nMyProc = " << nMyProc << endl; 

    sprintf(procChar, "%i", nMyProc);
    // strcat(pIntGen, "/processor");
    // strcat(pIntGen, procChar);
    strcat(pIntGen, "/pInt");

    label n=Pstream::nProcs();
    int ens = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;
    
    std::string ensString = std::to_string(ens);
    // std::string sampString = ensString + procString;
    // char sampChar[500];
    // strcpy(sampChar, sampString.c_str());
    char ensChar[500];
    strcpy(ensChar, ensString.c_str());

    // char* pIntPart_path = strcat(pIntGen, sampChar);
    char* pIntEns_path = strcat(pIntGen, ensChar);

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The file where pInt is created is " << pIntEns_path << endl;

    //* Configuration of streams to open the files : observation coordinates, sampling partioned file, member sampling file *
    std::ifstream obs_file;
    obs_file.open(obs_file_path);

    // std::ofstream UintPart_file;
    // UintPart_file.open(UIntPart_path, std::ios::out | std::ios::app);

    // string data = "";
    // int columns = 0;
    // int parameters = 1;

     if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "Creating the sampling files..." << endl;

    point pointCoord;
    int cellIndex[cwipiObsp] = {0};
    double pp[cwipiObsp] = {0};
    
    if (obs_file.is_open())
    {  
        std::string line;
        int countProbes = 0;
        while( std::getline(obs_file,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
            pointCoord.x() = stod(coordX);
            std::getline(ss, coordY,',');
            pointCoord.y() = stod(coordY);
            std::getline(ss, coordZ,','); 
            pointCoord.z() = stod(coordZ);

            label ownCell = mesh.findCell(pointCoord);
            if (ownCell != -1){
                pp[countProbes] = triangulateCellsp.interpolate(pointCoord, ownCell, -1);
                cellIndex[countProbes] = countProbes+1; //index in the list of observations, we do not use 0 because it corresponds to no cell in the domain
            }
            ++countProbes;
        }
    }
    else{
        Foam::Pout << " Not able to open the obs_coordinates.txt file " << endl;
        // printing the error message
        fprintf(stderr, "File does not exist\n");
    }
    obs_file.close();
                
    /* Gather all the values from the procs and then write the new values in a txt file */
    if (Pstream::master()){
        /* Gather the values */
        int tempCellIndex = 0;
        double tempp = 0; 

        for(int i=1; i<n; i++)
        {
            // Pout << "Receiving sampled data from proc " << i << endl;  
                // create the input stream from processor i
            IPstream streamFromMain(Pstream::commsTypes::blocking, i);
            for (int i = 0; i < cwipiObsp; ++i){
                streamFromMain >> tempCellIndex;
                streamFromMain >> tempp;
                cellIndex[i] = cellIndex[i] + tempCellIndex;
                p[i] = p[i] + tempp;
            }
        }

        /* Write them */
        remove(pIntEns_path);
        std::ofstream pintEns_file;
        pintEns_file.open(pIntEns_path, std::ios::out | std::ios::app);
        if (cwipiVerbose) Foam::Pout << "I am about to write the sampling file" << endl;
        if (pintEns_file.is_open()){
            // for (int i = 0; i < cwipiObsU; ++i){
            //     UintEns_file << cellIndex[i] << "\n";
            // }
            for (int i = 0; i < cwipiObsp; ++i){
                
                pintEns_file << p[i] << "\n";
            }
        }
        else{
            std::cerr << "Not able to open the sampling file" << std::endl;
        }
        pintEns_file.close();
    }
    else{
        // create the stream to send to the main proc
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        for (int i = 0; i < cwipiObsp; ++i){
            stream2Main << cellIndex[i];
            stream2Main << p[i];
        }
    }

        //========== Retrieve the cellIDs of my probes ==========//
    // if (myGlobalRank == 1){
    //     std::ofstream myfile_2;
    //     std::string cellSensors = pIntPath + "/results/cellIDs";

    //     if (remove(cellSensors.c_str()) != 0)
    //         perror( "Error deleting file with cell IDs of the sensors" );
    //     myfile_2.open(cellSensors, std::ios::out | std::ios::app);
    //         for (int i = 0; i < cwipiObsp; ++i){
    //             myfile_2 << cellID[i] << "\n";
    //         }
    // }
    if (cwipiVerbose) if (Pstream::master()) Pout << "Sampling files created " << endl;
}

/**
 * @brief 
 * 
 * @param U 
 * @param p 
 * @param mesh 
 * @param runTime 
 * @param cwipiObsU 
 * @param cwipiObsp 
 * @param nbParts 
 * @param cwipiVerbose 
 * @param globalPath 
 * @param UpIntPath 
 * @param filename 
 */
void UpInterpolation(volVectorField& U, volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsU, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string UpIntPath, const std::string filename)
{
    //=== Produce a file for each OF instance containing the sampled velocities and pressures 
    //(H.x term in the kalman gain calculation) === 
    interpolationCellPoint<vector> triangulateCellsU(U);
    interpolationCellPoint<double> triangulateCellsp(p);

    //* Path of the correct OF instance *
    std::string obs_file_path = globalPath + "/" + filename;
    if (cwipiVerbose) if (Pstream::master()) Pout << "The path where I have my observation coordinates is " << obs_file_path << endl;
    char UpIntGen[500];
    strcpy(UpIntGen, UpIntPath.c_str());

    int myGlobalRank;
    char procChar[50];
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

    label nMyProc=Pstream::myProcNo();
    std::string procString = std::to_string(nMyProc);
    Info << "nMyProc = " << nMyProc << endl; 

    sprintf(procChar, "%i", nMyProc);
    // strcat(UpIntGen, rankChar);
    // strcat(UpIntGen, "/UpInt");
    strcat(UpIntGen, "/UpInt");

    label n=Pstream::nProcs();
    int ens = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;

    std::string ensString = std::to_string(ens);
    // std::string sampString = ensString + procString;
    // char sampChar[500];
    // strcpy(sampChar, sampString.c_str());
    char ensChar[500];
    strcpy(ensChar, ensString.c_str());

    // char* UIntPart_path = strcat(UIntGen, sampChar);
    char* UpIntEns_path = strcat(UpIntGen, ensChar);

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "The file where UpInt is created is " << UpIntEns_path << endl;

    //* Configuration of streams to open the files : observation coordinates, sampling partioned file, member sampling file *
    std::ifstream obs_file;
    obs_file.open(obs_file_path);

    // std::ofstream UpintPart_file;
    // UintPart_file.open(UpIntPart_path, std::ios::out | std::ios::app);
    
    // string data = "";
    // int columns = 0;
    // int parameters = 3; // Always coordinates here 

    if (cwipiVerbose) if (Pstream::master()) Foam::Pout<< "Creating the sampling files..." << endl;

    point pointCoord;
    double Ux[cwipiObsU] = {0};
    double Uy[cwipiObsU] = {0};
    double Uz[cwipiObsU] = {0};
    double pInt[cwipiObsp] = {0};
    int cellIndex[cwipiObsU + cwipiObsp] = {0};

    // remove(UIntPart_path);

    if (obs_file.is_open())
    {
        std::string line;
        int countProbes = 0;
        while( std::getline(obs_file,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
            pointCoord.x() = stod(coordX);
            std::getline(ss, coordY,',');
            pointCoord.y() = stod(coordY);
            std::getline(ss, coordZ,','); 
            pointCoord.z() = stod(coordZ);

            label ownCell = mesh.findCell(pointCoord);
            if (ownCell != -1){
                if (countProbes < cwipiObsU) {
                    vector UatSP(vector::zero);
                    UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

                    Ux[countProbes] = UatSP.x();
                    Uy[countProbes] = UatSP.y();
                    Uz[countProbes] = UatSP.z();
                }
                else {
                    scalar patSP = triangulateCellsp.interpolate(pointCoord, ownCell, -1);
                    pInt[countProbes - cwipiObsU] = patSP;
                }
            }
            ++countProbes;
        }
    }
    else {
        Foam::Pout << " Not able to open the obs_coordinates.txt file " << endl;
        // printing the error message
        fprintf(stderr, "File does not exist\n");
    }
    obs_file.close();

    /* Gather all the values from the procs and then write the new values in a txt file */
    if (Pstream::master()){
        /* Gather the values */
        int tempCellIndex = 0;
        double tempUx = 0; 
        double tempUy = 0;
        double tempUz = 0;
        double tempp = 0;
 
        for(int i = 1; i < n; i++)
        {
            // Pout << "Receiving sampled data from proc " << i << endl;  
                // create the input stream from processor i
            IPstream streamFromMain(Pstream::commsTypes::blocking, i);
            for (int j = 0; j < cwipiObsU; ++j){
                streamFromMain >> tempCellIndex;
                streamFromMain >> tempUx;
                streamFromMain >> tempUy;
                streamFromMain >> tempUz;
                cellIndex[j] = cellIndex[j] + tempCellIndex;
                Ux[j] = Ux[j] + tempUx;
                Uy[j] = Uy[j] + tempUy;
                Uz[j] = Uz[j] + tempUz;
            }
            for (int j = cwipiObsU; j < (cwipiObsU + cwipiObsp); ++j){
                streamFromMain >> tempCellIndex;
                streamFromMain >> tempp;
                cellIndex[j] = cellIndex[j] + tempCellIndex;
                pInt[j - cwipiObsU] = pInt[j - cwipiObsU] + tempp;
            }
        }

        /* Write them */
        remove(UpIntEns_path);
        std::ofstream UpintEns_file;
        UpintEns_file.open(UpIntEns_path, std::ios::out | std::ios::app);
        if (UpintEns_file.is_open()){
            for (int i = 0; i < cwipiObsU; ++i){
                UpintEns_file << Ux[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                UpintEns_file << Uy[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                UpintEns_file << Uz[i] << "\n";
            }
            for (int i = 0; i < cwipiObsp; ++i){
                UpintEns_file << p[i] << "\n";
            }
        }
        else {
            std::cerr << "Not able to open the sampling file" << std::endl;
        }
        UpintEns_file.close();
    }
    else {       
    // create the stream to send to the main proc
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        for (int i = 0; i < cwipiObsU; ++i){
            stream2Main << cellIndex[i];
            stream2Main << Ux[i];
            stream2Main << Uy[i];
            stream2Main << Uz[i];
        }
        for (int i = cwipiObsU; i < (cwipiObsU + cwipiObsp); ++i){
            stream2Main << cellIndex[i];
            stream2Main << pInt[i - cwipiObsU];
        }
    }

    //========== Retrieve the cellIDs of my probes ==========//
    // if (myGlobalRank == 1){
    //     std::ofstream myfile_2;
    //     std::string cellSensors = UpIntPath + "/results/cellIDs";

    //     if (remove(cellSensors.c_str()) != 0)
    //         perror( "Error deleting file with cell IDs of the sensors" );
    //     myfile_2.open(cellSensors, std::ios::out | std::ios::app);
    //         for (int i = 0; i < (cwipiObsU + cwipiObsp); ++i){
    //             myfile_2 << cellID[i] << "\n";
    //         }
    // }

    if (cwipiVerbose) if (Pstream::master()) Pout << "Sampling files created " << endl;
}




void UCfInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, scalar dragCoeff, scalar liftCoeff, int nbParts, 
int cwipiStep, int cwipiVerbose, std::string globalPath, std::string UCfIntPath)
{
    //========== Produce a file for each OF instance containing the sampled velocities/global drag coefficient 
    //(H.x term in the kalman gain calculation) ==========
    interpolationCellPoint<vector> triangulateCellsU(U);
    
    //========== Path of the correct OF instance ==========
    std::string obs_file_path = globalPath + "/obs_coordinates.txt";
    if (cwipiVerbose == 1) Info<< "The path where I have my observation coordinates is " << obs_file_path << endl;
    char UCfIntGen[500];
    strcpy(UCfIntGen, UCfIntPath.c_str());
    strcat(UCfIntGen, "/UCfInt");

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    label nMyProc = Pstream::myProcNo();
    int ens = round((myGlobalRank - nMyProc - 1)/nbParts) + 1;
    std::string ensString = std::to_string(ens);
    char ensChar[500];
    strcpy(ensChar, ensString.c_str());
    char* UCfIntEns_path = strcat(UCfIntGen, ensChar);

    if (cwipiVerbose == 1) Info<< "The file where UCfInt is created is " << UCfIntEns_path << endl;

    // ========= Configuration of processors to open obs_coordinates.txt =========
    std::ifstream obs_file;
    obs_file.open(obs_file_path);

    if (cwipiVerbose == 1) Info<< "Creating the sampling files..." << endl;

    point pointCoord;
    double Ux[cwipiObsU] = {0}; 
    double Uy[cwipiObsU] = {0};
    double Uz[cwipiObsU] = {0};
    int cellID[cwipiObsU] = {0};
    int cellIndex[cwipiObsU] = {0};
    int flagCellID[cwipiObsU] = {0};

    if (obs_file.is_open())
    {
        if (cwipiVerbose == 1) Info<< "Inside the observation coordinates file for interpolation" << endl;
        std::string line;
        int countProbes = 0;
        while( std::getline(obs_file,line) )
        {
            std::stringstream ss(line);

            std::string coordX, coordY, coordZ;
            std::getline(ss, coordX,',');
            pointCoord.x() = std::stod(coordX);
            std::getline(ss, coordY,',');
            pointCoord.y() = std::stod(coordY);
            std::getline(ss, coordZ,','); 
            pointCoord.z() = std::stod(coordZ);

            label ownCell = mesh.findCell(pointCoord);
            if (ownCell != -1){
                vector UatSP(vector::zero);
                UatSP = triangulateCellsU.interpolate(pointCoord, ownCell, -1);

                Ux[countProbes] = UatSP.x();
                Uy[countProbes] = UatSP.y();
                Uz[countProbes] = UatSP.z();
                cellID[countProbes] = ownCell;
                
                // Index in the list of observations, we do not use 0 because it corresponds to 
                // no cell in the domain
                cellIndex[countProbes] = countProbes+1;
                flagCellID[countProbes] = 1;
            }
            ++countProbes;
        }
    }
    else{
        FatalErrorInFunction
            << "In UCdInterpolation() function, "
            << "obs_coordinates.txt files has not been found."
            << Foam::abort(FatalError);
    }
    obs_file.close();

    //========= Gather all the values from the procs and then write the new values in a txt file =========//

    if (Pstream::master()){
        /* Gather the values */
        int tempCellIndex = 0;
        int tempCellID = 0;
        int tempflagCellID = 0;
        double tempUx = 0; 
        double tempUy = 0;
        double tempUz = 0;

        for(int i = 1; i < nbParts; i++)
        {
            IPstream streamFromMain(Pstream::commsTypes::blocking, i);
            for (int j = 0; j < cwipiObsU; ++j)
            {
                streamFromMain >> tempCellIndex;
                streamFromMain >> tempCellID;
                streamFromMain >> tempflagCellID;
                streamFromMain >> tempUx;
                streamFromMain >> tempUy;
                streamFromMain >> tempUz;

                cellIndex[j] = cellIndex[j] + tempCellIndex;
                cellID[j] = cellID[j] + tempCellID;
                flagCellID[j] = flagCellID[j] + tempflagCellID;
                Ux[j] = Ux[j] + tempUx;
                Uy[j] = Uy[j] + tempUy;
                Uz[j] = Uz[j] + tempUz;

                //========== In case of probes located over processor Boundary Conditions, we are summing up
                // the value twice, so we divide by 2 to take it into account only once ==========
                if (flagCellID[j] > 1){
                    Ux[j] = Ux[j] / flagCellID[j];
                    Uy[j] = Uy[j] / flagCellID[j];
                    Uz[j] = Uz[j] / flagCellID[j];
                }  
            }

            Pout<< "Information from processor " << i << " retrieved correctly" << endl;
        }

        if (cwipiVerbose == 1) Info<< "Writing the sampling file in the interpolation function" << endl;
        /* Write them */
        remove(UCfIntEns_path);
        std::ofstream UCfIntEns_file;
        UCfIntEns_file.open(UCfIntEns_path, std::ios::out | std::ios::app);
        if (UCfIntEns_file.is_open()){

            for (int i = 0; i < cwipiObsU; ++i){
                
                UCfIntEns_file << Ux[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                
                UCfIntEns_file << Uy[i] << "\n";
            }
            for (int i = 0; i < cwipiObsU; ++i){
                
                UCfIntEns_file << Uz[i] << "\n";
            }
            UCfIntEns_file << dragCoeff << "\n";
            UCfIntEns_file << liftCoeff << "\n";
        }
        else{
            FatalErrorInFunction
                << "Not able to create/open the file UCfInt."
                << Foam::abort(FatalError);
        }
        UCfIntEns_file.close();
    }
    else{
        // create the stream to send to the main proc
        OPstream stream2Main(Pstream::commsTypes::blocking, 0);
        for (int i = 0; i < cwipiObsU; ++i){
            stream2Main << cellIndex[i];
            stream2Main << cellID[i];
            stream2Main << flagCellID[i];
            stream2Main << Ux[i];
            stream2Main << Uy[i];
            stream2Main << Uz[i];
        }
    }

    //========== Retrieve the cellIDs of my probes ==========//
    if (myGlobalRank == 1 && (runTime.value() - runTime.deltaTValue()*cwipiStep) < 1e-4){
        std::ofstream myfile_2;
        std::string cellSensors = globalPath + "/results/cellIDs";

        remove(cellSensors.c_str());
        myfile_2.open(cellSensors, std::ios::out | std::ios::app);
        for (int i = 0; i < cwipiObsU; ++i){
            myfile_2 << cellID[i] << "\n";
        }
        myfile_2.close();
    }

    if (cwipiVerbose == 1) Info<< "Sampling files created " << endl;
}

}