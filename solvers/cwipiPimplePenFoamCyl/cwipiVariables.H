if (Pstream::master()) Foam::Pout << "Read cwipi variables" << endl;

Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");

std::ifstream cFile ("cwipiConfig");

int values[35]={0};
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
        // auto name = line.substr(0, delimiterPos);
        auto value = line.substr(delimiterPos + 1);

        values[k] = stod(value);
        // std::cout << values[k] << '\n';
        k = k+1;
    }
}
else 
{
    std::cerr << "Couldn't open config file for reading.\n";
}

int cwipiStep = values[0];
int nbParts = values[2];
int cwipiObsU = values[3];
int cwipiObsp = values[4];
int cwipiParams = values[6];
double geom_tol = values[7];
float cwipiVerbose = values[12];
int cwipiParamsObs = values[15];

if (cwipiVerbose) if (Pstream::master())  Foam::Pout << "Information about decomposition read" << endl;

scalar cwipiTimestep = 1;
scalar numberCwipiPhaseDouble = (runTime.endTime().value()/runTime.deltaTValue())/cwipiStep;
int numberCwipiPhase = round(numberCwipiPhaseDouble);
int cwipiIteration = 1;
scalar cwipiPhaseCheck = 0;

double* pointCoords = new double[3*mesh.nPoints()];

int* face_index = new int[mesh.nCells()+1];
int* face_connectivity_index = new int[mesh.nFaces()+1];

int c2fconnec_size = 0;
forAll(mesh.cells(),i){
    c2fconnec_size = c2fconnec_size + mesh.cells()[i].size();
}
int* cell_to_face_connectivity = new int[c2fconnec_size];

int fconnec_size = 0;
forAll(mesh.faces(),i){
    fconnec_size = fconnec_size + mesh.faces()[i].size();
}
int* face_connectivity = new int[fconnec_size];

//========== Definition of all MPI routines and directories ============//
int myGlobalRank;
MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

std::string globalCaseName = args.globalCaseName();
std::string globalRootPath = args.rootPath(); // Root path
std::string globalCasePath = globalRootPath + "/" + globalCaseName; // For UInt

// Get the start time and convert to string (with a good format) for creating a directory
char initialTimeString[10];  // maximum expected length of the float
std::snprintf(initialTimeString, 10, "%.0f", runTime.startTime().value());
std::string str(initialTimeString);

std::string postProcPath = globalRootPath + "/" + globalCaseName + "/postProcessing"; // postProcessing file path
std::string fCoeffIBMPath = globalRootPath + "/" + globalCaseName + "/postProcessing/forceCoeffs_IBM";
std::string writefCoeffPath = fCoeffIBMPath + "/" + initialTimeString; // Path where the forceCoeffs.txt file should be written

std::string command = "mkdir -p " + writefCoeffPath;
Foam::system(command);

if (cwipiVerbose) if (Pstream::master()){
    Foam::Pout<< "Solving with rank " << myGlobalRank << " the testcase with path " << globalCasePath << endl;
}
