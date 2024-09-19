if (Pstream::master()) Foam::Pout << "Read cwipi variables" << endl;

Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");

//========== Parameters from config file ==========
IOdictionary conesDict
(
IOobject
(
    "conesDict",
    runTime.system(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
)
);

IOdictionary decomposeParDict
(
IOobject
(
    "decomposeParDict",
    runTime.system(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
)
);

// Sub dict declaration
Foam::dictionary observationSubDict = conesDict.subDict("observationSubDict");
Foam::dictionary inflationSubDict = conesDict.subDict("inflationSubDict");
Foam::dictionary localizationSubDict = conesDict.subDict("localizationSubDict");

// int cwipiStep = values[0];
int cwipiStep = round(conesDict.lookupOrDefault<scalar>("observationWindow", 1));                         // Number of time step between each DA phase 
// int nbParts = values[2];
int subdomains = round(decomposeParDict.lookup<scalar>("numberOfSubdomains"));                            // Number of subdomains    
// int cwipiObsU = values[3];
int cwipiObsU = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesVelocity", 0));          // Number of observation probes for velocity  
// int cwipiObsp = values[4];
int cwipiObsp = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesPressure", 0));          // Number of observation probes for pressure 
// int cwipiParams = values[6];
int cwipiParams = round(conesDict.lookup<scalar>("numberParameters"));                                    // Number of parameters to optimize with the DA
// double geom_tol = values[7];
double geometricTolerance = conesDict.lookupOrDefault<scalar>("geometricTolerance", 0.1);                 // Geometric tolerance of the cwipi coupling meshes
// float cwipiVerbose = values[12];
int cwipiVerbose = round(conesDict.lookupOrDefault<scalar>("verbosityLevel", 0));                         // Print all the debuging messages or not: 1 = printed, 0 = nothing
// int cwipiParamsObs = values[15];
Foam::word obsType = observationSubDict.lookupOrDefault<word>("obsType", "U");
Foam::word obsCoordFile = observationSubDict.lookupOrDefault<word>("obsCoordFile", "obs_coordinates.txt");  /*!< Name of the file containing the observation coordinates. */ 

if (cwipiVerbose) if (Pstream::master()) Foam::Pout << "Information about decomposition read" << endl;

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

if (cwipiVerbose) if (Pstream::master){
    Foam::Pout<< "Solving with rank " << myGlobalRank << " the testcase with path " << globalCasePath << endl;
}
