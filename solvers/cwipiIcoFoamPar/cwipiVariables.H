if (Pstream::master()) Foam::Pout << "Read cwipi variables" << endl;

Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");

std::ifstream cFile ("cwipiConfig");

int values[34]={0};
int k = 0;

// if (cFile.is_open())
// {
//     std::string line;
//     while(std::getline(cFile, line))
//     {
//         line.erase(std::remove_if(line.begin(), line.end(), ::isspace),line.end());
//         if( line.empty() || line[0] == '#' )
//         {
//             continue;
//         }
//         auto delimiterPos = line.find("=");
//         // auto name = line.substr(0, delimiterPos);
//         auto value = line.substr(delimiterPos + 1);

//         values[k] = stod(value);
//         // std::cout << values[k] << '\n';
//         k = k+1;
//     }
// }
// else 
// {
//     std::cerr << "Couldn't open config file for reading.\n";
// }

// int cwipiStep = values[0];
// int nbParts = values[2];
// int cwipiObsU = values[3];
// int cwipiObsp = values[4];
// int cwipiParams = values[6];
// double geom_tol = values[7];
// float cwipiVerbose = values[12];
// int cwipiParamsObs = values[15];

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

  int cwipiStep = round(conesDict.lookupOrDefault<scalar>("observationWindow", 1));         // Number of time step between each DA phase 
  int cwipiMembers = round(conesDict.lookupOrDefault<scalar>("ensemble", 2));
  int subdomains = round(decomposeParDict.lookup<scalar>("numberOfSubdomains"));
  int cwipiObsU = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesVelocity", 0));
  int cwipiObsp = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesPressure", 0));
  int cwipiObsCf = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesCf", 0));
  int cwipiParams = round(conesDict.lookup<scalar>("numberParameters"));
  double geometricTolerance = conesDict.lookupOrDefault<scalar>("geometricTolerance", 0.1);
  int cwipiOutputNb = round(conesDict.lookupOrDefault<scalar>("numberOutputs", 0));
  double sigmaUserU = observationSubDict.lookupOrDefault<scalar>("velocityObsNoise", 0.05);
  double sigmaUserp = observationSubDict.lookupOrDefault<scalar>("pressureObsNoise", 0.05);
  double sigmaUserCf = observationSubDict.lookupOrDefault<scalar>("cfObsNoise", 0.05);
  int cwipiVerbose = round(conesDict.lookupOrDefault<scalar>("verbosityLevel", 0));
  bool cwipiTimedObs = observationSubDict.lookupOrDefault<bool>("obsTimeDependency", false);
  double obsTimeStep = observationSubDict.lookup<scalar>("obsTimeStep");
  Foam::word obsType = observationSubDict.lookupOrDefault<word>("obsType", "U");
  double stateInfl = inflationSubDict.lookupOrDefault<scalar>("stateInflation", 0);
  double paramsInfl = inflationSubDict.lookupOrDefault<scalar>("parametersInflation", 0);
  Foam::word typeInfl = inflationSubDict.lookupOrDefault<word>("inflationType", "stochastic");
  bool clippingSwitch = localizationSubDict.lookupOrDefault<bool>("clippingSwitch" , false);
  bool localSwitch = localizationSubDict.lookupOrDefault<bool>("covarianceLocalizationSwitch", false);
  bool hyperlocSwitch = localizationSubDict.lookupOrDefault<bool>("hyperlocalizationSwitch", false);
  bool paramEstSwitch = conesDict.lookupOrDefault<bool>("paramEstSwitch", false);
  bool stateEstSwitch = conesDict.lookupOrDefault<bool>("stateEstSwitch", true);
  Foam::vector obsVelocityComponents = observationSubDict.lookupOrDefault<vector>("obsVelocityComponents", Foam::vector(1,1,1));
  int Ux = round(obsVelocityComponents.x());
  int Uy = round(obsVelocityComponents.y());
  int Uz = round(obsVelocityComponents.z());
  Foam::word typeInputs = observationSubDict.lookupOrDefault<word>("obsNoiseType", "absVel");
  Foam::vector correlationScale = localizationSubDict.lookup<vector>("correlationScale");
  double sigmaLocX = correlationScale.x();
  double sigmaLocY = correlationScale.y();
  double sigmaLocZ = correlationScale.z();
  double epsilon = conesDict.lookupOrDefault<scalar>("epsilon", 1e-10);

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

if (cwipiVerbose) if (Pstream::master()){
    Foam::Pout<< "Solving with rank " << myGlobalRank << " the testcase with path " << globalCasePath << endl;
}
