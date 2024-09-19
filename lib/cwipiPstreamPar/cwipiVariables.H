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
 * @file cwipiVariables.H
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Variables used for CWIPI.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
if (Pstream::master()) Foam::Pout << "Read cwipi variables" << endl;

Switch cwipiSwitch = runTime.controlDict().lookup("cwipiSwitch");

//========== Parameters from config file ==========
  /**
   * @brief conesDict dictionnary declaration.
   * 
   */
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

  /**
   * @brief decomposeParDict dictionnary declaration.
   * 
   */
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
  Foam::dictionary observationSubDict = conesDict.subDict("observationSubDict"); /*!< Observation sub-dictionnary declaration from conesDict dictionnary.*/
  Foam::dictionary inflationSubDict = conesDict.subDict("inflationSubDict"); /*!< Inflation sub-dictionnary declaration from conesDict dictionnary.*/
  Foam::dictionary geometricToleranceSubDict = conesDict.subDict("geometricToleranceSubDict"); /*!< Geometric tolerance sub-dictionnary declaration from conesDict dictionnary.*/
  

  int cwipiStep = round(conesDict.lookupOrDefault<scalar>("observationWindow", 1)); /*!< Number of time step between each DA phase */
  int subdomains = round(decomposeParDict.lookup<scalar>("numberOfSubdomains")); /*!< Number of decomposition subdomains from decomposeParDict dictionnary. */
  int cwipiObsU = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesVelocity", 0)); /*!< Number of velocity observation probes. */
  int cwipiObsp = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesPressure", 0)); /*!< Number of pressure observation probes. */
  int cwipiParams = round(conesDict.lookup<scalar>("numberParameters")); /*!< Number of parameters. */
  int cwipiMembers = round(conesDict.lookup<scalar>("ensemble")); /*!< Number of coarse ensemble members.*/
  int fineCwipiMembers = round(conesDict.lookup<scalar>("fineEnsemble")); /*!< Number of fine ensemble members. */
  double geometricToleranceFine = geometricToleranceSubDict.lookupOrDefault<scalar>("geometricToleranceFine", 0.1); /*!< Coarse geometric tolerance.*/
  double geometricToleranceCoarse = geometricToleranceSubDict.lookupOrDefault<scalar>("geometricToleranceCoarse", 0.9); /*!< Fine geometric tolerance. */
  int cwipiVerbose = round(conesDict.lookupOrDefault<scalar>("verbosityLevel", 0)); /*!< Verbosity level. */
  Foam::word obsType = observationSubDict.lookupOrDefault<word>("obsType", "U"); /*!< Observation type : velocity, pressure or both. */
  Foam::word typeInputs = observationSubDict.lookupOrDefault<word>("obsNoiseType", "absVel"); /*!< Type of noise : absolute or relative. */
  Foam::word typeInfl = inflationSubDict.lookupOrDefault<word>("inflationType", "stochastic"); /*!< Inflation type : stochastic or deterministic. */
  Foam::word obsCoordFile = observationSubDict.lookupOrDefault<word>("obsCoordFile", "obs_coordinates.txt");  /*!< Name of the file containing the observation coordinates. */ 

if (cwipiVerbose) if (Pstream::master())  Foam::Pout << "Information about decomposition read" << endl;

int principalControlMembers = fineCwipiMembers + fineCwipiMembers; /*!< Principal + control members.*/
int allMembers = cwipiMembers + fineCwipiMembers + fineCwipiMembers; /*!< Principal + control + ancillary members.*/

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
int myGlobalRank; /*!< Global rank. */
MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);

int myGlobalSize; /*!< Global size. */
MPI_Comm_size(MPI_COMM_WORLD, &myGlobalSize);

std::string globalCaseName = args.globalCaseName(); /*!< Current case name. */
std::string globalRootPath = args.rootPath(); /*!< Current root path.*/
std::string globalCasePath = globalRootPath + "/" + globalCaseName; /*!< Current case path.*/

/*
Foam::Info << "globalCaseName : " << globalCaseName << Foam::endl;
Foam::Info << "globalRootPath : " << globalRootPath << Foam::endl;
Foam::Info << "globalCasePath : " << globalCasePath << Foam::endl;
*/

if (cwipiVerbose) if (Pstream::master()){
    Foam::Pout<< "Solving with rank " << myGlobalRank << " the testcase with path " << globalCasePath << endl;
}

