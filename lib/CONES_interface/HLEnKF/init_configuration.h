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
 * @file init_configuration.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Configuration initialization.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
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


int cwipiStep = round(conesDict.lookupOrDefault<scalar>("observationWindow", 1));                           // Number of time step between each DA phase 
int cwipiMembers = round(conesDict.lookupOrDefault<scalar>("ensemble", 2));                                 // Number of members in the ensemble
int subdomains = round(decomposeParDict.lookup<scalar>("numberOfSubdomains"));                              // Number of subdomains    
int cwipiObsU = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesVelocity", 0));            // Number of observation probes for velocity       
int cwipiObsp = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesPressure", 0));            // Number of observation probes for pressure       
int cwipiObsCf = round(observationSubDict.lookupOrDefault<scalar>("numberObsProbesCf", 0));                 // Number of observation probes for friction coefficient (default = 1) 
int cwipiParams = round(conesDict.lookup<scalar>("numberParameters"));                                      // Number of parameters to optimize with the DA     
double geometricTolerance = conesDict.lookupOrDefault<scalar>("geometricTolerance", 0.1);                   // Geometric tolerance of the cwipi coupling meshes
int cwipiOutputNb = round(conesDict.lookupOrDefault<scalar>("numberOutputs", 0));                           // Number of txt file written beginning from the last one
Foam::vector sigmaUserU = observationSubDict.lookupOrDefault<vector>("velocityObsNoise", Foam::vector(0.05,0.05,0.05));                 // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for velocity)
double sigmaUserp = observationSubDict.lookupOrDefault<scalar>("pressureObsNoise", 0.05);                   // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for pressure)
double sigmaUserCf = observationSubDict.lookupOrDefault<scalar>("cfObsNoise", 0.05);                        // sigma of the EnKF (pertubation of the obs and diagonal of R matrix for friction coefficient)
int cwipiVerbose = round(conesDict.lookupOrDefault<scalar>("verbosityLevel", 0));                           // Print all the debuging messages or not: 1 = printed, 0 = nothing
bool cwipiTimedObs = observationSubDict.lookupOrDefault<bool>("obsTimeDependency", false);                  // Switch to for the obs: 1 = obs depends on time, 0 = obs does not depend on time
double obsTimeStep = observationSubDict.lookup<scalar>("obsTimeStep");                                      // The time step of the observations if the case is unsteady  
double obsStartTime = observationSubDict.lookup<scalar>("obsStartTime");                                    // The start time for the observation data 
Foam::word obsType = observationSubDict.lookupOrDefault<word>("obsType", "U");                              // Definition of observation parameter : 0 = vel, 1 = pres, 2 = both, 3 = vel+cf   
double stateInfl = inflationSubDict.lookupOrDefault<scalar>("stateInflation", 0);                           // State inflation (default = 0)  
double preStateInfl = inflationSubDict.lookupOrDefault<scalar>("preStateInflation", 0);                     // State inflation applied before the analysis, considered as model error (default = 0)  
double paramsInfl = inflationSubDict.lookupOrDefault<scalar>("parametersInflation", 0);                     // Parameters inflation (default = 0)      
Foam::word typeInfl = inflationSubDict.lookupOrDefault<word>("inflationType", "stochastic");                // Definition of inflation (0 = stochastic, 1 = deterministic)
bool clippingSwitch = localizationSubDict.lookupOrDefault<bool>("clippingSwitch" , false);                  // Switch to use the clipping file  
Foam::word clippingFile = localizationSubDict.lookupOrDefault<word>("clippingFile", "clippingCells.txt");   // Name of the file containing the clippings IDs
bool localSwitch = localizationSubDict.lookupOrDefault<bool>("covarianceLocalizationSwitch", false);        // Switch for the localisation (basic clipping or hyperlocalisation activated)
bool paramEstSwitch = conesDict.lookupOrDefault<bool>("paramEstSwitch", false);                             // Switch for parameter estimation
bool stateEstSwitch = conesDict.lookupOrDefault<bool>("stateEstSwitch", true);                              // Switch for state estimation           
Foam::vector obsVelocityComponents = observationSubDict.lookupOrDefault<vector>("obsVelocityComponents", Foam::vector(1,1,1));
Foam::word typeInputs = observationSubDict.lookupOrDefault<word>("obsNoiseType", "absVal");                 // Inputs for R are given in absolute values (0), percentage (1) or potential function (2) (default = 0)
Foam::vector correlationScale = localizationSubDict.lookup<vector>("correlationScale");
double sigmaLocX = correlationScale.x();                                                                    // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in X direction)
double sigmaLocY = correlationScale.y();                                                                    // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in X direction)
double sigmaLocZ = correlationScale.z();                                                                    // eta of the EnKF (pertubation of the Kalman gain to take into consideration the localization in X direction)    
double epsilon = conesDict.lookupOrDefault<scalar>("epsilon", 1e-10);                                       // Value used to calculate the NSRMD and comparing with 0
Foam::word obsCoordFile = observationSubDict.lookupOrDefault<word>("obsCoordFile", "obs_coordinates.txt");  // Name of the file containing the observation coordinates
Foam::word obsDataFile = observationSubDict.lookupOrDefault<word>("obsDataFile", "obs_field.txt");          // Name of the file containing the observation data