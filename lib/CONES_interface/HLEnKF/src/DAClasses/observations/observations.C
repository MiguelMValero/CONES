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
 * @dir ./lib/CONES_interface/HLEnKF/src/DAClasses
 * @brief **Data Assimilation HLEnKF classes declarations and definitions.**
 */

/**
 * @dir ./lib/CONES_interface/HLEnKF/src/DAClasses/observations
 * @brief **Observations class declarations and definitions.**
 */

/**
 * @file observations.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Observation class definition.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#include "observations.h"

/*** Observations Class ***/
/** Member functions **/
/* Constructors */
observations::observations(/* Default */)
{
}

observations::observations(int nObsProbesVel, int nObsProbesPres, int nObsProbesCf, Foam::vector obsVelocityComponents, Foam::vector noiseVel, double noisep, double noiseCf, std::string obsNoiseType, bool timeDependency, double obsTimeStep, double obsStartTime, std::string obsType, const std::string dir_in, const std::string coordFilename, const std::string dataFilename):
nObsProbesVel_(nObsProbesVel),
nObsProbesPres_(nObsProbesPres),
nObsProbesCf_(nObsProbesCf),
nObsProbes_(nObsProbesVel+nObsProbesPres+nObsProbesCf),
obsType_(obsType),

std_u_(noiseVel.x()),
std_v_(noiseVel.y()),
std_w_(noiseVel.z()),
std_p_(noisep),
std_Cf_(noiseCf),
obsNoiseType_(obsNoiseType),

timeDependency_(timeDependency),
obsTimeStep_ (obsTimeStep),
obsStartTime_(obsStartTime)
{
  //Define the list of probes inside the constructor
  init_velocityCaseSwitch(obsVelocityComponents);
  std::cout << " velocityCaseSwitch is " << velocityCaseSwitch_ << std::endl << "\n";

  init_nObs();
  std::cout << " nObs is " << nObs_ << std::endl << "\n";

  std::vector<std::vector<double>> obs_coordinates = read_obs_coordinates(dir_in, coordFilename);
  std::vector<std::vector<double>> raw_data = read_obs_data(dir_in, dataFilename);

  set_obsCoord(obs_coordinates);
  set_obsData(raw_data);

  init_obsProbes(obs_coordinates, raw_data);

  std::cout << " Initialization of observation probes done " << std::endl << "\n";
}


/* Destructor */
observations::~observations()
{
}


/* Setter methods */
void observations::init_nObs(){
	if (obsType_ == "U")
	{
		switch (velocityCaseSwitch_)
		{
		case 1:
		case 2:
		case 3:
			nObs_ = nObsProbesVel_;
		break;
		case 4:
		case 5:
		case 6:
			nObs_ = 2 * nObsProbesVel_;
		break;
		case 7:
			nObs_ = 3 * nObsProbesVel_;
		break;
		}
	}
	else if (obsType_ == "p")
		nObs_ = nObsProbesPres_;
	else if (obsType_ == "Up")
	{
		switch (velocityCaseSwitch_)
		{
		case 1:
		case 2:
		case 3:
			nObs_ = nObsProbesVel_ + nObsProbesPres_;
		break;
		case 4:
		case 5:
		case 6:
			nObs_ = 2 * nObsProbesVel_ + nObsProbesPres_;
		break;
		case 7:
			nObs_ = 3 * nObsProbesVel_ + nObsProbesPres_;
		break;
		}
	}
	else if (obsType_ == "UCf")
	{
		switch (velocityCaseSwitch_)
		{
		case 1:
		case 2:
		case 3:
			nObs_ = nObsProbesVel_ + nObsProbesCf_;
		break;
		case 4:
		case 5:
		case 6:
			nObs_ = 2 * nObsProbesVel_ + nObsProbesCf_;
		break;
		case 7:
			nObs_ = 3 * nObsProbesVel_ + nObsProbesCf_;
		break;
		}
	}
}

void observations::init_velocityCaseSwitch(const Foam::vector obsVelocityComponents){
  int Ux = round(obsVelocityComponents.x());           // Specification if Ux is read or not (cwipiParamsObs needs to be either 0, 2 or 3)
  int Uy = round(obsVelocityComponents.y());           // Specification if Uy is read or not (cwipiParamsObs needs to be either 0, 2 or 3)
  int Uz = round(obsVelocityComponents.z());           // Specification if Uz is read or not (cwipiParamsObs needs to be either 0, 2 or 3)  
  
  if (obsType_ != "p")
  {
    if (Ux == 1 && Uy == 0 && Uz == 0)
      velocityCaseSwitch_ = 1;
    else if (Ux == 0 && Uy == 1 && Uz == 0)
      velocityCaseSwitch_ = 2;
    else if (Ux == 0 && Uy == 0 && Uz == 1)
      velocityCaseSwitch_ = 3;
    else if (Ux == 1 && Uy == 1 && Uz == 0)
      velocityCaseSwitch_ = 4;
    else if (Ux == 1 && Uy == 0 && Uz == 1)
      velocityCaseSwitch_ = 5;
    else if (Ux == 0 && Uy == 1 && Uz == 1)
      velocityCaseSwitch_ = 6;
    else if (Ux == 1 && Uy == 1 && Uz == 1)
      velocityCaseSwitch_ = 7;
    else
    {
      std::cerr << "Check your inputs. Either obsType or the components of the velocity are not properly defined.\n";
    }
  }
}

void observations::init_obsProbes(const std::vector<std::vector<double>> obs_coordinates,const std::vector<std::vector<double>> raw_data)
{
	if (obsType_ == "U"){
		switch (velocityCaseSwitch_){
			case 1 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "u", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 2 :
      	for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "v", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 3 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "w", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 4 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uv", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 5 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uw", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 6 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
         			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "vw", obs_coordinates[i], obsData_temp));
				}
        		break;
			case 7 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(3);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
          			obsData_temp[2] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(3, "uvw", obs_coordinates[i], obsData_temp));
				}
        		break;
		}
	}
	else if (obsType_ == "p"){
      	for (int i=0 ; i<nObsProbesPres_; i++)
			{
				std::vector<std::vector<double>> obsData_temp(1);
				obsData_temp[0] = raw_data[i];
				obsProbes_.push_back(probe(1, "p", obs_coordinates[i], obsData_temp));
			}
	}
	else if (obsType_ == "Up"){
		switch (velocityCaseSwitch_){
			case 1 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "u", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
        		break;
			case 2 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "v", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
        		break;
			case 3 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "w", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
				break;
			case 4 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uv", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 5 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uw", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 6 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "vw", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "p", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 7 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(3);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
          			obsData_temp[2] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(3, "uvw", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesPres_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
          			obsProbes_.push_back(probe(1, "p", obs_coordinates[i+3*nObsProbesVel_], obsData_temp));
				}
        		break;
		}
	}
	else if (obsType_ == "UCf"){
		switch (velocityCaseSwitch_){
			case 1 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "u", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
          			obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
        		break;
			case 2 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "v", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
          			obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
        		break;
			case 3 :
      			for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i];
					obsProbes_.push_back(probe(1, "w", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+nObsProbesVel_];
          			obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+nObsProbesVel_], obsData_temp));
				}
				break;
			case 4 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uv", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 5 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "uw", obs_coordinates[i], obsData_temp));
				}
        for (int i=0; i<nObsProbesCf_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 6 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(2);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
					obsProbes_.push_back(probe(2, "vw", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+2*nObsProbesVel_], obsData_temp));
				}
        		break;
			case 7 :
        		for (int i=0 ; i<nObsProbesVel_; i++)
				{
					std::vector<std::vector<double>> obsData_temp(3);
					obsData_temp[0] = raw_data[i];
          			obsData_temp[1] = raw_data[i+nObsProbesVel_];
          			obsData_temp[2] = raw_data[i+2*nObsProbesVel_];
					obsProbes_.push_back(probe(3, "uvw", obs_coordinates[i], obsData_temp));
				}
        		for (int i=0; i<nObsProbesCf_; i++)
				{
          			std::vector<std::vector<double>> obsData_temp(1);
					obsData_temp[0] = raw_data[i+3*nObsProbesVel_];
					obsProbes_.push_back(probe(1, "Cf", obs_coordinates[i+3*nObsProbesVel_], obsData_temp));
				}
        		break;
		}
	}
}

void observations::set_obsCoord(const std::vector<std::vector<double>> coord){
	obs_coordinates_ = coord;
}

void observations::set_obsData(const std::vector<std::vector<double>> data){
	obs_raw_data_ = data;
}


//*Not used atm*

// void observations::set_obsType(const std::string obsType){
// 	obsType_ = obsType;
// }

// void observations::set_std(const double stdVal, const std::string obsType){
// 	if(obsType == "u"){
// 		std_u = stdVal;
// 	}else if(obsType == "v"){
// 		std_v = stdVal;
// 	}else if(obsType == "w"){
// 		std_w = stdVal;
// 	}else if(obsType == "p"){
// 		std_p = stdVal;
// 	}else if(obsType == "Cf"){
// 		std_Cf =stdVal;
// 	}else{
// 		std::cout << "Variable " << obsType << " has not been coded yet!\n";
// 	}
// }

// void observations::set_stdType(const std::string stdType){
// 	if(stdType == "absVal"){
// 		obsNoiseType_ = stdType;
// 	}else if(stdType == "relVal"){
// 		obsNoiseType_ = stdType;
// 	}else{
// 		std::cout << "Standard deviation type " << stdType << " has not been coded yet!\n";
// 	}
// }


/* Getter methods */
// See header file

/* Getter methods for probes info */
int observations::get_probeVar(int probeIndex)
{
  probe probe_temp = get_probe(probeIndex);
  return probe_temp.get_nVar();
}

std::string observations::get_probeVarLabels(int probeIndex)
{
  probe probe_temp = get_probe(probeIndex);
  return probe_temp.get_varLabels();
}

std::vector<double> observations::get_probeCoord(int probeIndex)
{
  probe probe_temp = get_probe(probeIndex);
  return probe_temp.get_probeCoordinates();
}

std::vector<std::vector<double>> observations::get_probeData(int probeIndex)
{
  probe probe_temp = get_probe(probeIndex);
  return probe_temp.get_probeData();
}


/* Info methods */
void observations::info_obsCoord(){
	for(unsigned int i = 0; i < obs_coordinates_.size() ; i++){
		std::cout << "\tObs coord " << i << "\n";
		for(unsigned int j = 0; j < obs_coordinates_[i].size(); j++){
			std::cout << "\t\t" << obs_coordinates_[i][j] << "\t";
		std::cout << "\n" ;
		}
	}
}

void observations::info_obsData(){
	for(unsigned int i = 0; i < obs_raw_data_.size() ; i++){
		std::cout << "\tObs raw data " << i << "\n";
		for(unsigned int j = 0; j < obs_raw_data_[i].size(); j++){
			std::cout << "\t\t" << obs_raw_data_[i][j] << "\t";
		std::cout << "\n" ;
		}
	}
}

void observations::info_std(){
  if (obsType_ == "U")
  {
    switch (velocityCaseSwitch_)
    {
    case 1:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "\n";
      break;
    case 2:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "\n";
      break;
    case 3:
      std::cout << obsType_ << " standard deviation set to " << std_w_ << "\n";
      break;
    case 4:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "\n";
      break;
    case 5:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_w_ << "\n";
      break;
    case 6:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "and " << std_w_ << "\n";
      break;
    case 7:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "and " << std_w_ << "\n";
      break;
    }
  }
  else if (obsType_ == "p")
    std::cout << obsType_ << " standard deviation set to " << std_p_ << "\n";
  else if (obsType_ == "Up")
  {
    switch (velocityCaseSwitch_)
    {
    case 1:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_p_ << "\n";
      break;
    case 2:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "and " << std_p_ << "\n";
      break;
    case 3:
      std::cout << obsType_ << " standard deviation set to " << std_w_ << "and " << std_p_ << "\n";
      break;
    case 4:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "and " << std_p_ << "\n";
      break;
    case 5:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_w_ << "and " << std_p_ << "\n";
      break;
    case 6:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "and " << std_w_ << "and " << std_p_ << "\n";
      break;
    case 7:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "and " << std_w_ << "and " << std_p_ << "\n";
      break;
    }
  }
  else if (obsType_ == "UCf")
  {
    switch (velocityCaseSwitch_)
    {
    case 1:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_Cf_ << "\n";
      break;
    case 2:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "and " << std_Cf_ << "\n";
      break;
    case 3:
      std::cout << obsType_ << " standard deviation set to " << std_w_ << "and " << std_Cf_ << "\n";
      break;
    case 4:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "and " << std_Cf_ << "\n";
      break;
    case 5:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_w_ << "and " << std_Cf_ << "\n";
      break;
    case 6:
      std::cout << obsType_ << " standard deviation set to " << std_v_ << "and " << std_w_ << "and " << std_Cf_ << "\n";
      break;
    case 7:
      std::cout << obsType_ << " standard deviation set to " << std_u_ << "and " << std_v_ << "and " << std_w_ << "and " << std_Cf_ << "\n";
      break;
    }
  }
}

void observations::info_timeDependency(){std::cout << "Time dependency set to " << timeDependency_ << "\n";}

void observations::info_stdType(){std::cout << "Standard deviation type set to " << obsNoiseType_ << "\n";}


/* Other Methods */
std::vector<std::vector<double>> observations::read_obs_coordinates(const std::string dir_in, const std::string filename){
	std::string concString = dir_in + '/' + filename;
	std::cout << "\n Reading Observations coordinates file" <<  concString << "\n";
	std::vector<std::vector<double>> obs_coordinates;
	int num_obs = 0;
	// Number of observations points
	std::ifstream f(concString);
	std::string line;
	int i;
	std::fstream ifs;
	std::string token;
	double cord;
	for (i = 0; std::getline(f, line); ++i){
		std::vector<double> c1;
		for( int j = 0 ; j < 3 ; j++){
			token = line.substr(0, line.find(","));
			line.erase(0, line.find(",") + 1);
			cord = std::stod(token);
			c1.push_back(cord);
		}
		obs_coordinates.push_back(c1);
	};
	num_obs = i;
	std::cout << "\n\tNumber of observation points " << num_obs << "\n";
	// for(i = 0 ; i < num_obs ; i++){
	// 	for(int j = 0; j < 3 ; j++){
	// 		std::cout << "\t\t"  << obs_coordinates[i][j];
	// 	}
	// 	std::cout << "\n" ;
	// }
	std::cout << std::endl;

	// set_obsCoord(obs_coordinates);
	// set_nObs(num_obs, obsType_);

	return obs_coordinates;
}

std::vector<std::vector<double>> observations::read_obs_data(const std::string dir_in, const std::string filename){
	//========== Depending on time : Read the observations data from a file named "obs_field.txt" only and return the values 
    // of the specific time step in a eigen array. The txt file has to contain velocity values with X Y and Z in the same column
    // (first all the X, then Y then Z). Number of column = number of DA phases ==========
    
    //=========== Extract observation data from file ===========
    std::string obs_file= dir_in + "/" + filename;
	std::cout << "\n Reading Observations data file" <<  obs_file << "\n";
 
    std::vector<std::vector<double>> obs_content;
    std::vector<double> obs_row;
    std::string obs_line, obs_word;
 
    std::fstream obs_stream (obs_file, std::ios::in);
    if(obs_stream.is_open())
    {
        while(getline(obs_stream, obs_line))
        {
            obs_row.clear();
 
            std::stringstream obs_str(obs_line);
 
            while(getline(obs_str, obs_word, ','))
                obs_row.push_back(std::stod(obs_word));
            obs_content.push_back(obs_row);
        }
    }
    else {
        std::cerr << "Couldn't open observation file.\n";
    }
    
    std::cout << "\t Obs raw data ok " << "\n";
    // for (unsigned int i=0; i<obs_content.size(); i++)
    // {
	// 	std::cout << "\t\t";
	// 	for (unsigned int j=0; j<obs_content[i].size(); j++)
	// 	{
    //     	std::cout << obs_content[i][j];
	// 	}
	// 	std::cout << "\n";
    // }
    std::cout << std::endl << "\n";

    return obs_content;
}

// std::vector<std::vector<double>> observations::read_obs_data(const std::string dir_in, const std::string filename){
// 	std::string concString = dir_in + '/' + filename;
// 	std::cout << "\n Reading Observations data file" <<  concString << "\n";
// 	std::vector<std::vector<double>> raw_data;
// 	int num_obs = nObsProbesVel_;
// 	if(obsType_ == "U"){
// 		num_obs = nObsProbesVel_;
// 	}else{
// 		std::cout << "obsType not defined!\n";
// 	}
// 	// Number of observations points
// 	std::ifstream f(concString);
// 	std::string line;
// 	int i;
// 	std::fstream ifs;
// 	std::string token;
// 	double data;
// 	int j = 0;
// 	if(!timeDependency_){
// 		std::cout << "\n\tSteady observations\n" ;
// 		if(i%num_obs != 0) std::cout << "Number of data is inconsistent with the number of observations" << "\n"; 
// 		std::vector<double> c1 = {0, 0, 0};
// 		for (i = 0; std::getline(f, line); i++){
// 			data = std::stod(line);
// 			c1[j] =data;
// 			std::cout << "\t" << c1[j]  ;
// 			j++;
// 			if(j==3){
// 				std::cout << "\n" ;
// 				j=0;
// 				raw_data.push_back(c1);
// 			}
// 		}
// 		std::cout << "\n\t Number of observation data " << i << "\n";
// 	}else{
// 		std::cout << "\t Unsteady observations\n" ; 
// 	}
// 	// int nVar_per_obs = i/num_obs;
// 	// set_nVar(nVar_per_obs);
// 	// std::cout << "\t" << nVar_per_obs << " variables per observations point\n" ;
// 	//for(i = 0 ; i < num_obs ; i++){
// 	//	for(int j = 0; j < 3 ; j++){
// 	//		std::cout << "\t\t"  << obs_coordinates[i][j] ;
// 	//	}
// 	//	std::cout << "\n" ;
// 	//}
// 	// set_obsData(raw_data);
// 	return raw_data;
// }

