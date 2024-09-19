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
 * @file observations.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Observation class declaration.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#ifndef OBSERVATIONS_H    //include guards to prevent the "redefinition" of class error
#define OBSERVATIONS_H

#include "fvCFD.H"
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <vector>
#include <numeric>  
#include <algorithm>
#include "../probe/probe.h"

class observations
{
private:
    /* data */
	// Number of observations per type
	int nObsProbesVel_;                  // Number of observation probes for velocity
	int nObsProbesPres_;                 // Number of observation probes for pressure
	int nObsProbesCf_;                   // Number of observation probes for friction coeff
	int nObsProbes_;                     // Total Number of observation probes
	std::string obsType_;                // Definition of observation parameter :  U (velocity only), p (pressure only), Up (velocity and pressure), UCf (velocity and friction coefficient)
	// std::vector<bool> velComponents_; // Specify velocity components in the observation
	int velocityCaseSwitch_;             // Switch to define the velocity case depending on which variable we consider
	int nObs_;                           // Total Number of observations

	// Noise configuration
	double std_u_;                        // Observation noise / confidence for velocity
	double std_v_;                        // Observation noise / confidence for velocity
	double std_w_;                        // Observation noise / confidence for velocity
	double std_p_;                        // Observation noise / confidence for pressure
	double std_Cf_;                       // Observation noise / confidence for friction coeff
	std::string obsNoiseType_;           // Inputs for R are given in absolute values (absVal), percentage (relVal)
	
	// Time configuration
	bool timeDependency_;                // Type of the observations in terms of Steady/Unsteady
	double obsTimeStep_;                 // The time step of the observations if the case is unsteady  
	double obsStartTime_;                // The start time for the observation data  

	// Raw information
	std::vector<std::vector<double>> obs_coordinates_;
	std::vector<std::vector<double>> obs_raw_data_;
	
	// Define a list of probes
	std::vector<probe> obsProbes_;        // [probe_Data(nVar, varLabels, raw_data[[ii, (ii+nProbes)], :]) for ii in range(nProbes)]

public:
    observations(/* Default */);
	observations(int nObsProbesVel, int nObsProbesPres, int nObsProbesCf, Foam::vector obsVelocityComponents, Foam::vector noiseVel, double noisep, double noiseCf, std::string obsNoiseType, bool timeDependency, double obsTimeStep, double obsStartTime, std::string obsType, const std::string dir_in, const std::string coordFilename, const std::string dataFilename);
    ~observations();
		
		/* Set methods */
		void init_nObs();
		void init_velocityCaseSwitch(const Foam::vector obsVelocityComponents);
		void init_obsProbes(const std::vector<std::vector<double>> obs_coordinates,const std::vector<std::vector<double>> raw_data);
		
		void set_obsCoord(const std::vector<std::vector<double>> coords);
		void set_obsData(const std::vector<std::vector<double>> data);

		//*Not used atm*
		// void set_obsType(const std::string obsType);
		// void set_std(const double stdVal, const std::string obsType);
		// void set_stdType(const std::string stdType);

		/* Getter methods */
		int get_nObsProbesVel() const{return nObsProbesVel_;};
		int get_nObsProbesPres() const{return nObsProbesPres_;};
		int get_nObsProbesCf() const{return nObsProbesCf_;};
		int get_nObsProbes() const{return nObsProbes_;};
		std::string get_obsType() const{return obsType_;};
		int get_velocityCaseSwitch() const{return velocityCaseSwitch_;};
		int get_nObs() const{return nObs_;};
		double get_std_u() const{return std_u_;};
		double get_std_v() const{return std_v_;};
		double get_std_w() const{return std_w_;};
		double get_std_p() const{return std_p_;};
		double get_std_Cf() const{return std_Cf_;};
		std::string get_obsNoiseType() const{return obsNoiseType_;};
		bool get_timeDependency() const{return timeDependency_;};
		double get_obsTimeStep() const{return obsTimeStep_;};
		double get_obsStartTime() const{return obsStartTime_;};
		std::vector<std::vector<double>> get_obs_coordinates() const{return obs_coordinates_;};
		std::vector<std::vector<double>> get_obs_raw_data() const{return obs_raw_data_;};

		/* Getter methods for probes info */
		probe get_probe(int probeIndex) const{return obsProbes_[probeIndex];};
		int get_probeVar(int probeIndex);
		std::string get_probeVarLabels(int probeIndex);
		std::vector<double> get_probeCoord(int probeIndex);
		std::vector<std::vector<double>> get_probeData(int probeIndex);


		/* info methods */
		void info_obsCoord();
		void info_obsData();
		void info_std();
		void info_timeDependency();
		void info_stdType();

		/* Other Methods */
		std::vector<std::vector<double>> read_obs_coordinates(const std::string dir_in, const std::string filename);
		std::vector<std::vector<double>> read_obs_data(const std::string dir_in, const std::string filename);
};

#endif // OBSERVATIONS_H


