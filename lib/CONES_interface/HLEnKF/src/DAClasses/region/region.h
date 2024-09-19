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
 * @file region.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Region class declaration.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#ifndef REGION_H    //include guards to prevent the "redefinition" of class error
#define REGION_H

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <vector>
#include <numeric>  
#include <algorithm>

#include "../observations/observations.h"

class region
{
private:
    /* data */
    int nClips_;                 // number of clips in the DA region
    Eigen::ArrayXi clipIDs_;     // IDs of the clips constituting the region
    int nCells_;                 // number of cells in the DA region
    Eigen::ArrayXi cellIDs_;     // IDs of the cells in the region
    int nProbes_;                // number of probes in the DA region
    Eigen::ArrayXi probeIDs_;    // IDs of the probes in DA region

    int nObsProbesVel_;                  // Number of observation probes for velocity in the region
	int nObsProbesPres_;                 // Number of observation probes for pressure in the region
	int nObsProbesCf_;                   // Number of observation probes for friction coeff in the region
	int nObs_;                           // Total Number of observations in the region

    Eigen::MatrixXd obsCoords_;           // Coordinates of all the observations values (not the probes) so some coordinates are repeated
    // Eigen::MatrixXd obsRawData_;          // Raw data with all the times so that we can access it faster during the calcutation
    
    
    //= TO DO =
    // Add the two probes attributes, for now it is not usefull since probeID = clipID and nProbes = nClips

public:
    /** Member functions **/
    /* Constructors */
    region();                    // Default constructor needed to declare non specified regions in a std::vector 
    region(Eigen::ArrayXi clipIDs, Eigen::ArrayXi cellIDs);
    region(Eigen::ArrayXi clipIDs, Eigen::ArrayXi cellIDs, Eigen::ArrayXi probeIDs);

    /* Destructor */
    ~region();

    /* Getters */
    int get_nClips() const{return nClips_;};
    Eigen::ArrayXi get_clipIDs() const{return clipIDs_;};
    int get_nCells() const{return nCells_;};
    Eigen::ArrayXi get_cellIDs() const{return cellIDs_;};
    int get_nProbes() const{return nProbes_;};
    Eigen::ArrayXi get_probesIDs() const{return probeIDs_;};
    Eigen::MatrixXd get_obsCoords() const{return obsCoords_;};
    // Eigen::MatrixXd get_obsRawdata(int timeIndex) const{return obsRawData_.col(timeIndex);};

    int get_nObsProbesVel() const{return nObsProbesVel_;};
	int get_nObsProbesPres() const{return nObsProbesPres_;};
	int get_nObsProbesCf() const{return nObsProbesCf_;};
	int get_nObs() const{return nObs_;};

    /* Setters */
    void set_nObsProbes(observations obsData);
    void set_nObsProbesNoClip(observations obsData);
};

/* Outside class but related functions */
Eigen::ArrayXi removeDuplicatesSet(const Eigen::Ref<const Eigen::ArrayXi>& arr);

Eigen::ArrayXi read_clipFile_info(std::string stringRootPath, std::string filename, float cwipiVerbose);

std::vector<int> calculate_nCells_clips(const Eigen::Ref<const Eigen::ArrayXi>& clippingFile, float cwipiVerbose);

std::vector<Eigen::ArrayXi> calculate_cellIDs_clips(const Eigen::Ref<const Eigen::ArrayXi>& clippingFile, float cwipiVerbose);

std::vector<region> clips_to_DA_regions(std::vector<int> nCellsClips, std::vector<Eigen::ArrayXi> clipsCells, float cwipiVerbose);

std::vector<region> initialize_DA_regions(const std::string stringRootPath, const std::string filename, float cwipiVerbose, bool clippingSwitch, int nb_cells, const observations obsData);


#endif // REGION_H