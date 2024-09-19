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
 * @file probe.h
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Probe class declaration.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#ifndef PROBE_H    //include guards to prevent the "redefinition" of class error
#define PROBE_H

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <vector>
#include <numeric>  
#include <algorithm>

class probe
{
private:
    /* data */
    int nVar_;                        // Number of variables per probe (U,V,W,P = 4)
    std::string varLabels_;     // Labels of the variables (U,V,W,P)
    std::vector<double> probeCoordinates_;     // Labels of the variables (U,V,W,P)
    std::vector<std::vector<double>> probeData_;         // Data associated with the probe, rows = components (U,V,W,P), columns = time if time dependency

public:
    /** Member functions **/
    /* Constructors */
    probe(/* Default */);
    probe(int nVar, std::string varLabels, std::vector<double> probeCoordinates, std::vector<std::vector<double>> probeData);

    /* Destructor */
    ~probe();

    /* Getters */
    int get_nVar() const{return nVar_;};
    std::string get_varLabels() const{return varLabels_;};
    std::vector<double> get_probeCoordinates() const{return probeCoordinates_;}; 
    std::vector<std::vector<double>> get_probeData() const{return probeData_;};

    /* Setters */
    // No setters at the moment
};


#endif // PROBE_H


