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
 * @dir ./lib/CONES_interface/HLEnKF/src/DAClasses/probe
 * @brief **Probes class declarations and definitions.**
 */

/**
 * @file probe.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Probe class definition.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#include "probe.h"

/*** Region Class ***/
/** Member functions **/
/* Constructors */
probe::probe(/* Default */)
{}

probe::probe(int nVar, std::string varLabels, std::vector<double> probeCoordinates, std::vector<std::vector<double>> probeData):
nVar_(nVar),
varLabels_(varLabels),
probeCoordinates_(probeCoordinates),
probeData_(probeData)
{}

/* Destructor */
probe::~probe()
{}
