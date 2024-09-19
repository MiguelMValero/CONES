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
 * @file cwipiPstreamPar.H
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Defines all functions needed to employ CWIPI routines
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "volFields.H"
#include "interpolationCellPointWallModified.H"
#include "kinematicMomentumTransportModel.H"
// #include "interpolationCellPointFace.H"
// #include <cwipi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void tic(int mode=0);

void toc();

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime);

void deleteControlParams();


/* COUPLING */
void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, 
int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol);

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int nbParts, float cwipiVerbose);

void cwipiCouplingMGEnKF(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol_fine, double geom_tol_coarse, int subdomains, int cwipiMembers);

void cwipiCouplingProjection(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol);

void cwipiCouplingMFEnKF(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol_fine, double geom_tol_coarse, int cwipiMembers);


/* INTERPOLATION */
void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int nbParts, float cwipiVerbose, std::string globalPath, std::string UIntPath, const std::string filename);

void pInterpolation(volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string pIntPath, const std::string filename);

void UpInterpolation(volVectorField& U, volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsU, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string UpIntPath, const std::string filename);

void UCfInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, scalar dragCoeff, scalar liftCoeff, int nbParts, int cwipiStep, int cwipiVerbose, std::string globalPath, std::string UIntPath);


/* SEND */
void cwipiSendVolVectorField(const fvMesh& mesh, const volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiSendVolVectorFieldMGEnKF(const fvMesh& mesh, const volVectorField& stateMatricesDiff, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiSendVolVectorFieldMFEnKF(const fvMesh& mesh, const volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsAirfoil(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsMGEnKF(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParams_sin(const fvMesh& mesh, const volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsChannel(const fvMesh& mesh, const volScalarField& Ck, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParams_OFRsin(const fvMesh& mesh, const volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsCylinder_IBM(const fvMesh& mesh, int cwipiParams, int nbParts, int cwipiVerbose, std::string parPath);


/* RECEIVE */

void cwipiRecvVolVectorField(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiRecvVolVectorFieldMGEnKF(const fvMesh& mesh, volVectorField& U, volVectorField& stateMatricesDiff, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiRecvVolVectorFieldMFEnKF(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsAirfoil(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsMGEnKF(const fvMesh& mesh, volVectorField& stateMatrices, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParams_sin(fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsChannel(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsCylinder_IBM(int cwipiParams, double* paramsToRecv, int cwipiVerbose, std::string parPath);





// Foam::vector cwipiRecvParamsAvg(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

}



// ************************************************************************* //
