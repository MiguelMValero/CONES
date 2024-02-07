/**
 * \file cwipiPstreamPar.H
 * \brief Declares all functions needed to employ CWIPI routines
 */

#include "volFields.H"
#include "interpolationCellPointWallModified.H"
#include "kinematicMomentumTransportModel.H"
// #include "interpolationCellPointFace.H"
// #include <cwipi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime);

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, 
int c2fconnec_size, int fconnec_size, int nbParts, float cwipiVerbose, double geom_tol);

void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int nbParts, float cwipiVerbose, std::string globalPath, std::string UIntPath);

void pInterpolation(volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string pIntPath);

void UpInterpolation(volVectorField& U, volScalarField& p, fvMesh& mesh, Time& runTime, int cwipiObsU, int cwipiObsp, int nbParts, float cwipiVerbose, std::string globalPath, std::string UpIntPath);

void cwipiSend(const fvMesh& mesh, const volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParams_sin(const fvMesh& mesh, const volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsChannel(const fvMesh& mesh, const volScalarField& Ck, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiSendParams_OFRsin(const fvMesh& mesh, const volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, int nbParts, float cwipiVerbose);

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParams_sin(fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsChannel(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecvParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, float cwipiVerbose, std::string globalRootPath);

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int nbParts, float cwipiVerbose);
}

// ************************************************************************* //
