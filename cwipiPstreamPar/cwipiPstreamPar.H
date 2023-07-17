/**
 * \file cwipiPstreamPar.H
 * \brief Declares all functions needed to employ CWIPI routines
 */

#include "volFields.H"
#include "interpolationCellPointWallModified.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime);

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, float cwipiVerbose, double geom_tol);

void UInterpolation(volVectorField& U, fvMesh& mesh, Time& runTime, int cwipiObsU, int mainsubDomain, interpolationCellPointWallModified<vector> triangulateCellsU,
float cwipiVerbose, std::string globalPath, std::string UIntPath);

void cwipiSend(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, float cwipiVerbose);

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, float cwipiVerbose);

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, float cwipiVerbose);

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, float cwipiVerbose);
}

// ************************************************************************* //
