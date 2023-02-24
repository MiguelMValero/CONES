/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class fvMesh;

void addControlParams(int numberCwipiPhase, double deltaT, double currentTime, int nbParts, int partsReparty);

void cwipiCoupling(const fvMesh& mesh, double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int c2fconnec_size, int fconnec_size, int cwipiVerbose, double geom_tol);

void UInterpolation(Foam::volVectorField& U, Foam::fvMesh& mesh, int cwipiObsU, int cwipiVerbose, std::string UIntPath);

void pInterpolation(Foam::volScalarField& p, Foam::fvMesh& mesh, int cwipiVerbose, std::string pIntPath);

void UpInterpolation(Foam::volVectorField& U, Foam::volScalarField& p, Foam::fvMesh& mesh, int cwipiObsU, int cwipiObsp, int cwipiVerbose, std::string UpIntPath);

void UCfInterpolation(Foam::volVectorField& U, double ReTau, double deltaChannel, volScalarField& nu, Foam::fvMesh& mesh, int cwipiObsU, int cwipiVerbose, std::string UpIntPath);

void cwipiSend(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiVerbose);

void cwipiSendParams(const fvMesh& mesh, const volVectorField& vf, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

void cwipiSendParamsChannel(const fvMesh& mesh, const volScalarField& Ck, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

void cwipiSendParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

void cwipiSendParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, const Time& runTime, int cwipiIteration, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

void cwipiRecv(const fvMesh& mesh, volVectorField& U, const Time& runTime, int cwipiIteration, int cwipiVerbose);

void cwipiRecvParams(const fvMesh& mesh, volVectorField& U, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

void cwipiRecvParamsChannel(const fvMesh& mesh, volScalarField& Ck, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose);

std::fstream& GotoLine(std::fstream& file, unsigned int num);

void cwipiSendParamsChannel_IBM(const fvMesh& mesh, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsChannel_IBM(const fvMesh& mesh, volVectorField& U, volScalarField& viscosity, volVectorField& F, volVectorField& D, int cwipiParams, int Fourier_m, int Fourier_n, int Fourier_o, 
double lambdaX, double lambdaY, double lambdaZ, double deltaChannel, int DarcyForchheimer, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsChannelNOFourier(const fvMesh& mesh, volVectorField& U, volScalarField& viscosity, volScalarField& nu, volVectorField& F, volVectorField& D, int cwipiParams, double deltaChannel, int DarcyForchheimer, 
int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsKEps(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath);

void cwipiRecvParamsKOmegaSST(const fvMesh& mesh, incompressible::momentumTransportModel& turbulence, int cwipiParams, int nbParts, int partsReparty, int cwipiVerbose, std::string globalRootPath);

void cwipideleteCoupling(double* pointCoords, int* face_index, int* face_connectivity_index, int* cell_to_face_connectivity, int* face_connectivity, int cwipiVerbose);

}

// ************************************************************************* //
