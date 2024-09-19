/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "DASmagorinsky.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> DASmagorinsky<BasicMomentumTransportModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return volScalarField::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
    );
}


template<class BasicMomentumTransportModel>
void DASmagorinsky<BasicMomentumTransportModel>::correctNut()
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
DASmagorinsky<BasicMomentumTransportModel>::DASmagorinsky
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    Ck_
	(
		IOobject
		(
			"Ck",
			this->runTime_.timeName(),
			this->mesh_,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh_
	)

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool DASmagorinsky<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        Foam::Info << "Ck field updated in the turbulence " << Foam::endl;

        volScalarField Cktemp
        (
            IOobject
            (
                "Ck",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_
        );

        Ck_ = Cktemp;

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> DASmagorinsky<BasicMomentumTransportModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k*sqrt(k)/this->delta()
    );
}


template<class BasicMomentumTransportModel>
void DASmagorinsky<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
