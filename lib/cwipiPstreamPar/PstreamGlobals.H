/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for mpi

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/
/**
 * @file PstreamGlobals.H
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Global functions and variables for working with parallel streams, but principally for mpi. 
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef PstreamGlobals_H
#define PstreamGlobals_H

#include "DynamicList.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PstreamGlobals Declaration
\*---------------------------------------------------------------------------*/

namespace PstreamGlobals
{
    extern MPI_Comm MPI_COMM_FOAM;

    extern DynamicList<MPI_Request> outstandingRequests_;

    extern int nTags_;

    extern DynamicList<int> freedTags_;

    // Current communicators. First element will be MPI_COMM_FOAM
    extern DynamicList<MPI_Comm> MPICommunicators_;

    extern DynamicList<MPI_Group> MPIGroups_;

    void checkCommunicator(const label, const label procNo);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
