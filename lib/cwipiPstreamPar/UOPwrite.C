/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/
/**
 * @file UOPwrite.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Write primitive and binary block from OPstream.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "UOPstream.H"
#include "PstreamGlobals.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPstream::write
(
    const commsTypes commsType,
    const int toProcNo,
    const char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator
)
{
    if (debug)
    {
        Pout<< "UOPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << communicator << " size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "UOPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << communicator << " size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }


    PstreamGlobals::checkCommunicator(communicator, toProcNo);


    bool transferFailed = true;

    if (commsType == commsTypes::blocking)
    {
        transferFailed = MPI_Bsend
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   // procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[communicator]
        );

        if (debug)
        {
            Pout<< "UOPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == commsTypes::scheduled)
    {
        transferFailed = MPI_Send
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   // procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[communicator]
        );

        if (debug)
        {
            Pout<< "UOPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == commsTypes::nonBlocking)
    {
        MPI_Request request;

        transferFailed = MPI_Isend
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   // procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
            &request
        );

        if (debug)
        {
            Pout<< "UOPstream::write : started write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << Foam::abort(FatalError);
    }

    return !transferFailed;
}


// ************************************************************************* //
