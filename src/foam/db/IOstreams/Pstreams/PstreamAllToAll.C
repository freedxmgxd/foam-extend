/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamGlobals.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Pstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        if
        (
            MPI_Alltoall
            (
                const_cast<label*>(sendData.begin()),
                sizeof(label),
                MPI_BYTE,
                recvData.begin(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


void Foam::Pstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {
        if
        (
            MPI_Alltoallv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


// ************************************************************************* //
