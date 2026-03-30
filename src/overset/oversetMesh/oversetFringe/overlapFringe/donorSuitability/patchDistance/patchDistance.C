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

#include "patchDistance.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "patchWave.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(patchDistance, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    patchDistance,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::patchDistance::patchDistance
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict),
    magDistance_()
{
    // Get reference to fvMesh
    const fvMesh& mesh = oversetFringeAlgorithm.mesh();

    // Get distance patch names for master and donor regions
    wordList patchNames =
        coeffDict().lookup("distancePatches");

    // Insert patch IDs into hash sets
    labelHashSet masterPatchIDs(patchNames.size());

    forAll (patchNames, patchI)
    {
        polyPatchID pID(patchNames[patchI], mesh.boundaryMesh());

        if (pID.active())
        {
            masterPatchIDs.insert(pID.index());
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find distance patch named: "
                << patchNames[patchI]
                << "Available patch names: " << mesh.boundaryMesh().names()
                << abort(FatalError);
        }
    }


    // Calculate distance from specified patches and do not correct for accurate
    // near wall distance (=false paramater)
    patchWave masterDistance(mesh, masterPatchIDs, false);

    // Combine both master and donor distances into a single field
    magDistance_ = masterDistance.distance();

    magDistance_ /= gMax(magDistance_);

    Info<< "magDistance_ : " << min(magDistance_) << " " << max(magDistance_ ) << endl;
}
    

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::donorSuitability::patchDistance::value
(
    const label& cellID
) const
{
    return 1;
}


Foam::scalar Foam::donorSuitability::patchDistance::suitabilityFraction
(
    const donorAcceptor& daPair
) const
{
    // Check whether the donor is valid for this pair
    if (!daPair.donorFound())
    {
        // No donor: return zero
        return 0;
    }
    else
    {
        if (daPair.donorProcNo() != Pstream::myProcNo())
        {
            FatalErrorInFunction
                << "Donor on different processor: this cannot happen: "
                << "myProc = " << Pstream::myProcNo()
                << " donorProc = " << daPair.donorProcNo()
                << abort(FatalError);
        }

        // Get local donor suitability function
        return magDistance_[daPair.donorCell()];
    }
}


bool Foam::donorSuitability::patchDistance::isDonorSuitable
(
    const donorAcceptor& daPair
) const
{
    // Suitable if donor found
    return (suitabilityFraction(daPair)) > SMALL;
}

    
// ************************************************************************* //
