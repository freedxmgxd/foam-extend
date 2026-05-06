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

#include "cellBoundingBoxDiagDiff.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(cellBoundingBoxDiagDiff, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    cellBoundingBoxDiagDiff,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::cellBoundingBoxDiagDiff::cellBoundingBoxDiagDiff
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict),
    threshold_(readScalar(coeffDict().lookup("threshold"))),
    cellBBDiag_(oversetFringeAlgorithm.mesh().nCells())
{
    // Sanity check
    if (threshold_ < SMALL)
    {
        FatalIOErrorInFunction(coeffDict())
            << "Zero or negative threshold specified. This is not allowed"
            << abort(FatalIOError);
    }

    // Get reference to fvMesh
    const fvMesh& mesh = oversetFringeAlgorithm.mesh();

    // Get necessary mesh data
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    // Create local donor suitability function
    scalarField localDsf(mesh.nCells(), 0);

    // Loop through cells and calculate the bounding box diagonal of each cell
    forAll (cells, cellI)
    {
        const boundBox bb(cells[cellI].points(faces, points), false);
        cellBBDiag_[cellI] = bb.mag();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::donorSuitability::cellBoundingBoxDiagDiff::value
(
    const label& cellID
) const
{
    return cellBBDiag_[cellID];
}


Foam::scalar 
Foam::donorSuitability::cellBoundingBoxDiagDiff::suitabilityFraction
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
        // Return relative difference in donor and acceptor value

        const scalar dsfAcceptor = daPair.acceptorSuitability();

        const scalar dsfDonor = daPair.donorSuitability();

        // Calculate suitability from difference between donor and acceptor
        // min face area
        return
        (
            1 - mag(dsfAcceptor - dsfDonor)/
            (Foam::max(dsfAcceptor, dsfDonor) + SMALL)
        );
    }
}


bool Foam::donorSuitability::cellBoundingBoxDiagDiff::isDonorSuitable
(
    const donorAcceptor& daPair
) const
{
    return (suitabilityFraction(daPair)) > threshold();
}


// ************************************************************************* //
