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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvPatch.H"
#include "fvMesh.H"
#include "fvPatchFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "slicedSurfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGgiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicGgiFvPatch, polyPatch);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::cyclicGgiFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (cyclicGgiPolyPatch_.master())
    {
        vectorField n = nf();

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        scalarField nfc =
            mag
            (
                n & (cyclicGgiPolyPatch_.reconFaceCellCentres() - Cf())
            );

        w = nfc/(mag(n & (Cf() - Cn())) + nfc + SMALL);

        // Master bridging is not needed, as reconFaceCellCentres
        // has already been bridged.  HJ, 12/Dec/2022
    }
    else
    {
        // Slave side. Interpolate the master side weights, scale them for
        // partially covered faces and set weights for fully uncovered faces if
        // the bridge overlap is switched on. VV, 15/Feb/2018.

        // Pick up weights from the master side
        fvsPatchScalarField masterWeights
        (
            shadow(),
            w.dimensionedInternalField()
        );

        shadow().makeWeights(masterWeights);

        // Interpolate master weights to this side
        w = interpolate(masterWeights);

        // Slave bridging is necesary.  HJ, 12/Dec/2022
        if (bridgeOverlap())
        {
            // Weights for fully uncovered faces
            const scalarField uncoveredWeights(w.size(), 0.5);

            // Set weights for uncovered faces
            setUncoveredFaces(uncoveredWeights, w);

            // Cannot manipulate partially covered faces, as this breaks
            // symmetry and causes conservation errors.
            // HJ, 12/Dec/2022
        }

        // Finally construct these weights as 1 - master weights
        w = 1 - w;
    }
}


// Make patch face - neighbour cell distances
void Foam::cyclicGgiFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    if (cyclicGgiPolyPatch_.master())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.

        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        const vectorField d = delta();

        dc = 1.0/max(nf() & d, 0.05*mag(d));

        // Note: no need to bridge the overlap since delta already takes it into
        // account. VV, 18/Oct/2017.
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltaCoeffs for fully uncovered faces if the
        // bridge overlap is switched on. VV, 15/Feb/2018.

        fvsPatchScalarField masterDeltas
        (
            shadow(),
            dc.dimensionedInternalField()
        );

        shadow().makeDeltaCoeffs(masterDeltas);

        dc = interpolate(masterDeltas);

        // For bridging, we use a mirror: weights are 0.5, not 1
        // 12/Nov/2022
        if (bridgeOverlap())
        {
            // Bugfix: delta coeffs for double distance are halved
            // Function fvPatch::deltaCoeffs() cannot be called here:
            // calculation is not complete and it returns zero.
            // HJ, 4/Dec/2022
            const scalarField uncoveredDeltaCoeffs =
                1/(2*mag(fvPatch::delta()));

            // Set delta coeffs for uncovered faces
            setUncoveredFaces(uncoveredDeltaCoeffs, dc);

            // Cannot manipulate partially covered faces, as this breaks
            // symmetry and causes conservation errors.
            // HJ, 12/Dec/2022
        }
    }
}

// Make patch face long distance factors
void Foam::cyclicGgiFvPatch::makeMagLongDeltas(fvsPatchScalarField& mld) const
{
    if (cyclicGgiPolyPatch_.master())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.

        vectorField d = fvPatch::delta();

        // NOT stabilised for bad meshes.  HJ, 11/May/2020
        mld = (mag(Sf() & d) + mag(Sf() & (delta() - d)))/magSf();

        // Note: no need to bridge the overlap since delta already takes it into
        // account. VV, 18/Oct/2017.  Correct.  HJ, 12/Dec/2022
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltaCoeffs for fully uncovered faces if the
        // bridge overlap is switched on. VV, 15/Feb/2018.

        fvsPatchScalarField masterLongDeltas
        (
            shadow(),
            mld.dimensionedInternalField()
        );

        shadow().makeMagLongDeltas(masterLongDeltas);

        mld = interpolate(masterLongDeltas);

        if (bridgeOverlap())
        {
            // Bugfix: delta coeffs for double distance are halved
            // HJ, 4/Dec/2022
            const scalarField uncoveredMagDeltas =
                1/(2*mag(fvPatch::delta()));

            // Set delta coeffs for uncovered faces
            setUncoveredFaces(uncoveredMagDeltas, mld);

            // Cannot manipulate partially covered faces, as this breaks
            // symmetry and causes conservation errors.
            // HJ, 12/Dec/2022
        }
    }
}


// Make patch face non-orthogonality correction vectors
void Foam::cyclicGgiFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
{
    // Non-orthogonality correction on a ggi interface
    // MB, 7/April/2009

    // No non-orthogonal correction if the bridge overlap is switched on to
    // ensure conservative interpolation for partially overlapping faces.
    // VV??  This is wrong.  HJ, 2/Dec/2022

    // Calculate correction vectors on coupled patches
    const scalarField& patchDeltaCoeffs = deltaCoeffs();

    const vectorField patchDeltas = delta();
    const vectorField n = nf();

    // If non-orthogonality is over 90 deg, kill correction vector
    // HJ, 6/Jan/2011
    cv = pos(patchDeltas & n)*(n - patchDeltas*patchDeltaCoeffs);

    if (bridgeOverlap())
    {
        // On uncovered faces, corr vectors are zero
        const vectorField bridgeCorrVecs(size(), vector::zero);

        // Set delta coeffs for uncovered faces
        setUncoveredFaces(bridgeCorrVecs, cv);

        // Do nothing to partially overlapping faces: correction for the
        // uncovered part is zero
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::cyclicGgiFvPatch::delta() const
{
    if (cyclicGgiPolyPatch_.master())
    {
        // Master side. Note: scaling partially covered faces and setting deltas
        // to fully uncovered faces correctly taken into account in
        // reconFaceCellCentres function. VV, 15/Feb/2018.

        tmp<vectorField> tdelta =
            cyclicGgiPolyPatch_.reconFaceCellCentres() - Cn();

        // Master bridging is not needed, as reconFaceCellCentres
        // has already been bridged.  HJ, 12/Dec/2022

        return tdelta;
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltas for fully uncovered faces if the bridge
        // overlap is switched on. VV, 15/Feb/2018.

        tmp<vectorField> tdelta = interpolate
        (
            shadow().Cn() - cyclicGgiPolyPatch_.shadow().reconFaceCellCentres()
        );
        vectorField& delta = tdelta.ref();

        if (bridgeOverlap())
        {
            // Deltas for fully uncovered faces
            const vectorField uncoveredDeltas = 2*fvPatch::delta();

            // Set deltas for fully uncovered faces
            setUncoveredFaces(uncoveredDeltas, delta);

            // Cannot manipulate partially covered faces, as this breaks
            // symmetry and causes conservation errors.
            // HJ, 12/Dec/2022
        }

        return tdelta;
    }
}



const Foam::cyclicGgiFvPatch& Foam::cyclicGgiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[cyclicGgiPolyPatch_.shadowIndex()];
    return refCast<const cyclicGgiFvPatch>(p);
}


// ************************************************************************* //
