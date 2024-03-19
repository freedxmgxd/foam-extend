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

Description
    Region couple patch

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCoupleFvPatch.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "fvBoundaryMesh.H"
#include "mapDistribute.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupleFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupleFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::regionCoupleFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            vectorField n = nf();

            // Note: mag in the dot-product.
            // For all valid meshes, the non-orthogonality will be less than
            // 90 deg and the dot-product will be positive.  For invalid
            // meshes (d & s <= 0), this will stabilise the calculation
            // but the result will be poor.  HJ, 24/Aug/2011
            scalarField nfc =
                mag(n & (rcPolyPatch_.reconFaceCellCentres() - Cf()));

            w = nfc/(mag(n & (Cf() - Cn())) + nfc);

            // Master bridging is not needed, as reconFaceCellCentres
            // has already been bridged.  HJ, 12/Dec/2022

            // For partial faces, recon cell centre is adjusted for live part
            // HJ, 3/Jan/2022
        }
        else
        {
            // Pick up weights from the master side
            fvsPatchScalarField masterWeights
            (
                shadow(),
                w.dimensionedInternalField()
            );

            shadow().makeWeights(masterWeights);

            w = interpolate(1 - masterWeights);

            if (bridgeOverlap())
            {
                // Weights for fully uncovered faces
                const scalarField uncoveredWeights(w.size(), 0.5);

                // Set weights for uncovered faces
                setUncoveredFaces(uncoveredWeights, w);

                // Since weights are not scaled in partial cover faces, add
                addToPartialFaces(uncoveredWeights, w);
            }
        }
    }
    else
    {
        fvPatch::makeWeights(w);
    }
}


// Make patch face delta coefficients
void Foam::regionCoupleFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            // Stabilised form for bad meshes.  HJ, 24/Aug/2011
            vectorField d = delta();

            dc = 1.0/max(nf() & d, 0.05*mag(d));

            // Note: no need to bridge the overlap since delta already takes
            //  it into account. VV, 18/Oct/2017.  Correct

            // Partially covered faces will get delta from the live part
        }
        else
        {
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

                // For partial faces, scale delta to account for only the
                // covered part.  HJ, 3/Jan/2022
                scalePartialFaces(dc);
            }
        }
    }
    else
    {
        fvPatch::makeDeltaCoeffs(dc);
    }
}


// Make patch face long deltas
void Foam::regionCoupleFvPatch::makeMagLongDeltas
(
    fvsPatchScalarField& mld
) const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            // NOT stabilised for bad meshes.  HJ, 11/May/2020
            vectorField d = fvPatch::delta();

            mld = (mag(Sf() & d) + mag(Sf() & (delta() - d)))/magSf();

            // Note: no need to bridge the overlap since delta already takes
            // it into account. VV, 18/Oct/2017.  Correct.  HJ, 12/Dec/2022
        }
        else
        {
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

                // For partial faces, scale delta to account for only the
                // covered part.  HJ, 3/Jan/2022
                scalePartialFaces(mld);
            }
        }
    }
    else
    {
        fvPatch::makeMagLongDeltas(mld);
    }
}


// Make patch face non-orthogonality correction vectors
void Foam::regionCoupleFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
{
    if (rcPolyPatch_.coupled())
    {
        // Non-orthogonality correction in attached state identical to ggi
        // interface

        // Calculate correction vectors on coupled patches
        const scalarField& patchDeltaCoeffs = deltaCoeffs();

        vectorField patchDeltas = delta();
        vectorField n = nf();

        // If non-orthogonality is over 90 deg, kill correction vector
        // HJ, 6/Jan/2011
        cv = pos(patchDeltas & n)*(n - patchDeltas*patchDeltaCoeffs);

        if (bridgeOverlap())
        {
            // On uncovered faces, corr vectors are zero
            const vectorField bridgeCorrVecs(size(), vector::zero);

            // Set delta coeffs for uncovered faces
            setUncoveredFaces(bridgeCorrVecs, cv);

            // Kill correction on partially overlapping faces for consistency:
            // Non-orthogonal correction is not allowed on walls and symm planes
            setPartialFaces(bridgeCorrVecs, cv);
        }
    }
    else
    {
        // No correction in detached state.  HJ, 26/Jul/2011
        cv = vector::zero;
    }
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::regionCoupleFvPatch::delta() const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            tmp<vectorField> tDelta =
                rcPolyPatch_.reconFaceCellCentres() - Cn();

            // Master bridging is not needed, as reconFaceCellCentres
            // has already been bridged.  HJ, 12/Dec/2022

            return tDelta;
        }
        else
        {
            tmp<vectorField> tdelta = interpolate
            (
                shadow().Cn() - rcPolyPatch_.shadow().reconFaceCellCentres()
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
    else
    {
        return fvPatch::delta();
    }
}


bool Foam::regionCoupleFvPatch::coupled() const
{
    return rcPolyPatch_.coupled();
}


const Foam::fvMesh& Foam::regionCoupleFvPatch::shadowRegion() const
{
    return
        boundaryMesh().mesh().objectRegistry::parent().lookupObject<fvMesh>
        (
            rcPolyPatch_.shadowRegionName()
        );
}


const Foam::regionCoupleFvPatch& Foam::regionCoupleFvPatch::shadow() const
{
    const fvPatch& p =
        shadowRegion().boundary()[rcPolyPatch_.shadowIndex()];

    return refCast<const regionCoupleFvPatch>(p);
}


bool Foam::regionCoupleFvPatch::master() const
{
    return rcPolyPatch_.master();
}


bool Foam::regionCoupleFvPatch::fineLevel() const
{
    return true;
}


Foam::label Foam::regionCoupleFvPatch::shadowIndex() const
{
    return rcPolyPatch_.shadowIndex();
}


const Foam::ggiLduInterface&
Foam::regionCoupleFvPatch::shadowInterface() const
{
    const fvPatch& p =
        shadowRegion().boundary()[rcPolyPatch_.shadowIndex()];

    return refCast<const ggiLduInterface>(p);
}


Foam::label Foam::regionCoupleFvPatch::interfaceSize() const
{
    return rcPolyPatch_.size();
}


Foam::label Foam::regionCoupleFvPatch::zoneSize() const
{
    return rcPolyPatch_.zone().size();
}


const Foam::labelList& Foam::regionCoupleFvPatch::zoneAddressing() const
{
    return rcPolyPatch_.zoneAddressing();
}


const Foam::labelListList& Foam::regionCoupleFvPatch::ggiAddressing() const
{
    if (rcPolyPatch_.master())
    {
        return rcPolyPatch_.patchToPatch().masterAddr();
    }
    else
    {
        return rcPolyPatch_.patchToPatch().slaveAddr();
    }
}


bool Foam::regionCoupleFvPatch::localParallel() const
{
    return rcPolyPatch_.localParallel();
}


const Foam::mapDistribute& Foam::regionCoupleFvPatch::map() const
{
    return rcPolyPatch_.map();
}


const Foam::scalarListList& Foam::regionCoupleFvPatch::ggiWeights() const
{
    if (rcPolyPatch_.master())
    {
        return rcPolyPatch_.patchToPatch().masterWeights();
    }
    else
    {
        return rcPolyPatch_.patchToPatch().slaveWeights();
    }
}


void Foam::regionCoupleFvPatch::expandAddrToZone(labelField& lf) const
{
    // Missing code.  Activate for AMG solvers across regionCoupleFvPatch
    notImplemented
    (
        "void regionCoupleFvPatch::expandAddrToZone(labelField& lf) const"
    );
}


void Foam::regionCoupleFvPatch::expandCrMatrixToZone(crMatrix&) const
{
    // Missing code.  Activate for SAMG solvers across regionCoupleFvPatch
    notImplemented
    (
        "void regionCoupleFvPatch::expandCrMatrixToZone(crMatrix&) const"
    );
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::regionCoupleFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{

    return shadow().labelTransferBuffer();
}


void Foam::regionCoupleFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return shadow().labelTransferBuffer();
}


// ************************************************************************* //
