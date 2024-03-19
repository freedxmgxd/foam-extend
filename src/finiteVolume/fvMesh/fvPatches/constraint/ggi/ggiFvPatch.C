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
    Generalized grid interface (GGI) patch, providing coupling
    between arbitrary patches which belong to the same fvMesh

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ggiFvPatch.H"
#include "fvPatchFields.H"
#include "fvBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, ggiFvPatch, polyPatch);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::ggiFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (ggiPolyPatch_.master())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.
        const vectorField n = nf();

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        const scalarField nfc =
            mag(n & (ggiPolyPatch_.reconFaceCellCentres() - Cf()));

        w = nfc/(mag(n & (Cf() - Cn())) + nfc);

        // Master bridging is not needed, as reconFaceCellCentres
        // has already been bridged.  HJ, 12/Dec/2022

        // For partial faces, recon cell centre is adjusted for the live part
        // HJ, 3/Jan/2022
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
        w = interpolate(1 - masterWeights);

        // Slave bridging is necesary.  HJ, 12/Dec/2022
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


// Make patch face delta coefficients
void Foam::ggiFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    if (ggiPolyPatch_.master())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.

        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        const vectorField d = delta();

        dc = 1.0/max(nf() & d, 0.05*mag(d));

        // Note: no need to bridge the overlap since delta already takes it into
        // account. VV, 18/Oct/2017.  Correct

        // Partially covered faces will get delta from the live part
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

            // For partial faces, scale delta to account for only the
            // covered part.  HJ, 3/Jan/2022
            scalePartialFaces(dc);
        }
    }
}


// Make patch face long distance factors
void Foam::ggiFvPatch::makeMagLongDeltas(fvsPatchScalarField& mld) const
{
    if (ggiPolyPatch_.master())
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

            // For partial faces, scale delta to account for only the
            // covered part.  HJ, 3/Jan/2022
            scalePartialFaces(mld);
        }
    }
}


// Make patch face non-orthogonality correction vectors
void Foam::ggiFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
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

        // Kill correction on partially overlapping faces for consistency:
        // Non-orthogonal correction is not allowed on walls and symm planes
        setPartialFaces(bridgeCorrVecs, cv);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::ggiFvPatch::delta() const
{
    if (ggiPolyPatch_.master())
    {
        // Master side. Note: scaling partially covered faces and setting deltas
        // to fully uncovered faces correctly taken into account in
        // reconFaceCellCentres function. VV, 15/Feb/2018.

        tmp<vectorField> tdelta = ggiPolyPatch_.reconFaceCellCentres() - Cn();

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
            shadow().Cn() - ggiPolyPatch_.shadow().reconFaceCellCentres()
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


const Foam::ggiFvPatch& Foam::ggiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiFvPatch>(p);
}


bool Foam::ggiFvPatch::master() const
{
    return ggiPolyPatch_.master();
}


bool Foam::ggiFvPatch::fineLevel() const
{
    return true;
}


Foam::label Foam::ggiFvPatch::shadowIndex() const
{
    return ggiPolyPatch_.shadowIndex();
}


const Foam::ggiLduInterface& Foam::ggiFvPatch::shadowInterface() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiLduInterface>(p);
}


Foam::label Foam::ggiFvPatch::interfaceSize() const
{
    return ggiPolyPatch_.size();
}


Foam::label Foam::ggiFvPatch::zoneSize() const
{
    return ggiPolyPatch_.zone().size();
}


const Foam::labelList& Foam::ggiFvPatch::zoneAddressing() const
{
    return ggiPolyPatch_.zoneAddressing();
}


const Foam::labelListList& Foam::ggiFvPatch::ggiAddressing() const
{
    if (ggiPolyPatch_.master())
    {
        return ggiPolyPatch_.patchToPatch().masterAddr();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveAddr();
    }
}


bool Foam::ggiFvPatch::localParallel() const
{
    return ggiPolyPatch_.localParallel();
}


const Foam::mapDistribute& Foam::ggiFvPatch::map() const
{
    return ggiPolyPatch_.map();
}


const Foam::scalarListList& Foam::ggiFvPatch::ggiWeights() const
{
    if (ggiPolyPatch_.master())
    {
        return ggiPolyPatch_.patchToPatch().masterWeights();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveWeights();
    }
}


void Foam::ggiFvPatch::expandAddrToZone(labelField& lf) const
{
    lf = ggiPolyPatch_.fastExpand(lf);
}


void Foam::ggiFvPatch::expandCrMatrixToZone(crMatrix& patchP) const
{
    if (!localParallel())
    {
        // Split the crMatrix into rows and expand it
        const crAddressing& patchCrAddr = patchP.crAddr();
        const labelList& patchRowStart = patchCrAddr.rowStart();
        const labelList& patchCol = patchCrAddr.column();
        const scalarField& patchCoeff = patchP.coeffs();

        List<labelField> cols(patchCrAddr.nRows());
        List<scalarField> coeffs(patchCrAddr.nRows());

        for (label faceI = 0; faceI < patchCrAddr.nRows(); faceI++)
        {
            // Unpack row
            const label rowStart = patchRowStart[faceI];
            const label rowLength = patchRowStart[faceI + 1] - rowStart;

            cols[faceI].setSize(rowLength);
            labelField& curCols = cols[faceI];

            coeffs[faceI].setSize(rowLength);
            scalarField& curCoeffs = coeffs[faceI];

            for (label coeffI = 0; coeffI < rowLength; coeffI++)
            {
                curCols[coeffI] = patchCol[rowStart + coeffI];
                curCoeffs[coeffI] = patchCoeff[rowStart + coeffI];
            }
        }

        // Expand to zone size
        List<labelField> zoneColsFF = ggiPolyPatch_.fastExpand(cols);
        List<scalarField> zoneCoeffsFF = ggiPolyPatch_.fastExpand(coeffs);

        scalar nZoneEntries = 0;

        forAll (zoneColsFF, zfI)
        {
            nZoneEntries += zoneColsFF[zfI].size();
        }

        // Reconstruct matrix
        labelList zoneRowStart(zoneSize() + 1);
        labelList zoneCols(nZoneEntries);
        scalarField zoneCoeffs(nZoneEntries);

        zoneRowStart[0] = 0;

        // Reset nZoneEntries for use as a counter
        nZoneEntries = 0;

        forAll(zoneColsFF, zfI)
        {
            const labelField& curCols = zoneColsFF[zfI];
            const scalarField& corCoeffs = zoneCoeffsFF[zfI];

            zoneRowStart[zfI + 1] = zoneRowStart[zfI] + curCols.size();

            forAll (curCols, coeffI)
            {
                zoneCols[nZoneEntries] = curCols[coeffI];
                zoneCoeffs[nZoneEntries] = corCoeffs[coeffI];
                nZoneEntries++;
            }
        }
        patchP = crMatrix
        (
            zoneSize(),
            patchCrAddr.nCols(),
            zoneRowStart,
            zoneCols
        );

        // Set coeffs
        patchP.coeffs() = zoneCoeffs;
    }
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::ggiFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::transfer
(
    const Pstream::commsTypes,
    const labelUList& interfaceData
) const
{
    return this->shadow().labelTransferBuffer();
}


void Foam::ggiFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    // Label transfer is local without global reduction
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    return shadow().labelTransferBuffer();
}


void Foam::ggiFvPatch::initProlongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    // crMatrix transfer is local without global reduction
    crMatrixTransferBuffer_ = filteredP;
}


Foam::autoPtr<Foam::crMatrix> Foam::ggiFvPatch::prolongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    autoPtr<crMatrix> tnbrP(new crMatrix(shadow().crMatrixTransferBuffer()));

    return tnbrP;
}


// ************************************************************************* //
