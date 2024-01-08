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

#include "cyclicFvPatch.H"
#include "fvMesh.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "slicedSurfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Make mesh cell centres.  Moved from fvMeshGeometry
void Foam::cyclicFvPatch::makeC(slicedSurfaceVectorField& C) const
{
    C.boundaryField()[index()].UList<vector>::shallowCopy
    (
        patchSlice(cyclicPolyPatch_.boundaryMesh().mesh().faceCentres())
    );
}


// Make patch weighting factors
void Foam::cyclicFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    const scalarField& magFa = magSf();

    // Note: mag in the dot-product.
    // For all valid meshes, the non-orthogonality will be less than
    // 90 deg and the dot-product will be positive.  For invalid
    // meshes (d & s <= 0), this will stabilise the calculation
    // but the result will be poor.  HJ, 24/Aug/2011
    scalarField deltas = mag(nf() & fvPatch::delta());
    label sizeby2 = deltas.size()/2;

    scalar maxMatchError = 0;
    label errorFace = -1;

    for (label faceI = 0; faceI < sizeby2; faceI++)
    {
        scalar avFa = (magFa[faceI] + magFa[faceI + sizeby2])/2.0;

        if
        (
            mag(magFa[faceI] - magFa[faceI + sizeby2])/avFa
          > polyPatch::matchTol_()
        )
        {
            // Found error.  Look for largest matching error
            maxMatchError =
                Foam::max
                (
                    maxMatchError,
                    mag(magFa[faceI] - magFa[faceI + sizeby2])/avFa
                );

            errorFace = faceI;
        }

        scalar di = deltas[faceI];
        scalar dni = deltas[faceI + sizeby2];

        w[faceI] = dni/(di + dni);
        w[faceI + sizeby2] = 1 - w[faceI];
    }

    // Check for error in matching
    if (maxMatchError > polyPatch::matchTol_())
    {
        scalar avFa = (magFa[errorFace] + magFa[errorFace + sizeby2])/2.0;

        FatalErrorIn("cyclicFvPatch::makeWeights(fvsPatchScalarField& w) const")
            << "face " << errorFace << " and " << errorFace + sizeby2
            <<  " areas do not match by "
            << 100*mag(magFa[errorFace] - magFa[errorFace + sizeby2])/avFa
            << "% -- possible face ordering problem." << nl
            << "Cyclic area match tolerance = "
            << polyPatch::matchTol_() << " patch: " << name()
            << abort(FatalError);
    }
}


// Make patch face - neighbour cell distances
void Foam::cyclicFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    vectorField d = delta();
    vectorField n = nf();
    label sizeby2 = d.size()/2;

    for (label faceI = 0; faceI < sizeby2; faceI++)
    {
        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        dc[faceI] = 1.0/max(n[faceI] & d[faceI], 0.05*mag(d[faceI]));
        dc[faceI + sizeby2] = dc[faceI];
    }
}


// Make patch face - neighbour cell distances
void Foam::cyclicFvPatch::makeMagLongDeltas(fvsPatchScalarField& mld) const
{
    vectorField d = fvPatch::delta();
    vectorField patchDelta = this->delta();
    vectorField pSf = Sf();
    scalarField pMagSf = magSf();
    label sizeby2 = d.size()/2;

    for (label faceI = 0; faceI < sizeby2; faceI++)
    {
        // NOT stabilised for bad meshes.  HJ, 11/May/2020
        mld[faceI] =
        (
            mag(pSf[faceI] & d[faceI])
          + mag(pSf[faceI] & (patchDelta[faceI] - d[faceI]))
        )/pMagSf[faceI];

        mld[faceI + sizeby2] = mld[faceI];
    }
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::cyclicFvPatch::delta() const
{
    vectorField patchD = fvPatch::delta();
    label sizeby2 = patchD.size()/2;

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv.ref();

    // To the transformation if necessary
    if (parallel())
    {
        for (label faceI = 0; faceI < sizeby2; faceI++)
        {
            vector ddi = patchD[faceI];
            vector dni = patchD[faceI + sizeby2];

            pdv[faceI] = ddi - dni;
            pdv[faceI + sizeby2] = -pdv[faceI];
        }
    }
    else
    {
        for (label faceI = 0; faceI < sizeby2; faceI++)
        {
            vector ddi = patchD[faceI];
            vector dni = patchD[faceI + sizeby2];

            pdv[faceI] = ddi - transform(forwardT()[0], dni);
            pdv[faceI + sizeby2] = -transform(reverseT()[0], pdv[faceI]);
        }
    }

    return tpdv;
}


Foam::tmp<Foam::labelField> Foam::cyclicFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicFvPatch::transfer
(
    const Pstream::commsTypes,
    const labelUList& interfaceData
) const
{
    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf.ref();

    label sizeby2 = this->size()/2;

    for (label faceI = 0; faceI < sizeby2; faceI++)
    {
        pnf[faceI] = interfaceData[faceI + sizeby2];
        pnf[faceI + sizeby2] = interfaceData[faceI];
    }

    return tpnf;
}


Foam::tmp<Foam::labelField> Foam::cyclicFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    const labelUList& faceCells = this->patch().faceCells();

    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf.ref();

    label sizeby2 = this->size()/2;

    for (label faceI = 0; faceI < sizeby2; faceI++)
    {
        pnf[faceI] = iF[faceCells[faceI + sizeby2]];
        pnf[faceI + sizeby2] = iF[faceCells[faceI]];
    }

    return tpnf;
}


// ************************************************************************* //
