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
    Hrvoje Jasak, Wikki Ltd.

Note on parallelisation
    In order to handle parallelisation correctly, I need to rely on the fact
    that all patches that require a global gather-scatter come before
    processor patches.  In that case, the communication pattern
    will be correct without intervention.  HJ, 6/Aug/2009

\*---------------------------------------------------------------------------*/

#include "jumpRegionCouplingFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
jumpRegionCouplingFvPatchField<Type>::jumpRegionCouplingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    regionCouplingFvPatchField<Type>(p, iF)
{}


template<class Type>
jumpRegionCouplingFvPatchField<Type>::jumpRegionCouplingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
jumpRegionCouplingFvPatchField<Type>::jumpRegionCouplingFvPatchField
(
    const jumpRegionCouplingFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
jumpRegionCouplingFvPatchField<Type>::jumpRegionCouplingFvPatchField
(
    const jumpRegionCouplingFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    regionCouplingFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> >
jumpRegionCouplingFvPatchField<Type>::patchNeighbourField() const
{
    tmp<Field<Type> > tpnf
    (
        this->regionCouplePatch().interpolate
        (
            this->shadowPatchField().patchInternalField()
        )
    );

    Field<Type>& pnf = tpnf.ref();

    if (this->regionCouplePatch().bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        // Use mirrored neighbour field for interpolation
        // HJ, 21/Jan/2009
        Field<Type> mirrorField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

        // Set mirror values to fully uncovered faces
        this->regionCouplePatch().setUncoveredFaces(mirrorField, pnf);

        // For partially covered faces, add mirror that causes no flux
        this->regionCouplePatch().addToPartialFaces(mirrorField, pnf);
    }

    // Add jump
    pnf += jump();

    return tpnf;
}


template<class Type>
void jumpRegionCouplingFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (this->regionCouplePatch().coupled())
    {
        // Prepare local matrix update buffer for the remote side.
        // Note that only remote side will have access to its psiInternal
        // as they are on different regions

        // Since interpolation needs to happen on the shadow, and within the
        // init, prepare interpolation for the other side.
        this->setMatrixUpdateBuffer
        (
            this->shadowPatchField().regionCouplePatch().interpolate
            (
                this->patch().patchInternalField(psiInternal)
            )
            // Add jump
          + jump()
        );
    }
    else
    {
        FatalErrorInFunction
            << "init matrix update called in detached state"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
