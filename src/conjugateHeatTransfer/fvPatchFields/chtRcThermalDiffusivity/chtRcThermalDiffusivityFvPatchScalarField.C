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
    Code removed in entirety: Henrik Rusche, Wikki GmbH.
    This was a disastrous heap of junk
    Written anew by Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcThermalDiffusivityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::chtRcThermalDiffusivityFvPatchScalarField::k() const
{
    // For thermal diffusivity (solving T equation), return decoupled k
    return decoupledField();
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // In order to deal with wall functions, decoupled patch field values
    // need to be interpolated
    // HJ, 7/Jun/2024

    if (!regionCouplePatch().coupled())
    {
        // Store coupled field
        decoupledField();
    }
    else
    {
        // Do interpolation
    
        // Get local decoupled patch field values
        const scalarField& diffOwn = decoupledField();

        // Interpolate neighbour decoupled patch field values

        // Get neighbour k
        const chtRcThermalDiffusivityFvPatchScalarField& pCht =
            refCast<const chtRcThermalDiffusivityFvPatchScalarField>
            (
                shadowPatchField()
            );

        const scalarField diffNei =
            regionCouplePatch().interpolate(pCht.k());

        // Evaluate patch field by direct interpolation.
        // There is no need for distances, as two sets of data
        // are on top of each other.
        scalarField::operator=(2*diffOwn*diffNei/(diffOwn + diffNei));
    }
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::patchInterpolate
(
    surfaceScalarField& fField,
    const scalarField& pL
) const
{
    // Use patch value both in coupled and decoupled state
    // HJ, 22/Oct2024
    fField.boundaryField()[patch().index()] = *this;
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::patchInterpolate
(
    surfaceScalarField& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    // Use patch value both in coupled and decoupled state
    // HJ, 22/Oct2024
    fField.boundaryField()[patch().index()] = *this;
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        chtRcThermalDiffusivityFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
