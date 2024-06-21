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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "chtRcTemperatureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF),
    kName_("none"),
    radiation_(false)
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict),
    kName_(dict.lookup("K")),
    radiation_(readBool(dict.lookup("radiation")))
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_)
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtRcTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    regionCouplingFvPatchScalarField::autoMap(m);
}


void Foam::chtRcTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    regionCouplingFvPatchScalarField::rmap(ptf, addr);

    const chtRcTemperatureFvPatchScalarField& pf =
        refCast<const chtRcTemperatureFvPatchScalarField>(ptf);

    kName_ = pf.kName_;
    radiation_ = pf.radiation_;
}


const Foam::chtRcTemperatureFvPatchScalarField&
Foam::chtRcTemperatureFvPatchScalarField::shadowPatchRcTemperatureField() const
{
    if
    (
        !isA<chtRcTemperatureFvPatchScalarField>
        (
            regionCouplingFvPatchScalarField::shadowPatchField()
        )
    )
    {
        FatalErrorInFunction
            << "Incorrect shadow patch type for patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " Should be chtRcTemperatureFvPatchScalarField.  Actual type "
            << regionCouplingFvPatchScalarField::shadowPatchField().type()
            << abort(FatalError);
    }

    return dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::patchInternalT() const
{
    // For temperature, return patchInternalField
    return patchInternalField();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::patchNeighbourField() const
{
    // For temperature, return mapped patchInternalT from other side

    if (!regionCouplePatch().coupled())
    {
        FatalErrorInFunction
            << "Requested patchNeighbourField in decoupled state"
            << abort(FatalError);
    }

    // Get neighbour T, interpolate and add jump
    return
        regionCouplePatch().interpolate
        (
            this->shadowPatchRcTemperatureField().patchInternalT()
        );
}


void Foam::chtRcTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Check name of local field
    if (this->dimensionedInternalField().name() != "T")
    {
        FatalErrorInFunction
            << "chtRcTemperatureFvPatchScalarField boundary condition "
            << " for patch " << patch().name()
            << " works with field T only.  Field name: "
            << this->dimensionedInternalField().name()
            << abort(FatalError);
    }

    fvPatchScalarField::updateCoeffs();
}


void Foam::chtRcTemperatureFvPatchScalarField::manipulateMatrix
(
    fvScalarMatrix& matrix
)
{
    const fvPatch& p = patch();
    const scalarField& magSf = p.magSf();
    const labelList& faceCells = p.faceCells();
    scalarField& source = matrix.source();

    // Move Qr to solid
    if (radiation())
    {
        scalarField Qr =
            this->regionCouplePatch().interpolate
            (
                this->shadowPatchRcTemperatureField().
                lookupPatchField<volScalarField, scalar>("Qr")
            );

        forAll(faceCells, faceI)
        {
            source[faceCells[faceI]] += Qr[faceI]*magSf[faceI];
        }
    }
}


void Foam::chtRcTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("K") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("radiation") << radiation_ << token::END_STATEMENT << nl;
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    chtRcTemperatureFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
