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
    Henrik Rusche, Wikki GmbH.  All rights reserved

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
    jumpRegionCouplingFvPatchScalarField(p, iF),
    kName_("none"),
    radiation_(false),
    jump_(p.size(), scalar(0))
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpRegionCouplingFvPatchScalarField(p, iF, dict),
    kName_(dict.lookup("K")),
    radiation_(readBool(dict.lookup("radiation"))),
    jump_(p.size(), scalar(0))
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpRegionCouplingFvPatchScalarField(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_),
    jump_(ptf.jump_, mapper)
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpRegionCouplingFvPatchScalarField(ptf, iF),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::chtRcTemperatureFvPatchScalarField&
Foam::chtRcTemperatureFvPatchScalarField::shadowPatchField() const
{
    if
    (
        !isA<chtRcTemperatureFvPatchScalarField>
        (
            jumpRegionCouplingFvPatchScalarField::shadowPatchField()
        )
    )
    {
        FatalErrorInFunction
            << "Incorrect shadow patch type for patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " Should be chtRcTemperatureFvPatchScalarField.  Actual type "
            << jumpRegionCouplingFvPatchScalarField::shadowPatchField().type()
            << abort(FatalError);
    }

    return dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
    (
        jumpRegionCouplingFvPatchScalarField::shadowPatchField()
    );
}

Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::patchInternalT() const
{
    // For temperature, return patchInternalField
    return patchInternalField();
}


void Foam::chtRcTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchScalarField::autoMap(m);
    jump_.autoMap(m);
}


void Foam::chtRcTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchScalarField::rmap(ptf, addr);

    const chtRcTemperatureFvPatchScalarField& mptf =
        refCast<const chtRcTemperatureFvPatchScalarField>(ptf);

    jump_.rmap(mptf.jump_, addr);
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

    // Get neighbour T
    const chtRcTemperatureFvPatchScalarField& pCht =
        refCast<const chtRcTemperatureFvPatchScalarField>
        (
            shadowPatchField()
        );

    // Interpolate and add jump
    return regionCouplePatch().interpolate(pCht.patchInternalT()) + jump();
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

    if (radiation())
    {
        Info<< "Radiation active - updating radiative jump Qr for field "
            << this->dimensionedInternalField().name()
            << " on region " << this->dimensionedInternalField().mesh().name()
            << endl;

        // Get radiative heat flux. Qr [W/m^2]
        scalarField Qr = lookupPatchField<volScalarField, scalar>("Qr");

        // Get thermal conductivity
        const fvPatchScalarField& kpf =
            lookupPatchField<volScalarField, scalar>(kName_);

        // Update jump.  Units of jump_ are [K/m^2]
        jump_ = Qr/(kpf*patch().deltaCoeffs());
    }
    else
    {
        // No radiation
        jump_ = scalarField(patch().size(), scalar(0));
    }

    fvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::jump() const
{
    return jump_;
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
