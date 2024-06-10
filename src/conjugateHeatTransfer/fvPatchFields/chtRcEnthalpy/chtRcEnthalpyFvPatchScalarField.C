/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved
    Bugs fixed by Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "chtRcEnthalpyFvPatchScalarField.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(p, iF),
    Cp_(0),
    H0_(0)
{
    // Set coupling conditions from T.  Reconsider.
    // HJ, 10/Jun/2024
    const chtRcTemperatureFvPatchScalarField& T =
        refCast<const chtRcTemperatureFvPatchScalarField>
        (
            lookupPatchField<volScalarField, scalar>(remoteFieldName())
        );

    setRemoteFieldName(T.remoteFieldName());
    radiation() = T.radiation();
    kName() = T.kName();
}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRcTemperatureFvPatchScalarField(p, iF, dict),
    Cp_(0),
    H0_(0)
{}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const chtRcEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRcTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    Cp_(0),
    H0_(0)
{}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const chtRcEnthalpyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(ptf, iF),
    Cp_(ptf.Cp_),
    H0_(ptf.H0_)
{}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const chtRcEnthalpyFvPatchScalarField& ptf
)
:
    chtRcTemperatureFvPatchScalarField(ptf),
    Cp_(ptf.Cp_),
    H0_(ptf.H0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtRcEnthalpyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    chtRcTemperatureFvPatchScalarField::autoMap(m);

    // Clear Cp and H0
    Cp_.clear();
    H0_.clear();
}


void Foam::chtRcEnthalpyFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    chtRcTemperatureFvPatchScalarField::rmap(ptf, addr);

    // Clear Cp and H0
    Cp_.clear();
    H0_.clear();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcEnthalpyFvPatchScalarField::patchNeighbourField() const
{
    // Find thermo if it exists
    if (db().foundObject<basicThermo>("thermophysicalProperties"))
    {
        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        // Check name of shadow patch field
        if (shadowPatchField().dimensionedInternalField().name() != "T")
        {
            FatalErrorInFunction
                << "chtRcEnthalpyFvPatchScalarField boundary condition "
                    << " for patch " << patch().name()
            << " works with T shadow patch field only.  Field name: "
            << this->shadowPatchField().dimensionedInternalField().name()
            << abort(FatalError);
        }

        // Get shadow temperature
        scalarField shadowT =
            this->regionCouplePatch().interpolate
            (
                this->shadowPatchField().patchInternalField()
            );

        // Prepare return
        tmp<scalarField> tshadowH(new scalarField(patch().size(), scalar(0)));
        scalarField& shadowH = tshadowH.ref();

        // Pick correct form of h-T conversion based on field name
        // HJ, 16/Jun/2018
        if (dimensionedInternalField().name() == "h")
        {
            shadowH = thermo.h(shadowT, patch().faceCells());
        }
        else if (dimensionedInternalField().name() == "hs")
        {
            shadowH = thermo.hs(shadowT, patch().faceCells());
        }
        else
        {
            FatalErrorInFunction
                << "Cannot recognise enthalpy field type for h-T conversion"
                << abort(FatalError);
        }

        // Add jump
        shadowH += jump();

        return tshadowH;
    }
    else
    {
        // Dummy return if thermo cannot be found
        InfoInFunction
            << "Cannot find thermo.  Simple return"
                << endl;

        return *this;
    }
}


void Foam::chtRcEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Check name of local field
    if
    (
        this->dimensionedInternalField().name() != "h"
     && this->dimensionedInternalField().name() != "hs"
    )
    {
        FatalErrorInFunction
            << "chtRcEnthalpyFvPatchScalarField boundary condition "
            << " for patch " << patch().name()
            << " works with h or hs field only.  Field name: "
            << this->dimensionedInternalField().name()
            << abort(FatalError);
    }

    // Find thermo if it exists
    if (db().foundObject<basicThermo>("thermophysicalProperties"))
    {
        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        const scalarField patchT =
            patch().lookupPatchField<volScalarField, scalar>("T");

        Cp_ = thermo.Cp(patchT, patch().faceCells());
        H0_ = patchInternalField() - Cp_*patchT;
    }
    else
    {
        InfoInFunction
            << "Cannot find thermo.  Simple return"
            << endl;
    }

    fvPatchScalarField::updateCoeffs();
}


// Initialise neighbour processor internal cell data
void Foam::chtRcEnthalpyFvPatchScalarField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& matrix,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (regionCouplePatch().coupled())
    {
        // Prepare local matrix update buffer for the remote side.
        // Note that only remote side will have access to its psiInternal
        // as they are on different regions

        // Since interpolation needs to happen on the shadow, and within the
        // init, prepare interpolation for the other side.

        // Initialise Cp for this linear solution - Reset in initEvaluate.

        // Add void pointer cast to keep compiler happy when instantiated
        // for vector/tensor fields.  HJ, 4/Jun/2013
        if
        (
            reinterpret_cast<const void*>(&psiInternal)
         == reinterpret_cast<const void*>(&this->internalField())
        )
        {
            setMatrixUpdateBuffer
            (
                shadowPatchField().regionCouplePatch().interpolate
                (
                    (patch().patchInternalField(psiInternal) - H0_)/Cp_
                    // Add jump
                  + jump()
                )
            );
        }
        else
        {
            setMatrixUpdateBuffer
            (
                shadowPatchField().regionCouplePatch().interpolate
                (
                    patch().patchInternalField(psiInternal)/Cp_
                )
                // Add jump
              + jump()
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "init matrix update called in detached state"
            << abort(FatalError);
    }
}


// Return matrix product for coupled boundary
void Foam::chtRcEnthalpyFvPatchScalarField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& matrix,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (regionCouplePatch().coupled())
    {
        scalarField pnf = shadowPatchField().matrixUpdateBuffer();

        if
        (
            reinterpret_cast<const void*>(&psiInternal)
         == reinterpret_cast<const void*>(&this->internalField())
        )
        {
            pnf = H0_ + pnf*Cp_;
        }
        else
        {
            pnf *= Cp_;
        }

        // Multiply the field by coefficients and add into the result
        const labelUList& fc = regionCouplePatch().faceCells();

        if (switchToLhs)
        {
            forAll (fc, elemI)
            {
                result[fc[elemI]] += coeffs[elemI]*pnf[elemI];
            }
        }
        else
        {
            forAll (fc, elemI)
            {
                result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    chtRcEnthalpyFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
