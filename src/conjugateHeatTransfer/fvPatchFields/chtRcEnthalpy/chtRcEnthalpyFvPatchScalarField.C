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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::chtRcEnthalpyFvPatchScalarField::clearOut()
{
    deleteDemandDrivenData(CpBarPtr_);
    deleteDemandDrivenData(H0Ptr_);
}


void Foam::chtRcEnthalpyFvPatchScalarField::checkCpH0() const
{
    // Perform conversion from hs to T using stored H0 and Cp
    // linearisation
    if (!CpBarPtr_ || !H0Ptr_)
    {
        FatalErrorInFunction
            << "H0 and Cp not available.  Call updateCoeffs before evaluating"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(p, iF),
    CpBarPtr_(nullptr),
    H0Ptr_(nullptr)
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
    CpBarPtr_(nullptr),
    H0Ptr_(nullptr)
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
    CpBarPtr_(nullptr),
    H0Ptr_(nullptr)
{}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const chtRcEnthalpyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(ptf, iF),
    CpBarPtr_(nullptr),
    H0Ptr_(nullptr)
{}


Foam::chtRcEnthalpyFvPatchScalarField::chtRcEnthalpyFvPatchScalarField
(
    const chtRcEnthalpyFvPatchScalarField& ptf
)
:
    chtRcTemperatureFvPatchScalarField(ptf),
    CpBarPtr_(nullptr),
    H0Ptr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chtRcEnthalpyFvPatchScalarField::~chtRcEnthalpyFvPatchScalarField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtRcEnthalpyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    chtRcTemperatureFvPatchScalarField::autoMap(m);

    // Clear Cp and H0
    clearOut();
}


void Foam::chtRcEnthalpyFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    chtRcTemperatureFvPatchScalarField::rmap(ptf, addr);

    // Clear Cp and H0
    clearOut();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcEnthalpyFvPatchScalarField::patchInternalT() const
{
    // Check if Cp and H0 are available
    checkCpH0();

    const scalarField& CpBar = *CpBarPtr_;
    const scalarField& H0 = *H0Ptr_;
    
    return (patchInternalField() - H0)/CpBar;
}


Foam::tmp<Foam::scalarField>
Foam::chtRcEnthalpyFvPatchScalarField::patchNeighbourField() const
{
    // Check if Cp and H0 are available
    checkCpH0();

    const scalarField& CpBar = *CpBarPtr_;
    const scalarField& H0 = *H0Ptr_;
    
    // Get neighbour T
    const chtRcTemperatureFvPatchScalarField& pCht =
        refCast<const chtRcTemperatureFvPatchScalarField>
        (
            shadowPatchField()
        );

    // Interpolate and add jump
    return
        H0
      + regionCouplePatch().interpolate(pCht.patchInternalT())*CpBar
      + jump();
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

        // Allocate fields if needed
        if (!CpBarPtr_)
        {
            CpBarPtr_ = new scalarField(patch().size());
        }
        scalarField& CpBar = *CpBarPtr_;

        if (!H0Ptr_)
        {
            H0Ptr_ = new scalarField(patch().size());
        }
        scalarField& H0 = *H0Ptr_;

        // Collect CpBar and H0 for linearisation
        CpBar = thermo.Cp(patchT, patch().faceCells());

        H0 = patchInternalField() - CpBar*patchT;
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
        // Check if Cp and H0 are available
        checkCpH0();

        const scalarField& CpBar = *CpBarPtr_;
        const scalarField& H0 = *H0Ptr_;
    
        // Prepare local matrix update buffer for the remote side.
        // Note that only remote side will have access to its psiInternal
        // as they are on different regions

        // Since interpolation needs to happen on the shadow, and within the
        // init, prepare interpolation for the other side.

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
                    (patch().patchInternalField(psiInternal) - H0)/CpBar
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
                    patch().patchInternalField(psiInternal)/CpBar
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
        // Check if Cp and H0 are available
        checkCpH0();

        const scalarField& CpBar = *CpBarPtr_;
        const scalarField& H0 = *H0Ptr_;
    
        scalarField pnf = shadowPatchField().matrixUpdateBuffer();

        if
        (
            reinterpret_cast<const void*>(&psiInternal)
         == reinterpret_cast<const void*>(&this->internalField())
        )
        {
            pnf = H0 + pnf*CpBar;
        }
        else
        {
            pnf *= CpBar;
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
