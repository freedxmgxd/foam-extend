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

\*---------------------------------------------------------------------------*/

#include "incompressiblePorousBafflePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


#include "RASModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::incompressiblePorousBafflePressureFvPatchScalarField::calcJump() const
{
    // Get flux
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    // Calculate normal velocity
    scalarField Un =
        scalarField::subField(phip, size()/2)/
        scalarField::subField(patch().magSf(), size()/2);

    scalarField magUn = Foam::min(mag(Un), Umax_);

    // Store old jump for under-relaxation
    scalarField jumpOld = jump_;

    // Calculate jump
    if
    (
        phip.dimensionedInternalField().dimensions()
     == dimVelocity*dimArea
    )
    {
        // Incompressible flux.  Lookup incompressible turbulence model
        const incompressible::RASModel& turbModel =
            this->dimensionedInternalField().mesh().lookupObject
            <
                incompressible::RASModel
            >("RASProperties");


         // Get effective viscosity for the jump.  Note: half-sized field
         const scalarField nuEffw =
             scalarField::subField
             (
                 turbModel.nuEff()().boundaryField()[patch().index()],
                 size()/2
             );

         // Calculate incompressible pressure jump
         // jump_ = -sign(Un)*magUn*
         //     (
         //         inertialCoeff_*nuEffw
         //       + DarcyCoeff_*0.5*magUn
         //     );

	 // MH - following the equation delta_p = a*u + b*u^2
	 // a and b have to take into accound the thickness of the baffle
	 // and viscosity
         jump_ = -sign(Un)*magUn*(A_+ B_*magUn);
    }
    else
    {
        FatalErrorInFunction
            << "Cannot recognise flux field dimensions: "
            << phip.dimensionedInternalField().dimensions()
            << abort(FatalError);
    }

    // Look up under-relaxation factor from dictionary
    const scalar alpha =
        patch().boundaryMesh().mesh().solutionDict().fieldRelaxationFactor
        (
            dimensionedInternalField().name() + "Jump"
        );

    jump_ = alpha*jump_ + (1 - alpha)*jumpOld;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePorousBafflePressureFvPatchScalarField::
incompressiblePorousBafflePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchScalarField(p, iF),
    phiName_("undefined"),
    A_(0),
    B_(0),
    Umax_(0),
    jump_(p.size()/2, 0),
    curTimeIndex_(-1)
{}


Foam::incompressiblePorousBafflePressureFvPatchScalarField::
incompressiblePorousBafflePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchScalarField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    A_(readScalar(dict.lookup("A"))),
    B_(readScalar(dict.lookup("B"))),
    Umax_(readScalar(dict.lookup("Umax"))),
    jump_(p.size()/2, 0),
    curTimeIndex_(-1)
{}


Foam::incompressiblePorousBafflePressureFvPatchScalarField::
incompressiblePorousBafflePressureFvPatchScalarField
(
    const incompressiblePorousBafflePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    A_(ptf.A_),
    B_(ptf.B_),
    Umax_(ptf.Umax_),
    jump_(p.size()/2, 0),
    curTimeIndex_(-1)
{}


Foam::incompressiblePorousBafflePressureFvPatchScalarField::
incompressiblePorousBafflePressureFvPatchScalarField
(
    const incompressiblePorousBafflePressureFvPatchScalarField& ptf
)
:
    cyclicLduInterfaceField(ptf),
    jumpCyclicFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    A_(ptf.A_),
    B_(ptf.B_),
    Umax_(ptf.Umax_),
    jump_(ptf.size()/2, 0),
    curTimeIndex_(-1)
{}


Foam::incompressiblePorousBafflePressureFvPatchScalarField::
incompressiblePorousBafflePressureFvPatchScalarField
(
    const incompressiblePorousBafflePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    A_(ptf.A_),
    B_(ptf.B_),
    Umax_(ptf.Umax_),
    jump_(ptf.size()/2, 0),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::incompressiblePorousBafflePressureFvPatchScalarField::jump() const
{
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        calcJump();
    }

    return jump_;
}

void Foam::incompressiblePorousBafflePressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("B") << B_
        << token::END_STATEMENT << nl;
    os.writeKeyword("Umax") << Umax_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        incompressiblePorousBafflePressureFvPatchScalarField
    );
} // End namespace Foam


// ************************************************************************* //
