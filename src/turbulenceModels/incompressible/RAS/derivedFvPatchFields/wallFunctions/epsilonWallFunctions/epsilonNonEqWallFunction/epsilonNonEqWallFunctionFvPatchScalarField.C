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

#include "epsilonNonEqWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonNonEqWallFunctionFvPatchScalarField::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void epsilonNonEqWallFunctionFvPatchScalarField::
writeLocalEntries(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41))
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf, iF),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonNonEqWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    // Epsilon is now a refValue (derived from fixedInternalValueFvPatchField)
    scalarField& epsilon = refValue();

    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();

    const scalarField& nuw =
        rasModel.nu().boundaryField()[patch().index()];

    const scalarField& nutw =
        rasModel.nut()().boundaryField()[patch().index()];

    const fvPatchVectorField& Uw =
        rasModel.U().boundaryField()[patch().index()];

    // Get face cells
    const labelUList& fc = patch().faceCells();

    // Set epsilon
    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];

        // Viscous sublayer thickness
        const scalar yVis = nuw[faceI]*11.225/(Cmu25*sqrt(k[faceCellI]));

        if (y[faceI] > yVis)
        {
           epsilon[faceI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);
        }
        else
        {
           epsilon[faceI] = 2.0*nuw[faceI]*k[faceCellI]/sqr(y[faceI]);
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchScalarField::updateCoeffs();
}


void epsilonNonEqWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchScalarField::evaluate(commsType);
}


void epsilonNonEqWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchScalarField::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonNonEqWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
