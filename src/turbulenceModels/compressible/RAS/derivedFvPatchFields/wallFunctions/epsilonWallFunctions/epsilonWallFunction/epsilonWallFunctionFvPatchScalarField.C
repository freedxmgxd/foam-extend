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

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
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

void epsilonWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf, iF),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(rasModel.GName()))
    {
        InfoInFunction
            << "Cannot access " << rasModel.GName() << " field for patch "
            << patch().name() << ".  Evaluating as zeroGradient"
            << endl;

        fvPatchScalarField::updateCoeffs();
        zeroGradientFvPatchScalarField::evaluate();

        return;
    }

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu75 = pow(Cmu_, 0.75);
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);

    const scalarField& y = rasModel.y()[patch().index()];

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(rasModel.GName()));

    // Note: epsilon is now a refValue and set in
    // fixedInternalValueFvPatchField
    // HJ, 3/Aug/2011
    scalarField& epsilon = refValue();

    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();

    const scalarField& rhow =
        rasModel.thermo().rho()().boundaryField()[patch().index()];

    const scalarField& muw =
        rasModel.mu().boundaryField()[patch().index()];

    const scalarField& mutw =
        rasModel.mut()().boundaryField()[patch().index()];

    const fvPatchVectorField& Uw =
        rasModel.U().boundaryField()[patch().index()];

    const scalarField magGradUw = mag(Uw.snGrad());

    const labelList& faceCells = patch().faceCells();

    // Set epsilon and G
    forAll(mutw, faceI)
    {
        label faceCellI = faceCells[faceI];

        scalar yPlus =
            Cmu25*y[faceI]*sqrt(k[faceCellI])
           /(muw[faceI]/rhow[faceI]);

        // Note: epsilon is now a refValue and set in
        // fixedInternalValueFvPatchField
        // HJ, 3/Aug/2011
        epsilon[faceI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);

        if (yPlus > yPlusLam)
        {
            G[faceCellI] =
                (mutw[faceI] + muw[faceI])
               *magGradUw[faceI]
               *Cmu25*sqrt(k[faceCellI])
               /(kappa_*y[faceI]);
        }
        else
        {
            G[faceCellI] = 0.0;
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchScalarField::updateCoeffs();
}


void epsilonWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchScalarField::evaluate(commsType);
}


void epsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchScalarField::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
