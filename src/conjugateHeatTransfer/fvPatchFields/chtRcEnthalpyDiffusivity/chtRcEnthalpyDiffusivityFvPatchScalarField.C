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
    Removed in entirety: Henrik Rusche, Wikki GmbH.  All rights reserved
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcEnthalpyDiffusivityFvPatchScalarField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::
chtRcEnthalpyDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcThermalDiffusivityFvPatchScalarField(p, iF)
{}


Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::
chtRcEnthalpyDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRcThermalDiffusivityFvPatchScalarField(p, iF, dict)
{}


Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::
chtRcEnthalpyDiffusivityFvPatchScalarField
(
    const chtRcEnthalpyDiffusivityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRcThermalDiffusivityFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::
chtRcEnthalpyDiffusivityFvPatchScalarField
(
    const chtRcEnthalpyDiffusivityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcThermalDiffusivityFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::k() const
{
    return decoupledField()*Cp();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::Cp() const
{
    // For thermal diffusivity (solving h or hs equation), return decoupled k

    // Find thermo if it exists
    if (db().foundObject<basicThermo>("thermophysicalProperties"))
    {
        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        const scalarField patchT =
            patch().lookupPatchField<volScalarField, scalar>("T");

        return thermo.Cp(patchT, patch().faceCells());
    }
    else
    {
        InfoInFunction
            << "Cannot find.  Simple return"
            << endl;

        // Note: returning zero will blow up k calculation
        return tmp<scalarField>(new scalarField(patch().size(), scalar(0)));
    }
}


void Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::initEvaluate
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

        const scalarField patchCp = Cp();

        // Get local decoupled patch field values, convert to k
        // k = KEff*Cp
        // const scalarField diffOwn = decoupledField()*patchCp;
        const scalarField diffOwn = k();

        // Interpolate neighbour decoupled patch field values

        // Get neighbour k
        const chtRcThermalDiffusivityFvPatchScalarField& pCht =
            refCast<const chtRcThermalDiffusivityFvPatchScalarField>
            (
                shadowPatchField()
            );

        const scalarField diffNei =
            regionCouplePatch().interpolate(pCht.k());

        // Evaluate patch field by harmonic interpolation of k.
        // Divide by Cp to create thermal diffusivity for hs
        // There is no need for distances, as two sets of data
        // are on top of each other.
        scalarField interpolatedK = 2*diffOwn*diffNei/(diffOwn + diffNei);

        scalarField::operator=(interpolatedK/patchCp);
    }
}


void Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcEnthalpyDiffusivityFvPatchScalarField::write(Ostream& os) const
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
        chtRcEnthalpyDiffusivityFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
