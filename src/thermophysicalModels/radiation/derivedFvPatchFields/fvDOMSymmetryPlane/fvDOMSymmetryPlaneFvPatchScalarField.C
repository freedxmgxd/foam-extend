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

#include "fvDOMSymmetryPlaneFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvDOM.H"
#include "specularReflectionAddressing.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvDOMSymmetryPlaneFvPatchScalarField::
fvDOMSymmetryPlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefinedT")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


fvDOMSymmetryPlaneFvPatchScalarField::
fvDOMSymmetryPlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T"))
{
    if (dict.found("refValue"))
    {
        refValue() = scalarField("value", dict, p.size());

        fvPatchScalarField::operator=
        (
            refValue()
        );

        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        // Bugfix: Do not initialize from temperature because it is unavailable
        // when running, e.g. decomposePar and loading radiation as
        // shared library. Initialize to zero instead.
        // 26/Mar/2014 - DC

        refValue() = 0;

        refGrad() = 0;
        valueFraction() = 1;

        fvPatchScalarField::operator=(refValue());
    }
}


fvDOMSymmetryPlaneFvPatchScalarField::
fvDOMSymmetryPlaneFvPatchScalarField
(
    const fvDOMSymmetryPlaneFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_)
{}


fvDOMSymmetryPlaneFvPatchScalarField::
fvDOMSymmetryPlaneFvPatchScalarField
(
    const fvDOMSymmetryPlaneFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_)
{}


fvDOMSymmetryPlaneFvPatchScalarField::
fvDOMSymmetryPlaneFvPatchScalarField
(
    const fvDOMSymmetryPlaneFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fvDOMSymmetryPlaneFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If patch is empty, do nothing
    if (empty())
    {
        return;
    }
    
    const label patchI = patch().index();

    // Access radiation model
    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom = dynamic_cast<const fvDOM&>(radiation);

    // Get rayId and lambda Id for this ray
    label rayId = -1;
    label lambdaId = -1;

    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);

    // Make shortcut to ray belonging to this field
    const radiativeIntensityRay& ray = dom.IRay(rayId);

    // Get face normals
    const vectorField nHat = patch().nf();

    // Calculate cos of incoming angle of current ray with every face
    const scalarField incomingAngle = nHat & ray.dAve();

    // Set to zeroGradient (=0; incomingAngle > 0) for faces with incoming rays
    // and to fixedValue (=1; incomingAngle < 0) for outgoing rays
    valueFraction() = neg(incomingAngle);

    // Set refGrad to zero (outgoing ray)
    refGrad() = 0;

    // Get specular reflection addressing
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const specularReflectionAddressing& specReflection =
        specularReflectionAddressing::New(mesh);

    // Get reflected ray index
    const label reflectedRayId =
        specReflection.patchReflectionAddr(patchI, rayId);

    if (reflectedRayId != -1)
    {
        // Get outgoing reflected ray
        const radiativeIntensityRay& outRay = dom.IRay(reflectedRayId);

        // Calculate cosine of angle between face and ray
        scalarField outgoingAngle = nHat & outRay.dAve();

        const fvPatchScalarField& curPatchOutRay =
            outRay.ILambda(lambdaId).boundaryField()[patchI];

        // Calculate reflected value.  Note change of sign
        refValue() = -curPatchOutRay*outgoingAngle/incomingAngle;
    }
    else
    {
        // No reflected ray or outgoing ray.  Set value to zero
        refValue() = 0;
    }

    // Update boundary field now, so values for incoming and outgoing rays
    // are in balance
    scalarField::operator=
    (
        valueFraction()*refValue()
      + (1.0 - valueFraction())*
        (
            patchInternalField()
          + refGrad()/patch().deltaCoeffs()
        )
    );

    mixedFvPatchScalarField::updateCoeffs();
}


void fvDOMSymmetryPlaneFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("patchType")
        << symmetryFvPatch::typeName << token::END_STATEMENT << nl;

    writeEntryIfDifferent(os, "T", word("T"), TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fvDOMSymmetryPlaneFvPatchScalarField
);

} // End namespace radiation
} // End namespace Foam

// ************************************************************************* //
