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

#include "fluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    flux_(p.size(), pTraits<Type>::zero),
    reactivity_(p.size(), 0),
    gammaName_("gamma"),
    fieldBound_(pTraits<Type>::min, pTraits<Type>::max)
{
    this->gradient() = pTraits<Type>::zero;
}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    flux_("flux", dict, p.size()),
    reactivity_("reactivity", dict, p.size()),
    gammaName_(dict.lookup("gamma")),
    fieldBound_(dict.lookup("fieldBound"))
{
    // Set dummy gradient
    this->gradient() = pTraits<Type>::zero;

    // Read the value entry from the dictionary
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Essential entry 'value' missing"
            << abort(FatalIOError);
    }

    if (fieldBound_.first() > fieldBound_.second())
    {
        FatalIOErrorInFunction(dict)
            << "Bad field bound: " << fieldBound_
            << abort(FatalError);
    }
}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    flux_(ptf.flux_),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_),
    fieldBound_(ptf.fieldBound_)
{}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    flux_(ptf.flux_),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_),
    fieldBound_(ptf.fieldBound_)
{}


template<class Type>
fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    flux_(ptf.flux_),
    reactivity_(ptf.reactivity_),
    gammaName_(ptf.gammaName_),
    fieldBound_(ptf.fieldBound_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fluxFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchField<Type>::autoMap(m);
    flux_.autoMap(m);
    reactivity_.autoMap(m);
}


template<class Type>
void fluxFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const fluxFvPatchField<Type>& fptf =
        refCast<const fluxFvPatchField<Type> >(ptf);

    flux_.rmap(fptf.flux_, addr);
    reactivity_.rmap(fptf.reactivity_, addr);

    gammaName_ = fptf.gammaName_;
    fieldBound_ = fptf.fieldBound_;
}


template<class Type>
void fluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchScalarField& gammap =
        this->template lookupPatchField<volScalarField, scalar>(gammaName_);

    // Calculate minimum and maximum gradient
    Field<Type> minGrad =
        (fieldBound_.first() - this->patchInternalField())*
        this->patch().deltaCoeffs();

    Field<Type> maxGrad =
        (fieldBound_.second() - this->patchInternalField())*
        this->patch().deltaCoeffs();

    this->gradient() =
        max
        (
            minGrad,
            min
            (
                maxGrad,
                reactivity_*flux_/gammap
            )
        );

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void fluxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    flux_.writeEntry("flux", os);
    reactivity_.writeEntry("reactivity", os);
    os.writeKeyword("gamma") << gammaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldBound") << fieldBound_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
