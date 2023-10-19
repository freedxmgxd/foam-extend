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

#include "inletOutletIbFvPatchField.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void inletOutletIbFvPatchField<Type>::updateIbValues()
{
    // Interpolate the values from tri surface using nearest triangle
    const labelList& nt = this->ibPatch().ibPolyPatch().nearestTri();

    this->refValue() = Field<Type>(triInletValue_, nt);

    this->refGrad() = Field<Type>(this->ibPatch().size(), pTraits<Type>::zero);

    this->valueFraction() = scalarField(this->ibPatch().size(), scalar(0));

    // Resize the field and initialise to a reasonable value.
    // Reconsider.  HJ, 16/Jun/2022
    Field<Type>::operator=
    (
        Field<Type>(this->ibPatch().size(), pTraits<Type>::zero)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inletOutletIbFvPatchField<Type>::inletOutletIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>(p, false, pTraits<Type>::zero),
    triInletValue_(this->ibPatch().ibMesh().size(), pTraits<Type>::zero)
{}


template<class Type>
inletOutletIbFvPatchField<Type>::inletOutletIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>
    (
        p,
        Switch(dict.lookup("setDeadValue")),
        pTraits<Type>(dict.lookup("deadValue"))
    ),
    triInletValue_("triInletValue", dict, this->ibPatch().ibMesh().size())
{
    // Since patch does not read a dictionary, the patch type needs to be read
    // manually.  HJ, 6/Sep/2018
    this->readPatchType(dict);

    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();

    inletOutletFvPatchField<Type>::evaluate();
}


template<class Type>
inletOutletIbFvPatchField<Type>::inletOutletIbFvPatchField
(
    const inletOutletIbFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>
    (
        p,
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triInletValue_(ptf.triInletValue())
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    this->setPatchType(ptf);

    // Re-interpolate the data related to immersed boundary is not safe
    // during decompose/reconstruct.  Set a value to avoid uninitialised
    // memory
    // this->updateIbValues();

    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1;

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(pTraits<Type>::zero);
}


template<class Type>
inletOutletIbFvPatchField<Type>::inletOutletIbFvPatchField
(
    const inletOutletIbFvPatchField<Type>& ptf
)
:
    inletOutletFvPatchField<Type>(ptf),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triInletValue_(ptf.triInletValue())
{
    this->setPatchType(ptf);
}


template<class Type>
inletOutletIbFvPatchField<Type>::inletOutletIbFvPatchField
(
    const inletOutletIbFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(ptf, iF),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triInletValue_(ptf.triInletValue())
{
    this->setPatchType(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void inletOutletIbFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Base fields do not map: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void inletOutletIbFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap: re-interpolate them from tri data
    const inletOutletIbFvPatchField<Type>& mptf =
        refCast<const inletOutletIbFvPatchField<Type> >(ptf);

    // Set rmap tri data
    triInletValue_ = mptf.triInletValue_;

    this->updateIbValues();
}


template<class Type>
void inletOutletIbFvPatchField<Type>::updateOnMotion()
{
    if (this->size() != this->ibPatch().size())
    {
        this->updateIbValues();

        // Note:
        // At this point, flux field has not been updated and it has the
        // wrong size.  Therefore, inletOutlet updateCoeffs will give the
        // wrong size of valueFraction.
        // To avoid the problem, the field is dimensioned to the right size in
        // updateIbValues().  HJ, 9/Sep/2022
    }
}


template<class Type>
void inletOutletIbFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    this->updateIbValues();

    // Set dead value
    this->setDeadValues(*this);

    // Evaluate fixed value condition
    inletOutletFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::inletOutletIbFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    this->setDeadValues(matrix);
}


template<class Type>
void inletOutletIbFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    // Resolve post-processing issues.  HJ, 1/Dec/2017
    if (this->phiName() != "phi")
    {
        os.writeKeyword("phi")
            << this->phiName() << token::END_STATEMENT << nl;
    }

    triInletValue_.writeEntry("triInletValue", os);
    immersedBoundaryFieldBase<Type>::writeDeadData(os);

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    this->writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
