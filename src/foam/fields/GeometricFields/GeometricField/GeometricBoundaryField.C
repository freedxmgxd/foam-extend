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

#include "emptyPolyPatch.H"
#include "commSchedule.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (debug)
    {
        InfoInFunction
            << "copy single type"
            << endl;
    }

    forAll(bmesh_, patchI)
    {
        this->set
        (
            patchI,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchI],
                field
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const wordList& patchFieldTypes
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (debug)
    {
        InfoInFunction
            << "copy with list of types"
            << endl;
    }

    if (patchFieldTypes.size() != this->size())
    {
        FatalErrorInFunction
            << "Incorrect number of patch type specifications given" << nl
            << "    Number of patches in mesh = " << bmesh.size()
            << " number of patch type specifications = "
            << patchFieldTypes.size()
            << abort(FatalError);
    }

    forAll(bmesh_, patchI)
    {
        this->set
        (
            patchI,
            PatchField<Type>::New
            (
                patchFieldTypes[patchI],
                bmesh_[patchI],
                field
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const PtrList<PatchField<Type> >& ptfl
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (debug)
    {
        InfoInFunction
            << "copy boundary field"
            << endl;
    }

    forAll(bmesh_, patchI)
    {
        if (ptfl.set(patchI))
        {
            this->set(patchI, ptfl[patchI].clone(field));
        }
        else
        {
            FatalErrorInFunction
                << "ptfl not set for index " << patchI
                << abort(FatalError);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const DimensionedField<Type, GeoMesh>& field,
    const typename GeometricField<Type, PatchField, GeoMesh>::
    GeometricBoundaryField& btf
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    if (debug)
    {
        InfoInFunction
            << "copy with boundary field with new internal field"
            << endl;
    }

    forAll(bmesh_, patchI)
    {
        if (btf.set(patchI))
        {
            this->set(patchI, btf[patchI].clone(field));
        }
        else
        {
            FatalErrorInFunction
                << "btf not set for index " << patchI
                << abort(FatalError);
        }
    }
}


// Construct as copy
// Dangerous because Field may be set to a field which gets deleted.
// Need new type of GeometricBoundaryField, one which IS part of a geometric
// field for which snGrad etc. may be called and a free standing
// GeometricBoundaryField for which such operations are unavailable.
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const typename GeometricField<Type, PatchField, GeoMesh>::
    GeometricBoundaryField& btf
)
:
    FieldField<PatchField, Type>(btf),
    bmesh_(btf.bmesh_)
{
    if (debug)
    {
        InfoInFunction
            << "GeometricBoundaryField copy"
            << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const dictionary& dict
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (debug)
    {
        InfoInFunction
            << endl;
    }

    forAll(bmesh_, patchI)
    {
        if (bmesh_[patchI].type() != emptyPolyPatch::typeName)
        {
            this->set
            (
                patchI,
                PatchField<Type>::New
                (
                    bmesh_[patchI],
                    field,
                    dict.subDict(bmesh_[patchI].name())
                )
            );
        }
        else
        {
            this->set
            (
                patchI,
                PatchField<Type>::New
                (
                    emptyPolyPatch::typeName,
                    bmesh_[patchI],
                    field
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
updateCoeffs()
{
    if (debug)
    {
        InfoInFunction
            << "updateCoeffs"
            << endl;
    }

    forAll(*this, patchI)
    {
        this->operator[](patchI).updateCoeffs();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
evaluate()
{
    if (debug)
    {
        InfoInFunction
            << "evaluate"
            << endl;
    }

    if
    (
        Pstream::defaultComms() == Pstream::blocking
     || Pstream::defaultComms() == Pstream::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(*this, patchI)
        {
            this->operator[](patchI).initEvaluate
            (
                Pstream::defaultComms()
            );
        }

        // Block for any outstanding requests
        if (Pstream::defaultComms() == Pstream::nonBlocking)
        {
            Pstream::waitRequests(nReq);
        }

        forAll(*this, patchI)
        {
            this->operator[](patchI).evaluate
            (
                Pstream::defaultComms()
            );
        }
    }
    else if (Pstream::defaultComms() == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            if (patchSchedule[patchEvali].init)
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .initEvaluate(Pstream::scheduled);
            }
            else
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .evaluate(Pstream::scheduled);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultComms()]
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
updateCoupledPatchFields() const
{
    if (debug)
    {
        InfoInFunction
            << "updateCoupledPatchFields"
            << endl;
    }

    bool couplesUpdated = true;

    forAll (*this, patchI)
    {
        if (this->operator[](patchI).coupled())
        {
            couplesUpdated &= this->operator[](patchI).couplesUpdated();
        }
    }

    reduce(couplesUpdated, andOp<bool>());

    // Update couples if not updated
    if (!couplesUpdated)
    {
        if (debug)
        {
            InfoInFunction
                << "Updating couples"
                << endl;
        }

        // Note: casting away const in initEvaluate and evaluate, as this
        // only changes the cached patchNeighbourField value and not the
        // field itself.  HJ, 10/Sep/2021
        if
        (
            Pstream::defaultComms() == Pstream::blocking
         || Pstream::defaultComms() == Pstream::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll (*this, patchI)
            {
                if (this->operator[](patchI).coupled())
                {
                    const_cast<PatchField<Type>&>(this->operator[](patchI))
                        .initEvaluate(Pstream::defaultComms());
                }
            }

            // Block for any outstanding requests
            if (Pstream::defaultComms() == Pstream::nonBlocking)
            {
                Pstream::waitRequests(nReq);
            }

            forAll (*this, patchI)
            {
                if (this->operator[](patchI).coupled())
                {
                    const_cast<PatchField<Type>&>(this->operator[](patchI))
                        .evaluate(Pstream::defaultComms());
                }
            }
        }
        else if (Pstream::defaultComms() == Pstream::scheduled)
        {
            const lduSchedule& patchSchedule =
                bmesh_.mesh().globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                if (patchSchedule[patchEvali].init)
                {
                    if
                    (
                        this->operator[](patchSchedule[patchEvali].patch)
                        .coupled()
                    )
                    {
                        const_cast<PatchField<Type>&>
                        (
                            this->operator[](patchSchedule[patchEvali].patch)
                        ).initEvaluate(Pstream::scheduled);
                    }
                }
                else
                {
                    if
                    (
                        this->operator[](patchSchedule[patchEvali].patch)
                        .coupled()
                    )
                    {
                        const_cast<PatchField<Type>&>
                        (
                            this->operator[](patchSchedule[patchEvali].patch)
                        ).evaluate(Pstream::scheduled);
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultComms()]
                    << exit(FatalError);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
types() const
{
    const FieldField<PatchField, Type>& pff = *this;

    wordList Types(pff.size());

    forAll(pff, patchI)
    {
        Types[patchI] = pff[patchI].type();
    }

    return Types;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
boundaryInternalField() const
{
    typename GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField
        BoundaryInternalField(*this);

    forAll(BoundaryInternalField, patchI)
    {
        BoundaryInternalField[patchI] ==
            this->operator[](patchI).patchInternalField();
    }

    return BoundaryInternalField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::lduInterfaceFieldPtrsList
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
interfaces() const
{
    lduInterfaceFieldPtrsList interfaces(this->size());

    forAll (interfaces, patchI)
    {
        if (isA<lduInterfaceField>(this->operator[](patchI)))
        {
            interfaces.set
            (
                patchI,
                &refCast<const lduInterfaceField>(this->operator[](patchI))
            );
        }
    }

    return interfaces;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename Foam::BlockLduInterfaceFieldPtrsList<Type>::Type
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
blockInterfaces() const
{
    typename BlockLduInterfaceFieldPtrsList<Type>::Type interfaces
    (
        this->size()
    );

    forAll (interfaces, patchI)
    {
        if (isA<BlockLduInterfaceField<Type> >(this->operator[](patchI)))
        {
            interfaces.set
            (
                patchI,
                &refCast<const BlockLduInterfaceField<Type> >
                (
                    this->operator[](patchI)
                )
            );
        }
    }

    return interfaces;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
clearCaches() const
{
    // Clear caches on all boundary patches
    forAll(*this, patchI)
    {
        this->operator[](patchI).clearCaches();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
writeEntry(const word& keyword, Ostream& os) const
{
    os  << keyword << nl << token::BEGIN_BLOCK << incrIndent << nl;

    forAll(*this, patchI)
    {
        os  << indent << this->operator[](patchI).patch().name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << this->operator[](patchI) << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_BLOCK << endl;

    // Check state of IOstream
    os.check
    (
        "GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::"
        "writeEntry(const word& keyword, Ostream& os) const"
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator=
(
    const typename GeometricField<Type, PatchField, GeoMesh>::
    GeometricBoundaryField& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator=
(
    const FieldField<PatchField, Type>& ptff
)
{
    FieldField<PatchField, Type>::operator=(ptff);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator=
(
    const Type& t
)
{
    FieldField<PatchField, Type>::operator=(t);
}


// Forced assignments
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator==
(
    const typename GeometricField<Type, PatchField, GeoMesh>::
    GeometricBoundaryField& bf
)
{
    forAll((*this), patchI)
    {
        this->operator[](patchI) == bf[patchI];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator==
(
    const FieldField<PatchField, Type>& ptff
)
{
    forAll((*this), patchI)
    {
        this->operator[](patchI) == ptff[patchI];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField::
operator==
(
    const Type& t
)
{
    forAll((*this), patchI)
    {
        this->operator[](patchI) == t;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const typename GeometricField<Type, PatchField, GeoMesh>::
    GeometricBoundaryField& bf
)
{
    os << static_cast<const FieldField<PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //
