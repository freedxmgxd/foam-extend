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

Class
    labelCoeffField

Description

\*---------------------------------------------------------------------------*/

#include "labelCoeffField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const char* const Foam::CoeffField<Foam::label>::typeName("CoeffField");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockCoeffBase::activeLevel
Foam::CoeffField<Foam::label>::activeType() const
{
    return blockCoeffBase::SCALAR;
}


Foam::tmp<Foam::CoeffField<Foam::label> >
Foam::CoeffField<Foam::label>::transpose() const
{
    return tmp<CoeffField<label> >(new CoeffField<label>(*this));
}


const Foam::labelField&
Foam::CoeffField<Foam::label>::asScalar() const
{
    return *this;
}


Foam::labelField&
Foam::CoeffField<Foam::label>::asScalar()
{
    return *this;
}


const Foam::labelField&
Foam::CoeffField<Foam::label>::asLinear() const
{
    return *this;
}


Foam::labelField&
Foam::CoeffField<Foam::label>::asLinear()
{
    return *this;
}


const Foam::labelField&
Foam::CoeffField<Foam::label>::asSquare() const
{
    return *this;
}


Foam::labelField&
Foam::CoeffField<Foam::label>::asSquare()
{
    return *this;
}


Foam::BlockCoeff<Foam::label>
Foam::CoeffField<Foam::label>::getCoeff(const label index) const
{
    BlockCoeff<label> result;

    result.asScalar() = (operator[](index));

    return result;
}


void Foam::CoeffField<Foam::label>::setCoeff
(
    const label index,
    const BlockCoeff<label>& coeff
)
{
    operator[](index) = coeff.asScalar();
}


void Foam::CoeffField<Foam::label>::getSubset
(
    CoeffField<label>& f,
    const label start,
    const label size
) const
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorInFunction
            << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    const labelField& localF = *this;

    forAll (f, fI)
    {
        f[fI] = localF[start + fI];
    }
}


void Foam::CoeffField<Foam::label>::getSubset
(
    CoeffField<label>& f,
    const labelList& addr
) const
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorInFunction
            << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    const labelField& localF = *this;

    forAll (f, fI)
    {
        f[fI] = localF[addr[fI]];
    }
}


void Foam::CoeffField<Foam::label>::setSubset
(
    const CoeffField<label>& f,
    const label start,
    const label size
)
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorInFunction
            << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    labelField& localF = *this;

    forAll (f, fI)
    {
        localF[start + fI] = f[fI];
    }
}


void Foam::CoeffField<Foam::label>::setSubset
(
    const CoeffField<label>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorInFunction
            << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    labelField& localF = this->asScalar();

    forAll (f, fI)
    {
        localF[addr[fI]] = f[fI];
    }
}


void Foam::CoeffField<Foam::label>::zeroOutSubset
(
    const label start,
    const label size
)
{
    labelField& localF = *this;

    for (label ffI = 0; ffI < size; ffI++)
    {
        localF[start + ffI] = pTraits<label>::zero;
    }
}


void Foam::CoeffField<Foam::label>::zeroOutSubset
(
    const labelList& addr
)
{
    labelField& localF = *this;

    forAll (addr, ffI)
    {
        localF[addr[ffI]] = pTraits<label>::zero;
    }
}


void Foam::CoeffField<Foam::label>::addSubset
(
    const CoeffField<label>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorInFunction
            << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    labelField& localF = this->asScalar();

    forAll (f, fI)
    {
        localF[addr[fI]] += f[fI];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::label>::operator=(const CoeffField<label>& f)
{
    labelField::operator=(f.asScalar());
}


void Foam::CoeffField<Foam::label>::operator=(const labelField& f)
{
    labelField::operator=(f);
}


void Foam::CoeffField<Foam::label>::operator=(const tmp<labelField>& tf)
{
    labelField::operator=(tf);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CoeffField<label>& f)
{
    os << static_cast<const labelField&>(f);

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<label> >& tf
)
{
    os << tf();
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<>
Foam::tmp<Foam::CoeffField<Foam::label> >
Foam::inv(const CoeffField<label>& f)
{
    notImplemented("Foam::inv(const CoeffField<label>& f)");

    return f;
}


template<>
void Foam::negate
(
    CoeffField<label>& f,
    const CoeffField<label>& f1
)
{
    f = f1;
    f.negate();
}


template<>
void Foam::multiply
(
    labelField& f,
    const CoeffField<label>& f1,
    const label& f2
)
{
    const labelField& sf = f1;
    f = sf*f2;
}


template<>
void Foam::multiply
(
    labelField& f,
    const CoeffField<label>& f1,
    const labelField& f2
)
{
    const labelField& sf = f1;
    f = sf*f2;
}


template<>
void Foam::multiply
(
    labelField& f,
    const labelField& f1,
    const CoeffField<label>& f2
)
{
    const labelField& sf = f2;
    f = f1*sf;
}


// ************************************************************************* //
