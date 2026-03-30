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

#include "labelBlockCoeff.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::label>::BlockCoeff()
:
    scalarCoeff_(pTraits<label>::zero)
{}


Foam::BlockCoeff<Foam::label>::BlockCoeff(const BlockCoeff<label>& f)
:
    scalarCoeff_(f.scalarCoeff_)
{}


Foam::BlockCoeff<Foam::label>::BlockCoeff(Istream& is)
:
    scalarCoeff_(readLabel(is))
{}


Foam::BlockCoeff<Foam::label> Foam::BlockCoeff<Foam::label>::clone() const
{
    return BlockCoeff<label>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::label>::~BlockCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockCoeffBase::activeLevel
Foam::BlockCoeff<Foam::label>::activeType() const
{
    return blockCoeffBase::SCALAR;
}


Foam::label Foam::BlockCoeff<Foam::label>::component(const direction) const
{
    return scalarCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::BlockCoeff<Foam::label>::operator=(const BlockCoeff<label>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    scalarCoeff_ = f.scalarCoeff_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const BlockCoeff<label>& f)
{
    os << f.scalarCoeff_;

    return os;
}


// ************************************************************************* //
