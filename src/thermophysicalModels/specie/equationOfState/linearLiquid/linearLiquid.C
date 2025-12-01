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

Description
    Compressible liquid with constant compressibility equation of state.

\*---------------------------------------------------------------------------*/

#include "linearLiquid.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearLiquid::linearLiquid(Istream& is)
:
    specie(is),
    rho0_(readScalar(is)),
    p0_(readScalar(is)),
    psiP_(readScalar(is)),
    T0_(readScalar(is)),
    psiT_(readScalar(is))
{
    is.check("linearLiquid::linearLiquid(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const linearLiquid& ll)
{
    os  << static_cast<const specie&>(ll)
        << ll.rho0_ << token::SPACE
        << ll.p0_ << token::SPACE
        << ll.psiP_
        << ll.T0_ << token::SPACE
        << ll.psiT_;

    os.check("Ostream& operator<<(Ostream& os, const linearLiquid& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
