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
    Foam::perfectGas

Description
    Compressible liquid with constant compressibility equation of state.

    rho(p) = rho0 + psiP*(p - p0) + psiT*(T - T0)

SourceFiles
    linearLiquidI.H
    linearLiquid.C

\*---------------------------------------------------------------------------*/

#ifndef linearLiquid_H
#define linearLiquid_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class linearLiquid Declaration
\*---------------------------------------------------------------------------*/

class linearLiquid
:
    public specie
{

public:

    // Private data

        //- Reference density
        //  Note: input in [kg/m3], but internally uses [kg/m3/kmol]
        scalar rho0_;

        //- Reference pressure
        scalar p0_;

        //- Pressure compressibility
        scalar psiP_;

        //- Reference temperature
        scalar T0_;

        //- Temperature compressibility
        scalar psiT_;


    // Constructors

        //- Construct from components
        inline linearLiquid
        (
            const specie& sp,
            const scalar rho0,
            const scalar p0,
            const scalar psiP,
            const scalar T0,
            const scalar psiT
        );

        //- Construct from Istream
        linearLiquid(Istream&);

        //- Construct as named copy
        inline linearLiquid(const word& name, const linearLiquid&);

        //- Construct and return a clone
        inline autoPtr<linearLiquid> clone() const;

        // Selector from Istream
        inline static autoPtr<linearLiquid> New(Istream& is);


    // Member functions

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;


    // Member operators

        inline void operator+=(const linearLiquid&);
        inline void operator-=(const linearLiquid&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend linearLiquid operator+
        (
            const linearLiquid&,
            const linearLiquid&
        );

        inline friend linearLiquid operator-
        (
            const linearLiquid&,
            const linearLiquid&
        );

        inline friend linearLiquid operator*
        (
            const scalar s,
            const linearLiquid&
        );

        inline friend linearLiquid operator==
        (
            const linearLiquid&,
            const linearLiquid&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const linearLiquid&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "linearLiquidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
