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

#include "linearLiquid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from components
inline linearLiquid::linearLiquid
(
    const specie& sp,
    const scalar rho0,
    const scalar p0,
    const scalar psiP,
    const scalar T0,
    const scalar psiT
)
:
    specie(sp),
    rho0_(rho0),
    p0_(p0),
    psiP_(psiP),
    T0_(T0),
    psiT_(psiT)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy
inline linearLiquid::linearLiquid(const word& name, const linearLiquid& ll)
:
    specie(name, ll),
    rho0_(ll.rho0_),
    p0_(ll.p0_),
    psiP_(ll.psiP_),
    T0_(ll.T0_),
    psiT_(ll.psiT_)
{}


// Construct and return a clone
inline autoPtr<linearLiquid> linearLiquid::clone() const
{
    return autoPtr<linearLiquid>(new linearLiquid(*this));
}


// Selector from Istream
inline autoPtr<linearLiquid> linearLiquid::New(Istream& is)
{
    return autoPtr<linearLiquid>(new linearLiquid(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return density [kg/m^3]
inline scalar linearLiquid::rho(scalar p, scalar T) const
{

    return rho0_ + psiP_*(p - p0_) + psiT_*(T - T0_) ;// /this->W();
}

//- Return compressibility rho/p [s^2/m^2]
inline scalar linearLiquid::psi(scalar, scalar T) const
{
    return psiP_;
}

//- Return compression factor []
inline scalar linearLiquid::Z(scalar, scalar) const
{
    return 1.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void linearLiquid::operator+=(const linearLiquid& ll)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ll);

    molr1 /= this->nMoles();
    scalar molr2 = ll.nMoles()/this->nMoles();

    // p0_ and T0_ remain unchanged: use first reference
    
    // Adjust rho0 by aligning the ll rho0_ with current
    rho0_ = molr1*rho0_ + molr2*ll.rho(p0_, T0_);

    psiP_ = molr1*psiP_ + molr2*ll.psiP_;
    psiT_ = molr1*psiT_ + molr2*ll.psiT_;
}

inline void linearLiquid::operator-=(const linearLiquid& ll)
{
    scalar molr1 = this->nMoles();

    specie::operator-=(ll);
    molr1 /= this->nMoles();
    scalar molr2 = ll.nMoles()/this->nMoles();

    // p0_ remains unchanged: use first reference p
    
    // Adjust rho0 by aligning the ll rho0_ with current
    rho0_ = molr1*rho0_ - molr2*ll.rho(p0_, T0_);

    psiP_ = molr1*psiP_ - molr2*ll.psiP_;
    psiT_ = molr1*psiT_ - molr2*ll.psiT_;
}

inline void linearLiquid::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline linearLiquid operator+
(
    const linearLiquid& ll1,
    const linearLiquid& ll2
)
{
    scalar nMoles = ll1.nMoles() + ll2.nMoles();
    scalar molr1 = ll1.nMoles()/nMoles;
    scalar molr2 = ll2.nMoles()/nMoles;

    return linearLiquid
    (
        static_cast<const specie&>(ll1)
      + static_cast<const specie&>(ll2),
        molr1*ll1.rho0_ + molr2*ll2.rho(ll1.p0_, ll1.T0_),
        ll1.p0_,
        molr1*ll1.psiP_ + molr2*ll2.psiP_,
        ll1.T0_,
        molr1*ll1.psiT_ + molr2*ll2.psiT_
    );
}


inline linearLiquid operator-
(
    const linearLiquid& ll1,
    const linearLiquid& ll2
)
{
    scalar nMoles = ll1.nMoles() - ll2.nMoles();
    scalar molr1 = ll1.nMoles()/nMoles;
    scalar molr2 = ll2.nMoles()/nMoles;

    return linearLiquid
    (
        static_cast<const specie&>(ll1)
      - static_cast<const specie&>(ll2),
        molr1*ll1.rho0_ - molr2*ll2.rho(ll1.p0_, ll1.T0_),
        ll1.p0_,
        molr1*ll1.psiP_ - molr2*ll2.psiP_,
        ll1.T0_,
        molr1*ll1.psiT_ - molr2*ll2.psiT_
    );
}


inline linearLiquid operator*
(
    const scalar s,
    const linearLiquid& ll
)
{
    return linearLiquid
    (
        s*static_cast<const specie&>(ll),
        ll.rho0_,
        ll.p0_,
        ll.psiP_,
        ll.T0_,
        ll.psiT_
    );
}


inline linearLiquid operator==
(
    const linearLiquid& ll1,
    const linearLiquid& ll2
)
{
    return ll2 - ll1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
