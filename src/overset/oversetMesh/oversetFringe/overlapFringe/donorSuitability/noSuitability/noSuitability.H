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
    Foam::donorSuitability::noSuitability

Description
    No suitability check.

SourceFiles
    noSuitability.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef noSuitability_H
#define noSuitability_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                       Class noSuitability Declaration
\*---------------------------------------------------------------------------*/

class noSuitability
:
    public donorSuitability
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        noSuitability(const noSuitability&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const noSuitability&) = delete;


public:

    //- Runtime type information
    TypeName("noSuitability");


    // Constructors

        //- Construct from components
        noSuitability
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        ):
            donorSuitability(oversetFringeAlgorithm, dict)
        {}


    //- Destructor
    virtual ~noSuitability() = default;


    // Member Functions

        // Overlap minimisation virtual functions

            //- Return suitability parameter
            virtual scalar value(const label& cellID) const;

            //- Return suitability
            virtual scalar suitabilityFraction
            (
                const donorAcceptor& daPair
            ) const;

            //- Is this donor/acceptor pair suitable?
            virtual bool isDonorSuitable
            (
                const donorAcceptor& daPair
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noSuitability
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
