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
    Foam::donorSuitability::patchDistance

Description
    Patch distance suitability.  Acceptor should be as close to a given set of
    patches as possible.

SourceFiles
    patchDistance.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef patchDistance_H
#define patchDistance_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                       Class patchDistance Declaration
\*---------------------------------------------------------------------------*/

class patchDistance
:
    public donorSuitability
{
private:

    // Private data

        //- Local suitability criterion: min distance to patch
        scalarField magDistance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchDistance(const patchDistance&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const patchDistance&) = delete;


public:

    //- Runtime type information
    TypeName("patchDistance");


    // Constructors

        //- Construct from components
        patchDistance
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchDistance() = default;


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

} // End namespace patchDistance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
