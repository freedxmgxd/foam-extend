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
    Foam::donorSuitability::minFaceAreaDiff

Description
    Face area suitability. Donor suitability function corresponds to
    difference in a minimum face area of donor and acceptor cell

    Minimises the overlap by shifting it towards cells with similar face sizes.

SourceFiles
    minFaceAreaDiff.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef minFaceAreaDiff_H
#define minFaceAreaDiff_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                        Class minFaceAreaDiff Declaration
\*---------------------------------------------------------------------------*/

class minFaceAreaDiff
:
    public donorSuitability
{
private:
    
    // Private data

        //- Threshold criterion
        scalar threshold_;

        //- Local suitability criterion: min face area of cell
        scalarField minFaceArea_;

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        minFaceAreaDiff(const minFaceAreaDiff&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const minFaceAreaDiff&) = delete;


public:

    //- Runtime type information
    TypeName("minFaceAreaDiff");


    // Constructors

        //- Construct from components
        minFaceAreaDiff
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~minFaceAreaDiff() = default;


    // Member Functions

        // Access

            //- Return threshold
            inline const scalar& threshold() const
            {
                return threshold_;
            }


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

} // End namespace minFaceAreaDiff
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
