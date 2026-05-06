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
    Foam::donorSuitability::cellVolumeDiff

Description
    Cell volumes suitability. Donor suitability function is the minimum
    difference in cell volume between donor and acceptor cell

    Minimises the overlap by shifting it towards cells with similar volume.

SourceFiles
    cellVolumeDiff.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellVolumeDiff_H
#define cellVolumeDiff_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                       Class cellVolumeDiff Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeDiff
:
    public donorSuitability
{
private:

    // Private data

        //- Threshold criterion
        scalar threshold_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellVolumeDiff(const cellVolumeDiff&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const cellVolumeDiff&) = delete;


public:

    //- Runtime type information
    TypeName("cellVolumeDiff");


    // Constructors

        //- Construct from components
        cellVolumeDiff
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellVolumeDiff() = default;


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

} // End namespace cellVolumeDiff
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
