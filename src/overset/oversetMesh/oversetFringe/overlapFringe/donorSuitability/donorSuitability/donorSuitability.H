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
    Foam::donorSuitability::donorSuitability

Description
    An abstract base class for donor suitability algorithms. This class is
    used by oversetFringe algorithm to optimise overlap given donor suitability
    function.

    Donor suitability criteria is given in isDonorSuitable member function.

SourceFiles
    donorSuitability.C
    newAcceptorSet.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

Author, complete rewrite
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Josip Katavic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef donorSuitability_H
#define donorSuitability_H

#include "dictionary.H"
#include "labelList.H"
#include "donorAcceptorList.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class oversetFringe;

namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                         Class donorSuitability Declaration
\*---------------------------------------------------------------------------*/

class donorSuitability
{
    // Private data

        //- Reference to oversetFringe algorithm
        const oversetFringe& oversetFringe_;

        //- Subdictionary for acceptor set specification
        const dictionary coeffDict_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        donorSuitability(const donorSuitability&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const donorSuitability&) = delete;


public:

    //- Runtime type information
    TypeName("donorSuitability");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            donorSuitability,
            dictionary,
            (
                const oversetFringe& oversetFringeAlgorithm,
                const dictionary& dict
            ),
            (oversetFringeAlgorithm, dict)
        );


    // Constructors

        //- Construct from components
        donorSuitability
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    // Selectors

        //- Return an autoPtr to the selected donorSuitability
        static autoPtr<donorSuitability> New
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~donorSuitability() = default;


    // Member Functions

        // Access

            //- Return a const reference to oversetFringe algorithm
            inline const oversetFringe& oversetFringeAlgorithm() const
            {
                return oversetFringe_;
            }

            //- Return coefficient dictionary
            inline const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


        // Overlap minimisation via Donor Suitability Function

            //- Return suitability parameter
            //  Meaning depends on suitability function, eg. cell volume
            virtual scalar value(const label& cellID) const = 0;

            //- Return donor/acceptor pair suitability fraction between 0 and 1
            //  Higher value means better suitability
            virtual scalar suitabilityFraction
            (
                const donorAcceptor& daPair
            ) const = 0;

            //- Is this donor/acceptor pair suitable?
            virtual bool isDonorSuitable
            (
                const donorAcceptor& daPair
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace donorSuitability
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
