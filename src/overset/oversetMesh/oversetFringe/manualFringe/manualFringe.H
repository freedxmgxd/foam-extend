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
    manualFringe

Description
    Manual fringe selection algorithm

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Full rewrite of interfaces, Mar/2026

Contributor
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    manualFringe.C

\*---------------------------------------------------------------------------*/

#ifndef manualFringe_H
#define manualFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class manualFringe Declaration
\*---------------------------------------------------------------------------*/

class manualFringe
:
    public oversetFringe
{
    // Private data

        //- Cell set of deactivated (hole) cell
        word holesSetName_;

        //- Cell set of acceptor cell
        word acceptorsSetName_;

        //- Hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;

        //- Final donor/acceptor pairs for this region (fringe)
        mutable donorAcceptorList* finalDonorAcceptorsPtr_;

        //- Switch whether to update acceptors/donors for possible changes in
        //  cell sets (switched off by default)
        const Switch updateFringe_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        manualFringe(const manualFringe&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const manualFringe&) = delete;


        // Calculate hole and acceptor addressing
        void calcAddressing() const;

        // Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct from dictionary
        manualFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    //- Destructor
    virtual ~manualFringe();


    // Member Functions

        // Access to fringe data

            //- Return list of hole cells
            virtual const labelList& fringeHoles() const;

            //- Return list of acceptor cells
            virtual const labelList& candidateAcceptors() const;

            //- Return list of final donor acceptor pairs
            virtual const donorAcceptorList& finalDonorAcceptors() const;


        // Iteration control

            //- Initialise donor search, before acceptors are distributed
            virtual void initSearch
            (
                const labelList& candidateAcceptors,
                donorAcceptorList& donorAcceptorRegionData
            ) const;

            //- Set donor suitability
            virtual void setDonorSuitability
            (
                donorAcceptorList& donorAcceptorRegionData
            ) const;
    
            //- Update iteration
            virtual bool updateIteration
            (
                donorAcceptorList& donorAcceptorRegionData
            ) const;


        //- Update the fringe
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
