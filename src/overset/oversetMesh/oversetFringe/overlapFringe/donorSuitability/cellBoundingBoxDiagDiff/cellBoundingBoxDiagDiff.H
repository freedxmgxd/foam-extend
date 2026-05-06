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
    Foam::donorSuitability::cellBoundingBoxDiagDiff

Description
    Cell bounding box diagonal suitability. Donor suitability function
    corresponds to the diagonal of a bounding box encompassing a cell.

    Minimises the overlap by shifting it towards cells with similar bounding box
    diagonal (example: a tet cell with vertices on a hex cell faces (opposite
    faces in a hex, have the same bounding box). This might be suitable for
    unstructured tet meshes.

SourceFiles
    cellBoundingBoxDiagDiff.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef cellBoundingBoxDiagDiff_H
#define cellBoundingBoxDiagDiff_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                    Class cellBoundingBoxDiagDiff Declaration
\*---------------------------------------------------------------------------*/

class cellBoundingBoxDiagDiff
:
    public donorSuitability
{
private:

    // Private data

        //- Threshold criterion
        scalar threshold_;

        //- Local suitability criterion: min face area of cell
        scalarField cellBBDiag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellBoundingBoxDiagDiff
        (
            const cellBoundingBoxDiagDiff&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const cellBoundingBoxDiagDiff&) = delete;


public:

    //- Runtime type information
    TypeName("cellBoundingBoxDiagDiff");


    // Constructors

        //- Construct from components
        cellBoundingBoxDiagDiff
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellBoundingBoxDiagDiff() = default;


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

} // End namespace cellBoundingBoxDiagDiff
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
