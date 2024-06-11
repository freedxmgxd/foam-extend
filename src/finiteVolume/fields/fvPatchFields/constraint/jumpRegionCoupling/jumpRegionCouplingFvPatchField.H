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
    jumpRegionCouplingFvPatchField

Description
    Region coupling interface patch field with a jump condition

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    jumpRegionCouplingFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef jumpRegionCouplingFvPatchField_H
#define jumpRegionCouplingFvPatchField_H

#include "regionCouplingFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class jumpRegionCouplingFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jumpRegionCouplingFvPatchField
:
    public regionCouplingFvPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("jumpRegionCoupling");


    // Constructors

        //- Construct from patch and internal field
        jumpRegionCouplingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        jumpRegionCouplingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given jumpRegionCouplingFvPatchField
        //  onto a new patch
        jumpRegionCouplingFvPatchField
        (
            const jumpRegionCouplingFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        jumpRegionCouplingFvPatchField
        (
            const jumpRegionCouplingFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );


    // Member functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return regionCouplingFvPatchField<Type>::type();
            }

            //- Return the "jump" across the patch
            virtual tmp<Field<Type> > jump() const = 0;


        // Evaluation functions

            //- Return neighbour field given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Initialise neighbour matrix update
            //  Accounts for jump
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            // Update result field based on interface functionality
            // Note: Since the jump is added in initInterfaceMatrixUpdate
            // there is no need to over-ride updateInterfaceMatrix
            // HJ, 10/Jun/2024
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "jumpRegionCouplingFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
