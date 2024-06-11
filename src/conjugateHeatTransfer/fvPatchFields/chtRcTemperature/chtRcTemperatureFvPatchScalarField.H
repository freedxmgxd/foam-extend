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
    chtRcTemperatureFvPatchScalarField

Description
    Region couple patch field (solving for T), taking radiation into account
    Complete rewrite and functionality reimplementation of garbage provided
    by Henrik Rusche.  That was shameful shit with unbalanced fluxes across
    the interface, indicating total ignorance of CFD.  Programming-wise,
    the code was wrong, illegible and uncommented.

    The code uses region coupling interface jump conditions to account
    for radiative heat flux off the CHT interface.
    It achieves perfect energy balance on the interface with and without
    radiation
    HJ, 9/Jun/2024

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    chtRcTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcTemperatureFvPatchScalarField_H
#define chtRcTemperatureFvPatchScalarField_H

#include "jumpRegionCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class chtRcTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcTemperatureFvPatchScalarField
:
    public jumpRegionCouplingFvPatchScalarField
{
    // Private data

        //- Name of thermal diffusivity field
        word kName_;

        //- Take radiation into account for this patch
        Switch radiation_;

        //- Radiation jump
        scalarField jump_;


public:

    //- Runtime type information
    TypeName("chtRcTemperature");


    // Constructors

        //- Construct from patch and internal field
        chtRcTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRcTemperatureFvPatchScalarField
        //  onto a new patch
        chtRcTemperatureFvPatchScalarField
        (
            const chtRcTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcTemperatureFvPatchScalarField
        (
            const chtRcTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcTemperatureFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~chtRcTemperatureFvPatchScalarField() = default;


    // Member functions

        // Access functions

            //- Return kName
            const word& kName() const
            {
                return kName_;
            }

            //- Return access to kName
            word& kName()
            {
                return kName_;
            }

            //- Return radiation
            const Switch& radiation() const
            {
                return radiation_;
            }

            //- Return access to radiation
            Switch& radiation()
            {
                return radiation_;
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<scalarField> jump() const;

            //- Return shadow patch field
            virtual const chtRcTemperatureFvPatchScalarField&
                shadowPatchField() const;

            //- Return patch internal temperature
            virtual tmp<scalarField> patchInternalT() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

        // Evaluation functions

            //- Return neighbour field of internal field across coupled patch
            virtual tmp<scalarField> patchNeighbourField() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
