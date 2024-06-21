/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    chtRcEnthalpyFvPatchScalarField

Description
    Region couple patch field (solving for enthalpy),
    taking radiation into account.
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
    chtRcEnthalpyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcEnthalpyFvPatchScalarField_H
#define chtRcEnthalpyFvPatchScalarField_H

#include "chtRcTemperatureFvPatchScalarField.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class chtRcEnthalpyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcEnthalpyFvPatchScalarField
:
    public chtRcTemperatureFvPatchScalarField
{
    // Private data

        //- Linearisation for thermo: CpBar
        //  This is temporary storage used in matrix update
        mutable scalarField* CpBarPtr_;

        //- Linearisation for thermo: H0
        //  This is temporary storage used in matrix update
        mutable scalarField* H0Ptr_;


    // Private Member Functions

        //- Clear storage
        void clearOut();

        //- Check Cp and H0
        void checkCpH0() const;


public:

    //- Runtime type information
    TypeName("chtRcEnthalpy");


    // Constructors

        //- Construct from patch and internal field
        chtRcEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRcEnthalpyFvPatchScalarField
        //  onto a new patch
        chtRcEnthalpyFvPatchScalarField
        (
            const chtRcEnthalpyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcEnthalpyFvPatchScalarField
        (
            const chtRcEnthalpyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct as copy setting internal field reference
        chtRcEnthalpyFvPatchScalarField
        (
            const chtRcEnthalpyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcEnthalpyFvPatchScalarField(*this)
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
                new chtRcEnthalpyFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~chtRcEnthalpyFvPatchScalarField();


    // Member functions

        // Access functions

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

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Initialise neighbour matrix update
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

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
