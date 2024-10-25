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
    chtRcEnthalpyDiffusivityFvPatchScalarField

Description
    Diffusivity for the region couple patch field when solving for hs.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    chtRcEnthalpyDiffusivityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcEnthalpyDiffusivityFvPatchScalarField_H
#define chtRcEnthalpyDiffusivityFvPatchScalarField_H

#include "chtRcThermalDiffusivityFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class chtRcEnthalpyDiffusivityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcEnthalpyDiffusivityFvPatchScalarField
:
    public chtRcThermalDiffusivityFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("chtRcEnthalpyDiffusivity");


    // Constructors

        //- Construct from patch and internal field
        chtRcEnthalpyDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcEnthalpyDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  chtRcEnthalpyDiffusivityFvPatchScalarField onto a new patch
        chtRcEnthalpyDiffusivityFvPatchScalarField
        (
            const chtRcEnthalpyDiffusivityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcEnthalpyDiffusivityFvPatchScalarField
        (
            const chtRcEnthalpyDiffusivityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcEnthalpyDiffusivityFvPatchScalarField(*this)
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
                new chtRcEnthalpyDiffusivityFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~chtRcEnthalpyDiffusivityFvPatchScalarField() = default;


    // Member functions

        // Access functions

            //- Return thermal diffusivity for temperature from decoupled state
            virtual tmp<scalarField> k() const;

            //- Return patch Cp
            virtual tmp<scalarField> Cp() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);


        // Patch interpolation and patch flux

            //- Calculate patch face interpolate given weighting factors
            virtual void patchInterpolate
            (
                surfaceScalarField& fField,
                const scalarField& pL
            ) const;

            //- Calculate patch face interpolate given two weighting factors
            virtual void patchInterpolate
            (
                surfaceScalarField& fField,
                const scalarField& pL,
                const scalarField& pY
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
