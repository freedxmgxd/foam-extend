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
    chtRcThermalDiffusivityFvPatchScalarField

Description
    Diffusivity for the region couple patch field when solving for T.

Author
    Code removed in entirety: Henrik Rusche, Wikki GmbH.
    This was a disastrous heap of junk
    Written anew by Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    chtRcThermalDiffusivityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcThermalDiffusivityFvPatchScalarField_H
#define chtRcThermalDiffusivityFvPatchScalarField_H

#include "regionCouplingFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class chtRcThermalDiffusivityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcThermalDiffusivityFvPatchScalarField
:
    public regionCouplingFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("chtRcThermalDiffusivity");


    // Constructors

        //- Construct from patch and internal field
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  chtRcThermalDiffusivityFvPatchScalarField onto a new patch
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const chtRcThermalDiffusivityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const chtRcThermalDiffusivityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcThermalDiffusivityFvPatchScalarField(*this)
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
                new chtRcThermalDiffusivityFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~chtRcThermalDiffusivityFvPatchScalarField() = default;


    // Member functions

        // Access functions

            //- Remove thermal diffusivity for temperature from decoupled state
            tmp<scalarField> k() const;


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
