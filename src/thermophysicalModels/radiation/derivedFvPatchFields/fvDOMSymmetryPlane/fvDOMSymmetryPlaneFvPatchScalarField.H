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
    Foam::fvDOMSymmetryPlaneFvPatchScalarField

Description
    Specular and symmetry plane radiation boundary for the fvDOM model.

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    fvDOMSymmetryPlaneFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fvDOMSymmetryPlaneFvPatchScalarField_H
#define fvDOMSymmetryPlaneFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
           Class fvDOMSymmetryPlaneFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fvDOMSymmetryPlaneFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of temperature field
        word TName_;


public:

    //- Runtime type information
    TypeName("fvDOMSymmetryPlane");


    // Constructors

        //- Construct from patch and internal field
        fvDOMSymmetryPlaneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fvDOMSymmetryPlaneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  fvDOMSymmetryPlaneFvPatchScalarField onto a new patch
        fvDOMSymmetryPlaneFvPatchScalarField
        (
            const fvDOMSymmetryPlaneFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvDOMSymmetryPlaneFvPatchScalarField
        (
            const fvDOMSymmetryPlaneFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fvDOMSymmetryPlaneFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fvDOMSymmetryPlaneFvPatchScalarField
        (
            const fvDOMSymmetryPlaneFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fvDOMSymmetryPlaneFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~fvDOMSymmetryPlaneFvPatchScalarField() = default;


    // Member functions

        // Access

            //- Return the temperature field name
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //  adjustment
            word& TName()
            {
                return TName_;
            }

    
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
