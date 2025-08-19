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
    Foam::incompressiblePorousBafflePressureFvPatchScalarField

Description
        This boundary condition provides a jump condition, using the \cyclic
    condition as a base.
    The porous baffle introduces a pressure jump defined by:
        \f[
            \Delta p = -(A*|U| + B*|U|^2)
        \f]
    where
    \vartable
        p      | pressure [Pa]
        D      | Darcy coefficient
        I      | inertial coefficient
        Umax   | max normal velocity. Used in stabilisation
    \endvartable
    \heading Patch usage
    \table
        Property      | Description              | Required    | Default value
        patchType     | underlying patch type should be cyclic | yes |
        A             | Linear coefficient       | yes         |
        B             | inertial coefficient     | yes         |
        Umax          | max normal velocity      | yes         |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            incompressiblePorousBafflePressure;
        patchType       cyclic;
        jump            uniform 0;
        A               1000000;
        B               0.001;
        Umax            30;
        value           uniform 0;
    }
    \endverbatim

Note
     The underlying \c patchType should be set to \c cyclic

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Marko Horvat, Wikki Ltd.
    Porous baffle behaviour has been used from the code originally written
    by Sergio Ferraris.  Full re-implementation was necessary

SourceFiles
    incompressiblePorousBafflePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef incompressiblePorousBafflePressureFvPatchScalarField_H
#define incompressiblePorousBafflePressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "jumpCyclicFvPatchFields.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class incompressiblePorousBafflePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class incompressiblePorousBafflePressureFvPatchScalarField
:
    public jumpCyclicFvPatchScalarField
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Darcy pressure loss coefficient
        scalar A_;

        //- Inertia pressure lost coefficient
        scalar B_;

        //- Max normal velocity
        scalar Umax_;

        //- Current jump
        mutable scalarField jump_;

        //- Current time index used to recalculate jump_
        mutable label curTimeIndex_;


    // Private member functions

        //- Calculate jump
        void calcJump() const;


public:

    //- Runtime type information
    TypeName("incompressiblePorousBafflePressure");


    // Constructors

        //- Construct from patch and internal field
        incompressiblePorousBafflePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        incompressiblePorousBafflePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  incompressiblePorousBafflePressureFvPatchScalarField
        //  onto a new patch
        incompressiblePorousBafflePressureFvPatchScalarField
        (
            const incompressiblePorousBafflePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        incompressiblePorousBafflePressureFvPatchScalarField
        (
            const incompressiblePorousBafflePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new incompressiblePorousBafflePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        incompressiblePorousBafflePressureFvPatchScalarField
        (
            const incompressiblePorousBafflePressureFvPatchScalarField&,
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
                new incompressiblePorousBafflePressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~incompressiblePorousBafflePressureFvPatchScalarField() = default;


    // Member functions

        // Access

            //- Return Darcy pressure loss coefficient
            const scalar& A() const
            {
                return A_;
            }

            //- Return inertial pressure lost coefficient
            const scalar& B() const
            {
                return  B_;
            }


        //- Return the "jump" across the patch as a "half" field
        virtual tmp<scalarField> jump() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
