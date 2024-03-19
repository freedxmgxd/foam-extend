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
    Foam::fv::CoEulerDdtScheme

Description
    Courant number limited first-order Euler implicit/explicit ddt.

    The time-step is adjusted locally so that the local Courant number
    does not exceed the specified value.

    This is a transient first-order scheme with "consistent locally introduced
    inertial under-relaxation" for cells where the specified Co limit has been
    breached.  HJ, 22/Dec/2010

SourceFiles
    CoEulerDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CoEulerDdtScheme_H
#define CoEulerDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class CoEulerDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CoEulerDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Name of the flux field used to calculate the local time-step
        word phiName_;

        //- Name of the density field used to obtain the volumetric flux
        //  from the mass flux if required
        word rhoName_;

        //- Maximum local Courant number
        scalar maxCo_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CoEulerDdtScheme(const CoEulerDdtScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const CoEulerDdtScheme&) = delete;

        //- Return the reciprocal of the Courant-number limited time-step
        tmp<volScalarField> CorDeltaT() const;

        //- Return the reciprocal of the face-Courant-number limited time-step
        tmp<surfaceScalarField> CofrDeltaT() const;


public:

    //- Runtime type information
    TypeName("CoEuler");


    // Constructors

        //- Construct from mesh and Istream
        CoEulerDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            phiName_(is),
            rhoName_(is),
            maxCo_(readScalar(is))
        {}


    //- Destructor
    virtual ~CoEulerDdtScheme() = default;


    // Member Functions

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensioned<Type>&
        ) const;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<fvMatrix<Type> > fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<fvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<fvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        ) const;

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        ) const;


        // Member functions for the new time consistent formulation
        virtual tmp<fluxFieldType> fvcDdtConsistentPhiCorr
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const surfaceScalarField& rAUf
        ) const;

        //- Time scheme consistemt mesh flux
        virtual tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


template<>
tmp<surfaceScalarField> CoEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> CoEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> CoEulerDdtScheme<scalar>::fvcDdtConsistentPhiCorr
(
    const surfaceScalarField& faceU,
    const volScalarField& U,
    const surfaceScalarField& rAUf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CoEulerDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
