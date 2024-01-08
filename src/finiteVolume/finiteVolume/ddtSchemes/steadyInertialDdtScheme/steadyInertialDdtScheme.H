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
    Foam::fv::steadyInertialDdtScheme

Description
    Stabilised local time-step first-order Euler implicit/explicit ddt
    for intertially relaxed steady-state simulations.
    The time-step is adjusted locally so that an advective equations remains
    diagonally dominant.

    This scheme should only be used for steady-state computations
    using transient codes where local time-stepping is preferably to
    under-relaxation for transport consistency reasons.

    See also CoEulerDdtScheme.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    steadyInertialDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef steadyInertialDdtScheme_H
#define steadyInertialDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class steadyInertialDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class steadyInertialDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Name of the flux field used to calculate the local
        //  convective time-step
        word phiName_;

        //- Name of the diffusivity field used to calculate the local
        //  diffusive time-step
        word nuName_;

        //- Name of the density field used to obtain the volumetric flux
        //  from the mass flux if required
        word rhoName_;

        //- Maximum local Courant number
        scalar maxCo_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        steadyInertialDdtScheme(const steadyInertialDdtScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const steadyInertialDdtScheme&) = delete;

        //- Return the reciprocal of the Courant-number limited time-step
        tmp<volScalarField> CorDeltaT() const;

        //- Return the reciprocal of the face-Courant-number limited time-step
        //  for convective transport
        tmp<surfaceScalarField> convectionCofrDeltaT() const;

        //- Return the reciprocal of the face-Courant-number limited time-step
        //  for diffusive transport
        tmp<surfaceScalarField> diffusionCofrDeltaT() const;

        //- Return the reciprocal of the face-Courant-number limited time-step
        //  for diffusive transport ginev face diffusivity
        tmp<surfaceScalarField> diffusionCofrDeltaT
        (
            const surfaceScalarField& nuf
        ) const;


public:

    //- Runtime type information
    TypeName("steadyInertial");


    // Constructors

        //- Construct from mesh and Istream
        steadyInertialDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            phiName_(is),
            nuName_(is),
            rhoName_(is),
            maxCo_(readScalar(is))
        {}


    //- Destructor
    virtual ~steadyInertialDdtScheme() = default;


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
tmp<surfaceScalarField> steadyInertialDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> steadyInertialDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> steadyInertialDdtScheme<scalar>::fvcDdtConsistentPhiCorr
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
#   include "steadyInertialDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
