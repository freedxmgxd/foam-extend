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
    Foam::fv::steadyStateDdtScheme

Description
    SteadyState implicit/explicit ddt which returns 0.

SourceFiles
    steadyStateDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef steadyStateDdtScheme_H
#define steadyStateDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class steadyStateDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class steadyStateDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        steadyStateDdtScheme(const steadyStateDdtScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const steadyStateDdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("steadyState");


    // Constructors

        //- Construct from mesh
        steadyStateDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        steadyStateDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {}


    //- Destructor
    virtual ~steadyStateDdtScheme() = default;


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
tmp<surfaceScalarField> steadyStateDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> steadyStateDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> steadyStateDdtScheme<scalar>::fvcDdtConsistentPhiCorr
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
#   include "steadyStateDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
