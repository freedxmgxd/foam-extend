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
    Foam::fv::backwardDualDdtScheme

Author
    Oliver Borm
    Aleksandar Jemcov

SourceFiles
    backwardDualDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef backwardDualDdtScheme_H
#define backwardDualDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class backwardDualDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardDualDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Name of the physical time-step (field)
        word deltaTName_;

        //- Name of the pseudo time-step field
        word deltaTauName_;

        //- Name of old field
        word oldName_;

        //- Current physical time-step
        mutable scalar physDeltaT_;

        //- Previous physical time-step
        mutable scalar physDeltaT0_;

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the current physical time-step
        template<class GeoField>
        scalar deltaT_(const GeoField&) const;

        //- Return the previous physical time-step
        template<class GeoField>
        scalar deltaT0_(const GeoField&, const GeoField&) const;

//         scalar deltaT0_(IOList<label>) const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;


        //- Disallow default bitwise copy construct
        backwardDualDdtScheme(const backwardDualDdtScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const backwardDualDdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("backwardDual");


    // Constructors

        //- Construct from mesh
        backwardDualDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh),
            deltaTName_("physDeltaT"),
            deltaTauName_("CoDeltaT"),
            oldName_("Old"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}

        //- Construct from mesh and Istream
        backwardDualDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            deltaTName_(is),
            deltaTauName_(is),
            oldName_("Old"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}


    //- Destructor
    virtual ~backwardDualDdtScheme() = default;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

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
        tmp<fluxFieldType> fvcDdtConsistentPhiCorr
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const surfaceScalarField& rAUf
        ) const;


        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


template<>
tmp<surfaceScalarField> backwardDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> backwardDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
) const;


template<>
tmp<surfaceScalarField> backwardDualDdtScheme<scalar>::fvcDdtConsistentPhiCorr
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
#   include "backwardDualDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
