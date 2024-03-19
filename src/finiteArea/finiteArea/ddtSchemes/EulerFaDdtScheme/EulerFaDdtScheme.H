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
    EulerDdtScheme

Description
    Basic first-order Euler implicit/explicit ddt using only the current and
    previous time-step values.

SourceFiles
    EulerDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerFaDdtScheme_H
#define EulerFaDdtScheme_H

#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class EulerFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerFaDdtScheme
:
    public fa::faDdtScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        EulerFaDdtScheme(const EulerFaDdtScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const EulerFaDdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from mesh
        EulerFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        EulerFaDdtScheme(const faMesh& mesh, Istream& is)
        :
            faDdtScheme<Type>(mesh, is)
        {}


    //- Destructor
    virtual ~EulerFaDdtScheme() = default;


    // Member Functions

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensioned<Type>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensioned<Type>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faePatchField, edgeMesh> > facDdt0
        (
            const GeometricField<Type, faePatchField, edgeMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<faMatrix<Type> > famDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<faMatrix<Type> > famDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        virtual tmp<faMatrix<Type> > famDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EulerFaDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
