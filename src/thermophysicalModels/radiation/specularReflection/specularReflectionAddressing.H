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
    Foam::specularReflectionAddressing

Description
    Specular reflection addressing mesh object for fvDOM radiation model.

    On specular patches, incoming and outgoing rays need to be matched,
    depending on the angle of the reflected ray and best match with the
    available incoming ray.  This depends on the resolution (number of rays)
    of the fvDOM model.

    Addressing is held as a list of cartesian2DStencil for each face of the
    patch and created patch-by-patch on demand.  The P ray is the closest match
    of the reflected ray from available incoming rays.  N-S-E-W is the stencil
    of neighbouring rays, provided they are incoming.

    Specular boundary condition can either use the best-match direction
    (without interpolation), or perform its own interpolation as needed.

    Addressing is calculated pathc-by-patch on demand.

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    specularReflectionAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef specularReflectionAddressing_H
#define specularReflectionAddressing_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "fvDOM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                Class specularReflectionAddressing Declaration
\*---------------------------------------------------------------------------*/

class specularReflectionAddressing
:
    public MeshObject<fvMesh, specularReflectionAddressing>
{
    // Private data

        //- Flat patch enumeration
        enum flatPatch
        {
            UNDEFINED = 0,
            FLAT = 1,
            NOT_FLAT = 2
        };

        //- Reflection addressing for patches.
        //  For each incoming ray, provides outgoing (reflected) ray index
        //  If ray is outgoing, index will be -1
        //  First index = patch
        //  Second index = ray
        mutable labelListList* reflectionAddrPtr_;


    // Private member functions

        //- Check for flat patch
        flatPatch checkFlatPatch(const label& patchI) const;

        //- Check if addressing is allocated for patch
        bool checkPatchAddressing(const label& patchI) const;

        //- Construct reflection addressing for patch
        void makePatchAddressing(const label& patchI) const;

        //- Clear out
        void clearOut() const;


public:

    // Declare name of the class and its debug switch
    TypeName("specularReflectionAddressing");


    // Static data members


    // Constructors

        //- Construct given an fvMesh
        explicit specularReflectionAddressing(const fvMesh&);


    //- Destructor
    virtual ~specularReflectionAddressing();


    // Member functions

        //- Return reference to owner least square vectors
        label patchReflectionAddr
        (
            const label& patchI,
            const label& rayId
        ) const;


        //- Update after mesh motion: delete addressing
        //  Delete
        virtual bool movePoints() const;

        //- Update after topo change: delete addressing
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
