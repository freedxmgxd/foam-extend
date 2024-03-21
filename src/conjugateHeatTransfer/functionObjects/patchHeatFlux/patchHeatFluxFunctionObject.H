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
    patchHeatFluxFunctionObject

Description
    Function object calculating convective-diffusive flux through a
    specified patch.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    patchHeatFluxFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef patchHeatFluxFunctionObject_H
#define patchHeatFluxFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class patchHeatFluxFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchHeatFluxFunctionObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Name of the patch (if only one patch is used)
        word patchName_;

        //- Name of the variable
        word TName_;

        //- Name of the flux field (phi by default)
        word phiName_;

        //- Name of the diffusivity
        word gammaName_;

        //- Output file
        autoPtr<OFstream> ofPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchHeatFluxFunctionObject
        (
            const patchHeatFluxFunctionObject&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const patchHeatFluxFunctionObject&) = delete;

        //- Get indices of all uncoupled patches
        labelList patchIDList() const;


public:

    //- Runtime type information
    TypeName("patchHeatFlux");


    // Constructors

        //- Construct from components
        patchHeatFluxFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    ~patchHeatFluxFunctionObject() = default;


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField& points)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
