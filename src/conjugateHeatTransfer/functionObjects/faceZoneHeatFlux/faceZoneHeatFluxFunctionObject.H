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
    faceZoneHeatFluxFunctionObject

Description
    FunctionObject calculating flux through a face zone.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    faceZoneHeatFluxFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneHeatFluxFunctionObject_H
#define faceZoneHeatFluxFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class faceZoneHeatFluxFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class faceZoneHeatFluxFunctionObject
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

        //- Name of the face zone
        word faceZoneName_;

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
        faceZoneHeatFluxFunctionObject
        (
            const faceZoneHeatFluxFunctionObject&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const faceZoneHeatFluxFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("faceZoneHeatFlux");


    // Constructors

        //- Construct from components
        faceZoneHeatFluxFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    ~faceZoneHeatFluxFunctionObject() = default;


    // Member Functions

        //- start is called at the start of the time-loop
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
