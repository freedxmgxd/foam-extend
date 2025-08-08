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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "surfaceFields.H"
#include "faceZoneFluxFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneFluxFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        faceZoneFluxFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneFluxFunctionObject::
faceZoneFluxFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    faceZoneName_(dict.lookup("faceZone")),
    phiName_(dict.lookupOrDefault<word>("phiName", "phi")),
    ofPtr_()
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating faceZoneFluxFunctionObject for face zone"
        << faceZoneName_  << " on region " << regionName_ << endl;

    // Set stream pointer: only master writes
    if (Pstream::master())
    {
        // Write into case directory instead of processor directory
        if (Pstream::parRun())
        {
            ofPtr_.reset
            (
                new OFstream(time_.path()/".."/word(dict.lookup("file")))
            );
        }
        else
        {
            ofPtr_.reset
            (
                new OFstream(time_.path()/word(dict.lookup("file")))
            );
        }

        // Write header
        ofPtr_() << "# Time, face zone flux" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceZoneFluxFunctionObject::start()
{
    return true;
}


bool Foam::faceZoneFluxFunctionObject::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Find the face zone
    const label faceZoneID = mesh.faceZones().findZoneID(faceZoneName_);

    // Check the zone has been found
    if (faceZoneID < 0)
    {
        InfoInFunction
            << "Face zone: " << faceZoneName_ << " not found."
            << " Returning."
            << endl;

        return false;
    }

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        // Get flux field and separately internal and boundary fields
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        const scalarField& phiIn = phi.internalField();
        const surfaceScalarField::GeometricBoundaryField& phiBndry =
            phi.boundaryField();

        // Get the face zone
        const faceZone& zone = mesh.faceZones()[faceZoneID];

        // Get flip map from the face zone
        const boolList& faceFlips = zone.flipMap();

        // Loop through zone faces and calculate the total flux
        scalar faceZoneFlux = 0;
        forAll (zone, i)
        {
            // Get face index and its flip
            const label faceI = zone[i];
            const bool faceFlipI = faceFlips[i];

            if (mesh.isInternalFace(faceI))
            {
                if (faceFlipI)
                {
                    // Face needs to be flipped (sign needs to be reverted)
                    faceZoneFlux -= phiIn[faceI];
                }
                else
                {
                    faceZoneFlux += phiIn[faceI];
                }
            }
            else
            {
                // This face is on the boundary, get patch and patch face index
                const label patchI = mesh.boundaryMesh().whichPatch(faceI);
                const label patchFaceI =
                    mesh.boundaryMesh()[patchI].whichFace(faceI);

                if (patchI < 0)
                {
                    FatalErrorInFunction
                        << "Cannot find patch for boundary face: " << faceI
                        << abort(FatalError);
                }

                if (faceFlipI)
                {
                    // Face needs to be flipped (sign needs to be reverted)
                    faceZoneFlux -= phiBndry[patchI][patchFaceI];
                }
                else
                {
                    faceZoneFlux += phiBndry[patchI][patchFaceI];
                }
            }
        }

        // Reduce for parallel runs
        reduce(faceZoneFlux, sumOp<scalar>());

        if (ofPtr_.valid())
        {
            ofPtr_() << time_.value() << tab << faceZoneFlux << endl;
        }

        Info<< "Flux through " << faceZoneName_ << ": " << faceZoneFlux << endl;

        return true;
    }
    else
    {
        InfoInFunction
            << "Flux field " << phiName_ << " not found. Returning."
            << endl;

        return false;
    }

}


bool Foam::faceZoneFluxFunctionObject::read(const dictionary& dict)
{
    faceZoneName_ = word(dict.lookup("faceZone"));
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");

    return false;
}

// ************************************************************************* //
