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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "faceZoneHeatFluxFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneHeatFluxFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        faceZoneHeatFluxFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneHeatFluxFunctionObject::
faceZoneHeatFluxFunctionObject
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
    TName_(),
    phiName_(),
    gammaName_(),
    ofPtr_()
{
    // Read parameters
    read(dict);

    Info<< "Creating faceZoneHeatFluxFunctionObject for face zone"
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
        ofPtr_()
            << "# Time, face zone convective, diffusive, total flux"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceZoneHeatFluxFunctionObject::start()
{
    return true;
}


bool Foam::faceZoneHeatFluxFunctionObject::execute(const bool forceWrite)
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

    // Get the face zone
    const faceZone& zone = mesh.faceZones()[faceZoneID];

    // Get flip map from the face zone
    const boolList& faceFlips = zone.flipMap();

    // Loop through zone faces and calculate the flux
    scalar faceZoneConvectiveHeatFlux = 0;
    scalar faceZoneDiffusiveHeatFlux = 0;

    if (mesh.foundObject<volScalarField>(TName_))
    {
        const surfaceScalarField& magSf = mesh.magSf();
        const scalarField& magSfIn = magSf.internalField();
        const surfaceScalarField::GeometricBoundaryField& magSfBndry =
            magSf.boundaryField();

        // Get flux field and separately internal and boundary fields
        const volScalarField& T =
            mesh.lookupObject<volScalarField>(TName_);

        if (mesh.foundObject<surfaceScalarField>(phiName_))
        {
            const surfaceScalarField& phi =
                mesh.lookupObject<surfaceScalarField>(phiName_);

            // Dimension check and correction
            if (debug)
            {
                if (phi.dimensions() == dimVolume/dimTime)
                {
                    Info<< "Flux dimension = m^3/s" << endl;
                }
                else if (phi.dimensions() == dimMass/dimTime)
                {
                    Info<< "Flux dimension = kg/s" << endl;
                }
                else if
                (
                    phi.dimensions() == dimSpecificHeatCapacity*dimMass/dimTime
                )
                {
                    Info<< "Flux dimension = W/K" << endl;
                }
            }

            const scalarField& phiIn = phi.internalField();
            const surfaceScalarField::GeometricBoundaryField& phiBndry =
                phi.boundaryField();

            surfaceScalarField faceT = fvc::interpolate(T);
            const scalarField& faceTIn = faceT.internalField();
            const surfaceScalarField::GeometricBoundaryField& faceTBndry =
                faceT.boundaryField();

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
                        faceZoneConvectiveHeatFlux -=
                            phiIn[faceI]*faceTIn[faceI];
                    }
                    else
                    {
                        faceZoneConvectiveHeatFlux +=
                            phiIn[faceI]*faceTIn[faceI];
                    }
                }
                else
                {
                    // This face is on the boundary, get patch and face index
                    const label patchI = mesh.boundaryMesh().whichPatch(faceI);
                    const label patchFaceI =
                        mesh.boundaryMesh()[patchI].whichFace(faceI);

                    // Note:
                    // Processor patches not handled correctly
                    // HJ, 12/Dec/2022

                    if (patchI < 0)
                    {
                        FatalErrorInFunction
                            << "Cannot find patch for boundary face: " << faceI
                                << abort(FatalError);
                    }

                    if (faceFlipI)
                    {
                        // Face needs to be flipped (sign needs to be reverted)
                        faceZoneConvectiveHeatFlux -=
                            phiBndry[patchI][patchFaceI]*
                            faceTBndry[patchI][patchFaceI];
                    }
                    else
                    {
                        faceZoneConvectiveHeatFlux +=
                            phiBndry[patchI][patchFaceI]*
                            faceTBndry[patchI][patchFaceI];
                    }
                }
            }

            // Reduce for parallel runs
            reduce(faceZoneConvectiveHeatFlux, sumOp<scalar>());
        }

        if (mesh.foundObject<volScalarField>(gammaName_))
        {
            const volScalarField& gamma =
                mesh.lookupObject<volScalarField>(gammaName_);

            // Dimension check and correction
            if (debug)
            {
                if (gamma.dimensions() == sqr(dimLength)/dimTime)
                {
                    Info<< "Flux dimension = m^3/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                {
                    Info<< "Flux dimension = kg/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                {
                    Info<< "Flux dimension = W/K" << endl;
                }
            }

            surfaceScalarField faceGamma = fvc::interpolate(gamma);
            const scalarField& faceGammaIn = faceGamma.internalField();
            const surfaceScalarField::GeometricBoundaryField& faceGammaBndry =
                faceGamma.boundaryField();

            surfaceScalarField faceSnGradT = fvc::snGrad(T);
            const scalarField& faceSnGradTIn = faceSnGradT.internalField();
            const surfaceScalarField::GeometricBoundaryField& faceSnGradTBndry =
                faceSnGradT.boundaryField();

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
                        faceZoneDiffusiveHeatFlux -=
                            faceGammaIn[faceI]*faceSnGradTIn[faceI]*
                            magSfIn[faceI];
                    }
                    else
                    {
                        faceZoneDiffusiveHeatFlux +=
                            faceGammaIn[faceI]*faceSnGradTIn[faceI]*
                            magSfIn[faceI];
                    }
                }
                else
                {
                    // This face is on the boundary, get patch and face index
                    const label patchI = mesh.boundaryMesh().whichPatch(faceI);
                    const label patchFaceI =
                        mesh.boundaryMesh()[patchI].whichFace(faceI);

                    // Note:
                    // Processor patches not handled correctly
                    // HJ, 12/Dec/2022

                    if (patchI < 0)
                    {
                        FatalErrorInFunction
                            << "Cannot find patch for boundary face: " << faceI
                                << abort(FatalError);
                    }

                    if (faceFlipI)
                    {
                        // Face needs to be flipped (sign needs to be reverted)
                        faceZoneDiffusiveHeatFlux -=
                            faceGammaBndry[patchI][patchFaceI]*
                            faceSnGradTBndry[patchI][patchFaceI]*
                            magSfBndry[patchI][patchFaceI];
                    }
                    else
                    {
                        faceZoneDiffusiveHeatFlux +=
                            faceGammaBndry[patchI][patchFaceI]*
                            faceSnGradTBndry[patchI][patchFaceI]*
                            magSfBndry[patchI][patchFaceI];
                    }
                }
            }

            // Reduce for parallel runs
            reduce(faceZoneDiffusiveHeatFlux, sumOp<scalar>());
        }

        if (ofPtr_.valid())
        {
            ofPtr_()
                << time_.value() << tab
                << faceZoneConvectiveHeatFlux << tab
                << -faceZoneDiffusiveHeatFlux << tab
                << faceZoneConvectiveHeatFlux - faceZoneDiffusiveHeatFlux
                << endl;
        }

        Info<< "Flux " << phiName_ << " through " << faceZoneName_
            << " (convective, diffusive, total): "
            << faceZoneConvectiveHeatFlux << tab
            << -faceZoneDiffusiveHeatFlux << tab
            << faceZoneConvectiveHeatFlux - faceZoneDiffusiveHeatFlux
            << endl;

        return true;
    }
    else
    {
        InfoInFunction
            << "Temperature field " << TName_ << " not found. Returning."
            << endl;

        return false;
    }
}


bool Foam::faceZoneHeatFluxFunctionObject::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    faceZoneName_ = word(dict.lookup("faceZone"));

    TName_ = word(dict.lookup("T"));
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    gammaName_ = word(dict.lookup("gamma"));

    return false;
}

// ************************************************************************* //
