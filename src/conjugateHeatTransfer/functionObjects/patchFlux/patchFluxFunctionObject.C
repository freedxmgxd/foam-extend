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
#include "patchFluxFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchFluxFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchFluxFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchFluxFunctionObject::patchFluxFunctionObject
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
    patchName_(dict.lookup("patch")),
    phiName_(dict.lookupOrDefault<word>("phiName", "phi")),
    ofPtr_()
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating patchFluxFunctionObject for patch "
        << patchName_  << " on region " << regionName_ << endl;

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
        ofPtr_() << "# Time, patch flux" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchFluxFunctionObject::start()
{
    return true;
}


bool Foam::patchFluxFunctionObject::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Get patch index
    label patchID = mesh.boundaryMesh().findPatchID(patchName_);

    // Check whether the patch has been found
    if (patchID < 0)
    {
        InfoInFunction
            << "Patch named: " << patchName_ << " not found."
            << " Returning."
            << endl;

        return false;
    }

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        // Calculate the flux through the patch
        const scalar patchFlux = gSum(phi.boundaryField()[patchID]);

        if (ofPtr_.valid())
        {
            ofPtr_() << time_.value() << tab << patchFlux << endl;
        }

        Info<< "Flux " << phiName_ << " through " << patchName_
            << " region " << regionName_
            << ": " << patchFlux << endl;

        return true;
    }
    else
    {
        InfoInFunction
            << "Flux " << phiName_ << " for region "
            << " region " << regionName_ << " not found. Returning."
            << endl;

        return false;
    }
}


bool Foam::patchFluxFunctionObject::read(const dictionary& dict)
{
    patchName_ = word(dict.lookup("patch"));
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");

    return false;
}

// ************************************************************************* //
