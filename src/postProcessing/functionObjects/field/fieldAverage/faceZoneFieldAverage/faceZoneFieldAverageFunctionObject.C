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
    Marko Horvat, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "surfaceFields.H"
#include "faceZoneFieldAverageFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneFieldAverageFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        faceZoneFieldAverageFunctionObject,
        dictionary
    );
}

// * * * * * * * * * * * * * Private member functions* * * * * * * * * * * * //

template<class Type>
const Type
Foam::faceZoneFieldAverageFunctionObject::calcAverage
(
    const faceZone& zone,
    const GeometricField<Type, fvPatchField, volMesh> volField
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfaceField;

    const surfaceField surfField = fvc::interpolate(volField);

    const Field<Type>& fieldIn = surfField.internalField();

    const typename surfaceField::GeometricBoundaryField& fieldBndry =
        surfField.boundaryField();

    const fvMesh& mesh = volField.mesh();

    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    Type faceZoneFieldAverage = pTraits<Type>::zero;

    scalar area = 0;

    // Get flip map from the face zone
    const boolList& faceFlips = zone.flipMap();

    forAll (zone, i)
    {
        const label faceI = zone[i];
        const bool faceFlipI = faceFlips[i];
        scalar curArea = 0;

        if (mesh.isInternalFace(faceI))
        {
            curArea = magSfIn[faceI];

            if (faceFlipI)
            {
                faceZoneFieldAverage -= curArea*fieldIn[faceI];
            }
            else
            {
                faceZoneFieldAverage += curArea*fieldIn[faceI];
            }
        }
        else
        {
            const label patchI = mesh.boundaryMesh().whichPatch(faceI);

            if (patchI < 0)
            {
                // Face in zone is beyond live faces list
                continue;
            }

            const label patchFaceI =
                mesh.boundaryMesh()[patchI].whichFace(faceI);

            curArea = magSf.boundaryField()[patchI][patchFaceI];

            if (faceFlipI)
            {
                faceZoneFieldAverage -= curArea*fieldBndry[patchI][patchFaceI];
            }
            else
            {
                faceZoneFieldAverage += curArea*fieldBndry[patchI][patchFaceI];
            }
        }

        area += curArea;
    }

    // Reduce for parallel runs
    reduce(faceZoneFieldAverage, sumOp<Type>());
    reduce(area, sumOp<scalar>());

    if (area < SMALL)
    {
        WarningInFunction
            << "Face zone " << zone.name() << " has zero area" << endl;
    }
    else
    {
        faceZoneFieldAverage /= area;
    }

    return faceZoneFieldAverage;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneFieldAverageFunctionObject::
faceZoneFieldAverageFunctionObject
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
    fieldName_(dict.lookup("fieldName")),
    ofPtr_()
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating faceZoneFieldAverageFunctionObject for face zone "
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
        ofPtr_() << "# Time, face zone average" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceZoneFieldAverageFunctionObject::start()
{
    return true;
}


bool Foam::faceZoneFieldAverageFunctionObject::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const label faceZoneID = mesh.faceZones().findZoneID(faceZoneName_);
    if (faceZoneID < 0)
    {
        InfoInFunction
            << "Face zone: " << faceZoneName_ << " not found."
            << " Returning."
            << endl;
        return false;
    }

    const faceZone& zone = mesh.faceZones()[faceZoneID];

    if (mesh.foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& volField =
            mesh.lookupObject<volScalarField>(fieldName_);
        scalar faceZoneFieldAverage = calcAverage(zone, volField);
        if (ofPtr_.valid())
        {
            ofPtr_() << time_.value() << tab << faceZoneFieldAverage << endl;
        }
        Info<< "Average of the field " << fieldName_ << " on the zone "
            << faceZoneName_ << ": " << faceZoneFieldAverage << endl;

        return true;
    }
    else if (mesh.foundObject<vectorField>(fieldName_))
    {
        const volVectorField& volField =
            mesh.lookupObject<volVectorField>(fieldName_);

        vector faceZoneFieldAverage = calcAverage(zone, volField);
        if (ofPtr_.valid())
        {
            ofPtr_() << time_.value() << tab << faceZoneFieldAverage << endl;
        }
        Info<< "Average of the field " << fieldName_ << " on the zone "
            << faceZoneName_ << ": " << faceZoneFieldAverage << endl;

        return true;
    }
    else
    {
        InfoInFunction
            << "Field " << fieldName_ << " not found. Returning."
            << endl;
        return false;
    }

}


bool Foam::faceZoneFieldAverageFunctionObject::read(const dictionary& dict)
{
    faceZoneName_ = word(dict.lookup("faceZone"));
    fieldName_ = word(dict.lookup("fieldName"));

    return false;
}





// ************************************************************************* //
