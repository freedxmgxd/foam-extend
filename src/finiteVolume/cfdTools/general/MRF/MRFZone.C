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

\*---------------------------------------------------------------------------*/

#include "MRFZone.H"
#include "cylindricalCS.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "fvc.H"
#include "syncTools.H"
#include "faceSet.H"
#include "geometricOneField.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::MRFZone, 0);


const Foam::debug::tolerancesSwitch Foam::MRFZone::rotAngle_
(
    "MRFZoneRotAngle",
    10
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0: not in zone
    //  1: moving with frame
    //  2: other
    labelList faceType(mesh_.nFaces(), 0);

    // Determine faces in cell zone
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    if (cellZoneID_ != -1)
    {
        const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
        forAll (cellLabels, i)
        {
            zoneCell[cellLabels[i]] = true;
        }
    }

    label nZoneFaces = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (zoneCell[own[faceI]] || zoneCell[nei[faceI]])
        {
            // By default, set type "other" for faces in MRF zone
            faceType[faceI] = 1;
            nZoneFaces++;
        }
    }

    labelHashSet excludedPatches(excludedPatchLabels_);

    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Note: coupled patches should be rotated.
        // However, sometimes it is difficult to create a mesh with the
        // GGI sitting properly on a cylinder or a sphere.  In such cases
        // planar/cylindrical GGIs can be excluded, and cyclic GGIs can
        // be preserved.  I am not happy about this, but if left uncorrected
        // it causes a massive conservation issue on rotating GGI pairs
        // HJ, 16/Feb/2021
        if
        (
            (pp.isWall() && !excludedPatches.found(patchI))
        )
        {
            forAll (pp, i)
            {
                label faceI = pp.start() + i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 1;
                    nZoneFaces++;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll (pp, i)
            {
                label faceI = pp.start() + i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 2;
                    nZoneFaces++;
                }
            }
        }
    }

    // Now we have for faceType:
    //  0   : face not in cellZone
    //  1   : internal face or wall patch or empty patch
    //  2   : coupled patch face or excluded patch face

    // Sort into lists per patch.

    internalFaces_.setSize(mesh_.nFaces());
    label nInternal = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceType[faceI] == 1)
        {
            internalFaces_[nInternal++] = faceI;
        }
    }
    internalFaces_.setSize(nInternal);

    labelList nIncludedFaces(patches.size(), 0);
    labelList nExcludedFaces(patches.size(), 0);

    label faceI;

    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll (pp, patchFaceI)
        {
            faceI = pp.start() + patchFaceI;

            if (faceType[faceI] == 1)
            {
                nIncludedFaces[patchI]++;
            }
            else if (faceType[faceI] == 2)
            {
                nExcludedFaces[patchI]++;
            }
        }
        Info<< "Patch " << patchI << " name " << pp.name()
            << " included = " << nIncludedFaces[patchI]
            << " excluded = " << nExcludedFaces[patchI]
            << endl;
    }

    includedFaces_.setSize(patches.size());
    excludedFaces_.setSize(patches.size());

    forAll (nIncludedFaces, patchI)
    {
        includedFaces_[patchI].setSize(nIncludedFaces[patchI]);
        excludedFaces_[patchI].setSize(nExcludedFaces[patchI]);
    }

    // Reset counters
    nIncludedFaces = 0;
    nExcludedFaces = 0;

    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll (pp, patchFaceI)
        {
            faceI = pp.start() + patchFaceI;

            if (faceType[faceI] == 1)
            {
                includedFaces_[patchI][nIncludedFaces[patchI]++] = patchFaceI;
            }
            else if (faceType[faceI] == 2)
            {
                excludedFaces_[patchI][nExcludedFaces[patchI]++] = patchFaceI;
            }
        }
    }

    if (debug)
    {
        faceSet internalFaces
        (
            mesh_,
            "internalFaces",
            labelHashSet(internalFaces_)
        );
        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;
        internalFaces.write();

        faceSet MRFFaces(mesh_, "includedFaces", 100);

        forAll (includedFaces_, patchI)
        {
            forAll (includedFaces_[patchI], i)
            {
                label patchFaceI = includedFaces_[patchI][i];
                MRFFaces.insert(patches[patchI].start() + patchFaceI);
            }
        }
        Pout<< "Writing " << MRFFaces.size()
            << " patch faces in MRF zone to faceSet "
            << MRFFaces.name() << endl;
        MRFFaces.write();

        faceSet excludedFaces(mesh_, "excludedFaces", 100);
        forAll (excludedFaces_, patchI)
        {
            forAll (excludedFaces_[patchI], i)
            {
                label patchFaceI = excludedFaces_[patchI][i];
                excludedFaces.insert(patches[patchI].start() + patchFaceI);
            }
        }
        Pout<< "Writing " << excludedFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << excludedFaces.name() << endl;
        excludedFaces.write();
    }
}


void Foam::MRFZone::calcMRFMeshVelocity() const
{
    // Calculate mesh velocity from deformed mesh without executing mesh motion
    // HJ, 6/Jun/2017

    if (debug)
    {
        InfoInFunction
            << "Calculating mesh velocity for zone " << name_
            << endl;
    }

    if (MRFMeshVelocityPtr_)
    {
        FatalErrorInFunction
            << "Mesh velocity for zone " << name_
            << " already calculated"
            << abort(FatalError);
    }

    // Create the mesh velocity
    MRFMeshVelocityPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "meshVelocity" + name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity*dimArea, 0)
        );
    surfaceScalarField& meshVel = *MRFMeshVelocityPtr_;

    // Record the time of creation of mesh velocity
    meshVelTime_ = mesh_.time().value();

    // If there is no rotation, return
    if (mag(omega_.value()) < SMALL)
    {
        return;
    }

    // Create moving points mask
    scalarField movingPointsMask(mesh_.allPoints().size(), 0);

    const pointField& p = mesh_.allPoints();
    const faceList& f = mesh_.allFaces();
    const cellList& c = mesh_.cells();

    const labelList& cellAddr = mesh_.cellZones()[cellZoneID_];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    // Create a rotational coordinate system

    // Create direction by transposing axis
    vector dir;

    if (mag(axis_.value().x()) > SMALL || mag(axis_.value().y()) > SMALL)
    {
        dir = vector(axis_.value().y(), axis_.value().x(), axis_.value().z());
    }
    else if (mag(axis_.value().z()) > SMALL)
    {
        dir = vector(axis_.value().x(), axis_.value().z(), axis_.value().y());
    }
    else
    {
        FatalErrorInFunction
            << "Cannot transpose axis: " << axis_.value()
            << " for MRF zone " << name_
            << abort(FatalError);
    }

    cylindricalCS cs
    (
        "cs",
        origin_.value(),
        axis_.value(),
        dir,
        false            // inDegrees, for consistency
    );

    // Calculate time-step to be used in rotation from rotAngle_
    // Consider using ramped omega?  HJ, 6/Jun/2017
    const scalar deltaT =
        rotAngle_()*mathematicalConstant::pi/(180*omega_.value());

    if (debug)
    {
        InfoInFunction
            << "deltaT for zone " << name_ << " = " << deltaT
            << endl;
    }

    // mag(Omega()) loses direction information
    // Bug fix, GC, 23/Aug/2018

    // Calculate new points
    const vectorField newP =
        cs.globalPosition
        (
            cs.localPosition(p)
          + vector(0, mag(Omega())*deltaT*sign(omega_.value()), 0)*
            movingPointsMask
        );

    // Calculate mesh velocity for all moving faces
    label faceI, patchFaceI;

    scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];
        meshVelIn[faceI] = f[faceI].sweptVol(p, newP)/deltaT;
    }

    forAll (mesh_.boundary(), patchI)
    {
        const label patchStart = mesh_.boundaryMesh()[patchI].start();

        // Get reference to patch velocity
        scalarField& patchMeshVel = meshVel.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            patchMeshVel[patchFaceI] =
                f[patchStart + patchFaceI].sweptVol(p, newP)/deltaT;
        }

        // Excluded faces - same as included faces
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            patchMeshVel[patchFaceI] =
                f[patchStart + patchFaceI].sweptVol(p, newP)/deltaT;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone(const fvMesh& mesh, Istream& is)
:
    name_(is),
    mesh_(mesh),
    dict_(is),
    cellZoneID_(mesh_.cellZones().findZoneID(name_)),
    excludedPatchNames_
    (
        dict_.lookupOrDefault("nonRotatingPatches", wordList(0))
    ),
    origin_(dict_.lookup("origin")),
    axis_(dict_.lookup("axis")),
    omega_(dict_.lookup("omega")),
    rampTime_(dict_.lookupOrDefault<scalar>("rampTime", 0)),
    MRFMeshVelocityPtr_(nullptr),
    meshVelTime_(-1)
{
    if (dict_.found("patches"))
    {
        WarningInFunction
            << "Ignoring entry 'patches'\n"
            << "    By default all patches within the rotating region rotate."
            << nl << "    Optionally supply excluded patches using "
            << "'nonRotatingPatches'."
            << endl;
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (mag(axis_.value()) < SMALL)
    {
        FatalErrorInFunction
            << "Axis vector has zero magnitude: " << axis_
            << ".  This is not allowed"
            << abort(FatalError);
    }

    axis_ = axis_/mag(axis_);

    excludedPatchLabels_.setSize(excludedPatchNames_.size());

    forAll (excludedPatchNames_, i)
    {
        excludedPatchLabels_[i] = patches.findPatchID(excludedPatchNames_[i]);

        if (excludedPatchLabels_[i] == -1)
        {
            FatalErrorInFunction
               << "cannot find MRF patch " << excludedPatchNames_[i]
                << exit(FatalError);
        }
    }

    bool cellZoneFound = (cellZoneID_ != -1);
    reduce(cellZoneFound, orOp<bool>());

    if (!cellZoneFound)
    {
        FatalErrorInFunction
            << "cannot find MRF cellZone " << name_
            << exit(FatalError);
    }

    Info<< "Creating MRF for cell zone " << name_ << ". rpm = "
        << 30*omega_.value()/mathematicalConstant::pi
        << endl;

    setMRFFaces();
}


const Foam::surfaceScalarField& Foam::MRFZone::MRFMeshVelocity() const
{
    // Check time
    if
    (
        rampTime_ > SMALL
     && mesh_.time().value() < rampTime_
     && mag(meshVelTime_ - mesh_.time().value()) > SMALL
    )
    {
        // MRF is ramping.  Recalculate mesh velocity
        deleteDemandDrivenData(MRFMeshVelocityPtr_);
    }

    // Calculate mesh velocity
    if (!MRFMeshVelocityPtr_)
    {
        calcMRFMeshVelocity();
    }

    return *MRFMeshVelocityPtr_;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::MRFZone::~MRFZone()
{
    deleteDemandDrivenData(MRFMeshVelocityPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFZone::Omega() const
{
    if (rampTime_ < SMALL)
    {
        return omega_.value()*axis_.value();
    }
    else
    {
        // Ramping
        const scalar t = mesh_.time().value();
        const scalar ramp = sin(2*pi/(4*rampTime_)*Foam::min(rampTime_, t));
        Info<< "MRF ramp: " << ramp << " Omega: " << ramp*omega_.value()
            << endl;

        return ramp*omega_.value()*axis_.value();
    }
}


void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector rotVel = Omega();

    // Set Coriolis force
    forAll (MRFCells, i)
    {
        const label cellI = MRFCells[i];
        Usource[cellI] -= V[cellI]*(rotVel ^ U[cellI]);
    }
}


void Foam::MRFZone::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V().field();

    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector rotVel = Omega();

    const scalarField& rhoIn = rho.internalField();

    forAll (MRFCells, i)
    {
        label cellI = MRFCells[i];
        Usource[cellI] -= V[cellI]*rhoIn[cellI]*(rotVel ^ U[cellI]);
    }
}


void Foam::MRFZone::addOmega(volVectorField& omg) const
{
    const vector rotVel = Omega();

    // Set omega in all cells of the rotating cell zone
    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];

    forAll (MRFCells, i)
    {
        omg[MRFCells[i]] = rotVel;
    }
}


void Foam::MRFZone::relativeVelocity(volVectorField& U) const
{
    // Get mesh data
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];

    forAll (MRFCells, i)
    {
        label cellI = MRFCells[i];
        U[cellI] -= (rotVel ^ (C[cellI] - origin));
    }

    // Get the MRF mesh velocity
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    autoPtr<surfaceScalarField> meshPhi;

    if (mesh_.changing())
    {
        meshPhi.set(fvc::meshPhi(U).ptr());
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }
        
        // Get patch data
        const vectorField& patchC = C.boundaryField()[patchI];

        // Get patch normal
        vectorField n = mesh_.boundary()[patchI].nf();

        // Get patch mesh flux and adjust for mesh motion flux if present
        scalarField patchUn = meshVel.boundaryField()[patchI];

        // Add mesh motion flux if present
        if (meshPhi.valid())
        {
            patchUn += meshPhi().boundaryField()[patchI];
        }

        // Normalise by magSf to get normal velocity
        patchUn /= mesh_.magSf().boundaryField()[patchI];

        // Get patch velocity for adjustment
        vectorField& patchU = U.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            const label patchFaceI = includedFaces_[patchI][i];

            // Set rotation
             vector Up = (rotVel ^ (patchC[patchFaceI] - origin));

             // Adjust normal component for mesh motion
             Up +=
                n[patchFaceI]*
                (
                    patchUn[patchFaceI]
                  - (n[patchFaceI] & patchU[patchFaceI])
                );

             // Subtract rotational velocity
             patchU[patchFaceI] -= Up;
        }

        // Excluded faces: no normal flux adjustment
        forAll (excludedFaces_[patchI], i)
        {
            const label patchFaceI = excludedFaces_[patchI][i];

            patchU[patchFaceI] -= (rotVel ^ (patchC[patchFaceI] - origin));
        }
    }
}


void Foam::MRFZone::absoluteVelocity(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];

    forAll (MRFCells, i)
    {
        label cellI = MRFCells[i];
        U[cellI] += (rotVel ^ (C[cellI] - origin));
    }

    // For precision and consistency, normal component of the
    // absolute velocity must match the mesh flux exactly
    // HJ, 3/Jul/2020

    const surfaceScalarField& meshVel = MRFMeshVelocity();

    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    autoPtr<surfaceScalarField> meshPhi;

    if (mesh_.changing())
    {
        meshPhi.set(fvc::meshPhi(U).ptr());
    }

    // Included faces
    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch data
        const vectorField& patchC = C.boundaryField()[patchI];

        // Get patch normal
        vectorField n = mesh_.boundary()[patchI].nf();

        // Get patch mesh flux and adjust for mesh motion flux if present
        scalarField patchUn = meshVel.boundaryField()[patchI];

        // Add mesh motion flux if present
        if (meshPhi.valid())
        {
            patchUn += meshPhi().boundaryField()[patchI];
        }

        // Normalise by magSf to get normal velocity
        patchUn /= mesh_.magSf().boundaryField()[patchI];

        // Get patch velocity
        vectorField& patchU = U.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            label patchFaceI = includedFaces_[patchI][i];

             // Add rotational velocity
            patchU[patchFaceI] +=
                (rotVel ^ (patchC[patchFaceI] - origin));

             // Adjust normal component for mesh motion
             patchU[patchFaceI] -=
                 n[patchFaceI]*
                 (
                     patchUn[patchFaceI]
                   - (n[patchFaceI] & patchU[patchFaceI])
                 );
        }

        // Excluded faces
        forAll (excludedFaces_[patchI], i)
        {
            label patchFaceI = excludedFaces_[patchI][i];

            patchU[patchFaceI] +=
                (rotVel ^ (patchC[patchFaceI] - origin));
        }
    }
}


void Foam::MRFZone::relativeFlux
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    if (mesh_.changing())
    {
        fvc::makeRelative(phi, U);
    }

    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        const label faceI = internalFaces_[i];

        phiIn[faceI] -= meshVelIn[faceI];
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch MRF flux
        const scalarField& patchU = meshVel.boundaryField()[patchI];

        scalarField& patchPhi = phi.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            const label patchFaceI = includedFaces_[patchI][i];

            // Bugfix, HJ and HN, 3/Jul/2020
            // Note: this should be zero if velocity is correctly adjusted
            // Reconsider
            patchPhi[patchFaceI] -= patchU[patchFaceI];
        }

        // Excluded faces
        forAll (excludedFaces_[patchI], i)
        {
            const label patchFaceI = excludedFaces_[patchI][i];

            patchPhi[patchFaceI] -= patchU[patchFaceI];
        }
    }
}


void Foam::MRFZone::relativeFlux
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
) const
{
    surfaceScalarField faceRho = fvc::interpolate(rho);
    
    relativeFlux(phi, faceRho, U);
}


void Foam::MRFZone::relativeFlux
(
    surfaceScalarField& phi,
    const surfaceScalarField& rho,
    const volVectorField& U
) const
{
    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    if (mesh_.changing())
    {
        fvc::makeRelative(phi, rho, U);
    }

    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        const label faceI = internalFaces_[i];

        phiIn[faceI] -= rho[faceI]*meshVelIn[faceI];
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch rho
        const scalarField& patchRho = rho.boundaryField()[patchI];

        // Get patch MRF flux
        const scalarField& patchU = meshVel.boundaryField()[patchI];

        scalarField& patchPhi = phi.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            const label patchFaceI = includedFaces_[patchI][i];

            // Bugfix, HJ and HN, 3/Jul/2020
            // Note: this should be zero if velocity is correctly adjusted
            // Reconsider
            patchPhi[patchFaceI] -= patchRho[patchFaceI]*patchU[patchFaceI];
        }

        // Excluded faces
        forAll (excludedFaces_[patchI], i)
        {
            const label patchFaceI = excludedFaces_[patchI][i];

            patchPhi[patchFaceI] -= patchRho[patchFaceI]*patchU[patchFaceI];
        }
    }
}


void Foam::MRFZone::absoluteFlux
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
    if (mesh_.changing())
    {
        Info<< "Adjust for absolute mesh motion flux" << endl;
        fvc::makeAbsolute(phi, U);
    }

    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    label faceI, patchFaceI;

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];

        phiIn[faceI] += meshVelIn[faceI];
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch MRF flux
        const scalarField& patchU = meshVel.boundaryField()[patchI];

        scalarField& patchPhi = phi.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchU[patchFaceI];
        }

        // Excluded patches
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchU[patchFaceI];
        }
    }
}


void Foam::MRFZone::absoluteFlux
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
) const
{
    surfaceScalarField faceRho = fvc::interpolate(rho);
    
    absoluteFlux(phi, faceRho, U);
}


void Foam::MRFZone::absoluteFlux
(
    surfaceScalarField& phi,
    const surfaceScalarField& rho,
    const volVectorField& U
) const
{
    if (mesh_.changing())
    {
        Info<< "Adjust for absolute mesh motion flux" << endl;
        fvc::makeAbsolute(phi, rho, U);
    }

    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    label faceI, patchFaceI;

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];

        phiIn[faceI] += rho[faceI]*meshVelIn[faceI];
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch rho
        const scalarField& patchRho = rho.boundaryField()[patchI];

        // Get patch MRF flux
        const scalarField& patchU = meshVel.boundaryField()[patchI];

        scalarField& patchPhi = phi.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchRho[patchFaceI]*patchU[patchFaceI];
        }

        // Excluded patches
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchRho[patchFaceI]*patchU[patchFaceI];
        }
    }
}


void Foam::MRFZone::meshPhi
(
    surfaceScalarField& phi
) const
{
    // Account for background mesh motion
    if (mesh_.changing())
    {
        Info<< "Get mesh motion flux" << endl;
        phi = mesh_.phi();
    }
    else
    {
        phi == dimensionedScalar("zero", phi.dimensions(), scalar(0));
    }

    // Get mesh motion velocity
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        const label faceI = internalFaces_[i];

        phiIn[faceI] += meshVelIn[faceI];
    }

    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Get patch MRF flux
        const scalarField& patchU = meshVel.boundaryField()[patchI];

        // Get patch flux
        scalarField& patchPhi = phi.boundaryField()[patchI];

        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            const label patchFaceI = includedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchU[patchFaceI];
        }

        // Excluded faces
        forAll (excludedFaces_[patchI], i)
        {
            const label patchFaceI = excludedFaces_[patchI][i];

            patchPhi[patchFaceI] += patchU[patchFaceI];
        }
    }
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    // Get mesh motion velocity
    const surfaceScalarField& meshVel = MRFMeshVelocity();

    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    autoPtr<surfaceScalarField> meshPhi;

    if (mesh_.changing())
    {
        meshPhi.set(fvc::meshPhi(U).ptr());
    }

    // Included faces
    forAll (mesh_.boundary(), patchI)
    {
        if (includedFaces_[patchI].empty() && excludedFaces_[patchI].empty())
        {
            // Nothing to do
            continue;
        }

        // Correct velocity for non-coupled patches
        if (!mesh_.boundary()[patchI].coupled())
        {
            // Get patch data
            const vectorField& patchC = C.boundaryField()[patchI];

            // Get patch normal
            vectorField n = mesh_.boundary()[patchI].nf();

            // Get patch MRF flux and adjust for mesh motion flux if present
            scalarField patchUn = meshVel.boundaryField()[patchI];

            // Add mesh motion flux if present
            if (meshPhi.valid())
            {
                patchUn += meshPhi().boundaryField()[patchI];
            }

            // Normalise by magSf to get normal velocity
            patchUn /= mesh_.magSf().boundaryField()[patchI];

            // Copy original patch velocity
            vectorField patchU(U.boundaryField()[patchI]);

            // Get parallel velocity for moving meshes
            vectorField patchUp(patchU.size(), vector::zero);

            // If the mesh is moving, calculate patch-parallel velocity
            // from mesh motion for included faces.
            // For excluded faces, velocity will be as before
            if (mesh_.changing())
            {
                const polyPatch& pp = mesh_.boundaryMesh()[patchI];
                const pointField& oldPoints = mesh_.oldPoints();

                const vectorField& patchC = C.boundaryField()[patchI];

                // For included faces, calculate tangential velocity from
                // mesh motion.  See movingWallVelocityFvPatchVectorField
                forAll (includedFaces_[patchI], i)
                {
                    const label patchFaceI = includedFaces_[patchI][i];

                    // Calculate old face centre from old point positions
                    vector oldFc = pp[patchFaceI].centre(oldPoints);

                    // Get wall-parallel mesh motion velocity from geometry
                    patchUp[patchFaceI] =
                        (patchC[patchFaceI] - oldFc)/
                        mesh_.time().deltaT().value();
                }
            }

            // Adjust patch U from rotation
            forAll (includedFaces_[patchI], i)
            {
                const label patchFaceI = includedFaces_[patchI][i];

                // Set rotation
                patchU[patchFaceI] =
                    patchUp[patchFaceI]
                  + (rotVel ^ (patchC[patchFaceI] - origin));

                // Adjust normal component for mesh motion
                patchU[patchFaceI] +=
                    n[patchFaceI]*
                    (
                        patchUn[patchFaceI]
                      - (n[patchFaceI] & patchU[patchFaceI])
                    );
            }

            // Set final U onto patch
            U.boundaryField()[patchI] == patchU;
        }
    }
}


void Foam::MRFZone::Su
(
    const volScalarField& phi,
    const volVectorField& gradPhi,
    volScalarField& source
) const
{
    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    forAll (MRFCells, i)
    {
        source[MRFCells[i]] =
            (rotVel ^ (C[MRFCells[i]] - origin)) & gradPhi[MRFCells[i]];
    }
}


void Foam::MRFZone::Su
(
    const volVectorField& phi,
    const volTensorField& gradPhi,
    volVectorField& source
) const
{
    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    forAll (MRFCells, i)
    {
        source[MRFCells[i]] =
            ((rotVel ^ (C[MRFCells[i]] - origin)) & gradPhi[MRFCells[i]])
          - (rotVel ^ phi[MRFCells[i]]);
    }
}


void Foam::MRFZone::calcMagUTheta
(
    const volVectorField& U,
    volScalarField& cumulativeUTheta
)
{
    volScalarField includedExcludedMask
    (
        IOobject
        (
            "includedExcludedMask",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Mark internal cells
    const labelList& MRFCells = mesh_.cellZones()[cellZoneID_];

    forAll (MRFCells, i)
    {
        const label cellI = MRFCells[i];
        includedExcludedMask[cellI] = 1.0;
    }

    forAll (mesh_.boundary(), patchI)
    {
        // Included faces
        forAll (includedFaces_[patchI], i)
        {
            const label patchFaceI = includedFaces_[patchI][i];

            includedExcludedMask.boundaryField()[patchI][patchFaceI] = 1.0;
        }

        // Excluded faces - same as included faces
        forAll (excludedFaces_[patchI], i)
        {
            const label patchFaceI = excludedFaces_[patchI][i];

            includedExcludedMask.boundaryField()[patchI][patchFaceI] = 1.0;
        }
    }

    vector dir;

    if (mag(axis_.value().x()) > SMALL || mag(axis_.value().y()) > SMALL)
    {
        dir = vector(axis_.value().y(), axis_.value().x(), axis_.value().z());
    }
    else if (mag(axis_.value().z()) > SMALL)
    {
        dir = vector(axis_.value().x(), axis_.value().z(), axis_.value().y());
    }
    else
    {
        FatalErrorInFunction
            << "Cannot transpose axis: " << axis_.value()
            << " for MRF zone " << name_
            << abort(FatalError);
    }

    cylindricalCS cs
    (
        "cs",
        origin_.value(),
        axis_.value(),
        dir,
        false
    );

    volVectorField ULocalCS("ULocalCS", U);

    const vectorField newPin =
        cs.globalPosition
        (
            cs.localPosition(ULocalCS.internalField())
        );

    ULocalCS.internalField() = newPin;

    forAll (ULocalCS.boundaryField(), patchI)
    {
        const vectorField newPp =
            cs.globalPosition
            (
                cs.localPosition(ULocalCS.boundaryField()[patchI])
            );

        ULocalCS.boundaryField()[patchI] = newPp;
    }

    // Vector::Y is tangential component of velocity in cylindrical CS
    cumulativeUTheta ==
        cumulativeUTheta + ULocalCS.component(vector::Y)*includedExcludedMask;
}


// ************************************************************************* //
