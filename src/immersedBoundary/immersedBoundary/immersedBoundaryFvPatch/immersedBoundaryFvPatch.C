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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "immersedBoundaryFvPatch.H"
#include "emptyFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvPatch, 1);

    addToRunTimeSelectionTable(fvPatch, immersedBoundaryFvPatch, polyPatch);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::nonOrthogonalFactor_
(
    "immersedBoundaryNonOrthogonalFactor",
    0.1
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeCf(slicedSurfaceVectorField& Cf) const
{
    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    // Inserting only local data
    Cf.boundaryField()[index()].UList::shallowCopy
    (
        ibPolyPatch().ibPatch().faceCentres()
    );
}


void Foam::immersedBoundaryFvPatch::makeSf(slicedSurfaceVectorField& Sf) const
{
    // Insert the patch data for the immersed boundary
    // Note: use the corrected face areas from immersed boundary instead of
    // the stand-alone patch areas within the IB
    // HJ, 30/Nov/2017
    // Inserting only local data
    Sf.boundaryField()[index()].UList::shallowCopy
    (
        ibPolyPatch().correctedIbPatchFaceAreas()
    );
}


void Foam::immersedBoundaryFvPatch::makeC(slicedVolVectorField& C) const
{
    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    // Inserting only local data
    C.boundaryField()[index()].UList::shallowCopy
    (
        ibPolyPatch().ibPatch().faceCentres()
    );
}


void Foam::immersedBoundaryFvPatch::makeV(scalarField& V) const
{}


void Foam::immersedBoundaryFvPatch::updatePhi
(
    DimensionedField<scalar, volMesh>& V,
    DimensionedField<scalar, volMesh>& V0,
    surfaceScalarField& phi
) const
{
    // Correct face fluxes for cut area and insert the immersed patch fluxes

    const fvMesh& mesh = boundaryMesh().mesh();

    const polyBoundaryMesh& bm = boundaryMesh().mesh().boundaryMesh();

    scalar deltaT = mesh.time().deltaT().value();
    scalar rDeltaT = 1.0/deltaT;


    // Scaling of internal mesh flux field should be done only for the current
    // ib patch to avoid scaling multiple times in case of multiple Ib patches
    // present. (IG 3/Dec/2018)

    // Scale internalField
    scalarField& phiIn = phi.internalField();

    const labelList& deadFaces = ibPolyPatch_.deadFaces();
    forAll (deadFaces, dfI)
    {
        const label faceI = deadFaces[dfI];
        if (mesh.isInternalFace(faceI))
        {
            phiIn[faceI] = scalar(0);
        }
        else
        {
            // Boundary face
            const label patchID = bm.whichPatch(faceI);

            if (!isA<emptyFvPatch>(boundaryMesh()[patchID]))
            {
                const label faceID = bm[patchID].whichFace(faceI);

                phi.boundaryField()[patchID][faceID] = scalar(0);
            }
        }
    }

    // Multiply the raw mesh motion flux with the masking function

    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();

    const vectorField& faceAreas = mesh.faceAreas();

    const labelList& cutFaces = ibPolyPatch_.ibFaces();
    forAll (cutFaces, cfI)
    {
        const label faceI = cutFaces[cfI];

        const scalar ibAreaRatio =
            mag(faceAreas[faceI])/faces[faceI].mag(points);

        if (mesh.isInternalFace(faceI))
        {
            // Multiply by masking function
            phiIn[faceI] *= ibAreaRatio;
        }
        else
        {
            // Boundary face
            const label patchID = bm.whichPatch(faceI);

            if (!isA<emptyFvPatch>(boundaryMesh()[patchID]))
            {
                const label faceID = bm[patchID].whichFace(faceI);

                phi.boundaryField()[patchID][faceID] *= ibAreaRatio;
            }
        }
    }

    // Immersed boundary patch
    // Calculate the mesh motion flux from the old and new coordinate of
    // triangular face centres and the time step dotted with the new face area
    phi.boundaryField()[index()] =
    (
        ibPolyPatch_.motionDistance()
      & ibPolyPatch_.correctedIbPatchFaceAreas()
    )*rDeltaT;

    // Check and adjust the immersed boundary space conservation law
    // The mesh motion fluxes come from the actual mesh motion or the motion
    // of the immersed boundary
    // The new cell volumes come from the current mesh configuration
    // The space conservation law will be satisfied by adjusting either
    // the old or the new cell volume.  HJ, 15/Dec/2017

    // First sum up all the fluxes
    scalarField divPhi(mesh.nCells(), 0);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll (owner, faceI)
    {
        divPhi[owner[faceI]] += phiIn[faceI];
        divPhi[neighbour[faceI]] -= phiIn[faceI];
    }

    // Add the mesh motion fluxes from all patches including immersed boundary
    forAll (mesh.boundary(), patchI)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        const scalarField& pssf = phi.boundaryField()[patchI];

        // Check for size since uninitialised ib patches can have zero size at
        // this point (IG 7/Nov/2018)
        if (pssf.size() > 0)
        {
            forAll (pFaceCells, faceI)
            {
                divPhi[pFaceCells[faceI]] += pssf[faceI];
            }
        }
    }

    // Use corrected cell volume
    scalarField& newVols = V.field();
    scalarField& oldVols = V0.field();

    // Multiply by the time-step size and add new volume
    scalarField magDivPhi = mag((newVols - oldVols)*rDeltaT - divPhi);

    // Note:
    // The immersed boundary is now in the new position.  Therefore, some
    // cells that were cut are no longer in the contact with the IB, meaning
    // that ALL cells need to be checked and corrected
    // HJ, 22/Dec/2017
    forAll (magDivPhi, cellI)
    {
        // if (magDivPhi[cellI] > SMALL)
        if (magDivPhi[cellI] > 1e-40)
        {
            // Attempt to correct via old volume
            scalar corrOldVol = newVols[cellI] - divPhi[cellI]*deltaT;

            // Pout<< "Flux maneouvre for cell " << cellI << ": "
            //     << " error: " << magDivPhi[cellI]
            //     << " V: " << newVols[cellI]
            //     << " V0: " << oldVols[cellI]
            //     << " divPhi: " << divPhi[cellI];

            if (corrOldVol < SMALL)
            {
                // Update new volume because old volume cannot carry
                // the correction
                newVols[cellI] = oldVols[cellI] + divPhi[cellI]*deltaT;
            }
            else
            {
                oldVols[cellI] = corrOldVol;
            }

            // scalar corrDivMeshPhi =
            //     mag((newVols[cellI] - oldVols[cellI]) - divPhi[cellI]*deltaT);
            // Pout<< " Corrected: " << corrDivMeshPhi << endl;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeltaCoeffs
(
    fvsPatchScalarField& dc
) const
{
    const vectorField d = delta();

    dc = 1.0/max((nf() & d), 0.05*mag(d));
}


void Foam::immersedBoundaryFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
{
    // Set patch non-orthogonality correction to zero on the patch
    cv = vector::zero;

    // Kill correction vectors in dead cells
    // Potential problem: cannot kill correction vectors on coupled boundaries
    // because the are set later.  For the moment, only the internal
    // correction vectors are killed.
    // HJ, 3/May/2022

    vectorField& cvIn = const_cast<vectorField&>(cv.internalField());

    // Get dead faces
    const labelList& deadFaces = ibPolyPatch_.deadFaces();

    const fvMesh& mesh = boundaryMesh().mesh();

    forAll (deadFaces, dfI)
    {
        if (mesh.isInternalFace(deadFaces[dfI]))
        {
            cvIn[deadFaces[dfI]] = vector::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvPatch::immersedBoundaryFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    ibPolyPatch_(refCast<const immersedBoundaryPolyPatch>(patch)),
    mesh_(bm.mesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::immersedBoundaryFvPatch::size() const
{
    // Immersed boundary patch size equals to the number of intersected cells
    // HJ, 28/Nov/2017

    // Note: asking for patch size triggers the cutting which involves
    // parallel communication.  This should be avoided under read/write, ie
    // when the ibPolyPatch_ is not initialised.
    // Initialisation happens when the fvMesh is initialised, which should be
    // sufficient
    //  HJ, 12/Dec/2018
    // if (!ibPolyPatch_.active())
    // {
    //     return 0;
    // }

    return ibPolyPatch_.ibCells().size();
}


const Foam::labelUList&
Foam::immersedBoundaryFvPatch::faceCells() const
{
    return ibPolyPatch_.ibCells();
}


Foam::tmp<Foam::vectorField> Foam::immersedBoundaryFvPatch::nf() const
{
    // The algorithm has been changed because basic IB patch information
    // (nf and delta) is used in assembly of derived information
    // (eg. deltaCoeffs) and circular dependency needs to be avoided.
    // nf and delta vectors shall be calculated directly from the intersected
    // patch.  HJ, 21/Mar/2019

    return ibPolyPatch_.ibPatch().faceNormals();
}


Foam::tmp<Foam::vectorField> Foam::immersedBoundaryFvPatch::delta() const
{
    // Not strictly needed: this is for debug only.  HJ, 5/Apr/2019
    return ibPolyPatch_.ibPatch().faceCentres() - Cn();
}


// ************************************************************************* //
