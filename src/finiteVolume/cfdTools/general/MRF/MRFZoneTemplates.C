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
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::MRFZone::relativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    // If the mesh is changing, additional mesh motion flux needs to be added
    // HJ, 14/May/2025
    if (mesh_.changing())
    {
        Info<< "Adjust for relative mesh motion flux" << endl;
        fvc::makeRelative(rho, phi);
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


template<class RhoFieldType>
void Foam::MRFZone::absoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (mesh_.changing())
    {
        Info<< "Adjust for absolute mesh motion flux" << endl;
        fvc::makeAbsolute(rho, phi);
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


// ************************************************************************* //
