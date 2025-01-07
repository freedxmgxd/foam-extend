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

#include "specularReflectionAddressing.H"
#include "mapPolyMesh.H"
#include "radiationModel.H"
#include "fvDOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(specularReflectionAddressing, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::specularReflectionAddressing::flatPatch
Foam::specularReflectionAddressing::checkFlatPatch(const label& patchI) const
{
    const vectorField nf = mesh().boundary()[patchI].nf();

    scalar minDotProd = 1;

    if (!nf.empty())
    {
        // Check the dot product of the first nf agains the rest.
        // If smaller than the tolerance, declare non-flat
        minDotProd = min(nf[0] & nf);
    }

    reduce(minDotProd, minOp<scalar>());

    // Hard-coded flatness tolerance
    if (minDotProd < 0.95)
    {
        return NOT_FLAT;
    }
    else
    {
        return FLAT;
    }
}


void Foam::specularReflectionAddressing::clearOut() const
{
    deleteDemandDrivenData(reflectionAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::specularReflectionAddressing::specularReflectionAddressing
(
    const fvMesh& mesh
)
:
    MeshObject<fvMesh, specularReflectionAddressing>(mesh),
    reflectionAddrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::specularReflectionAddressing::~specularReflectionAddressing()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::specularReflectionAddressing::checkPatchAddressing
(
    const label& patchI
) const
{
    if (!reflectionAddrPtr_)
    {
        // If addressing is not allocated, return false
        return false;
    }
    else if
    (
        reflectionAddrPtr_->operator[](patchI).empty()
     && mesh().boundary()[patchI].size() != 0
    )
    {
        // Check if current patch is not empty and data is calculated
        return false;
    }
    else
    {
        // Already calculated
        return true;
    }
}


void Foam::specularReflectionAddressing::makePatchAddressing
(
    const label& patchI
) const
{
    if (debug)
    {
        InfoInFunction
            << "Constructing specular addressing for patch "
            << mesh().boundary()[patchI].name()
            << endl;
    }

    // If patch is empty, nothing to do
    if (mesh().boundary()[patchI].size() == 0)
    {
        return;
    }

    // Check radiation model
    if (!mesh().foundObject<radiation::radiationModel>("radiationProperties"))
    {
        FatalErrorInFunction
            << "Cannot find radiation model"
            << abort(FatalError);
    }

    // Access radiation model
    const radiation::radiationModel& radiation =
        db().lookupObject<radiation::radiationModel>("radiationProperties");

    // Check if the radiation model is fvDOM
    if (!isA<radiation::fvDOM>(radiation))
    {
        FatalErrorInFunction
            << "Radiation model is not fvDOM.  Cannot use specular reflection"
            << abort(FatalError);
    }

    // Create reference to fvDOM model
    const radiation::fvDOM& dom =
        dynamic_cast<const radiation::fvDOM&>(radiation);

    // Check if reflection is calculated
    if (!reflectionAddrPtr_)
    {
        reflectionAddrPtr_ =
            new labelListList(mesh().boundary().size());
    }

    if (!reflectionAddrPtr_->operator[](patchI).empty())
    {
        FatalErrorInFunction
            << "specular addressing for patch "
            << mesh().boundary()[patchI].name()
            << " already calculated."
            << abort(FatalError);
    }

    if (mesh().boundary()[patchI].coupled())
    {
        FatalErrorInFunction
            << "Requested specular addressing for coupled patch "
            << mesh().boundary()[patchI].name()
            << ".  This is not allowed: regular patches only"
            << abort(FatalError);
    }

    const fvPatch& p = mesh().boundary()[patchI];

    // Create addressing for patch
    labelList& curReflection = reflectionAddrPtr_->operator[](patchI);

    curReflection.setSize(dom.nRay(), -1);

    // Check if the patch is flat
    if (checkPatchAddressing(patchI) != FLAT)
    {
        FatalErrorInFunction
            << "Patch " << patchI << " named " << p.name()
            << " is not flat.  This is currently not supported"
            << abort(FatalError);
    }

    // Get face normals
    const vectorField nHat = p.nf();

    const vector patchNormal = nHat[0];

    // Go through all ray directions
    for (label rayI = 0; rayI < dom.nRay(); rayI++)
    {
        // Normalise ray vector
        vector rayNormal = dom.IRay(rayI).dAve();
        rayNormal /= mag(rayNormal);

        // Check if the ray is outward pointing
        if ((patchNormal & rayNormal) > SMALL)
        {
            // Calculate reflected ray
            const vector reflectedNormal =
                rayNormal - 2*patchNormal*(patchNormal & rayNormal);

            // Check if the reflected normal is inward-pointing
            if ((patchNormal & reflectedNormal) > -SMALL)
            {
                FatalErrorInFunction
                    << "Reflected normal is not inward pointing.  "
                    << "ray = " << rayNormal
                    << " n = " << patchNormal
                    << " reflected ray = " << reflectedNormal
                    << abort(FatalError);
            }

            // Check all ray directions to find best reflected direction
            scalar maxDotProduct = 0;
            label bestRay = -1;

            for (label rayJ = 0; rayJ < dom.nRay(); rayJ++)
            {
                // Get direction of reflection candidate
                vector jNormal = dom.IRay(rayJ).dAve();
                jNormal /= mag(jNormal);

                // Check if reflected candidate is inbound
                if ((patchNormal & reflectedNormal) < -SMALL)
                {
                    scalar curDotProduct = (reflectedNormal & jNormal);

                    if (curDotProduct > maxDotProduct)
                    {
                        maxDotProduct = curDotProduct;

                        bestRay = rayJ;
                    }
                }
            }

            // Check if direction has been found
            if (maxDotProduct > SMALL && bestRay > -1)
            {
                // Set reflected direction for ray
                // curReflection[rayI] = bestRay; HJ, HERE!!!
                curReflection[bestRay] = rayI;
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find best reflected ray for ray " << rayI
                    << " ray = " << rayNormal
                    << " n = " << patchNormal
                    << " reflected ray = " << reflectedNormal
                    << abort(FatalError);
            }
        }
    }

    Info<< "Reflection for patch " << patchI << ": "
        << curReflection
        << endl;

    // Check
    // forAll (curReflection, rayI)
    // {
    //     if (curReflection[rayI] > -1)
    //     {
    //         vector myRay = dom.IRay(rayI).dAve();
    //         myRay /= mag(myRay);

    //         // Check that my ray is outgoing
    //         if ((patchNormal & myRay) < -SMALL)
    //         {
    //             vector reflectedRay = dom.IRay(curReflection[rayI]).dAve();
    //             reflectedRay /= mag(reflectedRay);

    //             // Untangle reflected ray
    //             reflectedRay =
    //                 reflectedRay - 2*patchNormal*(patchNormal & reflectedRay);

    //             Info<< "R " << rayI << " ref " << curReflection[rayI]
    //                 << " my ray = " << myRay
    //                 << " corrected reflected = " << reflectedRay
    //                 << " DOT = " << (myRay & reflectedRay)
    //                 << endl;
    //         }
    //     }
    // }
}


Foam::label Foam::specularReflectionAddressing::patchReflectionAddr
(
    const label& patchI,
    const label& rayId
) const
{
    if (!checkPatchAddressing(patchI))
    {
        makePatchAddressing(patchI);
    }

    const labelListList& r = *reflectionAddrPtr_;

    return r[patchI][rayId];
}


bool Foam::specularReflectionAddressing::movePoints() const
{
    if (debug)
    {
        InfoInFunction
            << "Clearing specular addressing data" << endl;
    }

    clearOut();

    return true;
}


bool Foam::specularReflectionAddressing::updateMesh
(
    const mapPolyMesh& mpm
) const
{
    if (debug)
    {
        InfoInFunction
            << "Clearing specular addressing data" << endl;
    }

    clearOut();

    return true;
}


// ************************************************************************* //
