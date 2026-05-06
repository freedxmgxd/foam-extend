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

#include "overlapFringe.H"
#include "oversetRegion.H"
#include "oversetMesh.H"
#include "processorFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::overlapFringe::initSearch
(
    const labelList& candidateAcceptors,
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // Give all acceptors to suitability to set data
    forAll (donorAcceptorRegionData, aI)
    {
        donorAcceptor& daPair = donorAcceptorRegionData[aI];

        // Check processor ID
        if (daPair.acceptorProcNo() != Pstream::myProcNo())
        {
            FatalErrorInFunction
                << "Acceptor on different processor: this cannot happen: "
                << "acceptorCell = " << daPair.acceptorProcNo()
                << "myProc = " << Pstream::myProcNo()
                << " acceptorProcNo = " << daPair.acceptorProcNo()
                << abort(FatalError);
        }

        daPair.acceptorSuitability() =
            donorSuitability_->value
            (
                candidateAcceptors[daPair.acceptorCell()]
            );
    }
}


void Foam::overlapFringe::setDonorSuitability
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // Give all acceptors to suitability to set data
    forAll (donorAcceptorRegionData, aI)
    {
        donorAcceptor& daPair = donorAcceptorRegionData[aI];

        if (daPair.donorFound())
        {
            // Check processor ID
            if (daPair.donorProcNo() != Pstream::myProcNo())
            {
                FatalErrorInFunction
                    << "Donor on different processor: this cannot happen: "
                    << "donorCell = " << daPair.donorCell()
                    << " myProc = " << Pstream::myProcNo()
                    << " donorProcNo = " << daPair.donorProcNo()
                    << abort(FatalError);
            }

            daPair.donorSuitability() =
                donorSuitability_->value(daPair.donorCell());
        }
    }
}


bool Foam::overlapFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    if (!fringeHolesPtr_ || !acceptorsPtr_)
    {
        FatalErrorInFunction
            << "fringeHolesPtr_ or acceptorsPtr_ is not allocated. "
            << "Make sure you have called acceptors() or fringeHoles() to "
            << "calculate the initial set of donor/acceptors before "
            << "actually updating iteration."
            << abort(FatalError);
    }

    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorInFunction
            << "Called iteration update with finalDonorAcceptorsPtr_ "
            << "allocated. The final overlap has been achieved, "
            << "prohibiting further calls to updateIteration."
            << abort(FatalError);
    }

    // Increment iteration counter for output
    ++fringeIter_;

    // Allocate worker cumulative donor/acceptor list if it has not been
    // allocated yet (first iteration). Use largest possible size to prevent
    // any resizing
    if (!cumulativeDonorAcceptorsPtr_)
    {
        cumulativeDonorAcceptorsPtr_ = new donorAcceptorDynamicList
        (
            region().mesh().nCells()
        );
    }
    donorAcceptorDynamicList& cumDAPairs = *cumulativeDonorAcceptorsPtr_;

    // Create a list containing unsuitable donors
    donorAcceptorDynamicList unsuitableDAPairs(donorAcceptorRegionData.size());

    // Calculate and report donor suitability
    scalarField ds(donorAcceptorRegionData.size());

    forAll (donorAcceptorRegionData, daPairI)
    {
        ds[daPairI] =
            donorSuitability_->suitabilityFraction
            (
                donorAcceptorRegionData[daPairI]
            );
    }

    // Loop through donor/acceptor pairs and perform mark-up
    forAll (donorAcceptorRegionData, daPairI)
    {
        if
        (
            donorSuitability_->isDonorSuitable(donorAcceptorRegionData[daPairI])
        )
        {
            // Donor is suitable, add it directly to the cumulative list
            cumDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
        else
        {
            // Donor is not suitable, append it to the unsuitable list
            unsuitableDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
    }

    // Calculate the number of total suitable pairs found so far and the number
    // of total pairs
    const label nSuitablePairs =
        returnReduce<label>(cumDAPairs.size(), sumOp<label>());

    const label nTotalPairs = nSuitablePairs
      + returnReduce<label>(unsuitableDAPairs.size(), sumOp<label>());

    const scalar suitabilityFrac = scalar(nSuitablePairs)/scalar(nTotalPairs);

    // Print information
    Info<< "Overlap fringe iteration: " << fringeIter_
        << " for region: " << region().name()
        << nl
        << "Cumulative suitable pairs: " << nSuitablePairs
        << ", total number of pairs: " << nTotalPairs
        << " (" << suitabilityFrac*100 << "%)"
        << endl;

    // Check whether the criterion has been satisfied
    if (suitabilityFrac > minGlobalFraction_)
    {
        // Append unsuitable donors to the list as well
        cumDAPairs.append(unsuitableDAPairs);

        // Now that we have reached suitability criterion specified by the user,
        // we need to clean up a bit. Namely, it is possible that a certain
        // acceptor cell is completely surrounded by holes or other acceptor, so
        // this cell needs to become a hole as well. For easier parallel
        // processing, we will create an indicator field where hole and acceptor
        // cells are marked with 1 and all the other cells (live cells) are
        // marked with -1. We will then use this indicator field to determine
        // whether this acceptor needs to become a hole.

        // Get mesh
        const fvMesh& mesh = region().mesh();

        // Create the indicator field to transfer hole cells to the
        // other side
        volLabelField holeIndicator
        (
            IOobject
            (
                "holeIndicator_" + region().name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedLabel("minusOne", dimless, -1)
        );
        labelField& holeIndicatorIn = holeIndicator.internalField();

        // Transfer fringeHolesPtr into the dynamic list for efficiency. Note:
        // will be transfered back at the end of the scope.
        dynamicLabelList allFringeHoles(*fringeHolesPtr_);

        // Loop through all fringe holes and mark them
        forAll (allFringeHoles, hcI)
        {
            holeIndicatorIn[allFringeHoles[hcI]] = 1;
        }

        // Loop through all acceptors and mark them
        forAll (cumDAPairs, daPairI)
        {
            holeIndicatorIn[cumDAPairs[daPairI].acceptorCell()] = 1;
        }

        // Get boundary field
        volLabelField::GeometricBoundaryField& holeIndicatorb =
            holeIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        evaluateNonOversetBoundaries(holeIndicatorb);

        // Get necessary mesh data
        const cellList& meshCells = mesh.cells();
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        // List of acceptors to be converted to holes
        boolList accBecomingHoles(cumDAPairs.size(), false);

        // Loop through all donor/acceptor pairs collected so far
        forAll (cumDAPairs, daPairI)
        {
            // Get acceptor cell index
            const label& accI = cumDAPairs[daPairI].acceptorCell();

            // Get faces of this cell
            const cell& accFaces = meshCells[accI];

            // Create a bool whether this acceptor needs to be converted to hole
            bool convertToHole = true;

            // Loop through faces
            forAll (accFaces, faceI)
            {
                // Get global face index
                const label& gfI = accFaces[faceI];

                // Check whether this is an internal face or patch face
                if (mesh.isInternalFace(gfI))
                {
                    // Internal face, check whether I'm owner or neighbour
                    if (own[gfI] == accI)
                    {
                        // I'm owner, check whether the neighbour is live
                        if (holeIndicatorIn[nei[gfI]] < 0)
                        {
                            // This acceptor has a live cell for neighbour,
                            // update the flag and continue
                            convertToHole = false;

                            // Found live neighbour.  No need for further search
                            continue;
                        }
                    }
                    else
                    {
                        // I'm neighbour, check whether the owner is live
                        if (holeIndicatorIn[own[gfI]] < 0)
                        {
                            // This acceptor has a live cell for neighbour,
                            // update the flag and continue
                            convertToHole = false;

                            // Found live neighbour.  No need for further search
                            continue;
                        }
                    }
                }
                else
                {
                    // Get patch and face index
                    const label patchI = mesh.boundaryMesh().whichPatch(gfI);
                    const label pfI =
                        mesh.boundaryMesh()[patchI].whichFace(gfI);

                    // Only consider processor patches
                    if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
                    {
                        // Note: patch stores neighbour field after evaluation
                        if (holeIndicatorb[patchI][pfI] < 0)
                        {
                            // This acceptor has a live cell for neighbour on
                            // the other processor, update the flag and continue
                            convertToHole = false;

                            // Found live neighbour.  No need for further search
                            continue;
                        }
                    }
                }
            }

            // Mark whether this acceptor cell has to be converted to hole
            accBecomingHoles[daPairI] = convertToHole;
        }

        // Now we need to filter the data: append acceptors that need to be
        // converted to holes into allFringeHoles and insert all other acceptors
        // into finalDAPairs temporary container
        // Create another dynamic list to collect final donor/acceptor pairs
        donorAcceptorDynamicList finalDAPairs(cumDAPairs.size());

        // Count number of acceptor holes that need to be converted to holes
        label nAccToHoles = 0;

        // Loop all current donor/acceptor pairs
        forAll (cumDAPairs, daPairI)
        {
            if (accBecomingHoles[daPairI])
            {
                // Append the acceptor to list of holes
                allFringeHoles.append(cumDAPairs[daPairI].acceptorCell());
                ++nAccToHoles;
            }
            else
            {
                // Append the donor/acceptor pair to finalDAPairs list
                finalDAPairs.append(cumDAPairs[daPairI]);
            }
        }

        // Bugfix: Although we have found suitable overlap, we need to update
        // acceptors as well because eligible donors for acceptors of other
        // regions are calculated based on these acceptors (and holes)
        labelList& acceptors = *acceptorsPtr_;
        acceptors.setSize(finalDAPairs.size());

        forAll (acceptors, aI)
        {
            acceptors[aI] = finalDAPairs[aI].acceptorCell();
        }

        // Transfer final donor/acceptor list
        finalDonorAcceptorsPtr_ = new donorAcceptorList(finalDAPairs);

        // Tranfer back the allFringeHoles dynamic list into member data
        fringeHolesPtr_->transfer(allFringeHoles);

        // At least 100*minGlobalFraction_ % of suitable donor/acceptor pairs
        // have been found.
        Info<< "Converted " << nAccToHoles << " acceptors to holes."
            << nl
            << "Finished assembling overlap fringe. " << endl;

        // Set the flag to true
        updateSuitableOverlapFlag(true);
    }
    else
    {
        // A sufficient number of suitable donor/acceptors has not been
        // found. Go through unsuitable donor/acceptor pairs and find a new
        // batch of acceptors and holes for the next iteration

        // Get necessary mesh data
        const fvMesh& mesh = region().mesh();
        const labelListList& cc = mesh.cellCells();

        // Create the indicator field to transfer the unsuitable
        // acceptors to the other side
        volLabelField cellTypeIndicator
        (
            IOobject
            (
                "cellTypeIndicator",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedLabel("minusOne", dimless, -1)
        );
        labelField& cellTypeIndicatorIn = cellTypeIndicator.internalField();

        // Transfer fringeHolesPtr into the dynamic list for efficiency. Note:
        // will be transfered back at the end of the scope.
        dynamicLabelList cumFringeHoles(*fringeHolesPtr_);

        // Create mask to prevent wrong and duplicate entries (i.e. we cannot
        // search backwards through existing acceptors and holes)
        boolList freeCells(mesh.nCells(), true);

        // Mask all considered suitable acceptor cells so far
        forAll (cumDAPairs, cpI)
        {
            freeCells[cumDAPairs[cpI].acceptorCell()] = false;
        }

        // Mask all current unsuitable acceptor pairs as well
        forAll (unsuitableDAPairs, upI)
        {
            const label& accCellI = unsuitableDAPairs[upI].acceptorCell();

            freeCells[accCellI] = false;

            // Mark unsuitable pair for possible processor transfer
            cellTypeIndicatorIn[accCellI] = 1;
        }

        // Mask all fringe holes
        forAll (cumFringeHoles, cfhI)
        {
            const label& fhCellI = cumFringeHoles[cfhI];

            freeCells[fhCellI] = false;

            // Mark fringe hole for possible processor transfer
            cellTypeIndicatorIn[fhCellI] = 1;
        }

        // Create dynamic list to efficiently append new batch of
        // acceptors. Note: allocate enough storage.
        dynamicLabelList newAcceptors(10*unsuitableDAPairs.size());

        // Loop through unsuitable acceptors
        forAll (unsuitableDAPairs, upI)
        {
            // Get acceptor cell and its neighbours
            const label& accI = unsuitableDAPairs[upI].acceptorCell();
            const labelList& aNbrs = cc[accI];

            // Loop through neighbours of this acceptor cell
            forAll (aNbrs, nbrI)
            {
                // Check whether the neighbouring cell is free
                const label& nbrCellI = aNbrs[nbrI];

                if (freeCells[nbrCellI])
                {
                    // This cell is neither an old acceptor, fringe hole nor it
                    // has been considered previously. Append it to the
                    // newAcceptors list and mark it as visited
                    newAcceptors.append(nbrCellI);
                    freeCells[nbrCellI] = false;
                }
            }

            // Append this "old" acceptor cell into fringe holes list
            cumFringeHoles.append(accI);
        }

        // Transfer the fringe accross processor boundaries

        // Get boundary field
        volLabelField::GeometricBoundaryField& cellTypeIndicatorBf =
            cellTypeIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        evaluateNonOversetBoundaries(cellTypeIndicatorBf);

        // Loop through boundary field
        forAll (cellTypeIndicatorBf, patchI)
        {
            // Get patch field
            const fvPatchLabelField& chipf = cellTypeIndicatorBf[patchI];

            // Only perform acceptor search if this is a processor boundary
            if (isA<processorFvPatchLabelField>(chipf))
            {
                // Get neighbour field
                const labelField nbrProcIndicator = chipf.patchNeighbourField();

                // Get face cells
                const labelUList& fc = chipf.patch().faceCells();

                // Loop through neighbouring processor field
                forAll (nbrProcIndicator, pfaceI)
                {
                    if
                    (
                        nbrProcIndicator[pfaceI] > 0
                     && freeCells[fc[pfaceI]]
                    )
                    {
                        // The cell on the other side is a hole or acceptor,
                        // while the cell on this side has not been marked yet
                        // as an acceptor. Append the cell to new set of
                        // acceptors and mark it as ineligible in order to
                        // propage the fringe on this side
                        newAcceptors.append(fc[pfaceI]);
                        freeCells[fc[pfaceI]] = false;
                    }
                }
            }
        }

        if (returnReduce(newAcceptors.empty(), andOp<bool>()))
        {
            FatalErrorInFunction
                << "Did not find any new candidate acceptors."
                << nl
                << "Please review your overlap fringe assembly settings."
                << abort(FatalError);
        }

        // Transfer back cumulative fringe holes into the fringeHolesPtr_
        fringeHolesPtr_->transfer(cumFringeHoles);

        // Transfer new acceptors into the acceptors list
        acceptorsPtr_->transfer(newAcceptors);

        // Set the flag to false (suitable overlap not found)
        updateSuitableOverlapFlag(false);
    }

    return foundSuitableOverlap();
}


// ************************************************************************* //
