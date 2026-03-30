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

#include "oversetRegion.H"
#include "oversetMesh.H"
#include "oversetFringe.H"
#include "polyPatchID.H"
#include "triSurfaceTools.H"
#include "cellSet.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetRegion::calcCutHoleCells() const
{
    if (cutHoleCellsPtr_)
    {
        FatalErrorInFunction
            << "Cut hole cells already calculated"
            << abort(FatalError);
    }

    // Algorithm
    // - go through all regions apart from the current region
    // - get access to hole surface search
    // - mark as hole all cells that fall outside of hole search
    // - combine all searches into a single inside-outside list

    // Get local cell indices
    const labelList& rc = regionCells();

    // Prepare local cell centres for inside-outside search on all regions
    vectorField localC(rc.size());

    const vectorField& c = mesh().cellCentres();

    forAll (localC, i)
    {
        localC[i] = c[rc[i]];
    }

    // Prepare hole mask
    boolList holeMask(mesh().nCells(), false);

    // Mark all hole cells using their hole boundary patch inside search

    // Get regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Go through all regions apart from the current
    forAll (regions, regionI)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const oversetRegion& otherRegion = regions[regionI];

        // If there are no hole patches on other region, skip it
        if (!otherRegion.holePatchesPresent())
        {
            continue;
        }

        // Get reference to hole search
        const triSurfaceSearch& holeSearch = otherRegion.holeSearch();

        boolList regionInside = holeSearch.calcInside(localC);

        // Note: hole mask has the size of all mesh cells and regionInside
        // only of the size of local region
        forAll (regionInside, i)
        {
            holeMask[rc[i]] |= regionInside[i];
        }
    }

    // Count hole cells
    label nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            nHoleCells++;
        }
    }

    // Allocate hole cells storage
    cutHoleCellsPtr_ = new labelList(nHoleCells);
    labelList& ch = *cutHoleCellsPtr_;

    // Reset counter and collect hole cells
    nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            ch[nHoleCells] = rc[i];
            nHoleCells++;
        }
    }

    if (oversetMesh::debug)
    {
        Pout<< "Region " << name()
            << ": number of local holes = " << cutHoleCellsPtr_->size()
            << endl;
    }
}


void Foam::oversetRegion::calcHoleCells() const
{
    if (holeCellsPtr_)
    {
        FatalErrorInFunction
            << "Hole cells already calculated"
            << abort(FatalError);
    }

    // Combine cut hole cells and fringe cells into a single list

    // Prepare hole mask
    boolList holeMask(mesh().nCells(), false);

    // Mask all cut hole cells
    const labelList& cutHoleCells = cutHoles();

    forAll (cutHoleCells, i)
    {
        holeMask[cutHoleCells[i]] = true;
    }

    // Mask fringe hole cells
    const labelList& fringeHoleCells = fringePtr_->fringeHoles();

    forAll (fringeHoleCells, i)
    {
        holeMask[fringeHoleCells[i]] = true;
    }

    // Count hole cells in the region
    const labelList& rc = regionCells();

    label nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            ++nHoleCells;
        }
    }

    // Allocate hole cells storage
    holeCellsPtr_ = new labelList(nHoleCells);
    labelList& h = *holeCellsPtr_;

    // Reset counter and collect hole cells
    nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            h[nHoleCells] = rc[i];
            ++nHoleCells;
        }
    }
}


void Foam::oversetRegion::calcHoleTriMesh() const
{
    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Calculating hole tri mesh for region " << name()
            << endl;
    }

    if (holeTriMeshPtr_)
    {
        FatalErrorInFunction
            << "Hole tri mesh already calculated"
            << abort(FatalError);
    }

    // Create region mask to check if patch touches region
    boolList regionMask(mesh().nCells(), false);

    const labelList& rc = regionCells();

    forAll (rc, rcI)
    {
        regionMask[rc[rcI]] = true;
    }

    // Get hole patch names
    const wordList& holePatchNames = overset().holePatchNames();

    // Collect local hole faces
    labelHashSet holePatches;

    forAll (holePatchNames, nameI)
    {
        polyPatchID curHolePatch
        (
            holePatchNames[nameI],
            mesh().boundaryMesh()
        );

        if (curHolePatch.active())
        {
            // If the patch has zero size, do not insert it
            // Parallel cutting bug.  HJ, 17/Apr/2014
            if (!mesh().boundaryMesh()[curHolePatch.index()].empty())
            {
                // Check if the patch is touching the current region
                const labelList& faceCells =
                    mesh().boundary()[curHolePatch.index()].faceCells();

                label nFound = 0;

                forAll (faceCells, fcI)
                {
                    if (regionMask[faceCells[fcI]])
                    {
                        nFound++;
                    }
                }

                // Check if the complete patch belongs to current region
                if (nFound == faceCells.size())
                {
                    holePatches.insert(curHolePatch.index());
                }
                else if (nFound > 0)
                {
                    WarningInFunction
                        << "Patch " << holePatchNames[nameI]
                        << " seems to be split between multiple regions.  "
                        << "Please check overset region structure.  "
                        << "nFound: " << nFound
                        << " faceCells: " << faceCells.size()
                        << endl;
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Patch "  << holePatchNames[nameI]
                << " cannot be found.  Available patch names: "
                << mesh().boundaryMesh().names()
                << abort(FatalError);
        }
    }

    // Make and invert local triSurface
    triFaceList triFaces;
    pointField triPoints;

    // Memory management
    {
        triSurface ts = triSurfaceTools::triangulate
        (
            mesh().boundaryMesh(),
            holePatches
        );

        // Clean mutiple points and zero-sized triangles
        ts.cleanup(false);

        triFaces.setSize(ts.size());
        triPoints = ts.points();

        forAll (ts, tsI)
        {
            triFaces[tsI] = ts[tsI].reverseFace();
        }
    }

    if (Pstream::parRun())
    {
        // Combine all faces and points into a single list

        // Gather-scatter triPoints
        List<pointField> allTriPoints(Pstream::nProcs());

        {
            allTriPoints[Pstream::myProcNo()] = triPoints;

            Pstream::gatherList(allTriPoints);

            Pstream::scatterList(allTriPoints);

            // Block for all requests and remove storage
            Pstream::waitRequests();
        }


        // Gather-scatter triFaces
        List<triFaceList> allTriFaces(Pstream::nProcs());

        {
            allTriFaces[Pstream::myProcNo()] = triFaces;

            Pstream::gatherList(allTriFaces);
            Pstream::scatterList(allTriFaces);

            // Block for all requests and remove storage
            Pstream::waitRequests();
        }

        // Re-pack points and faces

        label nTris = 0;
        label nPoints = 0;

        forAll (allTriFaces, procI)
        {
            nTris += allTriFaces[procI].size();
            nPoints += allTriPoints[procI].size();
        }

        // Pack points
        triPoints.setSize(nPoints);

        // Prepare point renumbering array
        labelListList renumberPoints(Pstream::nProcs());

        nPoints = 0;

        forAll (allTriPoints, procI)
        {
            const pointField& ptp = allTriPoints[procI];

            renumberPoints[procI].setSize(ptp.size());

            labelList& procRenumberPoints = renumberPoints[procI];

            forAll (ptp, ptpI)
            {
                triPoints[nPoints] = ptp[ptpI];
                procRenumberPoints[ptpI] = nPoints;

                nPoints++;
            }
        }

        // Pack triangles and renumber into complete points on the fly
        triFaces.setSize(nTris);

        nTris = 0;

        forAll (allTriFaces, procI)
        {
            const triFaceList& ptf = allTriFaces[procI];

            const labelList& procRenumberPoints = renumberPoints[procI];

            forAll (ptf, ptfI)
            {
                const triFace& procFace = ptf[ptfI];

                triFace& renumberFace = triFaces[nTris];

                forAll (renumberFace, rfI)
                {
                    renumberFace[rfI] = procRenumberPoints[procFace[rfI]];
                }

                nTris++;
            }
        }
    }

    // Make a complete triSurface from local data
    holeTriMeshPtr_ = new triSurface
    (
        triFaces,
        triPoints
    );

    // Clean up duplicate points and zero sized triangles
    holeTriMeshPtr_->cleanup(false);

    Info<< "Overset region " << name() << ": "
        << holeTriMeshPtr_->size() << " triangles in hole cutting"
        << endl;

    // Debug: write holeTriMesh
    if (Pstream::master())
    {
        if (!holeTriMeshPtr_->empty())
        {
            holeTriMeshPtr_->write(word("holeTriSurface_") + name() + ".vtk");
        }
    }

    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Finished calculating hole tri mesh for region " << name()
            << endl;
    }
}


// ************************************************************************* //
