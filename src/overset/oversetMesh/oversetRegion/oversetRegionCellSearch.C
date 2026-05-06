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

void Foam::oversetRegion::calcBounds() const
{
    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Calculating oversetRegion bounds for region " << name()
            << endl;
    }

    if (localBoundsPtr_ || globalBoundsPtr_)
    {
        FatalErrorInFunction
            << "Bounds already calculated"
            << abort(FatalError);
    }

    // Make a global bounding box for this region
    boolList usedPoints(mesh_.nPoints(), false);

    // Get cells-points from mesh
    const labelListList& pc = mesh_.cellPoints();

    // Get region cell indices
    const labelList& rc = zone();

    forAll (rc, rcI)
    {
        // Get points of region cells
        const labelList& curPc = pc[rc[rcI]];

        forAll (curPc, i)
        {
            usedPoints[curPc[i]] = true;
        }
    }

    // Count used points
    label nUsedPoints = 0;

    forAll (usedPoints, pointI)
    {
        if (usedPoints[pointI])
        {
            nUsedPoints++;
        }
    }

    // Make a list of used points
    const pointField& points = mesh_.points();

    pointField regionPoints(nUsedPoints);

    // Reset point counter to zero
    nUsedPoints = 0;

    forAll (usedPoints, pointI)
    {
        if (usedPoints[pointI])
        {
            regionPoints[nUsedPoints] = points[pointI];
            nUsedPoints++;
        }
    }

    // Local (processor) bounding box is calculated without a reduce
    localBoundsPtr_ = new boundBox(regionPoints, false);

    // Global bounding box is calculated with a reduce
    globalBoundsPtr_ = new boundBox(regionPoints, true);

    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Finished calculating oversetRegion bounds for region " << name()
            << endl;
    }
}


void Foam::oversetRegion::calcCellSearch() const
{
    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Calculating oversetRegion cellSearch for region " << name()
            << endl;
    }

    if (cellSearchPtr_)
    {
        FatalErrorInFunction
            << "Cell tree already calculated"
            << abort(FatalError);
    }

    // Create the octree search for this region.  It will be used by other
    // regions when searching for donor cells

    // Bounding box containing only local region cells
    treeBoundBox overallBb(localBounds());
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1e-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // Search
    cellSearchPtr_ = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,  // Cache bb.  Reconsider for moving mesh cases
            mesh_,
            eligibleDonors()
        ),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );

    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Finished calculating oversetRegion cellSearch for region "
            << name()
            << endl;
    }
}


void Foam::oversetRegion::calcProcBoundBoxes() const
{
    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Calculating oversetRegion procBoundBoxes for region " << name()
            << endl;
    }

    if (procBoundBoxesPtr_)
    {
        FatalErrorInFunction
            << "Processor bounding boxes already calculated"
            << abort(FatalError);
    }

    // Create the list
    procBoundBoxesPtr_ = new List<List<boundBox> >(Pstream::nProcs());
    List<List<boundBox> >& procBoundBoxes = *procBoundBoxesPtr_;

    // Get pointer list of bounding boxes for this processor
    List<boundBox>& localBoundBoxes = procBoundBoxes[Pstream::myProcNo()];

    // Get all regions that are present on this processor
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Set the size for this processor
    localBoundBoxes.setSize(regions.size());

    // Loop through overset regions and populate the list
    forAll (regions, orI)
    {
        // Use local bounds to optimise sending of acceptors
        localBoundBoxes[orI] = regions[orI].localBounds();
    }

    // Now that each processor has filled in its own part, combine the data
    {
        Pstream::gatherList(procBoundBoxes);
        Pstream::scatterList(procBoundBoxes);

        // Block for all requests and remove storage
        Pstream::waitRequests();
    }

    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Finished calculating oversetRegion procBoundBoxes for region "
            << name()
            << endl;
    }
}


// ************************************************************************* //
