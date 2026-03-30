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

void Foam::oversetRegion::clearOut() const
{
    if (oversetMesh::debug)
    {
        InfoInFunction
            << "Called clearOut for region " << name()
            << endl;
    }

    deleteDemandDrivenData(donorRegionsPtr_);
    deleteDemandDrivenData(acceptorRegionsPtr_);

    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(cutHoleCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);
    deleteDemandDrivenData(eligibleDonorCellsPtr_);

    deleteDemandDrivenData(holeTriMeshPtr_);
    deleteDemandDrivenData(holeSearchPtr_);

    deleteDemandDrivenData(localBoundsPtr_);
    deleteDemandDrivenData(globalBoundsPtr_);
    deleteDemandDrivenData(cellSearchPtr_);

    deleteDemandDrivenData(procBoundBoxesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetRegion::oversetRegion
(
    const word& name,
    const label index,
    const fvMesh& mesh,
    const oversetMesh& oversetMesh,
    const dictionary& dict
)
:
    name_(name),
    index_(index),
    mesh_(mesh),
    oversetMesh_(oversetMesh),
    zoneIndex_(mesh_.cellZones().findZoneID(name_)),
    donorRegionNames_(dict.lookup("donorRegions")),
    fringePtr_(),
    donorRegionsPtr_(nullptr),
    acceptorRegionsPtr_(nullptr),

    acceptorCellsPtr_(nullptr),
    donorCellsPtr_(nullptr),
    cutHoleCellsPtr_(nullptr),
    holeCellsPtr_(nullptr),
    eligibleDonorCellsPtr_(nullptr),

    holeTriMeshPtr_(nullptr),
    holeSearchPtr_(nullptr),

    localBoundsPtr_(nullptr),
    globalBoundsPtr_(nullptr),
    cellSearchPtr_(nullptr),
    procBoundBoxesPtr_(nullptr)
{
    // Check zone index
    if (zoneIndex_ < 0)
    {
        FatalErrorInFunction
            << "Cannot find cell zone for region " << name << nl
            << "Available cell zones: " << mesh_.cellZones().names()
            << abort(FatalError);
    }

    fringePtr_ = oversetFringe::New
    (
        mesh,
        *this,
        dict.subDict("fringe")
    );

    calcBounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetRegion::~oversetRegion()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::oversetFringe& Foam::oversetRegion::fringe() const
{
    if (fringePtr_.empty())
    {
        FatalErrorInFunction
            << "Fringe pointer not allocated. It should have been initialized"
            << " properly at construction. Something went wrong..."
            << abort(FatalError);
    }

    return fringePtr_();
}


const Foam::labelList& Foam::oversetRegion::donorRegions() const
{
    if (!donorRegionsPtr_)
    {
        calcDonorRegions();
    }

    return *donorRegionsPtr_;
}


const Foam::labelList& Foam::oversetRegion::acceptorRegions() const
{
    if (!acceptorRegionsPtr_)
    {
        calcAcceptorRegions();
    }

    return *acceptorRegionsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::acceptors() const
{
    if (!acceptorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *acceptorCellsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::donors() const
{
    if (!donorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *donorCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::cutHoles() const
{
    if (!cutHoleCellsPtr_)
    {
        calcCutHoleCells();
    }

    return *cutHoleCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::holes() const
{
    if (!holeCellsPtr_)
    {
        calcHoleCells();
    }

    return *holeCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::eligibleDonors() const
{
    if (!eligibleDonorCellsPtr_)
    {
        calcEligibleDonorCells();
    }

    return *eligibleDonorCellsPtr_;
}


bool Foam::oversetRegion::holePatchesPresent() const
{
    return !holeTriMesh().empty();
}


const Foam::triSurface& Foam::oversetRegion::holeTriMesh() const
{
    if (!holeTriMeshPtr_)
    {
        calcHoleTriMesh();
    }

    return *holeTriMeshPtr_;
}


const Foam::triSurfaceSearch& Foam::oversetRegion::holeSearch() const
{
    if (!holeSearchPtr_)
    {
        holeSearchPtr_ = new triSurfaceSearch
        (
            holeTriMesh()
        );
    }

    return *holeSearchPtr_;
}


const Foam::boundBox& Foam::oversetRegion::localBounds() const
{
    if (!localBoundsPtr_)
    {
        calcBounds();
    }

    return *localBoundsPtr_;
}


const Foam::boundBox& Foam::oversetRegion::globalBounds() const
{
    if (!globalBoundsPtr_)
    {
        calcBounds();
    }

    return *globalBoundsPtr_;
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::oversetRegion::cellSearch() const
{
    if (!cellSearchPtr_)
    {
        calcCellSearch();
    }

    return *cellSearchPtr_;
}


const Foam::List<Foam::List<Foam::boundBox> >&
Foam::oversetRegion::procBoundBoxes() const
{
    if (!procBoundBoxesPtr_)
    {
        calcProcBoundBoxes();
    }

    return *procBoundBoxesPtr_;
}


void Foam::oversetRegion::update() const
{
    Info<< "oversetRegion " << name() << " update" << endl;

    fringePtr_->update();

    clearOut();
}


// ************************************************************************* //
