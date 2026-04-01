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

#include "layeredOverlapFringe.H"
#include "oversetRegion.H"
#include "oversetMesh.H"
#include "polyPatchID.H"
#include "processorFvPatchFields.H"
#include "oversetFvPatchFields.H"
#include "typeInfo.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(layeredOverlapFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, layeredOverlapFringe, dictionary);
}


// * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * //

void Foam::layeredOverlapFringe::evaluateNonOversetBoundaries
(
    volLabelField::GeometricBoundaryField& psib
)
{
    // Code practically copy/pasted from
    // GeometricBoundaryField::updateCoupledPatchFields
    // GeometricBoundaryField should be redesigned to accomodate for such needs
    if
    (
        Pstream::defaultComms() == Pstream::blocking
     || Pstream::defaultComms() == Pstream::nonBlocking
    )
    {
        forAll(psib, patchI)
        {
            // Get fvPatchField
            fvPatchLabelField& psip = psib[patchI];

            if (psip.coupled() && !isA<oversetFvPatchLabelField>(psip))
            {
                psip.initEvaluate(Pstream::defaultComms());
            }
        }

        // Block for any outstanding requests
        if (Pstream::defaultComms() == Pstream::nonBlocking)
        {
            IPstream::waitRequests();
            OPstream::waitRequests();
        }

        forAll(psib, patchI)
        {
            // Get fvPatchField
            fvPatchLabelField& psip = psib[patchI];

            if (psip.coupled() && !isA<oversetFvPatchLabelField>(psip))
            {
                psip.evaluate(Pstream::defaultComms());
            }
        }
    }
    else if (Pstream::defaultComms() == Pstream::scheduled)
    {
        // Get the mesh by looking at first fvPatchField
        const lduSchedule& patchSchedule =
            psib[0].dimensionedInternalField().mesh().globalData().
            patchSchedule();

        forAll(patchSchedule, patchEvalI)
        {
            if (patchSchedule[patchEvalI].init)
            {
                // Get fvPatchField
                fvPatchLabelField psip = psib[patchSchedule[patchEvalI].patch];

                if (psip.coupled() && !isA<oversetFvPatchLabelField>(psip))
                {
                    psip.initEvaluate(Pstream::scheduled);
                }
            }
            else
            {
                // Get fvPatchField
                fvPatchLabelField psip = psib[patchSchedule[patchEvalI].patch];

                if (psip.coupled() && !isA<oversetFvPatchLabelField>(psip))
                {
                    psip.evaluate(Pstream::scheduled);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultComms()]
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layeredOverlapFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorInFunction
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get initial holes and surrounding acceptors
    // Algorithm:
    //    - Create indicator field for correct data exchange accross processor
    //      boundaries
    //    - Get holes from overset region (and optionally from specified set)
    //      and mark immediate neighbours of holes as acceptors

    // Get necessary mesh data
    const fvMesh& mesh = region().mesh();
    const labelListList& cc = mesh.cellCells();

    // Get cut holes from overset region
    const labelList& cutHoles = region().cutHoles();

    // Debug
    if (oversetMesh::debug && cutHoles.empty())
    {
        Pout<< "Did not find any holes to initialise the overlap fringe "
            << "assembly. Proceeding to specified holes..."
            << endl;
    }

    // Initialise mask field for eligible acceptors (cells that are not
    // holes)
    boolList eligibleAcceptors(mesh.nCells(), true);

    // Read user specified holes into allHoles list. Note: if the cell set
    // is not found, the list will be empty
    cellSet allHoles
    (
        mesh,
        holesSetName_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    // Extend allHoles with cutHoles
    allHoles.insert(cutHoles);

    // Note: Because we cannot assume anything about parallel decomposition and
    // we use neighbourhood walk algorithm, there is no easy way to go across
    // processor boundaries. We will create an indicator field marking all
    // possible cut cells and all possible face cells of given patches. Then, we
    // will use this indicator field to transfer the search to the other side.

    // Create the indicator field
    volLabelField processorIndicator
    (
        IOobject
        (
            "processorIndicator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedLabel("minusOne", dimless, -1)
    );
    labelField& processorIndicatorIn = processorIndicator.internalField();

    // Mark all holes
    forAllConstIter (cellSet, allHoles, iter)
    {
        const label& holeCellI = iter.key();

        // Mask eligible acceptors
        eligibleAcceptors[holeCellI] = false;

        // Mark cut hole cell in processor indicator field
        processorIndicatorIn[holeCellI] = 1;
    }

    // Hash set for storing acceptors
    labelHashSet candidateAcceptors(mesh.nCells());

    // Extend nLayers away from initial set. Note postfix increment
    while (layerIter_++ != nLayers_)
    {
        // Append candidate acceptors to holes for this iteration and reset
        // candidate acceptors
        allHoles += candidateAcceptors;
        candidateAcceptors.clear();

        // Loop through all holes and find acceptor candidates
        forAllConstIter (cellSet, allHoles, iter)
        {
            const label& holeCellI = iter.key();

            // Get neighbours of this hole cell
            const labelList& hNbrs = cc[holeCellI];

            // Loop through neighbours of this hole cell
            forAll (hNbrs, nbrI)
            {
                // Check whether the neighbouring cell is eligible
                const label& nbrCellI = hNbrs[nbrI];

                if (eligibleAcceptors[nbrCellI])
                {
                    // Append the cell and mask it to avoid duplicate entries
                    candidateAcceptors.insert(nbrCellI);
                    eligibleAcceptors[nbrCellI] = false;
                    processorIndicatorIn[nbrCellI] = 1;
                }
            }
        }

        // Get boundary field
        volLabelField::GeometricBoundaryField& processorIndicatorBf =
            processorIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        evaluateNonOversetBoundaries(processorIndicatorBf);

        // Loop through boundary field
        forAll (processorIndicatorBf, patchI)
        {
            // Get patch field
            const fvPatchLabelField& chipf = processorIndicatorBf[patchI];

            // Only perform acceptor search if this is a processor boundary
            if (isA<processorFvPatchLabelField>(chipf))
            {
                // Get neighbour field
                const labelField nbrProcIndicator =
                    chipf.patchNeighbourField();

                // Get face cells
                const labelUList& fc = chipf.patch().faceCells();

                // Loop through neighbouring processor field
                forAll (nbrProcIndicator, pfaceI)
                {
                    if
                    (
                        nbrProcIndicator[pfaceI] > 0
                     && eligibleAcceptors[fc[pfaceI]]
                    )
                    {
                        // The cell on the other side is a hole or acceptor from
                        // face cells of a given patch, while the cell on this
                        // side has not been marked yet neither as an acceptor
                        // or as a hole. Append the cell to candidate acceptors
                        // and mark it as ineligible in order to propage the
                        // fringe on this side
                        candidateAcceptors.insert(fc[pfaceI]);
                        eligibleAcceptors[fc[pfaceI]] = false;
                        processorIndicatorIn[fc[pfaceI]] = 1;
                    }
                }
            }
        }
    }

    // Issue an error if no acceptors have been found for initial guess
    if (returnReduce(candidateAcceptors.size(), sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "Did not find any acceptors to begin with." << nl
            << "Check definition of holes for layeredOverlapFringe "
            << "in oversetMeshDict for region: " << this->region().name() << nl
            << "More specifically, check definition of:" << nl
            << "1. holePatches (mandatory entry)" << nl
            << "2. holes (optional entry)" << nl
            << abort(FatalError);
    }

    // Now we have a nice layout of acceptors that are nLayers away from
    // specified holes. Copy into members
    acceptorsPtr_ = new labelList(candidateAcceptors.toc());
    fringeHolesPtr_ = new labelList(allHoles.toc());
}


void Foam::layeredOverlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::layeredOverlapFringe::layeredOverlapFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    fringeHolesPtr_(nullptr),
    acceptorsPtr_(nullptr),
    finalDonorAcceptorsPtr_(nullptr),

    nLayers_(readLabel(dict.lookup("nLayers"))),
    layerIter_(0),

    holesSetName_(dict.lookupOrDefault<word>("holes", word()))
{
    // Sanity check
    if (nLayers_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid number of layers specified. "
            << nl
            << "Please specify value above 1."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::layeredOverlapFringe::~layeredOverlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::layeredOverlapFringe::initSearch
(
    const labelList& candidateAcceptors,
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // Check only
    forAll (donorAcceptorRegionData, aI)
    {
        donorAcceptor& daPair = donorAcceptorRegionData[aI];

        // Check processor ID
        if (daPair.acceptorProcNo() != Pstream::myProcNo())
        {
            FatalErrorInFunction
                << "Donor on different processor: this cannot happen: "
                << "myProc = " << Pstream::myProcNo()
                << " acceptorProcNo = " << daPair.acceptorProcNo()
                << abort(FatalError);
        }
    }
}


void Foam::layeredOverlapFringe::setDonorSuitability
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // Check only
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
        }
    }
}


bool Foam::layeredOverlapFringe::updateIteration
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
            << "allocated. This means that the final overlap has been "
            << "achieved, prohibiting calls to updateIteration."
            << abort(FatalError);
    }

    // Allocate the list by reusing the argument list. Note: all the work has
    // been done in calcAddressing actually since this is not an iterative
    // process by definition
    finalDonorAcceptorsPtr_ = new donorAcceptorList
    (
        donorAcceptorRegionData,
        true
    );

    // Set the flag to true and return
    updateSuitableOverlapFlag(true);

    return foundSuitableOverlap();
}


const Foam::labelList& Foam::layeredOverlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::layeredOverlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


const Foam::donorAcceptorList&
Foam::layeredOverlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorInFunction
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called layeredOverlapFringe::updateIteration() before "
            << "asking for final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorInFunction
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::layeredOverlapFringe::update() const
{
    Info<< "layeredOverlapFringe::update() const" << endl;

    // Clear everything, including acceptors and fringe holes
    clearAddressing();

    // Reset counter to zero
    layerIter_ = 0;

    // Set flag to false and clear final donor/acceptors only
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
