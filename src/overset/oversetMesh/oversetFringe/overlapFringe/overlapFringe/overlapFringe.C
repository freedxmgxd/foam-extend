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
#include "polyPatchID.H"
#include "processorFvPatchFields.H"
#include "oversetFvPatchFields.H"
#include "typeInfo.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, overlapFringe, dictionary);
}


// * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * //

void Foam::overlapFringe::evaluateNonOversetBoundaries
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
        forAll (psib, patchI)
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
            Pstream::waitRequests();
        }

        forAll (psib, patchI)
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

        forAll (patchSchedule, patchEvalI)
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

void Foam::overlapFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorInFunction
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get initial guess for holes and acceptors.
    // Algorithm:
    //    - Create indicator field for correct data exchange accross processor
    //      boundaries
    //    - Get holes from overset region (and optionally from specified set)
    //      and mark immediate neighbours of holes as acceptors
    //    - Loop through (optionally) user specified patches for
    //      initialising the overlap fringe assembly, marking face cells

    // Get necessary mesh data
    const fvMesh& mesh = region().mesh();
    const labelListList& cc = mesh.cellCells();

    // Note: Because we cannot assume anything about parallel decomposition and
    // we use neighbourhood walk algorithm, there is no easy way to go across
    // processor boundaries. We will create an indicator field marking all
    // possible cut cells and all possible face cells of given patches. Then, we
    // will use this indicator field to transfer the search to the other side.

    // Meaning of indicator values
    // -1 = unset
    //  1 = ineligible acceptor: hole

    // Create the indicator field
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

    // Get cut holes from overset region
    const labelList& cutHoles = region().cutHoles();

    // Debug
    if (oversetMesh::debug && cutHoles.empty())
    {
        Pout<< "Did not find any holes to initialise the overlap fringe "
            << "assembly. Proceeding to patches..."
            << endl;
    }

    // Initialise mask field for eligible acceptors (cells that are not
    // holes)
    boolList eligibleAcceptors(mesh.nCells(), true);

    // Read user specified holes into allHoles list. Note: use cellZone rather
    // than cellSet to have correct behaviour on dynamic mesh simulations
    // We will silently proceed if the zone is not found since this option is
    // not mandatory but is useful in certain cases

    // Create a hash set for allHoles
    labelHashSet allHoles;

    // Collect holes that are specified as a cell zone
    // Scoping
    {
        // Get zone index
        const label holeZoneID = mesh.cellZones().findZoneID(holesZoneName_);

        if (holeZoneID > -1)
        {
            // Get the zone for holes and append them to set
            const labelList& specifiedHoles = mesh.cellZones()[holeZoneID];

            allHoles.insert(specifiedHoles);
        }
    }
    // else silently proceed without user-specified holes

    // Extend allHoles with cutHoles
    forAll (cutHoles, chI)
    {
        // Note: duplicated are removed by hash set
        allHoles.insert(cutHoles[chI]);
    }

    // Mark all holes as ineligible acceptors
    forAllConstIter (labelHashSet, allHoles, iter)
    {
        const label& holeCellI = iter.key();

        // Mask eligible acceptors
        eligibleAcceptors[holeCellI] = false;

        // Mark cut hole cell in indicator field
        cellTypeIndicatorIn[holeCellI] = 1;
    }


    // Dynamic list for storing acceptors.
    // Note: inserting duplicates is avoided by updating eligibleAcceptors
    // mask
    dynamicLabelList candidateAcceptors(Foam::max(50, mesh.nCells()/5));

    // Loop through all holes and find acceptor candidates
    forAllConstIter (labelHashSet, allHoles, iter)
    {
        // Get neighbours of this hole cell
        const labelList& hNbrs = cc[iter.key()];

        // Loop through neighbours of this hole cell
        forAll (hNbrs, nbrI)
        {
            // Check whether the neighbouring cell is eligible
            const label& nbrCellI = hNbrs[nbrI];

            if (eligibleAcceptors[nbrCellI])
            {
                // Append the cell and mask it to avoid duplicate entries
                candidateAcceptors.append(nbrCellI);
                eligibleAcceptors[nbrCellI] = false;
            }
        }
    }

    // Debug
    if (oversetMesh::debug() && initPatchNames_.empty())
    {
        Pout<< "Did not find any specified patches to initialise the "
            << "overlap fringe assembly."
            << endl;
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Loop through patches and mark face cells as eligible acceptors
    forAll (initPatchNames_, nameI)
    {
        const polyPatchID curPatch
        (
            initPatchNames_[nameI],
            mesh.boundaryMesh()
        );

        if (!curPatch.active())
        {
            FatalErrorInFunction
                << "Patch specified for fringe initialisation "
                << initPatchNames_[nameI] << " cannot be found"
                << abort(FatalError);
        }

        const labelUList& curFaceCells =
            mesh.boundaryMesh()[curPatch.index()].faceCells();

        // Loop through face cells and mark candidate acceptors if
        // eligible
        forAll (curFaceCells, fcI)
        {
            // Get cell index
            const label& cellI = curFaceCells[fcI];

            // Mark acceptor face cell in indicator field
            cellTypeIndicatorIn[cellI] = 1;

            // Check if the cell is eligible and if it is in region zone
            // (Note: the second check is costly)
            if
            (
                eligibleAcceptors[cellI]
             && rcz.whichCell(cellI) > -1
            )
            {
                candidateAcceptors.append(cellI);
                eligibleAcceptors[cellI] = false;
            }
        }
    }

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
                    // face cells of a given patch, while the cell on this side
                    // has not been marked yet neither as an acceptor or as a
                    // hole. Append the cell to candidate acceptors and mark it
                    // as ineligible in order to propage the fringe on this side
                    candidateAcceptors.append(fc[pfaceI]);
                    eligibleAcceptors[fc[pfaceI]] = false;
                }
            }
        }
    }

    // Issue an error if no acceptors have been found for initial guess
    if (returnReduce(candidateAcceptors.size(), sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "Did not find any acceptors to begin with." << nl
            << "Check definition of holes for overlapFringe in oversetMeshDict"
            << " for region: " << this->region().name() << nl
            << "More specifically, check definition of:" << nl
            << "1. holePatches (mandatory entry)" << nl
            << "2. holes (optional entry)" << nl
            << "3. initPatchNames (optional entry)"
            << abort(FatalError);
    }


    // Now we have a decent first guess for acceptors that will be used as
    // an initial condition for the iterative overlap assembly
    // process.
    // Transfer the acceptor list and allocate empty fringeHoles list, which
    // may be populated in updateIteration member function
    acceptorsPtr_ = new labelList(candidateAcceptors.shrink());
    fringeHolesPtr_ = new labelList(allHoles.toc());
}


void Foam::overlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::overlapFringe::overlapFringe
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

    holesZoneName_(dict.lookupOrDefault<word>("holes", word())),
    initPatchNames_
    (
        dict.lookupOrDefault<wordList>("initPatchNames", wordList())
    ),

    donorSuitability_
    (
        donorSuitability::donorSuitability::New(*this, dict)
    ),
    minGlobalFraction_
    (
        readScalar(dict.lookup("suitablePairFraction"))
    ),
    cumulativeDonorAcceptorsPtr_(nullptr),
    cacheFringe_(dict.lookupOrDefault<Switch>("cacheFringe", false)),
    fringeIter_(0)
{
    // Sanity check
    if (minGlobalFraction_ < SMALL || minGlobalFraction_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid suitablePairFraction found while reading the overlap "
            << "fringe dictionary."
            << nl
            << "Please specify value between 0 and 1."
            << exit(FatalIOError);
    }

    Info<< "initPatchNames = " << initPatchNames_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapFringe::~overlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::overlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    // Debug, write fringe holes as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor fringe holes into a cell set." << endl;

        cellSet holesSet
        (
            mesh(),
            "fringeHolesProc" + name(Pstream::myProcNo()) + region().name(),
            labelHashSet(*fringeHolesPtr_)
        );

        holesSet.write();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::overlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    // Debug, write candidate acceptors as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor candidate acceptors into a cell set." << endl;

        cellSet candidateAcceptorsSet
        (
            mesh(),
            "candidateAcceptorsProc" + name(Pstream::myProcNo())
          + region().name(),
            labelHashSet(*acceptorsPtr_)
        );

        candidateAcceptorsSet.write();
    }

    return *acceptorsPtr_;
}


const Foam::donorAcceptorList& Foam::overlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorInFunction
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called overlapFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorInFunction
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    if (cacheFringe_)
    {
        // Get reference to final donor/acceptor pairs
        const donorAcceptorList& finalDAPairs = *finalDonorAcceptorsPtr_;

        // Clear acceptors
        deleteDemandDrivenData(acceptorsPtr_);

        // Now allocate with the correct size. Note: since it is expected that
        // acceptorsPtr_->size() (before destruction) is smaller than
        // finalDonorAcceptorsPtr_.size(), this destruction and initialization
        // should not represent an overhead.
        acceptorsPtr_ = new labelList(finalDAPairs.size());
        labelList& acceptors = *acceptorsPtr_;

        // Set acceptors for the next fringe assembly process
        forAll (finalDAPairs, daPairI)
        {
            acceptors[daPairI] = finalDAPairs[daPairI].acceptorCell();
        }

        // Note: fringe holes actually hold the complete list, there's nothing
        // to do
        Info<< "Cached "
            << returnReduce<label>(acceptorsPtr_->size(), sumOp<label>())
            << " acceptors and "
            << returnReduce<label>(fringeHolesPtr_->size(), sumOp<label>())
            << " fringe holes for region: " << region().name()
            << endl;
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::overlapFringe::update() const
{
    Info<< "overlapFringe::update() const" << endl;

    if (cacheFringe_)
    {
        // If the cache is switched on, simply do not clear acceptorsPtr_ and
        // fringeHolesPtr_ which now hold all final acceptors and holes that
        // will be used to start the next iteration

        // Now clear final and cumulative donor acceptors
        deleteDemandDrivenData(finalDonorAcceptorsPtr_);
        deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
    }
    else
    {
        // Clear everything, including acceptors and fringe holes
        clearAddressing();
    }

    // Reset iteration counter
    fringeIter_ = 0;

    // Set flag to false
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
