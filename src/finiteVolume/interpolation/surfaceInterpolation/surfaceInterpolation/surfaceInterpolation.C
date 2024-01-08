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
#include "demandDrivenData.H"
#include "coupledFvPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceInterpolation, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::surfaceInterpolation::clearOut()
{
    deleteDemandDrivenData(weightingFactorsPtr_);
    deleteDemandDrivenData(deltaCoeffsPtr_);
    deleteDemandDrivenData(magLongDeltasPtr_);

    orthogonal_ = false;
    deleteDemandDrivenData(correctionVectorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::surfaceInterpolation(const fvMesh& fvm)
:
    mesh_(fvm),
    schemesDict_(fvm),
    solutionDict_(fvm),
    weightingFactorsPtr_(nullptr),
    deltaCoeffsPtr_(nullptr),
    magLongDeltasPtr_(nullptr),
    orthogonal_(false),
    correctionVectorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::~surfaceInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::surfaceScalarField& Foam::surfaceInterpolation::weights() const
{
    if (!weightingFactorsPtr_)
    {
        makeWeights();
    }

    return *weightingFactorsPtr_;
}


const Foam::surfaceScalarField& Foam::surfaceInterpolation::deltaCoeffs() const
{
    if (!deltaCoeffsPtr_)
    {
        makeDeltaCoeffs();
    }

    return *deltaCoeffsPtr_;
}


const Foam::surfaceScalarField&
Foam::surfaceInterpolation::magLongDeltas() const
{
    if (!magLongDeltasPtr_)
    {
        makeMagLongDeltas();
    }

    return *magLongDeltasPtr_;
}


bool Foam::surfaceInterpolation::orthogonal() const
{
    if (orthogonal_ == false && !correctionVectorsPtr_)
    {
        makeCorrectionVectors();
    }

    return orthogonal_;
}


const Foam::surfaceVectorField&
Foam::surfaceInterpolation::correctionVectors() const
{
    if (orthogonal())
    {
        FatalErrorInFunction
            << "cannot return correctionVectors; mesh is orthogonal"
            << abort(FatalError);
    }

    return *correctionVectorsPtr_;
}


bool Foam::surfaceInterpolation::movePoints()
{
    clearOut();

    return true;
}


void Foam::surfaceInterpolation::makeWeights() const
{
    if (weightingFactorsPtr_)
    {
        FatalErrorInFunction
            << "Weighting factors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoInFunction
            << "Constructing weighting factors for face interpolation"
            << endl;
    }

    weightingFactorsPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "weightingFactors",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless
    );
    surfaceScalarField& weightingFactors = *weightingFactorsPtr_;


    // Set local references to mesh data
    // (note that we should not use fvMesh sliced fields at this point yet
    //  since this causes a loop when generating weighting factors in
    //  coupledFvPatchField evaluation phase)
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // ... and reference to the internal field of the weighting factors
    scalarField& w = weightingFactors.internalField();

    forAll (owner, faceI)
    {
        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.
        scalar SfdOwn = mag(Sf[faceI] & (Cf[faceI] - C[owner[faceI]]));
        scalar SfdNei = mag(Sf[faceI] & (C[neighbour[faceI]] - Cf[faceI]));
        w[faceI] = SfdNei/(SfdOwn + SfdNei);
    }

    forAll (mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].makeWeights
        (
            weightingFactors.boundaryField()[patchi]
        );
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
}


void Foam::surfaceInterpolation::makeDeltaCoeffs() const
{
    if (deltaCoeffsPtr_)
    {
        FatalErrorInFunction
            << "Delta coefficients already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoInFunction
            << "Constructing delta coefficients"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    deltaCoeffsPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "deltaCoeffs",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless/dimLength
    );
    surfaceScalarField& DeltaCoeffs = *deltaCoeffsPtr_;


    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    forAll (owner, faceI)
    {
        vector delta = C[neighbour[faceI]] - C[owner[faceI]];
        vector unitArea = Sf[faceI]/magSf[faceI];

        // Standard cell-centre distance form
        //DeltaCoeffs[faceI] = (unitArea & delta)/magSqr(delta);

        // Slightly under-relaxed form
        //DeltaCoeffs[faceI] = 1.0/mag(delta);

        // More under-relaxed form
        //DeltaCoeffs[faceI] = 1.0/(mag(unitArea & delta) + VSMALL);

        // Stabilised form for bad meshes
        DeltaCoeffs[faceI] = 1.0/max(unitArea & delta, 0.05*mag(delta));
    }

    forAll (DeltaCoeffs.boundaryField(), patchi)
    {
        mesh_.boundary()[patchi].makeDeltaCoeffs
        (
            DeltaCoeffs.boundaryField()[patchi]
        );
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing delta coefficients"
            << endl;
    }
}


void Foam::surfaceInterpolation::makeMagLongDeltas() const
{
    if (magLongDeltasPtr_)
    {
        FatalErrorInFunction
            << "Long deltas already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoInFunction
            << "Constructing long deltas"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    magLongDeltasPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "magLongDeltas",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless/dimLength
    );
    surfaceScalarField& magLongDeltas = *magLongDeltasPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& C = mesh_.C().internalField();
    const vectorField& Cf = mesh_.Cf().internalField();
    const vectorField& Sf = mesh_.Sf().internalField();
    const scalarField& magSf = mesh_.magSf().internalField();

    scalarField& mldIn = magLongDeltas.internalField();

    forAll (owner, faceI)
    {
        // This must be the same as in surfaceInterpolation.C - but it is not!
        // Check.  HJ, 11/May/2020
        scalar SfdOwn = mag(Sf[faceI] & (Cf[faceI] - C[owner[faceI]]));
        scalar SfdNei = mag(Sf[faceI] & (C[neighbour[faceI]] - Cf[faceI]));

        mldIn[faceI] = (SfdOwn + SfdNei)/magSf[faceI];
    }

    forAll (magLongDeltas.boundaryField(), patchi)
    {
        mesh_.boundary()[patchi].makeMagLongDeltas
        (
            magLongDeltas.boundaryField()[patchi]
        );
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing long deltas"
            << endl;
    }
}


void Foam::surfaceInterpolation::makeCorrectionVectors() const
{
    if (correctionVectorsPtr_)
    {
        FatalErrorInFunction
            << "Non-orthogonal correction vectors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoInFunction
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    correctionVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "correctionVectors",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless
    );
    surfaceVectorField& corrVecs = *correctionVectorsPtr_;

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& DeltaCoeffs = deltaCoeffs();

    forAll (owner, faceI)
    {
        vector unitArea = Sf[faceI]/magSf[faceI];
        vector delta = C[neighbour[faceI]] - C[owner[faceI]];

        // If non-orthogonality is over 90 deg, kill correction vector
        // HJ, 27/Feb/2011
        corrVecs[faceI] = pos(unitArea & delta)*
            (unitArea - delta*DeltaCoeffs[faceI]);
    }

    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    forAll (corrVecs.boundaryField(), patchI)
    {
        mesh_.boundary()[patchI].makeCorrVecs
        (
            corrVecs.boundaryField()[patchI]
        );
    }

    // Set max non-orthogonality to zero
    scalar MaxNonOrthog = 0.0;

    // Calculate the non-orthogonality for meshes with 1 face or more
    if (!corrVecs.internalField().empty())
    {
        MaxNonOrthog = max(mag(corrVecs.internalField()));

        forAll (corrVecs.boundaryField(), patchI)
        {
            if (!corrVecs.boundaryField()[patchI].empty())
            {
                MaxNonOrthog =
                    max
                    (
                        MaxNonOrthog,
                        max(mag(corrVecs.boundaryField()[patchI]))
                    );
            }
        }
    }

    reduce(MaxNonOrthog, maxOp<scalar>());

    // Convert to angle
    MaxNonOrthog =
        asin(min(MaxNonOrthog, scalar(1)))*
        180.0/mathematicalConstant::pi;

    if (debug)
    {
        InfoInFunction
            << "maximum non-orthogonality = " << MaxNonOrthog << " deg."
            << endl;
    }

    if (MaxNonOrthog < 5)
    {
        orthogonal_ = true;
        deleteDemandDrivenData(correctionVectorsPtr_);
    }
    else
    {
        orthogonal_ = false;
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
}


// ************************************************************************* //
