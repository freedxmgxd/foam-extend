/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "patchHeatFluxFunctionObject.H"
#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchHeatFluxFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchHeatFluxFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::patchHeatFluxFunctionObject::patchIDList() const
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    labelList patchIDs(mesh.boundary().size(), -1);

    label nPatches = 0;

    forAll (mesh.boundary(), patchI)
    {
        if (!mesh.boundary()[patchI].coupled())
        {
            patchIDs[nPatches] = patchI;
            nPatches++;
        }
    }

    // Resize the list
    patchIDs.setSize(nPatches);

    return patchIDs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchHeatFluxFunctionObject::patchHeatFluxFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    patchName_(dict.lookup("patch")),
    TName_(),
    phiName_(),
    gammaName_(),
    ofPtr_()
{
    // Read parameters
    read(dict);

    Info<< "Creating patchHeatFluxFunctionObject for patch "
        << patchName_ << " on region " << regionName_ << endl;

    // Set stream pointer: only master writes
    if (Pstream::master())
    {
        // Write into case directory instead of processor directory
        if (Pstream::parRun())
        {
            ofPtr_.reset
            (
                new OFstream(time_.path()/".."/word(dict.lookup("file")))
            );
        }
        else
        {
            ofPtr_.reset
            (
                new OFstream(time_.path()/word(dict.lookup("file")))
            );
        }

        // Write header
        ofPtr_()
            << "# Time, patch convective, diffusive, radiative, total flux"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchHeatFluxFunctionObject::start()
{
    return true;
}


bool Foam::patchHeatFluxFunctionObject::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Get patch index
    label patchID = mesh.boundaryMesh().findPatchID(patchName_);

    // Check whether the patch has been found
    if (patchID < 0)
    {
        InfoInFunction
            << "Patch named: " << patchName_ << " not found."
            << " Returning."
            << endl;

        return false;
    }

    if (mesh.foundObject<volScalarField>(TName_))
    {
        const surfaceScalarField& magSf = mesh.magSf();

        const volScalarField& T =
            mesh.lookupObject<volScalarField>(TName_);

        const fvPatchScalarField& patchT = T.boundaryField()[patchID];

        // Calculate the flux through the patch
        scalar convectiveFlux = 0;

        if (mesh.foundObject<surfaceScalarField>(phiName_))
        {
            const surfaceScalarField& phi =
                mesh.lookupObject<surfaceScalarField>(phiName_);

            const scalarField& patchPhi = phi.boundaryField()[patchID];

            // Dimension check and correction
            if (debug)
            {
                if (phi.dimensions() == dimVolume/dimTime)
                {
                    Info<< "Convective flux dimension = m^3/s" << endl;
                }
                else if (phi.dimensions() == dimMass/dimTime)
                {
                    Info<< "Convective flux dimension = kg/s" << endl;
                }
                else if
                (
                    phi.dimensions() == dimSpecificHeatCapacity*dimMass/dimTime
                )
                {
                    Info<< "Convective flux dimension = W/K" << endl;
                }
            }

            if (isA<cyclicPolyPatch>(mesh.boundaryMesh()[patchID]))
            {
                convectiveFlux =
                    gSum
                    (
                        SubField<scalar>(patchPhi, patchPhi.size()/2)*
                        SubField<scalar>(patchT, patchT.size()/2)
                    );
            }
            else
            {
                convectiveFlux = gSum(patchPhi*patchT);
            }
        }

        // Calculate diffusive flux
        scalar diffusiveFlux = 0;

        // Note: diffusive flux can be a scalar or a symmTensor
        if (mesh.foundObject<volScalarField>(gammaName_))
        {
            const volScalarField& gamma =
                mesh.lookupObject<volScalarField>(gammaName_);

            const scalarField& patchGamma = gamma.boundaryField()[patchID];

            // Dimension check and correction
            if (debug)
            {
                if (gamma.dimensions() == sqr(dimLength)/dimTime)
                {
                    Info<< "Flux dimension = m^3/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                {
                    Info<< "Flux dimension = kg/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                    Info<< "Flux dimension = W/K" << endl;
            }

            if (isA<cyclicPolyPatch>(mesh.boundaryMesh()[patchID]))
            {
                diffusiveFlux =
                    gSum
                    (
                        SubField<scalar>(patchGamma, patchGamma.size()/2)*
                        SubField<scalar>(patchT.snGrad(), patchT.size()/2)*
                        SubField<scalar>
                        (
                            magSf.boundaryField()[patchID], patchT.size()/2
                        )
                    );
            }
            else
            {
                diffusiveFlux =
                    gSum
                    (
                        patchGamma*patchT.snGrad()*
                        magSf.boundaryField()[patchID]
                    );
            }
        }
        else if (mesh.foundObject<volSymmTensorField>(gammaName_))
        {
            const volSymmTensorField& gamma =
                mesh.lookupObject<volSymmTensorField>(gammaName_);

            const symmTensorField& patchGamma = gamma.boundaryField()[patchID];

            // Dimension check and correction
            if (debug)
            {
                if (gamma.dimensions() == sqr(dimLength)/dimTime)
                {
                    Info<< "Flux dimension = m^3/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                {
                    Info<< "Flux dimension = kg/s" << endl;
                }
                else if
                (
                    gamma.dimensions() == dimDensity*sqr(dimLength)/dimTime
                )
                    Info<< "Flux dimension = W/K" << endl;
            }

            vectorField Sn =
                mesh.Sf().boundaryField()[patchID]/
                magSf.boundaryField()[patchID];

            if (isA<cyclicPolyPatch>(mesh.boundaryMesh()[patchID]))
            {
                diffusiveFlux =
                    gSum
                    (
                        (
                            (
                                SubField<vector>
                                (
                                    mesh.Sf().boundaryField()[patchID],
                                    patchT.size()/2
                                )
                              & SubField<symmTensor>
                                (
                                    patchGamma,
                                    patchT.size()/2
                                )
                            ) & SubField<vector>(Sn, Sn.size()/2)
                        )*SubField<scalar>(patchT.snGrad(), patchT.size()/2)
                    );
            }
            else
            {
                diffusiveFlux =
                    gSum
                    (
                        (
                            (
                                mesh.Sf().boundaryField()[patchID] & patchGamma
                            ) & Sn
                        )*patchT.snGrad()
                    );
            }
        }

        // Calculate radiative flux
        scalar radiativeFlux = 0;

        // Radiative heat flux on normal patches
        if (mesh.foundObject<volScalarField>("Qr"))
        {
            const volScalarField& Qr =
                mesh.lookupObject<volScalarField>("Qr");

            const scalarField& patchQr = Qr.boundaryField()[patchID];

            const scalarField& patchMagSf =
                mesh.magSf().boundaryField()[patchID];

            if (isA<cyclicPolyPatch>(mesh.boundaryMesh()[patchID]))
            {
                radiativeFlux =
                    gSum
                    (
                        SubField<scalar>(patchQr, patchQr.size()/2)*
                        patchMagSf
                    );
            }
            else
            {
                radiativeFlux =
                    gSum
                    (
                        SubField<scalar>(patchQr, patchQr.size()/2)*
                        SubField<scalar>(patchMagSf, patchMagSf.size()/2)
                    );
            }
        }

        if (ofPtr_.valid())
        {
            ofPtr_()
                << time_.value() << tab
                << convectiveFlux << tab
                << -diffusiveFlux << tab
                << radiativeFlux << tab
                << convectiveFlux - diffusiveFlux + radiativeFlux
                << endl;
        }

        Info<< "Heat flux through " << patchName_
            << " region " << regionName_
            << " (convective, diffusive, radiative, total): "
            << convectiveFlux << tab
            << -diffusiveFlux << tab
            << radiativeFlux << tab
            << convectiveFlux - diffusiveFlux + radiativeFlux
            << endl;

        return true;
    }
    else
    {
        InfoInFunction
            << "Defining field "
            << TName_ << " not found. Returning."
            << endl;

        return false;
    }
}


bool Foam::patchHeatFluxFunctionObject::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    patchName_ = word(dict.lookup("patch"));

    TName_ = word(dict.lookup("T"));
    phiName_ = word(dict.lookupOrDefault<word>("phi", "phi"));
    gammaName_ = word(dict.lookup("gamma"));

    return false;
}

// ************************************************************************* //
