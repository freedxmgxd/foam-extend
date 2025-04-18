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

#include "centreOfPressure.H"
#include "volFields.H"
#include "dictionary.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(centreOfPressure, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::centreOfPressure::centreOfPressure
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    patchSet_(),
    pName_(word::null),
    centreOfPressureFilePtr_(nullptr)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        WarningInFunction
            << "No fvMesh available, deactivating."
            << endl;

        active_ = false;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::centreOfPressure::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        patchSet_ =
            mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));

        // Optional entry for p
        pName_ = dict.lookupOrDefault<word>("pName", "p");

        // Check whether pName exists,
        // if not deactivate centreOfPressure
        if (!obr_.foundObject<volScalarField>(pName_))
        {
            active_ = false;
            
            WarningInFunction
                << "Could not find " << pName_
                << " in database." << nl << "   De-activating centreOfPressure."
                << endl;
        }
    }
}


void Foam::centreOfPressure::makeFile()
{
    // Create the centreOfPressure file if not already created
    if (centreOfPressureFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating centreOfPressure file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName centreOfPressureDir;

            word startTimeName =
                obr_.time().timeName(obr_.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                centreOfPressureDir =
                    obr_.time().path()/".."/name_/startTimeName;
            }
            else
            {
                centreOfPressureDir = obr_.time().path()/name_/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(centreOfPressureDir);

            // Open new file at start up
            centreOfPressureFilePtr_.reset
            (
                new OFstream(centreOfPressureDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::centreOfPressure::writeFileHeader()
{
    if (centreOfPressureFilePtr_.valid())
    {
        centreOfPressureFilePtr_()
            << "# Time" << tab << "centreOfPressure" << tab << "pressure"
            << endl;
    }
}


void Foam::centreOfPressure::execute()
{
    // Do nothing - only valid on write
}


void Foam::centreOfPressure::end()
{
    // Do nothing - only valid on write
}


void Foam::centreOfPressure::write()
{
    if (active_)
    {
        // Create the centreOfPressure file if not already created
        makeFile();

        copTuple cop = calcCentreOfPressure();

        // Calculate centre of pressure from reduced values
        cop.first() = cop.first()/stabilise(cop.second(), SMALL);

        if (Pstream::master())
        {
            centreOfPressureFilePtr_()
                << obr_.time().value() << tab << cop << endl;

            if (log_)
            {
                Info<< "centreOfPressure = " << cop.first() << tab 
                    << "pressure = " << cop.second()
                    << endl;
            }
        }
    }
}


Foam::centreOfPressure::copTuple
Foam::centreOfPressure::calcCentreOfPressure() const
{
    copTuple cop(vector::zero, scalar(0));

    // Get pressure field
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    // Get mesh
    const fvMesh& mesh = p.mesh();
    
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        scalarField pf =
            magSf.boundaryField()[patchi]*p.boundaryField()[patchi];

        cop.first() += sum(Cf.boundaryField()[patchi]*pf);
        cop.second() += sum(pf);
    }

    reduce(cop, sumOp());

    return cop;
}


// ************************************************************************* //
