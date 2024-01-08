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

Class
    thermalSourceModel

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "thermalSourceModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::thermalSourceModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalSourceModel::thermalSourceModel
(
    const dictionary& dict,
    const volScalarField& T
)
:
    dict_(dict),
    T_(T)
{
    // If thermal dictionary is found, read it
    if (dict_.found("sources"))
    {
        PtrList<entry> entries(dict_.lookup("sources"));
        sources_.setSize(entries.size());

        forAll (sources_, sourceI)
        {
            sources_.set
            (
                sourceI,
                thermalSource::New
                (
                    entries[sourceI].keyword(),
                    T,
                    entries[sourceI].dict()
                )
            );
        }
    }

    Info<< "SOURCE SIZE: " << sources_.size() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermalSourceModel::S() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "heatSource",
                T_.time().timeName(),
                T_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar
            (
                "zero",
                // Watt/m^3
                dimEnergy/dimTime/dimVolume,
                scalar(0)
            )
        )
    );
    volScalarField& source = tsource.ref();

    forAll(sources_, sourceI)
    {
        sources_[sourceI].addSource(source);
    }

    return tsource;
}


// ************************************************************************* //
