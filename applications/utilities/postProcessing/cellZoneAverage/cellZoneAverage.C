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

Application
    cellZoneAverage

Description
    Calculates the average of the specified field over the specified cellZone

Author
    Marko Horvat, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("cellZoneName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    word fieldName(args.additionalArgs()[0]);
    word cellZoneName(args.additionalArgs()[1]);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            const label cellZoneID = mesh.cellZones().findZoneID(cellZoneName);
            if (cellZoneID < 0)
            {
                FatalError
                    << "Unable to find cellZone " << cellZoneName << nl
                    << exit(FatalError);
            }

            if (fieldHeader.headerClassName() == "volScalarField")
            {
                Info<< "    Reading volScalarField " << fieldName << endl;
                volScalarField field(fieldHeader, mesh);

                scalar volume = 0;
                scalar sumField = 0;
                scalarField volumes = mesh.V();
                const cellZone& zone = mesh.cellZones()[cellZoneID];

                forAll (zone, i)
                {
                    const label cellI = zone[i];
                    sumField += volumes[cellI] * field[cellI];
                    volume += volumes[cellI];
                }

                Info<< "    Average of " << fieldName << " over cellZone "
                    << cellZoneName << " = " << sumField/volume << endl;
            }
            else
            {
                FatalError
                    << "Only possible to average volScalarFields "
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
