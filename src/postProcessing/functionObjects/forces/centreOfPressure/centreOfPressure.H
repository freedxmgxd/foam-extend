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
    Foam::centreOfPressure

Description
    Calculates the centre of pressure over a given list of patches.

    Member function calcCentreOfPressure() calculates and returns the centreOfPressure and
    moments.

    Member function centreOfPressure::write() calls calcCentreOfPressure() and writes the
    centreOfPressure and moments into the file \<timeDir\>/centreOfPressure.dat

SourceFiles
    centreOfPressure.C
    IOcentreOfPressure.H

\*---------------------------------------------------------------------------*/

#ifndef centreOfPressure_H
#define centreOfPressure_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class centreOfPressure Declaration
\*---------------------------------------------------------------------------*/

class centreOfPressure
{
public:

    // Tuple which holds the pressure*x (.first()) and pressure (.second)
    // for centreOfPressure
    typedef Tuple2<vector, scalar> copTuple;

    //- Sum operation class to accumulate the pressure, viscous centreOfPressure
    //  and moments
    class sumOp
    {
    public:

        copTuple operator()
        (
            const copTuple& cop1,
            const copTuple& cop2
        ) const
        {
            return copTuple
            (
                cop1.first() + cop2.first(),
                cop1.second() + cop2.second()
            );
        }
    };


protected:

    // Private data

        //- Name of this set of centreOfPressure,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        // Read from dictionary

            //- Patches to integrate centreOfPressure over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;


        //- CentreOfPressure/moment file ptr
        autoPtr<OFstream> centreOfPressureFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        centreOfPressure(const centreOfPressure&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const centreOfPressure&) = delete;


        //- If the centreOfPressure file has not been created create it
        void makeFile();

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("centreOfPressure");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        centreOfPressure
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~centreOfPressure() = default;


    // Member Functions

        //- Return name of the set of centreOfPressure
        virtual const word& name() const
        {
            return name_;
        }

        //- Calculate and return centreOfPressure and moment
        copTuple calcCentreOfPressure() const;

        //- Read the centreOfPressure data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet()
        {}

        //- Write the centreOfPressure
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
