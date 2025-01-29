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

#include "rotation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotation::calcTransformation
(
    const scalar t
) const
{
    // Rotational motion in rad
    vector eulerAngles = 2*pi*axis_*rpm_/60*t;

    const   quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
    const   septernion TR(septernion(origin_)*R*septernion(-origin_));

    return TR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotation::
rotation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotation::
transformation() const
{
    scalar t = time_.value();

    const septernion TR = calcTransformation(t);

    Info<< "solidBodyMotionFunctions::rotation::"
        << "transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::rotation::velocity() const
{
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

Info<< "solidBodyMotionFunctions::rotation::velocity" << endl;
    const septernion velocity
    (
        (calcTransformation(t).t() - calcTransformation(t - dt).t())/dt,
        (calcTransformation(t).r()/calcTransformation(t - dt).r())/dt
    );

    return velocity;
}


bool Foam::solidBodyMotionFunctions::rotation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("centreOfRotation") >> origin_;
    SBMFCoeffs_.lookup("axisOfRotation") >> axis_;
    SBMFCoeffs_.lookup("rpm") >> rpm_;

    if (mag(axis_) < SMALL)
    {
        FatalErrorInFunction
            << "Axis of rotation has zero length: " << axis_
            << abort(FatalError);
    }

    // Normalise axis
    axis_ /= mag(axis_);

    return true;
}


// ************************************************************************* //
