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

#include "objectRegistry.H"
#include "finiteRotation.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::vector Foam::finiteRotation::rotVector(const tensor& rotT)
{
    vector ur = - *( inv(I + rotT) & (I - rotT) );

    // Scaling to a unit vector. Problems with round-off. HJ, 4/Aug/2008
    if (mag(ur) > SMALL)
    {
        return ur/(mag(ur) + SMALL);
    }
    else
    {
        // Rotation vector is undertermined at zero rotation
        // Returning arbitrary unit vector. HJ, 4/Mar/2015
        return vector(0, 0, 1);
    }
}


Foam::scalar Foam::finiteRotation::rotAngle(const tensor& rotT)
{
    // Alternative formulation: Daniel Schmode, 15/Feb/2009
    scalar x = rotT.zy() - rotT.yz();
    scalar y = rotT.xz() - rotT.zx();
    scalar z = rotT.yx() - rotT.xy();

    scalar r = hypot(x, hypot(y, z));
    scalar t = tr(rotT);

    return atan2(r, t - 1);
}


Foam::vector Foam::finiteRotation::eulerAngles(const tensor& rotT)
{
    // Create a vector containing euler angles (x = roll, y = pitch, z = yaw)
    vector eulerAngles;

    scalar& rollAngle = eulerAngles.x();
    scalar& pitchAngle = eulerAngles.y();
    scalar& yawAngle = eulerAngles.z();

    // Calculate pitch angle
    scalar xz = rotT.xz();
    if (xz >= 1.0)
    {
        pitchAngle = M_PI / 2.0; // set to the closest value within the domain
    }
    else if (xz <= -1.0)
    {
        pitchAngle = -M_PI / 2.0; // set to the closest value within the domain
    }
    else
    {
        pitchAngle = asin(xz);
    }

    // Calculate roll and yaw angles
    const scalar cosPitch = cos(pitchAngle);
    if (mag(cosPitch) > SMALL) // Check if not in gimbal lock
    {
        rollAngle = atan2(-rotT.yz(), rotT.zz());
        yawAngle = atan2(-rotT.xy(), rotT.xx());
    }
    else
    {
        // Gimbal lock occurred, use an alternative representation
        rollAngle = atan2(rotT.zy(), rotT.yy());
        // Set yaw angle to an arbitrary value (e.g., 0.0) since it is
        // undetermined in gimbal lock
        yawAngle = 0.0;
    }

    return eulerAngles;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteRotation::finiteRotation(const HamiltonRodriguezRot& rot)
:
    eInitial_(rot),
    rotTensor_(rot.R()),
    rotIncrementTensor_(tensor::zero)
{}


Foam::finiteRotation::finiteRotation
(
    const vector& r,
    const scalar& angle
)
:
    eInitial_(HamiltonRodriguezRot(r, angle)),
    rotTensor_(eInitial_.R()),
    rotIncrementTensor_(tensor::zero)
{}


Foam::finiteRotation::finiteRotation(const tensor& R)
:
    eInitial_(R),
    rotTensor_(R),
    rotIncrementTensor_(tensor::zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteRotation::updateRotation(const tensor& R)
{
    rotIncrementTensor_ = (R & rotTensor_.T());
    rotTensor_ = R;
}


void Foam::finiteRotation::updateRotation(const HamiltonRodriguezRot& rot)
{
    updateRotation(rot.R());
}


void Foam::finiteRotation::updateRotationWithIncrement(const tensor& incR)
{
    rotIncrementTensor_ = incR;
    rotTensor_ = (incR & rotTensor_);
}


const Foam::HamiltonRodriguezRot& Foam::finiteRotation::eInitial() const
{
    return eInitial_;
}


Foam::HamiltonRodriguezRot Foam::finiteRotation::eCurrent() const
{
    return HamiltonRodriguezRot(rotVector(), rotAngle());
}


const Foam::tensor& Foam::finiteRotation::rotTensor() const
{
    return rotTensor_;
}


Foam::vector Foam::finiteRotation::rotVector() const
{
    return rotVector(rotTensor_);
}


Foam::scalar Foam::finiteRotation::rotAngle() const
{
    return rotAngle(rotTensor_);
}


Foam::vector Foam::finiteRotation::eulerAngles() const
{
    return eulerAngles(rotTensor_);
}


const Foam::tensor& Foam::finiteRotation::rotIncrementTensor() const
{
    return rotIncrementTensor_;
}


Foam::vector Foam::finiteRotation::rotIncrementVector() const
{
    return rotVector(rotIncrementTensor_);
}


Foam::scalar Foam::finiteRotation::rotIncrementAngle() const
{
    return rotAngle(rotIncrementTensor_);
}


Foam::vector Foam::finiteRotation::omegaAverage(const scalar deltaT) const
{
    return (rotIncrementAngle()/deltaT)*rotIncrementVector();
}


Foam::tensor Foam::finiteRotation::rotTensorAverage() const
{
    return
    (
        HamiltonRodriguezRot
        (
            rotIncrementVector(),
            0.5*rotIncrementAngle()
        ).R().T()
      & rotTensor_
    );
}


Foam::vector Foam::finiteRotation::rotVectorAverage() const
{
    return rotVector(rotTensorAverage());
}


Foam::vector Foam::finiteRotation::omegaAverageAbsolute
(
    const scalar deltaT
) const
{
    return (rotTensorAverage().T() & omegaAverage(deltaT));
}


// ************************************************************************* //
