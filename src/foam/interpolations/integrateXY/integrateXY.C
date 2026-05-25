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

#include "integrateXY.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Type integrateXY
(
    const scalar xLow,
    const scalar xHigh,
    const scalarField& x,
    const Field<Type>& y
)
{
    if (xLow > xHigh)
    {
        FatalErrorInFunction
            << "Bad integration bounds: " << xLow << " - " << xHigh
            << abort(FatalError);
    }

    const label n = x.size();

#   ifdef FULLDEBUG
    // Check ascending order of x
    scalarField diffs(n - 1);

    forAll (diffs, i)
    {
        diffs[i] = x[i + 1] - x[i];
    }

    if (min(x) < VSMALL)
    {
        FatalErrorInFunction
            << "x values not in strict ascending order.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Find the low location

    // Locate intersection be looking for change of sign
    scalarField xDiff = x - xLow;

    label lo = -1;

    for (label i = 0; i < n - 1; i++)
    {
        if (xDiff[i]*xDiff[i + 1] <= 0)
        {
            // Found straddling difference
            lo = i;
            break;
        }
    }

    // If lo is not found, check bounds
    if (lo == -1)
    {
        if (xLow < x[0])
        {
            lo = 0;
        }
        else if (xLow > x[n - 1])
        {
            lo = n - 1;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find lower bound"
                << abort(FatalError);
        }
    }

    // Find the high location
    xDiff = x - xHigh;

    label hi = -1;

    for (label i = 0; i < n - 1; i++)
    {
        if (xDiff[i]*xDiff[i + 1] <= 0)
        {
            // Found straddling difference
            hi = i;
            break;
        }
    }

    // If hi is not found, check bounds
    if (hi == -1)
    {
        if (xHigh < x[0])
        {
            hi = 0;
        }
        else if (xHigh > x[n - 1])
        {
            hi = n - 1;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find upper bound"
                << abort(FatalError);
        }
    }
    Info<< "lo: " << lo << " hi: " << hi << endl;
    Type integral = pTraits<Type>::zero;

    // Calculate integral of first bracket
    {
        Type yLow = y[lo]
            + ((xLow - x[lo])/(x[lo + 1] - x[lo]))*(y[lo + 1] - y[lo]);

        integral += 0.5*(y[lo + 1] + yLow)*(x[lo + 1] - xLow);
    }
    Info<< "First: " << integral << endl;
    // Sum up regular intervals
    for (label i = lo + 1; i < hi; i++)
    {
        integral += 0.5*(y[i + 1] + y[i])*(x[i + 1] - x[i]);
        Info<< i << " Sum: " << integral << endl;
    }

    // Calculate integral of the last bracket
    {
        Type yHigh = y[hi]
            + ((xHigh - x[hi])/(x[hi + 1] - x[hi]))*(y[hi + 1] - y[hi]);
        Info<< "yHigh: " << yHigh << " int: " << 0.5*(yHigh + y[hi])*(xHigh - x[hi]) << endl;
        integral += 0.5*(yHigh + y[hi])*(xHigh - x[hi]);
    }

    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
