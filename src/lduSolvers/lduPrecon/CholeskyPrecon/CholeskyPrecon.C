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
    CholeskyPrecon

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "CholeskyPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CholeskyPrecon, 0);

    lduPreconditioner::
        addsymMatrixConstructorToTable<CholeskyPrecon>
        addCholeskyPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CholeskyPrecon::calcPreconDiag()
{
    // Precondition the diagonal
    if (matrix_.symmetric())
    {
        const labelUList& upperAddr = matrix_.lduAddr().upperAddr();
        const labelUList& lowerAddr = matrix_.lduAddr().lowerAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();

        forAll (upper, coeffI)
        {
            preconDiag_[upperAddr[coeffI]] -=
                sqr(upper[coeffI])/preconDiag_[lowerAddr[coeffI]];
        }
    }

    // Invert the diagonal for future use
    forAll (preconDiag_, i)
    {
        preconDiag_[i] = 1.0/preconDiag_[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CholeskyPrecon::CholeskyPrecon
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    preconDiag_(matrix_.diag())
{
    calcPreconDiag();
}


Foam::CholeskyPrecon::CholeskyPrecon
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    preconDiag_(matrix_.diag())
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CholeskyPrecon::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    if (matrix_.asymmetric())
    {
        FatalErrorInFunction
            << "Calling CholeskyPrecon on an assymetric matrix.  "
            << "Please use ILU0 instead"
            << abort(FatalError);
    }

    // Note: coupled boundary updated is not needed because x is zero
    // HJ and VV, 19/Jun/2017

    // Diagonal block
    // Note: multiplication over-write x: no need to initialise
    // HJ, and VV, 19/Jun/2017
    x = b*preconDiag_;

    if (matrix_.symmetric())
    {
        // Addressing
        const labelUList& u = matrix_.lduAddr().upperAddr();

        const labelUList& l = matrix_.lduAddr().lowerAddr();

        // Coeffs
        const scalarField& upper = matrix_.upper();

        const label nCoeffs = upper.size();

        // Forward sweep
        for (label coeffI = 0; coeffI < nCoeffs; coeffI++)
        {
            x[u[coeffI]] -= preconDiag_[u[coeffI]]*upper[coeffI]*x[l[coeffI]];
        }

        // Reverse sweep
        for (label coeffI = nCoeffs - 1; coeffI >= 0; coeffI--)
        {
            x[l[coeffI]] -= preconDiag_[l[coeffI]]*upper[coeffI]*x[u[coeffI]];
        }
    }
}


// ************************************************************************* //
