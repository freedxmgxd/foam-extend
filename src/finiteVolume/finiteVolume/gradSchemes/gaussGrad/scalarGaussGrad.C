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

Description
    Specialisation of gaussGrad for scalars. Needed for implicit fvmGrad
    operator for block coupled systems.

\*---------------------------------------------------------------------------*/

#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, vector> > gaussGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    surfaceScalarField weights = this->tinterpScheme_().weights(vf);
    const scalarField& wIn = weights.internalField();

    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs.ref();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    const vectorField& SfIn = mesh.Sf().internalField();

    l = -wIn*SfIn;
    u = (1 - wIn)*SfIn;
    bs.negSumDiag();

    // Boundary contributions
    forAll (vf.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = patch.Sf();
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));

        // Diag contribution
        forAll (pf, faceI)
        {
            d[fc[faceI]] += internalCoeffs[faceI]*pSf[faceI];
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const vectorField pcl = -pw*pSf;
            const vectorField pcu = (1 - pw)*pSf;

            // Coupling  contributions.  Note change of sign because of
            // boundary coeffients are on the other side.  HJ, 3/Dec/2024
            pcoupleLower = -pcl;
            pcoupleUpper = -pcu;
        }
        else
        {
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Boundary contribution
            forAll (pf, faceI)
            {
                source[fc[faceI]] -= boundaryCoeffs[faceI]*pSf[faceI];
            }
        }
    }

    // Interpolation schemes with corrections not accounted for

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
