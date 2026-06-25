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

#include "gaussDivScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template<class Type>
// tmp
// <
//     GeometricField
//     <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
// >
// gaussDivScheme<Type>::fvcDiv
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// ) const
// {
//     tmp
//     <
//         GeometricField
//         <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
//     > tDiv
//     (
//         fvc::surfaceIntegrate
//         (
//             this->mesh_.Sf() & this->tinterpScheme_().interpolate(vf)
//         )
//     );

//     tDiv.ref().rename("div(" + vf.name() + ')');

//     return tDiv;
// }


// NEW FORMULATION: deltas.  See CJ Marooney OFW17.  HJ, 8/Dec/2022
template<class Type>
tmp
<
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
>
gaussDivScheme<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef typename innerProduct<vector, Type>::type DivType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<DivType, fvPatchField, volMesh> > tgDiv
    (
        new GeometricField<DivType, fvPatchField, volMesh>
        (
            IOobject
            (
                word("div(" + vf.name() + ')'),
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<DivType>
            (
                "zero",
                vf.dimensions()/dimLength,
                pTraits<DivType>::zero
            ),
            extrapolatedCalculatedFvPatchField<DivType>::typeName
        )
    );
    GeometricField<DivType, fvPatchField, volMesh>& gDiv = tgDiv.ref();

    // Get weights
    surfaceScalarField w = this->tinterpScheme_().weights(vf);

    const surfaceVectorField& Sf = mesh.Sf();
    
    // updateCoupledPatchFields for patchNeighbourField update
    // HJ, 10/Sep/2021
    vf.boundaryField().updateCoupledPatchFields();

    // Get access to internal fields

    const Field<Type>& vfIn = vf.internalField();

    Field<DivType>& gDivIn = gDiv.internalField();

    const scalarField& wIn = w.internalField();
    const vectorField& SfIn = Sf.internalField();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    label ownFaceI, neiFaceI;

    forAll (own, faceI)
    {
        ownFaceI = own[faceI];
        neiFaceI = nei[faceI];

        Type deltaVf = vfIn[neiFaceI] - vfIn[ownFaceI];

        // Both Sf and own-nei swap values: sign remains the same
        // HJ, 8/Dec/2022
        gDivIn[ownFaceI] += (1 - wIn[faceI])*(SfIn[faceI] & deltaVf);
        gDivIn[neiFaceI] += wIn[faceI]*(SfIn[faceI] & deltaVf);
    }

    // Boundary faces
    forAll (vf.boundaryField(), patchI)
    {
        const scalarField& patchW = w.boundaryField()[patchI];

        const vectorField& patchSf = Sf.boundaryField()[patchI];

        const labelUList& faceCells =
            gDiv.boundaryField()[patchI].patch().faceCells();

        if (vf.boundaryField()[patchI].coupled())
        {
            Field<Type> neiVf =
                vf.boundaryField()[patchI].patchNeighbourField();

            forAll (neiVf, patchFaceI)
            {
                gDiv[faceCells[patchFaceI]] +=
                    patchW[patchFaceI]*
                    (
                        patchSf[patchFaceI]
                      & (neiVf[patchFaceI] - vfIn[faceCells[patchFaceI]])
                    );
            }
        }
        else
        {
            const fvPatchField<Type>& patchVf = vf.boundaryField()[patchI];

            forAll (patchVf, patchFaceI)
            {
                gDiv[faceCells[patchFaceI]] +=
                    patchW[patchFaceI]*
                    (
                        patchSf[patchFaceI]
                      & (patchVf[patchFaceI] - vfIn[faceCells[patchFaceI]])
                    );
            }
        }
    }

    gDivIn /= mesh.V();

    gDiv.correctBoundaryConditions();

    return tgDiv;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> gaussDivScheme<Type>::fvmUDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorInFunction
        << "Implicit div operator defined only for vector."
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> gaussDivScheme<Type>::fvmUDiv
(
    const surfaceScalarField& flux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorInFunction
        << "Implicit div operator defined only for vector."
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
