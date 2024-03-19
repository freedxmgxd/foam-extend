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

#include "EulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt.ref().internalField() =
            rDeltaT.value()*dt.value()*(1.0 - this->mesh().V0()/
            this->mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*this->mesh().V0()/
                    this->mesh().V()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*this->mesh().V0()/this->mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*this->mesh().V0()/
                    this->mesh().V()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
EulerDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/this->mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*this->mesh().V();

    if (this->mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*this->mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*this->mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/this->mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.value()*this->mesh().V();

    if (this->mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*this->mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*this->mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/this->mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.internalField()*this->mesh().V();

    if (this->mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*this->mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*this->mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phiAbs.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    tmp<fluxFieldType> phiCorr =
        phiAbs.oldTime() - (fvc::interpolate(U.oldTime()) & this->mesh().Sf());

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime(), phiCorr())
           *fvc::interpolate(rDeltaT*rA)*phiCorr
        )
    );
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
) const
{
    dimensionedScalar rDeltaT = 1.0/this->mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phiAbs.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
               *(
                   fvc::interpolate(rA*rho.oldTime())*phiAbs.oldTime()
                 - (fvc::interpolate(rA*rho.oldTime()*U.oldTime())
                  & this->mesh().Sf())
                )
            )
        );
    }
    else if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                )
               *(
                   fvc::interpolate(rA*rho.oldTime())
                  *phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                 - (
                       fvc::interpolate
                       (
                           rA*rho.oldTime()*U.oldTime()
                       ) & this->mesh().Sf()
                   )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    rho.oldTime(),
                    U.oldTime(),
                    phiAbs.oldTime()
                )
               *(
                   fvc::interpolate(rA)*phiAbs.oldTime()
                 - (fvc::interpolate(rA*U.oldTime()) & this->mesh().Sf())
                )
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phiAbs are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
) const
{
    tmp<fluxFieldType> toldTimeFlux =
        (this->mesh().Sf() & faceU.oldTime())*rAUf/this->mesh().time().deltaT();

    if (this->mesh().moving())
    {
        // Mesh is moving, need to take into account the ratio between old and
        // current cell volumes
        volScalarField V0ByV
        (
            IOobject
            (
                "V0ByV",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );
        V0ByV.internalField() = this->mesh().V0()/this->mesh().V();
        V0ByV.correctBoundaryConditions();

        // Correct the flux with interpolated volume ratio
        toldTimeFlux.ref() *= fvc::interpolate(V0ByV);
    }

    return toldTimeFlux;
}


template<class Type>
tmp<surfaceScalarField> EulerDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return this->mesh().phi();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
