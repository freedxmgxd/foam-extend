
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

#include "SLTSDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void SLTSDdtScheme<Type>::relaxedDiag
(
    scalarField& rD,
    const surfaceScalarField& phi
) const
{
    const labelUList& owner = this->mesh().owner();
    const labelUList& neighbour = this->mesh().neighbour();
    scalarField diag(rD.size(), 0.0);

    forAll (owner, faceI)
    {
        if (phi[faceI] > 0.0)
        {
            diag[owner[faceI]] += phi[faceI];
            rD[neighbour[faceI]] += phi[faceI];
        }
        else
        {
            diag[neighbour[faceI]] -= phi[faceI];
            rD[owner[faceI]] -= phi[faceI];
        }
    }

    forAll (phi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
        const labelUList& faceCells = pphi.patch().patch().faceCells();

        forAll (pphi, patchFacei)
        {
            if (pphi[patchFacei] > 0.0)
            {
                diag[faceCells[patchFacei]] += pphi[patchFacei];
            }
            else
            {
                rD[faceCells[patchFacei]] -= pphi[patchFacei];
            }
        }
    }

    rD += (1.0/alpha_ - 2.0)*diag;
}


template<class Type>
tmp<volScalarField> SLTSDdtScheme<Type>::SLrDeltaT() const
{
    const surfaceScalarField& phi =
        this->mesh().objectRegistry::template
        lookupObject<surfaceScalarField>(phiName_);

    const dimensionedScalar& deltaT = this->mesh().time().deltaT();

    tmp<volScalarField> trDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT",
                phi.instance(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("rDeltaT", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& rDeltaT = trDeltaT.ref();

    relaxedDiag(rDeltaT, phi);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        rDeltaT.internalField() = max
        (
            rDeltaT.internalField()/this->mesh().V(),
            scalar(1)/deltaT.value()
        );
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            this->mesh().objectRegistry::template
            lookupObject<volScalarField>(rhoName_)
           .oldTime();

        rDeltaT.internalField() = max
        (
            rDeltaT.internalField()/(rho.internalField()*this->mesh().V()),
            scalar(1)/deltaT.value()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);
    }

    rDeltaT.correctBoundaryConditions();

    return trDeltaT;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
) const
{
    volScalarField rDeltaT = SLrDeltaT();

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
            rDeltaT.internalField()*dt.value()*(1.0 - this->mesh().V0()/
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
SLTSDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    volScalarField rDeltaT = SLrDeltaT();

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
                rDeltaT.internalField()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*this->mesh().V0()/
                    this->mesh().V()
                ),
                rDeltaT.boundaryField()*
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
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    volScalarField rDeltaT = SLrDeltaT();

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
                rDeltaT.internalField()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*this->mesh().V0()/
                    this->mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
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
SLTSDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    volScalarField rDeltaT = SLrDeltaT();

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
                rDeltaT.internalField()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*this->mesh().V0()/
                    this->mesh().V()
                ),
                rDeltaT.boundaryField()*
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
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT = SLrDeltaT()().internalField();

    Info<< "max/min rDeltaT " << gMax(rDeltaT) << " " << gMin(rDeltaT) << endl;

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
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT = SLrDeltaT()().internalField();

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
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT = SLrDeltaT()().internalField();

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
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
) const
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        volScalarField rDeltaT = SLrDeltaT();

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                (
                    fvc::interpolate(rDeltaT*rA)*phi.oldTime()
                  - (
                        fvc::interpolate(rDeltaT*rA*U.oldTime())
                      & this->mesh().Sf()
                    )
                )
            )
        );
    }
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
) const
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    if (this->mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        volScalarField rDeltaT = SLrDeltaT();

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
                   *(
                        fvc::interpolate(rDeltaT*rA*rho.oldTime())*phi.oldTime()
                      - (fvc::interpolate(rDeltaT*rA*rho.oldTime()*U.oldTime())
                      & this->mesh().Sf())
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (
                        U.oldTime(),
                        phi.oldTime()/fvc::interpolate(rho.oldTime())
                    )
                   *(
                        fvc::interpolate(rDeltaT*rA*rho.oldTime())
                       *phi.oldTime()/fvc::interpolate(rho.oldTime())
                      - (
                            fvc::interpolate
                            (
                                rDeltaT*rA*rho.oldTime()*U.oldTime()
                            ) & this->mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimDensity*dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (
                        rho.oldTime(),
                        U.oldTime(),
                        phi.oldTime()
                    )
                   *(
                        fvc::interpolate(rDeltaT*rA)*phi.oldTime()
                      - (
                            fvc::interpolate(rDeltaT*rA*U.oldTime())
                          & this->mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
) const
{
    tmp<fluxFieldType> toldTimeFlux =
        (this->mesh().Sf() & faceU.oldTime())*rAUf*
        fvc::interpolate(SLrDeltaT());

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
tmp<surfaceScalarField> SLTSDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh()
            ,
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
