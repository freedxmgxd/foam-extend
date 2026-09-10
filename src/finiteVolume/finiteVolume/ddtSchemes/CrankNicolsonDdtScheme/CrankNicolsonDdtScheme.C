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

#include "CrankNicolsonDdtScheme.H"
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
template<class GeoField>
CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh
)
:
    GeoField(io, mesh),
    startTimeIndex_(-2) // This field is for a restart and thus correct so set
                        // the start time-index to corespond to a previous run
{
    // Set the time-index to the beginning of the run to ensure the field
    // is updated during the first time-step
    this->timeIndex() = mesh.time().startTimeIndex();
}


template<class Type>
template<class GeoField>
CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh,
    const dimensioned<typename GeoField::value_type>& dimType
)
:
    GeoField(io, mesh, dimType),
    startTimeIndex_(mesh.time().timeIndex())
{}


template<class Type>
template<class GeoField>
label CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
startTimeIndex() const
{
    return startTimeIndex_;
}


template<class Type>
template<class GeoField>
GeoField& CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
operator()()
{
    return *this;
}


template<class Type>
template<class GeoField>
void CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
operator=(const GeoField& gf)
{
    GeoField::operator=(gf);
}


template<class Type>
template<class GeoField>
typename CrankNicolsonDdtScheme<Type>::template DDt0Field<GeoField>&
CrankNicolsonDdtScheme<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
) const
{
    if (!this->mesh().objectRegistry::template foundObject<GeoField>(name))
    {
        const Time& runTime = this->mesh().time();
        word startTimeName = runTime.timeName(runTime.startTime().value());

        if
        (
            (
                runTime.timeIndex() == runTime.startTimeIndex()
             || runTime.timeIndex() == runTime.startTimeIndex() + 1
            )
         && IOobject
            (
                name,
                startTimeName,
                this->mesh()
            ).headerOk()
        )
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        startTimeName,
                        this->mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh()
                )
            );
        }
        else
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensioned<typename GeoField::value_type>
                    (
                        "0",
                        dims/dimTime,
                        pTraits<typename GeoField::value_type>::zero
                    )
                )
            );
        }
    }

    DDt0Field<GeoField>& ddt0 = static_cast<DDt0Field<GeoField>&>
    (
        const_cast<GeoField&>
        (
            this->mesh().objectRegistry::template lookupObject<GeoField>(name)
        )
    );

    return ddt0;
}


template<class Type>
template<class GeoField>
bool CrankNicolsonDdtScheme<Type>::evaluate
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return ddt0.timeIndex() != this->mesh().time().timeIndex();
}


template<class Type>
template<class GeoField>
scalar CrankNicolsonDdtScheme<Type>::coef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (this->mesh().time().timeIndex() - ddt0.startTimeIndex() > 0)
    {
        return 1.0 + ocCoeff_;
    }
    else
    {
        return 1.0;
    }
}


template<class Type>
template<class GeoField>
scalar CrankNicolsonDdtScheme<Type>::coef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (this->mesh().time().timeIndex() - ddt0.startTimeIndex() > 1)
    {
        return 1.0 + ocCoeff_;
    }
    else
    {
        return 1.0;
    }
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicolsonDdtScheme<Type>::rDtCoef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef_(ddt0)/this->mesh().time().deltaT();
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicolsonDdtScheme<Type>::rDtCoef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef0_(ddt0)/this->mesh().time().deltaT0();
}


template<class Type>
template<class GeoField>
tmp<GeoField> CrankNicolsonDdtScheme<Type>::offCentre_
(
    const GeoField& ddt0
) const
{
    if (ocCoeff_ < 1.0)
    {
        return ocCoeff_*ddt0;
    }
    else
    {
        return ddt0;
    }
}


template<class Type>
const FieldField<fvPatchField, Type>& ff
(
    const FieldField<fvPatchField, Type>& bf
)
{
    return bf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + dt.name() + ')',
            dt.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

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

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            dimensionedScalar rDtCoef0 = rDtCoef0_(ddt0);

            ddt0.dimensionedInternalField() =
            (
                (rDtCoef0*dt)*(this->mesh().V0() - this->mesh().V00())
              - this->mesh().V00()*offCentre_(ddt0.dimensionedInternalField())
            )/this->mesh().V0();
        }

        tdtdt.ref().dimensionedInternalField() =
        (
            (rDtCoef*dt)*(this->mesh().V() - this->mesh().V0())
          - this->mesh().V0()*offCentre_(ddt0.dimensionedInternalField())
        )/this->mesh().V();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    this->mesh().V0()*vf.oldTime().internalField()
                  - this->mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDtCoef.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        this->mesh().V()*vf.internalField()
                      - this->mesh().V0()*vf.oldTime().internalField()
                    ) - this->mesh().V0()*offCentre_(ddt0.internalField())
                )/this->mesh().V(),
                rDtCoef.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(vf - vf.oldTime()) - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*rho.value()*
                (
                    this->mesh().V0()*vf.oldTime().internalField()
                  - this->mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*rho.value()*
                    (
                        this->mesh().V()*vf.internalField()
                      - this->mesh().V0()*vf.oldTime().internalField()
                    ) - this->mesh().V0()*offCentre_(ddt0.internalField())
                )/this->mesh().V(),
                rDtCoef.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*rho*(vf - vf.oldTime()) - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    this->mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - this->mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                this->mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        this->mesh().V()*rho.internalField()*vf.internalField()
                      - this->mesh().V0()*rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - this->mesh().V00()*offCentre_(ddt0.internalField())
                )/this->mesh().V(),
                rDtCoef.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()*vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(rho*vf - rho.oldTime()*vf.oldTime())
              - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDtCoef = rDtCoef_(ddt0).value();

    fvm.diag() = rDtCoef*this->mesh().V();

    vf.oldTime().oldTime();

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    this->mesh().V0()*vf.oldTime().internalField()
                  - this->mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.value()*this->mesh().V();

    vf.oldTime().oldTime();

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*rho.value()*
                (
                    this->mesh().V0()*vf.oldTime().internalField()
                  - this->mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.internalField()*this->mesh().V();

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (this->mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    this->mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - this->mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - this->mesh().V00()*offCentre_(ddt0.internalField())
            )/this->mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.oldTime().internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*rho.oldTime().internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*this->mesh().V();
    }

    return tfvm;
}



template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& dUdt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            U.dimensions()
        );

    DDt0Field<fluxFieldType>& dphidt0 =
        ddt0_<fluxFieldType>
        (
            "ddt0(" + phi.name() + ')',
            phi.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(dUdt0);

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
                    rDtCoef.dimensions()*rA.dimensions()*phi.dimensions(),
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        if (evaluate(dUdt0))
        {
            dUdt0 =
                rDtCoef0_(dUdt0)*(U.oldTime() - U.oldTime().oldTime())
              - offCentre_(dUdt0());
        }

        if (evaluate(dphidt0))
        {
            dphidt0 =
                rDtCoef0_(dphidt0)*(phi.oldTime() - phi.oldTime().oldTime())
              - offCentre_(dphidt0());
        }

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
               *fvc::interpolate(rA)
               *(
                    (rDtCoef*phi.oldTime() + offCentre_(dphidt0()))
                  - (
                        fvc::interpolate
                        (
                            (rDtCoef*U.oldTime() + offCentre_(dUdt0()))
                        ) & this->mesh().Sf()
                    )
                )
            )
        );
    }
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
) const
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& dUdt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            U.dimensions()
        );

    DDt0Field<fluxFieldType>& dphidt0 =
        ddt0_<fluxFieldType>
        (
            "ddt0(" + phi.name() + ')',
            U.dimensions()*dimArea
        );

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        this->mesh().time().timeName(),
        this->mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(dUdt0);

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
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            if (evaluate(dUdt0))
            {
                dUdt0 = rDtCoef0_(dUdt0)*
                (
                    U.oldTime() - U.oldTime().oldTime()
                ) - offCentre_(dUdt0());
            }

            if (evaluate(dphidt0))
            {
                dphidt0 = rDtCoef0_(dphidt0)*
                (
                    phi.oldTime()
                  - fvc::interpolate(rho.oldTime().oldTime()/rho.oldTime())
                   *phi.oldTime().oldTime()
                ) - offCentre_(dphidt0());
            }

            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                    (
                        fvc::interpolate(rA*rho.oldTime())
                       *(rDtCoef*phi.oldTime() + offCentre_(dphidt0()))
                      - (
                            fvc::interpolate
                            (
                                rA*rho.oldTime()
                               *(rDtCoef*U.oldTime() + offCentre_(dUdt0()))
                            ) & this->mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            if (evaluate(dUdt0))
            {
                dUdt0 = rDtCoef0_(dUdt0)*
                (
                    U.oldTime() - U.oldTime().oldTime()
                ) - offCentre_(dUdt0());
            }

            if (evaluate(dphidt0))
            {
                dphidt0 = rDtCoef0_(dphidt0)*
                (
                    phi.oldTime()
                   /fvc::interpolate(rho.oldTime())
                  - phi.oldTime().oldTime()
                   /fvc::interpolate(rho.oldTime().oldTime())
                ) - offCentre_(dphidt0());
            }

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
                        fvc::interpolate(rA*rho.oldTime())
                       *(
                           rDtCoef*phi.oldTime()/fvc::interpolate(rho.oldTime())
                         + offCentre_(dphidt0())
                        )
                      - (
                            fvc::interpolate
                            (
                                rA*rho.oldTime()
                               *(rDtCoef*U.oldTime() + offCentre_(dUdt0()))
                            ) & this->mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == rho.dimensions()*dimVelocity
         && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
        )
        {
            if (evaluate(dUdt0))
            {
                dUdt0 = rDtCoef0_(dUdt0)*
                (
                    U.oldTime() - U.oldTime().oldTime()
                ) - offCentre_(dUdt0());
            }

            if (evaluate(dphidt0))
            {
                dphidt0 = rDtCoef0_(dphidt0)*
                (
                    phi.oldTime() - phi.oldTime().oldTime()
                ) - offCentre_(dphidt0());
            }

            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime())*
                    (
                        fvc::interpolate(rA)
                       *(rDtCoef*phi.oldTime() + offCentre_(dphidt0()))
                      - (
                            fvc::interpolate
                            (
                                rA*(rDtCoef*U.oldTime() + offCentre_(dUdt0()))
                            ) & this->mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
) const
{
    // Store old ddt field for faceU, necessary for consistent flux evaluation
    DDt0Field<GeometricField<Type, fvsPatchField, surfaceMesh> >& faceUDdt0 =
        ddt0_<GeometricField<Type, fvsPatchField, surfaceMesh> >
        (
            "ddt0(" + faceU.name() + ')',
            faceU.dimensions()
        );

    // Trigger storage of faceU old time values
    faceU.oldTime().oldTime();

    // Evaluate faceU ddt if necessary
    if (evaluate(faceUDdt0))
    {
        // Update ddt(faceU)
        faceUDdt0 =
        (
            rDtCoef0_(faceUDdt0)*
            (
                faceU.oldTime()
              - faceU.oldTime().oldTime()
            )
          - offCentre_(faceUDdt0())
        );
    }

    // Calculate old time flux
    fluxFieldType oldTimeFlux =
        rAUf*rDtCoef_(faceUDdt0)*(this->mesh().Sf() & faceU.oldTime());

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
        oldTimeFlux *= fvc::interpolate(V0ByV);
    }

    return
        oldTimeFlux
      + rAUf*rDtCoef_(faceUDdt0)*(this->mesh().Sf() & offCentre_(faceUDdt0()));
}


template<class Type>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DDt0Field<surfaceScalarField>& meshPhi0 = ddt0_<surfaceScalarField>
    (
        "meshPhiCN_0",
        dimVolume
    );

    if (evaluate(meshPhi0))
    {
        meshPhi0 =
            coef0_(meshPhi0)*this->mesh().phi().oldTime() - offCentre_(meshPhi0());
    }

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                this->mesh().phi().name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            coef_(meshPhi0)*this->mesh().phi() - offCentre_(meshPhi0())
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
