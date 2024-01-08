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

#ifndef GeometricFieldReuseFunctions_H
#define GeometricFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
(
    const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initRet = false
)
{
    GeometricField<TypeR, PatchField, GeoMesh>& gf1 =
        const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf1());

    if (reusable(tgf1))
    {
        gf1.rename(name);
        gf1.dimensions().reset(dimensions);
        return tgf1;
    }
    else
    {
        tmp<GeometricField<TypeR, PatchField, GeoMesh>> rtgf
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );

        if (initRet)
        {
            rtgf.ref() == tgf1();
        }

        return rtgf;
    }
}


template
<
    class TypeR,
    class Type1,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
class reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh> >& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, PatchField, GeoMesh>& gf1 =
            const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf1());

        if (tgf1.isTmp())
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            gf1.clearCaches();

            return tgf1;
        }
        else
        {
            return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
            (
                new GeometricField<TypeR, PatchField, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        gf1.instance(),
                        gf1.db()
                    ),
                    gf1.mesh(),
                    dimensions
                )
            );
        }
    }
};


template
<
    class TypeR,
    class Type1,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


template
<
    class TypeR,
    class Type1,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpTmpGeometricField<TypeR, Type1, TypeR, PatchField, GeoMesh>
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh> >& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();
        GeometricField<TypeR, PatchField, GeoMesh>& gf2 =
            const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf2());

        if (tgf2.isTmp())
        {
            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            gf2.clearCaches();

            return tgf2;
        }
        else
        {
            return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
            (
                new GeometricField<TypeR, PatchField, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        gf1.instance(),
                        gf1.db()
                    ),
                    gf1.mesh(),
                    dimensions
                )
            );
        }
    }
};


template
<
    class TypeR,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpTmpGeometricField<TypeR, TypeR, Type2, PatchField, GeoMesh>
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh> >& tgf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, PatchField, GeoMesh>& gf1 =
            const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf1());

        if (tgf1.isTmp())
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            gf1.clearCaches();

            return tgf1;
        }
        else
        {
            return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
            (
                new GeometricField<TypeR, PatchField, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        gf1.instance(),
                        gf1.db()
                    ),
                    gf1.mesh(),
                    dimensions
                )
            );
        }
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
class reuseTmpTmpGeometricField<TypeR, TypeR, TypeR, PatchField, GeoMesh>
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh> > New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh> >& tgf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh> >& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, PatchField, GeoMesh>& gf1 =
            const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf1());
        GeometricField<TypeR, PatchField, GeoMesh>& gf2 =
            const_cast<GeometricField<TypeR, PatchField, GeoMesh>& >(tgf2());

        if (tgf1.isTmp())
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            gf1.clearCaches();

            return tgf1;
        }
        else if (tgf2.isTmp())
        {
            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            gf2.clearCaches();

            return tgf2;
        }
        else
        {
            return tmp<GeometricField<TypeR, PatchField, GeoMesh> >
            (
                new GeometricField<TypeR, PatchField, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        gf1.instance(),
                        gf1.db()
                    ),
                    gf1.mesh(),
                    dimensions
                )
            );
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
