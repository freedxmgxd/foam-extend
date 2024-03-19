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

#ifndef FieldReuseFunctions_H
#define FieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR>
tmp<Field<TypeR>> New
(
    const tmp<Field<TypeR>>& tf1,
    const bool initRet = false
)
{
    if (tf1.isTmp())
    {
        return tf1;
    }
    else
    {
        tmp<Field<TypeR>> rtf(new Field<TypeR>(tf1().size()));

        if (initRet)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
}


template<class TypeR, class Type1>
class reuseTmp
{
public:

    static tmp<Field<TypeR> > New(const tmp<Field<Type1> >& tf1)
    {
        return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
    }
};


template<class TypeR>
class reuseTmp<TypeR, TypeR>
{
public:

    static tmp<Field<TypeR> > New(const tmp<Field<TypeR> >& tf1)
    {
        if (tf1.isTmp())
        {
            // Clear caches before reuse
            tf1().clearCaches();

            return tf1;
        }
        else
        {
            return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
        }
    }
};


template<class TypeR, class Type1, class Type12, class Type2>
class reuseTmpTmp
{
public:

    static tmp<Field<TypeR> > New
    (
        const tmp<Field<Type1> >& tf1,
        const tmp<Field<Type2> >& tf2
    )
    {
        return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
    }
};


template<class TypeR, class Type1, class Type12>
class reuseTmpTmp<TypeR, Type1, Type12, TypeR>
{
public:

    static tmp<Field<TypeR> > New
    (
        const tmp<Field<Type1> >& tf1,
        const tmp<Field<TypeR> >& tf2
    )
    {
        if (tf2.isTmp())
        {
            // Clear caches before reuse
            tf2().clearCaches();

            return tf2;
        }
        else
        {
            return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
        }
    }
};


template<class TypeR, class Type2>
class reuseTmpTmp<TypeR, TypeR, TypeR, Type2>
{
public:

    static tmp<Field<TypeR> > New
    (
        const tmp<Field<TypeR> >& tf1,
        const tmp<Field<Type2> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            // Clear caches before reuse
            tf1().clearCaches();

            return tf1;
        }
        else
        {
            return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
        }
    }
};


template<class TypeR>
class reuseTmpTmp<TypeR, TypeR, TypeR, TypeR>
{
public:

    static tmp<Field<TypeR> > New
    (
        const tmp<Field<TypeR> >& tf1,
        const tmp<Field<TypeR> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            // Clear caches before reuse
            tf1().clearCaches();

            return tf1;
        }
        else if (tf2.isTmp())
        {
            // Clear caches before reuse
            tf2().clearCaches();

            return tf2;
        }
        else
        {
            return tmp<Field<TypeR> >(new Field<TypeR>(tf1().size()));
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
