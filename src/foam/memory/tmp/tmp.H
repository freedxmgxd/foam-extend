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
    Foam::tmp

Description
    A class for managing temporary objects

SourceFiles
    tmpI.H

\*---------------------------------------------------------------------------*/

#ifndef tmp_H
#define tmp_H

#include "refCount.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class tmp Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmp
{
    // Private Data

        //- Temporary Object Types
        enum typeOfTmp
        {
            REUSABLE_TMP,
            NON_REUSABLE_TMP,
            CONST_REF
        };

        //- Type of object
        typeOfTmp tmpType_;

        //- Pointer to object
        mutable T* ptr_;


    // Private Member Operators

        //- Return true if the object is a reusable or non-reusable
        inline bool isAnyTmp() const;

        //- Increase count with check: maximum of 2 tmp objects pointing
        //  to same tmp data
        inline void operator++();


public:

    // Constructors

        //- Store object pointer
        inline explicit tmp(T* = 0, bool nonReusable = false);

        //- Store object const reference
        inline tmp(const T&);

        //- Construct copy and increment reference count
        inline tmp(const tmp<T>&);

        // //- Construct copy moving content, does not increment reference count
        // inline tmp(const tmp<T>&&);

        //- Construct copy transferring content of temporary if required
        inline tmp(const tmp<T>&, bool allowTransfer);


    //- Destructor, delete object when reference count is 0
    inline ~tmp();


    // Member Functions

        // Access

            //- Return true if this is really a temporary object
            inline bool isTmp() const;

            //- Return true if this temporary object is empty,
            //  ie, a temporary without allocation
            inline bool empty() const;

            //- Is this temporary object valid,
            //  ie, it is a reference or a temporary that has been allocated
            inline bool valid() const;

            //- Return the type name of the tmp
            //  constructed from the type name of T
            inline word typeName() const;


        // Edit

            //- Return non-const reference or generate a fatal error
            //  if the object is const.
            inline T& ref() const;

            //- Return tmp pointer for reuse
            inline T* ptr() const;

            //- If object pointer points to valid object:
            //  delete object and set pointer to nullptr
            inline void clear() const;


    // Member Operators

        //- Const dereference operator
        inline const T& operator()() const;

        //- Const cast to the underlying type reference
        inline operator const T&() const;

        //- Return object pointer
        inline T* operator->();

        //- Return const object pointer
        inline const T* operator->() const;

        //- Assignment to pointer changing this tmp to a temporary T
        inline void operator=(T*);

        //- Assignment operator
        inline void operator=(const tmp<T>&);

        //- Move assignment transferring the temporary T to this tmp
        inline void operator=(const tmp<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
