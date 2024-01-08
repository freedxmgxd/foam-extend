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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * Private Member Operators  * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isAnyTmp() const
{
    return tmpType_ == REUSABLE_TMP || tmpType_ == NON_REUSABLE_TMP;
}


template<class T>
inline void Foam::tmp<T>::operator++()
{
    ptr_->operator++();

    if (ptr_->count() > 1)
    {
        FatalErrorInFunction
            << "Attempt to create more than 2 tmps referring to same object.  "
            << "This is not allowed: incorrect use of tmp detected for type "
            << this->typeName()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::tmp<T>::tmp(T* tPtr, bool nonReusable)
:
    tmpType_(nonReusable ? NON_REUSABLE_TMP : REUSABLE_TMP),
    ptr_(tPtr)
{
    if (tPtr && !tPtr->unique())
    {
        FatalErrorInFunction
            << "Attempted construction of a " << this->typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const T& tRef)
:
    tmpType_(CONST_REF),
    ptr_(const_cast<T*>(&tRef))
{}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t)
:
    tmpType_(t.tmpType_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        if (ptr_)
        {
            operator++();
        }
        else
        {
            FatalErrorInFunction
                << "attempted copy of a deallocated temporary"
                << " of type " << this->typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>&& t)
:
    tmpType_(t.tmpType_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        t.ptr_ = 0;
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool allowTransfer)
:
    tmpType_(t.tmpType_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        if (ptr_)
        {
            if (allowTransfer && tmpType_ == REUSABLE_TMP)
            {
                t.ptr_ = 0;
            }
            else
            {
                operator++();
            }
        }
        else
        {
            FatalErrorInFunction
                << "attempted copy of a deallocated temporary"
                << " of type " << this->typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::~tmp()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isTmp() const
{
    return tmpType_ == REUSABLE_TMP;
}


template<class T>
inline bool Foam::tmp<T>::empty() const
{
    return (isAnyTmp() && !ptr_);
}


template<class T>
inline bool Foam::tmp<T>::valid() const
{
    return (!isAnyTmp() || (isAnyTmp() && ptr_));
}


template<class T>
inline Foam::word Foam::tmp<T>::typeName() const
{
    return "tmp<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T& Foam::tmp<T>::ref() const
{
    if (isAnyTmp())
    {
         if (!ptr_)
         {
             FatalErrorInFunction
                 << "Temporary of type " << this->typeName() << " deallocated"
                 << abort(FatalError);
         }
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to acquire non-const reference to const object"
            << " from a " << this->typeName()
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "Temporary of type " << this->typeName() << " deallocated"
                << abort(FatalError);
        }

        if (!ptr_->unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type " << this->typeName()
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = 0;

        return ptr;
    }
    else
    {
        return ptr_->clone().ptr();
    }
}


template<class T>
inline void Foam::tmp<T>::clear() const
{
    if (isAnyTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
            ptr_ = 0;
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    if (isAnyTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "temporary of type " << this->typeName() << " deallocated"
                << abort(FatalError);
        }
    }

    return *ptr_;
}


template<class T>
inline Foam::tmp<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (isAnyTmp())
    {
         if (!ptr_)
         {
             FatalErrorInFunction
                 << "Temporary of type " << this->typeName() << " deallocated"
                 << abort(FatalError);
         }
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const for type "
            << this->typeName()
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::tmp<T>::operator->() const
{
    if (isAnyTmp() && !ptr_)
    {
        FatalErrorInFunction
            << "Temporary of type " << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmp<T>::operator=(T* tPtr)
{
    clear();

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated tmp of type "
                << this->typeName()
            << abort(FatalError);
    }

    if (tPtr && !tPtr->unique())
    {
        FatalErrorInFunction
            << "Attempted assignment of a tmp of type" << this->typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    tmpType_ = REUSABLE_TMP;
    ptr_ = tPtr;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    clear();

    if (t.isAnyTmp())
    {
        tmpType_ = t.tmpType_;
        ptr_ = t.ptr_;

        t.ptr_ = 0;
    }
    else
    {
        FatalErrorInFunction
            << "Attempted to assign to a const reference to an object"
            << " of type " << this->typeName()
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>&& t)
{
    clear();

    tmpType_ = t.tmpType_;
    ptr_ = t.ptr_;

    if (isAnyTmp())
    {
        t.ptr_ = 0;
    }
}


// ************************************************************************* //
