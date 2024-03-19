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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "transformField.H"
#include "coeffFields.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    sendBuf_(),
    receiveBuf_(),
    scalarSendBuf_(),
    scalarReceiveBuf_(),
    pnf_()
{
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Clear patch neighbour field
    clearCaches();

    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::processorFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    // Clear patch neighbour field
    clearCaches();

    fvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFvPatchField<Type>::patchNeighbourField() const
{
    if (Pstream::parRun())
    {
        if (pnf_.empty())
        {
            FatalErrorInFunction
                << "Processor patchNeighbourField not ready"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempting to access processor patchNeighbourField "
            << "in a serial run"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Tracing processorFvPatchField: reusing cache : "
            << this->dimensionedInternalField().name()
            << " on patch: " << this->patch().name()
            << endl;
    }

    return pnf_;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        // Set field size for receive
        pnf_.setSize(this->size());

        // Collect data into send buffer
        sendBuf_ = this->patchInternalField();

        if (commsType == Pstream::nonBlocking)
        {
            // Fast path. Receive into pnf_.  HJ, 9/Sep/2021

            outstandingRecvRequest_ = Pstream::nRequests();

            IPstream::read
            (
                Pstream::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<char*>(pnf_.begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );

            outstandingSendRequest_ = Pstream::nRequests();

            OPstream::write
            (
                Pstream::nonBlocking,
                procPatch_.neighbProcNo(),
                reinterpret_cast<const char*>(sendBuf_.begin()),
                this->byteSize(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        else
        {
            procPatch_.send(commsType, sendBuf_);
        }
    }

    // Signal completion not needed
    fvPatchField<Type>::initEvaluate(commsType);
}


template<class Type>
void Foam::processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "Tracing processorFvPatchField: updating cache: "
            << this->dimensionedInternalField().name()
            << " on patch: " << this->patch().name()
            << endl;
    }

    if (Pstream::parRun())
    {
        if (commsType == Pstream::nonBlocking)
        {
            // Fast path. Received into *this

            if
            (
                outstandingRecvRequest_ >= 0
             && outstandingRecvRequest_ < Pstream::nRequests()
            )
            {
                Pstream::waitRequest(outstandingRecvRequest_);
            }
            outstandingSendRequest_ = -1;
            outstandingRecvRequest_ = -1;
        }
        else
        {
            // Receive into pnf_.  HJ, 9/Sep/2021
            procPatch_.receive<Type>(commsType, pnf_);
        }

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), pnf_);
        }
    }

    // Evaluate patch field to carry interpolated values
    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1 - this->patch().weights())*pnf_
    );

    // Signal completion
    fvPatchField<Type>::evaluate(commsType);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFvPatchField<Type>::snGrad() const
{
    // Note: fixed patchNeighbourField storage problem
    return this->patch().deltaCoeffs()*
        (this->patchNeighbourField() - this->patchInternalField());
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    scalarSendBuf_ = this->patch().patchInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        // Check buffer size
        scalarReceiveBuf_.setSize(this->size());

        outstandingRecvRequest_ = Pstream::nRequests();

        IPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();

        OPstream::write
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.send(commsType, scalarSendBuf_);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            Pstream::waitRequest(outstandingRecvRequest_);
        }

        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;
    }
    else
    {
        // Check size
        scalarReceiveBuf_.setSize(this->size());

        procPatch_.receive<scalar>(commsType, scalarReceiveBuf_);
    }

    // The data is now in scalarReceiveBuf_ for both cases

    // Transform according to the transformation tensor
    transformCoupleField(scalarReceiveBuf_, cmpt);

    // Multiply the field by coefficients and add into the result

    const labelUList& faceCells = this->patch().faceCells();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*scalarReceiveBuf_[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*scalarReceiveBuf_[elemI];
        }
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    sendBuf_ = this->patch().patchInternalField(psiInternal);

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        // Check buffer size
        receiveBuf_.setSize(this->size());

        outstandingRecvRequest_ = Pstream::nRequests();

        IPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(receiveBuf_.begin()),
            receiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = Pstream::nRequests();

        OPstream::write
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.begin()),
            sendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.send
        (
            commsType,
            this->patch().patchInternalField(psiInternal)()
        );
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if (commsType == Pstream::nonBlocking)
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            Pstream::waitRequest(outstandingRecvRequest_);
        }

        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;
    }
    else
    {
        // Check size
        receiveBuf_.setSize(this->size());

        procPatch_.receive<Type>(commsType, receiveBuf_);
    }

    // The data is now in receiveBuf_ for both cases

    // Multiply neighbour field with coeffs and re-use buffer for result
    // of multiplication
    multiply(receiveBuf_, coeffs, receiveBuf_);

    const labelUList& faceCells = this->patch().faceCells();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += receiveBuf_[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= receiveBuf_[elemI];
        }
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
bool Foam::processorFvPatchField<Type>::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = Pstream::finishedRequest(outstandingSendRequest_);

        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = Pstream::finishedRequest(outstandingRecvRequest_);

        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}

template<class Type>
void Foam::processorFvPatchField<Type>::clearCaches() const
{
    if (debug && !pnf_.empty())
    {
        Pout<< "Tracing processorFvPatchField: clearing cache: "
            << this->dimensionedInternalField().name()
            << " on patch: " << this->patch().name()
            << endl;
    }

    // Clear patch neighbour field cache
    pnf_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(ul);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator+=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator-=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator*=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator/=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator+=(tf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator-=(tf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator*=(tf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator/=(tf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator=
(
    const Type& t
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(t);
}


// Note: it is possible to manipulate pnf_ in operations with primitive
// types, but I am currently avoiding this because this assumes synchronised
// field functions across multiple processors
// HJ, 8/Sep/2021


template<class Type>
void Foam::processorFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    this->clearCaches();
    fvPatchField<Type>::operator+=(t);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    this->clearCaches();
    fvPatchField<Type>::operator-=(t);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    this->clearCaches();
    fvPatchField<Type>::operator*=(s);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    this->clearCaches();
    fvPatchField<Type>::operator/=(s);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(ptf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(tf);
}


template<class Type>
void Foam::processorFvPatchField<Type>::operator==
(
    const Type& t
)
{
    this->clearCaches();
    fvPatchField<Type>::operator=(t);
}


// ************************************************************************* //
