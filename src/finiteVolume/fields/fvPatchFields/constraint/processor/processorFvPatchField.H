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
    Foam::processorFvPatchField

Description
    This boundary condition enables processor communication across patches.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            processor;
    }
    \endverbatim

SourceFiles
    processorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorFvPatchField_H
#define processorFvPatchField_H

#include "coupledFvPatchField.H"
#include "processorLduInterfaceField.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFvPatchField
:
    public processorLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const processorFvPatch& procPatch_;

        // Sending and receiving

            //- Outstanding request
            mutable label outstandingSendRequest_;

            //- Outstanding request
            mutable label outstandingRecvRequest_;

            //- Send buffer.
            mutable Field<Type> sendBuf_;

            //- Receive buffer.
            mutable Field<Type> receiveBuf_;

            //- Scalar send buffer
            mutable Field<scalar> scalarSendBuf_;

            //- Scalar receive buffer
            mutable Field<scalar> scalarReceiveBuf_;


        //- Cached patch neighbour field
        mutable Field<Type> pnf_;


public:

    //- Runtime type information
    TypeName(processorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorFvPatchField onto a new patch
        processorFvPatchField
        (
            const processorFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        processorFvPatchField(const processorFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        processorFvPatchField
        (
            const processorFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~processorFvPatchField() = default;


    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //- Return true if this coupled patch fields are updated
            virtual bool couplesUpdated() const
            {
                return !pnf_.empty();
            }

            //- Return neighbour field of internal field
            //  Returning correct patch neighbour field data
            //  HJ, 10/Sep/2021
            virtual tmp<Field<Type> > patchNeighbourField() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            //  Note: fixed patchNeighbourField storage problem.  HJ, 8/Sep/2021
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            //  Note: fixed patchNeighbourField storage problem.  HJ, 8/Sep/2021
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                processorLduInterfaceField::transformCoupleField(f, cmpt);
            }


        // Matrix manipulation

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        // Block coupled interface functionality

//             //- Transform given patch field
//             virtual void transformCoupleField
//             (
//                 Field<Type>& f
//             ) const
//             {
//                 processorLduInterfaceField::transformCoupleField(f);
//             }

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        // Processor coupled interface functions

            //- Return processor number
            virtual int myProcNo() const
            {
                return procPatch_.myProcNo();
            }

            //- Return neigbour processor number
            virtual int neighbProcNo() const
            {
                return procPatch_.neighbProcNo();
            }


            // Communication support

                //- Is all data available
                virtual bool ready() const;

                //- Return communicator used for parallel communication
                virtual int comm() const
                {
                    return procPatch_.comm();
                }


            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(procPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procPatch_.forwardT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }

            //- Clear caches
            virtual void clearCaches() const;


    // Member operators

        // Note: all field operators need to handle cached patchNeighbourField
        // HJ, 8/Sep/2021

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const fvsPatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
