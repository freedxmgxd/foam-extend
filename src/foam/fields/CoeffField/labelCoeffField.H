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
    CoeffField<T>

Description
    Template specialisation for label coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    tensorCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef labelCoeffField_H
#define labelCoeffField_H

#include "CoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
Ostream& operator<<(Ostream&, const CoeffField<label>&);
Ostream& operator<<(Ostream&, const tmp<CoeffField<label> >&);


/*---------------------------------------------------------------------------*\
                     Class labelCoeffField Declaration
\*---------------------------------------------------------------------------*/

template<>
class CoeffField<label>
:
    public labelField
{
public:

    // Public data types

        //- Component type
        typedef label scalarType;
        typedef label linearType;
        typedef label squareType;


        //- Field type
        typedef labelField scalarTypeField;
        typedef labelField linearTypeField;
        typedef labelField squareTypeField;


    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct given size
        explicit CoeffField<label>(const label size)
        :
            labelField(size, pTraits<label>::zero)
        {}

        //- Construct as copy
        CoeffField(const CoeffField<label>& f)
        :
            labelField(f)
        {}

        //- Construct from Istream
        CoeffField(Istream& is)
        :
            labelField(is)
        {}


    // Member functions

        //- Return active type
        blockCoeffBase::activeLevel activeType() const;

        //- Return the field transpose
        tmp<CoeffField<label> > transpose() const;


        // Return as typed.  Fails when asked for the incorrect type

            //- Return as scalar field
            const labelField& asScalar() const;
            labelField& asScalar();

            //- Return as linear field.  Equal to scalar field
            const labelField& asLinear() const;
            labelField& asLinear();

            //- Return as square field.  Equal to scalar field
            const labelField& asSquare() const;
            labelField& asSquare();


        //- Return coefficient as a block
        BlockCoeff<label> getCoeff(const label index) const;

        //- Set coefficient from a block
        void setCoeff(const label index, const BlockCoeff<label>& coeff);


        // Subset operations

            //- Get subset with offset and size and store in given field
            void getSubset
            (
                CoeffField<label>& f,
                const label start,
                const label size
            ) const;

            //- Get subset with addressing and store in given field
            void getSubset
            (
                CoeffField<label>& f,
                const labelList& addr
            ) const;


            //- Set subset with offset and size from given field
            void setSubset
            (
                const CoeffField<label>& f,
                const label start,
                const label size
            );

            //- Get subset with addressing and store in target field
            void setSubset
            (
                const CoeffField<label>& f,
                const labelList& addr
            );

            //- Zero out subset with offset and size
            void zeroOutSubset
            (
                const label start,
                const label size
            );

            //- Zero out subset with addressing
            void zeroOutSubset
            (
                const labelList& addr
            );

            //- Add subset with addressing to field
            void addSubset
            (
                const CoeffField<label>& f,
                const labelList& addr
            );


    // Member operators

        void operator=(const CoeffField<label>&);
        void operator=(const labelField&);
        void operator=(const tmp<labelField>&);


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const CoeffField<label>&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const tmp<CoeffField<label> >&
        );
};


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

template<>
tmp<CoeffField<label> > inv(const CoeffField<label>& f);


template<>
void negate
(
    CoeffField<label>& f,
    const CoeffField<label>& f1
);


template<>
void multiply
(
    labelField& f,
    const CoeffField<label>& f1,
    const label& f2
);


template<>
void multiply
(
    labelField& f,
    const CoeffField<label>& f1,
    const labelField& f2
);


template<>
void multiply
(
    labelField& f,
    const labelField& f1,
    const CoeffField<label>& f2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
