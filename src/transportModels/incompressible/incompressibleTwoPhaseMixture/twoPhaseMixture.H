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
    Foam::twoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    twoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixture_H
#define twoPhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixture
:
    public transportModel
{
protected:

    // Protected data

        //- Name of phase 1
        word phase1Name_;

        //- Name of phase 1
        word phase2Name_;

        //- Viscosity model for phase 1
        autoPtr<viscosityModel> nuModel1_;

        //- Viscosity model for phase 1
        autoPtr<viscosityModel> nuModel2_;

        //- Density for phase 1
        dimensionedScalar rho1_;

        //- Density for phase 1
        dimensionedScalar rho2_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to phase flux field
        const surfaceScalarField& phi_;

        //- Reference to phase indicator field alpha
        const volScalarField& alpha1_;

        //- Mixture viscosity
        volScalarField nu_;


    // Private Member Functions

        //- Disallow copy construct
        twoPhaseMixture(const twoPhaseMixture&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const twoPhaseMixture&) = delete;


        //- Get phase name (backward compatibility)
        word getPhaseName(const word& key) const;

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        twoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1",
	    const dictionary& pyDict = dictionary()
        );


    // Destructor
    virtual ~twoPhaseMixture() = default;


    // Member Functions

        //- Return name for phase 1
        const word& phase1Name() const
        {
            return phase1Name_;
        }

        //- Return name for phase 2
        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the kinematic laminar viscosity for the mixture
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Return the dynamic laminar viscosity for the mixture
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        //  for the mixture
        tmp<surfaceScalarField> muf() const;


        //- Return the face-interpolated dynamic laminar viscosity
        //  for the mixture
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
