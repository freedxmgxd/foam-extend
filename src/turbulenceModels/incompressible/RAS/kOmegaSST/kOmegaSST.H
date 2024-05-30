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
    Foam::incompressible::RASModels::kOmegaSST

Description
    Implementation of the k-omega-SST turbulence model for incompressible
    flows.

    Turbulence model described in:
    @verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001.
    @endverbatim

    with updated coefficients from
    @verbatim
        Menter, F. R., Kuntz, M., and Langtry, R.,
        "Ten Years of Industrial Experience with the SST Turbulence Model",
        Turbulence, Heat and Mass Transfer 4, 2003,
        pp. 625 - 632.
    @endverbatim

    but with the consistent production terms from the 2001 paper as form in the
    2003 paper is a typo, see
    @verbatim
        http://turbmodels.larc.nasa.gov/sst.html
    @endverbatim

    and the addition of the optional F3 term for rough walls from
    \verbatim
        Hellsten, A.
        "Some Improvements in Menter’s k-omega-SST turbulence model"
        29th AIAA Fluid Dynamics Conference,
        AIAA-98-2554,
        June 1998.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficuients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been corrected.

    The default model coefficients correspond to the following:
    @verbatim
        kOmegaSSTCoeffs
        {
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;
        }
    @endverbatim

SourceFiles
    kOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSST_H
#define kOmegaSST_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmega Declaration
\*---------------------------------------------------------------------------*/

class kOmegaSST
:
    public RASModel
{
    // Private data

        // Model coefficients

            //- Coefficient alphaK1
            dimensionedScalar alphaK1_;

            //- Coefficient alphaK1
            dimensionedScalar alphaK2_;

            //- Coefficient alphaOmega1
            dimensionedScalar alphaOmega1_;

            //- Coefficient
            dimensionedScalar alphaOmega2_;

            //- Coefficient gamma1
            dimensionedScalar gamma1_;

            //- Coefficient gamma2
            dimensionedScalar gamma2_;

            //- Coefficient beta1
            dimensionedScalar beta1_;

            //- Coefficient beta2
            dimensionedScalar beta2_;

            //- Coefficient betaStar
            dimensionedScalar betaStar_;

            //- Coefficient a1
            dimensionedScalar a1_;

            //- Coefficient b1
            dimensionedScalar b1_;

            //- Coefficient c1
            dimensionedScalar c1_;

            //- F3 switch
            Switch F3_;


        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
        wallDist y_;


        // Fields

            //- Turbulence kinetic energy field
            volScalarField k_;

            //- Eddy turnover time field (omega = epsilon/(betaStar k)), betaStar = Cmu
            volScalarField omega_;

            //- Turbulence viscosity
            volScalarField nut_;


    // Private member functions

        //- Disallow default bitwise copy construct
        kOmegaSST(const kOmegaSST&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const kOmegaSST&) = delete;


        //- Function F1
        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Function F1
        tmp<volScalarField> F2() const;

        //- Function F1
        tmp<volScalarField> F3() const;

        //- Function F1
        tmp<volScalarField> F23() const;

        //- Blending function (ke-kw)
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        //- Return F1 alphaK blend field
        tmp<volScalarField> alphaK
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        //- Return F1 alphaOmega blend field
        tmp<volScalarField> alphaOmega
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        //- Return beta blend field
        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, beta1_, beta2_);
        }

        //- Return gamma blend field
        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, gamma1_, gamma2_);
        }


public:

    //- Runtime type information
    TypeName("kOmegaSST");


    // Constructors

        //- Construct from components
        kOmegaSST
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kOmegaSST() = default;


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK(F1)*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega(F1)*nut_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
