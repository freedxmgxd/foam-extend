{
    // Update boundary velocity for consistency with the flux
    mrfZones.correctBoundaryVelocity(U);

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // HJ, explicit MRF
    mrfZones.addCoriolis(UEqn);
    
    // Add MRF and porous sources implicitly.  HJ, 18/Nov/2017
    tmp<volTensorField> tTU;

    // if (addMRF || addPorosity)
    // {
    //     tTU = tmp<volTensorField>
    //     (
    //         new volTensorField
    //         (
    //             IOobject
    //             (
    //                 "TU",
    //                 runTime.timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::NO_WRITE
    //             ),
    //             mesh,
    //             dimensionedTensor("zero", dimless/dimTime, tensor::zero)
    //         )
    //     );
    //     volTensorField& TU = tTU.ref();

    //     // Add implicit MRF source as a Hodge dual of the rotational velocity
    //     TU += *mrfZones.omega();

    //     // Add implicit resistance
    //     pZones.addResistance(UEqn, TU);

    //     trTU = inv(TU + tensor(I)*UEqn.A());
    //     trTU.ref().rename("rAU");
    // }
    // else
    {
        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

#   include "addImplicitMRFPorous.H"
#   include "addBlockCoupledBC.H"
}
