    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    pZones.addResistance(UEqn.ref());

    UEqn.ref().relax();

    solve
    (
        UEqn()
      ==
        fvc::reconstruct
        (
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        )
    );
