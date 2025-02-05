    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - fvc::div
        (
            muEff*
            (
                mesh.Sf()
              & fvc::interpolate(dev2(T(fvc::grad(U))))
            )
        )
    );

    UEqn.relax();

    solve
    (
        UEqn
     ==
      - gh*fvc::grad(rho)
      - fvc::grad(pd)
    );
