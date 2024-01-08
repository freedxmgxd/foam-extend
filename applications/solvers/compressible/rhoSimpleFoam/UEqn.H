    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.ref().relax();

    solve(UEqn() == -fvc::grad(p));
