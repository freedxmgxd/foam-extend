    // Solve the Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax();

    volScalarField rUA = 1.0/UEqn.A();

    solve(UEqn == -fvc::grad(p));
