    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );
    mrfZones.addCoriolis(UEqn.ref());

    UEqn.ref().relax();

    solve
    (
        UEqn() == -fvc::grad(p)
    );
