    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn().relax();

    solve
    (
        UEqn()
     ==
      - fvc::reconstruct(ghf*fvc::snGrad(rho)*mesh.magSf())
      - fvc::grad(p_rgh)
    );
