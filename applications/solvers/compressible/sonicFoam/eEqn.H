{
    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - p*fvc::div(phi/fvc::interpolate(rho))
        // viscous heating?
    );

    eEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor(e.select(pimple.finalIter()))
    );

    eEqn.solve();

    thermo.correct();

    // Recalculate compressibility
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();

    // Recalculate density
    rho = thermo.rho();
    rho.correctBoundaryConditions();
}
