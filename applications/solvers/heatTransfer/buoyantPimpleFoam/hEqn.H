{
    // Storing alphaH for post-processing function objects
    alphaH = turbulence->alphaEff();

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(alphaH, h)
     ==
        DpDt
    );

    if (thermalModelFluid.active())
    {
        Info<< "Adding thermal source: sum = "
            << gSum(thermalModelFluid.S()*mesh.V()) << " max = "
            << gMax(thermalModelFluid.S()().internalField())
            << endl;

        hEqn -= thermalModelFluid.S();
    }

    hEqn.relax();
    hEqn.solve();

    thermo.correct();
    rho = thermo.rho();
}
