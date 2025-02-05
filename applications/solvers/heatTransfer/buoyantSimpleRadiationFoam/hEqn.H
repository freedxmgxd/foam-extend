{
    // Storing alphaH for post-processing function objects
    alphaH = turbulence->alphaEff();

    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(alphaH, h)
     ==
        fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p))
      - p*fvc::div(phi/fvc::interpolate(rho))
      + radiation->Sh(thermo)
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

    radiation->correct();
}
