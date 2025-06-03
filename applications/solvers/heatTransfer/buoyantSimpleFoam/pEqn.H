{
    volScalarField rAU("rAU", 1.0/UEqn().A());
    surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));

    U = rAU*UEqn().H();
    UEqn.clear();

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        )
    );

    phi = phiU - ghf*fvc::snGrad(rho)*rhorAUf*mesh.magSf();

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();
        }
    }

    // Explicitly relax pressure for momentum corrector
    p_rgh.relax();

    U += rAU*fvc::reconstruct((phi - phiU)/rhorAUf);
    U.correctBoundaryConditions();

    // Update pressure
    p = p_rgh + rho*gh;
    p.correctBoundaryConditions();

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity

    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();

        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }

#   include "continuityErrs.H"

    rho = thermo.rho();
    rho.relax();

    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}
