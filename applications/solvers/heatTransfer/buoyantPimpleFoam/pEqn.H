{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));

    U = rAU*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        )
    );

    phi = phiU - ghf*fvc::snGrad(rho)*rhorAUf*mesh.magSf();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::ddt(psi, p_rgh)
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.solve
        (
            mesh.solutionDict().solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi += p_rghEqn.flux();
        }
    }

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

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

#   include "rhoEqn.H"
#   include "compressibleContinuityErrs.H"

    rho = thermo.rho();

    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}
