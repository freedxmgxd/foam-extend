{
    p.boundaryField().updateCoeffs();

    // Prepare clean 1/a_p without time derivative and under-relaxation
    // contribution
    rAU = 1.0/HUEqn.A();

    // Calculate U from convection-diffusion matrix
    U = rAU*HUEqn.H();

    // Consistently calculate flux
    pimple.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

    // Global flux balance
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        p.storePrevIter();

        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAU)/pimple.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }

        // Relax within the non-orthogonal loop
        p.relax();
    }

    // Explicitly relax pressure for momentum corrector except for last
    // corrector
    // if (!pimple.finalIter())
    // {
    //     p.relax();
    // }

#   include "movingMeshContinuityErrs.H"

    // Do not update velocity - SIMPLE!
    
    // Consistently reconstruct velocity after pressure equation. Note: flux is
    // made relative inside the function
    pimple.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);

#   include "limitU.H"
}
