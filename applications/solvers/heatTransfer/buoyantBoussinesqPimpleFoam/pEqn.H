{
    volScalarField rUA("rUA", 1.0/UEqn.A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn.H();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p_rgh);

    surfaceScalarField buoyancyPhi(rUAf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    phi -= buoyancyPhi;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            U -= rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rUAf);
            U.correctBoundaryConditions();
        }
    }

#   include "continuityErrs.H"
}
