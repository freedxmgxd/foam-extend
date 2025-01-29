{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(U) & mesh.Sf()
    );

    adjustPhi(phiU, U, pd);

    phi = phiU
      + (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    // Explicitly relax pressure
    pd.relax();

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

#   include "limitU.H"
}
