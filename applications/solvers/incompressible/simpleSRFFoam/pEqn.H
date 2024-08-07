    p.boundaryField().updateCoeffs();

    volScalarField rAUrel = 1.0/UrelEqn().A();
    Urel = rAUrel*UrelEqn().H();
    UrelEqn.clear();

#   include "limitUrel.H"

    // Calculate under-relaxation consistent flux
    simple.calcSteadyConsistentFlux(phi, Urel);
    adjustPhi(phi, Urel, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAUrel)/simple.aCoeff(Urel.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
          ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    scalar purf = mesh.solutionDict().fieldRelaxationFactor(p.name());

    Info<< "purf = " << purf << endl;
    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Reconstruct consistent velocity after pressure equation
    simple.reconstructSteadyVelocity(Urel, rAUrel, p);

#   include "limitUrel.H"
