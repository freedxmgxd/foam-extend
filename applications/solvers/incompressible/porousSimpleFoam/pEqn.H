if (pressureImplicitPorosity)
{
    U = (trTU() & UEqn().H());
}
else
{
    U = trAU()*UEqn().H();
}

UEqn.clear();
phi = fvc::interpolate(U) & mesh.Sf();
adjustPhi(phi, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn.ref() = (fvm::laplacian(trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn.ref() = (fvm::laplacian(trAU(), p) == fvc::div(phi));
    }

    tpEqn.ref().setReference(pRefCell, pRefValue);

    tpEqn.ref().solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi -= tpEqn().flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();
