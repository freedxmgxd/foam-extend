{
    // Bound the velocity
    volScalarField magU = mag(U);
    dimensionedScalar U1 = max(magU);

    if (U1 > UMax)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        volScalarField Ulimiter = pos(magU - UMax)*UMax/(magU + smallU)
            + neg(magU - UMax);
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}
