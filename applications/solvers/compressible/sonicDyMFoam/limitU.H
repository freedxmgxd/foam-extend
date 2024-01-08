{
    // Bound the velocity
    volScalarField magU = mag(U);

    if (max(magU) > UMax)
    {
        Info<< "magU " << max(magU).value()
            << " (" << gMax(magU.internalField()) << "): "
            << "Limiting."
            << endl;

        volScalarField Ulimiter =
            pos(magU - UMax)*UMax/(magU + smallU)
          + neg(magU - UMax);

        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}
