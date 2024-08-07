{
    scalar limitMagUrel = readScalar(simple.dict().lookup("limitMagUrel"));

    volScalarField magUrel(mag(Urel));

    scalar maxMagUrel = max(magUrel).value();

    Info<< "mag(Urel): max: " << maxMagUrel
        << " avg: " << magUrel.weightedAverage(mesh.V()).value();

    if (maxMagUrel > limitMagUrel)
    {
        Urel.internalField() *=
            neg(magUrel.internalField() - limitMagUrel)
          + pos(magUrel.internalField() - limitMagUrel)*
            limitMagUrel/(magUrel.internalField() + SMALL);
        Urel.correctBoundaryConditions();
        Info << " ...limiting" << endl;
    }
    else
    {
        Info << endl;
    }
}
