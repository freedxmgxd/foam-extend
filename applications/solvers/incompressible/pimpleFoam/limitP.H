{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin || p2 > pMax)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p.max(pMin);
        p.min(pMax);
        p.correctBoundaryConditions();
    }
}
