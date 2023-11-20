{
    // Bound h
    h.internalField() = max
    (
        max
        (
            h.internalField(),
            fac::average(max(h, h0))().internalField()
            *pos(h0.value() - h.internalField())
        ),
        h0.value()
    );

    h.correctBoundaryConditions();
}
