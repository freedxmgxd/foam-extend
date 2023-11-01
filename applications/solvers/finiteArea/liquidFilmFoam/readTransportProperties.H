    // Read transport properties

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Surface tension coefficient
    dimensionedScalar sigma
    (
        transportProperties.lookup("sigma")
    );

    // Density
    dimensionedScalar rhol
    (
        transportProperties.lookup("rhol")
    );

    dimensionedScalar h0
    (
        transportProperties.lookup("h0")
    );
