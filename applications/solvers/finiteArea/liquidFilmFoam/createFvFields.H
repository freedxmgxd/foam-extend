   // Create fields for finite area post-processing

    volVectorField faUs
    (
        IOobject
        (
            "faUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    volScalarField faH
    (
        IOobject
        (
            "faH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimLength, 0)
    );

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh);
