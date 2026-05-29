Info<< "Reading/calculating face flux field phiB\n" << endl;

surfaceScalarField phiB
(
    IOobject
    (
        "phiB",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (fvc::interpolate(B) & mesh.Sf())
);
