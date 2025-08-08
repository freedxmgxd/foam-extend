    Info<< "Reading transportProperties\n" << endl;
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

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

#   include "readTransportProperties.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_rgh, simple.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p_rgh.name());

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - beta*(T - TRef)
    );

    // Create thermal model for fluid sources
    thermalSourceModel thermalModelFluid(laminarTransport, T);

    // Create porous zones
    porousZones pZones(mesh);

    // Thermal flux
    surfaceScalarField rhoCpPhi
    (
        "rhoCpPhi",
        rho*Cp*phi
    );

    // Create thermal diffusivity field (post-processing)
    volScalarField kappaEff
    (
        "kappaEff",
        rho*Cp*(turbulence->nu()/Pr + turbulence->nut()/Prt)
    );
