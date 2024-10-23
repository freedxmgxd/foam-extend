    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicRhoThermo> pThermo
    (
        basicRhoThermo::New(mesh)
    );
    basicRhoThermo& thermo = pThermo();

    // Correct thermo for consistency
    thermo.correct();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    rho.correctBoundaryConditions();

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();

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

#   include "compressibleCreatePhi.H"

    // Read field p_rgh
    // Note: since initialisation for p_rgh is much easier than for p
    // p_rgh shall be considered as the base pressure field
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

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    // Force p to be consistent with p_rgh
    p = p_rgh + rho*gh;
    p.correctBoundaryConditions();

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    mesh.schemesDict().setFluxRequired(p_rgh.name());

    // Record initial mass
    Switch closedVolume(simple.dict().lookup("closedVolume"));

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    // Create thermal model for fluid sources
    thermalSourceModel thermalModelFluid(thermo, thermo.T());

    // Storing alphaEff for post-processing function objects
    volScalarField alphaH
    (
        "alphaH",
        turbulence->alphaEff()
    );
