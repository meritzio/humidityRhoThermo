Info<< "Reading thermophysical properties\n" << endl;

autoPtr<humidityRhoThermo> pThermo(humidityRhoThermo::New(mesh));
humidityRhoThermo& thermo = pThermo();

thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermoThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermoThermophysicalTransportModel::New(turbulence(), thermo)
);


#include "readGravitationalAcceleration.H"

Info << "\nReading hRef" << endl;

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, 0)
);

Info << "Calculating field g.h\n" << endl;
dimensionedScalar ghRef(- mag(g)*hRef);
volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell(0);
scalar pRefValue(0.0);
setRefCell(p, p_rgh, simple.dict(), pRefCell, pRefValue);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh - pRefValue;

mesh.schemes().setFluxRequired(p_rgh.name());

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createRadiationModel.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
