Info<< "Reading thermophysical properties\n" << endl;

autoPtr<humidityRhoThermo> pThermo
(
    humidityRhoThermo::New(mesh)
);
humidityRhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

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

Info<< "\nReading pRef" << endl;

uniformDimensionedScalarField pRef
(
    IOobject
    (
        "pRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimPressure, 0)
);

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

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.schemes().setFluxRequired(p_rgh.name());

hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.name(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

if (!mesh.foundObject<IOdictionary>("radiationProperties"))
{
    typeIOobject<IOdictionary> radiationProperties
    (
        "radiationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (radiationProperties.headerOk())
    {
        const word modelType
        (
            IOdictionary(radiationProperties).lookup("radiationModel")
        );

        if (modelType != "none")
        {
            FatalErrorInFunction
                << "Radiation model " << modelType
                << " selected but not enabled in fvModels" << nl
                << "To enable radiation add " << nl << nl
                << "radiation" << nl
                << "{" << nl
                << "    type    radiation;" << nl
                << "    libs    (\"libradiationModels.so\");" << nl
                << "}" << nl << nl
                << "to " << radiationProperties.relativePath()/"fvModels"
                << exit(FatalError);
        }
    }
}

