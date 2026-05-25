rho = thermo.rho();
rho.relax();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

const volScalarField rAU("rAU", 1.0/UEqn.A());
const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
);

bool adjustMass = mesh.schemes().steady() && adjustPhi(phiHbyA, U, p_rgh);

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf);

fvc::makeRelative(phiHbyA, rho, U);

fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p_rgh)*phiHbyA/fvc::interpolate(rho);

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA) + fvm::div(phid, p_rgh)
     ==
        fvModels.sourceProxy(rho, p_rgh)
    );

    while (pimple.correctNonOrthogonal())
    {
        p_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAUf, p_rgh);

        // Relax the pressure equation to ensure diagonal-dominance
        p_rghEqn.relax();

        p_rghEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        p_rghEqn.solve();
    }
}
else
{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
     ==
        fvModels.sourceProxy(rho, p_rgh)
    );

    while (pimple.correctNonOrthogonal())
    {
        p_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAUf, p_rgh);

        p_rghEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        p_rghEqn.solve();
    }
}

phi = phiHbyA + p_rghEqn.flux();

// Explicitly relax pressure for momentum corrector
p_rgh.relax();

p = p_rgh + rho*gh + pRef;

if (mesh.schemes().steady())
{
    #include "incompressible/continuityErrs.H"
}
else
{
    const bool constrained = fvConstraints.constrain(p);

    // Thermodynamic density update
    thermo.correctRho(psi*p - psip0);

    if (constrained)
    {
        rho = thermo.rho();
    }

    #include "rhoEqn.H"
}

// Correct the momentum source with the pressure gradient flux
// calculated from the relaxed pressure
U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
U.correctBoundaryConditions();
fvConstraints.constrain(U);
K = 0.5*magSqr(U);

if (mesh.schemes().steady())
{
    fvConstraints.constrain(p);
}

// For steady closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass && !thermo.incompressible())
{
    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    p_rgh = p - rho*gh - pRef;
    p.correctBoundaryConditions();
}

if (mesh.schemes().steady() || pimple.simpleRho())
{
    rho = thermo.rho();
}

if (mesh.schemes().steady() || pimple.simpleRho())
{
    rho.relax();
}

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
