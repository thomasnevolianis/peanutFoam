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

// Mobility
dimensionedScalar M
(
    transportProperties.lookup("M")
);
// Interface thickness
dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);
// dynamic viscosity of fluid 1
dimensionedScalar mu1
(
    transportProperties.lookup("mu1")
);
// dynamic viscosity of fluid 2
dimensionedScalar mu2
(
    transportProperties.lookup("mu2")
);
// density of fluid 1
dimensionedScalar rho1
(
    transportProperties.lookup("rho1")
);
// density of fluid 2
dimensionedScalar rho2
(
    transportProperties.lookup("rho2")
);
// bending stiffness
dimensionedScalar k
(
    transportProperties.lookup("k")
);
// stretching coefficient
dimensionedScalar k_P
(
    transportProperties.lookup("k_P")
);

// order parameter
Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    (rho1*(1+C))/2 + (rho2*(1-C))/2
);
// dynamic viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    (mu1*(1+C))/2 + (mu2*(1-C))/2
);
// pressure
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// velocity
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

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

