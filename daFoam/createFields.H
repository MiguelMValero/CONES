Info << "Reading field T\n" << endl;
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

Info << "Reading transportProperties\n" << endl;
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

Info << "Reading diffusivity DT\n" << endl;
dimensionedScalar nu("nu",  dimensionSet(0, 2, -1, 0, 0, 0, 0), transportProperties);
dimensionedScalar DT("DT",  dimensionSet(0, 2, -1, 0, 0, 0, 0), transportProperties);

#include "createPhi.H"
