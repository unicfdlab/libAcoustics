Info<< "\nReading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar c0("c0", dimVelocity, transportProperties);

dimensionedScalar rho("rho", dimDensity, transportProperties);

scalar MaxCo =
    max(mesh.surfaceInterpolation::deltaCoeffs()*c0).value()
   *runTime.deltaT().value();

Info<< "Max acoustic Courant Number = " << MaxCo << endl;
