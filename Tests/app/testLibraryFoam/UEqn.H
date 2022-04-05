// Solve the Momentum equation

fvVectorMatrix UEqn
(
    rho*fvm::ddt(U)
    ==
    -fvc::grad(pa)
);

if (solvePrimaryRegion)
{
    UEqn.relax();
    UEqn.solve();
}
else
{
    U.correctBoundaryConditions();
}