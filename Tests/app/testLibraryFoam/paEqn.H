
fvScalarMatrix paEqn
(
    fvm::d2dt2(pa) - sqr(c0)*fvc::laplacian(pa)
);

if (solvePrimaryRegion)
{
    paEqn.relax();
    paEqn.solve();
}
else
{
    pa.correctBoundaryConditions();
}
