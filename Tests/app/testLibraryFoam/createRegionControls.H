fvSolution solutionDict(runTime);

const dictionary& pimpleDict = solutionDict.subDict("PIMPLE");

bool solvePrimaryRegion
(
    pimpleDict.getOrDefault("solvePrimaryRegion", true)
);
