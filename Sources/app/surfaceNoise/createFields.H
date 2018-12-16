    word dictName("surfaceNoiseDict");

    if (args.optionFound("dict"))
    {
        dictName = args["dict"];
    }

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // upper frequency of frequency band
    scalar maxFreq = dict.lookupOrDefault("maxFrequency", 10000.0);

    // output format
    word outputFormat = dict.lookupOrDefault<word>("outputFormat", "gmsh");

    // output directory (similar with surface name)
    fileName inputFileName = dict.lookup("inputFileName");

    //interpolation scheme
    //word interpolationScheme = dict.lookup("interpolationScheme");
