    dictionary fwhControlDict {};
    dictionary observersDict{};
    dictionary positionDict{};
    wordList surfaces;

    fwhControlDict.add<word>("type", "FfowcsWilliamsHawkings");
    fwhControlDict.add<word>("patches","()");
    fwhControlDict.add<word>("formulationType","Farassat1AFormulation");
    fwhControlDict.add<word>("interpolationScheme","cell");
    
    fwhControlDict.add<vector>("CofR",vector::zero);
    fwhControlDict.add<scalar>("timeStart",0.0);
    fwhControlDict.add<scalar>("timeEnd",1.0);
    fwhControlDict.add<scalar>("c0",100);
    fwhControlDict.add<scalar>("rhoInf",1.2);
    fwhControlDict.add<scalar>("pInf",0);
    
    positionDict.add<vector>("position",vector(1,0,0));
    observersDict.add("TestPoint1",positionDict);
    fwhControlDict.add("observers", observersDict);

    fwhControlDict.add("surfaces",surfaces, true);
