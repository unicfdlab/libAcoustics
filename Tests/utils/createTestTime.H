    dictionary testDict {};
    testDict.add<scalar>("deltaT", 0.1);
    testDict.add<scalar>("startTime", 0.0);
    testDict.add<scalar>("endTime", 1);
    testDict.add<scalar>("writeInterval", 100000); //no write 
    testDict.add<scalar>("writeFrequency",1);
    Foam::Time timeObj(testDict, ".", ".");
    