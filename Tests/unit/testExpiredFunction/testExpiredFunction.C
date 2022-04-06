/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    testExpiredFunction

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "FfowcsWilliamsHawkings.H"
#include "Farassat1AFormulation.H"
#include "fwhFormulation.H"

using namespace Foam;
using namespace functionObjects;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
TEST_CASE("Checks expired functionality in the libAcoustics", "[testExpiredFunction]")
{
    #include "createTestTime.H"
    #include "createTestMesh.H"
    #include "createTestFwhDict.H"

    FfowcsWilliamsHawkings fwh("fwhName",timeObj,fwhControlDict);
    Farassat1AFormulation fwhFormulation1(fwh);

    //Create test data
    typedef Pair<DynamicList<scalar> > pointTimeData;
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    surfaceTimeData data;

    data.resize(1);
    data[0].resize(1);
    data[0][0].resize(1);

    SECTION ("Test boundary value")
    {
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 10) == -1);
    
        for(label t=1; t<=1;t++)
        {
            data[0][0][0].first().append(0.1*t);
            data[0][0][0].second().append(t);
        }
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 10) == 0);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.05) == -1);
    }

    SECTION ("Test findExpired function")
    {
        for(label t=1; t<=10;t++)
        {
            data[0][0][0].first().append(0.1*t);
            data[0][0][0].second().append(t);
        }

        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.0) == -1);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.1) == -1);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.1001) == 0);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.2) == 0);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.25) == 1);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.5) == 3);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.55) == 4);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 0.95) == 8);
        REQUIRE (fwhFormulation1.findExpiredIndex(data[0][0][0], 1.1) == -1);
    }

    SECTION ("Test clear expired data")
    {
        for(label t=1; t<=10;t++)
        {
            data[0][0][0].first().append(0.1*t);
            data[0][0][0].second().append(t);
        }

        data[0][0][0].operator=(fwhFormulation1.getNewPointData(data[0][0][0], fwhFormulation1.findExpiredIndex(data[0][0][0], 0.1)));
        REQUIRE (data[0][0][0].first().size() == 10);
        REQUIRE (data[0][0][0].second().size() == 10);
        data[0][0][0].operator=(fwhFormulation1.getNewPointData(data[0][0][0], fwhFormulation1.findExpiredIndex(data[0][0][0], 0.55)));
        REQUIRE (data[0][0][0].first().size() == 5);
        REQUIRE (data[0][0][0].second().size() == 5);
        data[0][0][0].operator=(fwhFormulation1.getNewPointData(data[0][0][0], fwhFormulation1.findExpiredIndex(data[0][0][0], 1.95)));
        REQUIRE (data[0][0][0].first().size() == 5);
        REQUIRE (data[0][0][0].second().size() == 5);

        scalar mean = 0;
        forAll (data[0][0][0].second(), i)
        {
            mean += data[0][0][0].second()[i];
        }
        REQUIRE (mean/data[0][0][0].second().size() == 8);
    }
}

// ************************************************************************* //
