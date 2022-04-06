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
    testValueAt

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "FfowcsWilliamsHawkings.H"
#include "Farassat1AFormulation.H"
#include "fwhFormulation.H"

using namespace Foam;
using namespace functionObjects;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
TEST_CASE("Checks valueAt functionality in the libAcoustics", "[testValueAt]")
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
    for(label t=1; t<=10;t++)
    {
        data[0][0][0].first().append(0.1*t);
        data[0][0][0].second().append(t);
    }

    //Check
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,1.01) == 0);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,0.099) == 0);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,-1) == 0);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,1) == 10);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,0.1) == 1);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,0.15) == 1.5);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,0.55) == 5.5);
    REQUIRE (fwhFormulation1.valueAt(data,0,0,0,0.99) == 9.9);
}

// ************************************************************************* //
