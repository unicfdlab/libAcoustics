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
    testFormulationFunction

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "FfowcsWilliamsHawkings.H"
#include "Farassat1AFormulation.H"
#include "fwhFormulation.H"

using namespace Foam;
using namespace functionObjects;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
TEST_CASE("Checks expired functionality in the libAcoustics", "[testFormulationFunction]")
{
    #include "createTestTime.H"
    #include "createTestMesh.H"
    #include "createTestFwhDict.H"

    positionDict.add<vector>("position",vector(0,10,0), true);
    observersDict.add("TestPoint1",positionDict, true);
    fwhControlDict.add("observers", observersDict, true);
    surfaces.resize(1);
    surfaces[0] = "testSurface{type sampledTriSurfaceMesh; surface line1.stl; source cells;}";
    fwhControlDict.add("surfaces",surfaces, true);
    fwhControlDict.add("cleanFreq", 1);

    FfowcsWilliamsHawkings fwh("fwhName",timeObj,fwhControlDict);
    Farassat1AFormulation fwhFormulation1(fwh);

    vectorField Sf; 
    vectorField uS;
    vectorField Cf;
    scalarField pS;
    scalarField rhoS;
    scalarField magSf;
    
    label iSurf = 0;
    label iObs = 0;
    
    Sf = fwh.getSampledSurface()[iSurf].Sf();
    REQUIRE (Sf.size() == 100);

    magSf = fwh.getSampledSurface()[iSurf].magSf();
    Cf = fwh.getSampledSurface()[iSurf].Cf();
    uS.setSize(Sf.size());
    pS.setSize(Sf.size());
    rhoS.setSize(Sf.size());    
    
    SECTION ("Test without clear function")
    {
        while (timeObj.run())
        {
            timeObj++;
            scalar t = timeObj.value();
            Info << "TIME = " << t << endl;
            fwhFormulation1.update();

            forAll(Sf,iFace)
            {
                rhoS[iFace] = 1;
                pS[iFace] = 0;
                vector normals = Sf[iFace]/magSf[iFace];
                uS[iFace] = 0*normals;
            }
            fwhFormulation1.observerAcousticPressure(Sf, uS, rhoS, pS, iObs, iSurf, t);
        }
        REQUIRE (fwhFormulation1.getQdsData()[0][0][0].first().size() == 10);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][0].second().size() == 10);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][99].first().size() == 10);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][99].second().size() == 10);
    }

    SECTION ("Test with clear function")
    {
        timeObj.setTime(0, 0);
        while (timeObj.run())
        {
            timeObj++;
            scalar t = timeObj.value();
            Info << "TIME = " << t << endl;
            fwhFormulation1.update();
            forAll(Sf,iFace)
            {
                rhoS[iFace] = 1;
                pS[iFace] = 0;
                vector normals = Sf[iFace]/magSf[iFace];
                uS[iFace] = 0*normals;
            }
            fwhFormulation1.observerAcousticPressure(Sf, uS, rhoS, pS, iObs, iSurf, t);
            fwhFormulation1.clearExpiredData();
        }
        REQUIRE (fwhFormulation1.getQdsData()[0][0][0].first().size() == 3);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][0].second().size() == 3);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][99].first().size() == 3);
        REQUIRE (fwhFormulation1.getQdsData()[0][0][99].second().size() == 3);
    }

}

// ************************************************************************* //
