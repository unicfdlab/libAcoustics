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
TEST_CASE("Checks functionality in the libAcoustics", "[testMonopoleSource]")
{
    #include "createTestTime.H"
    #include "createTestMesh.H"
    #include "createTestFwhDict.H"

    timeObj.setDeltaT(1e-4, true);
    timeObj.setEndTime(0.2);
    positionDict.add<vector>("position",vector(1,0,0), true);
    observersDict.add("TestPoint1",positionDict, true);
    positionDict.add<vector>("position",vector(-1,0,0), true);
    observersDict.add("TestPoint2",positionDict);
    positionDict.add<vector>("position",vector(0,1,0), true);
    observersDict.add("TestPoint3",positionDict);
    positionDict.add<vector>("position",vector(0,-1,0), true);
    observersDict.add("TestPoint4",positionDict);
    positionDict.add<vector>("position",vector(0.7,0.7,0), true);
    observersDict.add("TestPoint5",positionDict);
    positionDict.add<vector>("position",vector(0.7,-0.7,0), true);
    observersDict.add("TestPoint6",positionDict);
    
    fwhControlDict.add("observers", observersDict, true);
    fwhControlDict.set("c0", 343);
    surfaces.resize(1);
    surfaces[0] = "testSurface{type sampledTriSurfaceMesh; surface sphere_r0.2.stl; source cells;}";
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
    
    scalar c0;
    fwhControlDict.lookup("c0") >> c0;

    scalar U0 = 1;
    scalar rho = 1.2;
    scalar freq = 100;
    scalar r0 = 0.2;
    scalar k = freq*2.0*constant::mathematical::pi/c0;
    
    Sf = fwh.getSampledSurface()[iSurf].Sf();
    REQUIRE (Sf.size() == 2760);

    magSf = fwh.getSampledSurface()[iSurf].magSf();
    Cf = fwh.getSampledSurface()[iSurf].Cf();
    uS.setSize(Sf.size());
    pS.setSize(Sf.size());
    rhoS.setSize(Sf.size()); 

    REQUIRE (fwh.getSoundObservers().size() == 6);
    List<scalarField> oap;
    List<scalarField> analyticResult;
    oap.setSize(fwh.getSoundObservers().size());
    analyticResult.setSize(fwh.getSoundObservers().size());
    timeObj.setTime(0, 0);
    
    SECTION ("Monopole")
    {
        while (timeObj.run())
        {
            timeObj++;
            scalar t = timeObj.value();
            forAll(fwh.getSoundObservers(), iObs)
            {
                scalar rSource = mag(fwh.getSoundObservers()[iObs].position());
                analyticResult[iObs].append(r0*r0*c0*k*U0*rho*(r0*k*Foam::cos(k*(r0-rSource+c0*t))-Foam::sin(k*(r0-rSource+c0*t))) / (rSource+r0*r0*k*k*rSource));
            }
        }
        //Output max/min pressure
        forAll(fwh.getSoundObservers(), iObs)
        {
            Info << "Max pressure for " << iObs << " observer = " << max(analyticResult[iObs]) << endl;
            Info << "Min pressure for " << iObs << " observer = " << min(analyticResult[iObs]) << endl;
        }
        timeObj.setTime(0, 0);
        while (timeObj.run())
        {
            timeObj++;
            scalar t = timeObj.value();
            Info << "TIME = " << t << endl;
            fwhFormulation1.update();            

            forAll(fwh.getSoundObservers(), iObs)
            {
                forAll(Sf,iFace)
                {
                    rhoS[iFace] = rho;
                    pS[iFace] = r0*r0*c0*k*U0*rho*(r0*k*Foam::cos(k*(r0-r0 + c0*t)) - Foam::sin(k*(r0-r0 + c0*t))) / (r0+r0*r0*k*k*r0);
                    vector normals = Sf[iFace]/magSf[iFace];
                    uS[iFace] = r0*r0*U0*((1.0 + k*k*r0*r0)*Foam::cos(k*(r0-r0+c0*t)) + k*(r0-r0)*Foam::sin(k*(r0-r0+c0*t))) / ((1.0 + r0*r0*k*k)*r0*r0)*normals;
                }
                oap[iObs].append(fwhFormulation1.observerAcousticPressure(Sf, uS, rhoS, pS, iObs, iSurf, t));
            }
        }
        forAll(fwh.getSoundObservers(), iObs)
        {
            Info << "For observer " << iObs << " MinOAP = " << min(analyticResult[iObs]) << " MaxOAP = " << max(analyticResult[iObs]) << endl; 
            forAll(oap[iObs], pI)
            {
                if (pI >= 1000)
                {
                    REQUIRE ((oap[iObs][pI]-analyticResult[iObs][pI]) <= 0.055);
                }
            }
        }
    }
}

// ************************************************************************* //
