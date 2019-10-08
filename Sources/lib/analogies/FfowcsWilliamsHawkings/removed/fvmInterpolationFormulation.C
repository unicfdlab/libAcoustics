/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvmInterpolationFormulation.H"
#include "FfowcsWilliamsHawkings.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IFstream.H"

/*

#include "addToRunTimeSelectionTable.H"

//sample surface
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"

#include "ListListOps.H"
#include "stringListOps.H"

#include "fvc.H"
#include "sampledPatch.H"
*/

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//namespace Foam
//{
//namespace functionObjects
//{
//    defineTypeNameAndDebug(FfowcsWilliamsHawkings, 0);
//    
//    addToRunTimeSelectionTable(functionObject, FfowcsWilliamsHawkings, dictionary);
//}
//}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fvmInterpolationFormulation::fvmInterpolationFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    intQdS_(0.0, fwh_.obr_.time().value()),
    intFdS_(0.0, fwh_.obr_.time().value())
{
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fvmInterpolationFormulation::~fvmInterpolationFormulation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::functionObjects::fvmInterpolationFormulation::initialize()
{
    
    intFdS_.resize(fwh_.observers_.size());
    intQdS_.resize(fwh_.observers_.size());
}

Foam::scalar Foam::functionObjects::fvmInterpolationFormulation::observerAcousticPressure(label iObs)
{
    intQdS_.value(iObs) = 0.0;
    intFdS_.value(iObs) = 0.0;

    scalar ct   = fwh_.obr_.time().value();
    scalar r  (0.0);
    scalar dS (0.0);
    vector n  (vector::zero);
    scalar Mr = 0.0;
    vector r_x (vector::zero);
    tensor Pf (tensor::zero);
    scalar oneByR2Mr(0.0);
    vector gradOneByRMr(vector::zero);
    scalar f1dotn(0.0);
    //scalar f2dotn(0.0);
    //vector divT(vector::zero);
    
    forAll(fwh_.controlSurfaces_, iSurf)
    {
        const sampledSurface& surf = fwh_.controlSurfaces_[iSurf];
        if (surf.interpolate())
        {
            Info<< "WARNING: Interpolation for surface " << surf.name() << " is on, turn it off"
                << endl;
        }

        const vectorField& Sf = surf.Sf();
        vectorField uS (fwh_.surfaceVelocity(surf)());
        scalarField rhoS (fwh_.surfaceDensity(surf)());
        scalarField pS (fwh_.surfacePressure(surf)() - fwh_.pInf_);
        //vectorField divTSurf (fwh_.surfaceStressDivergence(surf)());
        forAll(Sf, iFace)
        {
            //for Observer No. iObs
            {
                r  = magrobs_[iObs][iSurf][iFace];
                dS = mag(Sf[iFace]);
                n  = ni_[iSurf].value(iFace);
                r_x = robs_[iObs][iSurf][iFace]/r;
                Mr = mag((fwh_.vS_[iSurf][iFace]/fwh_.c0_) & r_x);
                qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
                qds_[iObs][iSurf][iFace].second().append
                (
                    (
                        (
                            (fwh_.rhoRef_*fwh_.vS_[iSurf][iFace] + rhoS[iFace]*(uS[iFace] - fwh_.vS_[iSurf][iFace]))
                            /
                            (r - r*Mr)
                        ) & n
                    )*dS
                );

                Pf = pS[iFace]*tensor::I + rhoS[iFace]*uS[iFace]*(uS[iFace] - fwh_.vS_[iSurf][iFace]);
                    
                oneByR2Mr = (1.0 / r / r)*(1.0/(1.0-Mr));
                
                gradOneByRMr = 
                (
                    -vector::one*oneByR2Mr
                    +oneByR2Mr*(fwh_.vS_[iSurf][iFace])/ fwh_.c0_ - (fwh_.vS_[iSurf][iFace] & r_x)*vector::one / fwh_.c0_
                );
                
                //divT = divTSurf[iFace];
                
                f1dotn = (Pf & gradOneByRMr) & n;
                //f2dotn = (divT * r / oneByR2Mr) & n;

                fds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
                fds_[iObs][iSurf][iFace].second().append //needs account for Doppler
                (
                    //(f1dotn + f2dotn) * dS
                    f1dotn * dS
                );
            } //for observers_
        } //for Sf

    } // for controlSurfaces_
    
    scalar ct1 = ct + fwh_.obr_.time().deltaT().value()*1.0e-6;
    //calculate acoustic pressure, zero if source didn't reached observer
    forAll(fwh_.controlSurfaces_, iSurf)
    {
        forAll(qds_[iObs][iSurf], iFace)
        {
            intQdS_.value(iObs) += 
                valueAt(qds_, iObs, iSurf, iFace, ct1);
            intFdS_.value(iObs) += 
                valueAt(fds_, iObs, iSurf, iFace, ct1);
        }
    }
    reduce (intQdS_.value(iObs), sumOp<scalar>());
    reduce (intFdS_.value(iObs), sumOp<scalar>());
    scalar coeff1 = 1. / 4. / Foam::constant::mathematical::pi;
    return (intQdS_.dot(ct,iObs) - intFdS_.value(iObs))*coeff1;
}

void Foam::functionObjects::fvmInterpolationFormulation::update()
{
    fwhFormulation::update();
}

void Foam::functionObjects::fvmInterpolationFormulation::clearExpiredData()
{
    fwhFormulation::clearExpiredData();
}


// ************************************************************************* //

//
//END OF FILE
//



