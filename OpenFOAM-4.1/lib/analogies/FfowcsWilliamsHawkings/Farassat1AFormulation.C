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

#include "Farassat1AFormulation.H"
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

//Foam::scalar Foam::functionObjects::FfowcsWilliamsHawkings::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Farassat1AFormulation::Farassat1AFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    Un_(0),
    Lr_(0),
    Mr_(0),

    intDotQdS_(0.0, fwh_.obr_.time().value()),
    intFdS_(0.0, fwh_.obr_.time().value())
{
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Farassat1AFormulation::~Farassat1AFormulation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::functionObjects::Farassat1AFormulation::initialize()
{
    intFdS_.resize(fwh_.observers_.size());
    intDotQdS_.resize(fwh_.observers_.size());
    
    Lr_.resize(fwh_.observers_.size());
    Mr_.resize(fwh_.observers_.size());
    Un_.resize(fwh_.observers_.size());

    forAll(Lr_, iObs)
    {
        Lr_[iObs].resize(fwh_.controlSurfaces_.size());
        Mr_[iObs].resize(fwh_.controlSurfaces_.size());
        Un_[iObs].resize(fwh_.controlSurfaces_.size());

        forAll(Lr_[iObs], iSurf)
        {
            Lr_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            Mr_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            Un_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
        }
    }
}

Foam::scalar Foam::functionObjects::Farassat1AFormulation::observerAcousticPressure(label iObs)
{
    scalar ct = fwh_.obr_.time().value();
    
        //Farassat 1A
    vector L (vector::zero);
    scalar lr (0.0);
    scalar lM (0.0);
    scalar dotlr (0.0);
    vector r  (vector::zero);
    vector rh (vector::zero);
    vector n  (vector::zero);
    scalar dS (0.0);
    scalar magr(0.0);
    vector M  (vector::zero);
    scalar magM (0.0);
    scalar Mr (0.0);
    scalar dotMr (0.0);
    tensor Pf (tensor::zero);
    scalar OneByOneMr(0.0);
    scalar OneByOneMrSq(0.0);
    
    scalar fpart1 (0.0);
    scalar fpart2 (0.0);
    scalar fpart3 (0.0);
    vector U(vector::zero);
    scalar Un(0.0);
    scalar dotUn(0.0);
    vector dotn(vector::zero);

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
        
        //Farassat 1A formulation
        forAll(Sf, iFace)
        {
        
            //For observe No iObs
            {
                r = robs_[iObs][iSurf][iFace];
                magr = magrobs_[iObs][iSurf][iFace];
                rh = r / magr;
                dS = mag(Sf[iFace]);
                n = ni_[iSurf].value(iFace);
    
                M = fwh_.vS_[iSurf][iFace] / fwh_.c0_;
                Mr = M & rh;
                magM = mag(M);

                U = (1.0 - rhoS[iFace] / fwh_.rhoRef_) * fwh_.vS_[iSurf][iFace]
            	    + rhoS[iFace] * uS[iFace] / fwh_.rhoRef_;
		Pf = pS[iFace]*tensor::I + rhoS[iFace]*uS[iFace]*(uS[iFace] - fwh_.vS_[iSurf][iFace]);

                L = Pf & n;
                lM = L & M; 
		lr = L & rh;
                Mr = M & rh;
                Un = U & n;
            	
            	OneByOneMr = 1.0 / (1.0 - Mr);
                OneByOneMrSq = OneByOneMr*OneByOneMr;
                
                Un_[iObs][iSurf].value(iFace) = Un;
                Lr_[iObs][iSurf].value(iFace) = lr;
                Mr_[iObs][iSurf].value(iFace) = Mr; 
            
                dotlr = Lr_[iObs][iSurf].dot(ct, iFace);
                dotMr = Mr_[iObs][iSurf].dot(ct, iFace);
                dotUn = Un_[iObs][iSurf].dot(ct, iFace);
                dotn = ni_[iSurf].dot(ct, iFace);
                
                qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
                qds_[iObs][iSurf][iFace].second().append
                (
                    (
                        fwh_.rhoRef_ * (dotUn + (U & dotn)) * OneByOneMrSq / magr
                        +
                        fwh_.rhoRef_ * Un * (magr * dotMr + fwh_.c0_ * (Mr - magM*magM)) *
                        OneByOneMrSq * OneByOneMr / magr / magr
                    )*dS
                );
                
                fpart1 = dotlr * (dS / magr / fwh_.c0_) * OneByOneMrSq;
                fpart2 = (lr - lM) * (dS / magr / magr) * OneByOneMrSq;
                fpart3 = lr * (dS / magr / magr / fwh_.c0_) * OneByOneMrSq * OneByOneMr * 
                        (magr * dotMr + fwh_.c0_ * Mr - fwh_.c0_ * magM * magM);
                    
                fds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
                fds_[iObs][iSurf][iFace].second().append
                (
                        fpart1 + fpart2 + fpart3
                );
                
            }//observer
        } //For Sf
    } // for controlSurfaces_
    
    scalar ct1 = ct+fwh_.obr_.time().deltaT().value()*1.0e-6;//slightly increase time to get inside of time step

    scalar retv = 0.0;
    intDotQdS_.value(iObs) = 0.0;
    intFdS_.value(iObs)    = 0.0;
    //calculate acoustic pressure, zero if source didn't reached observer
    forAll(fwh_.controlSurfaces_, iSurf)
    {
        forAll(qds_[iObs][iSurf], iFace)
        {
            retv = valueAt(qds_, iObs, iSurf, iFace, ct1);
                
            //Code to check bisection
            /*
            scalar retv2 = 0.0;
            {
                const pointTimeData& timeData = qds_[iObs][iSurf][iFace];
                if (timeData.first().size() < 1)
                {
                    retv2 = 0.0;
                }
                if (ct < timeData.first()[0])
                {
                    retv2 = 0.0;
                }
                if (ct > timeData.first()[timeData.first().size()-1])
                {
                    retv2 = 0.0;
                }
                for(label k=1; k<timeData.first().size(); k++)
                {
                    label kl = k-1;
                    if (ct == timeData.first()[kl])
                    {
                        retv2 = timeData.second()[kl];
                        break;
                    }
                    if (ct == timeData.first()[k])
                    {
                        retv2 = timeData.second()[k];
                        break;
                    }
                    if ( (ct > timeData.first()[kl]) && (ct < timeData.first()[k]) )
                    {
                        retv2 = timeData.second()[kl] + 
                        (
                            (timeData.second()[k] - timeData.second()[kl])
                            /
                            (timeData.first()[k] - timeData.first()[kl])
                        ) * (ct - timeData.first()[kl]);
                        break;
                    }
                }
                if (mag(retv - retv2) > 1.0e-8)
                {
                    Info << "Error in retv: " << timeData << endl;
                }
            }
            */
            
            intDotQdS_.value(iObs) += retv;
            retv = valueAt(fds_, iObs, iSurf, iFace, ct1);
            intFdS_.value(iObs) += retv;
        }
    }

    reduce (intDotQdS_.value(iObs), sumOp<scalar>());
    reduce (intFdS_.value(iObs), sumOp<scalar>());

    scalar coeff1 = 1. / 4. / Foam::constant::mathematical::pi;
    
    return (intDotQdS_.value(iObs) + intFdS_.value(iObs))*coeff1;
}

void Foam::functionObjects::Farassat1AFormulation::update()
{
    fwhFormulation::update();
}

void Foam::functionObjects::Farassat1AFormulation::clearExpiredData()
{
    fwhFormulation::clearExpiredData();
}

// ************************************************************************* //

//
//END OF FILE
//


