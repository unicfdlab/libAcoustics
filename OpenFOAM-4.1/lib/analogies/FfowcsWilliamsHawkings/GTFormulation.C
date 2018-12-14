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

#include "GTFormulation.H"
#include "FfowcsWilliamsHawkings.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::GTFormulation::GTFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    Qn_(0),
    Lr_(0),
    Mr_(0),

    intDotQdS_(0.0, fwh_.obr_.time().value()),
    intFdS_(0.0, fwh_.obr_.time().value())
{
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::GTFormulation::~GTFormulation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::functionObjects::GTFormulation::initialize()
{
    intFdS_.resize(fwh_.observers_.size());
    intDotQdS_.resize(fwh_.observers_.size());
    
    Lr_.resize(fwh_.observers_.size());
    Mr_.resize(fwh_.observers_.size());
    Qn_.resize(fwh_.observers_.size());

    forAll(Lr_, iObs)
    {
        Lr_[iObs].resize(fwh_.controlSurfaces_.size());
        Mr_[iObs].resize(fwh_.controlSurfaces_.size());
        Qn_[iObs].resize(fwh_.controlSurfaces_.size());

        forAll(Lr_[iObs], iSurf)
        {
            Lr_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            Mr_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            Qn_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
        }
    }
}

Foam::scalar Foam::functionObjects::GTFormulation::observerAcousticPressure(label iObs)
{
    scalar ct = fwh_.obr_.time().value();
    
    //GTFormulation
    scalar R (0.0);
    vector Rh (vector::zero);
    scalar M0 = mag(fwh_.U0_) / fwh_.c0_;
    scalar sqrBetta  = 1 - sqr(M0);
    //-
    vector L (vector::zero);
    scalar Lr (0.0);
    scalar LM (0.0);
    scalar dotLr (0.0);
    vector r  (vector::zero);
    vector n  (vector::zero);
    scalar dS (0.0);
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
    scalar Qn(0.0);
    scalar dotQn(0.0);
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
        
        //GT formulation
        forAll(Sf, iFace)
        {
        
            //For observe No iObs
            {
                r = robs_[iObs][iSurf][iFace];
                
                scalar R_ = sqrt( sqr(r[0]) + sqrBetta * (sqr(r[1]) + sqr(r[2])) );
                R = -M0 * r[0] / sqrBetta + R_ / sqrBetta;
            	Rh = vector
            	    (
            		(-M0 * R_ + r[0]) / sqrBetta / R,
            		r[1] / R,
            		r[2] / R
            	    );
	    
                dS = mag(Sf[iFace]);
                n = ni_[iSurf].value(iFace);
                M = fwh_.vS_[iSurf][iFace] / fwh_.c0_;
                Mr = M & Rh;
                magM = mag(M);
                

                U = (-fwh_.rhoRef_ * fwh_.U0_ + rhoS[iFace] * (uS[iFace] + fwh_.U0_));
		Pf = pS[iFace] * tensor::I + rhoS[iFace] * uS[iFace] * (uS[iFace] + fwh_.U0_);

                Qn = U & n;
                L = Pf & n;
                Lr = L & Rh;
                LM = L & M; 
	        Mr = M & Rh;
                	
            	OneByOneMr = 1.0 / (1.0 - Mr); 
                OneByOneMrSq = OneByOneMr * OneByOneMr;
                
                Qn_[iObs][iSurf].value(iFace) = Qn;
                Lr_[iObs][iSurf].value(iFace) = Lr;
                Mr_[iObs][iSurf].value(iFace) = Mr; 
            
                dotLr = Lr_[iObs][iSurf].dot(ct, iFace);
                dotMr = Mr_[iObs][iSurf].dot(ct, iFace);
                dotQn = Qn_[iObs][iSurf].dot(ct, iFace);
                
                qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
                qds_[iObs][iSurf][iFace].second().append
                (
                    (
                        (dotQn) * OneByOneMrSq / R 
                        +
                        Qn * fwh_.c0_ * (Mr - magM*magM) * OneByOneMrSq * OneByOneMr / R / R 
                    ) * dS
                );
                
                fpart1 = dotLr * (dS / R / fwh_.c0_) * OneByOneMrSq;
                fpart2 = (Lr - LM) * (dS / R / R) * OneByOneMrSq; 
                fpart3 = Lr * (dS / R / R) * OneByOneMrSq * OneByOneMr * (Mr - magM * magM);
                    
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

void Foam::functionObjects::GTFormulation::update()
{
    fwhFormulation::update();
}

void Foam::functionObjects::GTFormulation::clearExpiredData()
{
    fwhFormulation::clearExpiredData();
}

// ************************************************************************* //

//
//END OF FILE
//


