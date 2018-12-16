#include "FfowcsWilliamsHawkings.H"
#include "fwhFormulation.H"

Foam::functionObjects::fwhFormulation::fwhFormulation(const FfowcsWilliamsHawkings& fwh)
:
    fwh_(fwh),
    fwhProbeI_(0),
    qds_(0),
    fds_(0),
    tobs_(0),
    robs_(0),
    magrobs_(0),
    ni_(0),
    nl_(0),
    rMax_(0),
    tauMax_(0),
    tauMin_(0)
{
    this->initialize();
}

void Foam::functionObjects::fwhFormulation::initialize()
{
    //allocate qds_, fds_ and vds_
    qds_.resize(fwh_.observers_.size());
    fds_.resize(fwh_.observers_.size());
    forAll(fwh_.observers_, iObs)
    {
        qds_[iObs].resize(fwh_.controlSurfaces_.size());
        fds_[iObs].resize(fwh_.controlSurfaces_.size());
        forAll(fwh_.controlSurfaces_, iSurf)
        {
            qds_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            fds_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
        }
    }

    //allocate tobs
    tobs_.resize(fwh_.observers_.size());
    forAll(fwh_.observers_, iObs)
    {
        tobs_[iObs].resize(fwh_.controlSurfaces_.size());
        forAll(fwh_.controlSurfaces_, iSurf)
        {
            tobs_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
        }
    }
    
    tauMax_.resize(fwh_.observers_.size(), 0.0);
    rMax_.resize(fwh_.observers_.size(), 0.0);
    
    //allocate robs
    robs_.resize(fwh_.observers_.size());
    magrobs_.resize(fwh_.observers_.size());
    
    forAll(fwh_.observers_, iObs)
    {
        rMax_[iObs] = 0.0;
        const SoundObserver& obs = fwh_.observers_[iObs];
        robs_[iObs].resize(fwh_.controlSurfaces_.size());
        magrobs_[iObs].resize(fwh_.controlSurfaces_.size());
        forAll(fwh_.controlSurfaces_, iSurf)
        {
            robs_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            magrobs_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            const vectorField& Cf = fwh_.controlSurfaces_[iSurf].Cf();
            forAll(Cf, i)
            {
                robs_[iObs][iSurf][i] = obs.position() - Cf[i];
                vector r = robs_[iObs][iSurf][i];
                scalar R_ = sqrt
            	    (
            		sqr(r[0]) 
            		+ 
            		( 1 - sqr(mag(fwh_.U0_)/fwh_.c0_))
            		*
            		( sqr(r[1]) + sqr(r[2]) )
            	    );
            	
//                magrobs_[iObs][iSurf][i] = mag(robs_[iObs][iSurf][i]);
		magrobs_[iObs][iSurf][i] = 
		    (
			-(mag(fwh_.U0_)/fwh_.c0_) * r[0] + R_
		    ) /	(1 - sqr(mag(fwh_.U0_)/fwh_.c0_));
		     

                if (magrobs_[iObs][iSurf][i] > rMax_[iObs])
                {
                    rMax_[iObs] = magrobs_[iObs][iSurf][i];
                }
            }
        }
        reduce(rMax_[iObs], maxOp<scalar>());
        tauMax_[iObs] = rMax_[iObs] / fwh_.c0_;
        reduce(tauMax_[iObs], maxOp<scalar>());
    }
    
    if (fwh_.fixedResponseDelay_)
    {
        tauMin_.resize(tauMax_.size());
        rMin_.resize(rMax_.size());
        tauMin_ = max(tauMax_);
        rMin_ = max(rMax_);
        
        forAll(rMin_, iObs)
        {
            forAll(fwh_.controlSurfaces_, iSurf)
            {
                const vectorField& Cf = fwh_.controlSurfaces_[iSurf].Cf();
                forAll(Cf, i)
                {
                    if (magrobs_[iObs][iSurf][i] < rMin_[iObs])
                    {
                        rMin_[iObs] = magrobs_[iObs][iSurf][i];
                    }
                }
            }
            tauMin_[iObs] = rMin_[iObs] / fwh_.c0_;
            reduce(rMin_[iObs], minOp<scalar>());
            reduce(tauMin_[iObs], minOp<scalar>());
        }
    }

    //calculate normals
    List<vector> Cs(fwh_.controlSurfaces_.size());
    ni_.resize(fwh_.controlSurfaces_.size());
    nl_.resize(fwh_.controlSurfaces_.size());
    forAll(ni_, iSurf)
    {
        const vectorField& Sf = fwh_.controlSurfaces_[iSurf].Sf();
        const vectorField& Cf = fwh_.controlSurfaces_[iSurf].Cf();
        ni_[iSurf].resize(Cf.size());
        nl_[iSurf].resize(Cf.size());
        Cs[iSurf] = gSum(Cf);
        scalar surfSize = scalar(Cf.size());
        reduce (surfSize, sumOp<scalar>());
        Cs[iSurf] /= surfSize;
        
        scalar magSf = 0.0;
        forAll(ni_[iSurf], iFace)
        {
            magSf = mag(Sf[iFace]);
            ni_[iSurf].value(iFace) = Sf[iFace]/magSf;
            if ( ((Cf[iFace] - Cs[iSurf]) & ni_[iSurf].value(iFace)) > 0 )
            {
                nl_[iSurf][iFace] = 1.0;
            }
            else
            {
                nl_[iSurf][iFace] = -1.0;
            }
            ni_[iSurf].value(iFace) *= nl_[iSurf][iFace];
        }
    }
}

Foam::functionObjects::fwhFormulation::~fwhFormulation()
{
}

Foam::scalar Foam::functionObjects::fwhFormulation::observerAcousticPressure(label iObs)
{
    return 0.0;
}

void Foam::functionObjects::fwhFormulation::clearExpiredData()
{
    scalar ct   = fwh_.obr_.time().value();// - fwh_.obr_.time().deltaT().value()*1.0e-6;
    reduce(ct, minOp<scalar>());
    
    fwhProbeI_++;
    
    if ( mag(fwhProbeI_ % fwh_.cleanFreq_) > VSMALL  )
    {
        
    }
    else
    {
        fwhProbeI_ = 0;
        scalar expiredTime = 0.0;
        label  expiredIndex= -1;
        label  newsize     = 0;
        forAll(qds_, iObs)
        {
            forAll(qds_[iObs], iSurf)
            {
                forAll(qds_[iObs][iSurf], iFace)
                {
                    if (tauMin_.size())
                    {
                        expiredTime = ct - (tauMax_[iObs] - tauMin_[iObs] + fwh_.responseDelay_);
                    }
                    else
                    {
                        expiredTime = ct - tauMax_[iObs];
                    }
                    const pointTimeData& qdsOldPointData = qds_[iObs][iSurf][iFace];
                    expiredIndex= findExpiredIndex(qdsOldPointData, expiredTime);
                    
                    // -1 - if nothing found, from 0 to (size-1) for indices to remove
                    if (expiredIndex > -1)
                    {
                        newsize = qdsOldPointData.first().size() - (expiredIndex + 1);
                        
                        //clean qds
                        pointTimeData newPointData;
                        newPointData.first().resize(newsize);
                        newPointData.second().resize(newsize);
                        for(label iTime=expiredIndex+1; iTime<qdsOldPointData.first().size(); iTime++)
                        {
                            newPointData.first() [iTime-(expiredIndex+1)] = qdsOldPointData.first()[iTime];
                            newPointData.second()[iTime-(expiredIndex+1)] = qdsOldPointData.second()[iTime];
                        }
                        qds_[iObs][iSurf][iFace].first().operator=(newPointData.first());
                        qds_[iObs][iSurf][iFace].second().operator=(newPointData.second());
                        
                        //clean fds
                        const pointTimeData& fdsOldPointData = fds_[iObs][iSurf][iFace];
                        for(label iTime=expiredIndex+1; iTime<fdsOldPointData.first().size(); iTime++)
                        {
                            newPointData.first()[iTime-(expiredIndex+1)] = fdsOldPointData.first()[iTime];
                            newPointData.second()[iTime-(expiredIndex+1)] = fdsOldPointData.second()[iTime];
                        }
                        
                        fds_[iObs][iSurf][iFace].first().operator=(newPointData.first());
                        fds_[iObs][iSurf][iFace].second().operator=(newPointData.second());
                    }
                }
            }
        }
    }
}

void Foam::functionObjects::fwhFormulation::update()
{
    scalar ct   = fwh_.obr_.time().value();
    
    if (mag(fwh_.Ufwh_) > SMALL || fwh_.nonUniformSurfaceMotion_)
    {
        forAll(fwh_.observers_, iObs)
        {
            rMax_[iObs] = 0.0;
            tauMax_[iObs] = 0.0;
            if (rMin_.size())
            {
                rMin_[iObs] = GREAT;
            }
            forAll(fwh_.controlSurfaces_, iSurf)
            {
                const vectorField& Cf = fwh_.controlSurfaces_[iSurf].Cf();
                forAll(Cf, i)
                {
                    robs_[iObs][iSurf][i] = fwh_.observers_[iObs].position() - Cf[i];
                    magrobs_[iObs][iSurf][i] = mag(robs_[iObs][iSurf][i]);
                    if (magrobs_[iObs][iSurf][i] > rMax_[iObs])
                    {
                        rMax_[iObs] = magrobs_[iObs][iSurf][i];
                    }
                    
                    if (rMin_.size() && (magrobs_[iObs][iSurf][i] < rMin_[iObs]))
                    {
                        rMin_[iObs] = magrobs_[iObs][iSurf][i];
                    }
                }
            }
            tauMax_[iObs] = rMax_[iObs] / fwh_.c0_;
            reduce(tauMax_[iObs], maxOp<scalar>());
            reduce(rMax_[iObs], maxOp<scalar>());
            tauMin_ = max(tauMax_);
            if (tauMin_.size())
            {
                tauMin_[iObs] = rMin_[iObs] / fwh_.c0_;
                reduce(tauMin_[iObs], minOp<scalar>());
                reduce(rMin_[iObs], minOp<scalar>());
            }
        }
        
        forAll(ni_, iSurf)
        {
            const vectorField& Sf = fwh_.controlSurfaces_[iSurf].Sf();
            scalar magSf = 0.0;
            forAll(ni_[iSurf], iFace)
            {
                magSf = mag(Sf[iFace]);
                ni_[iSurf].value(iFace) = nl_[iSurf][iFace]*Sf[iFace]/magSf;
            }
        }
    }
    
    forAll(fwh_.observers_, iObs)
    {
    //    Info << "iObs = " << iObs << " tauMax = " << tauMax_[iObs] << endl;
    //    Info << "iObs = " << iObs << " tauMin = " << tauMin_[iObs] << endl;
    //    Info << "iObs = " << iObs << " Delta  = " << (tauMax_[iObs] - tauMin_[iObs]) << endl;
        forAll(fwh_.controlSurfaces_, iSurf)
        {
            const sampledSurface& surf = fwh_.controlSurfaces_[iSurf];
            const vectorField& Cf = surf.Cf();
            forAll(Cf,i)
            {
                tobs_[iObs][iSurf][i] = 
                        ct + magrobs_[iObs][iSurf][i]/fwh_.c0_;
                if (tauMin_.size())
                {
                    tobs_[iObs][iSurf][i] -= (tauMin_[iObs]-fwh_.responseDelay_);
                }
            }
        }
    }
}


Foam::label Foam::functionObjects::fwhFormulation::findExpiredIndex
(
    const pointTimeData& timeData,
    scalar expiredTime
)
{
    label ui = timeData.first().size();
    label li = 0;
    label mi = (ui + li) / 2;
    if (ui < 1)
    {
        return -1;
    }
    if (mi < 1)
    {
        if (timeData.first()[mi] < expiredTime)
        {
            return 0;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        if (timeData.first()[li] >= expiredTime)
        {
            return -1;
        }
        do
        {
            if (timeData.first()[mi] < expiredTime)
            {
                li = mi;
                mi = (ui + li) / 2;
            }
            else
            {
                ui = mi;
                mi = (ui + li) / 2;
            }
            if (expiredTime == timeData.first()[mi])
            {
                return (mi - 1);
            }
        }
        while
        (
            ((ui - li) != 1)
        );
    }
    
    if (timeData.first()[ui] < expiredTime)
    {
        return ui;
    }
    
    return li;
}

Foam::scalar Foam::functionObjects::fwhFormulation::valueAt
(
    const surfaceTimeData& data,
    label iObs,
    label iSurf,
    label iFace,
    scalar tau
)
{
    const pointTimeData& timeData = data[iObs][iSurf][iFace];
    
    label ui = timeData.first().size();
    label li = 0;
    label mi = (ui + li) / 2;
    //exit if size is zero
    if (ui < 1)
    {
        return 0.0;
    }
    //exit if tau out of bounds
    if (
        (tau < timeData.first()[li])
        ||
        (tau > timeData.first()[ui-1])
       )
    {
        return 0.0;
    }
    //exit if tau is at bounds
    if (tau == timeData.first()[li])
    {
        return timeData.second()[li];
    }
    if (tau == timeData.first()[ui-1])
    {
        return timeData.second()[ui-1];
    }
    //check that size is sufficient
    if (mi >= 1)
    {
        do
        {
            if (tau < timeData.first()[mi])
            {
                ui = mi;
                mi = (ui + li) / 2;
            }
            else if (tau > timeData.first()[mi])
            {
                li = mi;
                mi = (ui + li) / 2;
            }
            else //(tau == timeData.first()[mi])
            {
                return timeData.second()[mi];
            }
        }
        while
        (
            ((ui - li) != 1)
        );
    }
    //return value between li and ui
    return timeData.second()[li] + 
    (
        (timeData.second()[ui] - timeData.second()[li])
        / 
        (timeData.first()[ui] - timeData.first()[li])
    ) * (tau - timeData.first()[li]);
}

//
//END-OF-FILE
//



