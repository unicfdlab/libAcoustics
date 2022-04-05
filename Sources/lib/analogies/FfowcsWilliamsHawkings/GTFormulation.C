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
    Q_(0),
    L_(0),

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

    L_.resize(fwh_.observers_.size());
    Q_.resize(fwh_.observers_.size());

    forAll(L_, iObs)
    {
        L_[iObs].resize(fwh_.controlSurfaces_.size());
        Q_[iObs].resize(fwh_.controlSurfaces_.size());

        forAll(L_[iObs], iSurf)
        {
            L_[iObs][iSurf].setTimeStepVariable(fwh_.obr_.time().isAdjustTimeStep());
            L_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
            Q_[iObs][iSurf].setTimeStepVariable(fwh_.obr_.time().isAdjustTimeStep());
            Q_[iObs][iSurf].resize(fwh_.controlSurfaces_[iSurf].Cf().size());
        }
    }
}

void Foam::functionObjects::GTFormulation::calculateAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    vector r  (vector::zero);
    vector n  (vector::zero);
    scalar dS (0.0);

    //GTFormulation
    scalar R (0.0);
    vector Rh (vector::zero);
    scalar M0 = mag(fwh_.U0_) / fwh_.c0_;
    scalar sqrBetta  = 1 - sqr(M0);
    //-
    tensor Lij (tensor::zero);
    vector Qi(vector::zero);
    //-
    scalar Qn(0.0);
    vector Li (vector::zero);
    scalar LR (0.0);
    scalar LM (0.0);

    vector M  (vector::zero);
    scalar magM (0.0);
    scalar MR (0.0);

    scalar OneByOneMr(0.0);
    scalar OneByOneMrSq(0.0);

    scalar fpart1 (0.0);
    scalar fpart2 (0.0);
    scalar fpart3 (0.0);

    scalar dotLR (0.0);
    scalar dotQn(0.0);

    //GT formulation
    forAll(Sf, iFace)
    {
        //For observe No iObs
        {
            r = robs_[iObs][iSurf][iFace];

            scalar Rz = sqrt( sqr(r[0]) + sqrBetta * (sqr(r[1]) + sqr(r[2])) );
            R = -M0 * r[0] / sqrBetta + Rz / sqrBetta;
            Rh = vector
                (
                    (-M0 * Rz + r[0]) / sqrBetta / R,
                    r[1] / R,
                    r[2] / R
                );

            dS = mag(Sf[iFace]);
            n = ni_[iSurf].value(iFace);
            M = fwh_.vS_[iSurf][iFace] / fwh_.c0_;
            magM = mag(M);

            Qi = (-fwh_.rhoRef_ * fwh_.U0_ + rhoS[iFace] * (uS[iFace] + fwh_.U0_));
            Lij = pS[iFace] * tensor::I + rhoS[iFace] * uS[iFace] * (uS[iFace] + fwh_.U0_);
            Qn = Qi & n;
            Li = Lij & n;

            LR = Li & Rh;
            LM = Li & M; 
            MR = M & Rh;
            OneByOneMr = 1.0 / (1.0 - MR); 
            OneByOneMrSq = OneByOneMr * OneByOneMr;

            Q_[iObs][iSurf].value(iFace) = Qi;
            L_[iObs][iSurf].value(iFace) = Li;

            dotQn = Q_[iObs][iSurf].dot(ct, iFace) & n;
            dotLR = L_[iObs][iSurf].dot(ct, iFace) & Rh;

            qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
            qds_[iObs][iSurf][iFace].second().append
            (
                (
                    (dotQn) * OneByOneMrSq / R 
                    +
                    Qn * fwh_.c0_ * (MR - magM*magM) * OneByOneMrSq * OneByOneMr / R / R 
                ) * dS
            );

            fpart1 = (dotLR) / R / fwh_.c0_ * OneByOneMrSq * dS;
            fpart2 = (LR - LM) / R / R * OneByOneMrSq * dS; 
            fpart3 = LR * (MR - magM * magM) / R / R * OneByOneMrSq * OneByOneMr * dS;

            fds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
            fds_[iObs][iSurf][iFace].second().append
            (
                fpart1 + fpart2 + fpart3
            );
        }//observer
    } //For Sf
}

Foam::scalar Foam::functionObjects::GTFormulation::observerAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    calculateAcousticPressure(Sf,uS,rhoS,pS,iObs,iSurf,ct);

    //slightly increase time to get inside of time step
    scalar ct1 = ct+fwh_.obr_.time().deltaT().value()*1.0e-6;

    scalar retv = 0.0;
    intDotQdS_.value(iObs) = 0.0;
    intFdS_.value(iObs)    = 0.0;
    //calculate acoustic pressure, zero if source didn't reached observer
    forAll(fwh_.controlSurfaces_, iSurf)
    {
        forAll(qds_[iObs][iSurf], iFace)
        {
            retv = valueAt(qds_, iObs, iSurf, iFace, ct1);
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