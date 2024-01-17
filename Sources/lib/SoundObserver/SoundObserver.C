/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SoundObserver.H"
#include "mathematicalConstants.H"

Foam::SoundObserver::SoundObserver()
:
    name_(Foam::word::null),
    position_(vector::zero),
    pref_(1.0e-5),
    apressure_(0.0),
    p_(0),
    fftFreq_(1024)
{
}

Foam::SoundObserver::SoundObserver(word name, vector pos, scalar pref, label fftFreq)
:
    name_(name),
    position_(pos),
    pref_(pref),
    apressure_(0.0),
    p_(0),
    fftFreq_(fftFreq)
{
}


Foam::SoundObserver::SoundObserver(const SoundObserver& so)
:
    name_(so.name_),
    position_(so.position_),
    pref_(so.pref_),
    apressure_(so.apressure_),
    p_(so.p_),
    fftFreq_(so.fftFreq_)
{
}

Foam::SoundObserver::SoundObserver(SoundObserver&& so)
:
    name_(so.name_),
    position_(so.position_),
    pref_(so.pref_),
    apressure_(so.apressure_),
    p_(so.p_),
    fftFreq_(so.fftFreq_)
{
}

const Foam::SoundObserver &Foam::SoundObserver::operator = (const SoundObserver& so)
{
    name_ = so.name_;
    position_ = so.position_;
    pref_ = so.pref_;
    apressure_ = so.apressure_;
    p_ = so.p_;
    fftFreq_ = so.fftFreq_;
    return *this;
}

const Foam::word& Foam::SoundObserver::name() const
{
    return name_;
}

const Foam::vector& Foam::SoundObserver::position() const
{
    return position_;
}

const Foam::scalar& Foam::SoundObserver::apressure() const
{
    return apressure_;
}

void Foam::SoundObserver::name(word name)
{
    name_ = name;
}

void Foam::SoundObserver::position(vector position)
{
    position_ = position;
}

void Foam::SoundObserver::apressure(scalar apressure)
{
    apressure_ = apressure;
    p_.append(apressure);
}

Foam::autoPtr<Foam::List<Foam::List<Foam::scalar> > > Foam::SoundObserver::fft(scalar tau) const
{

    List<List<scalar> > fft_res(3);
    forAll (fft_res, i)
    {
        fft_res[i].resize(0);
    }
    
    if ( (p_.size() > 0) && (p_.size() % fftFreq_ == 0) )
    {
        FoamFftwDriver fftw (p_, tau);
        
        autoPtr<Pair<List<scalar> > > pfft = fftw.simpleScalarForwardTransform();
        
        fft_res[0].resize(pfft().first().size());
        fft_res[1].resize(pfft().first().size());
        fft_res[2].resize(pfft().first().size());
        
        forAll (pfft().first(), k)
        {
            fft_res[0][k] = pfft().first()[k]; //Frequency, Hz
            fft_res[1][k] = pfft().second()[k]; //pressure amplitude, Pa
            if (fft_res[1][k] > SMALL)
            {
                fft_res[2][k] = 20*log10(fft_res[1][k] / pref_); //SPL, dB
            }
            else
            {
                fft_res[2][k] = 0.0;
            }
        }
    }
    
    return autoPtr<List<List<scalar> > >
    (
        new List<List<scalar> >
        (
            fft_res
        )
    );
}

//
//END-OF-FILE
//

