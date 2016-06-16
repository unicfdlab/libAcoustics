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

#include "IFstream.H"
#include "DynamicList.H"
#include "fft.H"
#include "SubField.H"
#include "mathematicalConstants.H"

Foam::SoundObserver::SoundObserver()
:
    name_(Foam::word::null),
    position_(vector::zero),
    pref_(1.0e-5),
    apressure_(0.0),
    atime_(0.0),
    p_(0),
    time_(0),
    fftFreq_(1024)
{
}

Foam::SoundObserver::SoundObserver(word name, vector pos, scalar pref, label fftFreq)
:
    name_(name),
    position_(pos),
    pref_(pref),
    apressure_(0.0),
    atime_(0.0),
    p_(0),
    time_(0),
    fftFreq_(fftFreq)
{
}


Foam::SoundObserver::SoundObserver(const SoundObserver& so)
:
    name_(so.name_),
    position_(so.position_),
    pref_(so.pref_),
    apressure_(so.apressure_),
    atime_(so.atime_),
    p_(so.p_),
    time_(so.time_),
    fftFreq_(so.fftFreq_)
{
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

const Foam::scalar& Foam::SoundObserver::atime() const
{
    return atime_;
}

const Foam::List<Foam::scalar>& Foam::SoundObserver::p() const
{
    return p_;
}

const Foam::List<Foam::scalar>& Foam::SoundObserver::time() const
{
    return time_;
}

const Foam::scalar& Foam::SoundObserver::pref() const
{
    return pref_;
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

void Foam::SoundObserver::atime(scalar atime)
{
    atime_ = atime;
    time_.append(atime);
}

bool Foam::SoundObserver::checkOrder(scalar size) const
{
  double reminder;
  reminder = log10(size)/log10(2.0);
  Info<<"Checking order rem="<<reminder<<endl;
  if (reminder > floor(reminder)){
    Info<<"False"<<endl;
    return false;
  }
  else{
    return true;
  };
}

Foam::autoPtr<Foam::List<Foam::List<Foam::scalar> > > Foam::SoundObserver::fft(scalar tau) const
{
    List<List<scalar> > fft_res(3);

    forAll (fft_res, i)
    {
	fft_res[i].resize(0);
    }
    
    if ( (p_.size() > 0) && (p_.size() % fftFreq_ == 0) && checkOrder(p_.size()) )
    {
        tmp<scalarField> tPn2
	(
	   mag
	   (
	    fft::reverseTransform
	    (
	     ReComplexField(p_),
	     labelList(1, p_.size())
	    )
	   )
	 );
   
	tmp<scalarField> tPn
	(
	 new scalarField
	 (
	  scalarField::subField(tPn2(), tPn2().size()/2)
	 )
	);
	
	scalarField& Pn = tPn();
	Pn *= 2.0/sqrt(scalar(tPn2().size()));
	Pn[0] /= 2.0;
	
	scalar N = p_.size();	
	scalarField f(N/2);
	scalar deltaf = 1.0/(N*tau);
	forAll(f, i)
	{
	  f[i] = i*deltaf;
	}

	 fft_res[0].resize(tPn().size());
	 fft_res[1].resize(tPn().size());
	 fft_res[2].resize(tPn().size());
	
	 forAll (tPn(), k)
	 {
	     fft_res[0][k] = f[k]; //Frequency, Hz
	     fft_res[1][k] = tPn()[k]; //pressure amplitude, Pa
	     fft_res[2][k] = 20*log10(fft_res[1][k] / pref_); //SLP, dB
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

