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

Class
    Foam::SoundObserver

Description
    
    Stores observed pressure fluctuations


SourceFiles
    SoundObserver.C

\*---------------------------------------------------------------------------*/

#ifndef SoundObserver_H
#define SoundObserver_H

#include "OFstream.H"
#include "vector.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SoundObserver Declaration
\*---------------------------------------------------------------------------*/

class SoundObserver
{
    //-
    word name_;

    //-
    vector position_;
    
    //-
    scalar pref_;
    
    //-
    scalar apressure_;

    //-
    scalar atime_;

    //-
    List<scalar> p_;
        
    //-
    List<scalar> time_; 
   
    //-
    label fftFreq_;

    bool checkPressureSize() const;

public:

    /// Constructors
    
    //-
    SoundObserver();
    
    //-
    SoundObserver(word name, vector pos, scalar pref, label fftFreq);
    
    //-
    SoundObserver(const SoundObserver&);

    /// Access
    
    //-
    const word& name() const;
    
    //-
    const vector& position() const;
    
    //-
    const scalar& apressure() const;

//noiseFFT stuff
    const scalar& atime() const;
//-
    const scalar& pref() const; 
    
    //-
    void name(word name);
    
    //-
    void position(vector position);
    
    //-
    void apressure(scalar apressure);

//noiseFFT stuff
    void atime(scalar atime);
//
    const List<scalar>& p() const;
    
//-
    const List<scalar>& time() const;
    
    //-
    autoPtr<List<List<scalar> > > fft (scalar tau) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
