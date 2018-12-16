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

#include "Curle.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Curle, 0);
    
    addToRunTimeSelectionTable(functionObject, Curle, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Curle::Curle
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    AcousticAnalogy
    (
        name,
        runTime,
        dict
    ),
    c_(vector::zero),
    F_(vector::zero, obr_.time().value()),
    dF_(vector::zero, obr_.time().value()),
    FF_(0, vector::zero),
    iter_(0)
{
    this->read(dict);
    F_.resize(1);
}

Foam::functionObjects::Curle::Curle
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    AcousticAnalogy
    (
        name,
        obr,
        dict
    ),
    c_(vector::zero),
    F_(vector::zero, obr_.time().value()),
    dF_(vector::zero, obr_.time().value()),
    FF_(0, vector::zero),
    iter_(0)
{
    this->read(dict);
    F_.resize(1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Curle::~Curle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Curle::read(const dictionary& dict)
{
    if (!AcousticAnalogy::read(dict))
    {
        return false;
    }
    
    calcDistances();

    return true;
}

void Foam::functionObjects::Curle::calcDistances()
{

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    vectorField ci(0);
    scalar ni = 0;
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        
        ci.append(mesh.boundary()[patchi].Cf());
        ni += scalar(ci.size());
    }
    reduce (ni, sumOp<scalar>());
    c_ = gSum(ci) / ni;
}

bool Foam::functionObjects::Curle::execute()
{
    return true;
}

bool Foam::functionObjects::Curle::write()
{
    //use forces library to calculate forces acting on patch
    if (!AcousticAnalogy::write())
    {
        return false;
    }
    
    return true;
}

void Foam::functionObjects::Curle::correct()
{
    iter_ += 1;
    F_.value(0) = forceEff();
    vector F1_= vector::zero;
    
    FF_.setSize(iter_);
    FF_[iter_ - 1] = F_.value(0);
    forAll(FF_, I)
    {
	F1_ += FF_[I];
    }
    vector Fav_ = F1_ / FF_.size();

    dF_.value(0) = F_.value(0) - Fav_; 
    
    vector dotF = F_.dot(obr_.time().value(), 0);
    vector dotdF = dF_.dot(obr_.time().value(), 0);
    
    if (Pstream::master() || !Pstream::parRun())
    {
        scalar coeff1_3d = 1. / 4. / Foam::constant::mathematical::pi / c0_;
        
        scalar coeff1_2d = 1. / 2.828427 / Foam::constant::mathematical::pi / sqrt(c0_);
        scalar coeff2_2d = 1. / 2. / Foam::constant::mathematical::pi;
        scalar coeff3_2d = sqrt(c0_) / 5.656854 / Foam::constant::mathematical::pi;
        
        forAll (observers_, iObs)
        {
            SoundObserver& obs = observers_[iObs];
            vector l = obs.position() - c_;
            scalar r = mag(l);
            scalar oap = 0;
            
            if (dRef_ > 0.0)
            {
        	scalar A1 = l & (dotF) * coeff1_2d / sqrt(r) / r;
        	scalar B1 = l & (Fav_) * coeff2_2d / r / r;
        	scalar C1 = l & (dotdF) * coeff3_2d / sqrt(r) / r / r;
        	oap = A1 + B1 + C1;
        	oap /= dRef_;
    	    }
    	    else
    	    {
    		oap = l & (dotF + c0_ * F_.value(0) / r) * coeff1_3d / r / r;
            }
            obs.apressure(oap);
        }
    }
}



// ************************************************************************* //
