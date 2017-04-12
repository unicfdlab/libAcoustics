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
    F_(vector::zero, obr_.time().value())
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
    F_(vector::zero, obr_.time().value())
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
    F_.value(0) = forceEff();
    vector dotF = F_.dot(obr_.time().value(), 0);
    
    if (Pstream::master() || !Pstream::parRun())
    {
        scalar coeff1 = 1. / 4. / Foam::constant::mathematical::pi / c0_;
        
        forAll (observers_, iObs)
        {
            SoundObserver& obs = observers_[iObs];
            vector l = obs.position() - c_;
            scalar r = mag(l);
            scalar oap = l & (dotF + c0_ * F_.value(0) / r) * coeff1 / r / r;
            if (dRef_ > 0.0)
            {
                oap /= dRef_;
            }
            obs.apressure(oap);
        }
    }
}



// ************************************************************************* //
