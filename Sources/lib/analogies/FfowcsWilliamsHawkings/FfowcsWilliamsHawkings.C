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

#include "FfowcsWilliamsHawkings.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

//sample surface
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"

#include "ListListOps.H"
#include "stringListOps.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "sampledPatch.H"

#include "fwhFormulation.H"
#include "Farassat1AFormulation.H"
#include "GTFormulation.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(FfowcsWilliamsHawkings, 0);
    
    addToRunTimeSelectionTable(functionObject, FfowcsWilliamsHawkings, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::FfowcsWilliamsHawkings::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
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
    formulationType_(word::null),
    fixedResponseDelay_(true),
    responseDelay_(0.0),
    Ufwh_(vector::zero),
    U0_(vector::zero),
    nonUniformSurfaceMotion_(false),
    Cf0_(0),
    vS_(0),
    pInf_(0.0),
    fwhFormulationPtr_(nullptr),
    cleanFreq_(100),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    mergeList_(0)
{
    this->read(dict);
    this->update();
    this->initialize();
}

Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
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
    formulationType_(word::null),
    fixedResponseDelay_(true),
    responseDelay_(0.0),
    Ufwh_(vector::zero),
    U0_(vector::zero),
    nonUniformSurfaceMotion_(false),
    Cf0_(0),
    vS_(0),
    pInf_(0.0),
    fwhFormulationPtr_(nullptr),
    cleanFreq_(100),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    mergeList_(0)
{
    this->read(dict);
    this->update();
    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::FfowcsWilliamsHawkings::~FfowcsWilliamsHawkings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::FfowcsWilliamsHawkings::initialize()
{
    if (nonUniformSurfaceMotion_)
    {
        Cf0_.resize(controlSurfaces_.size());
        forAll(controlSurfaces_, iSurf)
        {
            Cf0_[iSurf].resize(controlSurfaces_[iSurf].Cf().size());
            forAll(Cf0_[iSurf], iF)
            {
                Cf0_[iSurf][iF] = controlSurfaces_[iSurf].Cf()[iF];
            }
        }
    }
    
    vS_.resize(controlSurfaces_.size());
    forAll(controlSurfaces_, iSurf)
    {
        vS_[iSurf].resize(controlSurfaces_[iSurf].Cf().size());
        vS_[iSurf] = vector::zero;
    }
    
    //Allocate pointer to FWH formulation
    if ((formulationType_ == "Farassat1AFormulation") or (formulationType_ == "GTFormulation"))
    {
        if (formulationType_ == "Farassat1AFormulation")
        {
    	    fwhFormulationPtr_.set
    	    (
        	new Farassat1AFormulation(*this)
    	    );
    	}
    	else
    	{
    	    fwhFormulationPtr_.set 
    	    (
	    	new GTFormulation(*this)
    	    );
    	}
    }
    else
    {
        Info << "Wrong formulation type: " << formulationType_ << endl
        << "Please, select one of: " << endl
        << "1) Farassat1AFormulation " << endl
        << "2) GTFormulation " << endl;
    }

}


bool Foam::functionObjects::FfowcsWilliamsHawkings::read(const dictionary& dict)
{
    if (!AcousticAnalogy::read(dict))
    {
        return false;
    }
    
    dict.lookup("formulationType") >> formulationType_;
    
    dict.lookup("Ufwh") >> Ufwh_;
    
    dict.lookup("U0") >> U0_;
    
    dict.lookup("pInf") >> pInf_;
    
    dict.lookup("interpolationScheme") >> interpolationScheme_;
    
    if (dict.found("nonUniformSurfaceMotion"))
    {
        dict.lookup("nonUniformSurfaceMotion") >> nonUniformSurfaceMotion_;
    }
    
    if (dict.found("fixedResponseDelay"))
    {
        dict.lookup("fixedResponseDelay") >> fixedResponseDelay_;
    }
    
    if (fixedResponseDelay_ && dict.found("responseDelay"))
    {
        dict.lookup("responseDelay") >> responseDelay_;
    }
    else
    {
        fixedResponseDelay_ = false;
    }
    
    dict.lookup("cleanFreq") >> cleanFreq_;
    
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh)
    );
    
    controlSurfaces_.transfer(newList);
    if (Pstream::parRun())
    {
        mergeList_.setSize(controlSurfaces_.size());
    }

    // Ensure all surfaces and merge information are expired
    expire();
    
    if (controlSurfaces_.size())
    {
        Log << "Function object "<< name()<<":" << nl;
        Log << "    Reading FwocsWilliams-Hawkings analogy control surface description:" << nl;
        forAll(controlSurfaces_, surfI)
        {
            Log << "        " <<  controlSurfaces_.operator[](surfI).name() << nl;        
        }
        Log << endl;
    }
    
    if (Pstream::master() && debug)
    {
        Pout<< "FWH control surfaces additional info:" << nl << "(" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Pout<< "  " << controlSurfaces_.operator[](surfI) << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::execute()
{
    return true;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::write()
{
    //use forces library to calculate forces acting on patch
    if (!AcousticAnalogy::write())
    {
        return false;
    }
    
    //store old faces
    forAll(controlSurfaces_, iSurf)
    {
        if (nonUniformSurfaceMotion_)
        {
            forAll(Cf0_[iSurf], iF)
            {
                Cf0_[iSurf][iF] = controlSurfaces_[iSurf].Cf()[iF];
            }
        }
    }
    
    return true;
}

void Foam::functionObjects::FfowcsWilliamsHawkings::correct()
{
    if (nonUniformSurfaceMotion_)
    {
        this->expire();
    }
    this->update();

    scalar dt   = obr_.time().deltaT().value();
    
    //update surface velocities (if needed)
    forAll(controlSurfaces_, iSurf)
    {
        if (nonUniformSurfaceMotion_)
        {
            vS_[iSurf] = (controlSurfaces_[iSurf].Cf() - Cf0_[iSurf])/dt;
        }
        else
        {
            vS_[iSurf] = Ufwh_;
        }
    }
    
    //update formulation-specific data
    fwhFormulationPtr_->update();
    
    if (fwhFormulationPtr_.valid())
    {
        forAll(observers_, iObs)
        {
            scalar oap = fwhFormulationPtr_->observerAcousticPressure(iObs);
            if (Pstream::master())
            {
                SoundObserver& obs = observers_[iObs];
            
                if (dRef_ > 0.0)
                {
                    oap /= dRef_;
                }
                obs.apressure(oap); //appends new calculated acoustic pressure
                Log<<"OAP = "<< oap <<nl;
            }
        }
    }
    else
    {
        if (Pstream::master())
        {
            forAll(observers_, iObs)
            {
                SoundObserver& obs = observers_[iObs];
                obs.apressure(0.0);
            }
        }
    }

    //Remove old data if needed
    if (fwhFormulationPtr_.valid())
    {
        fwhFormulationPtr_->clearExpiredData();
    }
}

//bool Foam::functionObjects::FfowcsWilliamsHawkings::signalReachedObserver(const surfaceTimeData& data, label iObs)
//{
//    scalar ct = obr_.time().value();
//    label lasttau = 0;
//    label signalReached = 1;
//    Info << "in signalReachedObserver" << endl;
//    forAll(data[iObs], iSurf)
//    {
//        forAll(data[iObs][iSurf],iFace)
//        {
//            lasttau = data[iObs][iSurf][iFace].first().size() - 1;
//            if (lasttau < 0)
//            {
//                signalReached = 0;
//                break;
//            }
//            if (ct < data[iObs][iSurf][iFace].first()[lasttau])
//            {
//                signalReached = 0;
//                Info << "ct = " << ct << " lasttm = " << data[iObs][iSurf][iFace].first()[lasttau] << endl;
//                break;
//            }
//        }
//    }
//    
//    reduce (signalReached, minOp<label>());
//    
//    if (signalReached)
//    {
//        Log << "Obs " << iObs << " is listening at time " << ct << endl;
//    }
//    
//    return signalReached;
//}

Foam::tmp<Foam::scalarField> Foam::functionObjects::FfowcsWilliamsHawkings::surfaceDensity(const sampledSurface& surface) const
{
    tmp<Field<scalar> > rhoSampled
    (
        sampleOrInterpolate<scalar>(this->rho()(), surface)
    );

    return rhoSampled;
}

Foam::tmp<Foam::vectorField> Foam::functionObjects::FfowcsWilliamsHawkings::surfaceVelocity(const sampledSurface& surface) const
{
    const volVectorField& U = obr_.lookupObject<volVectorField>("U");
 
    tmp<Field<vector> > USampled;
    
    USampled = sampleOrInterpolate<vector>(U , surface);

    return USampled;
}

Foam::tmp<Foam::vectorField> Foam::functionObjects::FfowcsWilliamsHawkings::surfaceStressDivergence(const sampledSurface& surface) const
{
    tmp<Field<vector>> divTSampled;
    
    volScalarField p (obr_.lookupObject<volScalarField>("p"));
    volVectorField U (obr_.lookupObject<volVectorField>("U"));
    
    volVectorField gradp = fvc::grad(p);
    
    //Interpolate FWH surface velocity to CFD mesh
    
    volVectorField Ufwh("Ufwh", U*0.0);
    forAll(controlSurfaces_, iSurf)
    {
        if (isA<sampledPatch>(controlSurfaces_[iSurf]))
        {
            //Write data directly to Ufwh
            const sampledPatch& pSurf = refCast<const sampledPatch>(controlSurfaces_[iSurf]);
            
            sampledPatchAccess * spa = new sampledPatchAccess(pSurf);

            labelList pIds  = spa->aPatchIDs();
            labelList pAddr = spa->aPatchFaceLabels();
            labelList pStart= spa->aPatchStart();
            const vectorField& vp = vS_[iSurf];
            
            label patchId = -1;
            label startFace = -1;
            label lastFace = -2;
            label iF = -1;
            forAll(pIds, iPatch)
            {
                patchId = pIds[iPatch];
                startFace = pStart[iPatch];
                if (iPatch >= pIds.size() - 1)
                {
                    lastFace = pAddr.size() - 1;
                }
                else
                {
                    lastFace = pStart[iPatch+1] -1;
                }
                fvPatchField<vector>& pUfwh = Ufwh.boundaryFieldRef()[patchId];
                for (label k=startFace; k<=lastFace; k++)
                {
                    iF = pAddr[k];
                    pUfwh[iF] = vp[k];
                }
            }
        }
        else
        {
            //use interpolation procedure
        }
    }
    
    volTensorField momConv = this->rho()()*U*(U-Ufwh);
    volVectorField divMomConv = fvc::div(momConv);
    
    divTSampled = sampleOrInterpolate<vector>(gradp,surface);
    if (p.dimensions() != dimPressure)
    {
        divTSampled.ref() *= rhoRef_;
    }
    
    divTSampled.ref() += sampleOrInterpolate<vector>(divMomConv,surface);

    return divTSampled;
}

Foam::tmp<Foam::scalarField> Foam::functionObjects::FfowcsWilliamsHawkings::surfacePressure(const sampledSurface& surface) const
{
    tmp<Field<scalar> > pSampled;
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    
    pSampled = sampleOrInterpolate<scalar>(p , surface);
    
    if (p.dimensions() != dimPressure)
    {
        pSampled.ref() *= rhoRef_;
    }
    
    //Info << pSampled() << endl;
    
    return pSampled;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::expire()
{
    bool justExpired = false;

    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).expire())
        {
            justExpired = true;
        }

        //Clear merge information
        if (Pstream::parRun())
        {
          mergeList_[surfI].clear();
        }
    }

    // true if any surfaces just expired
    return justExpired;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::needsUpdate() const
{
    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // Serial: quick and easy, no merging required
    // Just like sampledSurfaces
    if (!Pstream::parRun())
    {
        forAll(controlSurfaces_, surfI)
        {
            if (controlSurfaces_.operator[](surfI).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
   
    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh.bounds().mag();

    if (Pstream::master() && debug)
    {
      Pout<< nl << "Merging all points within "
          << mergeDim << " metre" << endl;
    }

    forAll(controlSurfaces_, surfI)
    {
        sampledSurface& s = controlSurfaces_.operator[](surfI);

        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }

        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(s.faces(), s.faces().size()),
                s.points()
            ),
            mergeList_[surfI].points,
            mergeList_[surfI].faces,
            mergeList_[surfI].pointsMap
        );
    }
    


    return updated;
}


// ************************************************************************* //
