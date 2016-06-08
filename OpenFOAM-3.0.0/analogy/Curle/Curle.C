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

#include "Curle.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
//#include "wordReList.H"
//#include "PtrList.H"
//#include "PtrListIO.C"

#include "RASModel.H"
#include "LESModel.H"

#include "basicThermo.H"

//sampledSurfaces stuff
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"

#include "ListListOps.H"
#include "stringListOps.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

    defineTypeNameAndDebug(Curle, 0);

}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::Curle::normalStress(const sampledSurface& surface) const
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
 
    tmp<Field<scalar> > pSampled;
    
    pSampled = sampleOrInterpolate<scalar>(p, surface);

    if (p.dimensions() == dimPressure)
    {
	//return tmp<scalarField>
	//(
	//    new scalarField(pPatch)
	//);
    }
    else
    {
	if (rhoRef_ < 0) //density in volScalarField
	{
	  volScalarField pRho = obr_.lookupObject<volScalarField>(rhoName_);

	  tmp<Field<scalar> > rhoSampled;
	  rhoSampled = sampleOrInterpolate<scalar>(pRho, surface);

	  pSampled() *= rhoSampled();
	}
	else //density is constant
	{
	  pSampled() *= rhoRef_;
	}
    }

    return pSampled;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Curle::Curle
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    probeFreq_(1),
    log_(false),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    pName_(word::null),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    rhoName_(word::null),
    rhoRef_(1.0),
    c_(vector::zero),
    CurleFilePtr_(NULL),
    FOldPtr_(NULL),
    FOldOldPtr_(NULL),
    probeI_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "Foam::Curle::Curle"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Curle::~Curle()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Curle::read(const dictionary& dict)
{
    if (!active_)
    {
	return;
    }

    log_ = dict.lookupOrDefault<Switch>("log", false);
    
    if (!log_)
    {
	Info << "Direct logging to stdio disabled" << endl
	    << " to enable, please insert string:" << endl
	    << "log\t\t true;" << endl
	    << "in dictionary" << endl;
    }
    
    dict.lookup("probeFrequency") >> probeFreq_;

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh_)
    );

    controlSurfaces_.transfer(newList);
    //Turn on if developing parallel
    // if (Pstream::parRun())
    // {
    //     mergeList_.setSize(size());
    // }


    // Ensure all surfaces and merge information are expired
    expire();

    if (controlSurfaces_.size())
    {
        Info<< "Function object "<< name_<<":" << nl;
        Info<< "    Reading Curle analogy control surface description:" << nl;
        forAll(controlSurfaces_, surfI)
        {
	    Info<< "        " <<  controlSurfaces_.operator[](surfI).name() << nl;        
	}
        Info<< endl;
    }

    // if (Pstream::master() && debug)
    // {
    //     Pout<< "sample fields:" << fieldSelection_ << nl
    //         << "sample surfaces:" << nl << "(" << nl;

    //     forAll(*this, surfI)
    //     {
    //         Pout<< "  " << operator[](surfI) << endl;
    //     }
    //     Pout<< ")" << endl;
    // }
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("c0") >> c0_;
    
    dict.lookup("dRef") >> dRef_;

    dict.lookup("pName") >> pName_;
    
    dict.lookup("rhoName") >> rhoName_;
    
    dict.lookup("rhoRef") >> rhoRef_;

    //read observers
    {
	const dictionary& obsDict = dict.subDict("observers");
	wordList obsNames = obsDict.toc();
	forAll (obsNames, obsI)
	{
	    word oname = obsNames[obsI];
	    vector opos (vector::zero);
	    obsDict.subDict(oname).lookup("position") >> opos;
	    scalar pref = 1.0e-5;
	    obsDict.subDict(oname).lookup("pRef") >> pref;
	    label fftFreq = 1024;
	    obsDict.subDict(oname).lookup("fftFreq") >> fftFreq;
	    
	    observers_.append
	    (
		SoundObserver
		(
		    oname,
		    opos,
		    pref,
		    fftFreq
		)
	    );
	}
    }
    
    calcDistances();   
}

void Foam::Curle::correct()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    //sign '-' needed to calculate force, which exerts fluid by solid
    vector F	(0.0, 0.0, 0.0);
    vector dFdT (0.0, 0.0, 0.0);
    scalar deltaT = mesh.time().deltaT().value();

    //calling a function to update all sampledSurfaces
    //without it everything related will be empty
    update();
    //working with sampled surfaces
    forAll(controlSurfaces_, surfI)
    {
      sampledSurface& s = controlSurfaces_.operator[](surfI);
      scalarField sampledPressure = normalStress(s);
      F -= gSum (sampledPressure * s.Sf());
      Info<<s.name()<<" , sampled Force = "<<F<<nl;
    }

    if (Pstream::master() || !Pstream::parRun())
    {
	//calculate dFdT and store old values
	
	if (FOldPtr_.empty())
	{
	    FOldPtr_.set
	    (
		new vector(F)
	    );
	}
	else
	{
	    if (FOldOldPtr_.empty())
	    {
		//first order scheme
		dFdT = (F - FOldPtr_()) / deltaT;
		
		FOldOldPtr_.set
		(
		    FOldPtr_.ptr()
		);
		
		FOldPtr_.reset
		(
		    new vector(F)
		);
	    }
	    else
	    {
		//second order scheme (BDF)
		dFdT = (3.0*F - 4.0*FOldPtr_() + FOldOldPtr_()) / 2.0 / deltaT;
		
		FOldOldPtr_.reset
		(
		    FOldPtr_.ptr()
		);
		
		FOldPtr_.reset
		(
		    new vector(F)
		);
	    }
	}
	
	scalar coeff1 = 1. / 4. / Foam::constant::mathematical::pi / c0_;
	
	forAll (observers_, iObs)
	{
	    SoundObserver& obs = observers_[iObs];
	    //Vector from observer to center
	    vector l = obs.position() - c_;
	    //Calculate distance
	    scalar r = mag(l);
	    //Calculate ObservedAcousticPressure
	    scalar oap = l & (dFdT + c0_ * F / r) * coeff1 / r / r;
	    if (dRef_ > 0.0)
	    {
		oap /= dRef_;
	    }
	    obs.apressure(oap); //appends new calculated acoustic pressure
	    
	    //noiseFFT addition
	    obs.atime(mesh.time().value());
	}
	
    }
}

void Foam::Curle::makeFile()
{
    fileName CurleDir;

    if (Pstream::master() && Pstream::parRun())
    {
	CurleDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
	mkDir(CurleDir);
    }
    else if (!Pstream::parRun())
    {
	CurleDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
	mkDir(CurleDir);
    }
    else
    {
    }
    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
	// Create the Curle file if not already created
	if (CurleFilePtr_.empty())
	{
	    // Open new file at start up
	    CurleFilePtr_.reset
	    (
		new OFstream
		(
		    CurleDir + "/" + (name_ + "-time.dat")
		)
	    );
	    
	    writeFileHeader();
	}
    }
}


void Foam::Curle::writeFileHeader()
{
    if (CurleFilePtr_.valid())
    {
        CurleFilePtr_()
            << "Time" << " ";
	
        forAll(observers_, iObserver)
        {
	    CurleFilePtr_() << observers_[iObserver].name() << "_pFluct ";
        }

        CurleFilePtr_()<< endl;
    }
}

void Foam::Curle::calcDistances()
{
    if (!active_)
    {
	return;
    }

    update();

    vectorField ci;
    scalar ni;

    forAll(controlSurfaces_, surfI)
    {
      sampledSurface& s = controlSurfaces_.operator[](surfI);
      ci = s.Cf();
      ni = scalar(ci.size());
    }

    reduce (ni, sumOp<scalar>());
    c_ = gSum(ci) / ni;
}

void Foam::Curle::writeFft()
{
    fileName CurleDir;

    if (Pstream::master() && Pstream::parRun())
    {
	CurleDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
    }
    else if (!Pstream::parRun())
    {
	CurleDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
    }
    
    if (Pstream::master() || !Pstream::parRun())
    {
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	//Save timestep for FFT transformation in tau
	scalar tau = probeFreq_*mesh.time().deltaT().value();
        Info << "Executing fft for obs: " << name_ << endl;
	forAll(observers_, iObserver)
	{
	    SoundObserver& obs = observers_[iObserver];
	    
	    autoPtr<List<List<scalar> > > obsFftPtr (obs.fft(tau));
	    
	    List<List<scalar> >& obsFft = obsFftPtr();
	    
  	    if (obsFft[0].size() > 0)
	    {

		fileName fftFile = CurleDir + "/fft-" + name_ + "-" + obs.name() + ".dat";
		
		OFstream fftStream (fftFile);
		fftStream << "Freq p\' spl" << endl;
		
		forAll(obsFft[0], k)
		{
		    fftStream << obsFft[0][k] << " " << obsFft[1][k] << " " << obsFft[2][k] << endl;
		}
		
		fftStream.flush();
	    }
	}
    }
}

void Foam::Curle::execute()
{
    if (!active_)
    {
	return;
    }

    // Create the Curle file if not already created
    makeFile();

    scalar cTime = obr_.time().value();
    
    probeI_++;
    
    if ( mag(probeI_ % probeFreq_) > VSMALL  )
    {
	return;
    }
    else
    {
	if (log_)
	{
	    Info << "Starting acoustics probe" << endl;
	}
	probeI_ = 0.0;
    }
    
    if ( (cTime < timeStart_) || (cTime > timeEnd_))
    {
	return;
    }
    
    correct();
    
    if (Pstream::master() || !Pstream::parRun())
    {
	// time history output
	CurleFilePtr_() << (cTime - timeStart_) << " ";
	
	forAll(observers_, iObserver)
	{
	    const SoundObserver& obs = observers_[iObserver];
	    CurleFilePtr_() << obs.apressure() << " ";
	}
	
	CurleFilePtr_() << endl;
	
	//fft output
	writeFft();
	
	//output to stdio
	if (log_)
	{
	    Info << "Curle acoustic pressure" << endl;
	    forAll(observers_, iObserver)
	    {
		const SoundObserver& obs = observers_[iObserver];
		Info << "Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
	    }
	    Info << endl;
	}
    }
}

bool Foam::Curle::expire()
{
    bool justExpired = false;

    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).expire())
        {
            justExpired = true;
        }

        // Clear merge information
        // if (Pstream::parRun())
        // {
        //     mergeList_[surfI].clear();
        // }
    }

    // true if any surfaces just expired
    return justExpired;
}

bool Foam::Curle::needsUpdate() const
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

bool Foam::Curle::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // Serial: quick and easy, no merging required
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

    // Dimension as fraction of mesh bounding box
    // scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    // if (Pstream::master() && debug)
    // {
    //     Pout<< nl << "Merging all points within "
    //         << mergeDim << " metre" << endl;
    // }

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

        // PatchTools::gatherAndMerge
        // (
        //     mergeDim,
        //     primitivePatch
        //     (
        //         SubList<face>(s.faces(), s.faces().size()),
        //         s.points()
        //     ),
        //     mergeList_[surfI].points,
        //     mergeList_[surfI].faces,
        //     mergeList_[surfI].pointsMap
        // );
    }

    return updated;
}

void Foam::Curle::end()
{
    // Do nothing - only valid on execute
}

void Foam::Curle::timeSet()
{
    // Do nothing - only valid on write
}

void Foam::Curle::write()
{
    // Do nothing - only valid on execute
}

// ************************************************************************* //
