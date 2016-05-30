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

#include "FWH.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

#include "basicThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

    defineTypeNameAndDebug(FWH, 0);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FWH::FWH
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
/*
    fwhSampledSurface_(
		"FWHControlSurface", 
		obr, 
		dict.subDict("controlSurface"),
		loadFromFiles
		),
*/
    fwhSurfaceType_(word::null),
    fwhSurfaceName_(word::null),
//    patchNames_(word::null),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    pName_(word::null),
    UName_(word::null),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    rhoName_(word::null),
    rhoRef_(1.0),
    pRef_(1e5),
    c_(vector::zero),
    FWHFilePtr_(NULL),
    FFOldPtr_(0),
    FFOldOldPtr_(0),
    probeI_(0),
    graphFormat_("raw")
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "Foam::FWH::FWH"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << endl;
    }
    
    Info << "FWH Construct done!!" << endl;

    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FWH::~FWH()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FWH::read(const dictionary& dict)
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

    dict.subDict("controlSurface").lookup("type") >> fwhSurfaceType_;
       
    if (fwhSurfaceType_=="faceSet")
    {
    	dict.subDict("controlSurface").lookup("faceSetName") >> fwhSurfaceName_;
    	
    	Info<<"FwocsWilliams-Hawkings control surface type is FACESET"<<nl;
    }

/*    
    if (fwhSurfaceType_=="triSurface")
    {
    	Info <<"FwocsWilliams-Hawkings control surface type is SAMPLED triSurface"<<nl;
    	
    	fwhSampledSurface_.update();
    	fwhSampledSurface_.print(Info);
    	
    	Info << endl;
    	
    }
*/    

    dict.lookup("probeFrequency") >> probeFreq_;

//    dict.lookup("patchNames") >> patchNames_;
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("c0") >> c0_;
    
    dict.lookup("dRef") >> dRef_;
    
    dict.lookup("pName") >> pName_;
    
    dict.lookup("UName") >> UName_;
    
    dict.lookup("rhoName") >> rhoName_;
    
    dict.lookup("rhoRef") >> rhoRef_;
    
    dict.lookup("pRef") >> pRef_;
    
    dict.lookup("graphFormat") >> graphFormat_;

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
	    scalar Uref = 1.0;
	    obsDict.subDict(oname).lookup("URef") >> Uref;
	    scalar lref = 1.0;
	    obsDict.subDict(oname).lookup("lRef") >> lref;
	    label fftFreq = 1024;
	    obsDict.subDict(oname).lookup("fftFreq") >> fftFreq;
	    observers_.append
	    (
		SoundObserver
		(
		    oname,
		    opos,
		    pref,
		    Uref,
		    lref,
		    fftFreq
		)
	    );
	
	}
	FFOldPtr_.setSize(obsNames.size());
	FFOldOldPtr_.setSize(obsNames.size());
    }    
	Info <<" FWH read done" << endl;

}

void Foam::FWH::correctFaceSetFWH()
{
//    Info << "Correct FaceSet" << endl;

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    
    
    const surfaceScalarField pSf = fvc::interpolate(p);
    const surfaceVectorField uSf = fvc::interpolate(U);    
    
    scalar deltaT = mesh.time().deltaT().value();

    scalar fZcount_ = 0;
    double pfl_ = 0;
    double Uj_ = 0;
    double Un_ = 0;
    scalar r_ = 0;
    
    vector faceNormal_ (0, 0, 0);
    vector l_ (0, 0, 0);
    double FF_ = 0;
    double dFFdT_ = 0;
    double F1_ = 0;
    double Q_ = 0;
    double F2_ = 0;
//    const double psmall = 0.0000001;
    
    const label faceZoneID = mesh.faceZones().findZoneID(fwhSurfaceName_);
    
    if (faceZoneID < 0)
    {
	FatalError
	    << "Don't find faces zone " << fwhSurfaceName_ << endl
	    << exit(FatalError);
    }
    
    const faceZone& facesZone = mesh.faceZones()[faceZoneID];
    
    fZcount_ = facesZone.size();

    reduce(fZcount_, sumOp<scalar>() );
	
    	scalar count = 0;
    	List<scalar>  fc;
    /*
    	forAll(facesZone, I)
        {

		bool goodFace = true;
		label faceID = facesZone[I]; 
		forAll (mesh.boundaryMesh(), Bp)
		{
		    
		    if (isA<processorPolyPatch>(mesh.boundaryMesh()[Bp]))
		    {
//			label nfaces = mesh.boundaryMesh()[Bp].size();
			const scalar startFaces = mesh.boundaryMesh()[Bp].start();
			scalar faceI = startFaces;
//			Pout << "nfaces = " << nfaces << " startFaces = " << startFaces << endl;
//			Pout << "meshBoundary " << mesh.boundaryMesh()[Bp].name() << mesh.boundaryMesh()[Bp].whichFace<< endl;
			forAll(mesh.boundaryMesh()[Bp], Fp)
			{
			    if (faceI == faceID)
			    {
			
			    const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh.boundaryMesh()[Bp] );
			
				Pout << "polyPatch Name "<< mesh.boundaryMesh()[Bp].name() <<" POLYPATCH face "<< faceI << " faceID = " << faceID << endl;
				fc.setSize((count+1)*4);
				fc[0+count*4] = faceID;
			    	fc[1+count*4] = pp.myProcNo();
			    	fc[2+count*4] = pp.neighbProcNo(); 
			    	fc[3+count*4] = 0;
			
			    	count++;
			    	goodFace = false;
			    }
			    faceI++;
			}
		    }
		}	
	}

	List< List<scalar> > gatherFace(Pstream::nProcs());
	gatherFace[Pstream::myProcNo()] = fc;
	Pstream::gatherList(gatherFace);
    
    reduce(count , sumOp<scalar>() );
	
    scalar cc = count;
    Info << cc << endl;
    */
    forAll (observers_, iObs)
    {
    
	deltaT = mesh.time().deltaT().value();
	pfl_ = 0;
	Uj_ = 0;
	Un_ = 0;
	FF_ = 0;
	dFFdT_ = 0;
	F1_ = 0;
	Q_ = 0;
	F2_ = 0;
	count = 0;
	SoundObserver& obs = observers_[iObs];

//	Pout <<"facezone size = " << facesZone.size() << endl;
	
	forAll(facesZone, I)
        {

		bool goodFace = true;
		label faceID = facesZone[I]; 

		forAll (mesh.boundaryMesh(), Bp)
		{
		    
		    if (isA<processorPolyPatch>(mesh.boundaryMesh()[Bp]))
		    {
//			label nfaces = mesh.boundaryMesh()[Bp].size();
			const scalar startFaces = mesh.boundaryMesh()[Bp].start();
			scalar faceI = startFaces;

			forAll(mesh.boundaryMesh()[Bp], Fp)
			{
			    if ( (faceI == faceID) )// &(count <= cc/2)  )
			    {
			    	//Pout << "polyPatch Name "<< mesh.boundaryMesh()[Bp].name() <<" POLYPATCH face "<< faceI << " faceID = " << faceID << endl;
			    	count++;
			    	goodFace = false;
			    }
			    faceI++;
			}
		    }
		}	


		if (goodFace == true)
		{
			pfl_ = (pSf[faceID] - pRef_)*rhoRef_;
		        Uj_ = uSf[faceID].component(1);

		        faceNormal_ = mesh.Sf()[faceID] / mesh.magSf()[faceID];
		    	Un_ = uSf[faceID].component(0) * faceNormal_.component(0);
			l_ = obs.position() - mesh.Cf()[faceID];
			r_ = mag(l_);
			
			F1_ += ((pfl_* faceNormal_.component(1) * l_.component(1) + rhoRef_ * Uj_ * Un_ * l_.component(1))/(r_));
			Q_  += ((rhoRef_ * Un_ / r_));
			F2_ += ((pfl_* faceNormal_.component(1) * l_.component(1) + rhoRef_ * Uj_ * Un_ * l_.component(1))/ (r_ * r_));
		}
        }
        
	reduce(count , sumOp<scalar>() );
	
	reduce(F1_ , sumOp<scalar>() );
	reduce(F2_ , sumOp<scalar>() );
	reduce(Q_ , sumOp<scalar>() );

//	Info << "count duplicateMesh = " << count << endl;
//	Pout << "Reduce done" <<  endl;

	if (Pstream::master() || !Pstream::parRun())
	{
	    FF_ = Q_ + F1_ / c0_;
/*
	    Info << "sumF1_ = " << F1_  << endl;
	    Info << "sumQ_ = " << Q_  << endl;
	    Info << "sumF2_ = " << F2_  << endl;
	    Info << "FaceZoneCount_ = " << fZcount_  << endl;
	    Info << "FF_ = " << FF_ << endl;
*/		

	    if (FFOldPtr_[iObs].empty())
	    {
		FFOldPtr_[iObs].set
		(
		    new scalar(FF_)
		);
	    }
	    else
	    {
//		Info << "deltat_ = " << deltaT  << endl;
		dFFdT_ = (FF_ - FFOldPtr_[iObs]()) / deltaT;
		if (FFOldOldPtr_[iObs].empty())
		{
		    //first order scheme
		    dFFdT_ = (FF_ - FFOldPtr_[iObs]()) / deltaT;
		
		    FFOldOldPtr_[iObs].set
		    (
			FFOldPtr_[iObs].ptr()
		    );
		
		    FFOldPtr_[iObs].reset
		    (
			new scalar(FF_)
		    );
		}
	        else
		{
	    	    //second order scheme (BDF)
	    	    dFFdT_ = (3.0*FF_ - 4.0*FFOldPtr_[iObs]() + FFOldOldPtr_[iObs]()) / 2.0 / deltaT;
	    
	            FFOldOldPtr_[iObs].reset
		    (
		    	FFOldPtr_[iObs].ptr()
		    );
		    FFOldPtr_[iObs].reset
		    (
			new scalar(FF_)
		    );
	        }
	    }
	}
//	Info << "dFdT = " << dFFdT_ << endl;
	scalar oap = (dFFdT_ + F2_)/(4 * Foam::constant::mathematical::pi);
	    
	if (dRef_ > 0.0)
	{
	    oap /= dRef_;
	}
	
	obs.apressure(oap);
	obs.atime(mesh.time().value());
	
    }
    
}

/*
void Foam::FWH::correctSampledFWH()
{
    Info << "Correct SAMPLED" << endl; 

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    //sign '-' needed to calculate force, which exerts fluid by solid
//    vector F	(0.0, 0.0, 0.0);
//    vector dFdT (0.0, 0.0, 0.0);
    scalar deltaT = mesh.time().deltaT().value();

    float pfl_;
    float Uj_;
    float Un_;
    scalar r_;
    vector faceNormal_;
    vector l_;
    scalar F1_;
    scalar F2_;
    scalar FF_;
    scalar dFFdT_ = 0;
    scalar Q_;

     const pointField& FWHFaces = fwhSampledSurface_.points();
     Info << "FWH size " << FWHFaces.size() << endl;
//    tmp<scalarField> ps = fwhSampledSurface_.sample(p);

    tmp<scalarField>	ps = fwhSampledSurface_.sample(p);
     Info << "size ps " << ps().size() << endl;

     tmp<vectorField>	Us = fwhSampledSurface_.sample(U);
     Info << "size Us " << Us().size() << endl;
//    tmp<vectorField> Us = fwhSampledSurface_.sample(U);
        
 

    forAll (observers_, iObs)
    {

	SoundObserver& obs = observers_[iObs];
	F1_ = 0;
	F2_ = 0;
	Q_ = 0;
	FF_ = 0;
	Info << "obs" << obs.position() << endl;

        forAll(FWHFaces, I)
        {
		faceNormal_ = fwhSampledSurface_.Sf()[I] / fwhSampledSurface_.magSf()[I];
//		Info << faceNormal_ << endl;
		pfl_ = (ps()[I] - pRef_)*rhoRef_;
		Uj_ = Us()[I].component(1);
		Un_ = Us()[I].component(0) * faceNormal_.component(0);

	    	l_ = obs.position() - fwhSampledSurface_.Cf()[I];
		r_ = mag(l_);

//		Info << Un_ << endl;		    

	    	F1_ += (pfl_* faceNormal_.component(1) * l_.component(1) + rhoRef_ * Uj_ * Un_ * l_.component(1))/(r_);
	    	Q_  += rhoRef_ * Un_ / r_;
	    
	F2_ += (pfl_* faceNormal_.component(1) * l_.component(1) + rhoRef_ * Uj_ * Un_ * l_.component(1))/ (r_ * r_);

	}

	Info << "F1 = " << F1_ << endl;
	Info << "Q = " << Q_ << endl;
	Info << "F2 = " << F2_ << endl;
		
	FF_ = Q_ + F1_ / c0_;
	
	if (Pstream::master() || !Pstream::parRun())
	{
	//calculate dFdT and store old values
	
	    if (FFOldPtr_.empty())
	    {
		FFOldPtr_.set
		(
		    new scalar(FF_)
		);
	    }
	    else
	    {
		if (FFOldOldPtr_.empty())
		{
		    //first order scheme
		    dFFdT_ = (FF_ - FFOldPtr_()) / deltaT;
		
		    FFOldOldPtr_.set
		    (
			FFOldPtr_.ptr()
		    );
		
		    FFOldPtr_.reset
		    (
			new scalar(FF_)
		    );
		}
		else
		{
		    //second order scheme (BDF)
		    dFFdT_ = (3.0*FF_ - 4.0*FFOldPtr_() + FFOldOldPtr_()) / 2.0 / deltaT;
		
		    FFOldOldPtr_.reset
		    (
			FFOldPtr_.ptr()
		    );
		
		    FFOldPtr_.reset
		    (
			new scalar(FF_)
		    );
		}
	    }
	}

	
	scalar oap = (dFFdT_ + F2_)/(4 * Foam::constant::mathematical::pi);
	    
	if (dRef_ > 0.0)
	{
	    oap /= dRef_;
	}
	 obs.apressure(oap);
	 
    }
}

*/

void Foam::FWH::makeFile()
{

    fileName FWHDir;

    if (Pstream::master() && Pstream::parRun())
    {
	FWHDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
	mkDir(FWHDir);
    }
    else if (!Pstream::parRun())
    {
	FWHDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
	mkDir(FWHDir);
    }
    else
    {
    }

    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
	// Create the FWH file if not already created
	if (FWHFilePtr_.empty())
	{
	    // Open new file at start up
	    FWHFilePtr_.reset
	    (
		new OFstream
		(
		    FWHDir + "/" + (name_ + "-time.dat")
		)
	    );
	    
	    writeFileHeader();
	}
    }
}


void Foam::FWH::writeFileHeader()
{
    if (FWHFilePtr_.valid())
    {
        FWHFilePtr_()
            << "Time" << " ";
	
        forAll(observers_, iObserver)
        {
	    FWHFilePtr_() << observers_[iObserver].name() << "_pFluct ";
        }

        FWHFilePtr_()<< endl;
    }
}

void Foam::FWH::writeFft()
{
    fileName FWHDir;

    if (Pstream::master() && Pstream::parRun())
    {
	FWHDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
    }
    else if (!Pstream::parRun())
    {
	FWHDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
    }
    
    if (Pstream::master() || !Pstream::parRun())
    {
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	scalar tau = (mesh.time().value() - timeStart_);
	    
//	autoPtr<List<List<scalar> > > obsFftPtr (obs.fft(tau));
//	List<List<scalar> >& obsFft = obsFftPtr();

    	Info << "Executing fft for obs: " << name_ << endl;
	
	//    Info << "Creating noiseFFT" << endl;
	forAll(observers_, iObserver)
	{
	    SoundObserver& obs = observers_[iObserver];
	    fileName fftFile = FWHDir + "/fft-" + name_ + "-" + obs.name() + ".dat";
	    const scalarField p(obs.p());
	    const scalarField t(obs.time());

	    noiseFFTDriver nfft(probeFreq_*mesh.time().deltaT().value(), p);
		
	    if (nfft.size() >= 2)
	    {
	        nfft -= obs.pref();
	        scalar Nn = 0;
	        while (nfft.size() >= pow(2,Nn))
	        {
			Nn++;
		}
		scalar N = pow(2,Nn-1);
		
		graph rmsPf (nfft.RMSmeanPf(t, N, min(nfft.size()/N, 100)));
		rmsPf.write(FWHDir + "/"+ name_ + "-" + graph::wordify(rmsPf.title() + "-" + obs.name()), graphFormat_);

		graph Pt (nfft.pt(t));
		Pt.write(FWHDir + "/"+ name_ + "-" + graph::wordify(Pt.title() + "-" + obs.name()), graphFormat_);
		
		graph spl (nfft.spl(rmsPf));
		spl.write(FWHDir + "/"+ name_ + "-" + graph::wordify(spl.title() + "-" + obs.name()), graphFormat_);
	    }
	    else
	    {
	        Info << "Data size < 2" << endl;
	    }

/*
	    if (obsFft[0].size() > 0)
	    {
		Info << "Executing fft for obs: " << name_ << endl;
		fileName fftFile = FWHDir + "/fft-" + name_ + "-" + obs.name() + ".dat";
		
		OFstream fftStream (fftFile);
		fftStream << "Freq"<< tab << "p'" << tab << "SPL" << tab << "St" << endl;
		
		forAll(obsFft[0], k)
		{
		    fftStream << obsFft[0][k] << tab << obsFft[1][k] << tab << obsFft[2][k] << tab << obsFft[3][k] << endl;
		}
		
		fftStream.flush();
	    }
	    */
	}
    }
}

void Foam::FWH::execute()
{
    if (!active_)
    {
	return;
    }

    // Create the FWH file if not already created
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
	    Info << "Starting acoustics probe for FW-H analogy" << endl;
	}
	probeI_ = 0.0;
    }
    
    if ( (cTime < timeStart_) || (cTime > timeEnd_))
    {
	return;
    }
    
    
     if  (fwhSurfaceType_=="faceSet")
    {
	correctFaceSetFWH();
    }
    else if (fwhSurfaceType_=="triSurface")
    {
    	correctSampledFWH();
    }


    //Info << "correct Done!!!" << endl;

    //correct();
    
    if (Pstream::master() || !Pstream::parRun())
    {
	// time history output
	FWHFilePtr_() << (cTime - timeStart_) << " ";
	
	forAll(observers_, iObserver)
	{
	    const SoundObserver& obs = observers_[iObserver];
	    FWHFilePtr_() << obs.apressure() << " ";
	}
	
	FWHFilePtr_() << endl;
	
	//fft output
	writeFft();
	
	//output to stdio
	if (log_)
	{
	    Info << "FW-H acoustic pressure" << endl;
	    forAll(observers_, iObserver)
	    {
		const SoundObserver& obs = observers_[iObserver];
		Info << "Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
	    }
	    Info << endl;
	}
	
    }
}


void Foam::FWH::end()
{
    // Do nothing - only valid on execute
}

void Foam::FWH::timeSet()
{
    // Do nothing - only valid on execute
}


void Foam::FWH::write()
{
    // Do nothing - only valid on execute
}




// ************************************************************************* //
