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
#include "fvcDdt.H"
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

    defineTypeNameAndDebug(FWH, 0);
}

Foam::scalar Foam::FWH::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::scalarField> Foam::FWH::normalStress(const sampledSurface& surface) const
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    //Info<<"    Normal stress p-field was read"<<nl; 
    tmp<Field<scalar> > pSampled;
    
    //Sample pressure field and obtain difference with pressure in undisturbed flow
    pSampled = sampleOrInterpolate<scalar>(p , surface) - pInf_;
    //Info<<"    Normal stress p-field was sampled"<<nl; 
    if (p.dimensions() == dimPressure)
    {
	//return tmp<scalarField>
	//(
	//    new scalarField(pPatch)
	//);
    }
    else
    {
      tmp<Field<scalar> > rhoField = surfaceDensity(surface);
      pSampled() *= rhoField();
    }

    return pSampled;
}

Foam::tmp<Foam::vectorField> Foam::FWH::surfaceVelocity(const sampledSurface& surface) const
{
    const volVectorField& U = obr_.lookupObject<volVectorField>("U");
 
    tmp<Field<vector> > USampled;
    
    USampled = sampleOrInterpolate<vector>(U , surface);

    return USampled;
}

Foam::tmp<Foam::scalarField> Foam::FWH::surfaceDensity(const sampledSurface& surface) const
{
  tmp<Field<scalar> > rhoSampled;

  if (rhoRef_ < 0) //density in volScalarField
    {
      volScalarField pRho = obr_.lookupObject<volScalarField>(rhoName_);     
      rhoSampled = sampleOrInterpolate<scalar>(pRho, surface);
    }
  else //density is constant
    {
      Info<<"Incompressible crash"<<nl;
      const fvMesh& mesh_ = refCast<const fvMesh>(obr_);
      volScalarField rhoTemp
        (
	 IOobject
	 (
	  name_ + ":area",
	  mesh_.time().timeName(),
	  mesh_,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh_,
	 dimensionedScalar("temp", dimless, rhoRef_)
	 );
      
      rhoSampled = sampleOrInterpolate<scalar>(rhoTemp, surface); //= rhoRef_;
    }

  return rhoSampled;
}

//May be not need at all
//Trade off between this implementation of the derivative operation
//and using <type> dotProduct(type &)

Foam::tmp<Foam::scalarField> Foam::FWH::dotNormalStress(const sampledSurface& surface) const
{
    const volScalarField& dpdt = Foam::fvc::ddt(obr_.lookupObject<volScalarField>(pName_));
 
    tmp<Field<scalar> > pSampled;
    
    //Sample pressure field and obtain difference with pressure in undisturbed flow
    pSampled = sampleOrInterpolate<scalar>(dpdt , surface);

    if (rhoRef_ < 0) //density in volScalarField
      {
	tmp<Field<scalar> > rhoField = surfaceDensity(surface);
	pSampled() *= rhoField();
      }
    else //density is constant
      {
	pSampled() *= rhoRef_;
      }

    return pSampled;
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
    fwhSurfaceType_(word::null),
    fwhFaceSetName_(word::null),
    faceZoneID_(-1),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    pName_(word::null),
    pInf_(0),
    Ufwh_(vector::zero),
    c0_(340.0),
    dRef_(-1.0),
    observers_(0),
    rhoName_(word::null),
    rhoRef_(1.0),
    rhoInf_(1.0),
    c_(vector::zero),
    FWHFilePtr_(NULL),
    F1_(0.0),
    F2_(0.0),
    VOldPtr_(NULL),
    VOldOldPtr_(NULL),
    SOldPtr_(NULL),
    SOldOldPtr_(NULL),
    mergeList_(),
    probeI_(0)
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
    
    dict.lookup("probeFrequency") >> probeFreq_;

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);
  
    //Istream& is = dict.lookup("surfaces");
    //dictionary dc(is);
    //dc.readIfPresent("type", fwhSurfaceType_);

    dict.lookup("fwhSurfaceType") >> fwhSurfaceType_;
    
    if ( fwhSurfaceType_=="faceSet")
    {
      dict.lookup("faceSetName") >> fwhFaceSetName_;
      Info<< "Function object "<< name_<<":" << nl;
      Info<< "    FwocsWilliams-Hawkings analogy control surface is faceSet "
	  << fwhFaceSetName_ << nl << nl;

      const fvMesh& mesh = refCast<const fvMesh>(obr_);
      faceZoneID_ = mesh.faceZones().findZoneID(fwhFaceSetName_);
      
      if (faceZoneID_ < 0)
	{
	  FatalError
	    << "Foam::FWH::read(const dictionary& dict): can not find faceZone "
	    << fwhFaceSetName_ << endl
	    << exit(FatalError);
	}

    }
    else if( fwhSurfaceType_ == "sampled")
      {
	PtrList<sampledSurface> newList
	  (
	   dict.lookup("surfaces"),
	   sampledSurface::iNew(mesh_)
	   );
	
	controlSurfaces_.transfer(newList);
	//Parallel fix as it was implemented in sampledSuraces class
	if (Pstream::parRun())
	  {
	    mergeList_.setSize(controlSurfaces_.size());
	  }
	
	// Ensure all surfaces and merge information are expired
	expire();
	
	if (controlSurfaces_.size())
	  {
	    Info<< "Function object "<< name_<<":" << nl;
	    Info<< "    Reading FwocsWilliams-Hawkings analogy control surface description:" << nl;
	    forAll(controlSurfaces_, surfI)
	      {
		Info<< "        " <<  controlSurfaces_.operator[](surfI).name() << nl;        
	      }
	    Info<< endl;
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
	
      }
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("c0") >> c0_;
    
    dict.lookup("dRef") >> dRef_;

    dict.lookup("pName") >> pName_;

    dict.lookup("pInf") >> pInf_;

    dict.lookup("Ufwh") >> Ufwh_;
    
    dict.lookup("rhoName") >> rhoName_;
    
    dict.lookup("rhoRef") >> rhoRef_;

    dict.lookup("rhoInf") >> rhoInf_;
    
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

Foam::vector Foam::FWH::dotProduct(const vector& F)
{
    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    vector dFdT (0.0, 0.0, 0.0);

    scalar deltaT = probeFreq_*mesh_.time().deltaT().value();

    if (VOldPtr_.empty())
      {
	VOldPtr_.set
	  (
	   new vector(F)
	  );
      }
    else
      {
	if (VOldOldPtr_.empty())
	  {
	    //first order scheme
	    dFdT = (F - VOldPtr_()) / deltaT;
	    
	    VOldOldPtr_.set
	      (
	       VOldPtr_.ptr()
	       );
	    
	    VOldPtr_.reset
	      (
	       new vector(F)
	       );
	  }
	else
	  {
	    //second order scheme (BDF)
	    dFdT = (3.0*F - 4.0*VOldPtr_() + VOldOldPtr_()) / 2.0 / deltaT;
	    
	    VOldOldPtr_.reset
	      (
	       VOldPtr_.ptr()
	       );
	    
	    VOldPtr_.reset
	      (
	       new vector(F)
	       );
	  }
      }

    return dFdT;
}

Foam::scalar Foam::FWH::dotProduct(const scalar& S)
{
    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    scalar dSdT (0.0);

    scalar deltaT = probeFreq_*mesh_.time().deltaT().value();

    if (SOldPtr_.empty())
      {
	SOldPtr_.set
	  (
	   new scalar(S)
	  );
      }
    else
      {
	if (SOldOldPtr_.empty())
	  {
	    //first order scheme
	    dSdT = (S - SOldPtr_()) / deltaT;
	    
	    SOldOldPtr_.set
	      (
	       SOldPtr_.ptr()
	       );
	    
	    SOldPtr_.reset
	      (
	       new scalar(S)
	       );
	  }
	else
	  {
	    //second order scheme (BDF)
	    dSdT = (3.0*S - 4.0*SOldPtr_() + SOldOldPtr_()) / 2.0 / deltaT;
	    
	    SOldOldPtr_.reset
	      (
	       SOldPtr_.ptr()
	       );
	    
	    SOldPtr_.reset
	      (
	       new scalar(S)
	       );
	  }
	  }

    return dSdT;
}

Foam::scalar Foam::FWH::dotProduct(dScalar& var)
{
    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    scalar ddT (0.0);

    scalar deltaT = probeFreq_*mesh_.time().deltaT().value();
   
    if (var.old.empty())
      {	
	var.old.set
	  (
	   new scalar(var.value)
	  );
      }
    else
      {
	if (var.oldOld.empty())
	  {
	    //first order scheme
	    ddT = (var.value - var.old()) / deltaT;
	    var.oldOld.set
	      (
	       var.old.ptr()
	       );
	    
	    var.old.reset
	      (
	       new scalar(var.value)
	       );
	  }
	else
	  {
	    //second order scheme (BDF)
	    ddT = (3.0*var.value - 4.0*var.old() + var.oldOld()) / 2.0 / deltaT;
	    
	    var.oldOld.reset
	      (
	       var.old.ptr()
	       );
	    
	    var.old.reset
	      (
	       new scalar(var.value)
	       );
	  }
      }

    return ddT;
}

void Foam::FWH::correct()
{   
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    scalar dF1dT    (0.0);
    scalar dF2dT    (0.0);
    
    forAll (observers_, iObs)
      {
	SoundObserver& obs = observers_[iObs];
	
	if( fwhSurfaceType_ == "sampled")
	  {
	    //calling a function to update all sampledSurfaces
	    //without it everything related will be empty
	    update();
       
	    //working with sampled surfaces
	    Info<<"    Surfaces updated"<<nl;
       
	    forAll(controlSurfaces_, surfI)
	      {
		sampledSurface& s = controlSurfaces_.operator[](surfI);
		label surfaceSize = s.Cf().size();
		Info<<"    Surface OK"<<nl;
		scalarField pS = normalStress(s);
		Info<<"    Pressure OK"<<gSum(pS)<<nl;
		vectorField uS = surfaceVelocity(s);
		Info<<"    Velocity OK"<<gSum(uS)<<nl;
		scalarField rhoS = surfaceDensity(s);
		Info<<"    Density OK"<<gSum(rhoS)<<nl;

		forAll(s.Cf(), i)
		  {
		    //Vector from observer to center of each cells
		    vector n = s.Sf()[i];
		    scalar dS = mag(n);
		    vector n_i = n/dS;
		    vector x = s.Cf()[i] - obs.position();
		    vector y = obs.position() - s.Cf()[i];
		    //Calculate distance
		    scalar r = mag(x);
		    vector x_i = x/r;
		    vector y_i = y/r;
		    //Calculate Mr
		    scalar Mr = mag((Ufwh_/c0_/r)&x_i);

		    //parallel: integrate locally for each processor
		    //		    F1 += (
		    //	   )*dS;
		
		    //		    F2 += (
		    //	   )*dS; // rhoS - rhoInf_ changed to isentropic equation

		    F1_.value += (
				  (
				   pS[i]*n_i + (n_i&(rhoS[i]*uS[i]) )*(uS[i] - Ufwh_)
				   )&x_i
				  )
		      /(r - r*Mr)
		      *dS;
				      
		    F2_.value += (
				  (
				   ( rhoInf_*uS[i] + (1.4*pS[i]/(c0_*c0_))*(uS[i] - Ufwh_) )
				   /
				   (r - r*Mr)
				   )&n_i
				  )*dS; // rhoS - rhoInf_ changed to isentropic equation
		    
		    if(false)
		      {
			Info<<nl
			    <<"FWH face integration loop debug info >>"<<nl
			    <<"    n = "<< n<<nl
			    <<"    n_i = "<< n_i<<nl
			    <<"    x = "<< x <<", y = "<< y << nl
			    <<"    x_i = "<< x_i <<", y_i = "<< y_i << nl
			    <<"    r = "<< r << nl
			    <<"    Mr = "<< Mr << nl <<endl;
		      }
		  }
		
		//parallel: add integral from each processor
		reduce (F1_.value, sumOp<scalar>());
		reduce (F2_.value, sumOp<scalar>());
		//not using gSum
		
		Info<<s.name()<<", sampled integrals F1="<<F1_.value<<" F2="<<F2_.value<<nl;
	      }
	  }
	else if (fwhSurfaceType_ == "faceSet")
	  {
	    /*	    const volVectorField& U = obr_.lookupObject<volVectorField>("U");
	    const volScalarField& rho = obr_.lookupObject<volScalarField>("rho");
	    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
	    
	    const surfaceScalarField pSf = fvc::interpolate(p);
	    const surfaceVectorField uSf = fvc::interpolate(U);
	    const surfaceScalarField rhoSf = fvc::interpolate(rho);
	    
	    const faceZone& facesZone = mesh.faceZones()[faceZoneID_];
	    int zoneSize = facesZone.size();
	    
	    reduce(zoneSize, sumOp<scalar>() );

	    scalar counter = 0;
	    List<scalar>  faceCenters;

	    forAll(facesZone, idx)
	      {

		bool goodFace = true;
		label faceID = facesZone[idx]; 

		forAll (mesh.boundaryMesh(), patchIdx)
		  {
		    
		    if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchIdx]))
		      {
			//			label nFaces = mesh.boundaryMesh()[patchIdx].size();
			const scalar startFace = mesh.boundaryMesh()[patchIdx].start();
			scalar faceIdx = startFace;
			forAll(mesh.boundaryMesh()[patchIdx], patchFaceIdx)
			  {
			    if ( (faceIdx == faceID) )// &(count <= cc/2)  )
			      {
				//Pout << "polyPatch Name "<< mesh.boundaryMesh()[Bp].name()
				//<<" POLYPATCH face "<< faceI << " faceID = " << faceID << endl;
			    	counter++;
			    	goodFace = false;
			      }
			    faceIdx++;
			  }
		      }
		  }
		if (goodFace == true)
		  {
		    vector faceNormal_ = mesh.Sf()[faceID];
		    vector Uf = uSf[faceID];
		    scalar Pf = pSf[faceID];
		    scalar Rhof = rhoSf[faceID];
		    
		    F1 += Pf*faceNormal_ + 
		      (Rhof*Uf)*( (Uf - Ufwh_)&faceNormal_ );
		    
		    F2 += (rhoInf_*Uf + (Pf/(c0_*c0_))*(Uf - Ufwh_))&faceNormal_ ; //rhoS - rhoInf_
		    
		  }
	    */
	    F1_.value += 0.0; 
	
	    F2_.value += 0.0; //rhoS - rhoInf_
		
	  }

	if (Pstream::master() || !Pstream::parRun())
	  {
	    SoundObserver& obs = observers_[iObs];
	    //calculate dFdT and store old values
	    scalar coeff1 = 1. / 4. / Foam::constant::mathematical::pi;
	    
	    dF1dT = dotProduct(F1_);
	    dF2dT = dotProduct(F2_);
	    
	    //dFdT = F;
	    
	    scalar oap = (dF1dT/c0_ + dF2dT ) * coeff1;
	    
	    if (dRef_ > 0.0)
	      {
		oap /= dRef_;
	      }
	    
	    obs.apressure(oap); //appends new calculated acoustic pressure
	    
	    //noiseFFT addition
	    obs.atime(mesh.time().value());

	    F1_.value = 0.0;
	    F2_.value = 0.0;
	    dF1dT = 0.0;
	    dF2dT = 0.0;

	  }
      }
}
    /*forAll (observers_, iObs)
    {
      SoundObserver& obs = observers_[iObs];


    // Surface integral - loop over all patches
    forAll(controlSurfaces_, surfI)
    {
      sampledSurface& s = controlSurfaces_.operator[](surfI);
      //scalarField sampledPressure = normalStress(s);

      // Surface area vector and face center at patch
      vectorField Sf = s.Sf();
      vectorField Cf = s.Cf();

      // Normal direction vector 
      vectorField n = -Sf/mag(Sf);

      // Pressure field and time derivative at patch
      scalarField pp = normalStress(s);
      //scalarField dpdtp = dpdt.boundaryField()[patchI];

      // Lighthill tensor on patch
      //tensorField Tijp = Tij.boundaryField()[patchI];
      //tensorField dTijdtp = dTijdt.boundaryField()[patchI];

      // Distance surface-observer
      scalarField r = mag(obs.position() - Cf);
      vectorField l = (obs.position() - Cf) / r;

      // Calculate pressure fluctuations
      pPrime += coeff * gSum
	(
	 (
	  (l*n)
	  && 
	  (
	   (dpdtp*I - dTijdtp) / (cRef_*r)
	   + (pp*I - Tijp) / sqr(r)
	   )
	  )
	 * mag(Sf)
	 );
    }

	//obs.pPrime(pPrime);
	}*/


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

void Foam::FWH::calcDistances()
{
    if (!active_)
    {
	return;
    };
    
    vectorField ci;
    scalar ni;
    
    if( fwhSurfaceType_ == "sampled")
      {

	update();

	forAll(controlSurfaces_, surfI)
	  {
	    sampledSurface& s = controlSurfaces_.operator[](surfI);
	    ci = s.Cf();
	    ni = scalar(ci.size());
	  }
	
	reduce (ni, sumOp<scalar>());
	
	c_ = gSum(ci) / ni;
      }
    else if( fwhSurfaceType_ == "faceSet")
      {
       const fvMesh& mesh = refCast<const fvMesh>(obr_);
       const faceZone& facesZone = mesh.faceZones()[faceZoneID_];
       ni = facesZone.size();
       reduce(ni, sumOp<scalar>() );
       forAll(facesZone, idx)
	 {
	   label faceID = facesZone[idx];
	   ci += mesh.Cf()[faceID];
	 }

       c_ = gSum(ci)/ni;
       
      };
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

		fileName fftFile = FWHDir + "/fft-" + name_ + "-" + obs.name() + ".dat";
		
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
	    Info << "FWH acoustic pressure" << endl;
	    forAll(observers_, iObserver)
	    {
		const SoundObserver& obs = observers_[iObserver];
		Info << "Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
	    }
	    Info << endl;
	}
    }
}

bool Foam::FWH::expire()
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

bool Foam::FWH::needsUpdate() const
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

bool Foam::FWH::update()
{
  //Actually the update() function copy-pasted from libSampling
  
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

void Foam::FWH::end()
{
    // Do nothing - only valid on execute
}

void Foam::FWH::timeSet()
{
    // Do nothing - only valid on write
}

void Foam::FWH::write()
{
    // Do nothing - only valid on execute
}

Foam::scalar Foam::FWH::mergeTol()
{
    return mergeTol_;
}

Foam::scalar Foam::FWH::mergeTol(const scalar tol)
{
    scalar oldTol = mergeTol_;
    mergeTol_ = tol;
    return oldTol;
}

bool Foam::FWH::checkNormal( vector c, vector p, vector n) const
{
    vector cp_ (0.0, 0.0, 0.0);
    
    cp_ = c - p;
    cp_ /= mag(cp_);
    
    if ((cp_ & n) > 0)
    {
	return true;
    }
    else
    {
	return false;
    }
}

// ************************************************************************* //
