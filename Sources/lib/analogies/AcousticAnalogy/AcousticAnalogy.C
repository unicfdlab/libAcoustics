#include "AcousticAnalogy.H"

#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(AcousticAnalogy, 0);
//    addToRunTimeSelectionTable(functionObject, AcousticAnalogy, dictionary);
}
}

Foam::functionObjects::AcousticAnalogy::AcousticAnalogy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces
    (
        name,
        runTime,
        dict
    ),
    analogyOutPtr_(nullptr),
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    writeFft_(true),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    probeI_(0)
{
}

Foam::functionObjects::AcousticAnalogy::AcousticAnalogy
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    ),
    analogyOutPtr_(nullptr),
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    writeFft_(true),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    probeI_(0)
{
}

Foam::functionObjects::AcousticAnalogy::~AcousticAnalogy()
{}

void Foam::functionObjects::AcousticAnalogy::makeFile()
{
    if (Pstream::master())
    {
        if(analogyOutPtr_.valid())
        {
            return;
        }
    }
    
    fileName ResultsDir;
    
    if (Pstream::master() && Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
        mkDir(ResultsDir);
    }
    else if (!Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
        mkDir(ResultsDir);
    }
    else
    {
    }
    
    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
        
        analogyOutPtr_.set
        (
            new OFstream
            (
                ResultsDir + "/" + (name() + "-time.dat")
            )
        );
        
        analogyOutPtr_() << "Time" << " ";
        forAll(observers_, iObserver)
        {
            analogyOutPtr_() << observers_[iObserver].name() << "_pFluct ";
        }
        analogyOutPtr_() << endl;
    }
}

void Foam::functionObjects::AcousticAnalogy::writeFft()
{
    fileName ResultsDir;
        
    if (Pstream::master() && Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
    }
    else if (!Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
    }
    
    if (Pstream::master() || !Pstream::parRun())
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        scalar	tau;
        if (mesh.time().startTime().value() > timeStart_)
        {
    	    tau = (mesh.time().value() - mesh.time().startTime().value());
    	}
    	else
    	{
    	    tau = (mesh.time().value() - timeStart_);
    	}    
    	
        forAll(observers_, iObserver)
        {
            SoundObserver& obs = observers_[iObserver];
            autoPtr<List<List<scalar> > > obsFftPtr (obs.fft(tau));
            
            List<List<scalar> >& obsFft = obsFftPtr();
            
            if (obsFft[0].size() > 0)
            {
                Log << "Executing fft for obs: " << obs.name() << endl;
                fileName fftFile = ResultsDir + "/fft-" + name() + "-" + obs.name() + ".dat";
                
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

void Foam::functionObjects::AcousticAnalogy::correct()
{
}

bool Foam::functionObjects::AcousticAnalogy::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    Info << "Reading analogy settings" << endl;

    dict.lookup("probeFrequency") >> probeFreq_;
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("writeFft") >> writeFft_;
    
    dict.lookup("c0") >> c0_;
    
//    dict.lookup("U0") >> U0_;
    
    dict.lookup("dRef") >> dRef_;
    
    //Ask for rhoRef again, because libforces sometimes do not
    dict.lookup("rhoInf") >> rhoRef_;
    
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
    
    this->makeFile();
    
    return true;
}

bool Foam::functionObjects::AcousticAnalogy::execute()
{
    return true;
}

bool Foam::functionObjects::AcousticAnalogy::write()
{
    //use forces library to calculate forces acting on patch
    if (!forces::write())
    {
        return false;
    }
    
    scalar cTime = obr_.time().value();
    scalar stTime = obr_.time().startTime().value();
    
    probeI_++;
    
    if ( mag(probeI_ % probeFreq_) > VSMALL  )
    {
        return true;
    }
    else
    {
        Log << "Starting acoustics probe" << endl;
        probeI_ = 0.0;
    }
    
    if ( (cTime < timeStart_) || (cTime > timeEnd_))
    {
        return true;
    }
    
    correct();
    
    if (Pstream::master() || !Pstream::parRun())
    {
        // time history output
        if (stTime > timeStart_)
        {
            analogyOutPtr_() << (cTime - stTime) << " ";
        }
        else
        {
            analogyOutPtr_() << (cTime - timeStart_) << " ";
        }
        
        forAll(observers_, iObserver)
        {
            const SoundObserver& obs = observers_[iObserver];
            analogyOutPtr_() << obs.apressure() << " ";
        }
        
        analogyOutPtr_() << endl;
        
        //fft output
        if (writeFft_ == true)
        {
            writeFft();
        }
        
        //output to stdio
        Log << "Acoustic pressure" << endl;
        forAll(observers_, iObserver)
        {
            const SoundObserver& obs = observers_[iObserver];
            Info << iObserver << " Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
        }
    }
    
    return true;
}


//
//END-OF-FILE
//


