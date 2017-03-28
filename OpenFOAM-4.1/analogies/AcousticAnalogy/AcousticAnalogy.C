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
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    probeI_(0)
{
    this->read(dict);
    this->makeFile();
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
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    c0_(300.0),
    dRef_(-1.0),
    observers_(0),
    probeI_(0)
{
    this->read(dict);
    this->makeFile();
}

void Foam::functionObjects::AcousticAnalogy::makeFile()
{
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
        files().resize(3);
        
        files().set
        (
            2,
            new OFstream
            (
                ResultsDir + "/" + (name() + "-time.dat")
            )
        );
        
        file(2) << "Time" << " ";
        forAll(observers_, iObserver)
        {
            file(2) << observers_[iObserver].name() << "_pFluct ";
        }
        file(2) << endl;
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
        scalar tau = (mesh.time().value() - timeStart_);
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
    
    dict.lookup("probeFrequency") >> probeFreq_;
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("c0") >> c0_;
    
    dict.lookup("dRef") >> dRef_;
    
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
        file(2) << (cTime - timeStart_) << " ";
        
        forAll(observers_, iObserver)
        {
            const SoundObserver& obs = observers_[iObserver];
            file(2) << obs.apressure() << " ";
        }
        
        file(2) << endl;
        
        //fft output
        writeFft();
        
        //output to stdio
        //Log << "Acoustic pressure" << endl;
        //forAll(observers_, iObserver)
        //{
        //    const SoundObserver& obs = observers_[iObserver];
        //    Info << "Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
        //}
    }

    return true;
}


//
//END-OF-FILE
//


