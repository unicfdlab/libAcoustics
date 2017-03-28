#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "SoundObserver.H"

#ifndef AcousticAnalogy_H
#define AcousticAnalogy_H

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;


namespace functionObjects
{

class AcousticAnalogy
:
    public forces
{

protected:

        //-
        label probeFreq_;
        
        //-
        scalar timeStart_;
        
        //-
        scalar timeEnd_;
        
        //- Sonic speed, m/s
        scalar c0_;
        
        //- Domain depth, m for 2D simulations
        scalar dRef_;
        
        //- List of observers
        List<SoundObserver> observers_;
        
        //-
        label probeI_;


protected:

        //-
        void makeFile();
        
        //-
        void writeFft();
        
        //-
        virtual void correct();
        

public:

        //- Runtime type information
        TypeName("AcousticAnalogy");
    
        // Constructors
    
        //- Construct from given Time and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const Time&,
            const dictionary&
        );
    
        //- Construct from given objectRegistry and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );
    
        //- Destructor
        virtual ~AcousticAnalogy();

public:

        //- Read the Curle data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool write();

};

} //namespace functionObjects

} //namespace Foam

#endif

//
//END-OF-FILE
//


