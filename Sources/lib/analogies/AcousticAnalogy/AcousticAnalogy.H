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
        autoPtr<OFstream> analogyOutPtr_;

        //-
        label probeFreq_;

        //-
        scalar timeStart_;

        //-
        scalar timeEnd_;

        //-
        bool writeFft_;

        //- Sonic speed, m/s
        scalar c0_;

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
        
        //-
        virtual void updateSurfaceCf();

public:

    template<class T>
    class tdPrimitive
    {

    private:

        //-
        bool varTimeStep_;

        //-
        scalar prevt_;

        //-
        scalar deltat_;

        //-
        List<T> new_;

        //-
        List<T> old_;

        //-
        List<T> oldOld_;

    public:

        //-
        tdPrimitive(T initv, scalar tprev)
        :
            varTimeStep_(false),
            prevt_(tprev),
            deltat_(tprev),
            new_(1, initv),
            old_(0),
            oldOld_(0)
        {
        };

        //-
    
        tdPrimitive()
        :
            varTimeStep_(false),
            prevt_(0.0),
            deltat_(GREAT),
            new_(1, pTraits<T>::zero),
            old_(0),
            oldOld_(0)
        {
        };
        
        ~tdPrimitive()
        {
        };

        void resize(label newsz)
        {
            new_.resize(newsz);
        }
        
        void setTimeStepVariable(bool adjustTimeStep)
        {
            varTimeStep_ = adjustTimeStep;
        }

        label size() const
        {
            return new_.size();
        }

        const T& value(label idx) const
        {
            return new_[idx];
        }

        T& value(label idx)
        {
            return new_[idx];
        }

        T dot(scalar newt, label idx)
        {
            T dvdt (pTraits<T>::zero);

            if (newt > prevt_)
            {
                deltat_ = newt - prevt_;
                prevt_ = newt;
            }

            if (old_.size() <= idx)
            {
                old_.resize(idx+1);
                old_[idx] = new_[idx];
                dvdt = new_[idx] / deltat_;
            }
            else
            {
                if ((oldOld_.size() <= idx) || (varTimeStep_))
                {
                    dvdt = (new_[idx] - old_[idx])/deltat_;
                    oldOld_.resize
                    (
                        idx+1
                    );
                    oldOld_[idx] = old_[idx];
                    old_[idx] = new_[idx];

                }
                else
                {
                    dvdt = (3.0*new_[idx] - 4.0*old_[idx] + oldOld_[idx]) / 2.0 / deltat_;
                    oldOld_[idx] = old_[idx];
                    old_[idx] = new_[idx];
                }
            }
            return dvdt;
        }

    };

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
        //- Get observers list
        List<SoundObserver> getSoundObservers();

        //- Read the Curle data
        virtual bool read(const dictionary&);

        //-
        virtual bool write();

        //-
        virtual bool execute();
};

} //namespace functionObjects

} //namespace Foam

#endif

//
//END-OF-FILE
//



