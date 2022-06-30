#ifndef Farassat1AFormulation_H
#define Farassat1AFormulation_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class Farassat1AFormulation : public fwhFormulation
{

private:

    //types for internal usage
    
    typedef Pair<DynamicList<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    
    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

private:

    //-
    List<List<tdPrimitive<scalar> > > Un_;
    
    //-
    List<List<tdPrimitive<scalar> > > Lr_;
    
    //-
    List<List<tdPrimitive<scalar> > > Mr_;

    //-
    tdPrimitive<scalar> intDotQdS_;
    
    //-
    tdPrimitive<scalar> intFdS_;
    
private:
    
    //- Disable default constructor
    Farassat1AFormulation();
    
    //- Disable copy constructor
    Farassat1AFormulation(const Farassat1AFormulation&);

public:
    
    //-
    Farassat1AFormulation(const FfowcsWilliamsHawkings& fwh);
    
    //-
    ~Farassat1AFormulation();
    
    //-
    scalar observerAcousticPressure(label iObs);
    
    //-
    void clearExpiredData();
    
    //-
    void update();
    
    //-
    void initialize();
};

}

} // namespace Foam

#endif

