#ifndef GTFormulation_H
#define GTFormulation_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class GTFormulation : public fwhFormulation
{

private:

    //types for internal usage
    
    typedef Pair<DynamicList<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    
    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

private:

    //-
    List<List<tdPrimitive<scalar> > > Qn_;
    
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
    GTFormulation();
    
    //- Disable copy constructor
    GTFormulation(const GTFormulation&);

public:
    
    //-
    GTFormulation(const FfowcsWilliamsHawkings& fwh);
    
    //-
    ~GTFormulation();
    
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

