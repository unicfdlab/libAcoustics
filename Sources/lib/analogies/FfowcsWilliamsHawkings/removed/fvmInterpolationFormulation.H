#ifndef fvmInterpolationFormulation_H
#define fvmInterpolationFormulation_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class fvmInterpolationFormulation : public fwhFormulation
{

private:

    //types for internal usage
    
    typedef Pair<DynamicList<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    
    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

private:

    //-
    tdPrimitive<scalar> intQdS_;
    
    //-
    tdPrimitive<scalar> intFdS_;
    
private:
    
    //- Disable default constructor
    fvmInterpolationFormulation();
    
    //- Disable copy constructor
    fvmInterpolationFormulation(const fvmInterpolationFormulation&);

public:
    
    //-
    fvmInterpolationFormulation(const FfowcsWilliamsHawkings& fwh);
    
    //-
    ~fvmInterpolationFormulation();
    
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

