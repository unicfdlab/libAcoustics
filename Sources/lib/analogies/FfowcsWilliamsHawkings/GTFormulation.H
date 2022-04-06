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
    List<List<tdPrimitive<vector> > > Q_;

    //-
    List<List<tdPrimitive<vector> > > L_;

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
    scalar observerAcousticPressure (const vectorField& Sf, const vectorField& uS, const scalarField& rhoS, const scalarField& pS, label iObs, label iSurf, scalar ct);
    
    //-
    inline void calculateAcousticPressure (const vectorField& Sf, const vectorField& uS, const scalarField& rhoS, const scalarField& pS, label iObs, label iSurf, scalar ct);
    
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

// ************************************************************************* //

//
//END OF FILE
//