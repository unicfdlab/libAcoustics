#ifndef fwhFormulation_H
#define fwhFormulation_H

#include "label.H"
#include "scalar.H"
#include "DynamicList.H"
#include "runTimeSelectionTables.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class FfowcsWilliamsHawkings;

class fwhFormulation
{

private:

    //types for internal usage
    typedef Pair<DynamicList<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;
    
    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

protected:
    
    //-
    const FfowcsWilliamsHawkings& fwh_;
    
    //-
    label fwhProbeI_;

    //-
    surfaceTimeData qds_;
    
    //-
    surfaceTimeData fds_;
    
    //- Keys for identifying observer time
    List<List<List<scalar> > > tobs_;
    
    //- Distances from observers to surfaces elements
    List<List<List<vector> > > robs_;
    
    //- Mag of distances from observers to surfaces elements
    List<List<List<scalar> > > magrobs_;
    
    //-
    List<tdPrimitive<vector> > ni_;
    
    //-
    List<List<scalar> > nl_;

    //-
    List<scalar> rMax_;
    
    //-
    List<scalar> tauMax_;
    
    //-
    List<scalar> rMin_;
    
    //-
    List<scalar> tauMin_;
    

protected:
    
    //-
    void initialize();



    
private:
    //- forbid default constructor
    fwhFormulation();

    //- forbid copy constructor
    fwhFormulation(const fwhFormulation&);

public:
    
    //- constructs using given analogy object
    fwhFormulation(const FfowcsWilliamsHawkings& analogy);
    
    //-
    virtual ~fwhFormulation();
    
    //-
    scalar valueAt (const surfaceTimeData& data, label iObs, label iSurf, label iFace, scalar tau);
    
    //-
    label findExpiredIndex (const pointTimeData& timeData, scalar expiredTime);
    
    //-
    pointTimeData getNewPointData (const pointTimeData& timeData, label expiredIndex);
    
    //-
    surfaceTimeData getQdsData();
    
    //-
    surfaceTimeData getFdsData();
    
    //-
    virtual  scalar observerAcousticPressure (const vectorField& Sf, const vectorField& uS, const scalarField& rhoS, const scalarField& pS, label iObs, label iSurf, scalar ct);
    
    //-
    virtual void update();
    
    //-
    virtual void clearExpiredData();
};

}; //namespace functionObjects

}; // namespace Foam

#endif

//
//END-OF-FILE
//


