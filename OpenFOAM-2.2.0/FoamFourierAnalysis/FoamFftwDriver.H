#ifndef FoamFftwDirver_H
#define FoamFftwDriver_H
#include "fftw3.h"
#include "List.H"
#include "scalar.H"
#include "Pair.H"
#include "autoPtr.H"

namespace Foam
{

class FoamFftwDriver
{

private:
    
    //-
    List<scalar> in_;
    
    //-
    scalar Tau_;

public:

    //-
    FoamFftwDriver(const List<scalar>& values, scalar Tau);

    //-
    autoPtr<Pair<List<scalar> > > simpleForwardTransform() const;
    
    //-
    ~FoamFftwDriver();
};

};

#endif
//END_OF_FILE

