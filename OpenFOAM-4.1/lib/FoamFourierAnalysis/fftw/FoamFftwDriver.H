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
    
    //- list of input data
    List<scalar> in_;
    
    //- full time interval
    scalar Tau_;

public:

    //- constructor
    FoamFftwDriver(const List<scalar>& values, scalar Tau);

    //- return result of fft: freq/Re/Im
    autoPtr<List<List<scalar> > > simpleForwardTransform() const;
    
    //- destructor
    ~FoamFftwDriver();
};

};

#endif
//END_OF_FILE

