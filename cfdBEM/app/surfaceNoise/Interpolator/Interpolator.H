#ifndef Interpolator_H
#define Interpolator_H

#include "List.H"
#include "scalar.H"
#include "label.H"
#include "word.H"
#include "autoPtr.H"


namespace Foam
{

class Interpolator
{

private:

    //- linear interpolation
    scalar linearInterpolation(const scalar x, const scalar x1, const scalar y1, const scalar x2, const scalar y2);

public:

    //- constructor
    Interpolator();
    
    //- destructor
    ~Interpolator();

    //- interpolate
    autoPtr<List<scalar> > interpolate(const List<scalar> &tSampled, const List<scalar> &data, const label &nTimePoints, const scalar &deltaIntT);



}; // end of fileInterface

}; // end of namespace Foam

#endif

//END_OF_FILE

