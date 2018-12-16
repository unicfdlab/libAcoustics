#include "Interpolator.H"

//constructor & destructor

namespace Foam
{

Interpolator::Interpolator()
//:
    //name_(name),
    //pData_(null)
{
}

Interpolator::~Interpolator()
{
}

// private functions

scalar Interpolator::linearInterpolation(const scalar x, const scalar x1, const scalar y1, const scalar x2, const scalar y2)
{
    scalar a = (y2 - y1)/(x2 - x1);
    scalar b = y1 - a*x1;

    scalar y = a*x + b;

    return y;
}

//public functions

autoPtr<List<scalar> > Interpolator::interpolate(const List<scalar> &tSampled, const List<scalar> &data, const label &nTimePoints, const scalar &deltaIntT)
{
    List<scalar> dataInterpolated(nTimePoints,0.0);

    // interpolate pressure data
    scalar cTimeInt = deltaIntT;
    label tI = 1;

    dataInterpolated[0] = data[0];

    for (label j = 1; j < nTimePoints; ++j)
    {
        while (tSampled[tI] < cTimeInt)
        {
            tI++;
            continue;
        }

        dataInterpolated[j] = linearInterpolation(cTimeInt,tSampled[tI],data[tI],tSampled[tI-1],data[tI-1]);
            
        cTimeInt += deltaIntT;
    }

    return autoPtr<List<scalar> >
    (
        new List<scalar>
        (
            dataInterpolated
        )
    );
}

} // end namespace Foam




//END_OF_FILE

