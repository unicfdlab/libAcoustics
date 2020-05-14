#include "FoamFftwDriver.H"
#include <fftw3.h>

Foam::FoamFftwDriver::FoamFftwDriver(const List<scalar>& inValues, scalar Tau)
:
    in_(inValues),
    Tau_(Tau)
{
}

Foam::autoPtr<Foam::List<Foam::complex> > Foam::FoamFftwDriver::simpleComplexForwardTransform() const
{
    if (in_.size() <= 0)
    {
        Info << "Input array not allocated!" << endl;
    
        return autoPtr<List<complex > >
        (
            nullptr
        );
    }
    
    label N = in_.size();
    
    fftw_complex* in = (fftw_complex*) fftw_malloc (N*sizeof(fftw_complex));
    fftw_complex* out= (fftw_complex*) fftw_malloc (N*sizeof(fftw_complex));
    
    forAll(in_, k)
    {
        in[k][0] = in_[k];
        in[k][1] = 0.0;
    }
    
    fftw_plan p;
    
    p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    
    fftw_execute(p);
    
    autoPtr<List<complex> > out_resPtr;
    
    out_resPtr.set
    (
        new List<complex>(N, pTraits<complex>::zero)
    );
    
    List<complex>& out_res = out_resPtr();
    
    scalar normCoeff = 2.0/scalar(N); //2.0/Foam::sqrt(scalar(N));
    
    forAll(out_res, k)
    {
        out_res[k].Re() = out[k][0] * normCoeff;    // Re(u)
        out_res[k].Im() = out[k][1] * normCoeff;    // Im(u)
    }
    
    out_res[0].Re() *= 0.5;    // norm coeff for 0 freq = 1/N
    out_res[0].Im() *= 0.5;    // norm coeff for 0 freq = 1/N 
    
    fftw_destroy_plan(p);
    
    if (in)
    {
        fftw_free(in);
        in = 0;
    }
    
    if (out)
    {
        fftw_free(out);
        out = 0;
    }
    
    return out_resPtr;
}

Foam::autoPtr<Foam::Pair<Foam::List<Foam::scalar> > > Foam::FoamFftwDriver::simpleScalarForwardTransform() const
{
    if (in_.size() <= 0)
    {
        Info << "Input array not allocated!" << endl;
        return autoPtr<Pair<List<scalar> > >
        (
            nullptr
        );
    }
    
    label N = in_.size();

    fftw_complex* in = (fftw_complex*) fftw_malloc (N*sizeof(fftw_complex));
    fftw_complex* out= (fftw_complex*) fftw_malloc (N*sizeof(fftw_complex));
    forAll(in_, k)
    {
        in[k][0] = in_[k];
        in[k][1] = 0.0;
    }
    
    fftw_plan p;
    
    p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    
    fftw_execute(p);
    
    Pair<List<scalar> > out_res;
    out_res.first().resize(N);
    out_res.second().resize(N);
    forAll (out_res.first(), k)
    {
        out_res.first()[k] = k / Tau_;
        out_res.second()[k] = 
                        2*sqrt
                        (
                                (out[k][0])*(out[k][0])
                                +
                                (out[k][1])*(out[k][1])
                        )/N;
    }
    
    fftw_destroy_plan(p);

    if (in)
    {
        fftw_free(in);
        in = 0;
    }

    if (out)
    {
        fftw_free(out);
        out = 0;
    }
    
    return autoPtr<Pair<List<scalar> > >
    (
        new Pair<List<scalar> > (out_res)
    );
}

Foam::FoamFftwDriver::~FoamFftwDriver()
{

}

//END_OF_FILE

