#include "FoamFftwDriver.H"

Foam::FoamFftwDriver::FoamFftwDriver(const List<scalar>& inValues, scalar Tau)
:
    in_(inValues),
    Tau_(Tau)
{
}

Foam::autoPtr<Foam::List<Foam::List<Foam::scalar> > > Foam::FoamFftwDriver::simpleForwardTransform() const
{
    if (in_.size() <= 0)
    {
        Info << "Input array not allocated!" << endl;
        
        return autoPtr<List<List<scalar> > >
        (
            0
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

    
    List<List<scalar> > out_res;

    out_res.resize(2);

    forAll(out_res,i)
    {
        out_res[i].resize(N);
    }
    
    scalar normCoeff = 2.0/scalar(N); //2.0/Foam::sqrt(scalar(N));
    
    forAll(out_res[0], k)
    {
        out_res[0][k] = out[k][0] * normCoeff;    // Re(u)
        out_res[1][k] = out[k][1] * normCoeff;    // Im(u)
    }

    out_res[0][0] *= 0.5;    // norm coeff for 0 freq = 1/N
    out_res[1][0] *= 0.5;    // norm coeff for 0 freq = 1/N 
    
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
    
    return autoPtr<List<List<scalar> > >
    (
        new List<List<scalar> > (out_res)
    );
}

Foam::FoamFftwDriver::~FoamFftwDriver()
{

}

//END_OF_FILE

