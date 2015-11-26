#include "FoamFftwDriver.H"

Foam::FoamFftwDriver::FoamFftwDriver(const List<scalar>& inValues, scalar Tau)
:
    in_(inValues),
    Tau_(Tau)
{
}

Foam::autoPtr<Foam::Pair<Foam::List<Foam::scalar> > > Foam::FoamFftwDriver::simpleForwardTransform() const
{
    if (in_.size() <= 0)
    {
	Info << "Input array not allocated!" << endl;
	return autoPtr<Pair<List<scalar> > >
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
    
    Pair<List<scalar> > out_res;
    out_res.first().resize(N);
    out_res.second().resize(N);
    forAll (out_res.first(), k)
    {
	out_res.first()[k] = k / Tau_;
	out_res.second()[k] = 
			2*sqrt
			(
				(out[k][0]/N)*(out[k][0]/N)
				+
				(out[k][1]/N)*(out[k][1]/N)
			);
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

