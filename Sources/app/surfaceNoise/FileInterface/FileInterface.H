#ifndef FileInterface_H
#define FileInterface_H


#include "List.H"
#include "scalar.H"
#include "label.H"
#include "word.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "complex.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "Pstream.H"

namespace Foam
{

class FileInterface
{

private:
    
    //- name of input file
    fileName name_;

    //- output format
    word outputFormat_;

    //- write complex numbers in the form Re(z) + Im(z) j
    void writeComplexNumber(const complex& number, autoPtr<OFstream>& os);


public:

    //- constructor
    FileInterface(const fileName& name);

    //- read file
    autoPtr<List<List<scalar> > > read();

    //- write data to file
    void write(const scalar frequency, const List<complex >& data);
    
    //- destructor
    ~FileInterface();

}; // end of FileInterface

}; // end of namespace Foam

#endif

//END_OF_FILE

