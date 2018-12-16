#include "FileInterface.H"

//constructor & destructor

namespace Foam
{

FileInterface::FileInterface(const fileName& name)
:
    name_(name),
    outputFormat_("gmsh")
{
}

FileInterface::~FileInterface()
{
}

//private functions

void FileInterface::writeComplexNumber(const complex& number, autoPtr<OFstream>& os)
{
    os()    << number.Re();

    label sgn = sign(number.Im());

    if (sgn == 1)
    {
        os() << "+";
    }
    else
    {
        os() << "-";
    }
    
    os() << fabs(number.Im());
    os() << "j";
}

//public functions

autoPtr<List<List<scalar > > > FileInterface::read()
{
    autoPtr<IFstream> inFilePtr;

    if (inFilePtr.empty())
    {
        inFilePtr.reset
        (
            new IFstream
            (
                name_
            )
        );
    }

    List<List<scalar> > data;

    label probeI = 0;

    scalar num = 0;

    string line;

    if (inFilePtr.valid())
    {
        //Info << inFilePtr().name() << endl;

        while(inFilePtr().getLine(line))
        {
            probeI ++;

            data.resize(probeI);

            //Info << "resize ok" << endl;

            label pointI = 0;

            IStringStream lineStream(line);

            while(lineStream >> num)
            {  
                pointI++;

                data[probeI-1].resize(pointI);

                data[probeI-1][pointI-1] = num;
            }

            //inFilePtr() >> num ;
            //Info << num << nl;

            //Info << "///////////////////////////////" << nl;
            //istr++;
        }

        //check if last string has equal size with previous strings
        label refSize = data[0].size();
        label lastStrSize = data[probeI-1].size();

        //Info << "ref size = " << refSize << ", size of last string = " << lastStrSize << endl;

        if (lastStrSize < refSize) // short last string shows the lost part of data
        {
            Info << "Cut tail of damaged data"  << endl;
            data.resize(probeI-1);
        }

        Info << "OK" << nl;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::autoPtr<Foam::List<Foam::List<scalar > > > read()"
        )   << "IFstream is not valid " << nl
            << exit(FatalError);
    }

    return autoPtr<List<List<scalar> > >
    (
        new List<List<scalar> >
        (
            data
        )
    );

}

void FileInterface::write(const scalar frequency, const List<complex >& data)
{

    if (Pstream::master() || !Pstream::parRun())
    {
        if (name_.ext() == "msh")
        {
            autoPtr<OFstream> filePtr;

            if (filePtr.empty())
            {
            // Open new file at start up
                filePtr.reset
                (
                    new OFstream
                    (
                       name_
                   )
               );
            }

            filePtr()   << "$NodeData" << nl;

            //string tags
            filePtr()   << "1" << nl;           // number of string tags
            filePtr()   << "node_data" << nl;   // name of data

            //real tags
            filePtr()   << "1" << nl            // number of real tags
                        << frequency << nl;     // should be output time (or frequency)
            
            //integer tags
            filePtr()   << "4" << nl            // number of int tags
                        << "0" << nl            // time step index
                        << "1" << nl            // how many field components
                        << data.size() << nl    // number of entities (nodes or elements)
                        << "0" << nl;

            forAll(data, probeI)
            {
                filePtr() << probeI + 1 << ' ';
                writeComplexNumber(data[probeI],filePtr);
                filePtr() << nl;       
            }

            filePtr() << "$EndNodeData" << nl;
        }
        else
        {
            FatalErrorIn
            (
                "Foam::FileInterface::write(const List<complex >& data)"
            )   << "Unsupported output file format: " <<  name_.ext() << nl
                << "Use file extension .msh" << nl
                << exit(FatalError);
        }
    }
}

} // end namespace Foam






//END_OF_FILE

