/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gmshSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
//    makeSurfaceWriterType(gmshSurfaceWriter);

    defineTypeName(gmshSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, gmshSurfaceWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, gmshSurfaceWriter, wordDict);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceWriters::gmshSurfaceWriter::writeGeometry
(
    Ostream& os,
    const pointField& points,
    const faceList& faces
)
{
    // Header
    os
        << "$MeshFormat" << nl
        << "2.2 0 " << int(sizeof(double)) << nl
        << "$EndMeshFormat" << nl;

    // Write vertex coords
    os  
        << "$Nodes " << nl
        << points.size() << nl;
    forAll(points, pointI)
    {
        const point& pt = points[pointI];
        os  << pointI + 1 << ' '
            << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
    os  << "$EndNodes" << nl;


    // Write faces
    label nNodes = 0;
    forAll(faces, faceI)
    {
        nNodes += faces[faceI].size();
    }

    os  
        << "$Elements" << nl
        << faces.size() << nl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        os  << faceI + 1 << " 2 2 0 0 ";

        forAll(f, fp)
        {
            os << f[fp] + 1 << ' ';
        }
        os  << nl;
    }

    os  << "$EndElements" << nl;
}

void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<scalar>& values
)
{
    os  
        << "1" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;

    forAll(values, elemI)
    {
        os  << elemI + 1 << ' ' 
            << float(values[elemI]) 
            << nl;
    }
}

void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<label>& values
)
{
    os  
        << "1" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;

    forAll(values, elemI)
    {
        os  << elemI + 1 << ' ' 
            << label(values[elemI]) 
            << nl;
    }
}

void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<vector>& values
)
{
    os  
        << "3" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;

    forAll(values, elemI)
    {
        const vector& v = values[elemI];
        os  << elemI + 1 << ' '
            << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
            << nl;
    }
}


void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<sphericalTensor>& values
)
{
    os  
        << "1" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;

    forAll(values, elemI)
    {
        const sphericalTensor& v = values[elemI];
        os  << elemI + 1 << ' '
            << float(v[0]) 
            << nl;
    }
}


void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<symmTensor>& values
)
{
    os  
        << "6" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;
   
    forAll(values, elemI)
    {
        const symmTensor& v = values[elemI];
        os  << elemI + 1 << ' '
            << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
            << ' '
            << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
            << nl;
    }
}



void Foam::surfaceWriters::gmshSurfaceWriter::writeData
(
    Ostream& os,
    const Field<tensor>& values
)
{
    os  
        << "9" << nl            // how many field components
        << values.size() << nl // number of entities (nodes or elements)
        << "0" << nl;

    forAll(values, elemI)
    {
        const tensor& v = values[elemI];
        os  << elemI + 1 << ' '
            << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
            << ' '
            << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
            << ' '
            << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
            << nl;
    }
}


template<class Type>
Foam::fileName Foam::surfaceWriters::gmshSurfaceWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& values
) const
{
    const fileName outputDir   = this->outputPath_;
    const fileName surfaceName = this->outputPath_.name();
    const meshedSurf& surf     = this->surf_;
    const bool isNodeValues    = this->isPointData();
    
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".msh");

    if (verbose())
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    // Write geometry
    writeGeometry(os, surf.points(), surf.faces());


    // start writing data
    if (isNodeValues)
    {
        os  << "$NodeData" << nl;
    }
    else
    {
        os  << "$ElementData" << nl;
    }

    //string tags
    os  << "1" << nl;   // how much
    isNodeValues ? os << "node_data" << nl : os << "element_data" << nl; //name of data

    //real tags
    os  << "1" << nl    // how much
        << "0" << nl;    // should be output time
        
    //integer tags
    os  << "4" << nl    // how much
        << "0" << nl;    // time step index

    // Write data
    writeData(os, values);

    // end writing data
    if (isNodeValues)
    {
        os  << "$EndNodeData" << nl;
    }
    else
    {
        os  << "$EndElementData" << nl;
    }

    return os.name();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::gmshSurfaceWriter::gmshSurfaceWriter()
:
    surfaceWriter()
{}

Foam::surfaceWriters::gmshSurfaceWriter::gmshSurfaceWriter(const dictionary& dict)
:
    surfaceWriter(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriters::gmshSurfaceWriter::~gmshSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::gmshSurfaceWriter::write()
{
    Info << "Checking surface" << endl;
    checkOpen();
    
    fileName outputDir = outputPath_;
    fileName outputFile= outputPath_.path() / timeName() / outputPath_.name();
    outputFile.ext("msh");
    
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }
        
        OFstream os
        (
            outputFile
        );

        if (verbose())
        {
            Info<< "Writing geometry to " << os.name() << endl;
        }

        writeGeometry(os, surf.points(), surf.faces());
    }
    
    wroteGeom_ = true;
    return outputFile;
}

// create write methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::gmshSurfaceWriter);


// ************************************************************************* //
