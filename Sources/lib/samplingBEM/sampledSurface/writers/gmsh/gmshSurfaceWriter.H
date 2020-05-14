/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::gmshSurfaceWriter

Description
    A surfaceWriter for gmsh legacy format.

SourceFiles
    gmshSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef gmshSurfaceWriter_H
#define gmshSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                      Class gmshSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class gmshSurfaceWriter
:
    public surfaceWriter
{
    // Private Member Functions

        static void writeGeometry(Ostream&, const pointField&, const faceList&);
        
        static void writeData(Ostream& os, const Field<label>&);
        static void writeData(Ostream& os, const Field<scalar>&);
        static void writeData(Ostream& os, const Field<vector>&);
        static void writeData(Ostream& os, const Field<sphericalTensor>&);
        static void writeData(Ostream& os, const Field<symmTensor>&);
        static void writeData(Ostream& os, const Field<tensor>&);

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const Field<Type>& values
        ) const;

public:

    //- Runtime type information
    TypeName("gmsh");


    // Constructors

        //- Construct null
        gmshSurfaceWriter();
        
        //- Construct from dictionary
        gmshSurfaceWriter(const dictionary& );


        //- Destructor
        virtual ~gmshSurfaceWriter();


    // Member Functions
        
        //- Write surface geometry to file.
        virtual fileName write(); // override
                
        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
