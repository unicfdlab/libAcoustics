/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::Curle

Group
    grpCurleFunctionObjects

Description


SeeAlso


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef CurleAnalogy_H
#define CurleAnalogy_H

#include "AcousticAnalogy.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class CurleAnalogy Declaration
\*---------------------------------------------------------------------------*/

class CurleAnalogy
:
    public AcousticAnalogy
{

        //- Center of surface for sampling
        vector c_;

        //- Force acting on patch
        tdPrimitive<vector> F_;

protected:

        //- Disallow default bitwise copy construct
        CurleAnalogy(const CurleAnalogy&);

        //- Disallow default bitwise assignment
        void operator=(const CurleAnalogy&);

        //-
        void calcDistances();

        //-
        void correct();


public:

    //- Runtime type information
    TypeName("CurleAnalogy");


    // Constructors

        //- Construct from given Time and dictionary.
        CurleAnalogy
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        CurleAnalogy
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~CurleAnalogy();


    // Member Functions

        //- Read the CurleAnalogy data
        virtual bool read(const dictionary&);

        //-
        virtual bool execute();

        //- Write the CurleAnalogy
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
