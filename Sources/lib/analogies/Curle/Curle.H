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

#ifndef Curle_H
#define Curle_H

#include "AcousticAnalogy.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class Curle Declaration
\*---------------------------------------------------------------------------*/

class Curle
:
    public AcousticAnalogy
{

        //- Center of surface for sampling
        vector c_;
        
        //- Force acting on patch
        tdPrimitive<vector> F_;
        
        //- Force acting on patch
        tdPrimitive<vector> dF_;
        
        //
        List<vector> FF_;
        
        //
        scalar iter_;
        

protected:

        //- Disallow default bitwise copy construct
        Curle(const Curle&);

        //- Disallow default bitwise assignment
        void operator=(const Curle&);
        
        //-
        void calcDistances();
        
        //-
        void correct();


public:

    //- Runtime type information
    TypeName("Curle");


    // Constructors

        //- Construct from given Time and dictionary.
        Curle
        (
            const word& name,
            const Time&,
            const dictionary&
        );
        
        //- Construct from given objectRegistry and dictionary.
        Curle
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~Curle();


    // Member Functions

        //- Read the Curle data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool execute();
        
        //- Write the Curle
        virtual bool write();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
