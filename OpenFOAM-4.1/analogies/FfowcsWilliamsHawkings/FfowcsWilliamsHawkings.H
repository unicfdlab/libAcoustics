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
    Foam::FfowcsWilliamsHawkings

Group
    grpFfowcsWilliamsHawkingsFunctionObjects

Description
    

SeeAlso
    

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef FfowcsWilliamsHawkings_H
#define FfowcsWilliamsHawkings_H

#include "AcousticAnalogy.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class FfowcsWilliamsHawkings Declaration
\*---------------------------------------------------------------------------*/

class FfowcsWilliamsHawkings
:
    public AcousticAnalogy
{

        //- Center of surface for sampling
        vector c_;
        
        //- Old value of patch force
        autoPtr<vector> FOldPtr_;
        
        //- Old-old value of patch force
        autoPtr<vector> FOldOldPtr_;
        

protected:

        //- Disallow default bitwise copy construct
        FfowcsWilliamsHawkings(const FfowcsWilliamsHawkings&);

        //- Disallow default bitwise assignment
        void operator=(const FfowcsWilliamsHawkings&);
        
        //-
        void calcDistances();
        
        //-
        void correct();


public:

    //- Runtime type information
    TypeName("FfowcsWilliamsHawkings");


    // Constructors

        //- Construct from given Time and dictionary.
        FfowcsWilliamsHawkings
        (
            const word& name,
            const Time&,
            const dictionary&
        );
        
        //- Construct from given objectRegistry and dictionary.
        FfowcsWilliamsHawkings
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~FfowcsWilliamsHawkings();


    // Member Functions

        //- Read the FfowcsWilliamsHawkings data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool execute();
        
        //- Write the FfowcsWilliamsHawkings
        virtual bool write();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
