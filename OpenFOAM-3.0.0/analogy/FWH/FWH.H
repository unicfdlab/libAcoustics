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
    Foam::FWH

Description
    
    Calculatates acoustic pressure using FW-H analogy
    
    Writes data to specified file


SourceFiles
    FWH.C
    IOFWH.H

\*---------------------------------------------------------------------------*/

#ifndef FWH_H
#define FWH_H

#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"
#include "SoundObserverM.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "volFields.H"
#include "FWHSampledTriSurfaceMesh.H"
#include "writer.H"
#include "noiseFFTDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class FW-H Declaration
\*---------------------------------------------------------------------------*/

class FWH
{
public:

protected:

    // Private data

        //- Name of this set of FW-H,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;
	
	//- on/off switch
        bool active_;
        
        //- probe frequency
        label probeFreq_;

        //- Switch to send output to Info as well as to file
        Switch log_;
        
        bool bigdate_;

        //- fwh cells
//	FWHSampledTriSurfaceMesh fwhSampledSurface_;

	//- Patch to integrate FW-H over

	word	fwhSurfaceType_;

	word	fwhSurfaceName_;

//	word 	patchNames_;
        //-
        scalar timeStart_;
        
        //-
        scalar timeEnd_;

        //- Name of pressure field
        word pName_;
        
        //
        word UName_;
        
        //- Sonic speed, m/s
        scalar c0_;
        
        //- Domain depth, m
        scalar dRef_;
        
        //-
        List<SoundObserver> observers_;

        //- Name of density field 
        word rhoName_;

	//- Reference density needed for incompressible calculations
	scalar rhoRef_;
	
	//pref
	scalar pRef_;

        //-
        vector c_;

        //- FWH/moment file ptr
        autoPtr<OFstream> FWHFilePtr_;
        
        //- Old value of patch force
//        autoPtr<vector> FOldPtr_;
        
        //- Old-old value of patch force
//        autoPtr<vector> FOldOldPtr_;
        
        //- Old value of patch force
//        autoPtr<scalar> FFOldPtr_;
	List<autoPtr<scalar> > FFOldPtr_;

        
        //- Old-old value of patch force
        List<autoPtr<scalar> >FFOldOldPtr_;

//	List<scalar> FFOldOld_;
        //-
        label probeI_;
        
        //-
        word graphFormat_;


    // Private Member Functions

        //- If the FWH file has not been created create it
        void makeFile();

        //- Output file header information
        virtual void writeFileHeader();

        //-
        
        void correctFaceSetFWH();
        
        void correctSampledFWH();
        //-
        
        void correct();
        
        //-
        void writeFft();
//
//        //- Return the effective viscous stress (laminar + turbulent).
//        tmp<volSymmTensorField> devRhoReff() const;
//
//        //- Return rho if rhoName is specified otherwise rhoRef
//        tmp<volScalarField> rho(const volScalarField& p) const;

	//- Returns normal stresses using pressure (and optionally density) field
//	tmp<scalarField> normalStress(const word& patchName) const;

        //- Disallow default bitwise copy construct
        FWH(const FWH&);

        //- Disallow default bitwise assignment
        void operator=(const FWH&);


public:

    //- Runtime type information
    TypeName("FW-H");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        FWH
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

	//- Destructor
	virtual ~FWH();


    // Member Functions

        //- Return name of the set of FWH
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the FWH data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();
        
        //-for 2.3
        virtual void timeSet();

        //- Write the FWH
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const Foam::mapPolyMesh&)
        {};

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {};
        

    //- Update for changes of mesh
    //    virtual void updateMesh(const mapPolyMesh&)
    //    {};

        //- Update for changes of mesh
    //    virtual void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
