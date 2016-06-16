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
    
    Calculatates acoustic pressure using FWH analogy
    
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
#include "SoundObserver.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"

//control surfaces stuff
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "surfaceFieldsFwd.H"
#include "wordReList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class FWH Declaration
\*---------------------------------------------------------------------------*/

class FWH
{
public:

protected:
  
  // Private data
  
  //- Name of this set of FWH,
  //  Also used as the name of the probes directory.
  word name_;

  const objectRegistry& obr_;
  
  //- on/off switch
  bool active_;
  
        //- probe frequency
  label probeFreq_;
  
  //- Switch to send output to Info as well as to file
  Switch log_;
  
  // Read from dictionary
  
  //- Interpolation scheme to use
  word interpolationScheme_;
  
  PtrList<sampledSurface> controlSurfaces_;
  
  //-
  scalar timeStart_;
  
  //-
  scalar timeEnd_;
  
  //- Name of pressure field
  word pName_;
  
  //- Pressure value at infinity
  scalar pInf_;

  vector Ufwh_;
  
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
  
        //-
  vector c_;

  //- FWH/moment file ptr
  autoPtr<OFstream> FWHFilePtr_;
  
  //- Old value of vector integrand
        autoPtr<vector> VOldPtr_;
  
  //- Old-old value of vector integrand
  autoPtr<vector> VOldOldPtr_;
        
  //- Old value of scalar integrand
  autoPtr<scalar> SOldPtr_;
  
  //- Old-old value of scalar integrand
  autoPtr<scalar> SOldOldPtr_;
  
  
  //-
  label probeI_;
  
  // Private Member Functions

  //- If the FWH file has not been created create it
  void makeFile();

  //- Output file header information
  virtual void writeFileHeader();
  
  //-
  void calcDistances();
  
  //derivative function
  Foam::vector dotProduct(const vector& F);

  //derivative function
  Foam::scalar dotProduct(const scalar& S);
        
  //-
  void correct();
  
  //-
  void writeFft();
  
  //- Returns Normal Stresses using pressure (and optionally density) field
  tmp<scalarField> normalStress(const sampledSurface& surface) const;
  //Same for Velocity
  tmp<vectorField> surfaceVelocity(const sampledSurface& surface) const;
  //Same for Density
  tmp<scalarField> surfaceDensity(const sampledSurface& surface) const;

  //- Returns normal stresses time derivative using pressure (and optionally density) field
  tmp<scalarField> dotNormalStress(const sampledSurface& surface) const;
  
  //- Disallow default bitwise copy construct
  FWH(const FWH&);
  
  //- Disallow default bitwise assignment
  void operator=(const FWH&);
  
  //- Sample and write a particular volume field
  template<class Type>
  tmp<Field<Type> > sampleOrInterpolate
  (
   const GeometricField<Type, fvPatchField, volMesh>&,
   const sampledSurface& surface
   ) const;


public:

    //- Runtime type information
    TypeName("FfowcsWilliams-Hawkings");

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

        //- Mark the control surfaces as needing an update.
        //  The function is copy-pasted from sampledSurfaces.C
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        virtual bool needsUpdate() const; 

        virtual bool update();        

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the FWH 
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {};

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {};
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FWHTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
