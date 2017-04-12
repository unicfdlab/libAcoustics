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

Class
    Foam::soundPressureSampler

Description

    Return sound amplitudes of sound pressure in each point of STL-surface.
    Data are sorted by frequencies.

    Input data:
    -   path to triSurface;
    -   reference pressure.

    Output:
    -   file with pressure history.

SourceFiles
    soundPressureSampler.C

\*---------------------------------------------------------------------------*/

#ifndef soundPressureSampler_H
#define soundPressureSampler_H

#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

#include "polyMesh.H"
#include "mapPolyMesh.H"

#include "volFields.H"

//control surfaces stuff (!)
#include "functionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "surfaceFieldsFwd.H"
#include "wordReList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                    Class soundPressureSampler Declaration
\*---------------------------------------------------------------------------*/

class soundPressureSampler
:
    public functionObject
{


private:

//- Class used for surface merging information
    class mergeInfo
    {

    public:

	pointField points;
	faceList faces;
	labelList pointsMap;

	//- Clear all storage
	void clear()
	{
	    points.clear();
	    faces.clear();
    	    pointsMap.clear();
	}
    };
    
    //- Output verbosity
    static bool verbose_;
    
    //- Tolerance for merging points (fraction of mesh bounding box)
    static scalar mergeTol_;
    
protected:

    // Private data
    
	//- Const reference to fvMesh
	const fvMesh& mesh_;
	        
	//- Load fields from files (not from objectRegistry)
	const bool loadFromFiles_;
	                        
        
        //- Name of this set of soundPressureSampler,
        // Also used as the name of the probes directory.
        word name_;
        
        //- Standard OpenFOAM object registry
        //const objectRegistry& obr_;
        
        //- Switch to send output to Info as well as to file
        Switch log_;


        //- Interpolation scheme to use for sampling
        word interpolationScheme_;

        //- List of control surfaces
        PtrList<sampledSurface> controlSurfaces_;

        //- Output format for surface geometry writing
        word outputFormat_;

        //- Writer for surface geometry
        autoPtr<surfaceWriter> formatter_;

        //- Name of pressure field
        word pName_;

        //- Pressure value at infinity 
        //scalar pRef_;
               
        //- Smart pointer to file output stream
        List< autoPtr <OFstream> > soundPressureSamplerFilePtr_;

        //- Information for merging surfaces
        List<mergeInfo> mergeList_;

        //- Number of probe
        label probeI_;

private:    

        //- Disallow default bitwise copy construct
        soundPressureSampler(const soundPressureSampler&);

        //- Disallow default bitwise assignment
        void operator=(const soundPressureSampler&);


        //- Write geometry only
        void writeGeometry() const;


        //- If the soundPressureSampler file has not been created create it
        void makeFile();


        //- Sample a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            const label surfI,
            autoPtr <OFstream>& os,
            scalar cTime
        ) const;

        //- Collect data and write to file
        template<class Type>
        void writeSurface
        (
            const Field<Type>& values,
            const label surfI,
            autoPtr <OFstream>& os,
            scalar cTime
        ) const;

public:
    
    //- Runtime type information
    TypeName("soundPressureSampler");

    // Constructors
    
        //- Construct from Time and dictionary
        soundPressureSampler
        (
	    const word& name,
            const Time& time,
            const dictionary& dict
        );
                                                                    
	//- Construct for given objectRegistry and dictionary
	//  allow the possibility to load fields from files
        soundPressureSampler
	(
    	    const word& name,
    	    const objectRegistry&,
    	    const dictionary&,
    	    const bool loadFromFiles = false
	);	
                                                                                                                                                            

    //- Destructor
    virtual ~soundPressureSampler();

    // Member Functions

        //- Return name of the set of soundPressureSampler
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the soundPressureSampler data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        //virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the soundPressureSampler
        virtual bool write();

    // Copy-pasted functions from libSampling for correct working of surface sampling

        //- Mark the control surfaces as needing an update.
        // The function is copy-pasted from sampledSurfaces.C
        // May also free up unneeded data.
        // Return false if all surfaces were already marked as expired.
        virtual bool expire();

        virtual bool needsUpdate() const;

        virtual bool update();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {};

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {};

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set tolerance (and return old tolerance)
        static scalar mergeTol(const scalar);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "soundPressureSamplerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
