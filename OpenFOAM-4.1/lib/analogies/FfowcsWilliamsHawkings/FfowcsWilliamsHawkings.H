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
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "surfaceFieldsFwd.H"
#include "sampledPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

class fwhFormulation;
class Farassat1AFormulation;


/*---------------------------------------------------------------------------*\
                           Class FfowcsWilliamsHawkings Declaration
\*---------------------------------------------------------------------------*/

class FfowcsWilliamsHawkings
:
    public AcousticAnalogy
{
//friends declation

private:

    friend class fwhFormulation;
    friend class Farassat1AFormulation;
    friend class fvmInterpolationFormulation;

protected:

    friend class fwhFormulation;
    friend class Farassat1AFormulation;
    friend class fvmInterpolationFormulation;

private:

    //- Class used for surface merging information
    class mergeInfo
    {
    public:

        pointField points;
        faceList   faces;
        labelList  pointsMap;
    
        //- Clear all storage
        void clear()
        {
            points.clear();
            faces.clear();
            pointsMap.clear();
        }
    };

    //- Tolerance for merging points (fraction of mesh bounding box)
    static scalar mergeTol_;
    
    //types for internal usage
    
    typedef Pair<List<scalar> > pointTimeData;
    
    typedef List<List<List<pointTimeData > > > surfaceTimeData;


private:

        //-
        word formulationType_;
        
        //-
        Switch fixedResponseDelay_;
        
        //-
        scalar responseDelay_;

        //-
        vector Ufwh_;
        
        //-
        bool nonUniformSurfaceMotion_;
        
        //-
        List<vectorField> Cf0_;
        
        //-
        List<vectorField> vS_;
        
        //- undisturbed pressure, in Pa
        scalar pInf_;
        
        //-
        autoPtr<fwhFormulation> fwhFormulationPtr_;
        
        //
        label cleanFreq_;

        //- Interpolation scheme to use for sampling
        word interpolationScheme_;
        
        //- surfaces for sampling
        PtrList<sampledSurface> controlSurfaces_;

        //- Information for merging surfaces
        List<mergeInfo> mergeList_;

protected:
    
        class sampledPatchAccess : public sampledPatch
        {
        public:
        
            //-
            sampledPatchAccess(const sampledPatch& sp)
            :
                sampledPatch
                (
                    sp
                )
            {
            };
            
            ~sampledPatchAccess()
            {
            };
        
            //-
            const labelList& aPatchIDs() const
            {
                return patchIDs();
            };
            
            //-
            const labelList& aPatchStart() const
            {
                return patchStart();
            };

            //-
            const labelList& aPatchFaceLabels() const
            {
                return patchFaceLabels();
            };
        };

protected:

        //- Disallow default bitwise copy construct
        FfowcsWilliamsHawkings(const FfowcsWilliamsHawkings&);

        //- Disallow default bitwise assignment
        void operator=(const FfowcsWilliamsHawkings&);
        
        //-
        void correct();

        //- Sample and write a particular volume field
        template<class Type>
        tmp<Field<Type> > sampleOrInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const sampledSurface& surface
        ) const;

        //- Returns sampled density field
        tmp<scalarField> surfaceDensity(const sampledSurface& surface) const;

        //- Returns sampled velocity field
        tmp<vectorField> surfaceVelocity(const sampledSurface& surface) const;

        //- Returns sampled pressure field
        tmp<scalarField> surfacePressure(const sampledSurface& surface) const;
        
        //
        tmp<vectorField> surfaceStressDivergence(const sampledSurface& surface) const;
        
        //-
        void initialize();
        
        //- Not needed at current time
        //bool signalReachedObserver(const surfaceTimeData& data, label iObs);


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

        //- Mark the control surfaces as needing an update.
        //  The function is copy-pasted from sampledSurfaces.C
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();
        
        //-
        virtual bool needsUpdate() const; 

        //-
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FWHTemplates.C"

#endif

// ************************************************************************* //
