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

\*---------------------------------------------------------------------------*/

#include "soundPressureSampler.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

template<class Type>
void Foam::soundPressureSampler::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField,
    const label surfI,
    autoPtr <OFstream>& os,
    scalar cTime
) const
{

    // interpolator for this field
    autoPtr<interpolation<Type> > interpolatorPtr;
    Field<Type> values;

    const sampledSurface& surface = controlSurfaces_.operator[](surfI);


    
    if (surface.interpolate())
    {
        if (interpolatorPtr.empty())
        {
            interpolatorPtr = interpolation<Type>::New
            (
                interpolationScheme_,
                vField
            );
        }
        values = surface.interpolate(interpolatorPtr());
    }
    else
    {
        if (interpolatorPtr.empty())
        {
            interpolatorPtr = interpolation<Type>::New
            (
                "cell",
                vField
            );
        }
        values = surface.sample(interpolatorPtr());
    }

    writeSurface(values, surfI, os, cTime);
}


template<class Type>
void Foam::soundPressureSampler::writeSurface
(
    const Field<Type>& values,
    const label surfI,
    autoPtr <OFstream>& os,
    scalar cTime
) const
{
    if (Pstream::parRun())
    {
        List<Field<Type> > gatheredValues(Pstream::nProcs());
        
        gatheredValues[Pstream::myProcNo()] = values;
        
        Pstream::gatherList(gatheredValues);

        if (Pstream::master())
        {
            // Combine values into single field
            Field<Type> allValues
            (
                ListListOps::combine<Field<Type> >
                (
                    gatheredValues,
                    accessOp<Field<Type> >()
                )
            );

            // Renumber (point data) to correspond to merged points
            if (mergeList_[surfI].pointsMap.size() == allValues.size())
            {
                inplaceReorder(mergeList_[surfI].pointsMap, allValues);
                allValues.setSize(mergeList_[surfI].points.size());
            }
            
            // Write to defined directory
            if (mergeList_[surfI].faces.size())
            {
                if (os.valid())
                {
                    os() << cTime << ' ';
                    
                    forAll(allValues,i)
                    {                  
                        os() << allValues[i] << " ";
                    }
                    
                    os() << nl;
                    os().flush();
                }
            }
        }
    }
    else
    {
        Info << "size = " << values.size() << nl;

        if (os.valid())
        {
            os() << cTime << ' ';
            
            forAll(values,i)
            {
                os() << values[i] << ' ';
            }
            os() << nl;
            os().flush();
        }
    }
}
