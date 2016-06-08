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

#include "Curle.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::Curle::sampleOrInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vField, 
    const sampledSurface& surface
) const
{
    // interpolator for this field
    autoPtr<interpolation<Type> > interpolatorPtr;
   
    tmp<Field<Type> > values;

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
        values = surface.sample(vField);
    }

    return values;
}

// ************************************************************************* //
