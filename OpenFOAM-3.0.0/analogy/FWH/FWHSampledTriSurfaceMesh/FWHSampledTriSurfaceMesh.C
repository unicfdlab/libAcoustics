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

#include "FWHSampledTriSurfaceMesh.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "meshTools.H"

#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FWHSampledTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        FWHSampledTriSurfaceMesh,
        word
    );

    template<>
    const char* NamedEnum<FWHSampledTriSurfaceMesh::samplingSource, 3>::names[] =
    {
        "cells",
        "insideCells",
        "boundaryFaces"
    };

    const NamedEnum<FWHSampledTriSurfaceMesh::samplingSource, 3>
    FWHSampledTriSurfaceMesh::samplingSourceNames_;


    //- Private class for finding nearest
    //  Comprising:
    //  - global index
    //  - sqr(distance)
    typedef Tuple2<scalar, label> nearInfo;

    class nearestEqOp
    {

    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first() < x.first())
            {
                x = y;
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataFace>&
Foam::FWHSampledTriSurfaceMesh::nonCoupledboundaryTree() const
{
    // Variant of meshSearch::boundaryTree() that only does non-coupled
    // boundary faces.

    if (!boundaryTreePtr_.valid())
    {
        // all non-coupled boundary faces (not just walls)
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        labelList bndFaces(mesh().nFaces()-mesh().nInternalFaces());
        label bndI = 0;
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (!pp.coupled())
            {
                forAll(pp, i)
                {
                    bndFaces[bndI++] = pp.start()+i;
                }
            }
        }
        bndFaces.setSize(bndI);


        treeBoundBox overallBb(mesh().points());
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        boundaryTreePtr_.reset
        (
            new indexedOctree<treeDataFace>
            (
                treeDataFace    // all information needed to search faces
                (
                    false,                      // do not cache bb
                    mesh(),
                    bndFaces                    // boundary faces only
                ),
                overallBb,                      // overall search domain
                8,                              // maxLevel
                10,                             // leafsize
                3.0                             // duplicity
            )
        );
    }

    return boundaryTreePtr_();
}


bool Foam::FWHSampledTriSurfaceMesh::update(const meshSearch& meshSearcher)
{
    // Find the cells the triangles of the surface are in.
    // Does approximation by looking at the face centres only
    const pointField& fc = surface_.faceCentres();

    List<nearInfo> nearest(fc.size());

    // Global numbering for cells/faces - only used to uniquely identify local
    // elements
    globalIndex globalCells
    (
        (sampleSource_ == cells || sampleSource_ == insideCells)
      ? mesh().nCells()
      : mesh().nFaces()
    );

    forAll(nearest, i)
    {
        nearest[i].first() = GREAT;
        nearest[i].second() = labelMax;
    }

    if (sampleSource_ == cells)
    {
        // Search for nearest cell

        const indexedOctree<treeDataCell>& cellTree = meshSearcher.cellTree();

        forAll(fc, triI)
        {
            pointIndexHit nearInfo = cellTree.findNearest
            (
                fc[triI],
                sqr(GREAT)
            );
            if (nearInfo.hit())
            {
                nearest[triI].first() = magSqr(nearInfo.hitPoint()-fc[triI]);
                nearest[triI].second() = globalCells.toGlobal(nearInfo.index());
            }
        }
    }
    else if (sampleSource_ == insideCells)
    {
        // Search for cell containing point

        const indexedOctree<treeDataCell>& cellTree = meshSearcher.cellTree();

        forAll(fc, triI)
        {
            if (cellTree.bb().contains(fc[triI]))
            {
                label index = cellTree.findInside(fc[triI]);
                if (index != -1)
                {
                    nearest[triI].first() = 0.0;
                    nearest[triI].second() = globalCells.toGlobal(index);
                }
            }
        }
    }
    else
    {
        // Search for nearest boundaryFace
        ////- Search on all (including coupled) boundary faces
        //const indexedOctree<treeDataFace>& bTree = meshSearcher.boundaryTree()
        //- Search on all non-coupled boundary faces
        const indexedOctree<treeDataFace>& bTree = nonCoupledboundaryTree();

        forAll(fc, triI)
        {
            pointIndexHit nearInfo = bTree.findNearest
            (
                fc[triI],
                sqr(GREAT)
            );
            if (nearInfo.hit())
            {
                nearest[triI].first() = magSqr(nearInfo.hitPoint()-fc[triI]);
                nearest[triI].second() = globalCells.toGlobal
                (
                    bTree.shapes().faceLabels()[nearInfo.index()]
                );
            }
        }
    }


    // See which processor has the nearest. Mark and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    labelList cellOrFaceLabels(fc.size(), -1);

    label nFound = 0;
    forAll(nearest, triI)
    {
        if (nearest[triI].second() == labelMax)
        {
            // Not found on any processor. How to map?
        }
        else if (globalCells.isLocal(nearest[triI].second()))
        {
            cellOrFaceLabels[triI] = globalCells.toLocal
            (
                nearest[triI].second()
            );
            nFound++;
        }
    }


    if (debug)
    {
        Pout<< "Local out of faces:" << cellOrFaceLabels.size()
            << " keeping:" << nFound << endl;
    }
	
    reduce (cellOrFaceLabels, sumOp<labelList>());

    // Now subset the surface. Do not use triSurface::subsetMesh since requires
    // original surface to be in compact numbering.

    const triSurface& s = surface_;

    // Compact to original triangle
    labelList faceMap(s.size());
    
    //*************DOING THIS FOR faces() in parallel
    // Compact to original points
    labelList pointMap(s.points().size());
    // From original point to compact points
    labelList reversePointMap(s.points().size(), -1);

    {
        label newPointI = 0;
        label newFaceI = 0;

        forAll(s, faceI)
        {
            if (cellOrFaceLabels[faceI] != -1)
            {
                faceMap[newFaceI++] = faceI;

                const triSurface::FaceType& f = s[faceI];
                forAll(f, fp)
                {
                    if (reversePointMap[f[fp]] == -1)
                    {
                        pointMap[newPointI] = f[fp];
                        reversePointMap[f[fp]] = newPointI++;
                    }
                }
            }
        }
        faceMap.setSize(newFaceI);
        pointMap.setSize(newPointI);
    }


    // Subset cellOrFaceLabels
    cellOrFaceLabels = UIndirectList<label>(cellOrFaceLabels, faceMap)();

    // Store any face per point (without using pointFaces())
    labelList pointToFace(pointMap.size());

    // Create faces and points for subsetted surface
    faceList& faces = this->storedFaces();
    faces.setSize(faceMap.size());
    forAll(faceMap, i)
    {
        const triFace& f = s[faceMap[i]];
        triFace newF
        (
            reversePointMap[f[0]],
            reversePointMap[f[1]],
            reversePointMap[f[2]]
        );
        faces[i] = newF.triFaceFace();

        forAll(newF, fp)
        {
            pointToFace[newF[fp]] = i;
        }
    }
    
    //*****************OLD SHORT CODE BELOW
    
 /*   // Compact to original triangle
    labelList faceMap(s.size());
    // Compact to original points
	
    {
        label newFaceI = 0;

        forAll(s, faceI)
        {
            if (cellOrFaceLabels[faceI] != -1)
            {
                faceMap[newFaceI++] = faceI;
            }
        }
        faceMap.setSize(newFaceI);
    }


    // Subset cellOrFaceLabels
    cellOrFaceLabels = UIndirectList<label>(cellOrFaceLabels, faceMap)();

    faceList& faces = this->storedFaces();
    faces.setSize(faceMap.size());
    forAll(faceMap, i)
    {
        faces[i] = s[i];
    }
  */  
    this->storedPoints() = pointField(fc, faceMap);

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
   
    // Collect the samplePoints and sampleElements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (sampledSurface::interpolate())
    {
        samplePoints_.setSize(fc.size());
        sampleElements_.setSize(fc.size(), -1);

        if (sampleSource_ == cells)
        {
            // samplePoints_   : per surface point a location inside the cell
            // sampleElements_ : per surface point the cell

            forAll(fc, pointI)
            {
                const point& pt = fc[pointI];
                label cellI = cellOrFaceLabels[pointI];
                sampleElements_[pointI] = cellI;

                // Check if point inside cell
                if
                (
                    mesh().pointInCell
                    (
                        pt,
                        sampleElements_[pointI],
                        meshSearcher.decompMode()
                    )
                )
                {
                    samplePoints_[pointI] = pt;
                }
                else
                {
                    // Find nearest point on faces of cell
                    const cell& cFaces = mesh().cells()[cellI];

                    scalar minDistSqr = VGREAT;

                    forAll(cFaces, i)
                    {
                        const face& f = mesh().faces()[cFaces[i]];
                        pointHit info = f.nearestPoint(pt, mesh().points());
                        if (info.distance() < minDistSqr)
                        {
                            minDistSqr = info.distance();
                            samplePoints_[pointI] = info.rawPoint();
                        }
                    }
                }
            }
        }
        else if (sampleSource_ == insideCells)
        {
            // samplePoints_   : per surface point a location inside the cell
            // sampleElements_ : per surface point the cell

            forAll(fc, pointI)
            {
                const point& pt = fc[pointI];
                label cellI = cellOrFaceLabels[pointI];
                sampleElements_[pointI] = cellI;
                samplePoints_[pointI] = pt;
            }
        }
        else
        {
            // samplePoints_   : per surface point a location on the boundary
            // sampleElements_ : per surface point the boundary face containing
            //                   the location

            forAll(fc, pointI)
            {
                const point& pt = fc[pointI];
                label faceI = cellOrFaceLabels[pointI];
                sampleElements_[pointI] = faceI;
                samplePoints_[pointI] =  mesh().faces()[faceI].nearestPoint
                (
                    pt,
                    mesh().points()
                ).rawPoint();
            }
        }
    }
    else
    {
        // if sampleSource_ == cells:
        //      samplePoints_   : n/a
        //      sampleElements_ : per surface triangle the cell
        // if sampleSource_ == insideCells:
        //      samplePoints_   : n/a
        //      sampleElements_ : -1 or per surface triangle the cell
        // else:
        //      samplePoints_   : n/a
        //      sampleElements_ : per surface triangle the boundary face
        samplePoints_.clear();
        sampleElements_.transfer(cellOrFaceLabels);
    }


    if (debug)
    {
        this->clearOut();
        OFstream str(mesh().time().path()/"surfaceToMesh.obj");
        Info<< "Dumping correspondence from local surface (points or faces)"
            << " to mesh (cells or faces) to " << str.name() << endl;
        label vertI = 0;

        if (sampledSurface::interpolate())
        {
            if (sampleSource_ == cells || sampleSource_ == insideCells)
            {
                forAll(samplePoints_, pointI)
                {
                    meshTools::writeOBJ(str, samplePoints_[pointI]);
                    vertI++;

                    label cellI = sampleElements_[pointI];
                    meshTools::writeOBJ(str, mesh().cellCentres()[cellI]);
                    vertI++;
                    str << "l " << ' ' << vertI-1 << ' ' << vertI
                        << nl;
                }
            }
            else
            {
                forAll(samplePoints_, pointI)
                {
                    meshTools::writeOBJ(str, samplePoints_[pointI]);
                    vertI++;

                    label faceI = sampleElements_[pointI];
                    meshTools::writeOBJ(str, mesh().faceCentres()[faceI]);
                    vertI++;
                    str << "l " << ' ' << vertI-1 << ' ' << vertI
                        << nl;
                }
            }
        }
        else
        {
            if (sampleSource_ == cells || sampleSource_ == insideCells)
            {
                forAll(sampleElements_, triI)
                {
                    meshTools::writeOBJ(str, faceCentres()[triI]);
                    vertI++;

                    label cellI = sampleElements_[triI];
                    meshTools::writeOBJ(str, mesh().cellCentres()[cellI]);
                    vertI++;
                    str << "l " << ' ' << vertI << nl;
                }
            }
            else
            {
                forAll(sampleElements_, triI)
                {
                    meshTools::writeOBJ(str, faceCentres()[triI]);
                    vertI++;

                    label faceI = sampleElements_[triI];
                    meshTools::writeOBJ(str, mesh().faceCentres()[faceI]);
                    vertI++;
                    str << "l " << vertI-1 << ' ' << vertI << nl;
                }
            }
        }
    }

    needsUpdate_ = false;
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FWHSampledTriSurfaceMesh::FWHSampledTriSurfaceMesh
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    sampledSurface(name, refCast<const fvMesh>(obr), dict),
    obr_(obr),
    surface_
    (
        IOobject
        (
            dict.lookup("surface"),
            refCast<const fvMesh>(obr).time().constant(), // instance
            "triSurface",           // local
            refCast<const fvMesh>(obr),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    sampleSource_(samplingSourceNames_[dict.lookup("source")]),
    interpolationScheme_(dict.lookup("interpolationScheme")),
    needsUpdate_(true),
    sampleElements_(0),
    samplePoints_(0)
{
/*    if (!Pstream::parRun())
    {
    	update();
    	Info<<"Fwocs-Williams Hawkings control surface description:"<<endl;
    	print(Info);
    	Info<<endl;
    	//read(dict);
    }
    else{}*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FWHSampledTriSurfaceMesh::~FWHSampledTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const Foam::vectorField& Foam::FWHSampledTriSurfaceMesh::Sf() const
{
    return sampledSurface::Sf();
}


const Foam::scalarField& Foam::FWHSampledTriSurfaceMesh::magSf() const
{
    return sampledSurface::magSf();
}


const Foam::vectorField& Foam::FWHSampledTriSurfaceMesh::Cf() const
{
    return sampledSurface::Cf();
}


Foam::scalar Foam::FWHSampledTriSurfaceMesh::area() const
{
    return sampledSurface::area();
}


bool Foam::FWHSampledTriSurfaceMesh::needsUpdate() const
{
    return needsUpdate_;
}

void Foam::FWHSampledTriSurfaceMesh::read(const dictionary& dict)
{
}

bool Foam::FWHSampledTriSurfaceMesh::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    MeshStorage::clear();

    boundaryTreePtr_.clear();
    sampleElements_.clear();
    samplePoints_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::FWHSampledTriSurfaceMesh::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Mesh search engine, no triangulation of faces.
    meshSearch meshSearcher(mesh(), polyMesh::FACEPLANES);

    return update(meshSearcher);
}


bool Foam::FWHSampledTriSurfaceMesh::update(const treeBoundBox& bb)
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Mesh search engine on subset, no triangulation of faces.
    meshSearch meshSearcher(mesh(), bb, polyMesh::FACEPLANES);

    return update(meshSearcher);
}


Foam::tmp<Foam::scalarField> Foam::FWHSampledTriSurfaceMesh::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::FWHSampledTriSurfaceMesh::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField> Foam::FWHSampledTriSurfaceMesh::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::FWHSampledTriSurfaceMesh::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::FWHSampledTriSurfaceMesh::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::FWHSampledTriSurfaceMesh::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::FWHSampledTriSurfaceMesh::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::FWHSampledTriSurfaceMesh::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::FWHSampledTriSurfaceMesh::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::FWHSampledTriSurfaceMesh::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}

void Foam::FWHSampledTriSurfaceMesh::print(Ostream& os) const
{
    os  << "FWHSampledTriSurfaceMesh: " << name() << " :"
        << "  surface:" << surface_.objectRegistry::name()
        << "  faces:" << faces().size()
        << "  points:" << points().size()<<endl;
}

void Foam::FWHSampledTriSurfaceMesh::execute()
{
}

void Foam::FWHSampledTriSurfaceMesh::end()
{
    // Do nothing - only valid on write
}

void Foam::FWHSampledTriSurfaceMesh::timeSet()
{
    // Do nothing - only valid on write
}

void Foam::FWHSampledTriSurfaceMesh::write()
{
    /*Info<< "Calling calcFWHSourceTerms()"<<nl;
    calcFWHSourceTerms();
    
    if (size())
    {
        // finalize surfaces, merge points etc.
        update();

        const label nFields = classifyFields();

        if (Pstream::master())
        {
            if (debug)
            {
                Pout<< "Creating directory "
                    << outputPath_/mesh_.time().timeName() << nl << endl;

            }

            mkDir(outputPath_/mesh_.time().timeName());
        }

        // write geometry first if required,
        // or when no fields would otherwise be written
        if (nFields == 0 || formatter_->separateGeometry())
        {
            writeGeometry();
        }

        const IOobjectList objects(mesh_, mesh_.time().timeName());
//	Info<<"OBJECTS"<<objects<<nl;
        sampleAndWrite<volScalarField>(objects);
        sampleAndWrite<volVectorField>(objects);
        sampleAndWrite<volSphericalTensorField>(objects);
        sampleAndWrite<volSymmTensorField>(objects);
        sampleAndWrite<volTensorField>(objects);

        sampleAndWrite<surfaceScalarField>(objects);
        sampleAndWrite<surfaceVectorField>(objects);
        sampleAndWrite<surfaceSphericalTensorField>(objects);
        sampleAndWrite<surfaceSymmTensorField>(objects);
        sampleAndWrite<surfaceTensorField>(objects);
    }*/
}


void Foam::FWHSampledTriSurfaceMesh::updateMesh(const mapPolyMesh&)
{
    expire();

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::FWHSampledTriSurfaceMesh::movePoints(const polyMesh&)
{
    expire();
}


// ************************************************************************* //
