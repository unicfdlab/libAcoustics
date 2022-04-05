Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);
Foam::word regionName(Foam::polyMesh::defaultRegion);

List<vector> points =
{
{-5,-0.1,-0.1},
{-3,-0.1,-0.1},
{-1,-0.1,-0.1},
{1,-0.1,-0.1},
{3,-0.1,-0.1},
{5,-0.1,-0.1},
{-5,0.1,-0.1},
{-3,0.1,-0.1},
{-1,0.1,-0.1},
{1,0.1,-0.1},
{3,0.1,-0.1},
{5,0.1,-0.1},
{-5,-0.1,0.1},
{-3,-0.1,0.1},
{-1,-0.1,0.1},
{1,-0.1,0.1},
{3,-0.1,0.1},
{5,-0.1,0.1},
{-5,0.1,0.1},
{-3,0.1,0.1},
{-1,0.1,0.1},
{1,0.1,0.1},
{3,0.1,0.1},
{5,0.1,0.1}
};

faceList faces  = 
{
face{1,7,19,13},
face{2,8,20,14},
face{3,9,21,15},
face{4,10,22,16},
face{5,11,23,17},
face{0,12,18,6},
face{0,1,13,12},
face{1,2,14,13},
face{2,3,15,14},
face{3,4,16,15},
face{4,5,17,16},
face{12,13,19,18},
face{13,14,20,19},
face{14,15,21,20},
face{15,16,22,21},
face{16,17,23,22},
face{6,18,19,7},
face{7,19,20,8},
face{8,20,21,9},
face{9,21,22,10},
face{10,22,23,11},
face{0,6,7,1},
face{1,7,8,2},
face{2,8,9,3},
face{3,9,10,4},
face{4,10,11,5},
};
List<label> owner = 
{
0,
1,
2,
3,
4,
0,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4,
0,
1,
2,
3,
4
};

List<label> neighbour = {1,2,3,4};

pointField pf (points);
faceList& fl = faces;
labelList& own = owner;
labelList& nei = neighbour;

meshPtr.reset
( 
    new Foam::fvMesh
    (
        IOobject
        (
            regionName,
            timeObj.timeName(),
            timeObj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        std::move(pf),
        std::move(fl),
        std::move(own),
        std::move(nei)
    )
);

PtrList<polyPatch> patches (2);

patches.set
(
    0,
    new polyPatch("testPatch", 2, 4, 0, meshPtr().boundaryMesh(), "patch")
);
patches.set
(
    1,
    new polyPatch("empty", 20, 6, 1, meshPtr().boundaryMesh(), "patch")
);

meshPtr().addFvPatches(patches);

meshPtr().init(true);
