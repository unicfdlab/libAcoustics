/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2012                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dnl> -----------------------------------------------------------------
dnl> <STANDARD DEFINTIONS>
dnl>
changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
define(VCOUNT, 0)  dnl>
define(vlabel, [[// ]pt VCOUNT ($1) define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])  dnl>
dnl>
dnl> </STANDARD DEFINTIONS>
dnl> -----------------------------------------------------------------
dnl>
define(r1, 1) dnl>
define(r2, 2) dnl>
define(xr, calc(10*r1)) dnl>
define(h, calc(30*r1)) dnl>
define(l, calc(xr + 30*r1)) dnl>
define(w, 1) dnl>
define(Sin, 0.7071067812) dnl>   == sin(45) dnl>
dnl>
define(y0, 0) dnl>
define(y1, calc(h/2 + y0 - r2*Sin)) dnl>
define(y2, calc(h/2 + y0 - r1*Sin)) dnl>
define(y3, calc(h/2 + y0 + r1*Sin)) dnl>
define(y4, calc(h/2 + y0 + r2*Sin)) dnl>
define(y5, calc(y0 + h)) dnl>
dnl>
define(x0, 0) dnl>
define(x1, calc(xr + x0 - r2*Sin)) dnl>
define(x2, calc(xr + x0 - r1*Sin)) dnl>
define(x3, calc(xr + x0 + r1*Sin)) dnl>
define(x4, calc(xr + x0 + r2*Sin)) dnl>
define(x5, calc(x0 + l)) dnl>
dnl>
define(w0, -1) dnl>
define(w1, calc(w0 + w)) dnl>
dnl>
define(N, 20) dnl>
define(hex2D, hex ($1f $2f $3f $4f $1b $2b $3b $4b)) dnl>
dnl>
dnl>
dnl>
define(nScale, 4)
define(nx1, calc(nScale*30))  dnl>
define(nx2, calc(nScale*10))  dnl>
define(nx3, calc(nScale*100))  dnl>
define(ny1, calc(nScale*20))  dnl>
define(ny2, nx2)  dnl>
define(ny3, ny1)  dnl>
define(nr, calc(nScale*20))  dnl>
define(nz, 1)  dnl>
dnl>
dnl>
define(hex2D, hex ($1b $2b $3b $4b $1f $2f $3f $4f)) dnl>
define(quad2D, ($1f $1b $2b $2f))  dnl>
define(frontQuad, ($1f $2f $3f $4f)) dnl>
define(backQuad, ($4b $3b $2b $1b)) dnl>
dnl>
dnl> </STANDARD DEFINTIONS>
dnl> -----------------------------------------------------------------
dnl>
dnl>

scale   0.02;

vertices
(
    (x0 y0 w0)  vlabel(p0f)
    (x1 y0 w0)  vlabel(p1f)
    (x4 y0 w0)  vlabel(p2f)
    (x5 y0 w0)  vlabel(p3f)
    (x5 y1 w0)  vlabel(p4f)
    (x5 y4 w0)  vlabel(p5f)
    (x5 y5 w0)  vlabel(p6f)
    (x4 y5 w0)  vlabel(p7f)
    (x1 y5 w0)  vlabel(p8f)
    (x0 y5 w0)  vlabel(p9f)
    (x0 y4 w0)  vlabel(p10f)
    (x0 y1 w0)  vlabel(p11f)
    (x1 y1 w0)  vlabel(p12f)
    (x4 y1 w0)  vlabel(p13f)
    (x4 y4 w0)  vlabel(p14f)
    (x1 y4 w0)  vlabel(p15f)
    (x2 y2 w0)  vlabel(p16f)
    (x3 y2 w0)  vlabel(p17f)
    (x3 y3 w0)  vlabel(p18f)
    (x2 y3 w0)  vlabel(p19f)

    (x0 y0 w1)  vlabel(p0b)
    (x1 y0 w1)  vlabel(p1b)
    (x4 y0 w1)  vlabel(p2b)
    (x5 y0 w1)  vlabel(p3b)
    (x5 y1 w1)  vlabel(p4b)
    (x5 y4 w1)  vlabel(p5b)
    (x5 y5 w1)  vlabel(p6b)
    (x4 y5 w1)  vlabel(p7b)
    (x1 y5 w1)  vlabel(p8b)
    (x0 y5 w1)  vlabel(p9b)
    (x0 y4 w1)  vlabel(p10b)
    (x0 y1 w1)  vlabel(p11b)
    (x1 y1 w1)  vlabel(p12b)
    (x4 y1 w1)  vlabel(p13b)
    (x4 y4 w1)  vlabel(p14b)
    (x1 y4 w1)  vlabel(p15b)
    (x2 y2 w1)  vlabel(p16b)
    (x3 y2 w1)  vlabel(p17b)
    (x3 y3 w1)  vlabel(p18b)
    (x2 y3 w1)  vlabel(p19b)
);

blocks
(
    // 0
    hex2D(p0, p1, p12, p11)
    ( nx1 ny1 nz ) simpleGrading (1 0.2 1)

    // 1
    hex2D(p1, p2, p13, p12)
    ( nx2 ny1 nz ) simpleGrading (1 0.2 1)

    // 2
    hex2D(p2, p3, p4, p13)
    ( nx3 ny1 nz ) simpleGrading (1 0.2 1)

    // 3
    hex2D(p13, p4, p5, p14)
    ( nx3 ny2 nz ) simpleGrading (1 1 1)

    // 4
    hex2D(p14, p5, p6, p7)
    ( nx3 ny3 nz ) simpleGrading (1 5 1)

    // 5
    hex2D(p15, p14, p7, p8)
    ( nx2 ny3 nz ) simpleGrading (1 5 1)

    // 6
    hex2D(p10, p15, p8, p9)
    ( nx1 ny3 nz ) simpleGrading (1 5 1)

    // 7
    hex2D(p11, p12, p15, p10)
    ( nx1 ny2 nz ) simpleGrading (1 1 1)

    // 8
    hex2D(p12, p16, p19, p15)
    ( nr ny2 nz ) simpleGrading (0.025 1 1)

    // 9
    hex2D(p12, p13, p17, p16)
    ( nx2 nr nz ) simpleGrading (1 0.025 1)

    // 10
    hex2D(p17, p13, p14, p18)
    ( nr ny2 nz ) simpleGrading (40 1 1)

    // 11
    hex2D(p19, p18, p14, p15)
    ( nx2 nr nz ) simpleGrading (1 40 1)
);

edges
(
    // Inner circle
    arc  16 17 (calc(xr + x0) calc(h/2 + y0 - r1) w0)
    arc  17 18 (calc(xr + x0 + r1) calc(h/2 + y0) w0)
    arc  18 19 (calc(xr + x0) calc(h/2 + y0 + r1) w0)
    arc  19 16 (calc(xr + x0 - r1) calc(h/2 + y0) w0)

    arc  36 37 (calc(xr + x0) calc(h/2 + y0 - r1) w1)
    arc  37 38 (calc(xr + x0 + r1) calc(h/2 + y0) w1)
    arc  38 39 (calc(xr + x0) calc(h/2 + y0 + r1) w1)
    arc  39 36 (calc(xr + x0 - r1) calc(h/2 + y0) w1)

    // Outer circle
    arc  12 13 (calc(xr + x0) calc(h/2 + y0 - r2) w0)
    arc  13 14 (calc(xr + x0 + r2) calc(h/2 + y0) w0)
    arc  14 15 (calc(xr + x0) calc(h/2 + y0 + r2) w0)
    arc  15 12 (calc(xr + x0 - r2) calc(h/2 + y0) w0)

    arc  32 33 (calc(xr + x0) calc(h/2 + y0 - r2) w1)
    arc  33 34 (calc(xr + x0 + r2) calc(h/2 + y0) w1)
    arc  34 35 (calc(xr + x0) calc(h/2 + y0 + r2) w1)
    arc  35 32 (calc(xr + x0 - r2) calc(h/2 + y0) w1)
);

defaultPatch
{
    name    frontAndBack;
    type    empty;
}

boundary
(
    inlet
    {
        type        patch;
        faces
        (
            quad2D(p9, p10)
            quad2D(p10, p11)
            quad2D(p11, p0)
        );
    }
    outlet
    {
        type        patch;
        faces
        (
            quad2D(p3, p4)
            quad2D(p4, p5)
            quad2D(p5, p6)
        );
    }

    cylinder
    {
        type        wall;
        faces
        (
            quad2D(p16, p17)
            quad2D(p17, p18)
            quad2D(p18, p19)
            quad2D(p19, p16)
        );
    }

    top
    {
        type        patch;
        faces
        (
            quad2D(p6, p7)
            quad2D(p7, p8)
            quad2D(p8, p9)
        );
    }

    bottom
    {
        type        patch;
        faces
        (
            quad2D(p0, p1)
            quad2D(p1, p2)
            quad2D(p2, p3)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
