/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create 2D/extruded-2D meshes

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(pi, 3.14159265)

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(quad2D, ($1b $2b $2t $1t))
define(frontQuad, ($1t $2t $3t $4t))
define(backQuad, ($1b $4b $3b $2b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.1;

// Hub radius
define(r, 0.2)

// Impeller-tip radius
define(rb, 0.5)

// Baffle-tip radius
define(Rb, 0.7)

// Tank radius
define(R, 1)

// MRF region radius
define(ri, calc(0.5*(rb + Rb)))

// Thickness of 2D slab
define(z, 0.1)

// Base z
define(Zb, 0)

// Top z
define(Zt, calc(Zb + z))

// Number of cells radially between hub and impeller tip
define(Nr, 12)

// Number of cells radially in each of the two regions between
// impeller and baffle tips
define(Ni, 4)

// Number of cells radially between baffle tip and tank
define(NR, 12)

// Number of cells azimuthally in each of the 8 blocks
define(Na, 12)

// Number of cells in the thickness of the slab
define(Nz, 1)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

define(a0, 0)
define(a1, -45)
define(a2, -90)
define(a3, -135)
define(a4, 180)
define(a5, 135)
define(a6, 90)
define(a7, 45)

define(ea0, -22.5)
define(ea1, -67.5)
define(ea2, -112.5)
define(ea3, -157.5)
define(ea4, 157.5)
define(ea5, 112.5)
define(ea6, 67.5)
define(ea7, 22.5)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))
define(ca4, calc(cos((pi/180)*a4)))
define(ca5, calc(cos((pi/180)*a5)))
define(ca6, calc(cos((pi/180)*a6)))
define(ca7, calc(cos((pi/180)*a7)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))
define(sa4, calc(sin((pi/180)*a4)))
define(sa5, calc(sin((pi/180)*a5)))
define(sa6, calc(sin((pi/180)*a6)))
define(sa7, calc(sin((pi/180)*a7)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))
define(cea4, calc(cos((pi/180)*ea4)))
define(cea5, calc(cos((pi/180)*ea5)))
define(cea6, calc(cos((pi/180)*ea6)))
define(cea7, calc(cos((pi/180)*ea7)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))
define(sea4, calc(sin((pi/180)*ea4)))
define(sea5, calc(sin((pi/180)*ea5)))
define(sea6, calc(sin((pi/180)*ea6)))
define(sea7, calc(sin((pi/180)*ea7)))

define(x00, calc(r*ca0))
define(x01, calc(r*ca1))
define(x02, calc(r*ca2))
define(x03, calc(r*ca3))
define(x04, calc(r*ca4))
define(x05, calc(r*ca5))
define(x06, calc(r*ca6))
define(x07, calc(r*ca7))

define(x10, calc(rb*ca0))
define(x11, calc(rb*ca1))
define(x12, calc(rb*ca2))
define(x13, calc(rb*ca3))
define(x14, calc(rb*ca4))
define(x15, calc(rb*ca5))
define(x16, calc(rb*ca6))
define(x17, calc(rb*ca7))

define(x20, calc(ri*ca0))
define(x21, calc(ri*ca1))
define(x22, calc(ri*ca2))
define(x23, calc(ri*ca3))
define(x24, calc(ri*ca4))
define(x25, calc(ri*ca5))
define(x26, calc(ri*ca6))
define(x27, calc(ri*ca7))

define(x30, calc(Rb*ca0))
define(x31, calc(Rb*ca1))
define(x32, calc(Rb*ca2))
define(x33, calc(Rb*ca3))
define(x34, calc(Rb*ca4))
define(x35, calc(Rb*ca5))
define(x36, calc(Rb*ca6))
define(x37, calc(Rb*ca7))

define(x40, calc(R*ca0))
define(x41, calc(R*ca1))
define(x42, calc(R*ca2))
define(x43, calc(R*ca3))
define(x44, calc(R*ca4))
define(x45, calc(R*ca5))
define(x46, calc(R*ca6))
define(x47, calc(R*ca7))

define(y00, calc(r*sa0))
define(y01, calc(r*sa1))
define(y02, calc(r*sa2))
define(y03, calc(r*sa3))
define(y04, calc(r*sa4))
define(y05, calc(r*sa5))
define(y06, calc(r*sa6))
define(y07, calc(r*sa7))

define(y10, calc(rb*sa0))
define(y11, calc(rb*sa1))
define(y12, calc(rb*sa2))
define(y13, calc(rb*sa3))
define(y14, calc(rb*sa4))
define(y15, calc(rb*sa5))
define(y16, calc(rb*sa6))
define(y17, calc(rb*sa7))

define(y20, calc(ri*sa0))
define(y21, calc(ri*sa1))
define(y22, calc(ri*sa2))
define(y23, calc(ri*sa3))
define(y24, calc(ri*sa4))
define(y25, calc(ri*sa5))
define(y26, calc(ri*sa6))
define(y27, calc(ri*sa7))

define(y30, calc(Rb*sa0))
define(y31, calc(Rb*sa1))
define(y32, calc(Rb*sa2))
define(y33, calc(Rb*sa3))
define(y34, calc(Rb*sa4))
define(y35, calc(Rb*sa5))
define(y36, calc(Rb*sa6))
define(y37, calc(Rb*sa7))

define(y40, calc(R*sa0))
define(y41, calc(R*sa1))
define(y42, calc(R*sa2))
define(y43, calc(R*sa3))
define(y44, calc(R*sa4))
define(y45, calc(R*sa5))
define(y46, calc(R*sa6))
define(y47, calc(R*sa7))

define(ex00, calc(r*cea0))
define(ex01, calc(r*cea1))
define(ex02, calc(r*cea2))
define(ex03, calc(r*cea3))
define(ex04, calc(r*cea4))
define(ex05, calc(r*cea5))
define(ex06, calc(r*cea6))
define(ex07, calc(r*cea7))

define(ex10, calc(rb*cea0))
define(ex11, calc(rb*cea1))
define(ex12, calc(rb*cea2))
define(ex13, calc(rb*cea3))
define(ex14, calc(rb*cea4))
define(ex15, calc(rb*cea5))
define(ex16, calc(rb*cea6))
define(ex17, calc(rb*cea7))

define(ex20, calc(ri*cea0))
define(ex21, calc(ri*cea1))
define(ex22, calc(ri*cea2))
define(ex23, calc(ri*cea3))
define(ex24, calc(ri*cea4))
define(ex25, calc(ri*cea5))
define(ex26, calc(ri*cea6))
define(ex27, calc(ri*cea7))

define(ex30, calc(Rb*cea0))
define(ex31, calc(Rb*cea1))
define(ex32, calc(Rb*cea2))
define(ex33, calc(Rb*cea3))
define(ex34, calc(Rb*cea4))
define(ex35, calc(Rb*cea5))
define(ex36, calc(Rb*cea6))
define(ex37, calc(Rb*cea7))

define(ex40, calc(R*cea0))
define(ex41, calc(R*cea1))
define(ex42, calc(R*cea2))
define(ex43, calc(R*cea3))
define(ex44, calc(R*cea4))
define(ex45, calc(R*cea5))
define(ex46, calc(R*cea6))
define(ex47, calc(R*cea7))

define(ey00, calc(r*sea0))
define(ey01, calc(r*sea1))
define(ey02, calc(r*sea2))
define(ey03, calc(r*sea3))
define(ey04, calc(r*sea4))
define(ey05, calc(r*sea5))
define(ey06, calc(r*sea6))
define(ey07, calc(r*sea7))

define(ey10, calc(rb*sea0))
define(ey11, calc(rb*sea1))
define(ey12, calc(rb*sea2))
define(ey13, calc(rb*sea3))
define(ey14, calc(rb*sea4))
define(ey15, calc(rb*sea5))
define(ey16, calc(rb*sea6))
define(ey17, calc(rb*sea7))

define(ey20, calc(ri*sea0))
define(ey21, calc(ri*sea1))
define(ey22, calc(ri*sea2))
define(ey23, calc(ri*sea3))
define(ey24, calc(ri*sea4))
define(ey25, calc(ri*sea5))
define(ey26, calc(ri*sea6))
define(ey27, calc(ri*sea7))

define(ey30, calc(Rb*sea0))
define(ey31, calc(Rb*sea1))
define(ey32, calc(Rb*sea2))
define(ey33, calc(Rb*sea3))
define(ey34, calc(Rb*sea4))
define(ey35, calc(Rb*sea5))
define(ey36, calc(Rb*sea6))
define(ey37, calc(Rb*sea7))

define(ey40, calc(R*sea0))
define(ey41, calc(R*sea1))
define(ey42, calc(R*sea2))
define(ey43, calc(R*sea3))
define(ey44, calc(R*sea4))
define(ey45, calc(R*sea5))
define(ey46, calc(R*sea6))
define(ey47, calc(R*sea7))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(r0b)
    vert(0, 0, Zb) vlabel(r0sb)
    vert(0, 1, Zb) vlabel(r1b)
    vert(0, 2, Zb) vlabel(r2b)
    vert(0, 2, Zb) vlabel(r2sb)
    vert(0, 3, Zb) vlabel(r3b)
    vert(0, 4, Zb) vlabel(r4b)
    vert(0, 4, Zb) vlabel(r4sb)
    vert(0, 5, Zb) vlabel(r5b)
    vert(0, 6, Zb) vlabel(r6b)
    vert(0, 6, Zb) vlabel(r6sb)
    vert(0, 7, Zb) vlabel(r7b)

    vert(1, 0, Zb) vlabel(rb0b)
    vert(1, 1, Zb) vlabel(rb1b)
    vert(1, 2, Zb) vlabel(rb2b)
    vert(1, 3, Zb) vlabel(rb3b)
    vert(1, 4, Zb) vlabel(rb4b)
    vert(1, 5, Zb) vlabel(rb5b)
    vert(1, 6, Zb) vlabel(rb6b)
    vert(1, 7, Zb) vlabel(rb7b)

    vert(2, 0, Zb) vlabel(rii0b)
    vert(2, 0, Zb) vlabel(rie0b)
    vert(2, 1, Zb) vlabel(rii1b)
    vert(2, 1, Zb) vlabel(rie1b)
    vert(2, 2, Zb) vlabel(rii2b)
    vert(2, 2, Zb) vlabel(rie2b)
    vert(2, 3, Zb) vlabel(rii3b)
    vert(2, 3, Zb) vlabel(rie3b)
    vert(2, 4, Zb) vlabel(rii4b)
    vert(2, 4, Zb) vlabel(rie4b)
    vert(2, 5, Zb) vlabel(rii5b)
    vert(2, 5, Zb) vlabel(rie5b)
    vert(2, 6, Zb) vlabel(rii6b)
    vert(2, 6, Zb) vlabel(rie6b)
    vert(2, 7, Zb) vlabel(rii7b)
    vert(2, 7, Zb) vlabel(rie7b)

    vert(3, 0, Zb) vlabel(Rb0b)
    vert(3, 1, Zb) vlabel(Rb1b)
    vert(3, 2, Zb) vlabel(Rb2b)
    vert(3, 3, Zb) vlabel(Rb3b)
    vert(3, 4, Zb) vlabel(Rb4b)
    vert(3, 5, Zb) vlabel(Rb5b)
    vert(3, 6, Zb) vlabel(Rb6b)
    vert(3, 7, Zb) vlabel(Rb7b)

    vert(4, 0, Zb) vlabel(R0b)
    vert(4, 1, Zb) vlabel(R1b)
    vert(4, 1, Zb) vlabel(R1sb)
    vert(4, 2, Zb) vlabel(R2b)
    vert(4, 3, Zb) vlabel(R3b)
    vert(4, 3, Zb) vlabel(R3sb)
    vert(4, 4, Zb) vlabel(R4b)
    vert(4, 5, Zb) vlabel(R5b)
    vert(4, 5, Zb) vlabel(R5sb)
    vert(4, 6, Zb) vlabel(R6b)
    vert(4, 7, Zb) vlabel(R7b)
    vert(4, 7, Zb) vlabel(R7sb)

    vert(0, 0, Zt) vlabel(r0t)
    vert(0, 0, Zt) vlabel(r0st)
    vert(0, 1, Zt) vlabel(r1t)
    vert(0, 2, Zt) vlabel(r2t)
    vert(0, 2, Zt) vlabel(r2st)
    vert(0, 3, Zt) vlabel(r3t)
    vert(0, 4, Zt) vlabel(r4t)
    vert(0, 4, Zt) vlabel(r4st)
    vert(0, 5, Zt) vlabel(r5t)
    vert(0, 6, Zt) vlabel(r6t)
    vert(0, 6, Zt) vlabel(r6st)
    vert(0, 7, Zt) vlabel(r7t)

    vert(1, 0, Zt) vlabel(rb0t)
    vert(1, 1, Zt) vlabel(rb1t)
    vert(1, 2, Zt) vlabel(rb2t)
    vert(1, 3, Zt) vlabel(rb3t)
    vert(1, 4, Zt) vlabel(rb4t)
    vert(1, 5, Zt) vlabel(rb5t)
    vert(1, 6, Zt) vlabel(rb6t)
    vert(1, 7, Zt) vlabel(rb7t)

    vert(2, 0, Zt) vlabel(rii0t)
    vert(2, 0, Zt) vlabel(rie0t)
    vert(2, 1, Zt) vlabel(rii1t)
    vert(2, 1, Zt) vlabel(rie1t)
    vert(2, 2, Zt) vlabel(rii2t)
    vert(2, 2, Zt) vlabel(rie2t)
    vert(2, 3, Zt) vlabel(rii3t)
    vert(2, 3, Zt) vlabel(rie3t)
    vert(2, 4, Zt) vlabel(rii4t)
    vert(2, 4, Zt) vlabel(rie4t)
    vert(2, 5, Zt) vlabel(rii5t)
    vert(2, 5, Zt) vlabel(rie5t)
    vert(2, 6, Zt) vlabel(rii6t)
    vert(2, 6, Zt) vlabel(rie6t)
    vert(2, 7, Zt) vlabel(rii7t)
    vert(2, 7, Zt) vlabel(rie7t)

    vert(3, 0, Zt) vlabel(Rb0t)
    vert(3, 1, Zt) vlabel(Rb1t)
    vert(3, 2, Zt) vlabel(Rb2t)
    vert(3, 3, Zt) vlabel(Rb3t)
    vert(3, 4, Zt) vlabel(Rb4t)
    vert(3, 5, Zt) vlabel(Rb5t)
    vert(3, 6, Zt) vlabel(Rb6t)
    vert(3, 7, Zt) vlabel(Rb7t)

    vert(4, 0, Zt) vlabel(R0t)
    vert(4, 1, Zt) vlabel(R1t)
    vert(4, 1, Zt) vlabel(R1st)
    vert(4, 2, Zt) vlabel(R2t)
    vert(4, 3, Zt) vlabel(R3t)
    vert(4, 3, Zt) vlabel(R3st)
    vert(4, 4, Zt) vlabel(R4t)
    vert(4, 5, Zt) vlabel(R5t)
    vert(4, 5, Zt) vlabel(R5st)
    vert(4, 6, Zt) vlabel(R6t)
    vert(4, 7, Zt) vlabel(R7t)
    vert(4, 7, Zt) vlabel(R7st)
);

blocks
(
    // block0
    hex2D(r0, r1, rb1, rb0)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(r1, r2s, rb2, rb1)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(r2, r3, rb3, rb2)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(r3, r4s, rb4, rb3)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(r4, r5, rb5, rb4)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(r5, r6s, rb6, rb5)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(r6, r7, rb7, rb6)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(r7, r0s, rb0, rb7)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(rb0, rb1, rii1, rii0)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(rb1, rb2, rii2, rii1)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(rb2, rb3, rii3, rii2)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(rb3, rb4, rii4, rii3)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(rb4, rb5, rii5, rii4)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(rb5, rb6, rii6, rii5)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(rb6, rb7, rii7, rii6)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(rb7, rb0, rii0, rii7)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(rie0, rie1, Rb1, Rb0)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(rie1, rie2, Rb2, Rb1)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(rie2, rie3, Rb3, Rb2)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(rie3, rie4, Rb4, Rb3)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(rie4, rie5, Rb5, Rb4)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(rie5, rie6, Rb6, Rb5)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(rie6, rie7, Rb7, Rb6)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(rie7, rie0, Rb0, Rb7)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(Rb0, Rb1, R1s, R0)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(Rb1, Rb2, R2, R1)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(Rb2, Rb3, R3s, R2)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(Rb3, Rb4, R4, R3)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(Rb4, Rb5, R5s, R4)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(Rb5, Rb6, R6, R5)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(Rb6, Rb7, R7s, R6)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(Rb7, Rb0, R0, R7)
    (Na NR Nz)
    simpleGrading (1 1 1)
);

edges
(
    arc r0b r1b evert(0, 0, Zb)
    arc r1b r2sb evert(0, 1, Zb)
    arc r2b r3b evert(0, 2, Zb)
    arc r3b r4sb evert(0, 3, Zb)
    arc r4b r5b evert(0, 4, Zb)
    arc r5b r6sb evert(0, 5, Zb)
    arc r6b r7b evert(0, 6, Zb)
    arc r7b r0sb evert(0, 7, Zb)

    arc rb0b rb1b evert(1, 0, Zb)
    arc rb1b rb2b evert(1, 1, Zb)
    arc rb2b rb3b evert(1, 2, Zb)
    arc rb3b rb4b evert(1, 3, Zb)
    arc rb4b rb5b evert(1, 4, Zb)
    arc rb5b rb6b evert(1, 5, Zb)
    arc rb6b rb7b evert(1, 6, Zb)
    arc rb7b rb0b evert(1, 7, Zb)

    arc rii0b rii1b evert(2, 0, Zb)
    arc rie0b rie1b evert(2, 0, Zb)
    arc rii1b rii2b evert(2, 1, Zb)
    arc rie1b rie2b evert(2, 1, Zb)
    arc rii2b rii3b evert(2, 2, Zb)
    arc rie2b rie3b evert(2, 2, Zb)
    arc rii3b rii4b evert(2, 3, Zb)
    arc rie3b rie4b evert(2, 3, Zb)
    arc rii4b rii5b evert(2, 4, Zb)
    arc rie4b rie5b evert(2, 4, Zb)
    arc rii5b rii6b evert(2, 5, Zb)
    arc rie5b rie6b evert(2, 5, Zb)
    arc rii6b rii7b evert(2, 6, Zb)
    arc rie6b rie7b evert(2, 6, Zb)
    arc rii7b rii0b evert(2, 7, Zb)
    arc rie7b rie0b evert(2, 7, Zb)

    arc Rb0b Rb1b evert(3, 0, Zb)
    arc Rb1b Rb2b evert(3, 1, Zb)
    arc Rb2b Rb3b evert(3, 2, Zb)
    arc Rb3b Rb4b evert(3, 3, Zb)
    arc Rb4b Rb5b evert(3, 4, Zb)
    arc Rb5b Rb6b evert(3, 5, Zb)
    arc Rb6b Rb7b evert(3, 6, Zb)
    arc Rb7b Rb0b evert(3, 7, Zb)

    arc R0b R1sb evert(4, 0, Zb)
    arc R1b R2b evert(4, 1, Zb)
    arc R2b R3sb evert(4, 2, Zb)
    arc R3b R4b evert(4, 3, Zb)
    arc R4b R5sb evert(4, 4, Zb)
    arc R5b R6b evert(4, 5, Zb)
    arc R6b R7sb evert(4, 6, Zb)
    arc R7b R0b evert(4, 7, Zb)

    arc r0t r1t evert(0, 0, Zt)
    arc r1t r2st evert(0, 1, Zt)
    arc r2t r3t evert(0, 2, Zt)
    arc r3t r4st evert(0, 3, Zt)
    arc r4t r5t evert(0, 4, Zt)
    arc r5t r6st evert(0, 5, Zt)
    arc r6t r7t evert(0, 6, Zt)
    arc r7t r0st evert(0, 7, Zt)

    arc rb0t rb1t evert(1, 0, Zt)
    arc rb1t rb2t evert(1, 1, Zt)
    arc rb2t rb3t evert(1, 2, Zt)
    arc rb3t rb4t evert(1, 3, Zt)
    arc rb4t rb5t evert(1, 4, Zt)
    arc rb5t rb6t evert(1, 5, Zt)
    arc rb6t rb7t evert(1, 6, Zt)
    arc rb7t rb0t evert(1, 7, Zt)

    arc rii0t rii1t evert(2, 0, Zt)
    arc rie0t rie1t evert(2, 0, Zt)
    arc rii1t rii2t evert(2, 1, Zt)
    arc rie1t rie2t evert(2, 1, Zt)
    arc rii2t rii3t evert(2, 2, Zt)
    arc rie2t rie3t evert(2, 2, Zt)
    arc rii3t rii4t evert(2, 3, Zt)
    arc rie3t rie4t evert(2, 3, Zt)
    arc rii4t rii5t evert(2, 4, Zt)
    arc rie4t rie5t evert(2, 4, Zt)
    arc rii5t rii6t evert(2, 5, Zt)
    arc rie5t rie6t evert(2, 5, Zt)
    arc rii6t rii7t evert(2, 6, Zt)
    arc rie6t rie7t evert(2, 6, Zt)
    arc rii7t rii0t evert(2, 7, Zt)
    arc rie7t rie0t evert(2, 7, Zt)

    arc Rb0t Rb1t evert(3, 0, Zt)
    arc Rb1t Rb2t evert(3, 1, Zt)
    arc Rb2t Rb3t evert(3, 2, Zt)
    arc Rb3t Rb4t evert(3, 3, Zt)
    arc Rb4t Rb5t evert(3, 4, Zt)
    arc Rb5t Rb6t evert(3, 5, Zt)
    arc Rb6t Rb7t evert(3, 6, Zt)
    arc Rb7t Rb0t evert(3, 7, Zt)

    arc R0t R1st evert(4, 0, Zt)
    arc R1t R2t evert(4, 1, Zt)
    arc R2t R3st evert(4, 2, Zt)
    arc R3t R4t evert(4, 3, Zt)
    arc R4t R5st evert(4, 4, Zt)
    arc R5t R6t evert(4, 5, Zt)
    arc R6t R7st evert(4, 6, Zt)
    arc R7t R0t evert(4, 7, Zt)
);

patches
(
    wall rotor
    (
        quad2D(r0, r1)
        quad2D(r1, r2s)
        quad2D(r2, r3)
        quad2D(r3, r4s)
        quad2D(r4, r5)
        quad2D(r5, r6s)
        quad2D(r6, r7)
        quad2D(r7, r0s)

        quad2D(r0, rb0)
        quad2D(r0s, rb0)

        quad2D(r2, rb2)
        quad2D(r2s, rb2)

        quad2D(r4, rb4)
        quad2D(r4s, rb4)

        quad2D(r6, rb6)
        quad2D(r6s, rb6)
    )

    wall stator
    (
        quad2D(R0, R1s)
        quad2D(R1, R2)
        quad2D(R2, R3s)
        quad2D(R3, R4)
        quad2D(R4, R5s)
        quad2D(R5, R6)
        quad2D(R6, R7s)
        quad2D(R7, R0)

        quad2D(R1, Rb1)
        quad2D(R1s, Rb1)

        quad2D(R3, Rb3)
        quad2D(R3s, Rb3)

        quad2D(R5, Rb5)
        quad2D(R5s, Rb5)

        quad2D(R7, Rb7)
        quad2D(R7s, Rb7)
    )

    empty front
    (
        frontQuad(r0, r1, rb1, rb0)
        frontQuad(r1, r2s, rb2, rb1)
        frontQuad(r2, r3, rb3, rb2)
        frontQuad(r3, r4s, rb4, rb3)
        frontQuad(r4, r5, rb5, rb4)
        frontQuad(r5, r6s, rb6, rb5)
        frontQuad(r6, r7, rb7, rb6)
        frontQuad(r7, r0s, rb0, rb7)
        frontQuad(rb0, rb1, rii1, rii0)
        frontQuad(rb1, rb2, rii2, rii1)
        frontQuad(rb2, rb3, rii3, rii2)
        frontQuad(rb3, rb4, rii4, rii3)
        frontQuad(rb4, rb5, rii5, rii4)
        frontQuad(rb5, rb6, rii6, rii5)
        frontQuad(rb6, rb7, rii7, rii6)
        frontQuad(rb7, rb0, rii0, rii7)
        frontQuad(rie0, rie1, Rb1, Rb0)
        frontQuad(rie1, rie2, Rb2, Rb1)
        frontQuad(rie2, rie3, Rb3, Rb2)
        frontQuad(rie3, rie4, Rb4, Rb3)
        frontQuad(rie4, rie5, Rb5, Rb4)
        frontQuad(rie5, rie6, Rb6, Rb5)
        frontQuad(rie6, rie7, Rb7, Rb6)
        frontQuad(rie7, rie0, Rb0, Rb7)
        frontQuad(Rb0, Rb1, R1s, R0)
        frontQuad(Rb1, Rb2, R2, R1)
        frontQuad(Rb2, Rb3, R3s, R2)
        frontQuad(Rb3, Rb4, R4, R3)
        frontQuad(Rb4, Rb5, R5s, R4)
        frontQuad(Rb5, Rb6, R6, R5)
        frontQuad(Rb6, Rb7, R7s, R6)
        frontQuad(Rb7, Rb0, R0, R7)
    )

    empty back
    (
        backQuad(r0, r1, rb1, rb0)
        backQuad(r1, r2s, rb2, rb1)
        backQuad(r2, r3, rb3, rb2)
        backQuad(r3, r4s, rb4, rb3)
        backQuad(r4, r5, rb5, rb4)
        backQuad(r5, r6s, rb6, rb5)
        backQuad(r6, r7, rb7, rb6)
        backQuad(r7, r0s, rb0, rb7)
        backQuad(rb0, rb1, rii1, rii0)
        backQuad(rb1, rb2, rii2, rii1)
        backQuad(rb2, rb3, rii3, rii2)
        backQuad(rb3, rb4, rii4, rii3)
        backQuad(rb4, rb5, rii5, rii4)
        backQuad(rb5, rb6, rii6, rii5)
        backQuad(rb6, rb7, rii7, rii6)
        backQuad(rb7, rb0, rii0, rii7)
        backQuad(rie0, rie1, Rb1, Rb0)
        backQuad(rie1, rie2, Rb2, Rb1)
        backQuad(rie2, rie3, Rb3, Rb2)
        backQuad(rie3, rie4, Rb4, Rb3)
        backQuad(rie4, rie5, Rb5, Rb4)
        backQuad(rie5, rie6, Rb6, Rb5)
        backQuad(rie6, rie7, Rb7, Rb6)
        backQuad(rie7, rie0, Rb0, Rb7)
        backQuad(Rb0, Rb1, R1s, R0)
        backQuad(Rb1, Rb2, R2, R1)
        backQuad(Rb2, Rb3, R3s, R2)
        backQuad(Rb3, Rb4, R4, R3)
        backQuad(Rb4, Rb5, R5s, R4)
        backQuad(Rb5, Rb6, R6, R5)
        backQuad(Rb6, Rb7, R7s, R6)
        backQuad(Rb7, Rb0, R0, R7)
    )

    patch MixingPlaneInt
    (
        quad2D(rii0, rii1)
        quad2D(rii1, rii2)
        quad2D(rii2, rii3)
        quad2D(rii3, rii4)
        quad2D(rii4, rii5)
        quad2D(rii5, rii6)
        quad2D(rii6, rii7)
        quad2D(rii7, rii0)
    )

    patch MixingPlaneExt
    (
        quad2D(rie0, rie1)
        quad2D(rie1, rie2)
        quad2D(rie2, rie3)
        quad2D(rie3, rie4)
        quad2D(rie4, rie5)
        quad2D(rie5, rie6)
        quad2D(rie6, rie7)
        quad2D(rie7, rie0)
    )
);

// ************************************************************************* //
