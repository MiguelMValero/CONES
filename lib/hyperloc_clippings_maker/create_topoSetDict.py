#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 11:43:04 2023

@author: villanul
"""

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys
import csv
import math
from pylab import cm
import pandas as pd
csv.field_size_limit(sys.maxsize)

lref = 0.029
decmls = 4

filename = "obs_coordinates.txt"
data = pd.read_csv(filename, header=None)
nb_o = len(data.index)


with open('topoSetDict', 'w') as f:
    f.write('/*--------------------------------*- C++ -*----------------------------------*\\\n')
    f.write('| =========                 |                                                 |\n')
    f.write('| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
    f.write('|  \\\    /   O peration     | Version:  v2306                                 |\n')
    f.write('|   \\\  /    A nd           | Website:  www.openfoam.com                      |\n')
    f.write('|    \\\/     M anipulation  |                                                 |\n')
    f.write('\*---------------------------------------------------------------------------*/\n')
    f.write('FoamFile\n')
    f.write('{\n')
    f.write('    version     2.0;\n')
    f.write('    format      ascii;\n')
    f.write('    class       dictionary;\n')
    f.write('    object      topoSetDict;\n')
    f.write('}\n')
    f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
    f.write('\nactions\n')
    f.write('(\n')
    for i in range(nb_o):
        f.write('    {\n')
        f.write('        name    clippingCells'+ str(i) +'.txt;\n')
        f.write('        type    cellSet;\n')
        f.write('        action  new;\n')
        f.write('        source  boxToCell;\n')
        f.write('        sourceInfo\n')
        f.write('        {\n')
        f.write('            boxes\n')
        f.write('            (\n')
        obs_x = data.iloc[i,0]
        obs_y = data.iloc[i,1]
        obs_z = data.iloc[i,2]
        minx = '{:.{decmls}f}'.format(obs_x - lref , decmls=decmls)
        miny = '{:.{decmls}f}'.format(obs_y - lref , decmls=decmls)
        minz = '{:.{decmls}f}'.format(obs_z - lref , decmls=decmls)
        maxx = '{:.{decmls}f}'.format(obs_x + lref , decmls=decmls)
        maxy = '{:.{decmls}f}'.format(obs_y + lref , decmls=decmls)
        maxz = '{:.{decmls}f}'.format(obs_z + lref , decmls=decmls)
        f.write('                ('+minx+' '+miny+' '+minz+')('+maxx+' '+maxy+' '+maxz+')\n')
        f.write('            );\n')
        f.write('        }\n')
        f.write('    }\n')
        f.write('\n')
    f.write(');\n')
    f.write('\n// ************************************************************************* //\n')


    