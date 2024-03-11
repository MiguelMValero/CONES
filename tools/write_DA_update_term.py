#!/usr/bin/python3

#############################################################
#   Script computing and writing the update term
#             Written by Sarp - 2024
#############################################################

import numpy as np
import time
import pandas as pd
import os

import toolboxGeneral

#========================= FUNCTIONS =============================

def writeVectorFieldOpenFoam(varr, name_field, dir_out):
	# Writes the given vector field in the OpenFoam format
	
	flnm = "DA_Updt"

	with open(f"{dir_out}/{flnm}", 'w') as f:
		f.write('/*--------------------------------*- C++ -*----------------------------------*\\\n')
		f.write('|  =========                 |                                                 |\n')
		f.write('| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
		f.write('|  \\\    /   O peration     | Website:  https://openfoam.org                  |\n')
		f.write('|   \\\  /    A nd           | Version:  9                                     |\n')
		f.write('|    \\\/     M anipulation  |                                                 |\n')
		f.write('\*---------------------------------------------------------------------------*/\n')
		f.write('FoamFile\n')
		f.write('{\n')
		f.write('    format      ascii;\n')
		f.write('    class       volVectorField;\n')
		f.write(f'    location    "{time_fld:.3f}";\n')
		f.write(f'    object      {name_field};\n')
		f.write('}\n')
		f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
		f.write('\ndimensions      [0 1 -1 0 0 0 0];\n')
		f.write('\ninternalField   nonuniform List<vector>\n')
		f.write(f'{num_elem}\n')
		f.write('(\n')
		# Write the data
		for ii in range(num_elem):
			f.write(f"({varr[ii,0]:.8e} {varr[ii,1]:.8e} {varr[ii,2]:.8e})\n")
		f.write(');\n')
		# Write the boundary conditions
		f.write('\nboundaryField\n')
		f.write('{\n')
		f.write('\tBOTTOM\n')
		f.write('\t{\n')
		f.write('\t\ttype            symmetryPlane;\n')
		f.write('\t}\n')
		f.write('\tTOP\n')
		f.write('\t{\n')
		f.write('\t\ttype            symmetryPlane;\n')
		f.write('\t}\n')
		f.write('\tOUTLET\n')
		f.write('\t{\n')
		f.write('\t\ttype            symmetryPlane;\n')
		f.write('\t}\n')
		f.write('\tINLET\n')
		f.write('\t{\n')
		f.write('\t\ttype            fixedValue;\n')
		f.write('\t\tvalue           uniform (1 0 0);\n')
		f.write('\t}\n')

		f.write('\tFrontAndBack\n')
		f.write('\t{\n')
		f.write('\t\ttype            empty;\n')
		f.write('\t}\n')

		f.write('\tbaffleFacesSquare_master\n')
		f.write('\t{\n')
		f.write('\t\ttype            cyclic;\n')
		f.write('\t}\n')

		f.write('\tbaffleFacesSquare_slave\n')
		f.write('\t{\n')
		f.write('\t\ttype            cyclic;\n')
		f.write('\t}\n')
		f.write('}\n')

		f.write('// ************************************************************************* //')


	return

#==============================================================


#========================= INPUTS =============================

# Name of the case
case_nm     = "DA_IBM4_full_dom_newVer-hyperLocalization-firstSteps"

# Which member to work on?
memberID = 1

# Time
time_fld   = 300.016

# Field
field_label = "DA_Updt"

#==============================================================

flnm_U     = f"UMat{time_fld:.6f}"
flnm_U_upt = f"UMat_upt{time_fld:.6f}"

dir_in  = (f"/home/sarp/square_cylinder/{case_nm}/results")

dir_out = (f"/home/sarp/square_cylinder/{case_nm}/")


# Create directory for output figures ------------
dir_out = f'{dir_out}figs_postProc'
try:
    os.mkdir(dir_out)
except OSError as error:
    print(f'\n {dir_out} already exists..')
#-------------------------------------------------

print(f"\n Reading the velocity field from file '{dir_in}/{flnm_U}'")
# Read the data as NumPy array
data = np.loadtxt(f"{dir_in}/{flnm_U}")[:, (memberID-1)]

# There are num_elem*3 rows (U, V, W fields)
num_elem = int(data.shape[0]/3)

# Separate U, V and W data
u_fld = np.zeros((num_elem, 3))  # Initialize array for the three components of the velocity
u_fld[:,0] = data[0:num_elem]
u_fld[:,1] = data[num_elem:2*num_elem]
u_fld[:,2] = data[2*num_elem:3*num_elem]


#num = 10
#print(f"\n First {num} U, V, W values :")
#for ii in range(num):
#	print(f"{u_fld[ii,0]:.9f}, {u_fld[ii,1]:.9f}, {u_fld[ii,2]:.9f}")


print(f"\n Reading the Updated velocity field from file '{dir_in}/{flnm_U}'")
# Read the data as NumPy array
data = np.loadtxt(f"{dir_in}/{flnm_U_upt}")[:, (memberID-1)]

# There are num_elem*3 rows (U, V, W fields)
num_elem = int(data.shape[0]/3)

# Separate U, V and W data
u_fld_upt = np.zeros((num_elem, 3))  # Initialize array for the three components of the velocity
u_fld_upt[:,0] = data[0:num_elem]
u_fld_upt[:,1] = data[num_elem:2*num_elem]
u_fld_upt[:,2] = data[2*num_elem:3*num_elem]

# Compute the update term K*(x_f - H.y)
upt_term = u_fld_upt - u_fld

# Writing the output file
print(f"\n Writing the output file in the directory {dir_out}/.")
writeVectorFieldOpenFoam(upt_term, field_label, dir_out)


# Finish Program
print("\n Program finished!")
