#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2023  ONERA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

import mpi4py.MPI as MPI
import numpy as np
import sys

def runTest():
    """
    Run tests on Python interface of new API
    """

    # Initialize MPI
    comm = MPI.COMM_WORLD
    i_rank = comm.rank
    n_rank = comm.size

    # This test mimics the coupling between 2 codes running each
    # on one processor.
    if (n_rank != 2):
        if i_rank == 0:
            print("      Not executed : only available for 2 processes")
        return

    # Load Python CWIPI module :
    try:
        from pycwp import pycwp
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # Initialize CWIPI :
    # Use CWP_Init for code1 running on MPI rank 0 and code2
    # running on MPU rank 1.
    # ------------------------------------------------------- To fill in
    n_code = 1

    if i_rank == 0:
        code_name = ["code1"]
        I_am_code1 = True
        I_am_code2 = False

    if i_rank == 1:
        code_name = ["code2"]
        I_am_code1 = False
        I_am_code2 = True

    # ---------------------------------------------------- End To fill in

    # Create the coupling :
    # Use CWP_Cpl_create to couple code1 with code2 on a surface
    # interface. Operate the localization with the octree method.
    # In this tutorial, one coupling iteration is done on a mesh
    # partitionned over all processors. This induces that the mesh
    # is static.
    # ------------------------------------------------------- To fill in
    if I_am_code1:
        coupled_code_name = ["code2"]
    if I_am_code2:
        coupled_code_name = ["code1"]
    n_part = 1

    # ---------------------------------------------------- End To fill in

    # Set coupling visualisation:
    # Use CWP_Visu_set to output ASCII Ensight format files
    # at each iteration step.
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Set the mesh vertices coordinates :
    # Use CWP_Mesh_interf_vtx_set to set the mesh vertex coordinates,
    # no global numbering of the vertices will be given. In this
    # simple setting, there is only one partition per processor.
    # ------------------------------------------------------- To fill in
    coords = np.array([0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0, \
              3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0], dtype=np.double)

    # ---------------------------------------------------- End To fill in

    # Set the mesh polygons connectivity :
    # Use CWP_Mesh_interf_block_add to create a block of
    # of polygons. Choose the correct CWIPI function
    # to set a polygonal mesh, no need to give the elements
    # global numbering.
    # ------------------------------------------------------- To fill in

    connec_idx = np.array([0,3,7,11,16,21], dtype=np.int32)
    connec = np.array([1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8], dtype=np.int32)

    # ---------------------------------------------------- End To fill in

    # Finalize mesh :
    # Use the correct CWIPI function to generate the
    # mesh global numbering and the underlying mesh data structure
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Create and set the field values :
    # Use CWP_Field_create and CWP_Field_data_set to create and set
    # a field onto the mesh. code1 will send its field which code2
    # will receive. The field is located at the mesh nodes.
    # There is only one mesh partition in this tutorial. Activate
    # visualization for this field if you wish it to be in the
    # Ensight file.
    # Do not forget to begin the time step AFTER creating the fields,
    # but BEFORE setting the fields data!
    # ------------------------------------------------------- To fill in
    n_components = 1

    n_vtx = len(coords)//3

    send_field_data = np.arange(n_vtx*n_components, dtype=np.double)
    recv_field_data = np.arange(n_vtx*n_components, dtype=np.double)

    for i in range(n_vtx):
      send_field_data[i] = coords[3*i]

    # Create the field :
    # ...

    # Begin the time step :
    # ...

    # Set the field pointers :
    # for code1
    # if I_am_code1: # to uncomment

    # for code2
    # if I_am_code2: # to uncomment

    # ---------------------------------------------------- End To fill in

    # Compute interpolation weights :
    # Choose the two correct CWIPI functions to set the geometric
    # tolerance to 10% of an element size for point localisation
    # and to compute the interpolation weights.
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Exchange field values between codes :
    # Use the CWIPI exchange functions similar to the MPI ones
    # for code1 to send a field and code2 to receive that field.
    # The reason for the double if statements is that the exchange
    # operation is non-blocking. That means that work can be done
    # before calling the wait function where the field need to be
    # used or changed.
    # ------------------------------------------------------- To fill in

    # for code1
    # if I_am_code1: # to uncomment

    # for code2
    # if I_am_code2: # to uncomment

    # for code1
    # if Iam_code1: # to uncomment

    # for code2
    # if I_am_code2: # to uncomment

    # ---------------------------------------------------- End To fill in

    # End time step :
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Check interpolation :
    # For the receiving code, check the vertices for which the
    # interpolation has been unsuccessful.
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Delete field :
    # Is there something to do in Python?
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Delete Mesh :
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Delete the coupling :
    # Is there something to do in Python?
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Finalize CWIPI :
    # ------------------------------------------------------- To fill in

    # ---------------------------------------------------- End To fill in

    # Finalize MPI :
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
