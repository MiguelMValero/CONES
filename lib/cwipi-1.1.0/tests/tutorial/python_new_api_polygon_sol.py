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
    # Here 2 codes are coupled. code1 runs on the processor of
    # MPI rank 0 and code2 runs on the processor if MPI rank 1.
    # In this version of CWIPI several codes can execute on the
    # same MPI rank (here only one code per processor, so n_code = 1).
    # Therefore, an array of code names is given at initialization.
    # is_active_rank tells if current ranks will be used
    # in the CWIPI coupling computations.
    # intra_comm is an array of MPI communicators
    # giving the for each code on the processors the communicator
    # to communicate through the ranks of that code.
    n_code = 1

    if (i_rank == 0):
        code_name = ["code1"]

    if (i_rank == 1):
        code_name = ["code2"]

    is_active_rank = True

    intra_comm = pycwp.init(comm,
                            code_name,
                            is_active_rank)

    # Create the coupling :
    # One CWIPI context can hold several couplings. Let us set up the
    # coupling between code1 and code2. CWP_INTERFACE_SURFACE informs
    # that the geometrical interface of the meshes of the coupled
    # codes is a surface, still for CWIPI the coordinate system is 3D.
    # CWP_COMM_PAR_WITH_PART means that each mesh is partitionned
    # over the processors of its code. Here the mesh does not change
    # over the coupling, so CWP_DYNAMIC_MESH_STATIC is set.
    # CWP_TIME_EXCH_USER_CONTROLLED is not used yet.
    if (i_rank == 0):
        coupled_code_name = ["code2"]
    if (i_rank == 1):
        coupled_code_name = ["code1"]
    n_part = 1
    cpl = pycwp.Coupling(code_name[0],
                         "code1_code2",
                         coupled_code_name[0],
                         pycwp.INTERFACE_SURFACE,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         n_part,
                         pycwp.DYNAMIC_MESH_STATIC,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    # Set coupling visualisation:
    # Output files of the code fields will be written in Ensight ASCII
    # format (easily readable by Paraview) at each iteration.
    cpl.visu_set(1,
                 pycwp.VISU_FORMAT_ENSIGHT,
                 "text")

    # Create the field
    # It is possible to operate a bidirectional exchange (see c_new_vs_old_sendrecv).
    # For sake of simplicity, this example will only send the field
    # of code1 (CWP_FIELD_EXCH_SEND) to code2 (CWP_FIELD_EXCH_RECV).
    # On code1 there is a field (CWP_FIELD_MAP_SOURCE) located at
    # the vertices (CWP_DOF_LOCATION_NODE) with one component (n_components)
    # which is the x coordinate of the mesh in this test.
    n_components = 1

    # for code1
    if (i_rank == 0):
      field = cpl.field_create("a super fancy field",
                               pycwp.DOUBLE,
                               pycwp.FIELD_STORAGE_INTERLACED,
                               n_components,
                               pycwp.DOF_LOCATION_NODE,
                               pycwp.FIELD_EXCH_SEND,
                               pycwp.STATUS_ON)

    # for code2
    if (i_rank == 1):
      field = cpl.field_create("a super fancy field",
                               pycwp.DOUBLE,
                               pycwp.FIELD_STORAGE_INTERLACED,
                               n_components,
                               pycwp.DOF_LOCATION_NODE,
                               pycwp.FIELD_EXCH_RECV,
                               pycwp.STATUS_ON)

    # Begin time step :
    # In this example there is only one time step. It is mandatory to create the
    # coupling and the associated fields before starting the first time step.
    pycwp.time_step_beg(code_name[0],
                        0.0);

    # Set the mesh vertices coordinates :
    # The coordinate system in CWIPI is always 3D, so
    # we allocate an array of the time the number of vertices
    # (11 here) to set the coordinates in. The coordinates are
    # interlaced (x0, y0, z0, x1, y1, z1, ..., xn, yn, zn).
    # The None argument will be explained later.
    coords = np.array([0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0, \
              3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0], dtype=np.double)
    cpl.mesh_interf_vtx_set(0,
                            coords,
                            None)

    # Set the mesh polygons connectivity :
    # Let us set a mesh of 5 polygons (CWP_BLOCK_FACE_POLY).
    # An index array (connec_idx) of size n_elts+1 contains the
    # information of the number of vertices per polygon. The first
    # index is always 0, from there we add up the number of vertices
    # per element. Here one triangle, 2 quadrangles and 2 pentagons.
    # The connectivity between elements and vertices is an array of
    # size connec_idx(n_elts+1) (here 21).
    block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

    connec_idx = np.array([0,3,7,11,16,21], dtype=np.int32)
    connec = np.array([1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8], dtype=np.int32)
    cpl.mesh_interf_f_poly_block_set(0,
                                     block_id,
                                     connec_idx,
                                     connec,
                                     None)

    # Finalize mesh :
    # CWIPI hides the parallelism for users, so it is not
    # mandatory to give a global numbering for mesh data (the
    # None arguments earlier). If not given this numbering is
    # generated by CWIPI by the following function,
    # as well as the underlying mesh data structure
    cpl.mesh_interf_finalize()

    # Set the field values :
    # Note that the user has to allocate the array for the
    # field that will be received by code2 (CWP_FIELD_MAP_TARGET).
    n_vtx = len(coords)//3
    send_field_data = np.arange(n_vtx*n_components, dtype=np.double)
    recv_field_data = np.arange(n_vtx*n_components, dtype=np.double)

    for i in range(n_vtx):
      send_field_data[i] = coords[3*i]

    # for code1
    if (i_rank == 0):
      field.data_set(0,
                     pycwp.FIELD_MAP_SOURCE,
                     send_field_data)

    # for code2
    if (i_rank == 1):
      field.data_set(0,
                     pycwp.FIELD_MAP_TARGET,
                     recv_field_data)

    # Compute interpolation weights :
    # Set a geometric tolerance of 10% of an element size for
    # point localisation.
    cpl.spatial_interp_property_set("tolerance",
                                    pycwp.DOUBLE,
                                    "0.1")

    cpl.spatial_interp_weights_compute()

    # Exchange field values between codes :
    # The field exchange functions mimic the way the associated
    # MPI functions work, see MPI documentation for more information.

    # for code1
    if (i_rank == 0):
      field.issend()

    # for code2
    if (i_rank == 1):
      field.irecv()

    # for code1
    if (i_rank == 0):
      field.wait_issend()

    # for code2
    if (i_rank == 1):
      field.wait_irecv()

    # Check interpolation :
    # These functions allow to know how many and for which target
    # vertices the interpolation operation has been unsuccessful.
    if (i_rank == 1):
      n_uncomputed_tgts = field.n_uncomputed_tgts_get(0);
      uncomputed_tgts   = field.uncomputed_tgts_get(0);

    # End time step :
    pycwp.time_step_end(code_name[0])

    # Delete field :
    # del of the class Field is automatically called once there
    # are no more references to the field instance in cpl

    # Delete Mesh :
    cpl.mesh_interf_del()

    # Delete the coupling :
    # del of the class Coupling is automatically called once there
    # are no more references to cpl

    # Finalize CWIPI :
    pycwp.finalize()

    # Finalize MPI :
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
