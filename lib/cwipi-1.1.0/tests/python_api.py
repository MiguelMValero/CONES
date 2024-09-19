#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library. 
#
# Copyright (C) 2011  ONERA
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

f=None

def userInterp(entities_dim,
                 n_local_vertex,
                 n_local_element,
                 n_local_polhyedra,
                 n_distant_point,
                 local_coordinates,
                 local_connectivity_index,
                 local_connectivity,
                 local_polyhedra_face_index,
                 local_polyhedra_cell_to_face,
                 local_polyhedra_face_connectivity_index,
                 local_polyhedra_face_connectivity,
                 distant_points_coordinates,
                 distant_points_location,
                 distant_points_distance,
                 distant_points_barycentric_coordinates_index,
                 distant_points_barycentric_coordinates,
                 stride,
                 solver_type,
                 local_field,
                 distant_field):
    """
    User interpolation method
    """
    global f

    f.write('user interpolation')

    if (distant_field is not None):
        distant_field = 1.234
            

def runTest():
    """
    Run Python API test
    """
    global f
    commloc=MPI.Comm()
    commworld = MPI.COMM_WORLD

    rank = commworld.rank
    size = commworld.size

    if (rank == 0):
        print("\nSTART: python_api.py")

    if (size != 2):
        if rank == 0:
            print("      Not executed : only available for 2 processus")
        return

    applis = ["proc0","proc1"]
    fname = ["proc0","proc1"]

    try:
        from cwipi import cwipi
    except:
        if rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    #
    # Define a python file to write CWIPI outputs 

    if (rank == 0):
        print("        Output redirection")

    srank = '{0}'.format(rank)
    f=open("python_api_"+srank.zfill(4)+".txt",'w')
    cwipi.set_output_listing(f)
    comm_loc = cwipi.init(MPI.COMM_WORLD, applis[rank])

    #
    # Control parameters

    if (rank == 0):
        print("        Control parameters")

    if (rank == 0):
        cwipi.add_local_int_control_parameter("i_0_param_1", 1)
        cwipi.add_local_int_control_parameter("i_1_param_1", 1)
        cwipi.add_local_int_control_parameter("i_2_param_1", 1)
        cwipi.add_local_int_control_parameter("i_3_param_1", 1)
        cwipi.add_local_double_control_parameter("d_0_param_1", 1)
        cwipi.add_local_double_control_parameter("d_1_param_1", 1)
        cwipi.add_local_double_control_parameter("d_2_param_1", 1)
        cwipi.add_local_string_control_parameter("s_0_param_1", "1")
        cwipi.add_local_string_control_parameter("s_1_param_1",  "1")
        cwipi.add_local_string_control_parameter("s_2_param_1",  "1")
        cwipi.add_local_string_control_parameter("s_3_param_1",  "1")
        cwipi.add_local_string_control_parameter("s_4_param_1",  "1")
        cwipi.add_local_string_control_parameter("s_5_param_1",  "1")
    else :
        cwipi.add_local_int_control_parameter("i_0_param_2", 1)
        cwipi.add_local_int_control_parameter("i_1_param_2", 1)
        cwipi.add_local_int_control_parameter("i_2_param_2", 1)
        cwipi.add_local_int_control_parameter("i_3_param_2", 1)
        cwipi.add_local_int_control_parameter("i_4_param_2", 1)
        cwipi.add_local_double_control_parameter("d_0_param_2", 1)
        cwipi.add_local_double_control_parameter("d_1_param_2", 1)
        cwipi.add_local_string_control_parameter("s_0_param_2",  "1")
        cwipi.add_local_string_control_parameter("s_1_param_2",  "1")
        cwipi.add_local_string_control_parameter("s_2_param_2",  "1")
        cwipi.add_local_string_control_parameter("s_3_param_2",  "1")

    cwipi.synchronize_control_parameter(applis[(rank + 1) % 2])

    has_toto_0 = cwipi.has_int_parameter("proc0","toto")
    has_toto_1 = cwipi.has_int_parameter("proc1","toto")
    has_param1_0 = cwipi.has_double_parameter("proc0","d_0_param_1")
    has_param1_1 = cwipi.has_double_parameter("proc1","d_0_param_2")

#    f.write(cwipi.get_list_int_parameter("proc0"))
    f.write("  - has d_0_param_1 proc0 : {param}\n".format(param=has_param1_0))
    f.write("  - has d_0_param_2 proc1 : {param}\n".format(param=has_param1_1))
    f.write("  - list int param proc0 : {param}\n".format(param=cwipi.get_list_int_parameter("proc0")))
    f.write("  - list int param proc1 : {param}\n".format(param=cwipi.get_list_int_parameter("proc1")))
    f.write("  - list double param proc0 : {param}\n".format(param=cwipi.get_list_double_parameter("proc0")))
    f.write("  - list double param proc1 : {param}\n".format(param=cwipi.get_list_double_parameter("proc1")))
    f.write("  - list string param proc0 : {param}\n".format(param=cwipi.get_list_string_parameter("proc0")))
    f.write("  - list string param proc1 : {param}\n".format(param=cwipi.get_list_string_parameter("proc1")))

    f.flush()  # make sure f is flushed before cwipi writes its properties
    cwipi.dump_application_properties()

#    param_1 = cwipi.get_local_int_control_parameter("param_1")
#    param_2 = cwipi.get_distant_int_control_parameter(applis[(rank + 1) % 2],"param_1")
#    f.write('parametres: {param_1}, {param_2}'.format(param_1=param_1, param_2=param_2))

    #
    # Class coupling

    if (rank == 0):
        print("        Create coupling")

    # Constructor

    cpl = cwipi.Coupling("cpl",
                         cwipi.COUPLING_PARALLEL_WITH_PARTITIONING,
                         applis[(rank + 1) % 2] , 
                         2, 
                         0.1,  
                         cwipi.CYCLIC_MESH,
                         cwipi.SOLVER_CELL_VERTEX,
                         1, 
                         "Ensight", 
                         "txt")
                         

    # Mesh

    if (rank == 0):
        print("        Create mesh")

    coord = np.array([-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0], dtype=np.double)
    connec_idx = np.array([0, 4], dtype=np.int32)
    connec = np.array([1, 2, 3, 4], dtype=np.int32)

    cpl.define_mesh(4, 1, coord, connec_idx, connec)

    # Only send 

    if (rank == 0):
        print("        Exchange Proc 0 -> Proc 1")

    f.write("send :\n")

    sendField=np.array([0.1, 0.2, 0.3, 0.4], dtype=np.double)
    recvField=np.arange(4, dtype=np.double)

    if rank == 0:
        result = cpl.exchange("ech1", 
                              1, 
                              1, 
                              0.1, 
                              "field_s", sendField, 
                              "field_r", None)
    else:
        result = cpl.exchange("ech1", 
                              1, 
                              1, 
                              0.1, 
                              "field_s", None, 
                              "field_r", recvField)

    f.write('  - status : {param_1}\n'.format(param_1=result["status"]))
    if rank == 1:
        f.write("  - number of not located points : {param}\n".format(param=result["n_not_located_points"]))


    # Send and receive with user interpolation

    f.write("send receive :")

    if (rank == 0):
        print("        Define user interpolation")

    cpl.set_interpolation_function(userInterp) # User interpolation activation

    if (rank == 0):
        print("        Exchange Proc 0 <-> Proc 1 with user interpolation")

    result = cpl.exchange("ech2", 
                          1, 
                          1, 
                          0.1, 
                          "field_s2", sendField, 
                          "field_r2", recvField)

    f.write("  - status : {param_1}\n".format(param_1=result["status"]))
    f.write("  - number of not located points : {param}\n".format(param=result["n_not_located_points"]))

    # Delete coupling object

    if (rank == 0):
        print("        Delete coupling")

    #del cpl

    f.write("\nThe end!\n")

if __name__ == '__main__':
    runTest()

