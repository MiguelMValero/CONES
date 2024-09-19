#!/usr/bin/env python

import mpi4py.MPI as MPI

# Initialize MPI
comm = MPI.COMM_WORLD
i_rank = comm.rank
n_rank = comm.size

import numpy as np

# pycwp
try:
    from pycwp import pycwp
    if i_rank == 0:
        print("Yes, we have found the module :D")
except:
    if i_rank == 0:
        print("Oh no, couldn't find the module :'(")

# parse command line arguments to know wether to activate the bonus

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-b",
                    "--bonus",
                    help="Activate bonus",
                    action="store_true")

args = parser.parse_args()

bonus = args.bonus

n_code = 1
code_name = ["code2"]
is_active_rank = True
intra_comm = pycwp.init(comm,
                        code_name,
                        is_active_rank)

coupled_code_name = ["code1"]
n_part = 1

spatial_interp_algorithm = pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE
if (bonus):
    spatial_interp_algorithm = pycwp.SPATIAL_INTERP_FROM_INTERSECTION

cpl = pycwp.Coupling(code_name[0],
                     "code1_code2",
                     coupled_code_name[0],
                     pycwp.INTERFACE_SURFACE,
                     pycwp.COMM_PAR_WITH_PART,
                     spatial_interp_algorithm,
                     n_part,
                     pycwp.DYNAMIC_MESH_STATIC,
                     pycwp.TIME_EXCH_USER_CONTROLLED)

cpl.visu_set(1,
             pycwp.VISU_FORMAT_ENSIGHT,
             "text")

coords = np.array([0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0, \
          3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0], dtype=np.double)
cpl.mesh_interf_vtx_set(0,
                        coords,
                        None)

block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
connec_idx = np.array([0,3,7,11,16,21], dtype=np.int32)
connec = np.array([1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8], dtype=np.int32)
cpl.mesh_interf_f_poly_block_set(0,
                                 block_id,
                                 connec_idx,
                                 connec,
                                 None)

cpl.mesh_interf_finalize()

n_components = 1

location = pycwp.DOF_LOCATION_NODE
if (bonus):
  location = pycwp.DOF_LOCATION_CELL_CENTER

field = cpl.field_create("a super fancy field",
                           pycwp.DOUBLE,
                           pycwp.FIELD_STORAGE_INTERLACED,
                           n_components,
                           location,
                           pycwp.FIELD_EXCH_RECV,
                           pycwp.STATUS_ON)

n_vtx = len(coords)//3
size = n_vtx*n_components
if (bonus):
  n_elt = len(connec_idx)-1
  size = n_elt*n_components

recv_field_data = np.arange(size, dtype=np.double)

field.data_set(0,
               pycwp.FIELD_MAP_TARGET,
               recv_field_data)

pycwp.time_step_beg(code_name[0],
                    0.0);

cpl.spatial_interp_property_set("tolerance",
                                pycwp.DOUBLE,
                                "0.001")

cpl.spatial_interp_weights_compute()

field.irecv()

field.wait_irecv()

pycwp.time_step_end(code_name[0])

# Delete field

cpl.mesh_interf_del()

# Delete the coupling

pycwp.finalize()

MPI.Finalize()
