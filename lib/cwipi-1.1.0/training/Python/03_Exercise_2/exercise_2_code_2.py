#!/usr/bin/env python

import mpi4py.MPI as MPI

# Initialize MPI
comm = MPI.COMM_WORLD
i_rank = comm.rank
n_rank = comm.size

import numpy as np
import math

# pycwp
try:
    from pycwp import pycwp
    if i_rank == 0:
        print("Yes, we have found the module :D")
except:
    if i_rank == 0:
        print("Oh no, couldn't find the module :'(")

# pypdm
import Pypdm.Pypdm as PDM

n_code = 1
code_name = ["code2"]
is_active_rank = True
intra_comm = pycwp.init(comm,
                        code_name,
                        is_active_rank)

coupled_code_name = ["code1"]
n_part = 1
cpl = pycwp.Coupling(code_name[0],
                     "coupling",
                     coupled_code_name[0],
                     pycwp.INTERFACE_SURFACE,
                     pycwp.COMM_PAR_WITH_PART,
                     pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                     n_part,
                     pycwp.DYNAMIC_MESH_DEFORMABLE,
                     pycwp.TIME_EXCH_USER_CONTROLLED)

cpl.visu_set(1,
             pycwp.VISU_FORMAT_ENSIGHT,
             "text")

# Create mesh :
mesh = PDM.generate_mesh_rectangle_simplified(intra_comm[0], 10)

n_vtx       = mesh["n_vtx"]
n_elt       = mesh["n_elt"]
coords      = mesh["coords"]
elt_vtx_idx = mesh["elt_vtx_idx"]
elt_vtx     = mesh["elt_vtx"]

cpl.mesh_interf_vtx_set(0,
                      coords,
                      None)

block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

cpl.mesh_interf_f_poly_block_set(0,
                               block_id,
                               elt_vtx_idx,
                               elt_vtx,
                               None)

cpl.mesh_interf_finalize()

n_components = 1
field = cpl.field_create("a super fancy field",
                          pycwp.DOUBLE,
                          pycwp.FIELD_STORAGE_INTERLACED,
                          n_components,
                          pycwp.DOF_LOCATION_NODE,
                          pycwp.FIELD_EXCH_SEND,
                          pycwp.STATUS_ON)

n_vtx = len(coords)//3
send_field_data = np.arange(n_vtx*n_components, dtype=np.double)

for i in range(n_vtx):
  send_field_data[i] = coords[3*i]

field.data_set(0,
               pycwp.FIELD_MAP_SOURCE,
               send_field_data)

cpl.spatial_interp_property_set("tolerance",
                                pycwp.DOUBLE,
                                "0.001")

itdeb = 1
itend = 10
ttime = 0.0
dt    = 0.1

degrad = math.acos(-1.0)/180.
x = 0.0
y = 0.0
alpha = 2
alpha = alpha * degrad
sina = math.sin(alpha)
cosa = math.cos(alpha)

for it in range(itdeb, itend+1):

  ttime = (it-itdeb)*dt;

  pycwp.time_step_beg(code_name[0],
                      ttime);

  cpl.spatial_interp_weights_compute()

  field.issend()

  field.wait_issend()

  pycwp.time_step_end(code_name[0])


# Delete field

cpl.mesh_interf_del()

# Delete the coupling

pycwp.finalize()

MPI.Finalize()
