#!/usr/bin/env python
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2022-2023  ONERA
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
try:
  from pycwp import pycwp
except:
  print("Error : CWIPI module not found (update PYTHONPATH variable)")
  sys.exit(1)

def ijk_ordering(order):
  ijk = np.zeros(2*(order+1)**2, dtype=np.int32)
  ivtx = 0
  for j in range(order+1):
    for i in range(order+1):
      ijk[2*ivtx  ] = i
      ijk[2*ivtx+1] = j
      ivtx += 1
  return ijk

def generate_mesh():

  coord = np.zeros((9+16)*3, dtype=np.double)
  ivtx = 0
  for j in range(3):
    for i in range(3):
      coord[3*ivtx  ] = i/2
      coord[3*ivtx+1] = j/2
      ivtx += 1

  for j in range(4):
    for i in range(4):
      coord[3*ivtx  ] = i/3+1
      coord[3*ivtx+1] = j/3
      ivtx += 1

  blocks = {
    2 :     1+np.arange(9,  dtype=np.int32),
    3 : 9 + 1+np.arange(16, dtype=np.int32),
  }

  return {
    "coord" : coord,
    "blocks": blocks
  }


def run_test():
  comm = MPI.COMM_WORLD

  i_rank = comm.rank
  n_rank = comm.size

  assert(n_rank == 2)

  code_names = [f"code{i}" for i in range(2)]

  code_name = code_names[i_rank]


  # Initialize CWIPI
  n_code = 1
  is_active_rank = True
  intra_comm = pycwp.init(comm,
                          [code_name],
                          is_active_rank)

  # Create coupling
  cpl = pycwp.Coupling(code_name,
                       "python_new_api_multiblock_2d",
                       code_names[(i_rank+1)%2],
                       pycwp.INTERFACE_SURFACE,
                       pycwp.COMM_PAR_WITH_PART,
                       pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                       # pycwp.SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,
                       1,
                       pycwp.DYNAMIC_MESH_VARIABLE,
                       pycwp.TIME_EXCH_USER_CONTROLLED)

  # Visu does not work for high-order meshes
  # cpl.visu_set(1,
  #              pycwp.VISU_FORMAT_ENSIGHT,
  #              "text")


  # Generate mesh
  mesh = generate_mesh()

  print(mesh)

  # Set vertices
  cpl.mesh_interf_vtx_set(0,
                          mesh["coord"],
                          None)

  n_vtx = len(mesh["coord"])//3

  # Set blocks (FACE_QUADHO order = 2, 3)
  # block_ids = []
  ijk_grid = {}

  for order in mesh["blocks"]:
    block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_QUADHO)
    # block_ids.append(block_id)

    cpl.mesh_interf_block_ho_set(0,
                                 block_id,
                                 order,
                                 mesh["blocks"][order],
                                 None)

    ijk_grid[order] = ijk_ordering(order)
    cpl.mesh_interf_ho_ordering_from_IJK_set(pycwp.BLOCK_FACE_QUADHO,
                                             order,
                                             ijk_grid[order])
    # break

  # Finalize mesh
  cpl.mesh_interf_finalize()


  # Create field
  if i_rank == 0:
    field = cpl.field_create("field",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_SEND,
                             pycwp.STATUS_ON)

  else:
    field = cpl.field_create("field",
                             pycwp.DOUBLE,
                             pycwp.FIELD_STORAGE_INTERLACED,
                             1,
                             pycwp.DOF_LOCATION_NODE,
                             pycwp.FIELD_EXCH_RECV,
                             pycwp.STATUS_ON)

  pycwp.time_step_beg(code_name, 0.)

  if i_rank == 0:
    send_field = np.zeros(n_vtx, np.double)
    for i in range(n_vtx):
      send_field[i] = mesh["coord"][3*i]

    field.data_set(0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_field)
  else:
    recv_field = np.zeros(n_vtx, np.double)
    field.data_set(0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_field)


  # Spatial interpolation
  cpl.spatial_interp_property_set("tolerance", pycwp.DOUBLE, "1e-2")
  cpl.spatial_interp_weights_compute()


  # Exchange interpolated field
  if i_rank == 0:
    field.issend()
  else:
    field.irecv()

  if i_rank == 0:
    field.wait_issend()
  else:
    field.wait_irecv()

    computed_tgt = field.computed_tgts_get(0)

    # print("check recv_field :")
    for i, j in enumerate(computed_tgt):
      if abs(recv_field[i] - mesh['coord'][3*(j-1)]) > 1e-9:
        print(f"!! node {j-1} expected {mesh['coord'][3*(j-1)]} but recieved {recv_field[i]}")

  pycwp.time_step_end(code_name)

  # Finalize
  pycwp.finalize()

  if i_rank == 0:
    print("End", flush=True)

  MPI.Finalize()

if __name__ == '__main__':
    run_test()





