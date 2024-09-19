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

def generate_mesh_multiblock():

  coord = np.zeros((16*3), dtype=np.double)
  ivtx = 0
  for k in range(2):
    for j in range(2):
      for i in range(4):
        coord[3*ivtx  ] = i
        coord[3*ivtx+1] = j
        coord[3*ivtx+2] = k
        if i == 2:
          coord[3*ivtx+1] = 0.5 + 1.5*(coord[3*ivtx+1] - 0.5)
          coord[3*ivtx+2] = 0.5 + 1.5*(coord[3*ivtx+2] - 0.5)
        ivtx += 1

  blocks = {
    pycwp.BLOCK_CELL_HEXA8: {
      "connec": np.array([1, 2, 6, 5, 9, 10, 14, 13], dtype=np.int32),
    },
    pycwp.BLOCK_CELL_POLY: {
      "cell_face_idx": np.array([0, 10], dtype=np.int32),
      "cell_face"    : 1+np.arange(10, dtype=np.int32),
      "face_vtx_idx" : 4*np.arange(11, dtype=np.int32),
      "face_vtx"     : np.array([
                                2, 3, 11, 10,
                                3, 4, 12, 11,
                                6, 14, 15, 7,
                                7, 15, 16, 8,
                                2, 10, 14, 6,
                                4, 8, 16, 12,
                                2, 6, 7, 3,
                                3, 7, 8, 4,
                                10, 11, 15, 14,
                                11, 12, 16, 15
                                ], dtype=np.int32)
    }
  }

  return {
    "coord"  : coord,
    "blocks" : blocks
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
                       "python_new_api_multiblock_3d",
                       code_names[(i_rank+1)%2],
                       pycwp.INTERFACE_VOLUME,
                       pycwp.COMM_PAR_WITH_PART,
                       # pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                       pycwp.SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,
                       1,
                       pycwp.DYNAMIC_MESH_VARIABLE,
                       pycwp.TIME_EXCH_USER_CONTROLLED)

  cpl.visu_set(1,
               pycwp.VISU_FORMAT_ENSIGHT,
               "text")


  # Generate mesh
  mesh = generate_mesh_multiblock()

  print(mesh)

  # Set vertices
  cpl.mesh_interf_vtx_set(0,
                          mesh["coord"],
                          None)

  n_vtx = len(mesh["coord"])//3

  # Set blocks (CELL_HEXA_, CELL_POLY)
  # block_ids = []
  for elt_type in mesh["blocks"]:
    block_id = cpl.mesh_interf_block_add(elt_type)
    # block_ids.append(block_id)

    if elt_type == pycwp.BLOCK_CELL_POLY:
      # Generic polyhedra
      cpl.mesh_interf_c_poly_block_set(0,
                                       block_id,
                                       mesh["blocks"][elt_type]["face_vtx_idx"],
                                       mesh["blocks"][elt_type]["face_vtx"],
                                       mesh["blocks"][elt_type]["cell_face_idx"],
                                       mesh["blocks"][elt_type]["cell_face"],
                                       None)
    else:
      # Standard element (fixed number of vertices per element, no need for idx)
      cpl.mesh_interf_block_std_set(0,
                                    block_id,
                                    mesh["blocks"][elt_type]["connec"],
                                    None)

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

  pycwp.time_step_end(code_name)

  # Finalize
  pycwp.finalize()

  if i_rank == 0:
    print("End", flush=True)

  MPI.Finalize()

if __name__ == '__main__':
    run_test()





