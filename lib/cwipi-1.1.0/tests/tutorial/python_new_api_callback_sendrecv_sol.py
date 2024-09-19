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

def my_interpolation(field,
                     i_part,
                     buffer_in,
                     buffer_out):

  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  n_comp = field.n_components

  spatial_interp_algorithm = field.spatial_interp_algo

  if spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES:
    tgt_data       = field.tgt_data_properties_get(i_part)
    n_tgt          = tgt_data["n_tgt"]
    ref_tgt        = tgt_data["computed_tgt"]
    tgt_to_src_idx = tgt_data["tgt_to_src_idx"]

    distance2 = field.nearest_neighbors_distances_get(i_part)

    for i, jtgt in enumerate(ref_tgt):
      itgt = jtgt-1
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
      sum_w = 0
      for isrc in range(tgt_to_src_idx[i], tgt_to_src_idx[i+1]):
        w = 1./max(1.e-24, distance2[isrc])
        sum_w += w
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
        w*buffer_in[n_comp*isrc:n_comp*(isrc+1)]

      buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
      buffer_out[n_comp*itgt:n_comp*(itgt+1)] / sum_w

  elif spatial_interp_algorithm == pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
    src_data       = field.src_data_properties_get(i_part)
    n_src          = src_data["n_src"]
    src_to_tgt_idx = src_data["src_to_tgt_idx"]

    weight = field.location_weights_get(i_part)

    cell_data    = field.location_internal_cell_vtx_get(i_part)
    cell_vtx_idx = cell_data["cell_vtx_idx"]
    cell_vtx     = cell_data["cell_vtx"]

    idx = 0
    for isrc in range(n_src):
      for itgt in range(src_to_tgt_idx[isrc], src_to_tgt_idx[isrc+1]):
        buffer_out[n_comp*itgt:n_comp*(itgt+1)] = 0
        for i in range(cell_vtx_idx[isrc], cell_vtx_idx[isrc+1]):
          ivtx = cell_vtx[i] - 1

          buffer_out[n_comp*itgt:n_comp*(itgt+1)] = \
          buffer_out[n_comp*itgt:n_comp*(itgt+1)] + \
          buffer_in[n_comp*ivtx:n_comp*(ivtx+1)] * weight[idx]

          idx += 1

  else:
    print(f"      Error : wrong spatial_interp_algorithm ({spatial_interp_algorithm})")
    sys.exit(1)



def run_coupling():
  # Initialize MPI
  comm = MPI.COMM_WORLD
  i_rank = comm.rank
  n_rank = comm.size

  # Load Python PDM module
  try:
    from Pypdm import Pypdm
  except:
    if i_rank == 0:
      print("      Error : PDM module not found (update PYTHONPATH variable)")
      print(f"pdm : {Pypdm.__file__}")
      sys.exit(1)

  # Load Python CWIPI module
  try:
    from pycwp import pycwp
  except:
    if i_rank == 0:
      print("      Error : CWIPI module not found (update PYTHONPATH variable)")
      print(f"cwp : {pycwp.__file__}")
      sys.exit(1)

  print(f"Python : {i_rank}/{n_rank} je suis là")


  # Initialize CWIPI
  n_code = 1

  code_name = ["codePython"]

  is_active_rank = True

  intra_comm = pycwp.init(comm,
                          code_name,
                          is_active_rank)

  print(f"Python : {i_rank}/{n_rank} CWP_Init OK")

  n_part = 1

  # Generate mesh
  mesh = Pypdm.generate_mesh_rectangle_simplified(intra_comm[0],
                                                  5)
  print(f"Python : {i_rank}/{n_rank} generate_mesh_rectangle_simplified OK")
  # send_field_data = mesh["coords"][0::3] # does not work for some reason...
  send_field_data = np.zeros(mesh["n_vtx"], dtype=np.double)
  for i in range(mesh["n_vtx"]):
    send_field_data[i] = mesh["coords"][3*i]
  recv_field_data = np.zeros(mesh["n_vtx"], dtype=np.double)

  # Create first coupling C <-> Python
  cpl_CP = pycwp.Coupling(code_name[0],
                          "coupling_C_Python",
                          "codeC",
                          pycwp.INTERFACE_SURFACE,
                          pycwp.COMM_PAR_WITH_PART,
                          pycwp.SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,
                          n_part,
                          pycwp.DYNAMIC_MESH_STATIC,
                          pycwp.TIME_EXCH_USER_CONTROLLED)

  print(f"Python : {i_rank}/{n_rank} Coupling C OK")


  # Set coupling visualisation
  cpl_CP.visu_set(1,
                  pycwp.VISU_FORMAT_ENSIGHT,
                  "text")

  # Set the mesh vertices
  cpl_CP.mesh_interf_vtx_set(0,
                             mesh["coords"],
                             None)

  # Set the mesh elements
  block_id = cpl_CP.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  cpl_CP.mesh_interf_f_poly_block_set(0,
                                      block_id,
                                      mesh["elt_vtx_idx"],
                                      mesh["elt_vtx"],
                                      None)

  # Finalize mesh
  cpl_CP.mesh_interf_finalize()

  print(f"Python : {i_rank}/{n_rank} mesh_interf_finalize OK")

  # Define field
  field_name = "coord_x"

  fieldCP = cpl_CP.field_create(field_name,
                                pycwp.DOUBLE,
                                pycwp.FIELD_STORAGE_INTERLACED,
                                1,
                                pycwp.DOF_LOCATION_NODE,
                                pycwp.FIELD_EXCH_SENDRECV,
                                pycwp.STATUS_ON)

   # Create second coupling Python <-> Fortran
  cpl_PF = pycwp.Coupling(code_name[0],
                          "coupling_Python_Fortran",
                          "codeFortran",
                          pycwp.INTERFACE_SURFACE,
                          pycwp.COMM_PAR_WITH_PART,
                          pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                          n_part,
                          pycwp.DYNAMIC_MESH_STATIC,
                          pycwp.TIME_EXCH_USER_CONTROLLED)

  print(f"Python : {i_rank}/{n_rank} Coupling F OK")


  # Set coupling visualisation
  cpl_PF.visu_set(1,
                  pycwp.VISU_FORMAT_ENSIGHT,
                  "text")

  # Set the mesh vertices
  cpl_PF.mesh_interf_vtx_set(0,
                             mesh["coords"],
                             None)

  # Set the mesh elements
  block_id = cpl_PF.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)

  cpl_PF.mesh_interf_f_poly_block_set(0,
                                      block_id,
                                      mesh["elt_vtx_idx"],
                                      mesh["elt_vtx"],
                                      None)

  # Finalize mesh
  cpl_PF.mesh_interf_finalize()


  # Define field
  fieldPF = cpl_PF.field_create(field_name,
                                pycwp.DOUBLE,
                                pycwp.FIELD_STORAGE_INTERLACED,
                                1,
                                pycwp.DOF_LOCATION_NODE,
                                pycwp.FIELD_EXCH_SENDRECV,
                                pycwp.STATUS_ON)

  # Begin time step :
  pycwp.time_step_beg(code_name[0],
                      0.0);


  fieldCP.data_set(0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_field_data)

  fieldCP.data_set(0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_field_data)

  # Set user-defined interpolation function
  fieldCP.interp_function_set(my_interpolation)


  # Spatial interpolation
  cpl_CP.spatial_interp_property_set("n_neighbors",
                                     pycwp.INT,
                                     "3")

  cpl_CP.spatial_interp_weights_compute()

  # Exchange interpolated fields
  fieldCP.issend()
  fieldCP.irecv ()

  fieldCP.wait_issend()
  fieldCP.wait_irecv ()

  # Set field
  fieldPF.data_set(0,
                   pycwp.FIELD_MAP_SOURCE,
                   send_field_data)

  fieldPF.data_set(0,
                   pycwp.FIELD_MAP_TARGET,
                   recv_field_data)

  # Set user-defined interpolation function
  fieldPF.interp_function_set(my_interpolation)


  # Spatial interpolation
  cpl_PF.spatial_interp_property_set("tolerance",
                                     pycwp.DOUBLE,
                                     "0.1")

  cpl_PF.spatial_interp_weights_compute()

  # Exchange interpolated fields
  fieldPF.issend()
  fieldPF.irecv ()

  fieldPF.wait_issend()
  fieldPF.wait_irecv ()

  # End time step :
  pycwp.time_step_end(code_name[0])

  # Delete Mesh
  cpl_CP.mesh_interf_del()
  cpl_PF.mesh_interf_del()

  # Finalize CWIPI
  pycwp.finalize()

  print(f"Python rank {i_rank} FINISHED :D")

  # Finalize MPI
  MPI.Finalize()


if __name__ == '__main__':
  print("C'est parti!")
  run_coupling()

