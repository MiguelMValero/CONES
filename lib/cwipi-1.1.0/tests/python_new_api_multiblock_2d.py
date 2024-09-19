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

def generate_mesh(nx=5,
                  ny=None,
                  xmin=0,
                  ymin=0,
                  xmax=1,
                  ymax=1):
  if ny is None:
    ny = nx

  """
  Vertices
  """
  n_vtx1 = 2*nx
  n_vtx2 = nx + 1
  n_vtx3 = n_vtx1 + 2*n_vtx2

  n_vtx = (ny+1)*2*nx + (nx+1)*2*ny + 4

  stepx = (xmax - xmin) / float(3*nx)
  stepy = (ymax - ymin) / float(3*ny)

  n_vtx = 2*nx*(ny+1) + 2*(nx+1)*ny + 4

  coord = np.zeros(n_vtx*3, dtype=np.double)
  ivtx = 0
  for j in range(3*ny+1):
    y = ymin + j*stepy

    if j%3 == 0:
      for ii in range(nx):
        for i in range(2):
          coord[3*ivtx  ] = xmin + (1+i+3*ii)*stepx
          coord[3*ivtx+1] = y
          ivtx += 1
    else:
      for i in range(nx+1):
        coord[3*ivtx  ] = xmin + i*3*stepx
        coord[3*ivtx+1] = y
        ivtx += 1

  # corners
  coord[3*ivtx  ] = xmin
  coord[3*ivtx+1] = ymin
  ivtx += 1

  coord[3*ivtx  ] = xmax
  coord[3*ivtx+1] = ymin
  ivtx += 1

  coord[3*ivtx  ] = xmin
  coord[3*ivtx+1] = ymax
  ivtx += 1

  coord[3*ivtx  ] = xmax
  coord[3*ivtx+1] = ymax



  """
  Faces
  """
  n_tria = 4 + 2*(nx-1 + ny-1)
  n_quad = (nx-1)*(ny-1)
  n_octa = nx*ny

  n_face = n_tria + n_quad + n_octa

  connec_idx = np.zeros(n_face+1, dtype=np.int32)
  for i in range(n_tria):
    connec_idx[i+1] = connec_idx[i] + 3
  for i in range(n_quad):
    connec_idx[n_tria+i+1] = connec_idx[n_tria+i] + 4
  for i in range(n_octa):
    connec_idx[n_tria+n_quad+i+1] = connec_idx[n_tria+n_quad+i] + 8

  connec = np.zeros(3*n_tria + 4*n_quad + 8*n_octa, dtype=np.int32)
  connec[:] = 1

  idx = 0

  first_vtx_last_row = n_vtx - 4 - 2*nx + 1 # 1-based

  # Triangles
  connec[connec_idx[0]  ] = n_vtx - 3
  connec[connec_idx[0]+1] = 1
  connec[connec_idx[0]+2] = 2*nx + 1

  connec[connec_idx[1]  ] = 2*nx
  connec[connec_idx[1]+1] = n_vtx - 2
  connec[connec_idx[1]+2] = 3*nx + 1

  connec[connec_idx[2]  ] = n_vtx - 1
  connec[connec_idx[2]+1] = first_vtx_last_row - nx - 1
  connec[connec_idx[2]+2] = first_vtx_last_row

  connec[connec_idx[3]  ] = n_vtx - 4
  connec[connec_idx[3]+1] = n_vtx - 4 - 2*nx
  connec[connec_idx[3]+2] = n_vtx

  iface = 4
  for i in range(nx-1):
    connec[connec_idx[iface]  ] = 1 + 2*i+1
    connec[connec_idx[iface]+1] = 1 + 2*(i+1)
    connec[connec_idx[iface]+2] = 1 + 2*nx + i + 1
    iface += 1

  for i in range(nx-1):
    connec[connec_idx[iface]  ] = first_vtx_last_row + 2*i+1
    connec[connec_idx[iface]+1] = first_vtx_last_row - (nx+1) + i + 1
    connec[connec_idx[iface]+2] = first_vtx_last_row + 2*(i+1)
    iface += 1

  for j in range(ny-1):
    base = (j+1)*3*nx + j*(nx+2) + 2
    connec[connec_idx[iface]  ] = base
    connec[connec_idx[iface]+1] = base + nx + 1
    connec[connec_idx[iface]+2] = base + 3*nx + 1
    iface += 1

  for j in range(ny-1):
    base = (j+1)*3*nx + j*(nx+2) + 2
    connec[connec_idx[iface]  ] = base + nx
    connec[connec_idx[iface]+1] = base + 4*nx + 1
    connec[connec_idx[iface]+2] = base + 3*nx
    iface += 1

  # Quadrangles
  for j in range(ny-1):
    base = 1 + (j+1)*3*nx + j*(nx+2) + 2

    for i in range(nx-1):
      connec[connec_idx[iface]  ] = base + i
      connec[connec_idx[iface]+1] = base + nx + 2*i + 2
      connec[connec_idx[iface]+2] = base + 3*nx + i + 1
      connec[connec_idx[iface]+3] = connec[connec_idx[iface]+1] - 1
      iface += 1

  # Octagons
  for j in range(ny):
    base = (2*nx + 2*(nx+1))*j + 1
    for i in range(nx):
      connec[connec_idx[iface]  ] = base + 2*i
      connec[connec_idx[iface]+1] = connec[connec_idx[iface]] + 1
      connec[connec_idx[iface]+2] = base + 2*nx + i + 1
      connec[connec_idx[iface]+3] = connec[connec_idx[iface]+2] + nx + 1
      connec[connec_idx[iface]+4] = connec[connec_idx[iface]+1] + 2*nx + 2*(nx+1)
      connec[connec_idx[iface]+5] = connec[connec_idx[iface]+4] - 1
      connec[connec_idx[iface]+6] = connec[connec_idx[iface]+3] - 1
      connec[connec_idx[iface]+7] = connec[connec_idx[iface]+2] - 1
      iface += 1

  return {
    "coord"     : coord,
    "connec_idx": connec_idx,
    "connec"    : connec
  }


def generate_mesh_multiblock(nx=5,
                             ny=None,
                             xmin=0,
                             ymin=0,
                             xmax=1,
                             ymax=1):
  mesh = generate_mesh(nx,
                       ny,
                       xmin,
                       ymin,
                       xmax,
                       ymax)

  blocks = {
    pycwp.BLOCK_FACE_TRIA3: {"connec" : []},
    pycwp.BLOCK_FACE_QUAD4: {"connec" : []},
    pycwp.BLOCK_FACE_POLY : {"connec" : [], "connec_idx": [0]}
  }

  connec_idx = mesh["connec_idx"]
  connec     = mesh["connec"]
  n_elt      = len(connec_idx)-1
  for i in range(n_elt):
    elt_vtx_n = connec_idx[i+1] - connec_idx[i]
    if elt_vtx_n == 3:
      elt_type = pycwp.BLOCK_FACE_TRIA3
    elif elt_vtx_n == 4:
      elt_type = pycwp.BLOCK_FACE_QUAD4
    else:
      elt_type = pycwp.BLOCK_FACE_POLY
      blocks[elt_type]["connec_idx"].append(blocks[elt_type]["connec_idx"][-1] + elt_vtx_n)

    blocks[elt_type]["connec"].extend(connec[connec_idx[i]:connec_idx[i+1]])

  for elt_type in blocks:
    blocks[elt_type]["connec"] = np.array(blocks[elt_type]["connec"], dtype=np.int32)
    if elt_type == pycwp.BLOCK_FACE_POLY:
      blocks[elt_type]["connec_idx"] = np.array(blocks[elt_type]["connec_idx"], dtype=np.int32)

  return {
    "coord"  : mesh["coord"],
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
                       "python_new_api_multiblock_2d",
                       code_names[(i_rank+1)%2],
                       pycwp.INTERFACE_SURFACE,
                       pycwp.COMM_PAR_WITH_PART,
                       pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
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

  # Set blocks (FACE_TRIA3, FACE_QUAD4, FACE_POLY)
  # block_ids = []
  for elt_type in mesh["blocks"]:
    block_id = cpl.mesh_interf_block_add(elt_type)
    # block_ids.append(block_id)

    if elt_type == pycwp.BLOCK_FACE_POLY:
      # Generic polygons
      cpl.mesh_interf_f_poly_block_set(0,
                                       block_id,
                                       mesh["blocks"][elt_type]["connec_idx"],
                                       mesh["blocks"][elt_type]["connec"],
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





