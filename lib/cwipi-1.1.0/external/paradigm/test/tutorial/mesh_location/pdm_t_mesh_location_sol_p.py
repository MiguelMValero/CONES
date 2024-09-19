#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

import argparse


def run_tuto():
  comm = MPI.COMM_WORLD



  parser = argparse.ArgumentParser()

  parser.add_argument("-n",          "--n_vtx_seg",  type=int, default=10)
  parser.add_argument("-src_n_part", "--src_n_part", type=int, default=1)
  parser.add_argument("-tgt_n_part", "--tgt_n_part", type=int, default=1)
  parser.add_argument("-nodal",      "--nodal",      action="store_true")

  args = parser.parse_args()

  n          = args.n_vtx_seg
  src_n_part = args.src_n_part
  tgt_n_part = args.tgt_n_part
  nodal      = args.nodal

  """
  1) Localization
  """

  # Generate partitioned source mesh
  src_mesh = PDM.generate_mesh_rectangle_ngon(comm          = comm,
                                              elt_type      = PDM._PDM_MESH_NODAL_POLY_2D,
                                              xmin          = 0.,
                                              ymin          = 0.,
                                              zmin          = 0.,
                                              lengthx       = 1.,
                                              lengthy       = 1.,
                                              n_x           = n,
                                              n_y           = n,
                                              n_part        = src_n_part,
                                              part_method   = PDM._PDM_SPLIT_DUAL_WITH_PARMETIS,
                                              random_factor = 0.8)

  src_n_vtx        = src_mesh["pn_vtx"]
  src_n_face       = src_mesh["pn_face"]
  src_vtx_coord    = src_mesh["pvtx_coord"]
  src_face_vtx_idx = src_mesh["pface_edge_idx"]
  src_face_vtx     = src_mesh["pface_vtx"]
  src_face_edge    = src_mesh["pface_edge"]
  src_edge_vtx     = src_mesh["pedge_vtx"]

  src_vtx_ln_to_gn  = src_mesh["pvtx_ln_to_gn"]
  src_face_ln_to_gn = src_mesh["pface_ln_to_gn"]



  # Generate target source mesh
  tgt_mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                              elt_type    = PDM._PDM_MESH_NODAL_QUAD4,
                                              xmin        = 0.3,
                                              ymin        = 0.3,
                                              zmin        = 0,
                                              lengthx     = 1.,
                                              lengthy     = 1.,
                                              n_x         = n,
                                              n_y         = n,
                                              n_part      = tgt_n_part,
                                              part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT)

  tgt_n_vtx        = tgt_mesh["pn_vtx"]
  tgt_n_face       = tgt_mesh["pn_face"]
  tgt_vtx_coord    = tgt_mesh["pvtx_coord"]
  tgt_face_vtx_idx = tgt_mesh["pface_edge_idx"]
  tgt_face_vtx     = tgt_mesh["pface_vtx"]

  tgt_vtx_ln_to_gn  = tgt_mesh["pvtx_ln_to_gn"]
  tgt_face_ln_to_gn = tgt_mesh["pface_ln_to_gn"]



  # Create the MeshLocation object
  mesh_loc = PDM.MeshLocation(1,
                              comm)


  # Set the target point cloud
  mesh_loc.n_part_cloud_set(0,
                            tgt_n_part)

  for i_part in range(tgt_n_part):
    mesh_loc.cloud_set(0,
                       i_part,
                       tgt_vtx_coord   [i_part],
                       tgt_vtx_ln_to_gn[i_part])

  # Set the source mesh
  mesh_loc.mesh_n_part_set(src_n_part)

  for i_part in range(src_n_part):
    if nodal:
      mesh_loc.nodal_part_set_2d(i_part,
                                 src_face_vtx_idx [i_part],
                                 src_face_vtx     [i_part],
                                 src_face_ln_to_gn[i_part],
                                 src_vtx_coord    [i_part],
                                 src_vtx_ln_to_gn [i_part])
    else:
      mesh_loc.part_set_2d(i_part,
                           src_face_vtx_idx [i_part],
                           src_face_edge    [i_part],
                           src_face_ln_to_gn[i_part],
                           src_edge_vtx     [i_part],
                           src_vtx_coord    [i_part],
                           src_vtx_ln_to_gn [i_part])

  # Set the location preconditioning method (optional)
  mesh_loc.method = PDM.MeshLocation.OCTREE

  # Set the geometric tolerance (optional)
  mesh_loc.tolerance = 1e-6

  # Compute location
  mesh_loc.compute()

  # Dump elapsed and CPU times
  mesh_loc.dump_times()


  """
  2) Interpolation
  """

  # Get PartToPart object
  ptp = mesh_loc.part_to_part_get(0)

  # Initiate exchange of first field
  request1 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                       PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                       src_face_ln_to_gn,
                       part1_stride=1,
                       interlaced_str=True)

  # Interpolate second field (node-based)
  src_vtx_field2 = []
  for i_part in range(src_n_part):
    src_vtx_field2.append(np.array(src_vtx_coord[i_part][::3]))

  src_send_field2 = []

  for i_part in range(src_n_part):
    src_result = mesh_loc.points_in_elt_get(0, i_part)
    src_to_tgt_idx = src_result["elt_pts_inside_idx"]
    n_pts = src_to_tgt_idx[src_n_face[i_part]]

    # Interpolate second field (node-based)
    src_connect = mesh_loc.cell_vertex_get(i_part)
    src_cell_vtx_idx = src_connect["cell_vtx_idx"]
    src_cell_vtx     = src_connect["cell_vtx"]

    weights_idx = src_result["points_weights_idx"]
    weights     = src_result["points_weights"]

    field2 = np.zeros(n_pts, dtype=np.double)
    for i_elt in range(src_n_face[i_part]):
      for i_pt in range(src_to_tgt_idx[i_elt], src_to_tgt_idx[i_elt+1]):
        field2[i_pt] = 0

        elt_n_vtx = src_cell_vtx_idx[i_elt+1] - src_cell_vtx_idx[i_elt]
        assert(weights_idx[i_pt+1] - weights_idx[i_pt] == elt_n_vtx)
        for i_vtx in range(elt_n_vtx):
          vtx_id = src_cell_vtx[src_cell_vtx_idx[i_elt] + i_vtx] - 1
          field2[i_pt] += weights[weights_idx[i_pt] + i_vtx] * src_vtx_field2[i_part][vtx_id]

    src_send_field2.append(field2)

  # Initiate exchange of second field
  request2 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                       PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                       src_send_field2,
                       part1_stride=1,
                       interlaced_str=True)

  # Finalize both exchanges
  _, tgt_recv_field1 = ptp.wait(request1)
  _, tgt_recv_field2 = ptp.wait(request2)


  # Check received fields
  pis_located = []
  ptgt_field1 = []
  ptgt_field2 = []
  for i_part in range(tgt_n_part):
    located_tgt   = mesh_loc.located_get  (0, i_part)
    unlocated_tgt = mesh_loc.unlocated_get(0, i_part)

    is_located = np.empty(tgt_n_vtx[i_part], dtype=bool)
    tgt_field1 = np.empty(tgt_n_vtx[i_part], dtype=np.double)
    tgt_field2 = np.empty(tgt_n_vtx[i_part], dtype=np.double)

    for i, i_vtx in enumerate(unlocated_tgt):
      is_located[i_vtx-1] = False
      tgt_field1[i_vtx-1] = -1
      tgt_field2[i_vtx-1] = -1

    for i, i_vtx in enumerate(located_tgt):
      is_located[i_vtx-1] = True
      tgt_field1[i_vtx-1] = tgt_recv_field1[i_part][i]
      tgt_field2[i_vtx-1] = tgt_recv_field2[i_part][i]

      error = abs(tgt_recv_field2[i_part][i] - tgt_vtx_coord[i_part][3*(i_vtx-1)])
      if error > 1e-9:
        print(f"!! error vtx {tgt_vtx_ln_to_gn[i_part][i_vtx]} : {error}")

    pis_located.append(is_located)
    ptgt_field1.append(tgt_field1)
    ptgt_field2.append(tgt_field2)

  # Export for visualization
  PDM.writer_wrapper(comm,
                     "mesh_location_sol_p",
                     "src_mesh",
                     src_vtx_coord,
                     src_vtx_ln_to_gn,
                     src_face_vtx_idx,
                     src_face_vtx,
                     src_face_ln_to_gn,
                     elt_fields={
                     "field1" : src_face_ln_to_gn
                     },
                     vtx_fields={
                     "field2" : src_vtx_field2})

  PDM.writer_wrapper(comm,
                     "mesh_location_sol_p",
                     "tgt_mesh",
                     tgt_vtx_coord,
                     tgt_vtx_ln_to_gn,
                     tgt_face_vtx_idx,
                     tgt_face_vtx,
                     tgt_face_ln_to_gn,
                     vtx_fields={
                     "is_located" : pis_located,
                     "field1"     : ptgt_field1,
                     "field2"     : ptgt_field2})

  if comm.rank == 0:
    print("The End :D")

  MPI.Finalize()

if __name__ == '__main__':
  run_tuto()
