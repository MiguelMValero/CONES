#!/usr/bin/env python

def run_test():
  import numpy as np
  import argparse
  import mpi4py.MPI as MPI
  import Pypdm.Pypdm as PDM

  parser = argparse.ArgumentParser()

  parser.add_argument("-fe", "--nodal", action="store_true")

  args = parser.parse_args()

  fe = args.nodal

  # Initialize MPI environment
  comm   = MPI.COMM_WORLD
  n_rank = comm.size
  i_rank = comm.rank

  # Generate block-distributed parallelepided mesh
  n_x      = 10
  n_y      = 10
  n_z      = 10
  lengthx  = 1.
  xmin     = 0.
  ymin     = 0.
  zmin     = 0.
  elt_type = PDM._PDM_MESH_NODAL_TETRA4
  order    = 1 # call PDM_dcube_nodal_gen_ordering_set if order > 1
  dcube = PDM.DCubeNodalGenerator(n_x,
                                  n_y,
                                  n_z,
                                  lengthx,
                                  xmin,
                                  ymin,
                                  zmin,
                                  elt_type,
                                  order,
                                  comm)

  dcube.compute()

  dmn = dcube.get_dmesh_nodal()

  PDM.generate_distribution(dmn)

  # Create partitioning object
  n_domain = 1 # fixed
  n_part   = 1 # fixed
  i_part   = 0 # fixed
  i_domain = 0 # fixed
  part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT;
  mpart = PDM.MultiPart(n_domain,
                        np.array([n_part]).astype(np.intc),
                        0,
                        part_method,
                        1,
                        np.ones(1).astype(np.double),
                        comm)

  renum_cell = bytes("PDM_PART_RENUM_CELL_NONE", 'ascii')
  renum_face = bytes("PDM_PART_RENUM_FACE_NONE", 'ascii')
  mpart.reordering_set(-1, # i_domain
                       renum_cell,
                       None,
                       renum_face)

  mpart.dmesh_nodal_set(i_domain, dmn)

  mpart.compute()

  # Get mesh arrays in FE structure
  if fe:
    coords = mpart.vtx_coord_get(i_domain,
                                 i_part)

    pmn = mpart.part_mesh_nodal_get(i_domain)
    i_section = 0 # fixed
    output = PDM.part_mesh_nodal_get_sections(pmn,
                                              PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                                              i_part)
    elt_vtx      = output[i_section]["np_connec"]
    elt_ln_to_gn = output[i_section]["np_numabs"]
    n_elt        = len(elt_ln_to_gn)
    vtx_ln_to_gn = PDM.part_mesh_nodal_vtx_g_num_get(pmn, i_part)
    n_vtx        = len(vtx_ln_to_gn)

  # Get mesh arrays in FV structure
  else :
    vtx_ln_to_gn = mpart.ln_to_gn_get(i_domain,
                                      i_part,
                                      PDM._PDM_MESH_ENTITY_VTX)
    n_vtx = len(vtx_ln_to_gn)

    coords = mpart.vtx_coord_get(i_domain,
                                 i_part)

    # edge_ln_to_gn = mpart.ln_to_gn_get(i_domain,
    #                                    i_part,
    #                                    PDM._PDM_MESH_ENTITY_EDGE)
    # n_edge = len(edge_ln_to_gn)

    # _, edge_vtx = mpart.connectivity_get(i_domain,
    #                                      i_part,
    #                                      PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    face_ln_to_gn = mpart.ln_to_gn_get(i_domain,
                                       i_part,
                                       PDM._PDM_MESH_ENTITY_FACE)
    n_face = len(face_ln_to_gn)

    # face_edge_idx, face_edge = mpart.connectivity_get(i_domain,
    #                                                   i_part,
    #                                                   PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
    face_vtx_idx, face_vtx = mpart.connectivity_get(i_domain,
                                                    i_part,
                                                    PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    cell_ln_to_gn = mpart.ln_to_gn_get(i_domain,
                                       i_part,
                                       PDM._PDM_MESH_ENTITY_CELL)

    n_cell = len(cell_ln_to_gn)

    cell_face_idx, cell_face = mpart.connectivity_get(i_domain,
                                                      i_part,
                                                      PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    # BONUS
    if fe:
      print("The end (FE)", flush=True)
      return

    # step 1 : create
    extend_type = PDM._PDM_EXTEND_FROM_VTX
    depth       = 1
    part_ext = PDM.PartExtension(n_domain,
                                 np.array([n_part]).astype(np.intc),
                                 extend_type,
                                 depth,
                                 comm)

    # step 2 : set
    output = mpart.graph_comm_get(i_domain,
                                  i_part,
                                  PDM._PDM_MESH_ENTITY_VTX)

    vtx_part_bound_proc_idx = output["np_entity_part_bound_proc_idx"]
    vtx_part_bound_part_idx = output["np_entity_part_bound_part_idx"]
    vtx_part_bound          = output["np_entity_part_bound"]

    # part_ext.set_part(i_domain,
    #                   i_part,
    #                   n_cell,
    #                   n_face,
    #                   0, # n_face_part_bound
    #                   0, # n_face_group
    #                   n_edge,
    #                   n_vtx,
    #                   cell_face_idx,
    #                   cell_face,
    #                   None, # face_cell
    #                   face_edge_idx,
    #                   face_edge,
    #                   None, # face_vtx_idx
    #                   None, # face_vtx
    #                   edge_vtx,
    #                   None, # face_group_idx
    #                   None, # face_group
    #                   None, # face_join_idx
    #                   None, # face_join
    #                   None, # face_part_bound_proc_idx
    #                   None, # face_part_bound_part_idx
    #                   None, # face_part_bound
    #                   vtx_part_bound_proc_idx,
    #                   vtx_part_bound_part_idx,
    #                   vtx_part_bound,
    #                   cell_ln_to_gn,
    #                   face_ln_to_gn,
    #                   edge_ln_to_gn,
    #                   vtx_ln_to_gn,
    #                   None, # face_group_ln_to_gn
    #                   coords)

    part_ext.connectivity_set(i_domain,
                              i_part,
                              PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE,
                              cell_face_idx,
                              cell_face)

    part_ext.connectivity_set(i_domain,
                              i_part,
                              PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX,
                              face_vtx_idx,
                              face_vtx)

    part_ext.vtx_coord_set(i_domain,
                           i_part,
                           coords)

    part_ext.ln_to_gn_set(i_domain,
                          i_part,
                          PDM._PDM_MESH_ENTITY_CELL,
                          cell_ln_to_gn)

    part_ext.ln_to_gn_set(i_domain,
                          i_part,
                          PDM._PDM_MESH_ENTITY_FACE,
                          face_ln_to_gn)

    part_ext.ln_to_gn_set(i_domain,
                          i_part,
                          PDM._PDM_MESH_ENTITY_VTX,
                          vtx_ln_to_gn)

    part_ext.part_bound_graph_set(i_domain,
                                  i_part,
                                  PDM._PDM_MESH_ENTITY_VTX,
                                  vtx_part_bound_proc_idx,
                                  vtx_part_bound_part_idx,
                                  vtx_part_bound)

    # step 3 : compute
    part_ext.compute()

    # step 4 : get
    # Cell
    cell_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                              i_part,
                                              PDM._PDM_MESH_ENTITY_CELL)

    cell_face_ext_idx, cell_face_ext = part_ext.connectivity_get(i_domain,
                                                                 i_part,
                                                                 PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    # Face
    face_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                              i_part,
                                              PDM._PDM_MESH_ENTITY_FACE)

    # face_edge_ext_idx, face_edge_ext = part_ext.connectivity_get(i_domain,
    #                                                              i_part,
    #                                                              PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
    face_vtx_ext_idx, face_vtx_ext = part_ext.connectivity_get(i_domain,
                                                               i_part,
                                                               PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)

    # Edge
    # edge_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
    #                                           i_part,
    #                                           PDM._PDM_MESH_ENTITY_EDGE)

    # edge_vtx_ext_idx, edge_vtx_ext = part_ext.connectivity_get(i_domain,
    #                                                            i_part,
    #                                                            PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    # Vertices
    vtx_ext_ln_to_gn = part_ext.ln_to_gn_get(i_domain,
                                             i_part,
                                             PDM._PDM_MESH_ENTITY_VTX)

    vtx_coord_ext = part_ext.vtx_coord_get(i_domain,
                                           i_part)


    # dummy gets to prevent leaks
    interface_cell = part_ext.get_interface(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_CELL)

    interface_face = part_ext.get_interface(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_FACE)

    interface_vtx  = part_ext.get_interface(i_domain,
                                            i_part,
                                            PDM._PDM_MESH_ENTITY_VTX)

    surfaces = part_ext.group_get(i_domain,
                                  i_part,
                                  PDM._PDM_MESH_ENTITY_FACE)

    # part_ext.get_composed_interface()

    # step 5 : free (implicit in Python)

if __name__ == '__main__':
  run_test()
