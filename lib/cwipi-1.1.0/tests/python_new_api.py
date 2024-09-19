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

f=None

def my_interpolation(field,
                     i_part,
                     buffer_in,
                     buffer_out):

    src_data = field.src_data_properties_get(i_part)
    weights  = field.location_weights_get(i_part)
    connec   = field.location_internal_cell_vtx_get(i_part)

    buffer_out[:] = 0

    idx = 0
    for i in range(src_data["n_src"]):
        for j in range(src_data["src_to_tgt_idx"][i],src_data["src_to_tgt_idx"][i+1]):
            for k in range(connec["cell_vtx_idx"][i],connec["cell_vtx_idx"][i+1]):
                ivtx = connec["cell_vtx"][k]-1
                buffer_out[field.n_components*j:field.n_components*(j+1)] = \
                buffer_out[field.n_components*j:field.n_components*(j+1)] + \
                buffer_in[field.n_components*ivtx:field.n_components*(ivtx+1)] * weights[idx]
                idx += 1

def userInterp(local_code_name,
               cpl_id,
               field_id,
               i_part,
               buffer_in,
               buffer_out):
    """
    User interpolation method
    """
    global f

    f.write("in a user interpolation function !")

def runTest():
    """
    Run tests on Python interface of new API
    """
    global f
    comm = MPI.COMM_WORLD

    i_rank = comm.rank
    n_rank = comm.size

    if (i_rank == 0):
        print("\nSTART: python_new_api.py")

    if (n_rank != 2):
        if i_rank == 0:
            print("      Not executed : only available for 2 processes")
        return

    code_names = ["proc0","proc1"]

    # if (i_rank == 0):
    #     code_name = ["proc0"]

    # if (i_rank == 1):
    #     code_name = ["proc1"]
    code_name = code_names[i_rank]

    try:
        from pycwp import pycwp
    except:
        if i_rank == 0:
            print("      Error : CWIPI module not found (update PYTHONPATH variable)")
        sys.exit(1)

    # OUTPUT
    srank = '{0}'.format(i_rank)
    f=open("python_new_api_"+srank.zfill(4)+".txt",'w')
    pycwp.output_file_set(f)

    # INIT
    f.write("pycwp.init:\n")
    n_code = 1
    is_active_rank = True
    out = pycwp.init(comm,
                     [code_name],
                     is_active_rank)
    f.write("  - intra_comms : {param}\n".format(param=out[0]))

    # STATE UPDATE
    pycwp.state_update(code_name, pycwp.STATE_IN_PROGRESS)
    f.write("pycwp.state_get:\n")
    state = pycwp.state_get(code_name)
    f.write("  - state : {param}\n".format(param=state))
    pycwp.state_update(code_name, pycwp.STATE_END)
    f.flush()

    # PROPERTIES DUMP
    # f.write("pycwp.properties_dump:\n")
    # pycwp.properties_dump()

    # CODES
    f.write("pycwp.code:\n")
    n_code = pycwp.codes_nb_get()
    code   = pycwp.codes_list_get()
    n_loc_code = pycwp.loc_codes_nb_get()
    loc_code   = pycwp.loc_codes_list_get()
    f.write("  - n_code : {param}\n".format(param=n_code))
    for i in range(n_code):
        f.write("    --> {param}\n".format(param=code[i]))
    f.write("  - n_loc_code : {param}\n".format(param=n_loc_code))
    for i in range(n_loc_code):
        f.write("    --> {param}\n".format(param=loc_code[i]))
    f.flush()

    # PARAM
    f.write("param_add_dbl:\n")
    pycwp.param_lock(code_name)
    pycwp.param_add_dbl(code_name, "double", 0.5)
    pycwp.param_unlock(code_name)
    f.flush()

    pycwp.param_lock(code_name)
    pycwp.param_add_str(code_name, "str", "chat")
    pycwp.param_unlock(code_name)

    if (i_rank == 0):
        pycwp.param_lock(code_name)
        pycwp.param_add_int(code_name, "entier", 1)
        pycwp.param_unlock(code_name)

    if (i_rank == 1):
        pycwp.param_lock(code_name)
        pycwp.param_add_int(code_name, "entier", -1)
        pycwp.param_unlock(code_name)

    comm.Barrier()

    f.write("cwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_name, "double", pycwp.DOUBLE)
    f.write("  - value (0): {param}\n".format(param=value))

    pycwp.param_lock(code_name)
    pycwp.param_set_dbl(code_name, "double", 0.25)
    pycwp.param_unlock(code_name)

    pycwp.param_lock(code_name)
    pycwp.param_set_str(code_name, "str", "chien")
    pycwp.param_unlock(code_name)

    if (i_rank == 0):
        pycwp.param_lock(code_name)
        pycwp.param_set_int(code_name, "entier", 2)
        pycwp.param_unlock(code_name)

    comm.Barrier()

    f.write("pycwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_name, "double", pycwp.DOUBLE)
    f.write("  - value (1): {param}\n".format(param=value))

    pycwp.param_lock(code_name)
    pycwp.param_del(code_name, "str", pycwp.CHAR)
    pycwp.param_unlock(code_name)

    comm.Barrier()

    f.write("pycwp.param_n_get:\n")
    n_param_str = pycwp.param_n_get(code_name, pycwp.CHAR)
    n_param_int = pycwp.param_n_get(code_name, pycwp.INT)
    f.write("  - n_param_str: {param}\n".format(param=n_param_str))
    f.write("  - n_param_int: {param}\n".format(param=n_param_int))

    f.write("pycwp.param_list_get:\n")
    str_param = pycwp.param_list_get(code_name, pycwp.CHAR)
    for name in str_param:
        f.write(f"    --> str_param: {name}\n")

    f.write("pycwp.param_is:\n")
    exists = pycwp.param_is(code_name, "entier", pycwp.INT)
    f.write(f"  - exists 'entier': {exists}\n")
    exists = pycwp.param_is(code_name, "chapeau", pycwp.INT)
    f.write(f"  - exists 'chapeau': {exists}\n")

    comm.Barrier()

    f.write("pycwp.param_list_get:\n")
    int_param = pycwp.param_list_get(code_name, pycwp.INT)
    for name in int_param:
        f.write(f"    --> int_param: {name}\n")

    f.write("pycwp.param_get ({param}):\n".format(param=i_rank))
    value = pycwp.param_get(code_name, "entier", pycwp.INT)
    f.write("  - value int: {param}\n".format(param=value))

    comm.Barrier()

    f.write("pycwp.param_reduce:\n")
    result = pycwp.param_reduce(pycwp.OP_MIN, "entier",  pycwp.INT, code_names)
    f.write("  - result: {param}\n".format(param=result))
    f.flush()

    # Cpl
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl = pycwp.Coupling(code_name,
                         "test",
                         code_names[(i_rank+1)%2],
                         pycwp.INTERFACE_SURFACE, #INTERFACE_LINEAR,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         1,
                         pycwp.DYNAMIC_MESH_VARIABLE,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    cpl.visu_set(1,
                 pycwp.VISU_FORMAT_ENSIGHT,
                 "text")

    # ------------------------------------------------------------------
    # std or polygon
    polygon = True

    # MESH
    if (polygon):
        if (i_rank == 0):
            coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
            connec_idx = np.array([0, 3, 6], dtype=np.int32)
            connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

        if (i_rank == 1):
            coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
            connec_idx = np.array([0, 3, 6], dtype=np.int32)
            connec = np.array([1, 2, 4, 1, 4, 3], dtype=np.int32)

        f.write("cpl.mesh_interf_vtx_set:\n")
        f.flush()
        cpl.mesh_interf_vtx_set(0,
                                coord,
                                None)

        f.write("cpl.mesh_interf_block_add:\n")
        f.flush()
        block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
        f.write("cpl.mesh_interf_f_poly_block_set ({param}):\n".format(param=i_rank))
        f.flush()

        cpl.mesh_interf_f_poly_block_set(0,
                                         block_id,
                                         connec_idx,
                                         connec,
                                         None)

        f.write("cpl.mesh_interf_finalize:\n")
        f.flush()
        cpl.mesh_interf_finalize()

        f.write("cpl.mesh_interf_f_poly_block_get:\n")
        f.flush()
        out = cpl.mesh_interf_f_poly_block_get(0, block_id)
        f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
        f.write("  - connec_idx {param}\n".format(param=out["connec_idx"]))
        f.write("  - connec {param}\n".format(param=out["connec"]))
        f.write("  - global_num : {param}\n".format(param=out["global_num"]))
        f.flush()

        # FIELD

        f.write("before creating field\n")
        f.flush()
        if (i_rank == 0):
            field = cpl.field_create("champs",
                                     pycwp.DOUBLE,
                                     pycwp.FIELD_STORAGE_INTERLACED,
                                     1,
                                     pycwp.DOF_LOCATION_NODE,
                                     pycwp.FIELD_EXCH_SEND,
                                     pycwp.STATUS_ON)

        if (i_rank == 1):
            field = cpl.field_create("champs",
                                     pycwp.DOUBLE,
                                     pycwp.FIELD_STORAGE_INTERLACED,
                                     1,
                                     pycwp.DOF_LOCATION_NODE,
                                     pycwp.FIELD_EXCH_RECV,
                                     pycwp.STATUS_ON)

        pycwp.time_step_beg(code_name, 0.)

        # sendField = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.double)
        sendField = np.zeros(4, dtype=np.double)
        recvField = np.zeros(4, dtype=np.double)
        for i in range(4):
            sendField[i] = coord[3*i]

        if (i_rank == 0):
            field.data_set(0,
                           pycwp.FIELD_MAP_SOURCE,
                           sendField)

        if (i_rank == 1):
            field.data_set(0,
                           pycwp.FIELD_MAP_TARGET,
                           recvField)
        f.write("after creating field\n")
        f.flush()

        # USER INTERPOLATION
        f.write("field.interp_function_set:\n")
        f.flush()
        field.interp_function_set(my_interpolation)

        comm.Barrier()

        f.write("cpl.spatial_interp_property_set:\n")
        f.flush()
        cpl.spatial_interp_property_set("tolerance", pycwp.DOUBLE, "1e-2")

        comm.Barrier()

        f.write("cpl.spatial_interp_weights_compute:\n")
        f.flush()
        cpl.spatial_interp_weights_compute()

        if (i_rank == 0):
            f.write("field.issend (0):\n")
            f.flush()
            field.issend()

        if (i_rank == 1):
            f.write("field.irecv (1):\n")
            f.flush()
            field.irecv()

        if (i_rank == 0):
            f.write("field.wait_issend (0):\n")
            f.flush()
            field.wait_issend()

        if (i_rank == 1):
            f.write("field.wait_irecv (1):\n")
            f.flush()
            field.wait_irecv()

        comm.Barrier()

        pycwp.time_step_end(code_name)

        f.write("field del:\n")
        f.flush()

        del field

        # USER STRUCTURE
        class userClass:
            animal = "chat"
            aliment = "aligot"
            ville = "Toulouse"

        userObj = userClass()

        f.write("pycwp.user_structure_set:\n")
        f.flush()
        pycwp.user_structure_set(code_name,
                                 userObj)
        f.write("pycwp.user_structure_get:\n")
        f.flush()
        user_structure = pycwp.user_structure_get(code_name)

        print(user_structure.animal)
        print(user_structure.aliment)
        print(user_structure.ville)

        f.write("cpl.mesh_interf_del:\n")
        f.flush()
        cpl.mesh_interf_del()

    else:

        # STD MESH
        f.write("cpl.mesh_interf_block_add:\n")
        f.flush()
        block_id = cpl.mesh_interf_block_add(pycwp.BLOCK_FACE_TRIA3)
        f.write("cpl.mesh_interf_block_std_set:\n")
        f.flush()
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)
        cpl.mesh_interf_block_std_set(0,
                                      block_id,
                                      connec,
                                      None)

        f.write("cpl.mesh_interf_block_std_get:\n")
        f.flush()
        out = cpl.mesh_interf_block_std_get(0,
                                            block_id)
        f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
        f.write("  - connec {param}\n".format(param=out["connec"]))
        f.write("  - global_num : {param}\n".format(param=out["global_num"]))
        f.flush()

        f.write("cpl.mesh_interf_from_faceedge_set:\n")
        f.flush()
        face_edge_idx = np.array([0, 3, 6], dtype=np.int32)
        face_edge = np.array([1, 2, 3, 5, 4, 3], dtype=np.int32)
        edge_vtx_idx = np.array([0, 2, 4, 6, 8, 10], dtype=np.int32)
        edge_vtx = np.array([3, 1, 1, 2, 2, 3, 4, 2, 3, 4], dtype=np.int32)
        cpl.mesh_interf_from_faceedge_set(0,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx_idx,
                                          edge_vtx,
                                          None)

        f.write("cpl.mesh_interf_del:\n")
        f.flush()
        cpl.mesh_interf_del()

    # ------------------------------------------------------------------
    # Volumic Cpl
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl2 = pycwp.Coupling(code_name,
                         "test_vol",
                         code_names[(i_rank+1)%2],
                         pycwp.INTERFACE_VOLUME,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                         1,
                         pycwp.DYNAMIC_MESH_STATIC,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    block_id = cpl2.mesh_interf_block_add(pycwp.BLOCK_CELL_POLY)

    connec_faces_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    connec_faces = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    connec_cells_idx = np.array([0, 4], dtype=np.int32)
    connec_cells = np.array([1, 2, 3, 4], dtype=np.int32)

    f.write("cpl2.mesh_interf_c_poly_block_set:\n")
    f.flush()
    cpl2.mesh_interf_c_poly_block_set(0,
                                     block_id,
                                     connec_faces_idx,
                                     connec_faces,
                                     connec_cells_idx,
                                     connec_cells,
                                     None)

    f.write("cpl2.mesh_interf_c_poly_block_get:\n")
    f.flush()
    out = cpl2.mesh_interf_c_poly_block_get(0,
                                           block_id)
    f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
    f.write("  - n_faces : {param}\n".format(param=out["n_faces"]))
    f.write("  - connec_faces_idx {param}\n".format(param=out["connec_faces_idx"]))
    f.write("  - connec_faces {param}\n".format(param=out["connec_faces"]))
    f.write("  - connec_cells_idx {param}\n".format(param=out["connec_cells_idx"]))
    f.write("  - connec_cells {param}\n".format(param=out["connec_cells"]))

    face_vtx_idx = np.array([0, 3, 6, 9, 12], dtype=np.int32)
    face_vtx = np.array([1, 2, 3, 1, 2, 4, 2, 3, 4, 1, 3, 4], dtype=np.int32)
    cell_face_idx = np.array([0, 4], dtype=np.int32)
    cell_face = np.array([1, 2, 3, 4], dtype=np.int32)

    f.write("cpl2.mesh_interf_from_cellface_set:\n")
    f.flush()
    cpl2.mesh_interf_from_cellface_set(0,
                                       cell_face_idx,
                                       cell_face,
                                       face_vtx_idx,
                                       face_vtx,
                                       None)

    f.write("cpl2.mesh_interf_del:\n")
    f.flush()
    cpl2.mesh_interf_del()

    # ------------------------------------------------------------------
    # High-order
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl3 = pycwp.Coupling(code_name,
                          "test_ho",
                          code_names[(i_rank+1)%2],
                          pycwp.INTERFACE_SURFACE,
                          pycwp.COMM_PAR_WITH_PART,
                          pycwp.SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                          1,
                          pycwp.DYNAMIC_MESH_STATIC,
                          pycwp.TIME_EXCH_USER_CONTROLLED)

    block_id = cpl3.mesh_interf_block_add(pycwp.BLOCK_FACE_TRIAHO)

    vtx_coord = np.zeros(6*3, dtype=np.double)
    ijk_grid  = np.array([0, 0,
                         2, 0,
                         0, 2,
                         1, 0,
                         1, 1,
                         0, 1], dtype=np.int32)
    k = 0
    for j in range(3):
        for i in range(3):
            vtx_coord[k] = i
            k += 1
            vtx_coord[k] = j
            k += 1
    face_vtx = np.array([1, 3, 6, 2, 5, 4], dtype=np.int32)

    f.write("cpl3.mesh_interf_block_ho_set:\n")
    f.flush()
    cpl3.mesh_interf_block_ho_set(0,
                                  block_id,
                                  2,
                                  face_vtx,
                                  None)

    cpl3.mesh_interf_ho_ordering_from_IJK_set(pycwp.BLOCK_FACE_TRIAHO,
                                              2,
                                              ijk_grid)

    f.write("mesh_interf_block_ho_get:\n")
    f.flush()
    out = cpl3.mesh_interf_block_ho_get(0,
                                        block_id)
    f.write("  - n_elts : {param}\n".format(param=out["n_elts"]))
    f.write("  - order {param}\n".format(param=out["order"]))
    f.write("  - connec {param}\n".format(param=out["connec"]))

    f.write("cpl3.mesh_interf_del:\n")
    f.flush()
    cpl3.mesh_interf_del()

    # ------------------------------------------------------------------
    # user target coupling
    f.write("pycwp.Coupling:\n")
    f.flush()
    cpl4 = pycwp.Coupling(code_name,
                         "test_user_target",
                         code_names[(i_rank+1)%2],
                         pycwp.INTERFACE_SURFACE,
                         pycwp.COMM_PAR_WITH_PART,
                         pycwp.SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,
                         1,
                         pycwp.DYNAMIC_MESH_STATIC,
                         pycwp.TIME_EXCH_USER_CONTROLLED)

    # face poly mesh
    if (i_rank == 0):
        coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 3, 2, 4, 3], dtype=np.int32)

    if (i_rank == 1):
        coord = np.array([0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0], dtype=np.double)
        connec_idx = np.array([0, 3, 6], dtype=np.int32)
        connec = np.array([1, 2, 4, 1, 4, 3], dtype=np.int32)

    f.write("cpl4.mesh_interf_vtx_set:\n")
    f.flush()
    cpl4.mesh_interf_vtx_set(0,
                            coord,
                            None)

    f.write("cpl4.mesh_interf_block_add:\n")
    f.flush()
    block_id = cpl4.mesh_interf_block_add(pycwp.BLOCK_FACE_POLY)
    f.write("cpl4.mesh_interf_f_poly_block_set ({param}):\n".format(param=i_rank))
    f.flush()

    cpl4.mesh_interf_f_poly_block_set(0,
                                     block_id,
                                     connec_idx,
                                     connec,
                                     None)

    # user targets
    user_tgt_coords = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4], dtype=np.double)
    cpl4.user_tgt_pts_set(0,
                          4,
                          user_tgt_coords,
                          None)

    # finalize
    f.write("cpl4.mesh_interf_finalize:\n")
    f.flush()
    cpl4.mesh_interf_finalize()

    # delete
    f.write("cpl4.mesh_interf_del:\n")
    f.flush()
    cpl4.mesh_interf_del()

    # FINALIZE
    pycwp.finalize()

    # END
    f.write("\nEnd.\n")
    f.close()
    comm.Barrier()
    MPI.Finalize()

if __name__ == '__main__':
    runTest()
