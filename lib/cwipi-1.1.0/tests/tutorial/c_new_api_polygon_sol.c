/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "cwp.h"
#include "cwp_priv.h"

/*----------------------------------------------------------------------
 *
 * Main : base test : Polygon
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // This test mimics the coupling between 2 codes running each
  // on one processor.
  assert (n_rank == 2);

  // Initialize CWIPI :
  // In this version of CWIPI each MPI rank can execute multiple codes.
  // In this particular tutorial, 2 codes are coupled :
  // 'code1' runs on MPI rank 0 and 'code2' runs on MPI rank 1.
  //

  // Here 2 codes are coupled. code1 runs on the processor of
  // MPI rank 0 and code2 runs on the processor if MPI rank 1.
  // In this version of CWIPI several codes can execute on the
  // same MPI rank (here only one code per processor, so n_code = 1).
  // Therefore an array of code names is given at initialization.
  // is_active_rank tells if current ranks will be used in
  // the CWIPI coupling computations.
  // intra_comm is an array of MPI communicators
  // giving the for each code on the processors the communicator
  // to communicate through the ranks of that code.
  int n_code = 1;
  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t  is_active_rank = CWP_STATUS_ON;
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  int I_am_code1 = 0;
  int I_am_code2 = 0;

  if (i_rank == 0) {
    code_name[0] = "code1";
    I_am_code1   = 1;
  }

  if (i_rank == 1) {
    code_name[0] = "code2";
    I_am_code2   = 1;
  }
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  // Create the coupling :
  // One CWIPI context can hold several couplings. Let us set up the
  // coupling between code1 and code2. CWP_INTERFACE_SURFACE informs
  // that the geometrical interface of the meshes of the coupled
  // codes is a surface, still for CWIPI the coordinate system is 3D.
  // CWP_COMM_PAR_WITH_PART means that each mesh is partitionned
  // over the processors of its code. Here the mesh does not change
  // over the coupling, so CWP_DYNAMIC_MESH_STATIC is set.
  // CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE designates
  // which spatial interpolation algorithm is used to map data between
  // the two codes.
  // CWP_TIME_EXCH_USER_CONTROLLED is not used yet.
  int n_part = 1;
  const char  *coupling_name     = "code1_code2";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);

  if (I_am_code1) {
    coupled_code_name[0] = "code2";
  }

  if (I_am_code2) {
    coupled_code_name[0] = "code1";
  }
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // Set coupling visualisation:
  // Output files of the code fields will be written in Ensight ASCII
  // format (easily readable by Paraview) at each iteration.
  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  // Create the field
  // It is possible to operate a bidirectional exchange (see c_new_vs_old_sendrecv).
  // For sake of simplicity, this example will only send the field
  // of code1 (CWP_FIELD_EXCH_SEND) to code2 (CWP_FIELD_EXCH_RECV).
  // On code1 there is a field (CWP_FIELD_MAP_SOURCE) located at
  // the vertices (CWP_DOF_LOCATION_NODE) with one component (n_components)
  // which is the x coordinate of the mesh in this test.
  const char *field_name      = "a super fancy field";
  int         n_components    = 1;

  if (I_am_code1) {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     n_components,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     CWP_STATUS_ON);
  }

  if (I_am_code2) {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     n_components,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     CWP_STATUS_ON);
  }

  // Begin time step :
  // In this example there is only one time step. It is mandatory to create the
  // coupling and the associated fields before starting the first time step.
  CWP_Time_step_beg(code_name[0],
                    0.0);

  // Set the mesh vertices coordinates :
  // The coordinate system in CWIPI is always 3D, so
  // we allocate an array of the time the number of vertices
  // (11 here) to set the coordinates in. The coordinates are
  // interlaced (x0, y0, z0, x1, y1, z1, ..., xn, yn, zn).
  // The NULL argument will be explained later.
  int    n_vtx = 11;
  double coords[33] = {0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0,
                       3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0};
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

  // Set the mesh polygons connectivity :
  // Let us set a mesh of 5 polygons (CWP_BLOCK_FACE_POLY).
  // An index array (connec_idx) of size n_elts+1 contains the
  // information of the number of vertices per polygon. The first
  // index is always 0, from there we add up the number of vertices
  // per element. Here one triangle, 2 quadrangles and 2 pentagons.
  // The connectivity between elements and vertices is an array of
  // size connec_idx(n_elts+1) (here 21).
  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling_name,
                                           CWP_BLOCK_FACE_POLY);

  int n_elts = 5;
  int connec_idx[6] = {0,3,7,11,16,21};
  int connec[21]    = {1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8};
  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   n_elts,
                                   connec_idx,
                                   connec,
                                   NULL);

  // Finalize mesh :
  // CWIPI hides the parallelism for users, so it is not
  // mandatory to give a global numbering for mesh data (the
  // NULL arguments earlier). If not given this numbering is
  // generated by CWIPI by the following function,
  // as well as the underlying mesh data structure
  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  // Set the field values :
  // Note that the user has to allocate the array for the
  // field that will be received by code2 (CWP_FIELD_MAP_TARGET).
  double     *send_field_data = malloc(sizeof(double) * n_vtx * n_components);
  double     *recv_field_data = malloc(sizeof(double) * n_vtx * n_components);

  // for code1
  if (I_am_code1) {
    for (int i = 0; i < n_vtx; i++) {
      send_field_data[i] = coords[3*i];
    }
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_field_data);
  }

  // for code2
  if (I_am_code2) {
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_field_data);
  }

  // Compute interpolation weights :
  // Set a geometric tolerance of 10% of an element size for
  // point localisation.
  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  CWP_DOUBLE,
                                  "0.1");
  CWP_Spatial_interp_weights_compute(code_name[0],
                                     coupling_name);

  // Exchange field values between codes :
  // The field exchange functions mimic the way the associated
  // MPI functions work, see MPI documentation for more information.
  // for code1
  if (I_am_code1) {
    CWP_Field_issend(code_name[0],
                     coupling_name,
                     field_name);
  }

  // for code2
  if (I_am_code2) {
    CWP_Field_irecv(code_name[0],
                    coupling_name,
                    field_name);
  }

  // for code1
  if (I_am_code1) {
    CWP_Field_wait_issend(code_name[0],
                          coupling_name,
                          field_name);
  }

  // for code2
  if (I_am_code2) {
    CWP_Field_wait_irecv(code_name[0],
                         coupling_name,
                         field_name);
  }

  // Check interpolation :
  // These functions allow to know how many and for which target
  // vertices the interpolation operation has been unsuccessful.
  int        n_uncomputed_tgts = -1;
  const int *uncomputed_tgts   = NULL;
  if (I_am_code2) {
    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0],
                                                  coupling_name,
                                                  field_name,
                                                  0);

    uncomputed_tgts = CWP_Uncomputed_tgts_get(code_name[0],
                                              coupling_name,
                                              field_name,
                                              0);
  }

  // End time step :
  CWP_Time_step_end(code_name[0]);

  // Delete field :
  CWP_Field_del(code_name[0],
                coupling_name,
                field_name);

  // Delete Mesh :
  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling_name);

  // free
  CWP_UNUSED(n_uncomputed_tgts);
  CWP_UNUSED(uncomputed_tgts);
  if (send_field_data != NULL) free(send_field_data);
  if (recv_field_data != NULL) free(recv_field_data);
  free(code_name);
  free(intra_comm);
  free(coupled_code_name);

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}

