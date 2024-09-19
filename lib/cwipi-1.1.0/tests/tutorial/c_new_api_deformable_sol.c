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
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_generate_mesh.h"

/*----------------------------------------------------------------------
 *
 * Main : advanced test : Deformable (code 2)
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI :
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Initialize CWIPI :
  int n_code = 1;

  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t  is_active_rank = CWP_STATUS_ON;
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  code_name[0]      = "code2";

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  // Create the coupling :
  // CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  // to the mesh over the coupling steps.
  int n_part = 1;
  const char  *coupling_name     = "coupling";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "code1";
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_DEFORMABLE,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // Set coupling visualisation:
  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  // Create mesh :
  // It is the users responsability to free arrays from _simplified mesh generation functions.
  // PDM_MPI_mpi_2_pdm_mpi_comm is used when calling ParaDiGM library functions to
  // convert the MPI communicator into the expected format.
  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;
  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[0]),
                                         10,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx);

  int n_components = 1;

  // Create the field to be send:
  const char *send_field_name = "chinchilla";
  CWP_Field_create(code_name[0],
                   coupling_name,
                   send_field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_SEND,
                   CWP_STATUS_ON);

  // Create the field to be received:
  const char *recv_field_name = "girafe";
  CWP_Field_create(code_name[0],
                   coupling_name,
                   recv_field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_RECV,
                   CWP_STATUS_ON);

  // Interations :
  // At each iteration the mesh coordinates and the exchanged fields are modified.
  double     *send_field_data = malloc(sizeof(double) * n_vtx);
  double     *recv_field_data = malloc(sizeof(double) * n_vtx);

  const int    itdeb = 1;
  const int    itend = 10;
  const double freq  = 0.20;
  const double ampl  = 0.012;
  const double phi   = 0.1;
  double       ttime = 0.0;
  double       dt    = 0.1;

  double omega = 2.0*acos(-1.0)*freq;

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    // Start time step
    CWP_Time_step_beg(code_name[0],
                      ttime);

    for (int i = 0; i < n_vtx; i++) {
      coords[3 * i + 2]  = ampl * (coords[3 * i]*coords[3 * i]+coords[1 + 3 * i]*coords[1 + 3 * i])*cos(omega*ttime+phi);
      send_field_data[i] = coords[3 * i + 2];
    }


    if (it == itdeb) {

      // Set the mesh vertices coordinates :
      // If the global numbering array is available, it can be given instead
      // of the last NULL argument. If not given, CWIPI will compute it
      // for you in CWP_Mesh_interf_finalize.
      CWP_Mesh_interf_vtx_set(code_name[0],
                              coupling_name,
                              0,
                              n_vtx,
                              coords,
                              NULL);

      // Set the mesh polygons connectivity :
      // Since the mesh elements are triangles, CWP_BLOCK_FACE_TRIA3 could
      // be used instead of CWP_BLOCK_FACE_POLY.
      int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                               coupling_name,
                                               CWP_BLOCK_FACE_POLY);

      // CWP_Mesh_interf_from_faceedge_set is not used here since the
      // mesh is in the form of an face->vertex connectivity. If the mesh
      // is only available with face->edge and edge->vertex connectivity,
      // the faceedge function should be used. The input of one does
      // not distinguish per element type. CWIPI filters that later.
      CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                       coupling_name,
                                       0,
                                       block_id,
                                       n_elt,
                                       elt_vtx_idx,
                                       elt_vtx,
                                       NULL);

      // Finalize mesh :
      CWP_Mesh_interf_finalize(code_name[0],
                               coupling_name);

      // Set the values of the field to be send:
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         send_field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_field_data);

      // Set the values of the field to be received:
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         recv_field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_field_data);

      // Set interpolation property :
      CWP_Spatial_interp_property_set(code_name[0],
                                      coupling_name,
                                      "tolerance",
                                      CWP_DOUBLE,
                                      "0.1");

    }

    // Compute interpolation weights :
    CWP_Spatial_interp_weights_compute(code_name[0],
                                       coupling_name);

    // Exchange field values :
    // If the codes operate a cross exchange, a deadlock could occur if
    // CWP_Field_wait_issend/CWP_Field_wait_irecv calls are mixed with
    // CWP_Field_issend/CWP_Field_irecv calls.
    CWP_Field_issend(code_name[0],
                     coupling_name,
                     send_field_name);

    CWP_Field_irecv(code_name[0],
                    coupling_name,
                    recv_field_name);

    CWP_Field_wait_issend(code_name[0],
                          coupling_name,
                          send_field_name);

    CWP_Field_wait_irecv(code_name[0],
                         coupling_name,
                         recv_field_name);

    // Check interpolation :
    int n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0],
                                                      coupling_name,
                                                      recv_field_name,
                                                      0);

    const int *uncomputed_tgts = NULL;
    if (n_uncomputed_tgts != 0) {
      uncomputed_tgts = CWP_Uncomputed_tgts_get(code_name[0],
                                                coupling_name,
                                                recv_field_name,
                                                0);
    }

    PDM_UNUSED(n_uncomputed_tgts);
    PDM_UNUSED(uncomputed_tgts);

    // End time step
    CWP_Time_step_end(code_name[0]);

  } // end interations

  // Delete field :
  CWP_Field_del(code_name[0],
                coupling_name,
                send_field_name);

  CWP_Field_del(code_name[0],
                coupling_name,
                recv_field_name);

  // Delete Mesh :
  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling_name);

  // free
  free(intra_comm);
  free(code_name);
  free(coupled_code_name);
  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(send_field_data);
  free(recv_field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}
