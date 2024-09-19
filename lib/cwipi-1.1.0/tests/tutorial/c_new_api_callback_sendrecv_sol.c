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
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_logging.h"
#include "pdm_generate_mesh.h"

/*----------------------------------------------------------------------
 *
 * User Interpolation function (callback)
 *
 *---------------------------------------------------------------------*/

static void
my_interpolation
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
       int             i_part,
       double         *buffer_in,
       double         *buffer_out
 )
  {
   CWP_Spatial_interp_t spatial_interp_algorithm = CWP_Cpl_spatial_interp_algo_get(local_code_name,
                                                                                   cpl_id);
   CWP_Field_storage_t  storage                  = CWP_Field_storage_get(local_code_name,
                                                                         cpl_id,
                                                                         field_id);

   if (spatial_interp_algorithm == CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES) {
    if (storage == CWP_FIELD_STORAGE_INTERLACED) {

      // Get interpolation information :
      // Interpolation is done on the target side for this KNN algorithm
      int  n_elt_tgt        = -1;
      int  n_referenced_tgt = -1;
      int *referenced_tgt        = NULL;
      int *tgt_come_from_src_ids = NULL;
      CWP_Field_tgt_data_properties_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part,
                                        &n_elt_tgt,
                                        &n_referenced_tgt,
                                        &referenced_tgt,
                                        &tgt_come_from_src_ids);

      // Get point location information :
      // For each target, this array gives the squared distance to its nearest sources
      double *d2 = NULL;
      CWP_Field_nearest_neighbors_distances_get(local_code_name,
                                                cpl_id,
                                                field_id,
                                                i_part,
                                                &d2);

      for (int i = 0; i < n_referenced_tgt; i++) {

        int    j_nearest  = tgt_come_from_src_ids[i];
        double d2_nearest = d2[j_nearest];

        for (int j = tgt_come_from_src_ids[i]; j < tgt_come_from_src_ids[i+1]; j++) {

          if (d2[j] < d2_nearest) {
            j_nearest  = j;
            d2_nearest = d2[j];
          }

        }

        buffer_out[i] = buffer_in[j_nearest];

      }

    }
  }
}

/*----------------------------------------------------------------------
 *
 * Main : advanced test : Callback (codeC)
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
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

  code_name[0]      = "codeC";

  printf("C : %d/%d je suis là\n", i_rank, n_rank);
  fflush(stdout);

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  printf("C : %d/%d CWP_Init OK\n", i_rank, n_rank);
  fflush(stdout);

  // Create the coupling with code 2 :
  int n_part = 1;
  const char  *coupling_name     = "coupling_C_Python";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "codePython";
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  // MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("C : CWP_Cpl_create OK\n");
    fflush(stdout);
  }

  // Create mesh :
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

  // Set mesh :
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling_name,
                                           CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling_name,
                                   0,
                                   block_id,
                                   n_elt,
                                   elt_vtx_idx,
                                   elt_vtx,
                                   NULL);

  // Set user targets :
  // For simplicity we locate them at nodes here but they could be located anywhere
  CWP_User_tgt_pts_set(code_name[0],
                       coupling_name,
                       0,
                       n_vtx,
                       coords,
                       NULL);

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  // MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("C : CWP_Mesh_interf_finalize OK\n");
    fflush(stdout);
  }

  // Create and set field :
  const char *field_name   = "coord_x";
  int         n_components = 1;
  CWP_Field_create(code_name[0],
                   coupling_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_USER,
                   CWP_FIELD_EXCH_SENDRECV,
                   CWP_STATUS_ON);

  // Start time step
  CWP_Time_step_beg(code_name[0],
                    0.);

  double *send_field_data = malloc(sizeof(double) * n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    send_field_data[i] = coords[3*i];
  }
  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     send_field_data);

  double *recv_field_data = malloc(sizeof(double) * n_vtx);
  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_TARGET,
                     recv_field_data);

  // Set user interpolation function :
  CWP_Field_interp_function_set(code_name[0], coupling_name, field_name, my_interpolation);

  // Set interpolation property and compute weights
  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "n_neighbors",
                                  CWP_INT,
                                  "3");

  CWP_Spatial_interp_weights_compute(code_name[0],
                                     coupling_name);

  // Exchange
  CWP_Field_issend(code_name[0], coupling_name, field_name);
  CWP_Field_irecv(code_name[0], coupling_name, field_name);

  CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
  CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);

  // End time step
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
  free(code_name);
  free(intra_comm);
  free(coupled_code_name);
  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(send_field_data);
  free(recv_field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  printf("C rank %d FINISHED :D\n", i_rank);

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;

}
