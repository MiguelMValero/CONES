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

#include "cwipi.h"
#include "cwp.h"
#include "cwipi_config.h"
#include "cwp_priv.h"

#include "grid_mesh.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"

#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))


static void
_read_mesh
(
 const char    *filename,
       int     *n_vtx,
       double **vtx_coord,
       int     *n_face,
       int    **face_vtx_idx,
       int    **face_vtx
)
{
  FILE *f = fopen(filename, "r");

  assert(f != NULL);

  fscanf(f, "%d %d", n_vtx, n_face);

  *vtx_coord = malloc(sizeof(double) * (*n_vtx) * 3);
  for (int i = 0; i < 3*(*n_vtx); i++) {
    fscanf(f, "%lf", *vtx_coord + i);
  }

  *face_vtx_idx = malloc(sizeof(int) * (*n_face + 1));
  for (int i = 0; i <= *n_face; i++) {
    fscanf(f, "%d", *face_vtx_idx + i);
  }

  *face_vtx = malloc(sizeof(int) * (*face_vtx_idx)[*n_face]);
  for (int i = 0; i < (*face_vtx_idx)[*n_face]; i++) {
    fscanf(f, "%d", *face_vtx + i);
  }

  fclose(f);
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P0P0
 *
 *---------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank == 2);


  // Initialize CWIPI

  const char **code_name         = malloc(sizeof(char *));
  const char **coupled_code_name = malloc(sizeof(char *));
  CWP_Status_t is_active_rank    = CWP_STATUS_ON;

  int n_code = 1;
  int n_part = 1;
  const char *all_code_names[2] = {"code1", "code2"};
  int code_id = i_rank;

  code_name        [0] = all_code_names[code_id];
  coupled_code_name[0] = all_code_names[(code_id+1)%2];



  char filename[999];
  sprintf(filename, CWP_MESH_DIR"guinea_pig_%d.dat", code_id+1);

  int     n_vtx        = 0;
  double *vtx_coord    = NULL;
  int     n_face       = 0;
  int    *face_vtx_idx = NULL;
  int    *face_vtx     = NULL;
  _read_mesh(filename,
             &n_vtx,
             &vtx_coord,
             &n_face,
             &face_vtx_idx,
             &face_vtx);


  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm));

  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  // Create coupling
  const char *cpl_name = "c_new_api_guinea_pig";
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_DIST_CLOUD_SURF;

  CWP_Cpl_create(code_name[0],                                  // Code name
                 cpl_name,                                      // Coupling id
                 coupled_code_name[0],                          // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,                        // Coupling type
                 spatial_interp,
                 n_part,                                        // Number of partitions
                 CWP_DYNAMIC_MESH_STATIC,                       // Mesh displacement type
                 CWP_TIME_EXCH_USER_CONTROLLED);                // Postprocessing frequency


  CWP_Visu_set(code_name[0],            // Code name
               cpl_name,                // Coupling id
               1,                       // Postprocessing frequency
               CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
               "text");                 // Postprocessing option


  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           cpl_name,
                                           CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_vtx_set(code_name[0],
                          cpl_name,
                          0,
                          n_vtx,
                          vtx_coord,
                          NULL);


  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   cpl_name,
                                   0,
                                   block_id,
                                   n_face,
                                   face_vtx_idx,
                                   face_vtx,
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0], cpl_name);


  // Create and set fields
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name1 = "field1";

  double *send_val = malloc(sizeof(double) * n_vtx);
  double *recv_val = malloc(sizeof(double) * n_vtx);

  if (code_id == 0) {
    for (int i = 0; i < n_vtx; i++) {
      send_val[i] = cos(0.5*vtx_coord[3*i]);
    }
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val);
  }
  else {
    for (int i = 0; i < n_vtx; i++) {
      send_val[i] = cos(0.5*vtx_coord[3*i+1]);
    }

    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0],
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       send_val);
  }


  CWP_Spatial_interp_property_set(code_name[0], cpl_name, "tolerance", CWP_DOUBLE, "0.1");

  CWP_Spatial_interp_weights_compute(code_name[0], cpl_name);



  if (code_id == 0) {
    CWP_Field_issend(code_name[0], cpl_name, field_name1);
  }
  else {
    CWP_Field_irecv (code_name[0], cpl_name, field_name1);
  }


  if (code_id == 0) {
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name1);
  }
  else {
    CWP_Field_wait_irecv (code_name[0], cpl_name, field_name1);
  }

  CWP_Time_step_end(code_name[0]);

  CWP_Mesh_interf_del(code_name[0], cpl_name);

  CWP_Cpl_del(code_name[0], cpl_name);

  // const char   *field_name[]   = {field_name1};//, field_name2};
  // const double *field_value[2];
  // if (code_id == 0) {
  //   field_value[0] = send_val;
  //   // field_value[1] = recv_val;
  // }
  // else {
  //   field_value[0] = send_val;
  //   // field_value[1] = recv_val;
  // }

  // sprintf(filename, "result_guinea_pig_%d.vtk", code_id+1);
  // PDM_vtk_write_polydata_field(filename,
  //                              n_vtx,
  //                              vtx_coord,
  //                              NULL,
  //                              n_face,
  //                              face_vtx_idx,
  //                              face_vtx,
  //                              NULL,
  //                              NULL,
  //                              NULL,
  //                              field_name[0],
  //                              field_value[0]);



  free(vtx_coord);
  free(face_vtx_idx);
  free(face_vtx);

  free(send_val);
  free(recv_val);

  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

