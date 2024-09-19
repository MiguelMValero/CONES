#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm.h"
#include "pdm_generate_mesh.h"

int
main(int argc, char *argv[]) {

  // Initialize MPI
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

  code_name[0] = "code2";

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  // Create the coupling :
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

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  const char *field_name      = "a super fancy field";
  int         n_components    = 1;

  CWP_Field_create(code_name[0],
                   coupling_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_SEND,
                   CWP_STATUS_ON);

  double *field_data = malloc(sizeof(double) * n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    field_data[i] = coords[3 * i];
  }

  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     field_data);

  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  CWP_DOUBLE,
                                  "0.001");

  const int    itdeb = 1;
  const int    itend = 10;
  double       ttime = 0.0;
  double       dt    = 0.1;

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    // Start time step
    CWP_Time_step_beg(code_name[0],
                      ttime);

    CWP_Spatial_interp_weights_compute(code_name[0],
                                       coupling_name);

    CWP_Field_issend(code_name[0],
                     coupling_name,
                     field_name);

    CWP_Field_wait_issend(code_name[0],
                          coupling_name,
                          field_name);

    CWP_Time_step_end(code_name[0]);

  } // end interations

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
  free(intra_comm);
  free(code_name);
  free(coupled_code_name);
  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}

