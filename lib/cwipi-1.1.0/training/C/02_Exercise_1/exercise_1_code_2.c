#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cwp.h"
#include "cwp_priv.h"

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -b           Enable bonu \n\n"
         "  -h           This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *bonus
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-b") == 0) {
      *bonus = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

int
main(int argc, char *argv[]) {

  int bonus = 0;

  _read_args (argc,
              argv,
              &bonus);

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

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

  int n_part = 1;
  const char  *coupling_name     = "code1_code2";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);

  coupled_code_name[0] = "code1";

  CWP_Spatial_interp_t spatial_interp_algorithm = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  if (bonus == 1) {
    spatial_interp_algorithm = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  }

  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 spatial_interp_algorithm,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  int    n_vtx = 11;
  double coords[33] = {0,0,0,  1,0,0,  2,0,0,  3,0,0,  0,1,0,  2,1,0,
                       3,1,0,  1,2,0,  0,3,0,  2,3,0,  3,3,0};
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

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

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  const char *field_name      = "a super fancy field";
  int         n_components    = 1;

  CWP_Dof_location_t location = CWP_DOF_LOCATION_NODE;
  if (bonus == 1) {
    location = CWP_DOF_LOCATION_CELL_CENTER;
  }

  CWP_Field_create(code_name[0],
                   coupling_name,
                   field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   location,
                   CWP_FIELD_EXCH_RECV,
                   CWP_STATUS_ON);

  int size = n_vtx * n_components;
  if (bonus == 1) {
    size = n_elts * n_components;
  }

  double *recv_field_data = malloc(sizeof(double) * size);

  CWP_Field_data_set(code_name[0],
                     coupling_name,
                     field_name,
                     0,
                     CWP_FIELD_MAP_TARGET,
                     recv_field_data);

  CWP_Time_step_beg(code_name[0],
                    0.0);

  CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  CWP_DOUBLE,
                                  "0.001");
  CWP_Spatial_interp_weights_compute(code_name[0],
                                     coupling_name);

  CWP_Field_irecv(code_name[0],
                  coupling_name,
                  field_name);

  CWP_Field_wait_irecv(code_name[0],
                       coupling_name,
                       field_name);

  CWP_Time_step_end(code_name[0]);

  CWP_Field_del(code_name[0],
                coupling_name,
                field_name);

  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  CWP_Cpl_del(code_name[0],
              coupling_name);

  free(recv_field_data);
  free(code_name);
  free(intra_comm);
  free(coupled_code_name);

  CWP_Finalize();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
