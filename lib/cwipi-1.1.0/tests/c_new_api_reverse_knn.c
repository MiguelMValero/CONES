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

#include "cwp_priv.h"
#include "cwp.h"

#include "pdm.h"
#include "pdm_multipart.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_block_to_part.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_generate_mesh.h"


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -v              verbose.\n\n"
         "  -n_rank1        number of MPI ranks for code1.\n\n"
         "  -n_rank2        number of MPI ranks for code2.\n\n"
         "  -v              verbose.\n\n"
         "  -n1             square root of number of vertices for code1.\n\n"
         "  -n2             square root of number of vertices for code2.\n\n"
         "  -swap_codes     swap rank order of code1 and 2.\n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *verbose,
  int                   *swap_codes,
  PDM_g_num_t            all_n_vtx_seg[],
  int                    all_n_rank[],
  int                   *visu
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-swap_codes") == 0) {
      *swap_codes = 1;
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_vtx_seg[0] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_vtx_seg[1] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_rank[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_rank[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_my_interpolation
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
 int                   i_part,
 double               *buffer_in,
 double               *buffer_out
)
{
  int  n_elt_tgt             = 0;
  int  n_referenced_tgt      = 0;
  int *referenced_tgt        = NULL;
  int *tgt_come_from_src_idx = NULL;
  CWP_Field_tgt_data_properties_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part,
                                    &n_elt_tgt,
                                    &n_referenced_tgt,
                                    &referenced_tgt,
                                    &tgt_come_from_src_idx);

  int n_components = CWP_Field_n_components_get(local_code_name,
                                                cpl_id,
                                                field_id);

  for (int i = 0; i < n_elt_tgt * n_components; i++) {
    buffer_out[i] = 0.;
  }

  for (int i = 0; i < n_referenced_tgt; i++) {
    int id = referenced_tgt[i] - 1;
    for (int j = tgt_come_from_src_idx[i]; j < tgt_come_from_src_idx[i+1]; j++) {
      for (int k = 0; k < n_components; k++) {
        buffer_out[n_components*id + k] += buffer_in[n_components*j + k];
      }
    }
    // PDM_log_trace_array_double(buffer_in + n_components*tgt_come_from_src_idx[i],
    //                            n_components*(tgt_come_from_src_idx[i+1] - tgt_come_from_src_idx[i]),
    //                            "buffer_in  : ");
    // PDM_log_trace_array_double(buffer_out + n_components*id,
    //                            n_components,
    //                            "buffer_out : ");
    // log_trace("\n");
  }
}


static void _deform
(
       double *x,
       double *y,
       double *z,
 const double  t
 )
{
  CWP_UNUSED(y);
  *z = 0.2*cos(*x + t);
}



/*----------------------------------------------------------------------
 *
 * Main
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
  int         verbose          = 0;
  int         swap_codes       = 0;
  PDM_g_num_t all_n_vtx_seg[2] = {10, 5};
  int         all_n_rank   [2] = {-1, -1};
  int         visu             = 0;

  _read_args(argc,
             argv,
             &verbose,
             &swap_codes,
             all_n_vtx_seg,
             all_n_rank,
             &visu);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  for (int i = 0; i < 2; i++) {
    if (all_n_rank[i] <= 0) {
      all_n_rank[i] = n_rank;
    }
  }

  const char *all_code_names[2] = {"code1", "code2"};
  int has_code[2] = {0, 0};


  has_code[0] = i_rank <  all_n_rank[0];
  has_code[1] = i_rank >= n_rank - all_n_rank[1];

  if (swap_codes) {
    int tmp = has_code[0];
    has_code[0] = has_code[1];
    has_code[1] = tmp;
  }

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);
  PDM_g_num_t   *n_vtx_seg         = malloc(sizeof(PDM_g_num_t ) * n_code);

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_vtx_seg        [n_code] = all_n_vtx_seg[icode];

      if (verbose) {
        log_trace("Running %s, coupled with %s\n",
                  code_name[n_code], coupled_code_name[n_code]);
      }
      n_code++;
    }
  }

  // Set up
  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }

  /* Create coupling */
  const char *cpl_name = "c_new_api_reverse_knn";

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES,
                   1,
                   CWP_DYNAMIC_MESH_DEFORMABLE,//CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);
  }


  for (int icode = 0; icode < n_code; icode++) {
    CWP_Visu_set(code_name[icode],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }


  /* Define interface mesh */
  int     *n_vtx       = malloc(sizeof(int     ) * n_code);
  int     *n_elt       = malloc(sizeof(int     ) * n_code);
  double **vtx_coord   = malloc(sizeof(double *) * n_code);
  int    **elt_vtx_idx = malloc(sizeof(int    *) * n_code);
  int    **elt_vtx     = malloc(sizeof(int    *) * n_code);
  int      n_tgt       = 0;
  double  *tgt_coord   = NULL;

  for (int icode = 0; icode < n_code; icode++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]);

    PDM_generate_mesh_rectangle_simplified(mesh_comm,
                                           n_vtx_seg   [icode],
                                           &n_vtx      [icode],
                                           &n_elt      [icode],
                                           &vtx_coord  [icode],
                                           &elt_vtx_idx[icode],
                                           &elt_vtx    [icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);

    CWP_Mesh_interf_vtx_set(code_name[icode],
                            cpl_name,
                            0,
                            n_vtx    [icode],
                            vtx_coord[icode],
                            NULL);

    if (code_id[icode] == 2) {
      tgt_coord = malloc(sizeof(double) * (n_vtx[icode]/2) * 3);
      for (int i = 1; i < n_vtx[icode]; i += 2) {
        memcpy(&tgt_coord[3*n_tgt], &vtx_coord[icode][3*i], sizeof(double) * 3);
        n_tgt++;
      }

      CWP_User_tgt_pts_set(code_name[icode],
                           cpl_name,
                           0,
                           n_tgt,
                           tgt_coord,
                           NULL);
    }

    CWP_Mesh_interf_f_poly_block_set(code_name[icode],
                                     cpl_name,
                                     0,
                                     block_id,
                                     n_elt      [icode],
                                     elt_vtx_idx[icode],
                                     elt_vtx    [icode],
                                     NULL);

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }

  /* Field */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "my_field";

  double *field_ptr = NULL;
  double *send_val  = NULL;
  double *recv_val  = NULL;

  for (int icode = 0; icode < n_code; icode++) {

    CWP_Field_exch_t   exch_type;
    CWP_Field_map_t    map_type;
    CWP_Dof_location_t dof_location;
    if (code_id[icode] == 1) {
      send_val = malloc(sizeof(double *) * n_elt[icode]);
      for (int i = 0; i < n_elt[icode]; i++) {
        // send_val[i] = (double) rand() / (double) RAND_MAX;
        send_val[i] = 0;
        for (int j = elt_vtx_idx[icode][i]; j < elt_vtx_idx[icode][i+1]; j++) {
          int vtx_id = elt_vtx[icode][j] - 1;
          send_val[i] += vtx_coord[icode][3*vtx_id];
        }
        send_val[i] /= (elt_vtx_idx[icode][i+1] - elt_vtx_idx[icode][i]);
      }

      exch_type    = CWP_FIELD_EXCH_SEND;
      map_type     = CWP_FIELD_MAP_SOURCE;
      field_ptr    = send_val;
      dof_location = CWP_DOF_LOCATION_CELL_CENTER;
    }
    else {
      recv_val = malloc(sizeof(double *) * n_vtx[icode]);

      exch_type    = CWP_FIELD_EXCH_RECV;
      map_type     = CWP_FIELD_MAP_TARGET;
      field_ptr    = recv_val;
      dof_location = CWP_DOF_LOCATION_USER;//CWP_DOF_LOCATION_NODE;
    }

    CWP_Field_create(code_name[icode],
                     cpl_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     dof_location,
                     exch_type,
                     visu_status);

    CWP_Field_data_set(code_name[icode],
                       cpl_name,
                       field_name,
                       0,
                       map_type,
                       field_ptr);

    CWP_Field_interp_function_set(code_name[icode],
                                  cpl_name,
                                  field_name,
                                  _my_interpolation);
  }


  /* Exchange */
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "n_neighbors",
                                    CWP_INT,
                                    "1");
  }

  double recv_time = 0.;
  for (int step = 0; step < 15; step++) {

    for (int icode = 0; icode < n_code; icode++) {
      CWP_Time_step_beg(code_name[icode],
                        recv_time);
    }

    if (i_rank == 0) {
      printf("\n  Step %d\n", step);
    }
    if (verbose) {
      log_trace("\n  Step %d\n", step);
    }
    for (int icode = 0; icode < n_code; icode++) {

      for (int i = 0; i < n_vtx[icode]; i++) {
        _deform(&vtx_coord[icode][3*i  ],
                &vtx_coord[icode][3*i+1],
                &vtx_coord[icode][3*i+2],
                recv_time);
      }

      if (code_id[icode] == 2) {
        for (int i = 0; i < n_tgt; i++) {
          _deform(&tgt_coord[3*i  ],
                  &tgt_coord[3*i+1],
                  &tgt_coord[3*i+2],
                  recv_time);
        }
      }
    }

    // Separate loops to avoid deadlock if multiple codes on same MPI rank
    for (int icode = 0; icode < n_code; icode++) {
      CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
    }

    for (int icode = 0; icode < n_code; icode++) {
      if (code_id[icode] == 1) {
        CWP_Field_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_irecv (code_name[icode], cpl_name, field_name);
      }
    }


    for (int icode = 0; icode < n_code; icode++) {
      if (code_id[icode] == 1) {
        CWP_Field_wait_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name);
      }
    }

    recv_time += 0.2;


    /* Check conservation */
    double l_integral[2] = {0., 0.};

    for (int icode = 0; icode < n_code; icode++) {
      int n = 0;
      if (code_id[icode] == 1) {
        n = n_elt[icode];
        field_ptr = send_val;
      }
      else {
        n = n_vtx[icode];
        field_ptr = recv_val;
      }

      for (int i = 0; i < n; i++) {
        l_integral[code_id[icode]-1] += field_ptr[i];
      }
    }

    double g_integral[2];
    MPI_Allreduce(l_integral, g_integral, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (i_rank == 0) {
      printf("g_integral = %20.16e / %20.16e, relative diff = %e\n",
             g_integral[0], g_integral[1], fabs(g_integral[0] - g_integral[1])/fabs(g_integral[1]));
    }

    for (int icode = 0; icode < n_code; icode++) {
      CWP_Time_step_end(code_name[icode]);
    }

  }


  /* Finalize */
  for (int icode = 0; icode < n_code; icode++) {

    if (code_id[icode] == 1) {
      free(send_val);
    }
    else {
      free(recv_val);
      free(tgt_coord);
    }

    free(vtx_coord  [icode]);
    free(elt_vtx_idx[icode]);
    free(elt_vtx    [icode]);
  }

  free(n_vtx      );
  free(n_elt      );
  free(vtx_coord  );
  free(elt_vtx_idx);
  free(elt_vtx    );

  for (int icode = 0; icode < n_code; icode++) {

    CWP_Mesh_interf_del(code_name[icode], cpl_name);

    CWP_Cpl_del(code_name[icode], cpl_name);
  }
  free(code_id);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);
  free(n_vtx_seg);

  CWP_Finalize();

  MPI_Finalize();

  return EXIT_SUCCESS;


}
