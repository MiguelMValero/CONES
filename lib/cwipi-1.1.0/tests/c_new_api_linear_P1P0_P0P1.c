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
#include "cwipi_config.h"
#include "cwp_priv.h"

#include "pdm.h"
#include "pdm_multipart.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_block_to_part.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MAX(a, b) ((a) > (b) ?  (a) : (b))

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
         "  -algo           spatial interpolation algorithm.\n\n"
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
  PDM_g_num_t            all_gn_vtx[],
  int                    all_n_rank[],
  int                    all_n_part[],
  CWP_Spatial_interp_t  *spatial_interp_algo,
  double                *tolerance,
  int                   *n_neighbors,
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
        all_gn_vtx[0] = atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_gn_vtx[1] = atol(argv[i]);
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
    else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_n_part[1] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-algo") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *spatial_interp_algo = (CWP_Spatial_interp_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tol") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *tolerance = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_cls") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_neighbors = atoi(argv[i]);
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
_gen_mesh
(
 const PDM_MPI_Comm         comm,
 const PDM_g_num_t          gn_vtx,
 const int                  n_part,
 const PDM_split_dual_t     part_method,
 int                      **pn_elt,
 int                      **pn_vtx,
 int                     ***pelt_vtx,
 double                  ***pvtx_coord,
 PDM_g_num_t             ***pelt_ln_to_gn,
 PDM_g_num_t             ***pvtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t gn_elt = gn_vtx - 1;

  PDM_g_num_t *distrib_vtx = PDM_compute_uniform_entity_distribution(comm,
                                                                     gn_vtx);

  PDM_g_num_t *distrib_elt = PDM_compute_uniform_entity_distribution(comm,
                                                                     gn_elt);

  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  int dn_elt = (int) (distrib_elt[i_rank+1] - distrib_elt[i_rank]);


  double step = 6.28 / (double) (gn_vtx - 1);

  double *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  for (int i = 0; i < dn_vtx; i++) {
    PDM_g_num_t g = distrib_vtx[i_rank] + i;

    double r = 1 + 0.3*cos(5*g*step);
    double t = g*step + 0.1*sin(5*g*step);

    dvtx_coord[3*i  ] = r*cos(t);
    dvtx_coord[3*i+1] = r*sin(t);
    dvtx_coord[3*i+2] = 0.1*cos(5*g*step);
  }
  free(distrib_vtx);


  PDM_g_num_t *delt_vtx = malloc(sizeof(PDM_g_num_t) * dn_elt * 2);
  for (int i = 0; i < dn_elt; i++) {
    PDM_g_num_t g = distrib_elt[i_rank] + i;

    delt_vtx[2*i  ] = g+1;
    delt_vtx[2*i+1] = g+2;
  }
  free(distrib_elt);

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                        0,
                                        0,
                                        dn_elt,
                                        dn_vtx,
                                        comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord,
                          PDM_OWNERSHIP_KEEP);

  int *delt_vtx_idx = NULL;
  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             delt_vtx,
                             delt_vtx_idx,
                             PDM_OWNERSHIP_KEEP);

  int n_bnd = 2;
  int *dbnd_vtx_idx = PDM_array_zeros_int(n_bnd + 1);
  if (i_rank == 0) {
    dbnd_vtx_idx[1] = 1;
    dbnd_vtx_idx[2] = 2;
  }
  PDM_g_num_t *dbnd_vtx = malloc(sizeof(PDM_g_num_t) * dbnd_vtx_idx[n_bnd]);
  if (i_rank == 0) {
    dbnd_vtx[0] = 1;
    dbnd_vtx[1] = gn_vtx;
  }

  PDM_dmesh_bound_set(dmesh,
                      PDM_BOUND_TYPE_VTX,
                      n_bnd,
                      dbnd_vtx,
                      dbnd_vtx_idx,
                      PDM_OWNERSHIP_KEEP);


  int n_zone = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_set(mpart,
                               0,
                               dmesh);

  PDM_multipart_compute(mpart);


  *pn_vtx        = malloc(sizeof(int          ) * n_part);
  *pvtx_coord    = malloc(sizeof(double      *) * n_part);
  *pvtx_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_elt        = malloc(sizeof(int          ) * n_part);
  *pelt_vtx      = malloc(sizeof(int         *) * n_part);
  *pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);


  for (int ipart = 0; ipart < n_part; ipart++) {

    /* Vertices */
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
                                                       &(*pvtx_ln_to_gn)[ipart],
                                                       PDM_OWNERSHIP_USER);

    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     ipart,
                                     &(*pvtx_coord)[ipart],
                                     PDM_OWNERSHIP_USER);


    /* Elements */
    (*pn_elt)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_EDGE,
                                                       &(*pelt_ln_to_gn)[ipart],
                                                       PDM_OWNERSHIP_USER);

    int *elt_vtx_idx;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &elt_vtx_idx,
                                        &(*pelt_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
  }

  PDM_dmesh_free(dmesh);
  PDM_multipart_free(mpart);
}




/*----------------------------------------------------------------------
 *
 * Main : Linear coupling interface
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
  int                  verbose        = 0;
  int                  swap_codes     = 0;
  PDM_g_num_t          all_gn_vtx[2]  = {100, 50};
  int                  all_n_rank[2]  = {-1, -1};
  int                  all_n_part[2]  = {1, 1};
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  double               tolerance      = 1e-2;
  int                  n_neighbors    = 5;
  int                  visu           = 0;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  _read_args(argc,
             argv,
             &verbose,
             &swap_codes,
             all_gn_vtx,
             all_n_rank,
             all_n_part,
             &spatial_interp,
             &tolerance,
             &n_neighbors,
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

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);
  int            n_part[2];
  PDM_g_num_t    gn_vtx[2];

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      gn_vtx           [n_code] = all_gn_vtx[icode];
      n_part           [n_code] = all_n_part[icode];

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
  const char *cpl_name = "c_new_api_linear";

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_LINEAR,
                   CWP_COMM_PAR_WITH_PART,
                   spatial_interp,
                   1,
                   CWP_DYNAMIC_MESH_STATIC,
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
  int          **pn_vtx        = malloc(sizeof(int          *) * n_code);
  int          **pn_elt        = malloc(sizeof(int          *) * n_code);
  double      ***pvtx_coord    = malloc(sizeof(double      **) * n_code);
  int         ***pelt_vtx      = malloc(sizeof(int         **) * n_code);
  PDM_g_num_t ***pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]);

    _gen_mesh(mesh_comm,
              gn_vtx[icode],
              n_part[icode],
              part_method,
              &pn_elt       [icode],
              &pn_vtx       [icode],
              &pelt_vtx     [icode],
              &pvtx_coord   [icode],
              &pelt_ln_to_gn[icode],
              &pvtx_ln_to_gn[icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_EDGE2);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx       [icode][ipart],
                              pvtx_coord   [icode][ipart],
                              pvtx_ln_to_gn[icode][ipart]);

      CWP_Mesh_interf_block_std_set(code_name[icode],
                                    cpl_name,
                                    ipart,
                                    block_id,
                                    pn_elt       [icode][ipart],
                                    pelt_vtx     [icode][ipart],
                                    pelt_ln_to_gn[icode][ipart]);
    }

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  /* Define fields */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name1 = "all_coords";
  const char *field_name2 = "coordX";



  double ***field1_val = malloc(sizeof(double **) * n_code);
  double ***field2_val = malloc(sizeof(double **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {

    field1_val[icode] = malloc(sizeof(double *) * n_part[icode]);
    field2_val[icode] = malloc(sizeof(double *) * n_part[icode]);

    if (code_id[icode] == 1) {
      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        field1_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart] * 3);
        field2_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
        for (int i = 0; i < 3*pn_vtx[icode][ipart]; i++) {
          field1_val[icode][ipart][i] = pvtx_coord[icode][ipart][i];
        }

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name1,
                           ipart,
                           CWP_FIELD_MAP_SOURCE,
                           field1_val[icode][ipart]);

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name2,
                           ipart,
                           CWP_FIELD_MAP_TARGET,
                           field2_val[icode][ipart]);
      }
    }
    else {
      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       3,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        field1_val[icode][ipart] = malloc(sizeof(double) * pn_elt[icode][ipart] * 3);
        field2_val[icode][ipart] = malloc(sizeof(double) * pn_elt[icode][ipart]);
        for (int i = 0; i < pn_elt[icode][ipart]; i++) {
          int vtx_id0 = pelt_vtx[icode][ipart][2*i  ] - 1;
          int vtx_id1 = pelt_vtx[icode][ipart][2*i+1] - 1;
          field2_val[icode][ipart][i] = 0.5*(pvtx_coord[icode][ipart][3*vtx_id0] + pvtx_coord[icode][ipart][3*vtx_id1]);
        }

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name1,
                           ipart,
                           CWP_FIELD_MAP_TARGET,
                           field1_val[icode][ipart]);

        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name2,
                           ipart,
                           CWP_FIELD_MAP_SOURCE,
                           field2_val[icode][ipart]);
      }
    }
  }


  /* Exchange fields */
  for (int icode = 0; icode < n_code; icode++) {
    char char_param[99];
    sprintf(char_param, "%e", tolerance);
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "tolerance",
                                    CWP_DOUBLE,
                                    char_param);

    sprintf(char_param, "%d", n_neighbors);
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "n_neighbors",
                                    CWP_INT,
                                    char_param);

    CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Field_issend(code_name[icode], cpl_name, field_name1);
      CWP_Field_irecv (code_name[icode], cpl_name, field_name2);
    }
    else {
      CWP_Field_irecv (code_name[icode], cpl_name, field_name1);
      CWP_Field_issend(code_name[icode], cpl_name, field_name2);
    }
  }


  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Field_wait_issend(code_name[icode], cpl_name, field_name1);
      CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name2);
    }
    else {
      CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name1);
      CWP_Field_wait_issend(code_name[icode], cpl_name, field_name2);
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Exchange fields OK\n");
    fflush(stdout);
  }


  /* Check interpolated fields */
  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {

  //   }
  //   else {

  //   }
  // }



  /* Finalize */
  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      free(pvtx_coord   [icode][ipart]);
      free(pvtx_ln_to_gn[icode][ipart]);
      free(pelt_vtx     [icode][ipart]);
      free(pelt_ln_to_gn[icode][ipart]);
      free(field1_val   [icode][ipart]);
      free(field2_val   [icode][ipart]);
    }
    free(pn_vtx       [icode]);
    free(pn_elt       [icode]);
    free(pvtx_coord   [icode]);
    free(pvtx_ln_to_gn[icode]);
    free(pelt_vtx     [icode]);
    free(pelt_ln_to_gn[icode]);
    free(field1_val   [icode]);
    free(field2_val   [icode]);
  }
  free(pn_vtx       );
  free(pn_elt       );
  free(pvtx_coord   );
  free(pvtx_ln_to_gn);
  free(pelt_vtx     );
  free(pelt_ln_to_gn);
  free(field1_val   );
  free(field2_val   );


  for (int icode = 0; icode < n_code; icode++) {
    CWP_Time_step_end(code_name[icode]);
    CWP_Mesh_interf_del(code_name[icode], cpl_name);
    CWP_Cpl_del        (code_name[icode], cpl_name);
  }
  free(code_id);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  CWP_Finalize();

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("End\n");
    fflush(stdout);
  }

  MPI_Finalize();


  return EXIT_SUCCESS;
}
