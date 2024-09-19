/*
  This file is part of the CWIPI library.

  Copyright (C) 2023w  ONERA

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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cwp.h"
#include "cwipi_config.h"
#include "cwp_priv.h"
#include "client_server/client.h"

#include "pdm.h"
#include "pdm_multipart.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_printf.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_part_to_block.h"

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MAX(a, b) ((a) > (b) ?  (a) : (b))

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -mc1      <> Mesh file of code1.\n\n"
     "  -mc2      <> Mesh file of code2.\n\n"
     "  -n_rankc1 <> Number of ranks for code1.\n\n"
     "  -n_rankc2 <> Number of ranks for code2.\n\n"
     "  -visu        Activate visualization.\n\n"
     "  -h     This message.\n\n");

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
 int            argc,
 char         **argv,
 char          *mesh_filenames[],
 int            code_n_rank[],
 int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-mc1") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        mesh_filenames[0] = argv[i];
      }
    }

    else if (strcmp(argv[i], "-mc2") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        mesh_filenames[1] = argv[i];
      }
    }
    else if (strcmp(argv[i], "-n_rankc1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rankc2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[1] = atoi(argv[i]);
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

/*----------------------------------------------------------------------
 *
 * Field
 *
 *---------------------------------------------------------------------*/

static double
_Franke
(
 const double x,
 const double y,
 const double z
 )
{
  return 0.75*exp(-(pow(9*x-2,2) + pow(9*y-2,2) + pow(9*z-2,2))/4)
  +      0.75*exp(-(pow(9*x+1,2)/49 + (9*y+1)/10 + (9*z+1)/10))
  +      0.5 *exp(-(pow(9*x-7,2) + pow(9*y-3,2) + pow(9*y-5,2))/4)
  -      0.2 *exp(-(pow(9*x-4,2) + pow(9*y-7,2) + pow(9*z-5,2)));
}


static double
_paper
(
 const double x,
 const double y,
 const double z
 )
{
  return 0.78 + cos(10*(x+y+z));
}

/*----------------------------------------------------------------------
 *
 * Read mesh file and partition it
 *
 *---------------------------------------------------------------------*/

CWP_GCC_SUPPRESS_WARNING("-Wcast-qual")
static void
_gen_part_data
(
 const char                *filename,
 const PDM_MPI_Comm         comm,
 const int                  n_part,
 const PDM_split_dual_t     part_method,
       int                **pn_face,
       int                **pn_vtx,
       int               ***pface_vtx_idx,
       int               ***pface_vtx,
       double            ***pvtx_coord,
       PDM_g_num_t       ***pface_ln_to_gn,
       PDM_g_num_t       ***pvtx_ln_to_gn,
       int                 *n_vtx_field,
       char              ***vtx_field_name,
       PDM_data_t         **vtx_field_type,
       int                **vtx_field_stride,
       double           ****pvtx_field_value

 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  void       **dvtx_field_value;
  int          n_elt_field;
  char       **elt_field_name;
  PDM_data_t  *elt_field_type;
  int         *elt_field_stride;
  void       **delt_field_value;

  PDM_dmesh_nodal_t *dmn = PDM_vtk_read_to_dmesh_nodal(comm,
                                                       filename,
                                                       n_vtx_field,
                                                       vtx_field_name,
                                                       vtx_field_type,
                                                       vtx_field_stride,
                                                       &dvtx_field_value,
                                                       &n_elt_field,
                                                       &elt_field_name,
                                                       &elt_field_type,
                                                       &elt_field_stride,
                                                       &delt_field_value);

  int n_zone = 1;
  int n_part_zones = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                &n_part_zones,
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

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  PDM_multipart_compute(mpart);


  *pn_vtx         = malloc(sizeof(int          ) * n_part);
  *pvtx_coord     = malloc(sizeof(double      *) * n_part);
  *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_face        = malloc(sizeof(int          ) * n_part);
  *pface_vtx_idx  = malloc(sizeof(int         *) * n_part);
  *pface_vtx      = malloc(sizeof(int         *) * n_part);
  *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

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


    /* Faces */
    (*pn_face)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_FACE,
                                                       &(*pface_ln_to_gn)[ipart],
                                                       PDM_OWNERSHIP_USER);

    int *_face_vtx;
    int *_face_vtx_idx;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &_face_vtx_idx,
                                        &_face_vtx,
                                        PDM_OWNERSHIP_KEEP);

    if (_face_vtx != NULL) {
      (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      memcpy((*pface_vtx_idx)[ipart], _face_vtx_idx, sizeof(int) * ((*pn_face)[ipart]+1));

      (*pface_vtx)[ipart] = malloc(sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
      memcpy((*pface_vtx)[ipart], _face_vtx, sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
    }

    else {
      int *_face_edge;
      int *_face_edge_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &_face_edge_idx,
                                          &_face_edge,
                                          PDM_OWNERSHIP_KEEP);

      int *_edge_vtx;
      int *_edge_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &_edge_vtx_idx,
                                          &_edge_vtx,
                                          PDM_OWNERSHIP_KEEP);

      (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      memcpy((*pface_vtx_idx)[ipart], _face_edge_idx, sizeof(int) * ((*pn_face)[ipart]+1));


      PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &(*pface_vtx)[ipart]);
    }

  }
  PDM_multipart_free(mpart);

  /* Exchange fields to partitions */
  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);

  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);

  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_vtx,
                               (const PDM_g_num_t **) *pvtx_ln_to_gn,
                                                      *pn_vtx,
                                                      n_part,
                                                      comm);

  *pvtx_field_value = malloc(sizeof(double **) * (*n_vtx_field));
  for (int ifield = 0; ifield < *n_vtx_field; ifield++) {

    double *_dvtx_field_value = NULL;
    if ((*vtx_field_type)[ifield] == PDM_INT) {
      _dvtx_field_value = malloc(sizeof(double) * dn_vtx * (*vtx_field_stride)[ifield]);
      int *dfv = (int *) dvtx_field_value[ifield];
      for (int i = 0; i < dn_vtx * (*vtx_field_stride)[ifield]; i++) {
        _dvtx_field_value[i] = (double) dfv[i];
      }
    }
    else {
      _dvtx_field_value = (double *) dvtx_field_value[ifield];
    }

    PDM_block_to_part_exch(btp,
                           sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &(*vtx_field_stride)[ifield],
                (void   *) _dvtx_field_value,
                           NULL,
                (void ***) &(*pvtx_field_value)[ifield]);

    if ((*vtx_field_type)[ifield] == PDM_INT) {
      free(_dvtx_field_value);
    }

    free(dvtx_field_value[ifield]);
  }

  PDM_block_to_part_free(btp);

  PDM_DMesh_nodal_free(dmn);

  if (dvtx_field_value != NULL) {
    free(dvtx_field_value);
  }

  if (*n_vtx_field > 0) {
    for (int ifield = 0; ifield < *n_vtx_field; ifield++) {
      free((*vtx_field_name)[ifield]);
      for (int ipart = 0; ipart < n_part; ipart++) {
        free((*pvtx_field_value)[ifield][ipart]);
      }
      free((*pvtx_field_value)[ifield]);
    }
    free(*vtx_field_name  );
    free(*vtx_field_type  );
    free(*vtx_field_stride);
    free(*pvtx_field_value);
  }

  *n_vtx_field = 1;

  *vtx_field_type = malloc(sizeof(PDM_data_t) * 1);
  (*vtx_field_type)[0] = PDM_DOUBLE;

  *vtx_field_stride = malloc(sizeof(int) * 1);
  (*vtx_field_stride)[0] = 1;

  *vtx_field_name = malloc(sizeof(char *) * 1);
  (*vtx_field_name)[0] = malloc(sizeof(char) * 6);
  sprintf((*vtx_field_name)[0], "field");

  *pvtx_field_value = malloc(sizeof(double **) * 1);
  (*pvtx_field_value)[0] = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pvtx_field_value)[0][ipart] = malloc(sizeof(double) * (*pn_vtx)[ipart]);
    for (int i = 0; i < (*pn_vtx)[ipart]; i++) {
      if (0) {
        (*pvtx_field_value)[0][ipart][i] = _Franke((*pvtx_coord)[ipart][3*i  ],
                                                   (*pvtx_coord)[ipart][3*i+1],
                                                   (*pvtx_coord)[ipart][3*i+2]);
      }
      else {
        (*pvtx_field_value)[0][ipart][i] = _paper((*pvtx_coord)[ipart][3*i  ],
                                                  (*pvtx_coord)[ipart][3*i+1],
                                                  (*pvtx_coord)[ipart][3*i+2]);
      }
    }
  }
}

/*----------------------------------------------------------------------
 *
 * Main : preCICE wind turbine blade test case
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  char   *mesh_filenames[2] = {NULL};
  int     code_n_rank   [2] = {1, 1};
  int     visu              = 0;

  _read_args (argc,
              argv,
              mesh_filenames,
              code_n_rank,
              &visu);

  if (mesh_filenames[0] == NULL) {
    mesh_filenames[0] = (char *) CWP_MESH_DIR"blade_0.01.vtk";
  }

  if (mesh_filenames[1] == NULL) {
    mesh_filenames[1] = (char *) CWP_MESH_DIR"blade_0.006.vtk";
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(code_n_rank[0] + code_n_rank[1] == n_rank);

  // Configuration file
  char *config = NULL;
  if (i_rank == 0) {
    config = (char *) "client_new_api_wind_turbine_blade_o/code1/cwp_config_srv.txt";
  }
  else {
    config = (char *) "client_new_api_wind_turbine_blade_o/code2/cwp_config_srv.txt";
  }

  // Launch server
  if (i_rank == 0) {
    system("mkdir -p client_new_api_wind_turbine_blade_o/code1");
    system("mkdir -p client_new_api_wind_turbine_blade_o/code2");
    system("rm -f ./client_new_api_wind_turbine_blade_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_new_api_wind_turbine_blade_o/code2/cwp_config_srv.txt");

    char str[999];
    sprintf(str, "mpirun -n %d cwp_server -cn code1 -p %d %d -c \"client_new_api_wind_turbine_blade_o/code1/cwp_config_srv.txt\" \
                  : -n %d  cwp_server -cn code2 -p %d %d -c \"client_new_api_wind_turbine_blade_o/code2/cwp_config_srv.txt\" &",
                  code_n_rank[0], 54100, 54100 + code_n_rank[0] - 1, code_n_rank[1], 54100 + code_n_rank[0],  54100 + code_n_rank[0] + code_n_rank[1] - 1);
    system(str);
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(5);

  // CWP_Init
  int code1 = i_rank < code_n_rank[0];
  const char *code_name;

  if (code1) {
    code_name = "code1";
  } else {
    code_name = "code2";
  }

  CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

  // Client intra-communicator
  MPI_Comm intra_comm;
  MPI_Comm_split(comm, code1, i_rank, &intra_comm);

  CWP_client_Init(intra_comm,
                  config,
                  code_name,
                  is_coupled_rank);

  // EXIT_SUCCESS ?
  int exit_check = 0;

  // CWP_Cpl_create
  const char *cpl_name = "c_new_api_wind_turbine_blade";
  const char *coupled_code_name;

  if (code1) {
    coupled_code_name = "code2";
  } else {
    coupled_code_name = "code1";
  }

  int n_part = 1;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_client_Cpl_create(code_name,
                        cpl_name,
                        coupled_code_name,
                        CWP_INTERFACE_SURFACE,
                        CWP_COMM_PAR_WITH_PART,
                        spatial_interp,
                        n_part,
                        CWP_DYNAMIC_MESH_STATIC,
                        CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_client_Visu_set(code_name,
                      cpl_name,
                      1,
                      CWP_VISU_FORMAT_ENSIGHT,
                      "text");

  // Read mesh
  int           *pn_face          = NULL;
  int           *pn_vtx           = NULL;
  int          **pface_vtx_idx    = NULL;
  int          **pface_vtx        = NULL;
  double       **pvtx_coord       = NULL;
  PDM_g_num_t  **pface_ln_to_gn   = NULL;
  PDM_g_num_t  **pvtx_ln_to_gn    = NULL;
  int            n_vtx_field      = 0;
  char         **vtx_field_name   = NULL;
  PDM_data_t    *vtx_field_type   = NULL;
  int           *vtx_field_stride = NULL;
  double      ***pvtx_field_value = NULL;

  PDM_MPI_Comm pdm_intra_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm);
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  char *filename = NULL;

  if (code1) {
    filename = mesh_filenames[0];
  } else {
    filename = mesh_filenames[1];
  }

  _gen_part_data(filename,
                 pdm_intra_comm,
                 n_part,
                 part_method,
                 &pn_face,
                 &pn_vtx,
                 &pface_vtx_idx,
                 &pface_vtx,
                 &pvtx_coord,
                 &pface_ln_to_gn,
                 &pvtx_ln_to_gn,
                 &n_vtx_field,
                 &vtx_field_name,
                 &vtx_field_type,
                 &vtx_field_stride,
                 &pvtx_field_value);

  // Set mesh
  int block_id = CWP_client_Mesh_interf_block_add(code_name,
                                                  cpl_name,
                                                  CWP_BLOCK_FACE_POLY);

  for (int i_part = 0; i_part < n_part; i_part++) {
    CWP_client_Mesh_interf_vtx_set(code_name,
                                   cpl_name,
                                   i_part,
                                   pn_vtx       [i_part],
                                   pvtx_coord   [i_part],
                                   pvtx_ln_to_gn[i_part]);

    CWP_client_Mesh_interf_f_poly_block_set(code_name,
                                            cpl_name,
                                            i_part,
                                            block_id,
                                            pn_face       [i_part],
                                            pface_vtx_idx [i_part],
                                            pface_vtx     [i_part],
                                            pface_ln_to_gn[i_part]);
  }

  CWP_client_Mesh_interf_finalize(code_name,
                                  cpl_name);

  // Create field
  CWP_Field_exch_t exch_type;
  CWP_Field_map_t  map_type;
  double **field_ptr = NULL;
  double ** recv_val = NULL;

  if (code1) {
    exch_type = CWP_FIELD_EXCH_SEND;
    map_type  = CWP_FIELD_MAP_SOURCE;
    field_ptr = pvtx_field_value[0];
  } else {
    recv_val = malloc(sizeof(double *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      recv_val[i_part] = malloc(sizeof(double) * pn_vtx[i_part]);
    }

    exch_type = CWP_FIELD_EXCH_RECV;
    map_type  = CWP_FIELD_MAP_TARGET;
    field_ptr = recv_val;
  }

  CWP_client_Field_create(code_name,
                          cpl_name,
                          vtx_field_name[0],
                          CWP_DOUBLE,
                          CWP_FIELD_STORAGE_INTERLACED,
                          vtx_field_stride[0],
                          CWP_DOF_LOCATION_NODE,
                          exch_type,
                          CWP_STATUS_ON);

  for (int i_part = 0; i_part < n_part; i_part++) {
    CWP_client_Field_data_set(code_name,
                              cpl_name,
                              vtx_field_name[0],
                              i_part,
                              map_type,
                              field_ptr[i_part]);
  }

  // Compute weights
  CWP_client_Spatial_interp_property_set(code_name,
                                         cpl_name,
                                         "tolerance",
                                         CWP_DOUBLE,
                                         "0.1");

  CWP_client_Spatial_interp_property_set(code_name,
                                         cpl_name,
                                         "n_neighbors",
                                         CWP_INT,
                                         "1");

  CWP_client_Spatial_interp_weights_compute(code_name,
                                     cpl_name);

  // Exchange field
  if (code1) {
    CWP_client_Field_issend(code_name,
                            cpl_name,
                            vtx_field_name[0]);
  } else {
    CWP_client_Field_irecv (code_name,
                            cpl_name,
                            vtx_field_name[0]);
  }

  if (code1) {
    CWP_client_Field_wait_issend(code_name,
                                 cpl_name,
                                 vtx_field_name[0]);
  } else {
    CWP_client_Field_wait_irecv (code_name,
                                 cpl_name,
                                 vtx_field_name[0]);
  }

  // Check interpolation error
  double linf_error = 0.;
  double l2_error   = 0.;
  double **pvtx_error = NULL;
  if (!code1) {
    int i_rank_intra;
    int n_rank_intra;
    MPI_Comm_rank(intra_comm, &i_rank_intra);
    MPI_Comm_size(intra_comm, &n_rank_intra);

    pvtx_error = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      pvtx_error[ipart] = malloc(sizeof(double) * pn_vtx[ipart]);
      for (int i = 0; i < pn_vtx[ipart] * vtx_field_stride[0]; i++) {
        double err = ABS(pvtx_field_value[0][ipart][i] - recv_val[ipart][i]);
        pvtx_error[ipart][i] = err;
        linf_error = MAX(linf_error, err);
      }
    }

    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        pvtx_ln_to_gn,
                                                        NULL,
                                                        pn_vtx,
                                                        n_part,
                                                        PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm));

    double *dvtx_error = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
                           (void **) pvtx_error,
                           NULL,
                           (void **) &dvtx_error);

    int dn_vtx = PDM_part_to_block_n_elt_block_get(ptb);
    for (int i = 0; i < dn_vtx; i++) {
      l2_error += dvtx_error[i]*dvtx_error[i];
    }
    free(dvtx_error);

    double gl2_error = 0;
    MPI_Allreduce(&l2_error, &gl2_error, 1, MPI_DOUBLE, MPI_SUM, intra_comm);

    PDM_g_num_t *distrib_vtx = PDM_part_to_block_distrib_index_get(ptb);
    gl2_error = sqrt(gl2_error/distrib_vtx[n_rank_intra]);

    PDM_part_to_block_free(ptb);

    // --> check
    if (!(gl2_error < 0.1)) {
      exit_check = 1;
    }

  }


  double glinf_error = 0;
  MPI_Allreduce(&linf_error, &glinf_error, 1, MPI_DOUBLE, MPI_MAX, comm);


  // --> check
  if (!(glinf_error < 0.2)) {
    exit_check = 1;
  }

  // free
  CWP_client_Mesh_interf_del(code_name,
                             cpl_name);
  CWP_client_Field_del(code_name,
                       cpl_name,
                       vtx_field_name[0]);
  CWP_client_Cpl_del(code_name,
                     cpl_name);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pface_vtx_idx [ipart]);
    free(pface_vtx     [ipart]);
    free(pvtx_coord    [ipart]);
    free(pface_ln_to_gn[ipart]);
    free(pvtx_ln_to_gn [ipart]);
    if (!code1) {
      free(pvtx_error[ipart]);
      free(recv_val  [ipart]);
    }
  }
  if (!code1) {
    free(pvtx_error);
    free(recv_val  );
  }
  for (int ifield = 0; ifield < n_vtx_field; ifield++) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(pvtx_field_value[ifield][ipart]);
    }
    free(vtx_field_name  [ifield]);
    free(pvtx_field_value[ifield]);
  }
  free(pn_face         );
  free(pn_vtx          );
  free(pface_vtx_idx   );
  free(pface_vtx       );
  free(pvtx_coord      );
  free(pface_ln_to_gn  );
  free(pvtx_ln_to_gn   );
  free(vtx_field_name  );
  free(vtx_field_type  );
  free(vtx_field_stride);
  free(pvtx_field_value);

  // CWP_Finalize
  CWP_client_Finalize();

  // Finalize MPI
  MPI_Comm_free(&intra_comm);
  MPI_Finalize();

  return exit_check;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
