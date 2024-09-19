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
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "grid_mesh.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_io.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_generate_mesh.h"
#include "pdm_array.h"

#include "cwp_priv.h"

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n_rank1 \n\n"
     "  -n_rank2 \n\n"
     "  -h       This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int                    argc,
 char                 **argv,
 int                    code_n_rank[]
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);
    else if (strcmp(argv[i], "-n_rank1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[0] = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_rank2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        code_n_rank[1] = atoi(argv[i]);
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}


static void
_compute_field1
(
 const int     n_vtx,
       double *vtx_coord,
       double *field_val
 )
{
  for (int i = 0; i < n_vtx; i++) {
    field_val[i] = vtx_coord[3*i];
  }
}

static void
_compute_field2
(
 const int     n_elt,
       int    *elt_vtx_idx,
       int    *elt_vtx,
       double *vtx_coord,
       double *field_val
 )
{
  for (int i = 0; i < n_elt; i++) {
    field_val[i] = 0;
    for (int j = elt_vtx_idx[i]; j < elt_vtx_idx[i+1]; j++) {
      field_val[i] += vtx_coord[3*(elt_vtx[j]-1)];
    }
    field_val[i] /= (elt_vtx_idx[i+1] - elt_vtx_idx[i]);
  }
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
)
{
  // default
  int                   code_n_rank[2]    = {-1, -1};
  _read_args(argc,
             argv,
             code_n_rank);

  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

  if (code_n_rank[0] < 0) {
    if (code_n_rank[1] < 0) {
      code_n_rank[0] = n_rank / 2;
    }
    else {
      code_n_rank[0] = n_rank - code_n_rank[1];
    }
  }
  if (code_n_rank[1] < 0) {
    code_n_rank[1] = n_rank - code_n_rank[0];
  }

  int is_code1 = (i_rank < code_n_rank[0]);


  // Initialization
  const char   *code_name         = NULL;
  const char   *coupled_code_name = NULL;
  CWP_Status_t  is_active_rank    = CWP_STATUS_ON;
  int           n_part            = 1;

  if (is_code1) {
    code_name         = "code1";
    coupled_code_name = "code2";
  }
  else {
    code_name         = "code2";
    coupled_code_name = "code1";
  }
  assert(n_part > 0);

  MPI_Comm intra_comm;

  CWP_Init(MPI_COMM_WORLD,
           1,
           &code_name,
           is_active_rank,
           &intra_comm);

  // EXIT_SUCCESS ?
  int exit_check = 0;

  // Create coupling
  const char *cpl_name = "c_new_api_multiple";
  CWP_Cpl_create(code_name,
                 cpl_name,
                 coupled_code_name,
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_VARIABLE,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name, cpl_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");


  // Create fields
  const char *field_name1 = "coord_x";
  const char *field_name2 = "elt_type";

  CWP_Status_t visu_status = CWP_STATUS_ON;

  if (is_code1) {
    CWP_Field_create(code_name,
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Field_create(code_name,
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
  }
  else {
    CWP_Field_create(code_name,
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Field_create(code_name,
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
  }

  CWP_Time_step_beg(code_name,
                    0.0);

  double *send_val1 = NULL;
  double *recv_val1 = NULL;
  double *send_val2 = NULL;
  double *recv_val2 = NULL;

  // Create mesh with several blocks
  int n_vtx = 9;
  double coords[27] = {0., 0., 0.,   1., 0., 0.,   2., 0., 0.,   0., 1., 0.,
                       1., 1., 0.,   2., 1., 0.,   0., 2., 0.,   1., 2., 0.,   2., 2., 0.};
  CWP_Mesh_interf_vtx_set(code_name,
                          cpl_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

  int first_block_id = CWP_Mesh_interf_block_add(code_name,
                                                 cpl_name,
                                                 CWP_BLOCK_FACE_TRIA3);
  int n_first_elts = 4;
  int first_connec[12] = {1, 2, 5,   1, 5, 4,   5, 6, 9,   5, 9, 8};
  CWP_Mesh_interf_block_std_set(code_name,
                                cpl_name,
                                0,
                                first_block_id,
                                n_first_elts,
                                first_connec,
                                NULL);

  // int second_block_id = CWP_Mesh_interf_block_add(code_name,
  //                                                 cpl_name,
  //                                                 CWP_BLOCK_FACE_QUAD4);

  // int n_second_elts = 2;
  // int second_connec[8] = {2, 3, 6, 5,   4, 5, 8, 7};
  // CWP_Mesh_interf_block_std_set(code_name,
  //                               cpl_name,
  //                               0,
  //                               second_block_id,
  //                               n_second_elts,
  //                               second_connec,
  //                               NULL);
  // int          *pn_vtx         = NULL;
  // int          *pn_edge        = NULL;
  // int          *pn_face        = NULL;
  // double      **pvtx_coord     = NULL;
  // int         **pedge_vtx      = NULL;
  // int         **pface_edge_idx = NULL;
  // int         **pface_edge     = NULL;
  // PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  // PDM_g_num_t **pedge_ln_to_gn = NULL;
  // PDM_g_num_t **pface_ln_to_gn = NULL;
  // PDM_generate_mesh_rectangle_ngon(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm),
  //                                  PDM_MESH_NODAL_POLY_2D,
  //                                  0.,
  //                                  0.,
  //                                  0.,
  //                                  1.,
  //                                  1.,
  //                                  10,
  //                                  10,
  //                                  n_part,
  //                                  PDM_SPLIT_DUAL_WITH_HILBERT,
  //                                  &pn_vtx,
  //                                  &pn_edge,
  //                                  &pn_face,
  //                                  &pvtx_coord,
  //                                  &pedge_vtx,
  //                                  &pface_edge_idx,
  //                                  &pface_edge,
  //                                  &pvtx_ln_to_gn,
  //                                  &pedge_ln_to_gn,
  //                                  &pface_ln_to_gn);

  // CWP_Mesh_interf_vtx_set(code_name,
  //                         cpl_name,
  //                         0,
  //                         pn_vtx[0],
  //                         pvtx_coord[0],
  //                         NULL);

  // CWP_Mesh_interf_from_faceedge_set(code_name,
  //                                   cpl_name,
  //                                   0,
  //                                   pn_face[0],
  //                                   pface_edge_idx[0],
  //                                   pface_edge[0],
  //                                   pn_edge[0],
  //                                   pedge_vtx[0],
  //                                   NULL);

  CWP_Mesh_interf_finalize(code_name, cpl_name);

  // CWP_Field_del(code_name,
  //               cpl_name,
  //               field_name1);

  // CWP_Field_del(code_name,
  //               cpl_name,
  //               field_name2);


  // if (is_code1) {
  //   CWP_Field_create(code_name,
  //                    cpl_name,
  //                    field_name1,
  //                    CWP_DOUBLE,
  //                    CWP_FIELD_STORAGE_INTERLACED,
  //                    1,
  //                    CWP_DOF_LOCATION_NODE,
  //                    CWP_FIELD_EXCH_SEND,
  //                    visu_status);

  //   CWP_Field_create(code_name,
  //                    cpl_name,
  //                    field_name2,
  //                    CWP_DOUBLE,
  //                    CWP_FIELD_STORAGE_INTERLACED,
  //                    1,
  //                    CWP_DOF_LOCATION_CELL_CENTER,
  //                    CWP_FIELD_EXCH_SEND,
  //                    visu_status);
  // }
  // else {
  //   CWP_Field_create(code_name,
  //                    cpl_name,
  //                    field_name1,
  //                    CWP_DOUBLE,
  //                    CWP_FIELD_STORAGE_INTERLACED,
  //                    1,
  //                    CWP_DOF_LOCATION_NODE,
  //                    CWP_FIELD_EXCH_RECV,
  //                    visu_status);

  //   CWP_Field_create(code_name,
  //                    cpl_name,
  //                    field_name2,
  //                    CWP_DOUBLE,
  //                    CWP_FIELD_STORAGE_INTERLACED,
  //                    1,
  //                    CWP_DOF_LOCATION_CELL_CENTER,
  //                    CWP_FIELD_EXCH_RECV,
  //                    visu_status);
  // }

  if (is_code1) {
    send_val1 = malloc(sizeof(double) * n_vtx);
    _compute_field1(n_vtx,
                    coords,
                    send_val1);

    send_val2 = malloc(sizeof(double) * n_first_elts);
    int *connec_idx = PDM_array_new_idx_from_const_stride_int(3, n_first_elts);
    _compute_field2(n_first_elts,
                    connec_idx,
                    first_connec,
                    coords,
                    send_val2);
    free(connec_idx);

    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val1);
    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val2);
  }
  else {
    recv_val1 = malloc(sizeof(double) * n_vtx);
    recv_val2 = malloc(sizeof(double) * n_first_elts);

    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val1);
    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val2);
  }

  CWP_Spatial_interp_weights_compute(code_name, cpl_name);

  if (is_code1) {
    CWP_Field_issend(code_name, cpl_name, field_name1);
    CWP_Field_issend(code_name, cpl_name, field_name2);
  }
  else {
    CWP_Field_irecv(code_name, cpl_name, field_name1);
    CWP_Field_irecv(code_name, cpl_name, field_name2);
  }

  if (is_code1) {
    CWP_Field_wait_issend(code_name, cpl_name, field_name1);
    CWP_Field_wait_issend(code_name, cpl_name, field_name2);
  }
  else {
    CWP_Field_wait_irecv(code_name, cpl_name, field_name1);
    CWP_Field_wait_irecv(code_name, cpl_name, field_name2);
  }

  // Delete mesh
  CWP_Mesh_interf_del(code_name, cpl_name);

  // free(pvtx_coord    [0]);
  // free(pedge_vtx     [0]);
  // free(pface_edge_idx[0]);
  // free(pface_edge    [0]);
  // free(pvtx_ln_to_gn [0]);
  // free(pedge_ln_to_gn[0]);
  // free(pface_ln_to_gn[0]);

  // free(pn_vtx        );
  // free(pn_edge       );
  // free(pn_face       );
  // free(pvtx_coord    );
  // free(pedge_vtx     );
  // free(pface_edge_idx);
  // free(pface_edge    );
  // free(pvtx_ln_to_gn );
  // free(pedge_ln_to_gn);
  // free(pface_ln_to_gn);



  if (is_code1) {
    free(send_val1);
    free(send_val2);
  }
  else {
    free(recv_val1);
    free(recv_val2);
  }

  CWP_Time_step_end(code_name);

  CWP_Time_step_beg(code_name,
                    1.0);

  // Create new mesh
  int n_second_vtx = 6;
  double second_coords[18] = {0., 0., 0.,   1., 0., 0.,   2., 0., 0.,   0., 1., 0.,
                              1., 1., 0.,   2., 1., 0.};
  CWP_Mesh_interf_vtx_set(code_name,
                          cpl_name,
                          0,
                          n_second_vtx,
                          second_coords,
                          NULL);

  int third_block_id = CWP_Mesh_interf_block_add(code_name,
                                                 cpl_name,
                                                 CWP_BLOCK_FACE_QUAD4);

  int n_third_elts = 2;
  int third_connec[8] = {1, 2, 5, 4,  2, 3, 6, 5};
  CWP_Mesh_interf_block_std_set(code_name,
                                cpl_name,
                                0,
                                third_block_id,
                                n_third_elts,
                                third_connec,
                                NULL);

  CWP_Mesh_interf_finalize(code_name, cpl_name);


  if (is_code1) {
    send_val1 = malloc(sizeof(double) * n_second_vtx);
    _compute_field1(n_second_vtx,
                    second_coords,
                    send_val1);

    send_val2 = malloc(sizeof(double) * n_third_elts);
    int *connec_idx = PDM_array_new_idx_from_const_stride_int(4, n_third_elts);
    _compute_field2(n_third_elts,
                    connec_idx,
                    third_connec,
                    second_coords,
                    send_val2);
    free(connec_idx);

    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val1);
    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_val2);
  }
  else {
    recv_val1 = malloc(sizeof(double) * n_second_vtx);
    recv_val2 = malloc(sizeof(double) * n_third_elts);

    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name1,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val1);
    CWP_Field_data_set(code_name,
                       cpl_name,
                       field_name2,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_val2);
  }

  CWP_Spatial_interp_weights_compute(code_name, cpl_name);

  if (is_code1) {
    CWP_Field_issend(code_name, cpl_name, field_name1);
    CWP_Field_issend(code_name, cpl_name, field_name2);
  }
  else {
    CWP_Field_irecv(code_name, cpl_name, field_name1);
    CWP_Field_irecv(code_name, cpl_name, field_name2);
  }

  if (is_code1) {
    CWP_Field_wait_issend(code_name, cpl_name, field_name1);
    CWP_Field_wait_issend(code_name, cpl_name, field_name2);
  }
  else {
    CWP_Field_wait_irecv(code_name, cpl_name, field_name1);
    CWP_Field_wait_irecv(code_name, cpl_name, field_name2);
  }



  if (is_code1) {
    free(send_val1);
    free(send_val2);
  }
  else {
    free(recv_val1);
    free(recv_val2);
  }

  CWP_Time_step_end(code_name);

  // Delete mesh
  CWP_Mesh_interf_del(code_name, cpl_name);

  // Delete coupling
  CWP_Cpl_del(code_name, cpl_name);

  // Finalize
  CWP_Finalize();

  MPI_Finalize();

  return exit_check;
}
