/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "src/new/cwp.h"
#include "pdm_dcube_gen.h"
#include "pdm_part.h"
#include "pdm_priv.h"


static void fix_connectivity
        (
                int n_cells,
                PDM_l_num_t *cell_vtx,
                int inv11,
                int inv12,
                int inv21,
                int inv22
        ) {
  int bck1, bck2;
//  printf("Inverting idx %d with %d and %d with %d\n", inv11, inv12, inv21, inv22);
  for (int i = 0 ; i < n_cells ; ++i) {
    bck1 = cell_vtx[8 * i + inv11];
    bck2 = cell_vtx[8 * i + inv21];

    cell_vtx[8 * i + inv11] = cell_vtx[8 * i + inv12];
    cell_vtx[8 * i + inv21] = cell_vtx[8 * i + inv22];
    cell_vtx[8 * i + inv12] = bck1;
    cell_vtx[8 * i + inv22] = bck2;
  }
}

static int
binary_search
        (
                const PDM_l_num_t elem,
                const PDM_l_num_t array[],
                const PDM_l_num_t n,
                PDM_bool_t *in_array
        ) {
  int l = 0;
  int r = n;

  *in_array = PDM_FALSE;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l) / 2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] == elem) {
    *in_array = PDM_TRUE;
    return l;

  }
  else if (array[l] < elem)
    return l + 1;
  else
    return l;
}

static void compute_cell_vtx_connectivity
        (
                const PDM_l_num_t n_cell,
                const PDM_l_num_t *face_vtx_idx,
                const PDM_l_num_t *face_vtx,
                const PDM_l_num_t *cell_face_idx,
                const PDM_l_num_t *cell_face,
                PDM_l_num_t **cell_vtx_idx,
                PDM_l_num_t **cell_vtx
        ) {
  *cell_vtx_idx = malloc(sizeof(int) * (n_cell + 1));
  PDM_l_num_t *_cell_vtx_idx = *cell_vtx_idx;

  _cell_vtx_idx[0] = 0;

  size_t s_cell_vtx = 10 * n_cell;
  *cell_vtx = malloc(sizeof(PDM_l_num_t) * s_cell_vtx);

  PDM_bool_t already_in_cell;
  int pos, i;
  PDM_l_num_t icell, iface, ivtx, id_face, id_vtx;

  PDM_l_num_t n_vtx_cell;
  for (icell = 0 ; icell < n_cell ; icell++) {
    PDM_l_num_t *_cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
    n_vtx_cell = 0;

    for (iface = cell_face_idx[icell] ; iface < cell_face_idx[icell + 1] ; iface++) {
      id_face = PDM_ABS (cell_face[iface]) - 1;

      for (ivtx = face_vtx_idx[id_face] ; ivtx < face_vtx_idx[id_face + 1] ; ivtx++) {
        id_vtx = face_vtx[ivtx];

        pos = binary_search(id_vtx,
                            _cell_vtx,
                            n_vtx_cell,
                            &already_in_cell);
        if (already_in_cell == PDM_TRUE) continue;

        if (n_vtx_cell + _cell_vtx_idx[icell] >= (int) s_cell_vtx) {
          s_cell_vtx = PDM_MAX ((int) (2 * s_cell_vtx), n_vtx_cell + _cell_vtx_idx[icell]);
          *cell_vtx = realloc(*cell_vtx, sizeof(PDM_l_num_t) * s_cell_vtx);
          _cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
        }

        for (i = n_vtx_cell ; i > pos ; i--) {
          _cell_vtx[i] = _cell_vtx[i - 1];
        }
        _cell_vtx[pos] = id_vtx;
        n_vtx_cell++;
      }
    }

    _cell_vtx_idx[icell + 1] = _cell_vtx_idx[icell] + n_vtx_cell;
  }
  *cell_vtx = realloc(*cell_vtx, sizeof(PDM_l_num_t) * _cell_vtx_idx[n_cell]);
}

int main(int argc, char *argv[]) {
  // Init MPI
  MPI_Init(&argc, &argv);
  int rank;
  int comm_world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);
  assert(comm_world_size > 1);

  // Init CWIPI
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;
  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  if (rank < comm_world_size / 2) {
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }
  else {
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  CWP_Init(MPI_COMM_WORLD, n_code, (const char **) code_name, is_active_rank, intra_comm);
  printf("%d: CWIPI Init OK, %s\n", rank, code_name[0]);

  // Create CWIPI coupling
  int n_part = 1;
  const char *coupling_name = "c_new_api_dcube";
  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  CWP_Cpl_create(code_name[0], coupling_name, coupled_code_name[0], CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITH_PART, interp_method, n_part, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);

  // Setup visualisation
  CWP_Visu_set(code_name[0], coupling_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");

  // Create PDM cube mesh
  PDM_g_num_t n_vtx_seg = 4;
  double length = 1.;
  const double xmin = 0, ymin = 0., zmin = 0.;


  PDM_dcube_t *dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD, n_vtx_seg, length, xmin, ymin, zmin, PDM_OWNERSHIP_KEEP);

  // Create mesh partitions
  int d_n_cell, d_n_face, d_n_vertices, s_face_vtx, s_face_group, n_face_group;
  int *d_face_vertex_idx = NULL, *d_face_group_idx = NULL;
  PDM_g_num_t *d_face_cell = NULL, *d_face_vertex = NULL, *d_face_group = NULL;
  double *d_vertex_coord = NULL;
  PDM_dcube_gen_dim_get(dcube, &n_face_group, &d_n_cell, &d_n_face, &d_n_vertices, &s_face_vtx, &s_face_group);
  PDM_dcube_gen_data_get(dcube, &d_face_cell, &d_face_vertex_idx, &d_face_vertex, &d_vertex_coord, &d_face_group_idx, &d_face_group);

  PDM_part_t *ppart_id = 0;
  int *d_cell_part = (int *) malloc(sizeof(int) * d_n_cell);
  int have_dcell_part = 0;
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_split_t method = PDM_PART_SPLIT_HILBERT;
  ppart_id = PDM_part_create(PDM_MPI_COMM_WORLD, method, "PDM_PART_RENUM_CELL_NONE", "PDM_PART_RENUM_FACE_NONE",
                             n_property_cell, renum_properties_cell, n_property_face, renum_properties_face, n_part,
                             d_n_cell, d_n_face, d_n_vertices, n_face_group, NULL, NULL, NULL, NULL, have_dcell_part, d_cell_part,
                             d_face_cell, d_face_vertex_idx, d_face_vertex, NULL, d_vertex_coord, NULL, d_face_group_idx, d_face_group);
  // free(d_face_vertex_idx);
  // free(d_face_group_idx);
  // free(d_face_cell);
  // free(d_face_vertex);
  // free(d_face_group);
  // free(d_vertex_coord);
  free(d_cell_part);

  PDM_dcube_gen_free(dcube);

  // Get connectivity
  // Cell face connectivity
  int          *n_cells       = (int         * ) malloc(sizeof(int          ) * n_part);
  int         **cell_face_idx = (int         **) malloc(sizeof(int         *) * n_part);
  int         **cell_face     = (int         **) malloc(sizeof(int         *) * n_part);
  PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  // Face vertex connectivity
  int  *n_faces      = (int  *) malloc(sizeof(int  ) * n_part);
  int **face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
  int **face_vtx     = (int **) malloc(sizeof(int *) * n_part);

  // Cell vertex connectivity
  int **cell_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
  int **cell_vtx     = (int **) malloc(sizeof(int *) * n_part);

  // Vertices
  int          *n_vtx        = (int          *) malloc(sizeof(int          ) * n_part);
  double      **vtx_coord    = (double      **) malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0 ; i_part < n_part ; i_part++) {
    int _n_cells, _n_faces, _n_face_part_bound, _n_vtx, _n_proc, _n_total_part, _s_cell_face, _s_face_vtx, _s_face_group, _n_edge_group2;
    PDM_part_part_dim_get(ppart_id, i_part, &_n_cells, &_n_faces, &_n_face_part_bound, &_n_vtx, &_n_proc, &_n_total_part,
                          &_s_cell_face, &_s_face_vtx, &_s_face_group, &_n_edge_group2);

    int *_cell_tag, *_cell_face_idx, *_cell_face, *_face_tag, *_face_cell, *_face_vtx_idx, *_face_vtx,
            *_face_part_bound_proc_idx, *_face_part_bound_part_idx, *_face_part_bound, *_vtx_tag, *_face_group_idx, *_face_group;
    int *_cell_vtx_idx, *_cell_vtx;
    double *_vtx_coords;
    PDM_g_num_t *_cell_ln_to_gn, *_face_ln_to_gn, *_vtx_ln_to_gn, *_face_group_ln_to_gn;
    PDM_part_part_val_get(ppart_id, i_part, &_cell_tag, &_cell_face_idx, &_cell_face, &_cell_ln_to_gn,
                          &_face_tag, &_face_cell, &_face_vtx_idx, &_face_vtx, &_face_ln_to_gn,
                          &_face_part_bound_proc_idx, &_face_part_bound_part_idx, &_face_part_bound, &_vtx_tag,
                          &_vtx_coords, &_vtx_ln_to_gn, &_face_group_idx, &_face_group, &_face_group_ln_to_gn);

    compute_cell_vtx_connectivity(_n_cells, _face_vtx_idx, _face_vtx, _cell_face_idx, _cell_face, &_cell_vtx_idx, &_cell_vtx);

    // Cell face connectivity
    n_cells[i_part] = _n_cells;
    cell_face_idx[i_part] = (int *) malloc(sizeof(int) * (_n_cells + 1));
    cell_face[i_part] = (int *) malloc(sizeof(int) * _s_cell_face);
    cell_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_cells);

    memcpy(cell_face_idx[i_part], _cell_face_idx, (_n_cells + 1) * sizeof(int));
    memcpy(cell_face[i_part], _cell_face, _s_cell_face * sizeof(int));
    memcpy(cell_ln_to_gn[i_part], _cell_ln_to_gn, _n_cells * sizeof(PDM_g_num_t));

    // Face vertex connectivty
    n_faces[i_part] = _n_faces;
    face_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (_n_faces + 1));
    face_vtx[i_part] = (int *) malloc(sizeof(int) * _s_face_vtx);

    memcpy(face_vtx_idx[i_part], _face_vtx_idx, (_n_faces + 1) * sizeof(int));
    memcpy(face_vtx[i_part], _face_vtx, _s_face_vtx * sizeof(int));

    // Vertices
    n_vtx[i_part] = _n_vtx;
    vtx_coord[i_part] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
    vtx_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);

    memcpy(vtx_coord[i_part], _vtx_coords, 3 * _n_vtx * sizeof(double));
    memcpy(vtx_ln_to_gn[i_part], _vtx_ln_to_gn, _n_vtx * sizeof(PDM_g_num_t));

    // Cell vertex connectivity
    cell_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (_n_cells + 1));
    cell_vtx[i_part] = (int *) malloc(sizeof(int) * _n_cells * 8);

    memcpy(cell_vtx_idx[i_part], _cell_vtx_idx, (_n_cells + 1) * sizeof(int));
    memcpy(cell_vtx[i_part], _cell_vtx, _n_cells * 8 * sizeof(int));

    free(_cell_vtx_idx);
    free(_cell_vtx);
  }

  // Free
  PDM_part_free(ppart_id);

  fix_connectivity(n_cells[0], cell_vtx[0], 2, 3, 6, 7);

  CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, n_vtx[0], vtx_coord[0], NULL);

//  CWP_Mesh_interf_from_cellface_set(code_name[0], coupling_name, 0, n_cells[0], cell_face_idx[0], cell_face[0], n_faces[0], face_vtx_idx[0], face_vtx[0], cell_ln_to_gn[0]);
  int block_id = CWP_Mesh_interf_block_add(code_name[0], coupling_name, CWP_BLOCK_CELL_HEXA8);
  CWP_Mesh_interf_block_std_set(code_name[0], coupling_name, 0, block_id, n_cells[0], cell_vtx[0], NULL);
  CWP_Mesh_interf_finalize(code_name[0], coupling_name);

  // Set fields
  double *send_val = NULL;
  double *recv_val = NULL;
  const char *field_name = "cooX";

  if (strcmp(code_name[0], "code1") != 0) {
    send_val = (double *) malloc(sizeof(double) * n_vtx[0]);
    for (int i = 0 ; i < n_vtx[0] ; i++) send_val[i] = vtx_coord[0][3 * i];
  }
  else recv_val = (double *) malloc(sizeof(double) * n_vtx[0]);

  CWP_Status_t visu_status = CWP_STATUS_ON;
  MPI_Barrier(MPI_COMM_WORLD);

  if (strcmp(code_name[0], "code1") != 0) {
    printf("%s is sending data\n", code_name[0]);
    CWP_Field_create(code_name[0], coupling_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_INTERLEAVED, 1, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0], coupling_name, field_name, 0, CWP_FIELD_MAP_SOURCE, send_val);
  }
  else {
    CWP_Field_create(code_name[0], coupling_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_INTERLEAVED, 1, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_RECV, visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0], coupling_name, field_name, 0, CWP_FIELD_MAP_TARGET, recv_val);
  }

  // Compute weights
  CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

  if (strcmp(code_name[0], "code1") != 0) {
    printf("%s is sending data\n", code_name[0]);
    CWP_Field_issend(code_name[0], coupling_name, field_name);
    // CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
  }
  else {
    printf("%s is receiving data\n", code_name[0]);
    CWP_Field_irecv(code_name[0], coupling_name, field_name);
    // CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);

//    int n_computed_tgts = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name, 0);
//    int n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name, 0);
//    const int* computed_tgts = CWP_Computed_tgts_get(code_name[0], coupling_name, field_name, 0);
//
//    printf("n_computed_tgts = %d\n", n_computed_tgts);
//    printf("n_uncomputed_tgts = %d\n", n_uncomputed_tgts);
//
//    for (int i = 0 ; i < n_computed_tgts ; i++) {
//      printf("%12.5e %12.5e\n", recv_val[3 * i], vtx_coord[0][3 * (computed_tgts[i] - 1)]);
//    }
  }

  if (strcmp(code_name[0], "code1") != 0) {
    CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    free(send_val);
  } else {
    CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    //    int n_computed_tgts = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name, 0);
//    int n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name, 0);
//    const int* computed_tgts = CWP_Computed_tgts_get(code_name[0], coupling_name, field_name, 0);
//
//    printf("n_computed_tgts = %d\n", n_computed_tgts);
//    printf("n_uncomputed_tgts = %d\n", n_uncomputed_tgts);
//
//    for (int i = 0 ; i < n_computed_tgts ; i++) {
//      printf("%12.5e %12.5e\n", recv_val[3 * i], vtx_coord[0][3 * (computed_tgts[i] - 1)]);
//    }
    free(recv_val);
  }

  for (int i_part = 0 ; i_part < n_part ; i_part++) {
    free(cell_face_idx[i_part]);
    free(cell_face    [i_part]);
    free(face_vtx_idx [i_part]);
    free(face_vtx     [i_part]);
    free(cell_vtx_idx [i_part]);
    free(cell_vtx     [i_part]);
    free(vtx_coord    [i_part]);
    free(vtx_ln_to_gn [i_part]);
    free(cell_ln_to_gn[i_part]);
  }
  free(cell_face_idx);
  free(cell_face    );
  free(face_vtx_idx );
  free(face_vtx     );
  free(cell_vtx_idx );
  free(cell_vtx     );
  free(vtx_coord    );
  free(vtx_ln_to_gn );
  free(cell_ln_to_gn);
  free(n_cells      );
  free(n_faces      );
  free(n_vtx        );

  // Finalize
  CWP_Time_step_end(code_name[0]);
  CWP_Mesh_interf_del(code_name[0], coupling_name);
  CWP_Cpl_del(code_name[0], coupling_name);

  free(code_name);
  free(coupled_code_name);
  free(intra_comm);
  CWP_Finalize();
  MPI_Finalize();
}
