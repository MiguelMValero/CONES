/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MIN(a, b) ((a) < (b) ?  (a) : (b))
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
         "  -n           <> Number of vertices in band length.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
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
  PDM_g_num_t            all_n_vtx_seg[],
  int                    all_n_rank[],
  int                    all_n_part[],
  PDM_split_dual_t      *part_method,
  int                   *some_inactive_ranks,
  int                   *verbose,
  int                   *swap_codes,
  CWP_Spatial_interp_t  *spatial_interp_algo,
  PDM_Mesh_nodal_elt_t   all_elt_type[],
  int                   *rotate,
  int                   *randomize
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        all_n_vtx_seg[0] = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        all_n_vtx_seg[1] = (PDM_g_num_t) _n;
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
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-swap_codes") == 0) {
      *swap_codes = 1;
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
    else if (strcmp(argv[i], "-t1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_elt_type[0] = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_elt_type[1] = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rotate") == 0) {
      *rotate = 1;
    }
    else if (strcmp(argv[i], "-randomize") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-some_inactive_ranks") == 0) {
      *some_inactive_ranks = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static double R[3][3] =
{
  {-0.14547275709949994,  0.8415293589391187 , -0.5202557207618055 },
  { 0.9893622576902102 ,  0.12373586628506748, -0.07649678720582984},
  { 0.                 , -0.5258495730132333 , -0.8505775840931856 }
};

static void
_rotate(const int n_pts, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[j][0] * x + R[j][1] * y + R[j][2] * z;
    }
  }
}

// static void
// _unrotate(const int n_pts, double *coord) {
//   for (int i = 0 ; i < n_pts ; i++) {
//     double x = coord[3 * i];
//     double y = coord[3 * i + 1];
//     double z = coord[3 * i + 2];

//     for (int j = 0 ; j < 3 ; j++) {
//       coord[3 * i + j] = R[0][j] * x + R[1][j] * y + R[2][j] * z;
//     }
//   }
// }

static void
_gen_mesh
(
 const PDM_MPI_Comm            comm,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const PDM_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     rotate,
 const int                     randomize,
 const PDM_Mesh_nodal_elt_t    elt_type,
       int                   **pn_cell,
       int                   **pn_face,
       int                   **pn_vtx,
       int                  ***pcell_face_idx,
       int                  ***pcell_face,
       int                  ***pface_vtx_idx,
       int                  ***pface_vtx,
       double               ***pvtx_coord,
       PDM_g_num_t          ***pcell_ln_to_gn,
       PDM_g_num_t          ***pface_ln_to_gn,
       PDM_g_num_t          ***pvtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

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

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           length,
                                                           xmin,
                                                           ymin,
                                                           zmin,
                                                           elt_type,
                                                           1,
                                                           PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build (dcube);


    PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);


    if (randomize) {
      double orig[3] = {xmin, ymin, zmin};
      const PDM_g_num_t *distrib_vtx = PDM_DMesh_nodal_distrib_vtx_get(dmn);
      int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
      double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

      double step = 0.3 * length / (double) n_vtx_seg;

      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          double x = dvtx_coord[3*i+j];
          if (MIN(ABS(x - orig[j]), ABS(x - orig[j] - length)) > 1e-9) {
            dvtx_coord[3*i+j] += step * (2*(double) rand() / (double) RAND_MAX - 1);
          }
        }
      }
    }

    /*
     * Split mesh
     */
    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
    PDM_multipart_compute(mpart);

    PDM_dcube_nodal_gen_free(dcube);
  }

  else {
    // Polyhedral mesh
    PDM_g_num_t  ng_cell = 0;
    PDM_g_num_t  ng_face = 0;
    PDM_g_num_t  ng_vtx  = 0;
    int          dn_cell = 0;
    int          dn_face = 0;
    int          dn_edge = 0;
    int          dn_vtx  = 0;
    double      *dvtx_coord     = NULL;
    int         *dcell_face_idx = NULL;
    PDM_g_num_t *dcell_face     = NULL;
    int          n_face_group    = 0;
    PDM_g_num_t *dface_cell      = NULL;
    int         *dface_vtx_idx   = NULL;
    PDM_g_num_t *dface_vtx       = NULL;
    int         *dface_group_idx = NULL;
    PDM_g_num_t *dface_group     = NULL;

    PDM_poly_vol_gen (comm,
                      xmin,
                      ymin,
                      zmin,
                      length,
                      length,
                      length,
                      n_vtx_seg,
                      n_vtx_seg,
                      n_vtx_seg,
                      randomize,
                      0,
                      &ng_cell,
                      &ng_face,
                      &ng_vtx,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dcell_face_idx,
                      &dcell_face,
                      &dface_cell,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);

    /* Generate dmesh */
    PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                           dn_cell,
                                           dn_face,
                                           dn_edge,
                                           dn_vtx,
                                           comm);

    PDM_dmesh_vtx_coord_set(dmesh,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);


    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_VTX,
                               dface_vtx,
                               dface_vtx_idx,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_connectivity_set(dmesh,
                               PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               dface_cell,
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_bound_set(dmesh,
                        PDM_BOUND_TYPE_FACE,
                        n_face_group,
                        dface_group,
                        dface_group_idx,
                        PDM_OWNERSHIP_USER);

    PDM_multipart_dmesh_set (mpart, 0, dmesh);

    PDM_multipart_compute(mpart);

    PDM_dmesh_free(dmesh);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_cell);
    free(dface_vtx_idx);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);
  }


  if (rotate) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      double *_vtx_coord = NULL;
      int _n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                    0,
                                                    i_part,
                                                    &_vtx_coord,
                                                    PDM_OWNERSHIP_KEEP);

      _rotate(_n_vtx, _vtx_coord);
    }
  }


  /* Get partitioned mesh */
  *pn_cell        = malloc(sizeof(int          *) * n_part);
  *pn_face        = malloc(sizeof(int          *) * n_part);
  *pn_vtx         = malloc(sizeof(int          *) * n_part);
  *pcell_face_idx = malloc(sizeof(int         **) * n_part);
  *pcell_face     = malloc(sizeof(int         **) * n_part);
  *pface_vtx_idx  = malloc(sizeof(int         **) * n_part);
  *pface_vtx      = malloc(sizeof(int         **) * n_part);
  *pcell_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_part);
  *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = malloc(sizeof(double      **) * n_part);

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
                                        PDM_OWNERSHIP_USER);

    if (_face_vtx != NULL) {
      (*pface_vtx_idx)[ipart] = _face_vtx_idx;
      (*pface_vtx)    [ipart] = _face_vtx;
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
                                          PDM_OWNERSHIP_USER);

      int *_edge_vtx;
      int *_edge_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &_edge_vtx_idx,
                                          &_edge_vtx,
                                          PDM_OWNERSHIP_USER);

      (*pface_vtx_idx)[ipart] = _face_edge_idx;

      PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &(*pface_vtx)[ipart]);

      free(_edge_vtx);
      free(_face_edge);
    }


    /* Cells */
    (*pn_cell)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_MESH_ENTITY_CELL,
                                                        &(*pcell_ln_to_gn)[ipart],
                                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &(*pcell_face_idx)[ipart],
                                        &(*pcell_face)[ipart],
                                        PDM_OWNERSHIP_USER);
  }

  PDM_multipart_free(mpart);

}


/*----------------------------------------------------------------------
 *
 * Main : volume coupling test : c_new_api_vol
 *
 *---------------------------------------------------------------------*/
int main
(
 int   argc,
 char *argv[]
 )
{
  /* Set default values */
  PDM_g_num_t      all_n_vtx_seg[2] = {3, 3};
  int              all_n_rank[2]    = {-1, -1};
  int              all_n_part[2]    = {1, 1};
  PDM_split_dual_t part_method      = PDM_SPLIT_DUAL_WITH_HILBERT;
  int              verbose          = 0;
  int              swap_codes       = 0;
  int              rotate           = 0;
  int              randomize        = 0;
  // CWP_Spatial_interp_t spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_Spatial_interp_t spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  // CWP_Spatial_interp_t spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_IDENTITY;
  PDM_Mesh_nodal_elt_t all_elt_type[2] = {PDM_MESH_NODAL_HEXA8, PDM_MESH_NODAL_HEXA8};
  int                  some_inactive_ranks = 0;

  _read_args(argc,
             argv,
             all_n_vtx_seg,
             all_n_rank,
             all_n_part,
             &part_method,
             &some_inactive_ranks,
             &verbose,
             &swap_codes,
             &spatial_interp_algo,
             all_elt_type,
             &rotate,
             &randomize);

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  assert(n_rank > 1);

  /* Initialize CWIPI */
  const char *all_code_names[2] = {"code1", "code2"};

  for (int i = 0; i < 2; i++) {
    if (all_n_rank[i] <= 0) {
      all_n_rank[i] = n_rank;
    }

    // if (all_n_active_rank[i] <= 0 || all_n_active_rank[i] > all_n_rank[i]) {
    //   all_n_active_rank[i] = all_n_rank[i];
    // }
  }

  int has_code[2] = {0, 0};
  has_code[0] = i_rank <  all_n_rank[0];
  has_code[1] = i_rank >= n_rank - all_n_rank[1];

  if (swap_codes) {
    int tmp = has_code[0];
    has_code[0] = has_code[1];
    has_code[1] = tmp;
  }

  int n_code = has_code[0] + has_code[1];

  int                   *code_id           = malloc(sizeof(int                 ) * n_code);
  int                   *n_part            = malloc(sizeof(int                 ) * n_code);
  const char           **code_name         = malloc(sizeof(char               *) * n_code);
  const char           **coupled_code_name = malloc(sizeof(char               *) * n_code);
  MPI_Comm              *intra_comm        = malloc(sizeof(MPI_Comm            ) * n_code);
  PDM_Mesh_nodal_elt_t  *elt_type          = malloc(sizeof(PDM_Mesh_nodal_elt_t) * n_code);
  PDM_g_num_t           *n_vtx_seg         = malloc(sizeof(PDM_g_num_t         ) * n_code);
  CWP_Status_t           is_active_rank    = CWP_STATUS_ON;

  if (some_inactive_ranks && i_rank%2 != 0) {
    is_active_rank = CWP_STATUS_OFF;
  }

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_part           [n_code] = all_n_part    [icode];
      elt_type         [n_code] = all_elt_type  [icode];
      n_vtx_seg        [n_code] = all_n_vtx_seg [icode];
      // if (icode == 0) {
      //   is_active_rank = (i_rank >= all_n_rank[icode] - all_n_active_rank[icode]) ? CWP_STATUS_ON : CWP_STATUS_OFF;
      // }
      // else {
      //   is_active_rank = (i_rank < n_rank - all_n_rank[icode] + all_n_active_rank[icode]) ? CWP_STATUS_ON : CWP_STATUS_OFF;
      // }

      if (verbose) {
        log_trace("Running %s, coupled with %s, n_part = %d, active %d\n",
                  code_name[n_code], coupled_code_name[n_code], n_part[n_code], is_active_rank);
      }
      n_code++;
    }
  }

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }

  /* Create coupling */
  const char *cpl_name = "c_new_api_vol";



  for (int icode = 0; icode < n_code; icode++) {
    int i_rank_intra;
    int n_rank_intra;
    MPI_Comm_rank(intra_comm[icode], &i_rank_intra);
    MPI_Comm_size(intra_comm[icode], &n_rank_intra);

    if (is_active_rank == CWP_STATUS_ON) {
      CWP_Cpl_create(code_name[icode],
                     cpl_name,
                     coupled_code_name[icode],
                     CWP_INTERFACE_VOLUME,
                     CWP_COMM_PAR_WITH_PART,
                     spatial_interp_algo,
                     n_part[icode],
                     CWP_DYNAMIC_MESH_STATIC,
                     CWP_TIME_EXCH_USER_CONTROLLED);
    }
  }

  // !! this must be performed in 2 separate loops if the intra comms do overlap
  for (int icode = 0; icode < n_code; icode++) {
    if (is_active_rank == CWP_STATUS_ON) {
      CWP_Visu_set(code_name[icode],        // Code name
                   cpl_name,                // Coupling id
                   1,                       // Postprocessing frequency
                   CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                   "text");                 // Postprocessing option
    }
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_barrier(code_name[icode],
                    cpl_name);
  }

  if (i_rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }


  /* Define mesh */
  int          **pn_cell        = malloc(sizeof(int          *) * n_code);
  int          **pn_face        = malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = malloc(sizeof(int          *) * n_code);
  int         ***pcell_face_idx = malloc(sizeof(int         **) * n_code);
  int         ***pcell_face     = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx_idx  = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = malloc(sizeof(double      **) * n_code);
  PDM_g_num_t ***pcell_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pface_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    MPI_Comm _mesh_comm;
    MPI_Comm_split(intra_comm[icode],
                   is_active_rank,
                   i_rank,
                   &_mesh_comm);

    if (is_active_rank == CWP_STATUS_ON) {
      PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &_mesh_comm);//(void *) &intra_comm[icode]);

      srand(code_id[icode]);

      _gen_mesh(mesh_comm,
                n_part[icode],
                part_method,
                n_vtx_seg[icode],
                0.,
                0.,
                0.,
                1.,
                rotate,
                randomize,
                elt_type       [icode],
                &pn_cell       [icode],
                &pn_face       [icode],
                &pn_vtx        [icode],
                &pcell_face_idx[icode],
                &pcell_face    [icode],
                &pface_vtx_idx [icode],
                &pface_vtx     [icode],
                &pvtx_coord    [icode],
                &pcell_ln_to_gn[icode],
                &pface_ln_to_gn[icode],
                &pvtx_ln_to_gn [icode]);

      // int id_block = CWP_Mesh_interf_block_add(code_name[icode],
      //                                          cpl_name,
      //                                          CWP_BLOCK_CELL_POLY);

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        CWP_Mesh_interf_vtx_set(code_name[icode],
                                cpl_name,
                                ipart,
                                pn_vtx       [icode][ipart],
                                pvtx_coord   [icode][ipart],
                                pvtx_ln_to_gn[icode][ipart]);

        CWP_Mesh_interf_from_cellface_set(code_name[icode],
                                          cpl_name,
                                          ipart,
                                          pn_cell       [icode][ipart],
                                          pcell_face_idx[icode][ipart],
                                          pcell_face    [icode][ipart],
                                          pn_face       [icode][ipart],
                                          pface_vtx_idx [icode][ipart],
                                          pface_vtx     [icode][ipart],
                                          pcell_ln_to_gn[icode][ipart]);

        // CWP_Mesh_interf_c_poly_block_set(code_name[icode],
        //                                  cpl_name,
        //                                  ipart,id_block,
        //                                  pn_cell       [icode][ipart],
        //                                  pn_face       [icode][ipart],
        //                                  pface_vtx_idx [icode][ipart],
        //                                  pface_vtx     [icode][ipart],
        //                                  pcell_face_idx[icode][ipart],
        //                                  pcell_face    [icode][ipart],
        //                                  pcell_ln_to_gn[icode][ipart]);
      }

      CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  /* Create and set fields */
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "field";

  int stride = 1;
  double **send_val  = NULL;
  double **recv_val  = NULL;
  double **field_ptr = NULL;

  for (int icode = 0; icode < n_code; icode++) {
    if (is_active_rank == CWP_STATUS_ON) {
      // srand(code_id[icode]);

      CWP_Field_exch_t exch_type;
      CWP_Field_map_t  map_type;
      if (code_id[icode] == 1) {
        send_val = malloc(sizeof(double *) * n_part[icode]);
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {

          send_val[ipart] = malloc(sizeof(double) * pn_cell[icode][ipart] * stride);
          for (int i = 0; i < pn_cell[icode][ipart]; i++) {
            for (int j = 0; j < stride; j++) {
              send_val[ipart][stride*i+j] = (double) (j+1)*pcell_ln_to_gn[icode][ipart][i];
              // send_val[ipart][stride*i+j] = (double) rand() / (double) RAND_MAX;
              // send_val[ipart][stride*i+j] = 1;
            }
          }
        }
        exch_type = CWP_FIELD_EXCH_SEND;
        map_type  = CWP_FIELD_MAP_SOURCE;
        field_ptr = send_val;
      }
      else {
        recv_val = malloc(sizeof(double *) * n_part[icode]);
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {
          recv_val[ipart] = malloc(sizeof(double) * pn_cell[icode][ipart] * stride);
        }
        exch_type = CWP_FIELD_EXCH_RECV;
        map_type  = CWP_FIELD_MAP_TARGET;
        field_ptr = recv_val;
      }

      CWP_Field_create(code_name[icode],
                       cpl_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       stride,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       exch_type,
                       visu_status);

      CWP_Time_step_beg(code_name[icode],
                        0.0);

      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        CWP_Field_data_set(code_name[icode],
                           cpl_name,
                           field_name,
                           ipart,
                           map_type,
                           field_ptr[ipart]);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Create fields OK\n");
    fflush(stdout);
  }

  /* Map source to target */
  for (int icode = 0; icode < n_code; icode++) {
    // if (spatial_interp == CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES) {
    //   CWP_Spatial_interp_property_set(code_name[icode],
    //                                   cpl_name,
    //                                   "n_neighbors",
    //                                   CWP_INT,
    //                                   "1");
    // }
    if (is_active_rank == CWP_STATUS_ON) {
      CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Interpolation weights computation OK\n");
    fflush(stdout);
  }


  for (int icode = 0; icode < n_code; icode++) {
    if (is_active_rank == CWP_STATUS_ON) {
      if (code_id[icode] == 1) {
        CWP_Field_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_irecv (code_name[icode], cpl_name, field_name);
      }

      if (code_id[icode] == 1) {
        CWP_Field_wait_issend(code_name[icode], cpl_name, field_name);
      }
      else {
        CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Exchange fields OK\n");
    fflush(stdout);
  }


  /* Check recv field */
  int error = 0;



  if (spatial_interp_algo == CWP_SPATIAL_INTERP_FROM_INTERSECTION) {

    double l_mass  [2] = {0., 0.};
    double l_volume[2] = {0., 0.};

    for (int icode = 0; icode < n_code; icode++) {
      if (is_active_rank == CWP_STATUS_ON) {
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {

          double *field_val = NULL;
          if (code_id[icode] == 1) {
            field_val = send_val[ipart];
          }
          else {
            field_val = recv_val[ipart];
          }

          double *cell_volume = malloc(sizeof(double) * pn_cell[icode][ipart]);
          double *cell_center = malloc(sizeof(double) * pn_cell[icode][ipart] * 3);
          PDM_geom_elem_polyhedra_properties_triangulated(1,
                                                          pn_cell       [icode][ipart],
                                                          pn_face       [icode][ipart],
                                                          pface_vtx_idx [icode][ipart],
                                                          pface_vtx     [icode][ipart],
                                                          pcell_face_idx[icode][ipart],
                                                          pcell_face    [icode][ipart],
                                                          pn_vtx        [icode][ipart],
                                                          pvtx_coord    [icode][ipart],
                                                          cell_volume,
                                                          cell_center,
                                                          NULL,
                                                          NULL);

          for (int i = 0; i < pn_cell[icode][ipart]; i++) {
            if (cell_volume[i] < 0) {
              log_trace("!! code %d, cell "PDM_FMT_G_NUM" : volume = %e\n",
                        code_id[icode],
                        pcell_ln_to_gn[icode][ipart][i],
                        cell_volume[i]);
            }
            l_volume[code_id[icode]-1] += cell_volume[i];
          }

          for (int i = 0; i < pn_cell[icode][ipart]; i++) {
            l_mass[code_id[icode]-1] += field_val[i] * cell_volume[i];
          }


          free(cell_volume);
          free(cell_center);
        }
      }
    }


    double g_mass[2];
    MPI_Allreduce(l_mass, g_mass, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double g_volume[2];
    MPI_Allreduce(l_volume, g_volume, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (i_rank == 0) {
      printf("g_volume = %20.16e / %20.16e, relative diff = %e\n",
             g_volume[0], g_volume[1], ABS(g_volume[0] - g_volume[1])/ABS(g_volume[1]));

      printf("g_mass   = %20.16e / %20.16e, relative diff = %e\n",
             g_mass[0], g_mass[1], ABS(g_mass[0] - g_mass[1])/ABS(g_mass[1]));
    }
  }

  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {
  //     if (verbose) {
  //       for (int ipart = 0; ipart < n_part[icode]; ipart++) {
  //         log_trace("\n-- src part %d --\n", ipart);
  //         for (int i = 0; i < pn_face[icode][ipart]; i++) {
  //           log_trace(PDM_FMT_G_NUM" sends:\n", pface_ln_to_gn[icode][ipart][i]);
  //           for (int j = 0; j < stride; j++) {
  //             log_trace("  %f\n",
  //                       send_val[ipart][stride*i+j]);
  //           }
  //         }
  //       }
  //     }
  //   }
  //   else {
  //     for (int ipart = 0; ipart < n_part[icode]; ipart++) {
  //       if (verbose) {
  //         log_trace("\n-- tgt part %d --\n", ipart);
  //       }
  //       for (int i = 0; i < pn_face[icode][ipart]; i++) {
  //         if (verbose) {
  //           log_trace(PDM_FMT_G_NUM" received:\n", pface_ln_to_gn[icode][ipart][i]);
  //         }
  //         for (int j = 0; j < stride; j++) {
  //           double expected = (double) (j+1)*pface_ln_to_gn[icode][ipart][i];
  //           if (verbose) {
  //             log_trace("  %f (expected %f)\n",
  //                       recv_val[ipart][stride*i+j],
  //                       expected);
  //           }

  //           if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
  //             if (ABS(recv_val[ipart][stride*i+j] - expected) > 1e-12) {
  //               error = 1;
  //               printf("[%d] error for "PDM_FMT_G_NUM" : received %e, expected %e\n",
  //                      i_rank, pface_ln_to_gn[icode][ipart][i],
  //                      recv_val[ipart][stride*i+j], expected);
  //               fflush(stdout);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Check fields OK\n");
    fflush(stdout);
  }






  // /* Part data */
  // const char *part_data_name = "part_data";

  // for (int icode = 0; icode < n_code; icode++) {
  //   CWP_PartData_exch_t exch_type;
  //   if (code_id[icode] == 1) {
  //     exch_type = CWP_PARTDATA_SEND;
  //   }
  //   else {
  //     // reset recv data
  //     for (int ipart = 0; ipart < n_part[icode]; ipart++) {
  //       for (int i = 0; i < pn_face[icode][ipart] * stride; i++) {
  //         recv_val[ipart][i] = -1234;
  //       }
  //     }
  //     exch_type = CWP_PARTDATA_RECV;
  //   }

  //   CWP_Part_data_create(code_name[icode],
  //                        cpl_name,
  //                        part_data_name,
  //                        exch_type,
  //                        pface_ln_to_gn[icode],
  //                        pn_face       [icode],
  //                        n_part        [icode]);
  // }

  // MPI_Barrier(MPI_COMM_WORLD);
  // if (i_rank == 0) {
  //   printf("Create part data OK\n");
  //   fflush(stdout);
  // }

  // int request[2] = {-13};

  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {
  //     CWP_Part_data_issend(code_name[icode],
  //                          cpl_name,
  //                          part_data_name,
  //                          sizeof(double),
  //                          stride,
  //                (void **) send_val,
  //                          &request[icode]);
  //   }
  //   else {
  //     CWP_Part_data_irecv(code_name[icode],
  //                         cpl_name,
  //                         part_data_name,
  //                         sizeof(double),
  //                         stride,
  //               (void **) recv_val,
  //                         &request[icode]);
  //   }
  // }

  // PDM_log_trace_array_int(request, n_code, "request : ");


  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {
  //     CWP_Part_data_wait_issend(code_name[icode],
  //                               cpl_name,
  //                               part_data_name,
  //                               request[icode]);
  //   }
  //   else {
  //     CWP_Part_data_wait_irecv(code_name[icode],
  //                              cpl_name,
  //                              part_data_name,
  //                              request[icode]);
  //   }
  // }

  // /* Check recv part data */
  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 2) {
  //     if (verbose) {
  //       log_trace("\n--- PartData ---\n");
  //     }
  //     for (int ipart = 0; ipart < n_part[icode]; ipart++) {
  //       if (verbose) {
  //         log_trace("\n-- tgt part %d --\n", ipart);
  //       }
  //       for (int i = 0; i < pn_face[icode][ipart]; i++) {
  //         if (verbose) {
  //           log_trace(PDM_FMT_G_NUM" received:\n", pface_ln_to_gn[icode][ipart][i]);
  //         }
  //         for (int j = 0; j < stride; j++) {
  //           double expected = (double) (j+1)*pface_ln_to_gn[icode][ipart][i];
  //           if (verbose) {
  //             log_trace("  %f (expected %f)\n",
  //                       recv_val[ipart][stride*i+j],
  //                       expected);
  //           }

  //           if (ABS(recv_val[ipart][stride*i+j] - expected) > 0) {
  //             error = 1;
  //             printf("[%d] error for "PDM_FMT_G_NUM" : received %e, expected %e\n",
  //                    i_rank, pface_ln_to_gn[icode][ipart][i],
  //                    recv_val[ipart][stride*i+j], expected);
  //             fflush(stdout);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // MPI_Barrier(MPI_COMM_WORLD);
  // if (i_rank == 0) {
  //   printf("Check received part data OK\n");
  //   fflush(stdout);
  // }


  // for (int icode = 0; icode < n_code; icode++) {
  //   CWP_PartData_exch_t exch_type;
  //   if (code_id[icode] == 1) {
  //     exch_type = CWP_PARTDATA_SEND;
  //   }
  //   else {
  //     exch_type = CWP_PARTDATA_RECV;
  //   }

  //   CWP_Part_data_del(code_name[icode],
  //                     cpl_name,
  //                     part_data_name,
  //                     exch_type);
  // }


  // /* Global data */
  // const char *global_data_name = "global_data";
  // int global_stride   = 3;
  // int global_n_entity = 4;
  // int *global_data = malloc(sizeof(int) * global_stride * global_n_entity);

  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {
  //     CWP_Global_data_irecv(code_name[icode],
  //                           cpl_name,
  //                           global_data_name,
  //                           sizeof(int),
  //                           global_stride,
  //                           global_n_entity,
  //                           global_data);
  //   }
  //   else {
  //     for (int i = 0; i < global_n_entity; i++) {
  //       for (int j = 0; j < global_stride; j++) {
  //         global_data[global_stride*i + j] = (i+1) * (j+1);
  //       }
  //     }

  //     CWP_Global_data_issend(code_name[icode],
  //                            cpl_name,
  //                            global_data_name,
  //                            sizeof(int),
  //                            global_stride,
  //                            global_n_entity,
  //                            global_data);
  //   }
  // }

  // for (int icode = 0; icode < n_code; icode++) {
  //   if (code_id[icode] == 1) {
  //     if (verbose) {
  //       log_trace("\n--- GlobalData ---\n");
  //     }
  //     CWP_Global_data_wait_irecv(code_name[icode],
  //                                cpl_name,
  //                                global_data_name);
  //     for (int i = 0; i < global_n_entity; i++) {
  //       if (verbose) {
  //         log_trace("global entity %d received ", i);
  //         PDM_log_trace_array_int(global_data + global_stride*i,
  //                                 global_stride,
  //                                 "");
  //       }
  //       for (int j = 0; j < global_stride; j++) {
  //         int expected = (i+1) * (j+1);
  //         if (global_data[global_stride*i + j] != expected) {
  //           error = 1;
  //           printf("[%d] error global entity %d comp %d : received %d (expected %d)\n",
  //                  i_rank, i, j, global_data[global_stride*i + j], expected);
  //           fflush(stdout);
  //         }
  //       }
  //     }
  //   }
  //   else {
  //     CWP_Global_data_wait_issend(code_name[icode],
  //                                 cpl_name,
  //                                 global_data_name);
  //   }
  // }


  // MPI_Barrier(MPI_COMM_WORLD);
  // if (i_rank == 0) {
  //   printf("Check received global data OK\n");
  //   fflush(stdout);
  // }


  /* Free memory */
  for (int icode = 0; icode < n_code; icode++) {
    if (is_active_rank == CWP_STATUS_ON) {
      CWP_Time_step_end(code_name[icode]);

      CWP_Mesh_interf_del(code_name[icode], cpl_name);

      CWP_Cpl_del(code_name[icode], cpl_name);
    }
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (is_active_rank == CWP_STATUS_ON) {
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        free(pcell_face_idx[icode][ipart]);
        free(pcell_face    [icode][ipart]);
        free(pface_vtx_idx [icode][ipart]);
        free(pface_vtx     [icode][ipart]);
        free(pvtx_coord    [icode][ipart]);
        free(pcell_ln_to_gn[icode][ipart]);
        free(pface_ln_to_gn[icode][ipart]);
        free(pvtx_ln_to_gn [icode][ipart]);
      }
      free(pn_cell       [icode]);
      free(pn_face       [icode]);
      free(pn_vtx        [icode]);
      free(pcell_face_idx[icode]);
      free(pcell_face    [icode]);
      free(pface_vtx_idx [icode]);
      free(pface_vtx     [icode]);
      free(pvtx_coord    [icode]);
      free(pcell_ln_to_gn[icode]);
      free(pface_ln_to_gn[icode]);
      free(pvtx_ln_to_gn [icode]);

      if (code_id[icode] == 1) {
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {
          free(send_val[ipart]);
        }
        free(send_val);
      }
      else {
        for (int ipart = 0; ipart < n_part[icode]; ipart++) {
          free(recv_val[ipart]);
        }
        free(recv_val);
      }
    }
  }
  free(pn_cell       );
  free(pn_face       );
  free(pn_vtx        );
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pvtx_coord    );
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );

  free(code_id);
  free(n_part);
  free(elt_type);
  free(n_vtx_seg);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  /* Finalize CWIPI */
  CWP_Finalize();


  MPI_Barrier(MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("End\n");
    fflush(stdout);
  }

  /* Finalize MPI */
  MPI_Finalize();

  // return error;
  return error;
}
