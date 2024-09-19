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

#include "pdm_dmesh_nodal.h"
#include "pdm_block_to_part.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_part_to_block.h"

#define ABS(a)    ((a) < 0   ? -(a) : (a))
#define MAX(a, b) ((a) > (b) ?  (a) : (b))

static int polyfit_degree = 1;

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
         "  -f1             first filename.\n\n"
         "  -f2             second filename.\n\n"
         "  -v              verbose.\n\n"
         "  -n_rank1        number of MPI ranks for code1.\n\n"
         "  -n_rank2        number of MPI ranks for code2.\n\n"
         "  -tol            geometrical tolerance for mesh location.\n\n"
         "  -n_cls          number of nearest sources for interpolation.\n\n"
         "  -polyfit_degree polynomial degree for LS interpolation.\n\n"
         "  -algo           spatial interpolation algorithm.\n\n"
         "  -visu           visualize interpolation error.\n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}

/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *verbose,
  char                  *all_file_names[],
  int                    all_n_rank[],
  double                *tolerance,
  int                   *n_neighbors,
  CWP_Spatial_interp_t  *spatial_interp_algo,
  int                   *visu,
  int                   *point_interface
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
    else if (strcmp(argv[i], "-f1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_file_names[0] = argv[i];
      }
    }
    else if (strcmp(argv[i], "-f2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        all_file_names[1] = argv[i];
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
    else if (strcmp(argv[i], "-polyfit_degree") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        polyfit_degree = atoi(argv[i]);
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
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else if (strcmp(argv[i], "-point_interface") == 0) {
      *point_interface = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

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


static double
_polynomial
(
 const double x,
 const double y,
 const double z
 )
{
  double f = 1.;

  for (int d = 1; d <= polyfit_degree; d++) {
    f += (pow(x, d) + pow(y, d) + pow(z, d)) / pow(2, d);
  }

  return f;
}

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

  CWP_GCC_SUPPRESS_WARNING("-Wcast-qual")
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



  if (1) {
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
          if (1) {
            (*pvtx_field_value)[0][ipart][i] = _paper((*pvtx_coord)[ipart][3*i  ],
                                                      (*pvtx_coord)[ipart][3*i+1],
                                                      (*pvtx_coord)[ipart][3*i+2]);
          }
          else {
            (*pvtx_field_value)[0][ipart][i] = _polynomial((*pvtx_coord)[ipart][3*i  ],
                                                           (*pvtx_coord)[ipart][3*i+1],
                                                           (*pvtx_coord)[ipart][3*i+2]);
          }
        }
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

  int     verbose = 0;
  char   *all_file_names[2] = {NULL};
  int     all_n_rank    [2] = {-1, -1};
  double  tolerance         = 1e-2;
  int     n_neighbors       = 5;
  int     visu              = 0;
  int     point_interface   = 0;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  _read_args (argc,
              argv,
              &verbose,
              all_file_names,
              all_n_rank,
              &tolerance,
              &n_neighbors,
              &spatial_interp,
              &visu,
              &point_interface);

  if (all_file_names[0] == NULL) {
    all_file_names[0] = (char *) CWP_MESH_DIR"blade_0.01.vtk";
  }

  if (all_file_names[1] == NULL) {
    all_file_names[1] = (char *) CWP_MESH_DIR"blade_0.006.vtk";
  }

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
  int all_n_part[2] = {1, 1};
  int has_code[2] = {0, 0};


  has_code[0] = i_rank <  all_n_rank[0];
  has_code[1] = i_rank >= n_rank - all_n_rank[1];

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  int           *n_part            = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t   is_active_rank    = CWP_STATUS_ON;
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);
  const char   **file_name         = malloc(sizeof(char       *) * n_code);

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_part           [n_code] = all_n_part    [icode];
      file_name        [n_code] = all_file_names[icode];

      if (verbose) {
        log_trace("Running %s, coupled with %s, n_part = %d\n",
                  code_name[n_code], coupled_code_name[n_code], n_part[n_code]);
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
  const char *cpl_name = "c_new_api_wind_turbine_blade";


  CWP_Interface_t interface_dim = CWP_INTERFACE_SURFACE;
  if (point_interface) {
    assert(spatial_interp == CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES ||
           spatial_interp == CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES);
    interface_dim = CWP_INTERFACE_POINT;
  }


  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   interface_dim,
                   CWP_COMM_PAR_WITH_PART,
                   spatial_interp,
                   n_part[icode],
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

  // Mesh and field
  int           **pn_face          = malloc(sizeof(int           *) * n_code);
  int           **pn_vtx           = malloc(sizeof(int           *) * n_code);
  int          ***pface_vtx_idx    = malloc(sizeof(int          **) * n_code);
  int          ***pface_vtx        = malloc(sizeof(int          **) * n_code);
  double       ***pvtx_coord       = malloc(sizeof(double       **) * n_code);
  PDM_g_num_t  ***pface_ln_to_gn   = malloc(sizeof(PDM_g_num_t  **) * n_code);
  PDM_g_num_t  ***pvtx_ln_to_gn    = malloc(sizeof(PDM_g_num_t  **) * n_code);
  int            *n_vtx_field      = malloc(sizeof(int            ) * n_code);
  char         ***vtx_field_name   = malloc(sizeof(char         **) * n_code);
  PDM_data_t    **vtx_field_type   = malloc(sizeof(PDM_data_t    *) * n_code);
  int           **vtx_field_stride = malloc(sizeof(int           *) * n_code);
  double      ****pvtx_field_value = malloc(sizeof(double      ***) * n_code);

  int          ***point_connec     = NULL;
  if (point_interface) {
    point_connec = malloc(sizeof(int **) * n_code);
  }

  for (int icode = 0; icode < n_code; icode++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]);

    _gen_part_data(file_name[icode],
                   mesh_comm,
                   n_part[icode],
                   part_method,
                   &pn_face         [icode],
                   &pn_vtx          [icode],
                   &pface_vtx_idx   [icode],
                   &pface_vtx       [icode],
                   &pvtx_coord      [icode],
                   &pface_ln_to_gn  [icode],
                   &pvtx_ln_to_gn   [icode],
                   &n_vtx_field     [icode],
                   &vtx_field_name  [icode],
                   &vtx_field_type  [icode],
                   &vtx_field_stride[icode],
                   &pvtx_field_value[icode]);

    int block_id = -1;
    if (point_interface) {
      block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                           cpl_name,
                                           CWP_BLOCK_NODE);
      point_connec[icode] = malloc(sizeof(int *) * n_part[icode]);
    }
    else {
      block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                           cpl_name,
                                           CWP_BLOCK_FACE_POLY);
    }

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx       [icode][ipart],
                              pvtx_coord   [icode][ipart],
                              pvtx_ln_to_gn[icode][ipart]);

      if (point_interface) {
        point_connec[icode][ipart] = malloc(sizeof(int) * pn_vtx[icode][ipart]);
        for (int i = 0; i < pn_vtx[icode][ipart]; i++) {
          point_connec[icode][ipart][i] = i + 1;
        }
        CWP_Mesh_interf_block_std_set(code_name[icode],
                                      cpl_name,
                                      ipart,
                                      block_id,
                                      pn_vtx       [icode][ipart],
                                      point_connec [icode][ipart],
                                      pvtx_ln_to_gn[icode][ipart]);
      }
      else {
        CWP_Mesh_interf_f_poly_block_set(code_name[icode],
                                         cpl_name,
                                         ipart,
                                         block_id,
                                         pn_face       [icode][ipart],
                                         pface_vtx_idx [icode][ipart],
                                         pface_vtx     [icode][ipart],
                                         pface_ln_to_gn[icode][ipart]);
      }
    }

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  // Field

  if (n_vtx_field == 0) {
    PDM_error(__FILE__, __LINE__, 0, "There must be at least one vtx_field\n");
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;
  double **field_ptr = NULL;

  double ***recv_val = malloc(sizeof(double **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {

    recv_val[icode] = malloc(sizeof(double *) * n_part[icode]);
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      recv_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
    }

    CWP_Field_exch_t exch_type;
    CWP_Field_map_t  map_type;
    if (code_id[icode] == 1) {
      exch_type = CWP_FIELD_EXCH_SEND;
      map_type  = CWP_FIELD_MAP_SOURCE;
      field_ptr = pvtx_field_value[icode][0];
    }
    else {
      exch_type = CWP_FIELD_EXCH_RECV;
      map_type  = CWP_FIELD_MAP_TARGET;
      field_ptr = recv_val[icode];
    }



    CWP_Field_create(code_name[icode],
                     cpl_name,
                     vtx_field_name[icode][0],
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     vtx_field_stride[icode][0],
                     CWP_DOF_LOCATION_NODE,
                     exch_type,
                     visu_status);

    CWP_Time_step_beg(code_name[icode],
                      0.0);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Field_data_set(code_name[icode],
                         cpl_name,
                         vtx_field_name[icode][0],
                         ipart,
                         map_type,
                         field_ptr[ipart]);
    }
  }

  // Exchange
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

    sprintf(char_param, "%d", polyfit_degree);
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "polyfit_degree",
                                    CWP_INT,
                                    char_param);

    CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 1) {
      CWP_Field_issend(code_name[icode], cpl_name, vtx_field_name[icode][0]);
    }
    else {
      CWP_Field_irecv (code_name[icode], cpl_name, vtx_field_name[icode][0]);
    }

    if (code_id[icode] == 1) {
      CWP_Field_wait_issend(code_name[icode], cpl_name, vtx_field_name[icode][0]);
    }
    else {
      CWP_Field_wait_irecv (code_name[icode], cpl_name, vtx_field_name[icode][0]);
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Exchange fields OK\n");
    fflush(stdout);
  }


  // Check interpolation error
  double linf_error = 0.;
  double l2_error   = 0.;
  double **pvtx_error = NULL;
  for (int icode = 0; icode < n_code; icode++) {
    if (code_id[icode] == 2) {
      int i_rank_intra;
      int n_rank_intra;
      MPI_Comm_rank(intra_comm[icode], &i_rank_intra);
      MPI_Comm_size(intra_comm[icode], &n_rank_intra);

      pvtx_error = malloc(sizeof(double *) * n_part[icode]);
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        pvtx_error[ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
        for (int i = 0; i < pn_vtx[icode][ipart] * vtx_field_stride[icode][0]; i++) {
          double err = ABS(pvtx_field_value[icode][0][ipart][i] - recv_val[icode][ipart][i]);
          pvtx_error[ipart][i] = err;
          linf_error = MAX(linf_error, err);
        }

        if (visu) {
          char filename[999];
          sprintf(filename, "interp_error_%d.vtk", n_part[icode]*i_rank_intra + ipart);
          PDM_vtk_write_polydata_field(filename,
                                       pn_vtx[icode][ipart],
                                       pvtx_coord[icode][ipart],
                                       pvtx_ln_to_gn[icode][ipart],
                                       pn_face[icode][ipart],
                                       pface_vtx_idx[icode][ipart],
                                       pface_vtx[icode][ipart],
                                       pface_ln_to_gn[icode][ipart],
                                       NULL,
                                       NULL,
                                       "interp_error",
                                       pvtx_error[ipart]);
        }
      }

      PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          pvtx_ln_to_gn[icode],
                                                          NULL,
                                                          pn_vtx[icode],
                                                          n_part[icode],
                                                          PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]));

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
      MPI_Allreduce(&l2_error, &gl2_error, 1, MPI_DOUBLE, MPI_SUM, intra_comm[icode]);

      PDM_g_num_t *distrib_vtx = PDM_part_to_block_distrib_index_get(ptb);
      gl2_error = sqrt(gl2_error/distrib_vtx[n_rank_intra]);

      PDM_part_to_block_free(ptb);

      if (i_rank_intra == 0) {
        printf("l2_error   = %e\n", gl2_error);
        fflush(stdout);
      }
    }
  }


  double glinf_error = 0;
  MPI_Allreduce(&linf_error, &glinf_error, 1, MPI_DOUBLE, MPI_MAX, comm);


  if (i_rank == 0) {
    printf("linf_error = %e\n", glinf_error);
    fflush(stdout);
  }


  /* Free memory */
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Time_step_end(code_name[icode]);

    CWP_Mesh_interf_del(code_name[icode], cpl_name);

    CWP_Cpl_del(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      free(pface_vtx_idx [icode][ipart]);
      free(pface_vtx     [icode][ipart]);
      free(pvtx_coord    [icode][ipart]);
      free(pface_ln_to_gn[icode][ipart]);
      free(pvtx_ln_to_gn [icode][ipart]);
      free(recv_val      [icode][ipart]);
      if (code_id[icode] == 2) {
        free(pvtx_error[ipart]);
      }
      if (point_interface) {
        free(point_connec[icode][ipart]);
      }
    }
    if (code_id[icode] == 2) {
      free(pvtx_error);
    }
    if (point_interface) {
        free(point_connec[icode]);
      }
    for (int ifield = 0; ifield < n_vtx_field[icode]; ifield++) {
      for (int ipart = 0; ipart < n_part[icode]; ipart++) {
        free(pvtx_field_value[icode][ifield][ipart]);
      }
      free(vtx_field_name  [icode][ifield]);
      free(pvtx_field_value[icode][ifield]);
    }
    free(pn_face         [icode]);
    free(pn_vtx          [icode]);
    free(pface_vtx_idx   [icode]);
    free(pface_vtx       [icode]);
    free(pvtx_coord      [icode]);
    free(pface_ln_to_gn  [icode]);
    free(pvtx_ln_to_gn   [icode]);
    free(vtx_field_name  [icode]);
    free(vtx_field_type  [icode]);
    free(vtx_field_stride[icode]);
    free(pvtx_field_value[icode]);
    free(recv_val        [icode]);
  }
  free(pn_face         );
  free(pn_vtx          );
  free(pface_vtx_idx   );
  free(pface_vtx       );
  free(pvtx_coord      );
  free(pface_ln_to_gn  );
  free(pvtx_ln_to_gn   );
  free(n_vtx_field     );
  free(vtx_field_name  );
  free(vtx_field_type  );
  free(vtx_field_stride);
  free(pvtx_field_value);
  free(recv_val        );
  if (point_interface) {
    free(point_connec);
  }

  free(code_id);
  free(n_part);
  free(coupled_code_name);
  free(code_name);
  free(intra_comm);
  free(file_name);

  // Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
