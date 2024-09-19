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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "grid_mesh.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

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
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_vtx_seg1,
  int                   *n_vtx_seg2,
  int                   *n_part1,
  int                   *n_part2,
  PDM_split_dual_t      *part_method,
  double                *tolerance,
  int                   *randomize,
  int                   *verbose
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
     else if (strcmp(argv[i], "-n_part1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
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
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_gen_mesh
(
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
 const PDM_g_num_t         n_vtx_seg,
 const int                 randomize,
 const int                 random_seed,
 const int                 is_partitionned,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 PDM_g_num_t            ***pface_ln_to_gn,
 PDM_g_num_t            ***pvtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Generate a distributed polygonal mesh */
  double xmin   = 0.;
  double ymin   = 0.;
  double zmin   = 0.;
  double length = 1.;
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1,
                                                        length,
                                                        xmin,
                                                        ymin,
                                                        zmin,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build(dcube);
  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);


  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

  if (randomize) {
    double noise = 0.2*length/(double) (n_vtx_seg - 1);
    for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
      if (ABS(dvtx_coord[3*i_vtx  ] - xmin         ) > 1.e-9 &&
          ABS(dvtx_coord[3*i_vtx  ] - xmin - length) > 1.e-9 &&
          ABS(dvtx_coord[3*i_vtx+1] - ymin         ) > 1.e-9 &&
          ABS(dvtx_coord[3*i_vtx+1] - ymin - length) > 1.e-9) {
        srand(distrib_vtx[i_rank] + i_vtx + random_seed);
        for (int i = 0; i < 2; i++) {
          dvtx_coord[3*i_vtx+i] += noise*0.5*(2*rand()/(double) RAND_MAX - 1);
        }
      }
    }
  }


  *pn_face        = (int *)          malloc(sizeof(int *)          * n_part);
  *pn_vtx         = (int *)          malloc(sizeof(int *)          * n_part);
  *pface_vtx_idx  = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_vtx      = (int **)         malloc(sizeof(int **)         * n_part);
  *pface_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t **) * n_part);
  *pvtx_coord     = (double **)      malloc(sizeof(double **)      * n_part);

  if (is_partitionned) {
    /* Spit the mesh */
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

    PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

    /* Run */
    PDM_multipart_compute(mpart);



    /* Get partitioned mesh */
    for (int i_part = 0; i_part < n_part; i_part++) {
      int          n_face;
      int          n_vtx;
      int         *face_edge;
      int         *face_edge_idx;
      PDM_g_num_t *face_ln_to_gn;
      int         *edge_vtx_idx;
      int         *edge_vtx;
      double      *vtx;
      PDM_g_num_t *vtx_ln_to_gn;

      /* Faces */
      n_face = PDM_multipart_part_connectivity_get(mpart,
                                                   0,
                                                   i_part,
                                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                   &face_edge_idx,
                                                   &face_edge,
                                                   PDM_OWNERSHIP_KEEP);
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx_idx,
                                          &edge_vtx,
                                          PDM_OWNERSHIP_KEEP);


      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      (*pn_face)[i_part] = n_face;

      (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * (n_face + 1));
      memcpy((*pface_vtx_idx)[i_part], face_edge_idx, sizeof(int) * (n_face + 1));

      PDM_compute_face_vtx_from_face_and_edge_unsigned(n_face,
                                                       face_edge_idx,
                                                       face_edge,
                                                       edge_vtx,
                                                       &(*pface_vtx)[i_part]);

      (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face);
      memcpy((*pface_ln_to_gn)[i_part], face_ln_to_gn, sizeof(PDM_g_num_t) * n_face);


      /* Vertices */
      n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                               0,
                                               i_part,
                                               &vtx,
                                               PDM_OWNERSHIP_KEEP);
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      (*pn_vtx)[i_part] = n_vtx;

      (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * n_vtx);
      memcpy((*pvtx_coord)[i_part], vtx, sizeof(double) * 3 * n_vtx);

      (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_vtx);
      memcpy((*pvtx_ln_to_gn)[i_part], vtx_ln_to_gn, sizeof(PDM_g_num_t) * n_vtx);

    }
    PDM_multipart_free(mpart);
  }

  else {
    int i_part = 0;
    // Not partitionned => allgather
    int *recv_count = malloc(sizeof(int) * n_rank);
    int *recv_shift = malloc(sizeof(int) * (n_rank + 1));
    recv_shift[0] = 0;
    PDM_MPI_Allgather(&dn_vtx, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }


    /* Vertices */
    PDM_MPI_Datatype mpi_coord;
    PDM_MPI_Type_create_contiguous(3, PDM_MPI_DOUBLE, &mpi_coord);
    PDM_MPI_Type_commit(&mpi_coord);

    (*pn_vtx)[i_part] = recv_shift[n_rank];
    (*pvtx_coord)[i_part] = (double *) malloc(sizeof(double) * 3 * (*pn_vtx)[i_part]);
    PDM_MPI_Allgatherv(dvtx_coord, dn_vtx, mpi_coord,
                       (*pvtx_coord)[i_part], recv_count, recv_shift, mpi_coord, comm);

    (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*pn_vtx)[i_part]);
    for (int i = 0; i < (*pn_vtx)[i_part]; i++) {
      (*pvtx_ln_to_gn)[i_part][i] = i + 1;
    }


    /* Faces */
    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, PDM_GEOMETRY_KIND_SURFACIC);

    int id_section = sections_id[0];
    int                   dn_face      = PDM_DMesh_nodal_section_n_elt_get  (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t          *dface_vtx    = PDM_DMesh_nodal_section_std_get    (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_Mesh_nodal_elt_t  t_elt        = PDM_DMesh_nodal_section_type_get   (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int stride = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    PDM_MPI_Allgather(&dn_face, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    (*pn_face)[i_part] = recv_shift[n_rank];

    int s_dface_vtx = stride*dn_face;
    int *_dface_vtx  = malloc(sizeof(int) * s_dface_vtx);
    int *dface_vtx_n = PDM_array_const_int(dn_face, stride);
    for (int i = 0; i < s_dface_vtx; i++) {
      _dface_vtx[i] = (int) dface_vtx[i];
    }

    int *face_vtx_n = malloc(sizeof(int) * (*pn_face)[i_part]);
    PDM_MPI_Allgatherv(dface_vtx_n, dn_face, PDM_MPI_INT,
                       face_vtx_n, recv_count, recv_shift, PDM_MPI_INT, comm);
    free(dface_vtx_n);


    (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * ((*pn_face)[i_part] + 1));
    (*pface_vtx_idx)[i_part][0] = 0;
    for (int i = 0; i < (*pn_face)[i_part]; i++) {
      (*pface_vtx_idx)[i_part][i+1] = (*pface_vtx_idx)[i_part][i] + face_vtx_n[i];
    }
    free(face_vtx_n);

    PDM_MPI_Allgather(&s_dface_vtx, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * (*pface_vtx_idx)[i_part][(*pn_face)[i_part]]);
    PDM_MPI_Allgatherv(_dface_vtx, s_dface_vtx, PDM_MPI_INT,
                       (*pface_vtx)[i_part], recv_count, recv_shift, PDM_MPI_INT, comm);
    free(recv_count);
    free(recv_shift);
    free(_dface_vtx);


    (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*pn_face)[i_part]);
    for (int i = 0; i < (*pn_face)[i_part]; i++) {
      (*pface_ln_to_gn)[i_part][i] = i + 1;
    }
  }

  PDM_dcube_nodal_gen_free(dcube);
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P0P0
 *
 *---------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  // Read args from command line
  int    n_vtx_seg1            = 4;
  int    n_vtx_seg2            = 4;
  int    n_part1               = 1;
  int    n_part2               = 1;
  int    randomize             = 1;
  double tolerance             = 1e-2;
  int    verbose               = 0;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif


  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &n_part1,
             &n_part2,
             &part_method,
             &tolerance,
             &randomize,
             &verbose);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);

  // Initialize CWIPI

  const char **code_name         = malloc(sizeof(char *) * 2);
  const char **coupled_code_name = malloc(sizeof(char *) * 2);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;


  int has_code[2] = {0, 0};
  if (i_rank < (2*n_rank) / 3) {
    has_code[0] = 1;
  }
  if (i_rank >= n_rank / 3) {
    has_code[1] = 1;
  }

  const char *all_code_names[2] = {"code1", "code2"};
  int all_n_vtx_seg[2] = {n_vtx_seg1, n_vtx_seg2};
  int all_n_part   [2] = {n_part1,    n_part2};
  CWP_Comm_t all_comm_type[2] = {CWP_COMM_PAR_WITHOUT_PART, CWP_COMM_PAR_WITH_PART};

  int n_code = 0;
  int n_vtx_seg[2];
  int n_part   [2];
  int code_id  [2];
  CWP_Comm_t comm_type[2];
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_vtx_seg        [n_code] = all_n_vtx_seg [icode];
      n_part           [n_code] = all_n_part    [icode];
      comm_type        [n_code] = all_comm_type[icode];
      if (all_comm_type[icode] == CWP_COMM_PAR_WITHOUT_PART) {
        n_part[n_code] = 1;
      }
      // log_trace("%s\n", code_name[n_code]);
      n_code++;
    }
  }

  // log_trace("n_code = %d : ", n_code);
  // for (int i_code = 0; i_code < n_code; i_code++) {
  //   log_trace("%s ", code_name[i_code]);
  // }
  // log_trace("\n");

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *cpl_name = "c_new_api_surf_P0P0_without_part";
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  for (int i_code = 0; i_code < n_code; i_code++) {
    // log_trace(">> CWP_Cpl_create %s\n", code_name[i_code]);
    CWP_Cpl_create(code_name[i_code],                                     // Code name
                   cpl_name,                                              // Coupling id
                   coupled_code_name[i_code],                             // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   comm_type[i_code],                                     // Coupling type
                   spatial_interp,
                   n_part[i_code],                                        // Number of partitions
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh displacement type
                   CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency
  }



  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Visu_set(code_name[i_code],       // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }
  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
  }

  // Mesh definition
  int          **pn_face        = malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = malloc(sizeof(double      **) * n_code);
  PDM_g_num_t ***pface_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_code);


  for (int i_code = 0; i_code < n_code; i_code++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[i_code]);
    _gen_mesh(mesh_comm,
              n_part[i_code],
              part_method,
              n_vtx_seg[i_code],
              randomize,
              code_id        [i_code],
              (comm_type[i_code] == CWP_COMM_PAR_WITH_PART),
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);

    int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);
    for (int i = 0; i < n_part[i_code]; i++) {
      CWP_Mesh_interf_vtx_set(code_name[i_code],
                              cpl_name,
                              i,
                              pn_vtx       [i_code][i],
                              pvtx_coord   [i_code][i],
                              pvtx_ln_to_gn[i_code][i]);


      CWP_Mesh_interf_f_poly_block_set(code_name[i_code],
                                       cpl_name,
                                       i,
                                       block_id,
                                       pn_face       [i_code][i],
                                       pface_vtx_idx [i_code][i],
                                       pface_vtx     [i_code][i],
                                       pface_ln_to_gn[i_code][i]);
    }

    CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
  }
  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Set mesh OK\n");
  }


  // Create and set fields
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name1 = "field1";
  const char *field_name2 = "field2";

  double ***send_val = malloc(sizeof(double **) * n_code);
  double ***recv_val = malloc(sizeof(double **) * n_code);
  for (int i_code = 0; i_code < n_code; i_code++) {
    send_val[i_code] = malloc(sizeof(double *) * n_part[i_code]);
    recv_val[i_code] = malloc(sizeof(double *) * n_part[i_code]);
    for (int i = 0; i < n_part[i_code]; i++) {
      send_val[i_code][i] = malloc(sizeof(double) * pn_face[i_code][i]);
      recv_val[i_code][i] = malloc(sizeof(double) * pn_face[i_code][i]);
    }

    if (code_id[i_code] == 1) {
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {
        for (int i = 0 ; i < pn_face[i_code][ipart]; i++) {
          send_val[i_code][ipart][i] = (double) rand() / (double) RAND_MAX;
        }
      }

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(code_name[i_code],
                        0.0);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
      }

      CWP_Involved_srcs_bcast_enable(code_name[i_code],
                                     cpl_name,
                                     field_name1);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
      }

      CWP_Computed_tgts_bcast_enable(code_name[i_code],
                                     cpl_name,
                                     field_name2);
    }

    if (code_id[i_code] == 2) {
      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Field_create(code_name[i_code],
                       cpl_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       1,
                       CWP_DOF_LOCATION_CELL_CENTER,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(code_name[i_code],
                        0.0);

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
      }

      for (int i = 0; i < n_part[i_code]; i++) {
        CWP_Field_data_set(code_name[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
      }
    }
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Create fields OK\n");
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
  }
  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Interpolation weights computation OK\n");
  }

  MPI_Barrier(comm);

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (code_id[i_code] == 1) {
      CWP_Field_issend(code_name[i_code], cpl_name, field_name1);
      CWP_Field_irecv (code_name[i_code], cpl_name, field_name2);
    }
    else {
      CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
      CWP_Field_issend(code_name[i_code], cpl_name, field_name2);
    }

    if (verbose) {
      log_trace("\n\n--- %s ---\n", code_name[i_code]);
    }


    if (code_id[i_code] == 1) {
      CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
      CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name2);
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {

        int n_involved_src = CWP_N_involved_srcs_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name1,
                                                     ipart);
        const int *involved_src = CWP_Involved_srcs_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name1,
                                                        ipart);
        if (verbose) {
          log_trace("Field 1 :\n");
          log_trace("  part %d, ", ipart);
          PDM_log_trace_array_int(involved_src, n_involved_src, "involved_src : ");
        }

        int n_computed_tgt = CWP_N_computed_tgts_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name2,
                                                     ipart);
        const int *computed_tgt = CWP_Computed_tgts_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name2,
                                                        ipart);
        if (verbose) {
          log_trace("Field 2 :\n");
          log_trace("  part %d, ", ipart);
          PDM_log_trace_array_int(computed_tgt, n_computed_tgt, "computed_tgt : ");
        }
      }
    }
    else {
      CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
      CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name2);
      for (int ipart = 0; ipart < n_part[i_code]; ipart++) {
        int n_computed_tgt = CWP_N_computed_tgts_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name1,
                                                     ipart);
        const int *computed_tgt = CWP_Computed_tgts_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name1,
                                                        ipart);
        if (verbose) {
          log_trace("Field 1 :\n");
          log_trace("  part %d, ", ipart);
          PDM_log_trace_array_int(computed_tgt, n_computed_tgt, "computed_tgt : ");
        }
        int n_involved_src = CWP_N_involved_srcs_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name2,
                                                     ipart);
        const int *involved_src = CWP_Involved_srcs_get(code_name[i_code],
                                                        cpl_name,
                                                        field_name2,
                                                        ipart);
        if (verbose) {
          log_trace("Field 2 :\n");
          log_trace("  part %d, ", ipart);
          PDM_log_trace_array_int(involved_src, n_involved_src, "involved_src : ");
        }
      }
    }
  }

  if (i_rank == 0) {
    printf("Exchange fields OK\n");
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Time_step_end(code_name[i_code]);

    CWP_Mesh_interf_del(code_name[i_code], cpl_name);

    CWP_Cpl_del(code_name[i_code], cpl_name);
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0 ; i_part < n_part[i_code]; i_part++) {
      free(pface_vtx_idx [i_code][i_part]);
      free(pface_vtx     [i_code][i_part]);
      free(pvtx_coord    [i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn [i_code][i_part]);
      free(send_val      [i_code][i_part]);
      free(recv_val      [i_code][i_part]);
    }
    free(pn_face       [i_code]);
    free(pn_vtx        [i_code]);
    free(pface_vtx_idx [i_code]);
    free(pface_vtx     [i_code]);
    free(pvtx_coord    [i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn [i_code]);
    free(send_val      [i_code]);
    free(recv_val      [i_code]);
  }
  free(pn_face);
  free(pn_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);
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

