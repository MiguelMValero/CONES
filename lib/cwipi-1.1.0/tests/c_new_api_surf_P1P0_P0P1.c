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
#include <string.h>
#include <math.h>
#include <time.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "grid_mesh.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_distrib.h"
#include "pdm_geom_elem.h"

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
_display_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n     <level>  Number of vertices in band width.\n\n"
         "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n"
         "  -h             this message.\n\n");
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
_read_args(int argc, char **argv, int *nVertex, double *randLevel) {
  int i = 1;

  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {_display_usage(EXIT_SUCCESS);}
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*nVertex = atoi(argv[i]);}
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*randLevel = atof(argv[i]);}
    }
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
  // PDM_g_num_t  ng_face         = 0;
  // PDM_g_num_t  ng_vtx          = 0;
  // PDM_g_num_t  ng_edge         = 0;
  // int          dn_vtx          = 0;
  // double      *dvtx_coord      = NULL;
  // int          dn_face         = 0;
  // int         *dface_vtx_idx   = NULL;
  // PDM_g_num_t *dface_vtx       = NULL;
  // PDM_g_num_t *dface_edge      = NULL;
  // int          dn_edge         = 0;
  // PDM_g_num_t *dedge_vtx       = NULL;
  // PDM_g_num_t *dedge_face      = NULL;
  // int          n_edge_group    = 0;
  // int         *dedge_group_idx = NULL;
  // PDM_g_num_t *dedge_group     = NULL;

  // PDM_poly_surf_gen(comm,
  //                   0.,
  //                   0.1,
  //                   0.,
  //                   0.1,
  //                   randomize,
  //                   random_seed,
  //                   n_vtx_seg,
  //                   n_vtx_seg,
  //                   &ng_face,
  //                   &ng_vtx,
  //                   &ng_edge,
  //                   &dn_vtx,
  //                   &dvtx_coord,
  //                   &dn_face,
  //                   &dface_vtx_idx,
  //                   &dface_vtx,
  //                   &dface_edge,
  //                   &dn_edge,
  //                   &dedge_vtx,
  //                   &dedge_face,
  //                   &n_edge_group,
  //                   &dedge_group_idx,
  //                   &dedge_group);
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
    // int n_join = 0;
    // PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
    //                                       dn_face,
    //                                       dn_edge,
    //                                       -1,
    //                                       dn_vtx,
    //                                       n_edge_group,
    //                                       n_join,
    //                                       comm);

    // int *djoins_ids = malloc (sizeof(int) * n_join);
    // int *dedge_join_idx = malloc (sizeof(int) * (n_join + 1));
    // dedge_join_idx[0] = 0;
    // PDM_g_num_t *dedge_join = malloc (sizeof(PDM_g_num_t) * dedge_join_idx[n_join]);

    // int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

    // PDM_dmesh_set(dmesh,
    //               dvtx_coord,
    //               dedge_vtx_idx,
    //               dedge_vtx,
    //               dedge_face,
    //               dedge_group_idx,
    //               dedge_group,
    //               djoins_ids,
    //               dedge_join_idx,
    //               dedge_join);

    // PDM_multipart_register_block(mpart, 0, dmesh);

    // /* Connection between zones */
    // int n_total_joins = 0;
    // int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
    // PDM_multipart_register_joins(mpart, n_total_joins, join_to_opposite);

    /* Run */
    PDM_multipart_compute(mpart);

    // free(djoins_ids);
    // free(dedge_join_idx);
    // free(dedge_join);
    // free(join_to_opposite);



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
                                                   &face_edge,
                                                   &face_edge_idx,
                                                   PDM_OWNERSHIP_KEEP);
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx,
                                          &edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);


      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      (*pn_face)[i_part] = n_face;

      (*pface_vtx_idx)[i_part] = face_edge_idx;

      PDM_compute_face_vtx_from_face_and_edge_unsigned(n_face,
                                                       face_edge_idx,
                                                       face_edge,
                                                       edge_vtx,
                                                       &(*pface_vtx)[i_part]);

      (*pface_ln_to_gn)[i_part] = face_ln_to_gn;

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

      (*pn_vtx)[i_part]        = n_vtx;
      (*pvtx_coord)[i_part]    = vtx;
      (*pvtx_ln_to_gn)[i_part] = vtx_ln_to_gn;

    }
    PDM_multipart_free(mpart);
    // PDM_dmesh_free(dmesh);
    // free(dedge_vtx_idx);
  }

  else {
    int i_part = 0;
    // Not partitionned => allgather
    int *recv_count = malloc(sizeof(int) * n_rank);
    int *recv_shift = malloc(sizeof(int) * (n_rank + 1));
    recv_shift[0] = 0;
    // log_trace("allgather dn_vtx...\n");
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
    // log_trace("allgatherv dvtx_coord...\n");
    PDM_MPI_Allgatherv(dvtx_coord, dn_vtx, mpi_coord,
                       (*pvtx_coord)[i_part], recv_count, recv_shift, mpi_coord, comm);

    (*pvtx_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*pn_vtx)[i_part]);
    for (int i = 0; i < (*pn_vtx)[i_part]; i++) {
      (*pvtx_ln_to_gn)[i_part][i] = i + 1;
    }


    /* Faces */

    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, PDM_GEOMETRY_KIND_SURFACIC);
    // int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, PDM_GEOMETRY_KIND_SURFACIC);

    int id_section = sections_id[0];
    // const PDM_g_num_t    *distrib_face = PDM_DMesh_nodal_distrib_section_get(dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    int                   dn_face      = PDM_DMesh_nodal_section_n_elt_get  (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t          *dface_vtx    = PDM_DMesh_nodal_section_std_get    (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_Mesh_nodal_elt_t  t_elt        = PDM_DMesh_nodal_section_type_get   (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    int stride = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    // log_trace("allgather dn_face...\n");
    PDM_MPI_Allgather(&dn_face, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    (*pn_face)[i_part] = recv_shift[n_rank];

    int s_dface_vtx = stride*dn_face;//dface_vtx_idx[dn_face];
    int *_dface_vtx  = malloc(sizeof(int) * s_dface_vtx);
    int *dface_vtx_n = PDM_array_const_int(dn_face, stride);
    for (int i = 0; i < s_dface_vtx; i++) {
      _dface_vtx[i] = (int) dface_vtx[i];
    }

    int *face_vtx_n = malloc(sizeof(int) * (*pn_face)[i_part]);
    // log_trace("allgatherv dface_vtx_n...\n");
    PDM_MPI_Allgatherv(dface_vtx_n, dn_face, PDM_MPI_INT,
                       face_vtx_n, recv_count, recv_shift, PDM_MPI_INT, comm);
    free(dface_vtx_n);


    (*pface_vtx_idx)[i_part] = (int *) malloc(sizeof(int) * ((*pn_face)[i_part] + 1));
    (*pface_vtx_idx)[i_part][0] = 0;
    for (int i = 0; i < (*pn_face)[i_part]; i++) {
      (*pface_vtx_idx)[i_part][i+1] = (*pface_vtx_idx)[i_part][i] + face_vtx_n[i];
    }
    free(face_vtx_n);

    // log_trace("allgather s_dface_vtx...\n");
    PDM_MPI_Allgather(&s_dface_vtx, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    for (int i = 0; i < n_rank; i++) {
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    (*pface_vtx)[i_part] = (int *) malloc(sizeof(int) * (*pface_vtx_idx)[i_part][(*pn_face)[i_part]]);
    // log_trace("allgatherv dface_vtx...\n");
    PDM_MPI_Allgatherv(_dface_vtx, s_dface_vtx, PDM_MPI_INT,
                       (*pface_vtx)[i_part], recv_count, recv_shift, PDM_MPI_INT, comm);
    free(_dface_vtx);


    (*pface_ln_to_gn)[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*pn_face)[i_part]);
    for (int i = 0; i < (*pn_face)[i_part]; i++) {
      (*pface_ln_to_gn)[i_part][i] = i + 1;
    }
  }



  // free(dvtx_coord);
  // free(dface_vtx_idx);
  // free(dface_vtx);
  // free(dface_edge);
  // free(dedge_vtx);
  // free(dedge_face);
  // free(dedge_group_idx);
  // free(dedge_group);
  PDM_dcube_nodal_gen_free(dcube);
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  // Read args from command line
  int nVertexSeg = 10;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  // Initialization
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;
  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  MPI_Comm *connectableLocalComm = malloc(sizeof(MPI_Comm) * n_code);
  int *connectableLocalCommSize = malloc(sizeof(int) * n_code);

  if (rank % 2 == 0) {
    printf("%d - Working for code1\n", rank);
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }

  else {
    printf("%d - Working for code2\n", rank);
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1";
  int nb_part = 1;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  CWP_Cpl_create(code_name[0],                                          // Code name
                 cpl_name,                                              // Coupling id
                 coupled_code_name[0],                                  // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,                                // Coupling type
                 spatial_interp,
                 nb_part,                                               // Number of partitions
                 CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                 CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

  printf("%d - Set visu\n", rank);
  CWP_Visu_set(code_name[0],            // Code name
               cpl_name,                // Coupling id
               1,                       // Postprocessing frequency
               CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
               "text");                 // Postprocessing option

  // Mesh definition
  printf("%d - Create mesh\n", rank);

  int nVertex;
  double **coords = NULL;
  int nElts;
  int **eltsConnecPointer = NULL;
  int **eltsConnec = NULL;

  // Domain bounds
  const double xmin = -10;
  const double xmax = 10;
  const double ymin = -10;
  const double ymax = 10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords            = (double **) malloc(sizeof(double *) * n_code);
  eltsConnecPointer = (int    **) malloc(sizeof(int    *) * n_code);
  eltsConnec        = (int    **) malloc(sizeof(int    *) * n_code);


  randLevel = 0.4;
  int randSeed = 0;

  connectableLocalComm[0] = CWP_Connectable_comm_get((char *) code_name[0]);
  MPI_Comm_size(connectableLocalComm[0], &connectableLocalCommSize[0]);

  srand(0);//time(NULL));


  if (strcmp(code_name[0], "code2") == 0) {
    randLevel = 0.2;
    randSeed  = 1;
  }

  if (0) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

#ifdef PDM_HAVE_PARMETIS
    PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
    PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
    PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif

    int          *pn_face        = NULL;
    int          *pn_vtx         = NULL;
    int         **pface_vtx_idx  = NULL;
    int         **pface_vtx      = NULL;
    double      **pvtx_coord     = NULL;
    PDM_g_num_t **pface_ln_to_gn = NULL;
    PDM_g_num_t **pvtx_ln_to_gn  = NULL;
    _gen_mesh(mesh_comm,
              nb_part,
              part_method,
              nVertexSeg,
              1,
              randSeed,
              1,
              &pn_face,
              &pn_vtx,
              &pface_vtx_idx,
              &pface_vtx,
              &pvtx_coord,
              &pface_ln_to_gn,
              &pvtx_ln_to_gn);
    nVertex              = pn_vtx       [0];
    nElts                = pn_face      [0];
    coords           [0] = pvtx_coord   [0];
    eltsConnecPointer[0] = pface_vtx_idx[0];
    eltsConnec       [0] = pface_vtx    [0];
    free(pvtx_ln_to_gn [0]);
    free(pface_ln_to_gn[0]);
    free(pvtx_ln_to_gn );
    free(pface_ln_to_gn);
    free(pn_face       );
    free(pn_vtx        );
    free(pface_vtx_idx );
    free(pface_vtx     );
    free(pvtx_coord    );
  }
  else {
    coords[0]            = (double *) malloc(sizeof(double) * 3 * nVertex);
    eltsConnecPointer[0] = (int    *) malloc(sizeof(int   ) * (nElts + 1));
    eltsConnec[0]        = (int    *) malloc(sizeof(int   ) * 4 * nElts);

    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(connectableLocalCommSize[0]),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              connectableLocalComm[0]);
  }

  printf("%d - Number of vertex   : %d\n", rank, nVertex);
  printf("%d - Number of elements : %d\n", rank, nElts);

  CWP_Mesh_interf_vtx_set(code_name[0], cpl_name, 0, nVertex, coords[0], NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0], cpl_name, CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   cpl_name,
                                   0,
                                   block_id,
                                   nElts,
                                   eltsConnecPointer[0],
                                   eltsConnec[0],
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0], cpl_name);

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  double **sendValues = (double **) malloc(sizeof(double *) * n_code);
  double **recvValues = (double **) malloc(sizeof(double *) * n_code);

  if (strcmp(code_name[0], "code1") == 0) {
    sendValues[0] = (double *) malloc(sizeof(double) * nVertex);
    recvValues[0] = (double *) malloc(sizeof(double) * nElts);
    for (int i = 0 ; i < nVertex ; i++) {
      sendValues[0][i] = coords[0][3 * i];
    }
  }
  else {
    sendValues[0] = (double *) malloc(sizeof(double) * nElts);
    recvValues[0] = (double *) malloc(sizeof(double) * nVertex);
    for (int i = 0 ; i < nElts ; i++) {
      sendValues[0][i] = (double) rand() / (double) RAND_MAX;//rank;
    }
  }

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX";
  const char *field_name2 = "rank";

  CWP_Status_t visu_status = CWP_STATUS_ON;
  printf("%d - Defining fields\n", rank);
  if (strcmp(code_name[0], "code1") == 0) {
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
  }
  
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%d - Before compute\n", rank);
  // CWP_next_recv_time_set(code_name[0],
  //                        cpl_name,
  //                        0.);
  CWP_Spatial_interp_weights_compute(code_name[0], cpl_name);

  int n_uncomputed = 0;

  if (strcmp(code_name[0], "code1") == 0) {
    n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name2, 0);
  }

  else if (strcmp(code_name[0], "code2") == 0) {
    n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name1, 0);
  }

  printf("%d - After compute %d\n", rank, n_uncomputed);

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_Field_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_irecv(code_name[0], cpl_name, field_name1);
    CWP_Field_issend(code_name[0], cpl_name, field_name2);
  }

  if (strcmp(code_name[0], "code1") == 0) {
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name1);
    CWP_Field_wait_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    CWP_Field_wait_irecv(code_name[0], cpl_name, field_name1);
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name2);
  }


  if (spatial_interp == CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
    /* Check mass conservation */
    double *field_value;
    int icode;
    if (strcmp(code_name[0], "code1") == 0) {
      field_value = recvValues[0];
      icode = 0;
    }
    else if (strcmp(code_name[0], "code2") == 0) {
      field_value = sendValues[0];
      icode = 1;
    }

    double *surface_vector = malloc(sizeof(double) * nElts * 3);
    double *center         = malloc(sizeof(double) * nElts * 3);
    PDM_geom_elem_polygon_properties(nElts,
                                     eltsConnecPointer[0],
                                     eltsConnec[0],
                                     coords[0],
                                     surface_vector,
                                     center,
                                     NULL,
                                     NULL);

    double l_integral[2] = {0., 0.};
    for (int i = 0; i < nElts; i++) {

      double area = 0.;
      for (int j = 0; j < 3; j++) {
        area += surface_vector[3*i+j] * surface_vector[3*i+j];
      }
      area = sqrt(area);

      l_integral[icode] += field_value[i] * area;
    }
    free(surface_vector);
    free(center);

    double g_integral[2];
    MPI_Allreduce(l_integral, g_integral, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    if (rank == 0) {
      printf("g_integral = %20.16e / %20.16e, error = %f%%\n",
             g_integral[0], g_integral[1], 100.*ABS(g_integral[0] - g_integral[1])/ABS(g_integral[1]));
    }
  }

  CWP_Time_step_end(code_name[0]);

  printf("%d - Delete mesh\n", rank);
  CWP_Mesh_interf_del(code_name[0], cpl_name);

  printf("%d - Delete coupling\n", rank);
  CWP_Cpl_del(code_name[0], cpl_name);

  // Freeing memory
  free(coords[0]           );
  free(eltsConnecPointer[0]);
  free(eltsConnec[0]       );

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec       );

  free(sendValues[0]);
  free(recvValues[0]);
  free(sendValues);
  free(recvValues);

  free(code_name);
  free(coupled_code_name);
  free(intra_comm);
  free(connectableLocalComm);
  free(connectableLocalCommSize);


  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}
