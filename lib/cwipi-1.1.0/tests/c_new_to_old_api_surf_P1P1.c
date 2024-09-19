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

#include "pdm_poly_surf_gen.h"
#include "pdm_part.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#define ABS(a) ((a) <  0  ? -(a) : (a))

typedef enum {
    CWP_VERSION_OLD,
    CWP_VERSION_NEW
} CWP_Version_t;


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
         "  -n           <> Number of vertices in band width.\n\n"
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
_read_args(int argc, char **argv, CWP_Version_t *version, int *n_vtx_seg1, int *n_vtx_seg2,
           double *width, double *depth, int *rotation, int *randomize, int *nProcData,
           int *part_method, CWP_Spatial_interp_t* loc_method, char** output_filename, int *verbose) {
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-new") == 0) {
      *version = CWP_VERSION_NEW;
    }
    else if (strcmp(argv[i], "-old") == 0) {
      *version = CWP_VERSION_OLD;
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
    else if (strcmp(argv[i], "-width") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *width = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-depth") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *depth = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-output") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *output_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-rot") == 0) {
      *rotation = 1;
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nProcData = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }

    i++;
  }
}


/**
 *
 * \brief  Compute faceVtx connectivity
 *
 * \param [in]      ppartId  ppart identifier
 * \param [in]      n_part    Number of partitions
 *
 * \return          faceVtx connectivity for each partition of each mesh
 */

//static void
//_compute_faceVtx(PDM_part_t *ppartId, int n_part, int **nFace, int ***faceVtxIdx, int ***faceVtx,
//                 PDM_g_num_t ***faceLNToGN, int **nVtx, double ***vtxCoord,
//                 PDM_g_num_t ***vtxLNToGN) {
//
//  *nFace = (int *) malloc(sizeof(int) * n_part);
//  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
//  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
//  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
//
//  *nVtx = (int *) malloc(sizeof(int) * n_part);
//  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
//  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
//
//  PDM_part_t *id_ppart = ppartId;
//
//  for (int ipart = 0 ; ipart < n_part ; ipart++) {
//    int _nFace;
//    int _nEdge;
//    int _nEdgePartBound;
//    int _nVtx;
//    int _nProc;
//    int _nTPart;
//    int _sFaceEdge;
//    int _sEdgeVtx;
//    int _sEdgeGroup;
//    int _nEdgeGroup2;
//
//    PDM_part_part_dim_get(id_ppart,
//                          ipart,
//                          &_nFace,
//                          &_nEdge,
//                          &_nEdgePartBound,
//                          &_nVtx,
//                          &_nProc,
//                          &_nTPart,
//                          &_sFaceEdge,
//                          &_sEdgeVtx,
//                          &_sEdgeGroup,
//                          &_nEdgeGroup2);
//
//    int *_faceTag;
//    int *_faceEdgeIdx;
//    int *_faceEdge;
//    PDM_g_num_t *_faceLNToGN;
//    int *_edgeTag;
//    int *_edgeFace;
//    int *_edgeVtxIdx;
//    int *_edgeVtx;
//    PDM_g_num_t *_edgeLNToGN;
//    int *_edgePartBoundProcIdx;
//    int *_edgePartBoundPartIdx;
//    int *_edgePartBound;
//    int *_vtxTag;
//    double *_vtx;
//    PDM_g_num_t *_vtxLNToGN;
//    int *_edgeGroupIdx;
//    int *_edgeGroup;
//    PDM_g_num_t *_edgeGroupLNToGN;
//
//    PDM_part_part_val_get(id_ppart,
//                          ipart,
//                          &_faceTag,
//                          &_faceEdgeIdx,
//                          &_faceEdge,
//                          &_faceLNToGN,
//                          &_edgeTag,
//                          &_edgeFace,
//                          &_edgeVtxIdx,
//                          &_edgeVtx,
//                          &_edgeLNToGN,
//                          &_edgePartBoundProcIdx,
//                          &_edgePartBoundPartIdx,
//                          &_edgePartBound,
//                          &_vtxTag,
//                          &_vtx,
//                          &_vtxLNToGN,
//                          &_edgeGroupIdx,
//                          &_edgeGroup,
//                          &_edgeGroupLNToGN);
//
//    (*nFace)[ipart] = _nFace;
//    (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
//    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
//    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);
//
//    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
//    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));
//
//    (*nVtx)[ipart] = _nVtx;
//    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
//    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
//
//    memcpy ((*vtxCoord)[ipart], _vtx, 3 * _nVtx * sizeof(double));
//    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));
//
//    int *_faceVtx = (*faceVtx)[ipart];
//    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));
//
//    for (int i = 0 ; i < _nVtx + 1 ; i++) {
//      vtxEdgeIdx[i] = 0;
//    }
//
//    for (int i = 0 ; i < _nEdge ; i++) {
//      int ivtx1 = _edgeVtx[2 * i];
//      int ivtx2 = _edgeVtx[2 * i + 1];
//
//      vtxEdgeIdx[ivtx1] += 1;
//      vtxEdgeIdx[ivtx2] += 1;
//    }
//
//    for (int i = 1 ; i < _nVtx + 1 ; i++) {
//      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i - 1];
//    }
//
//    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
//    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
//    for (int i = 0 ; i < _nVtx ; i++) {
//      vtxEdgeN[i] = 0;
//    }
//
//    for (int i = 0 ; i < _nEdge ; i++) {
//      int ivtx1 = _edgeVtx[2 * i] - 1;
//      int ivtx2 = _edgeVtx[2 * i + 1] - 1;
//      int iedge = i + 1;
//
//      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
//      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
//      vtxEdgeN[ivtx1] += 1;
//      vtxEdgeN[ivtx2] += 1;
//    }
//    free(vtxEdgeN);
//
//    for (int i = 0 ; i < _nFace ; i++) {
//      int idx = _faceEdgeIdx[i];
//      int _nEdge2 = _faceEdgeIdx[i + 1] - idx;
//      int *_edges = _faceEdge + idx;
//      int *_vertices = _faceVtx + idx;
//
//      int edge_cur = _edges[0];
//      int vtx_deb = _edgeVtx[2 * (edge_cur - 1)];
//      _vertices[0] = vtx_deb;
//      int vtx_cur = _edgeVtx[2 * (edge_cur - 1) + 1];
//      int idxVtx = 0;
//
//      while (vtx_deb != vtx_cur) {
//        _vertices[++idxVtx] = vtx_cur;
//        int find_vtx = 0;
//
//        for (int j = vtxEdgeIdx[vtx_cur - 1] ; j < vtxEdgeIdx[vtx_cur] ; j++) {
//          for (int k = 0 ; k < _nEdge2 ; k++) {
//            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
//              edge_cur = _edges[k];
//              if (_edgeVtx[2 * (_edges[k] - 1)] == vtx_cur) {
//                vtx_cur = _edgeVtx[2 * (_edges[k] - 1) + 1];
//              }
//              else {
//                vtx_cur = _edgeVtx[2 * (_edges[k] - 1)];
//              }
//              find_vtx = 1;
//              break;
//            }
//          }
//          if (find_vtx) {
//            break;
//          }
//        }
//        if (!find_vtx) {
//          PDM_error(__FILE__, __LINE__, 0, "Error to compute vtxedge !!!!\n");
//          abort();
//        }
//      }
//    }
//
//    free(vtxEdge);
//    free(vtxEdgeIdx);
//  }
//}

static void
_get_connectivity(PDM_part_t *ppartId, int n_part, int **nFace, int ***faceEdgeIdx, int ***faceEdge,
                  int ***faceVtxIdx, int ***faceVtx, PDM_g_num_t ***faceLNToGN, int **nEdge,
                  int ***edgeVtxIdx, int ***edgeVtx, int **nVtx, double ***vtxCoord,
                  PDM_g_num_t ***vtxLNToGN) {
  *nFace = (int *) malloc(sizeof(int) * n_part);
  *faceEdgeIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceEdge = (int **) malloc(sizeof(int *) * n_part);
  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *nEdge = (int *) malloc(sizeof(int) * n_part);
  *edgeVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *edgeVtx = (int **) malloc(sizeof(int *) * n_part);

  *nVtx = (int *) malloc(sizeof(int) * n_part);
  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  PDM_part_t *id_ppart = ppartId;

  for (int ipart = 0 ; ipart < n_part ; ipart++) {
    int _nFace;
    int _nEdge;
    int _nEdgePartBound;
    int _nVtx;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get(id_ppart,
                          ipart,
                          &_nFace,
                          &_nEdge,
                          &_nEdgePartBound,
                          &_nVtx,
                          &_nProc,
                          &_nTPart,
                          &_sFaceEdge,
                          &_sEdgeVtx,
                          &_sEdgeGroup,
                          &_nEdgeGroup2);

    int *_faceTag;
    int *_faceEdgeIdx;
    int *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int *_edgeTag;
    int *_edgeFace;
    int *_edgeVtxIdx;
    int *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int *_edgePartBoundProcIdx;
    int *_edgePartBoundPartIdx;
    int *_edgePartBound;
    int *_vtxTag;
    double *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int *_edgeGroupIdx;
    int *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get(id_ppart,
                          ipart,
                          &_faceTag,
                          &_faceEdgeIdx,
                          &_faceEdge,
                          &_faceLNToGN,
                          &_edgeTag,
                          &_edgeFace,
                          &_edgeVtxIdx,
                          &_edgeVtx,
                          &_edgeLNToGN,
                          &_edgePartBoundProcIdx,
                          &_edgePartBoundPartIdx,
                          &_edgePartBound,
                          &_vtxTag,
                          &_vtx,
                          &_vtxLNToGN,
                          &_edgeGroupIdx,
                          &_edgeGroup,
                          &_edgeGroupLNToGN);

    // Faces
    (*nFace)[ipart] = _nFace;
    (*faceEdgeIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceEdge)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

    memcpy ((*faceEdgeIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceEdge)[ipart], _faceEdge, _sFaceEdge * sizeof(int));
    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    // Edges
    (*nEdge)[ipart] = _nEdge;
    (*edgeVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
    (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

    memcpy ((*edgeVtxIdx)[ipart], _edgeVtxIdx, (_nEdge + 1) * sizeof(int));
    memcpy ((*edgeVtx)[ipart], _edgeVtx, _sEdgeVtx * sizeof(int));

    // Vertices
    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 * _nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));

    // Compute face-vtx connectivity
    int *_faceVtx = (*faceVtx)[ipart];
    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));

    for (int i = 0 ; i < _nVtx + 1 ; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0 ; i < _nEdge ; i++) {
      int ivtx1 = _edgeVtx[2 * i];
      int ivtx2 = _edgeVtx[2 * i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1 ; i < _nVtx + 1 ; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i - 1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
    for (int i = 0 ; i < _nVtx ; i++) {
      vtxEdgeN[i] = 0;
    }

    for (int i = 0 ; i < _nEdge ; i++) {
      int ivtx1 = _edgeVtx[2 * i] - 1;
      int ivtx2 = _edgeVtx[2 * i + 1] - 1;
      int iedge = i + 1;

      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
      vtxEdgeN[ivtx1] += 1;
      vtxEdgeN[ivtx2] += 1;
    }
    free(vtxEdgeN);

    for (int i = 0 ; i < _nFace ; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i + 1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _faceVtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb = _edgeVtx[2 * (edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur = _edgeVtx[2 * (edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtxEdgeIdx[vtx_cur - 1] ; j < vtxEdgeIdx[vtx_cur] ; j++) {
          for (int k = 0 ; k < __nEdge ; k++) {
            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edgeVtx[2 * (_edges[k] - 1)] == vtx_cur) {
                vtx_cur = _edgeVtx[2 * (_edges[k] - 1) + 1];
              }
              else {
                vtx_cur = _edgeVtx[2 * (_edges[k] - 1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx) {
            break;
          }
        }
        if (!find_vtx) {
          PDM_error(__FILE__, __LINE__, 0, "Error to compute vtxedge !!!!\n");
          abort();
        }
      }
    }

    free(vtxEdge);
    free(vtxEdgeIdx);
  }
}


static void
_add_depth(const int n_pts, const double width, const double depth, double *coord) {
  for (int i = 0 ; i < n_pts ; i++) {
    double x = 2. * coord[3 * i] / width;
    double y = 2. * coord[3 * i + 1] / width;
    coord[3 * i + 2] = 0.5 * depth * (1. - (x * x + y * y));
  }
}


static void
_rotate(const int n_pts, double *coord) {
  double R[3][3] = {{0.9362934,  -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025,  -0.0978434},
                    {-0.1593451, 0.1537920,  0.9751703}};

  for (int i = 0 ; i < n_pts ; i++) {
    double x = coord[3 * i];
    double y = coord[3 * i + 1];
    double z = coord[3 * i + 2];

    for (int j = 0 ; j < 3 ; j++) {
      coord[3 * i + j] = R[j][0] * x + R[j][1] * y + R[j][2] * z;
    }
  }
}


static void
_create_split_mesh(int activeRank, PDM_MPI_Comm pdm_mpi_comm, double xmin, double ymin,
                   PDM_g_num_t nVtxSeg, double length, double depth, int rotation, int n_part,
                   PDM_part_split_t method, int haveRandom, int initRandom, PDM_g_num_t *nGFace,
                   PDM_g_num_t *nGVtx, int **nFace, int ***faceEdgeIdx, int ***faceEdge,
                   int ***faceVtxIdx, int ***faceVtx, PDM_g_num_t ***faceLNToGN, int **nEdge,
                   int ***edgeVtxIdx, int ***edgeVtx, int **nVtx, double ***vtxCoord,
                   PDM_g_num_t ***vtxLNToGN) {
  int i_rank;
  int numProcs;

  if (activeRank) {
    PDM_MPI_Comm_rank(pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size(pdm_mpi_comm, &numProcs);

    double xmax = xmin + length;
    double ymax = ymin + length;
    PDM_g_num_t nx = nVtxSeg;
    PDM_g_num_t ny = nVtxSeg;

    int dNFace;
    int dNVtx;
    int dNEdge;
    int *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    double *dVtxCoord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int nEdgeGroup;
    int *dEdgeGroupIdx;
    PDM_g_num_t *dEdgeGroup;

    // Create mesh
    PDM_g_num_t nGEdge;

    PDM_poly_surf_gen(pdm_mpi_comm,
                      xmin,
                      xmax,
                      ymin,
                      ymax,
                      haveRandom,
                      initRandom,
                      nx,
                      ny,
                      nGFace,
                      nGVtx,
                      &nGEdge,
                      &dNVtx,
                      &dVtxCoord,
                      &dNFace,
                      &dFaceVtxIdx,
                      &dFaceVtx,
                      &dFaceEdge,
                      &dNEdge,
                      &dEdgeVtx,
                      &dEdgeFace,
                      &nEdgeGroup,
                      &dEdgeGroupIdx,
                      &dEdgeGroup);

    _add_depth(dNVtx, length, depth, dVtxCoord);

    if (rotation) {
      _rotate(dNVtx, dVtxCoord);
    }
    // Create mesh partitions
    int have_dCellPart = 0;

    int *dCellPart = (int *) malloc(dNFace * sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc((dNEdge + 1) * sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0 ; i < dNEdge ; i++) {
      dEdgeVtxIdx[i + 1] = 2 + dEdgeVtxIdx[i];
    }

    // Split mesh
    PDM_part_t *ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    ppartId = PDM_part_create(pdm_mpi_comm,
                              method,
                              "PDM_PART_RENUM_CELL_NONE",
                              "PDM_PART_RENUM_FACE_NONE",
                              nPropertyCell,
                              renum_properties_cell,
                              nPropertyFace,
                              renum_properties_face,
                              n_part,
                              dNFace,
                              dNEdge,
                              dNVtx,
                              nEdgeGroup,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              have_dCellPart,
                              dCellPart,
                              dEdgeFace,
                              dEdgeVtxIdx,
                              dEdgeVtx,
                              NULL,
                              dVtxCoord,
                              NULL,
                              dEdgeGroupIdx,
                              dEdgeGroup);

    free(dCellPart);

    free(dVtxCoord);
    free(dFaceVtxIdx);
    free(dFaceVtx);
    free(dFaceEdge);
    free(dEdgeVtxIdx);
    free(dEdgeVtx);
    free(dEdgeFace);
    free(dEdgeGroupIdx);
    free(dEdgeGroup);

    _get_connectivity(ppartId,
                      n_part,
                      nFace,
                      faceEdgeIdx,
                      faceEdge,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nEdge,
                      edgeVtxIdx,
                      edgeVtx,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN);

    PDM_part_free(ppartId);
  }
  else {
    *nFace = (int *) malloc(sizeof(int) * n_part);
    *faceEdgeIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceEdge = (int **) malloc(sizeof(int *) * n_part);
    *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceVtx = (int **) malloc(sizeof(int *) * n_part);
    *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    *nEdge = (int *) malloc(sizeof(int) * n_part);
    *edgeVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *edgeVtx = (int **) malloc(sizeof(int *) * n_part);

    *nVtx = (int *) malloc(sizeof(int) * n_part);
    *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
    *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0 ; ipart < n_part ; ipart++) {
      int _nFace = 0;
      int _sFaceEdge = 0;
      (*nFace)[ipart] = _nFace;
      (*faceEdgeIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceEdgeIdx)[ipart][0] = 0;
      (*faceEdge)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceVtxIdx)[ipart][0] = 0;
      (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

      int _nEdge = 0;
      int _sEdgeVtx = 0;
      (*nEdge)[ipart] = _nEdge;
      (*edgeVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
      (*edgeVtxIdx)[ipart][0] = 0;
      (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

      int _nVtx = 0;
      (*nVtx)[ipart] = _nVtx;
      (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
      (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
    }
  }

  PDM_MPI_Bcast(nGFace, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast(nGVtx, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
}


static int
_set_rank_has_mesh(const MPI_Comm comm, const int nProcData, PDM_MPI_Comm *meshComm) {
  int current_rank_has_mesh = 1;
  int rank;
  int commSize;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commSize);

  PDM_MPI_Comm _comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &comm);

  if (nProcData > 0 && nProcData < commSize) {
    int rankInNode = PDM_io_mpi_node_rank(_comm);

    int nNode = 0;
    int iNode = -1;
    int masterRank = (rankInNode == 0);

    int *rankInNodes = malloc(sizeof(int) * commSize);

    MPI_Allreduce(&masterRank, &nNode, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allgather(&rankInNode, 1, MPI_INT, rankInNodes, 1, MPI_INT, comm);

    current_rank_has_mesh = 0;

    for (int i = 0 ; i < rank ; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && masterRank) {
        current_rank_has_mesh = 1;
      }
    }
    else {
      if (rankInNode < (nProcData / nNode)) {
        current_rank_has_mesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        current_rank_has_mesh = 1;
      }
    }

    PDM_MPI_Comm_split(_comm, current_rank_has_mesh, rank, meshComm);
    free(rankInNodes);
  }

  return current_rank_has_mesh;
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  // Read args from command line
  CWP_Version_t version = CWP_VERSION_NEW;
  int n_vtx_seg1 = 4;
  int n_vtx_seg2 = 4;
  int randomize = 1;
  int n_proc_data = -1;
  CWP_Spatial_interp_t loc_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t part_method = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t part_method = PDM_PART_SPLIT_PTSCOTCH;
#else
  PDM_part_split_t part_method = PDM_PART_SPLIT_HILBERT;
#endif
#endif
  int verbose = 0;

  double width = 20.;
  double depth = 1.;
  int rotation = 0;

  char* output_filename = NULL;
  int filedump = 0;

  _read_args(argc,
             argv,
             &version,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &width,
             &depth,
             &rotation,
             &randomize,
             &n_proc_data,
             (int *) &part_method,
             &loc_method,
             &output_filename,
             &verbose);

  if (output_filename !=NULL) {
    filedump = 1;    
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }

  // Initialize CWIPI
  int n_part = 1;
  int n_code = 1;
  int code_id;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (rank < comm_world_size / 2) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  if (version == CWP_VERSION_OLD) {
    cwipi_init(MPI_COMM_WORLD, code_name[0], intra_comm);
  }

  else {
    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);
  }

  if (verbose && rank == 0) {
    printf("CWIPI Init OK\n");
  }

  // Create coupling
  const char *coupling_name = "c_surf_cpl_P1P1";

  if (version == CWP_VERSION_OLD) {
    cwipi_create_coupling(coupling_name,                             // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          coupled_code_name[0],                      // Coupled application id
                          2,                                         // Geometric entities dimension
                          0.01,                                       // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          CWIPI_SOLVER_CELL_VERTEX,                  // Solver type
                          -1,                                        // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");
  }
  else {
    CWP_Cpl_create(code_name[0],
                   coupling_name,
                   coupled_code_name[0],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   loc_method,
                   n_part,
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

    // CWP_Visu_set (code_name[0],
    //               coupling_name,
    //               -1,
    //               CWP_VISU_FORMAT_ENSIGHT,
    //               "text");

  }

  if (verbose && rank == 0) {
    printf("Create coupling OK\n");
  }

  // Define mesh
  PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

  int _n_proc_data = n_proc_data;
  if (n_proc_data > 0) {
    if (code_id == 1) {
      _n_proc_data /= 2;
    }
    else {
      _n_proc_data -= n_proc_data / 2;
    }
  }
  int current_rank_has_mesh = _set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  int true_n_proc_data;
  MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("nb procs with mesh data = %d\n", true_n_proc_data);
  }

  const double xmin = -0.5 * width;
  const double ymin = -0.5 * width;
//  int init_random = (int) time(NULL);
  int init_random = 5;

  PDM_g_num_t nGFace;
  PDM_g_num_t nGVtx;
  int *nFace = NULL;
  PDM_g_num_t **faceLNToGN = NULL;
  int **faceEdgeIdx = NULL;
  int **faceEdge = NULL;
  int **faceVtxIdx = NULL;
  int **faceVtx = NULL;
  int *nEdge = NULL;
  int **edgeVtxIdx = NULL;
  int **edgeVtx = NULL;
  int *nVtx = NULL;
  double **vtxCoord = NULL;
  PDM_g_num_t **vtxLNToGN = NULL;

  PDM_MPI_Comm code_mesh_comm;
  PDM_MPI_Comm_split(mesh_comm, code_id, rank, &code_mesh_comm);

  if (code_id == 2) {
    init_random++;
  }
  _create_split_mesh(current_rank_has_mesh,
                     code_mesh_comm,
                     xmin,
                     ymin,
                     n_vtx_seg,
                     width,
                     depth,
                     rotation,
                     n_part,
                     part_method,
                     randomize,
                     init_random,
                     &nGFace,
                     &nGVtx,
                     &nFace,
                     &faceEdgeIdx,
                     &faceEdge,
                     &faceVtxIdx,
                     &faceVtx,
                     &faceLNToGN,
                     &nEdge,
                     &edgeVtxIdx,
                     &edgeVtx,
                     &nVtx,
                     &vtxCoord,
                     &vtxLNToGN);

  // Set interface mesh
  if (version == CWP_VERSION_OLD) {
    cwipi_define_mesh(coupling_name, nVtx[0], nFace[0], vtxCoord[0], faceVtxIdx[0], faceVtx[0]);
  }
  else {
    CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, nVtx[0], vtxCoord[0], vtxLNToGN[0]);

    CWP_Mesh_interf_from_faceedge_set(code_name[0],
                                      coupling_name,
                                      0,
                                      nFace[0],
                                      faceEdgeIdx[0],
                                      faceEdge[0],
                                      nEdge[0],
                                      edgeVtx[0],
                                      faceLNToGN[0]);

    CWP_Mesh_interf_finalize(code_name[0], coupling_name);
  }

  if (verbose && rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  double *send_val = NULL;
  double *recv_val = NULL;

  const char *field_name = "cooX";
  const char *field_name2 = "coocooY";

  if (code_id == 1) {
    send_val = (double *) malloc(sizeof(double) * nVtx[0]);

    for (int i = 0 ; i < nVtx[0] ; i++) {
      send_val[i] = vtxCoord[0][3 * i];
    }
  }
  else {
    recv_val = (double *) malloc(sizeof(double) * nVtx[0]);
  }

  if (version == CWP_VERSION_NEW) {
    CWP_Status_t visu_status = CWP_STATUS_OFF;
    MPI_Barrier(MPI_COMM_WORLD);

    if (code_id == 1) {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_SOURCE,
                         send_val);
    }
    else {
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name2,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         0,
                         CWP_FIELD_MAP_TARGET,
                         recv_val);
    }
  }

  if (verbose && rank == 0) {
    printf("Fields OK\n");
  }

  // Perform geometric algorithm
  PDM_timer_t *timer = PDM_timer_create();
  double t_start, t_end;

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_init(timer);

  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  // int n_unlocated = 0;
  int n_located = 0;
  const int *located = NULL;
  if (version == CWP_VERSION_OLD) {
    cwipi_locate(coupling_name);

    if (code_id != 1) {
      // n_unlocated = cwipi_get_n_not_located_points(coupling_name);
      n_located = cwipi_get_n_located_points(coupling_name);
      located = cwipi_get_located_points(coupling_name);
    }

  }
  else {
    PDM_part_to_block_global_statistic_reset();
    PDM_block_to_part_global_statistic_reset();

    CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", CWP_DOUBLE, "1e-2");
    CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);
  
    if (code_id != 1) {
      // n_unlocated = CWP_N_uncomputed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      n_located = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name2, 0);
      located =  CWP_Computed_tgts_get(code_name[0], coupling_name, field_name2, 0);
    }

    double min_elaps_create_ptb;
    double max_elaps_create_ptb;
    double min_cpu_create_ptb;
    double max_cpu_create_ptb;
    double min_elaps_create2_ptb;
    double max_elaps_create2_ptb;
    double min_cpu_create2_ptb;
    double max_cpu_create2_ptb;
    double min_elaps_exch_ptb;
    double max_elaps_exch_ptb;
    double min_cpu_exch_ptb;
    double max_cpu_exch_ptb;

    PDM_part_to_block_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_ptb,
                                       &max_elaps_create_ptb,
                                       &min_cpu_create_ptb,
                                       &max_cpu_create_ptb,
                                       &min_elaps_create2_ptb,
                                       &max_elaps_create2_ptb,
                                       &min_cpu_create2_ptb,
                                       &max_cpu_create2_ptb,
                                       &min_elaps_exch_ptb,
                                       &max_elaps_exch_ptb,
                                       &min_cpu_exch_ptb,
                                       &max_cpu_exch_ptb);

    double min_elaps_create_btp;
    double max_elaps_create_btp;
    double min_cpu_create_btp;
    double max_cpu_create_btp;
    double min_elaps_exch_btp;
    double max_elaps_exch_btp;
    double min_cpu_exch_btp;
    double max_cpu_exch_btp;

    PDM_block_to_part_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_btp,
                                       &max_elaps_create_btp,
                                       &min_cpu_create_btp,
                                       &max_cpu_create_btp,
                                       &min_elaps_exch_btp,
                                       &max_elaps_exch_btp,
                                       &min_cpu_exch_btp,
                                       &max_cpu_exch_btp);

    int min_exch_rank_send_ptb;
    int min_exch_rank_recv_ptb;
    int max_exch_rank_send_ptb;
    int max_exch_rank_recv_ptb;
    unsigned long long min_exch_data_send_ptb;
    unsigned long long min_exch_data_recv_ptb;
    unsigned long long max_exch_data_send_ptb;
    unsigned long long max_exch_data_recv_ptb;

    PDM_part_to_block_global_statistic_get(PDM_MPI_COMM_WORLD,
                                           &min_exch_rank_send_ptb,
                                           &min_exch_rank_recv_ptb,
                                           &max_exch_rank_send_ptb,
                                           &max_exch_rank_recv_ptb,
                                           &min_exch_data_send_ptb,
                                           &min_exch_data_recv_ptb,
                                           &max_exch_data_send_ptb,
                                           &max_exch_data_recv_ptb);

    int min_exch_rank_send_btp;
    int min_exch_rank_recv_btp;
    int max_exch_rank_send_btp;
    int max_exch_rank_recv_btp;
    unsigned long long min_exch_data_send_btp;
    unsigned long long min_exch_data_recv_btp;
    unsigned long long max_exch_data_send_btp;
    unsigned long long max_exch_data_recv_btp;

    PDM_block_to_part_global_statistic_get(PDM_MPI_COMM_WORLD,
                                           &min_exch_rank_send_btp,
                                           &min_exch_rank_recv_btp,
                                           &max_exch_rank_send_btp,
                                           &max_exch_rank_recv_btp,
                                           &min_exch_data_send_btp,
                                           &min_exch_data_recv_btp,
                                           &max_exch_data_send_btp,
                                           &max_exch_data_recv_btp);

    if (rank == 0) {
      printf("Global time in PDM_part_to_block : \n");
      printf("   - ptb min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_ptb,
             max_elaps_create_ptb);
      printf("   - ptb min max elaps create2 : %12.5e %12.5e\n",
             min_elaps_create2_ptb,
             max_elaps_create2_ptb);
      printf("   - ptb min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_ptb,
             max_elaps_exch_ptb);

      printf("   - ptb min max send rank     : %d %d\n",
             min_exch_rank_send_ptb,
             max_exch_rank_send_ptb);
      printf("   - ptb min max recv rank     : %d %d\n",
             min_exch_rank_recv_ptb,
             max_exch_rank_recv_ptb);
      printf("   - ptb min max send data     : %llu %llu\n",
             min_exch_data_send_ptb,
             max_exch_data_send_ptb);
      printf("   - ptb min max recv data     : %llu %llu\n",
             min_exch_data_recv_ptb,
             max_exch_data_recv_ptb);
      fflush(stdout);

      printf("Global time in PDM_block_to_part : \n");
      printf("   - btp min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_btp,
             max_elaps_create_btp);
      printf("   - btp min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_btp,
             max_elaps_exch_btp);

      printf("   - btp min max send rank     : %d %d\n",
             min_exch_rank_send_btp,
             max_exch_rank_send_btp);
      printf("   - btp min max recv rank     : %d %d\n",
             min_exch_rank_recv_btp,
             max_exch_rank_recv_btp);
      printf("   - btp min max send data     : %llu %llu\n",
             min_exch_data_send_btp,
             max_exch_data_send_btp);
      printf("   - btp min max recv data     : %llu %llu\n",
             min_exch_data_recv_btp,
             max_exch_data_recv_btp);
      fflush(stdout);

      char output_filename_ptbbtp[50];
      FILE* output_ptbbtp = NULL;

      if (filedump) {
        sprintf(output_filename_ptbbtp, "ptbbtp_%s", output_filename);
        output_ptbbtp = fopen(output_filename_ptbbtp, "w");

        fprintf(output_ptbbtp,
                "%d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
                comm_world_size,
                min_elaps_create_ptb,
                max_elaps_create_ptb,
                min_elaps_create2_ptb,
                max_elaps_create2_ptb,
                min_elaps_exch_ptb,
                max_elaps_exch_ptb);


//      fprintf(output_ptbbtp,
//              "%d %d %d %d %llu %llu %llu %llu ",
//              min_exch_rank_send_ptb,
//              max_exch_rank_send_ptb,
//              min_exch_rank_recv_ptb,
//              max_exch_rank_recv_ptb,
//              min_exch_data_send_ptb,
//              max_exch_data_send_ptb,
//              min_exch_data_recv_ptb,
//              max_exch_data_recv_ptb);

      fprintf(output_ptbbtp,
              "%12.5e %12.5e %12.5e %12.5e\n",
              min_elaps_create_btp,
              max_elaps_create_btp,
              min_elaps_exch_btp,
              max_elaps_exch_btp);

//      fprintf(output_ptbbtp,
//              "%d %d %d %d %llu %llu %llu %llu ",
//              min_exch_rank_send_btp,
//              max_exch_rank_send_btp,
//              min_exch_rank_recv_btp,
//              max_exch_rank_recv_btp,
//              min_exch_data_send_btp,
//              max_exch_data_send_btp,
//              min_exch_data_recv_btp,
//              max_exch_data_recv_btp);

        fclose(output_ptbbtp);
      }
    }
  }

  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);
  double geom_time = t_end - t_start;
  double max_geom_time;
  MPI_Reduce(&geom_time, &max_geom_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  FILE *output;
  if (filedump) {
    output = fopen(output_filename, "w");
  } 

  if (rank == 0) {
    printf("\n\nGeometric algorithm :%12.5es\n", max_geom_time);

    if (filedump) {
      fprintf(output, "%d %12.5e ", comm_world_size, max_geom_time);
    }
  }

  if (verbose && rank == 0) {
    printf("Geometric algorithm OK\n");

  }

  //  Exchange interpolated fields 1
  MPI_Barrier(MPI_COMM_WORLD);

  PDM_timer_hang_on(timer);
  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  int request;
  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 1, 1, 0.1, field_name, send_val, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 1, 1, 0.1, field_name, recv_val, &request);
    }
  }

  else {
    if (code_id == 1) {
      CWP_Field_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_irecv(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  double exch_time1 = t_end - t_start;
  double max_exch_time1;
  t_start = t_end;
  MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  PDM_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
      CWP_Field_wait_issend(code_name[0], coupling_name, field_name);
    }
    else {
      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  double exch_time = t_end - t_start;
  double max_exch_time;
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Exchange 1 fields issend/irecv: %12.5es\n", max_exch_time1);
    printf("Exchange 1 fields wait        : %12.5es\n", max_exch_time);
    printf("Total Exchange 1              : %12.5es\n", max_exch_time1 + max_exch_time);

    if (filedump) {
      fprintf(output, "%12.5e ", max_exch_time1 + max_exch_time);
    }
  }

  double redondance_geom = max_exch_time1;
  max_geom_time += max_exch_time1;

  //  Exchange interpolated fields 2
  PDM_timer_resume(timer);
  MPI_Barrier(MPI_COMM_WORLD);

  PDM_timer_hang_on(timer);
  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_issend(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, send_val, &request);
    }
    else {
      cwipi_irecv(coupling_name, "ech", 0, 1, 1, 0.1, field_name2, recv_val, &request);
    }
  }

  else {
    if (code_id == 1) {
//      CWP_Field_issend(code_name[0], coupling_name, field_name2);
    }
    else {
//      CWP_Field_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  exch_time1 = t_end - t_start;
  CWP_UNUSED (max_exch_time1);
  t_start = t_end;
  MPI_Reduce(&exch_time1, &max_exch_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  PDM_timer_resume(timer);

  redondance_geom += -max_exch_time1;

  if (version == CWP_VERSION_OLD) {
    if (code_id == 1) {
      cwipi_wait_issend(coupling_name, request);
    }
    else {
      cwipi_wait_irecv(coupling_name, request);
    }
  }
  else {
    if (code_id == 1) {
//      CWP_Field_wait_issend(code_name[0], coupling_name, field_name2);
    }
    else {
//      CWP_Field_wait_irecv(code_name[0], coupling_name, field_name2);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);
  exch_time = t_end - t_start;
  CWP_UNUSED (max_exch_time1);
  MPI_Reduce(&exch_time, &max_exch_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Exchange 2 fields issend/irecv :%12.5es\n", max_exch_time1);
    printf("Exchange 2  fields wait        :%12.5es\n", max_exch_time);
    printf("Total exchange 2               :%12.5es\n", max_exch_time1 + max_exch_time);

    printf("\n\nTemps geometrie                            : %12.5es\n", max_geom_time);
    printf("Temps geometrie escompte (sans redondance) : %12.5es\n",
           max_geom_time - redondance_geom);
    printf("Temps un Echange aux noeuds                : %12.5es\n",
           max_exch_time1 + max_exch_time);

    if (filedump) {
      fprintf(output, "%12.5e\n", max_exch_time1 + max_exch_time);
      fclose(output);
    }
  }

  //  Check

  if (1) {
    double max_err = 0.;
    if (code_id == 2) {
      for (int i = 0 ; i < n_located ; i++) {
        double err = ABS (recv_val[i] - vtxCoord[0][3 * (located[i] - 1)]);
        if (err > 1.e-5) {
          printf("[%d] !! vtx %ld %d err = %g (x = %f, recv = %f)\n",
          rank, vtxLNToGN[0][(located[i] - 1)], located[i], err, vtxCoord[0][3*(located[i] - 1)], recv_val[i]);
        }
        if (err > max_err) {
          max_err = err;
        }
      }
    }

    double global_max_err = 0.;
    MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      printf("Max error = %g\n", global_max_err);
    }
  }

  //  Delete interface mesh
  if (version == CWP_VERSION_NEW) {
    CWP_Mesh_interf_del(code_name[0], coupling_name);
  }
  //  Delete coupling
  if (version == CWP_VERSION_OLD) {
    cwipi_delete_coupling(coupling_name);
  }
  else {
    CWP_Cpl_del(code_name[0], coupling_name);
  }
  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(intra_comm);


  for (int ipart = 0; ipart < n_part; ipart++) {
    if (edgeVtxIdx[ipart] != NULL ) free (edgeVtxIdx[ipart]);
    if (edgeVtx[ipart] != NULL    ) free (edgeVtx[ipart]);
    if (faceEdgeIdx[ipart] != NULL) free (faceEdgeIdx[ipart]);
    if (faceEdge[ipart] != NULL   ) free (faceEdge[ipart]);
    if (faceVtxIdx[ipart] != NULL ) free (faceVtxIdx[ipart]);
    if (faceVtx[ipart] != NULL    ) free (faceVtx[ipart]);
    if (faceLNToGN[ipart] != NULL ) free (faceLNToGN[ipart]);
    if (vtxCoord[ipart] != NULL   ) free (vtxCoord[ipart]);
    if (vtxLNToGN[ipart] != NULL  ) free (vtxLNToGN[ipart]);
  }


  free(nFace);
  free(faceVtxIdx);
  free(faceVtx);
  free(faceEdgeIdx);
  free(faceEdge);
  free(faceLNToGN);
  free(nEdge);
  free(edgeVtxIdx);
  free(edgeVtx);
  free(nVtx);
  free(vtxCoord);
  free(vtxLNToGN);

  if (code_id == 1) {
    free(send_val);
  }
  else {
    free(recv_val);
  }
  PDM_timer_free(timer);

  //  Finalize CWIPI
  if (version == CWP_VERSION_OLD) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }

  if (version == CWP_VERSION_NEW) {
    double min_elaps_create_ptb;
    double max_elaps_create_ptb;
    double min_cpu_create_ptb;
    double max_cpu_create_ptb;
    double min_elaps_create2_ptb;
    double max_elaps_create2_ptb;
    double min_cpu_create2_ptb;
    double max_cpu_create2_ptb;
    double min_elaps_exch_ptb;
    double max_elaps_exch_ptb;
    double min_cpu_exch_ptb;
    double max_cpu_exch_ptb;

    PDM_part_to_block_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_ptb,
                                       &max_elaps_create_ptb,
                                       &min_cpu_create_ptb,
                                       &max_cpu_create_ptb,
                                       &min_elaps_create2_ptb,
                                       &max_elaps_create2_ptb,
                                       &min_cpu_create2_ptb,
                                       &max_cpu_create2_ptb,
                                       &min_elaps_exch_ptb,
                                       &max_elaps_exch_ptb,
                                       &min_cpu_exch_ptb,
                                       &max_cpu_exch_ptb);

    double min_elaps_create_btp;
    double max_elaps_create_btp;
    double min_cpu_create_btp;
    double max_cpu_create_btp;
    double min_elaps_exch_btp;
    double max_elaps_exch_btp;
    double min_cpu_exch_btp;
    double max_cpu_exch_btp;

    PDM_block_to_part_global_timer_get(PDM_MPI_COMM_WORLD,
                                       &min_elaps_create_btp,
                                       &max_elaps_create_btp,
                                       &min_cpu_create_btp,
                                       &max_cpu_create_btp,
                                       &min_elaps_exch_btp,
                                       &max_elaps_exch_btp,
                                       &min_cpu_exch_btp,
                                       &max_cpu_exch_btp);

    if (rank == 0) {
      printf("Global time in PDM_part_to_block : \n");
      printf("   - ptb min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_ptb,
             max_elaps_create_ptb);
      printf("   - ptb min max elaps create2 : %12.5e %12.5e\n",
             min_elaps_create2_ptb,
             max_elaps_create2_ptb);
      printf("   - ptb min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_ptb,
             max_elaps_exch_ptb);
      fflush(stdout);

      printf("Global time in PDM_block_to_part : \n");
      printf("   - btp min max elaps create  : %12.5e %12.5e\n",
             min_elaps_create_btp,
             max_elaps_create_btp);
      printf("   - btp min max elaps exch    : %12.5e %12.5e\n",
             min_elaps_exch_btp,
             max_elaps_exch_btp);
      fflush(stdout);
    }
  }

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
