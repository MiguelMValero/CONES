#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                         argc,
           char                      **argv,
           PDM_g_num_t                *n_vtx_seg,
           double                     *length,
           double                     *depth,
           int                        *rotation,
           double                     *tolerance,
           double                     *marge,
           int                        *n_part,
           PDM_g_num_t                *n_pts,
           int                        *n_proc_data,
           int                        *post,
           int                        *have_random,
           int                        *init_random,
           int                        *part_method,
           PDM_mesh_location_method_t *loc_method,
           int                        *use_mesh_vtx)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-d") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *depth = atof(argv[i]);
    }
    else if (strcmp (argv[i], "-rot") == 0) {
      *rotation = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_pts = atol(argv[i]);
        *n_pts = (PDM_g_num_t) _n_pts;
      }
    }

    else if (strcmp (argv[i], "-no_random") == 0) {
      *have_random = 0;
    }

    else if (strcmp (argv[i], "-use_vtx") == 0) {
      *use_mesh_vtx = 1;
    }

    else if (strcmp (argv[i], "-random_init") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *init_random = atoi(argv[i]);
      }
    }

    else if (strcmp (argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_proc_data = atoi(argv[i]);
      }
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
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_point_cloud_from_mesh_vtx
(
 PDM_MPI_Comm   pdm_mpi_comm,
 double         xmin,
 double         ymin,
 PDM_g_num_t    nVtxSeg,
 double         length,
 int            haveRandom,
 int            initRandom,
 double       **coord,
 int           *n_pts_l
 )
{
  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

  double       xmax = xmin + length;
  double       ymax = ymin + length;
  PDM_g_num_t  nx = nVtxSeg;
  PDM_g_num_t  ny = nVtxSeg;

  int          dNFace;
  int          dNVtx;
  int          dNEdge;
  int         *dFaceVtxIdx;
  PDM_g_num_t *dFaceVtx;
  double      *dVtxCoord;
  PDM_g_num_t *dFaceEdge;
  PDM_g_num_t *dEdgeVtx;
  PDM_g_num_t *dEdgeFace;
  int          nEdgeGroup;
  int         *dEdgeGroupIdx;
  PDM_g_num_t *dEdgeGroup;

  /*
   *  Create mesh
   */
  PDM_g_num_t nGFace;
  PDM_g_num_t nGEdge;
  PDM_g_num_t nGVtx;

  PDM_poly_surf_gen (pdm_mpi_comm,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     haveRandom,
                     initRandom,
                     nx,
                     ny,
                     &nGFace,
                     &nGVtx,
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

  *n_pts_l = dNVtx;
  *coord   = dVtxCoord;

  free (dFaceVtxIdx);
  free (dFaceVtx);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);
}


static void _add_depth (const int     n_pts,
                        const double  length,
                        const double  depth,
                        double       *coord)
{
  double inv_length = 1.;
  if (PDM_ABS (length) > 1e-15) inv_length /= length;

  for (int i = 0; i < n_pts; i++) {
    double x = 2.*coord[3*i]   * inv_length;
    double y = 2.*coord[3*i+1] * inv_length;
    coord[3*i+2] = 0.5*depth*(1. - (x*x + y*y));
  }
}


static void _rotate (const int  n_pts,
                     double    *coord)
{
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}


static int _set_rank_has_mesh
(
 const PDM_MPI_Comm  comm,
 const int           nProcData,
 PDM_MPI_Comm       *meshComm
 )
{
  int current_rank_has_mesh = 1;

  int rank;
  int commSize;

  PDM_MPI_Comm_rank (comm, &rank);
  PDM_MPI_Comm_size (comm, &commSize);

  if (nProcData > 0 && nProcData < commSize) {

    int rankInNode = PDM_io_mpi_node_rank (comm);

    int nNode = 0;
    int iNode = -1;
    int masterRank = (rankInNode == 0);

    int *rankInNodes = malloc(sizeof(int) * commSize);

    PDM_MPI_Allreduce (&masterRank, &nNode, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    PDM_MPI_Allgather (&rankInNode, 1, PDM_MPI_INT, rankInNodes, 1, PDM_MPI_INT, comm);

    current_rank_has_mesh = 0;

    for (int i = 0; i < rank; i++) {
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

    PDM_MPI_Comm_split (comm,
                        current_rank_has_mesh,
                        rank,
                        meshComm);
    free (rankInNodes);
  }

  return current_rank_has_mesh;
}


static void
_get_connectivity
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **nFace,
 int         ***faceEdgeIdx,
 int         ***faceEdge,
 int         ***faceVtxIdx,
 int         ***faceVtx,
 PDM_g_num_t ***faceLNToGN,
 int          **nEdge,
 int         ***edgeVtxIdx,
 int         ***edgeVtx,
 int          **nVtx,
 double      ***vtxCoord,
 PDM_g_num_t ***vtxLNToGN
 )
{
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

   // int id_ppart = ppartId;

  for (int ipart = 0; ipart < n_part; ipart++) {

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

    PDM_part_part_dim_get (ppart,
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

    int         *_faceTag;
    int         *_faceEdgeIdx;
    int         *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int         *_edgeTag;
    int         *_edgeFace;
    int         *_edgeVtxIdx;
    int         *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int         *_edgePartBoundProcIdx;
    int         *_edgePartBoundPartIdx;
    int         *_edgePartBound;
    int         *_vtxTag;
    double      *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int         *_edgeGroupIdx;
    int         *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
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

    /* Faces */
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

    /* Edges */
    (*nEdge)[ipart] = _nEdge;
    (*edgeVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
    (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

    memcpy ((*edgeVtxIdx)[ipart], _edgeVtxIdx, (_nEdge + 1) * sizeof(int));
    memcpy ((*edgeVtx)[ipart], _edgeVtx, _sEdgeVtx * sizeof(int));

    /* Vertices */
    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));


    /* Compute face-vtx connectivity */
    int *_faceVtx = (*faceVtx)[ipart];

    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));

    for (int i = 0; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i];
      int ivtx2 = _edgeVtx[2*i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
    for (int i = 0; i < _nVtx; i++) {
      vtxEdgeN[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i] - 1;
      int ivtx2 = _edgeVtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
      vtxEdgeN[ivtx1] += 1;
      vtxEdgeN[ivtx2] += 1;
    }
    free(vtxEdgeN);

    for (int i = 0; i < _nFace; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i+1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _faceVtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edgeVtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edgeVtx[2*(edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtxEdgeIdx[vtx_cur - 1]; j <  vtxEdgeIdx[vtx_cur]; j++) {
          for (int k = 0; k < __nEdge; k++) {
            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edgeVtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edgeVtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edgeVtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtxedge !!!!\n");
          abort();
        }
      }
    }

    free (vtxEdge);
    free (vtxEdgeIdx);

  }
}


static void
_create_split_mesh
(
 int                 activeRank,
 PDM_MPI_Comm        pdm_mpi_comm,
 double              xmin,
 double              ymin,
 PDM_g_num_t         nVtxSeg,
 double              length,
 double              depth,
 int                 rotation,
 int                 n_part,
 PDM_part_split_t    method,
 int                 haveRandom,
 int                 initRandom,
 PDM_g_num_t        *nGFace,
 PDM_g_num_t        *nGVtx,
 int               **nFace,
 int              ***faceEdgeIdx,
 int              ***faceEdge,
 int              ***faceVtxIdx,
 int              ***faceVtx,
 PDM_g_num_t      ***faceLNToGN,
 int               **nEdge,
 int              ***edgeVtxIdx,
 int              ***edgeVtx,
 int               **nVtx,
 double           ***vtxCoord,
 PDM_g_num_t      ***vtxLNToGN
 )
{
  int i_rank;
  int numProcs;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

    double       xmax = xmin + length;
    double       ymax = ymin + length;
    PDM_g_num_t  nx = nVtxSeg;
    PDM_g_num_t  ny = nVtxSeg;

    int          dNFace;
    int          dNVtx;
    int          dNEdge;
    int         *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    double      *dVtxCoord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int          nEdgeGroup;
    int         *dEdgeGroupIdx;
    PDM_g_num_t *dEdgeGroup;

    /*
     *  Create mesh
     */
    PDM_g_num_t nGEdge;

    PDM_poly_surf_gen (pdm_mpi_comm,
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

    _add_depth (dNVtx,
                length,
                depth,
                dVtxCoord);

    if (rotation) {
      _rotate (dNVtx,
               dVtxCoord);
    }

    /*
     *  Create mesh partitions
     */
    int have_dCellPart = 0;

    int *dCellPart   = (int *) malloc (dNFace * sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc ((dNEdge + 1) * sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh
     */
    // int ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
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

    free (dCellPart);

    free (dVtxCoord);
    free (dFaceVtxIdx);
    free (dFaceVtx);
    free (dFaceEdge);
    free (dEdgeVtxIdx);
    free (dEdgeVtx);
    free (dEdgeFace);
    free (dEdgeGroupIdx);
    free (dEdgeGroup);

    _get_connectivity (ppart,
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

    PDM_part_free (ppart);
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

    for (int ipart = 0; ipart < n_part; ipart++) {
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

  PDM_MPI_Bcast (nGFace, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (nGVtx, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
}


static inline double
_eval_field
(
 double *xyz
 )
{
  return 1 + 2*xyz[0] + 3*xyz[1] + 4*xyz[2];
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t n_vtx_seg   = 10;
  double      length      = 1.;
  double      depth       = 0.;
  int         rotation    = 0;
  double      tolerance   = 1e-6;
  double      marge       = 0.;
  int         n_part      = 1;
  int         post        = 0;
  int         have_random = 1;
  int         init_random = 0;
  int         n_proc_data = -1;
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  int use_vtx = 0;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &depth,
              &rotation,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &n_proc_data,
              &post,
              &have_random,
              &init_random,
              (int *) &part_method,
              &loc_method,
              &use_vtx);


  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank      : %d\n", n_rank);
    PDM_printf ("  - n_vtx_seg   : "PDM_FMT_G_NUM"\n", n_vtx_seg);
    PDM_printf ("  - n_pts       : "PDM_FMT_G_NUM"\n", n_pts);
    PDM_printf ("  - length      : %f\n", length);
    PDM_printf ("  - depth       : %f\n", depth);
    PDM_printf ("  - tolerance   : %f\n", tolerance);
    PDM_printf ("  - part_method : %d\n", (int) part_method);
    PDM_printf ("  - loc_method  : %d\n", (int) loc_method);
    PDM_printf ("  - n_proc_data : %d\n", n_proc_data);
  }

  /*
   *  Create partitionned surface mesh
   */

  if (i_rank == 0) {
    printf("-- Build surface mesh\n");
    fflush(stdout);
  }

  PDM_MPI_Comm mesh_comm = PDM_MPI_COMM_WORLD;
  int current_rank_has_mesh = _set_rank_has_mesh (PDM_MPI_COMM_WORLD,
                                                  n_proc_data,
                                                  &mesh_comm);

  const double xmin = -0.5*length;
  const double ymin = -0.5*length;

  PDM_g_num_t   nGFace;
  PDM_g_num_t   nGVtx;
  int          *nFace       = NULL;
  PDM_g_num_t **faceLNToGN  = NULL;
  int         **faceEdgeIdx = NULL;
  int         **faceEdge    = NULL;
  int         **faceVtxIdx  = NULL;
  int         **faceVtx     = NULL;
  int          *nEdge       = NULL;
  int         **edgeVtxIdx  = NULL;
  int         **edgeVtx     = NULL;
  int          *nVtx        = NULL;
  double      **vtxCoord    = NULL;
  PDM_g_num_t **vtxLNToGN   = NULL;

  _create_split_mesh (current_rank_has_mesh,
                      mesh_comm,
                      xmin,
                      ymin,
                      n_vtx_seg,
                      length,
                      depth,
                      rotation,
                      n_part,
                      part_method,
                      have_random,
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


  /************************
   *
   * Point cloud definition
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double *pts_coords    = NULL;
  PDM_g_num_t *pts_gnum = NULL;

  marge *= length;
  double _min = -0.5*length - marge;
  double _max = -_min;

  if (use_vtx) {
    _point_cloud_from_mesh_vtx (PDM_MPI_COMM_WORLD,
                                _min,
                                _min,
                                n_pts,
                                length + 2.*marge,
                                have_random,
                                init_random + 1,
                                &pts_coords,
                                &n_pts_l);
  }
  else {
    PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                                0, // seed
                                0, // geometric_g_num
                                n_pts,
                                _min,
                                _min,
                                0.,
                                _max,
                                _max,
                                0.,
                                &n_pts_l,
                                &pts_coords,
                                &pts_gnum);
  }

  _add_depth (n_pts_l,
              length,
              depth,
              pts_coords);


  /* Point cloud global numbering */
  if (use_vtx) {
    if (i_rank == 0) {
      printf("-- Point cloud g_num\n");
      fflush(stdout);
    }
  #if 1
    PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

    double *char_length = malloc(sizeof(double) * n_pts_l);

    for (int i = 0; i < n_pts_l; i++) {
      char_length[i] = length * 1.e-9;
    }

    PDM_gnum_set_from_coords (gen_gnum, 0, n_pts_l, pts_coords, char_length);

    if (i_rank == 0) {
      printf(">> PDM_gnum_compute\n");
      fflush(stdout);
    }
    PDM_gnum_compute (gen_gnum);
    if (i_rank == 0) {
      printf("<< PDM_gnum_compute\n");
      fflush(stdout);
    }

    pts_gnum = PDM_gnum_get(gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
    free (char_length);
  #else
  PDM_g_num_t *distrib = PDM_compute_entity_distribution (PDM_MPI_COMM_WORLD,
                                                          n_pts_l);
  PDM_g_num_t *pts_gnum = malloc (sizeof(PDM_g_num_t) * n_pts_l);
  for (int i = 0; i < n_pts_l; i++) {
    pts_gnum[i] = distrib[i_rank] + i + 1;
  }

  free (distrib);
  #endif
  }


  if (post) {
    char filename[999];
    sprintf(filename, "point_cloud_%3.3d.vtk", i_rank);

    PDM_vtk_write_point_cloud (filename,
                               n_pts_l,
                               pts_coords,
                               pts_gnum,
                               NULL);
  }


  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Create mesh loc\n");
    fflush(stdout);
  }

  PDM_mesh_location_t* mesh_loc = PDM_mesh_location_create (1,//const int n_point_cloud,
                                                            PDM_MPI_COMM_WORLD,
                                                            PDM_OWNERSHIP_KEEP);

  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (mesh_loc,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coords,
                               pts_gnum);

  PDM_mesh_location_mesh_n_part_set (mesh_loc,
                                          n_part);

  /* Set mesh */
  if (i_rank == 0) {
    printf("-- Set mesh\n");
    fflush(stdout);
  }
  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_part_set_2d (mesh_loc,
                                   ipart,
                                   nFace[ipart],
                                   faceEdgeIdx[ipart],
                                   faceEdge[ipart],
                                   faceLNToGN[ipart],
                                   nEdge[ipart],
                                   edgeVtx[ipart],
                                   nVtx[ipart],
                                   vtxCoord[ipart],
                                   vtxLNToGN[ipart]);
  }




  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);


  /*
   * Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  // PDM_mesh_location_compute (mesh_loc);
  PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);





  /*
   * Check results
   */
  if (i_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                   0,//i_point_cloud,
                                                   0);//i_part,

  int *located = PDM_mesh_location_located_get (mesh_loc,
                                                0,//i_point_cloud,
                                                0);

  int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                       0,//i_point_cloud,
                                                       0);

  int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                    0,//i_point_cloud,
                                                    0);

  PDM_g_num_t *p_location    = NULL;
  double      *p_dist2  = NULL;
  double      *p_proj_coord  = NULL;
  PDM_mesh_location_point_location_get (mesh_loc,
                                        0,//i_point_cloud,
                                        0,//i_part,
                                        &p_location,
                                        &p_dist2,
                                        &p_proj_coord);

  if (0) {
    printf("Unlocated %d :\n", n_unlocated);
    for (int k1 = 0; k1 < n_unlocated; k1++) {
      printf("%d\n", unlocated[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      printf("%d\n", located[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      printf(PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e",
        pts_gnum[ipt],  p_location[k1],
        pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
        p_dist2[k1],
        p_proj_coord[3*k1], p_proj_coord[3*k1+1], p_proj_coord[3*k1+2]);
      printf("\n");
    }
  }


  // /* int p_n_points; */
  // PDM_g_num_t *p_location    = NULL;
  // int         *p_weights_idx = NULL;
  // double      *p_weights     = NULL;
  // double      *p_proj_coord  = NULL;

  // PDM_mesh_location_get (mesh_loc,
  //                        0,//i_point_cloud,
  //                        0,//i_part,
  //                        &p_location,
  //                        &p_weights_idx,
  //                        &p_weights,
  //                        &p_proj_coord);

  if (0) {
    printf("Unlocated %d :\n", n_unlocated);
    for (int k1 = 0; k1 < n_unlocated; k1++) {
      printf("%d\n", unlocated[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      printf("%d\n", located[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      printf(PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e",
        pts_gnum[ipt],  p_location[k1],
        pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
        p_dist2[k1],
        p_proj_coord[3*k1], p_proj_coord[3*k1+1], p_proj_coord[3*k1+2]);
      printf("\n");
    }
  }

  //  if (0) {
  //   for (int ipt = 0; ipt < n_pts_l; ipt++) {
  //     printf("Point ("PDM_FMT_G_NUM") (%f %f %f), location = ("PDM_FMT_G_NUM"), proj = (%f %f %f), weights =",
  //            pts_gnum[ipt],
  //            pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
  //            p_location[ipt],
  //            p_proj_coord[3*ipt], p_proj_coord[3*ipt+1], p_proj_coord[3*ipt+2]);
  //     for (int i = p_weights_idx[ipt]; i < p_weights_idx[ipt+1]; i++) {
  //       printf(" %f", p_weights[i]);
  //     }
  //     printf("\n");
  //   }
  // }




  /*
   *  Check location (interpolation of an affine field)
   */
  double **src_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    src_field[ipart] = malloc(sizeof(double) * nVtx[ipart]);
    for (int i = 0; i < nVtx[ipart]; i++) {
      src_field[ipart][i] = _eval_field(&vtxCoord[ipart][3*i]);
    }
  }


  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);
  if (ptp == NULL) {
    int         **pelt_pts_idx  = malloc(sizeof(int         *) * n_part);
    PDM_g_num_t **pelt_pts_gnum = malloc(sizeof(PDM_g_num_t *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      double *elt_pts_coord      = NULL;
      double *elt_pts_uvw        = NULL;
      int    *elt_pts_weight_idx = NULL;
      double *elt_pts_weight     = NULL;
      double *elt_pts_dist2      = NULL;
      double *elt_pts_proj_coord = NULL;
      PDM_mesh_location_points_in_elt_get(mesh_loc,
                                          ipart,
                                          0, // i_point_cloud,
                                          &pelt_pts_idx [ipart],
                                          &pelt_pts_gnum[ipart],
                                          &elt_pts_coord,
                                          &elt_pts_uvw,
                                          &elt_pts_weight_idx,
                                          &elt_pts_weight,
                                          &elt_pts_dist2,
                                          &elt_pts_proj_coord);
    }

    ptp = PDM_part_to_part_create((const PDM_g_num_t **) faceLNToGN,
                                  (const int          *) nFace,
                                  n_part,
                                  (const PDM_g_num_t **) &pts_gnum,
                                  (const int          *) &n_pts_l,
                                  1,
                                  (const int         **) pelt_pts_idx,
                                  (const PDM_g_num_t **) pelt_pts_gnum,
                                  PDM_MPI_COMM_WORLD);

    free(pelt_pts_idx );
    free(pelt_pts_gnum);
  }


  double **send_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int         *elt_pts_idx        = NULL;
    PDM_g_num_t *elt_pts_gnum       = NULL;
    double      *elt_pts_coord      = NULL;
    double      *elt_pts_uvw        = NULL;
    int         *elt_pts_weight_idx = NULL;
    double      *elt_pts_weight     = NULL;
    double      *elt_pts_dist2      = NULL;
    double      *elt_pts_proj_coord = NULL;
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        ipart,
                                        0, // i_point_cloud,
                                        &elt_pts_idx,
                                        &elt_pts_gnum,
                                        &elt_pts_coord,
                                        &elt_pts_uvw,
                                        &elt_pts_weight_idx,
                                        &elt_pts_weight,
                                        &elt_pts_dist2,
                                        &elt_pts_proj_coord);

    send_field[ipart] = malloc(sizeof(double) * elt_pts_idx[nFace[ipart]]);
    for (int ielt = 0; ielt < nFace[ipart]; ielt++) {

      int *fv = faceVtx[ipart] + faceVtxIdx[ipart][ielt];

      for (int idx_pt = elt_pts_idx[ielt]; idx_pt < elt_pts_idx[ielt+1]; idx_pt++) {
        send_field[ipart][idx_pt] = 0.;
        int idx_vtx = 0;
        // double e[3] = {
        //   elt_pts_proj_coord[3*idx_pt  ],
        //   elt_pts_proj_coord[3*idx_pt+1],
        //   elt_pts_proj_coord[3*idx_pt+2]
        // };
        double e[3] = {
          elt_pts_coord[3*idx_pt  ],
          elt_pts_coord[3*idx_pt+1],
          elt_pts_coord[3*idx_pt+2]
        };

        for (int idx_w = elt_pts_weight_idx[idx_pt]; idx_w < elt_pts_weight_idx[idx_pt+1]; idx_w++) {
          int vtx_id = fv[idx_vtx++] - 1;
          send_field[ipart][idx_pt] += elt_pts_weight[idx_w] * src_field[ipart][vtx_id];
          for (int j = 0; j < 3; j++) {
            e[j] -= elt_pts_weight[idx_w] * vtxCoord[ipart][3*vtx_id+j];
          }
        }

        // log_trace("pt "PDM_FMT_G_NUM" (%f %f %f), in elt "PDM_FMT_G_NUM" : dist = %e\n",
        //           elt_pts_gnum[idx_pt],
        //           elt_pts_coord[3*idx_pt], elt_pts_coord[3*idx_pt+1], elt_pts_coord[3*idx_pt+2],
        //           faceLNToGN[ipart][ielt],
        //           PDM_MODULE(e));
      }

    }
  }


  double **recv_field = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) send_field,
                         NULL,
        (      void ***) &recv_field,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);

  double lmax_err = 0.;
  for (int i = 0; i < n_located; i++) {
    int pt_id = located[i] - 1;

    double f = _eval_field(&pts_coords[3*pt_id]);

    double err = PDM_ABS(recv_field[0][i] - f);
    lmax_err = PDM_MAX(lmax_err, err);

    if (err > 1.e-12) {
      log_trace("point "PDM_FMT_G_NUM" (%f %f %f) located in elt "PDM_FMT_G_NUM" : error = %e (%20.16f / %20.16f)\n",
                pts_gnum[pt_id],
                pts_coords[3*pt_id], pts_coords[3*pt_id+1], pts_coords[3*pt_id+2],
                p_location[i], err, recv_field[0][i], f);
    }
  }

  free(recv_field[0]);
  free(recv_field);


  double gmax_err;
  PDM_MPI_Allreduce(&lmax_err, &gmax_err, 1, PDM_MPI_DOUBLE,
                    PDM_MPI_MAX, PDM_MPI_COMM_WORLD);


  if (i_rank == 0) {
    printf("global max interpolation error = %e\n", gmax_err);
  }


  /*
   * Finalize
   */
  PDM_mesh_location_free(mesh_loc);
  PDM_part_to_part_free (ptp);
                          

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(faceEdgeIdx[ipart]);
    free(faceEdge[ipart]);
    free(faceVtxIdx[ipart]);
    free(faceVtx[ipart]);
    free(faceLNToGN[ipart]);
    free(edgeVtxIdx[ipart]);
    free(edgeVtx[ipart]);
    free(vtxCoord[ipart]);
    free(vtxLNToGN[ipart]);

    free(src_field[ipart]);
    free(send_field[ipart]);
  }
  free(faceVtxIdx);
  free(faceVtx);
  free(nFace);
  free(faceEdgeIdx);
  free(faceEdge);
  free(faceLNToGN);
  free(nEdge);
  free(edgeVtxIdx);
  free(edgeVtx);
  free(nVtx);
  free(vtxCoord);
  free(vtxLNToGN);

  free(src_field);
  free(send_field);
  /*PDM_part_free (ppart_id);


  free (dvtx_coord);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dFaceEdge);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);
  free (dEdgeVtxIdx);*/

  free (pts_coords);
  free (pts_gnum);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
