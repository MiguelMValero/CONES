#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part.h"

#include "pdm_writer.h"
#include "pdm_geom_elem.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_priv.h"

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA       <value>  Cube A - Number of vertices on the side.\n\n"
     "  -lA       <value>  Cube A - length.\n\n"
     "  -xminA    <value>  Cube A - Origin X.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -n_partA  <level>  Cube A - Number of partitions par process.\n\n"
     "  -nB       <value>  Cube A - Number of vertices on the side (default : nA+4).\n\n"
     "  -lB       <value>  Cube A - length.\n\n"
     "  -xminB    <value>  Cube A - Origin X.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -n_partB  <level>  Cube A - Number of partitions par process.\n\n"
     "  -post              Ensight output.\n\n"
     "  -no_random         No random to define points coordinates.\n\n"
     "  -random_time_init  Intialize random with the current time.\n\n"
     "  -parmetis          Call ParMETIS.\n\n"
     "  -pt-scocth         Call PT-Scotch.\n\n"
     "  -n_proc_data       Number of processses where there are data.\n\n"
     "  -h                 This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtxSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t  *n_vtxSegA,
 double        *lengthA,
 double        *xminA,
 double        *yminA,
 int           *n_partA,
 int           *post,
 int           *method,
 int           *haveRandom,
 int           *randomTimeInit,
 int           *randomMeshAInit,
 int           *nProcData
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-nA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtxSeg = atol (argv[i]);
        *n_vtxSegA = (PDM_g_num_t) _n_vtxSeg;
      }
    }
    else if (strcmp (argv[i], "-lA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *lengthA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-xminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *xminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-yminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *yminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-n_partA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_partA = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-randomMeshAInit") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *randomMeshAInit = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-no_random") == 0) {
      *haveRandom = 0;
    }
    else if (strcmp (argv[i], "-random_time_init") == 0) {
      *randomTimeInit = 1;
    }
    else if (strcmp (argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp (argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp (argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else if (strcmp (argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *nProcData = atoi (argv[i]);
      }
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Compute face_vtx connectivity
 *
 * \param [in]      ppartId  ppart identifier
 * \param [in]      n_part    Number of partitions
 *
 * \return          face_vtx connectivity for each partition of each mesh
 */

static void
_compute_face_vtx
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **n_face,
 int         ***face_vtx_idx,
 int         ***face_vtx,
 PDM_g_num_t ***face_ln_to_gn,
 int          **n_vtx,
 double      ***vtx_coord,
 PDM_g_num_t ***vtx_ln_to_gn
)
{

  *n_face = (int *) malloc(sizeof(int) * n_part);
  *face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
  *face_vtx = (int **) malloc(sizeof(int *) * n_part);
  *face_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *n_vtx = (int *) malloc(sizeof(int) * n_part);
  *vtx_coord = (double **) malloc(sizeof(double *) * n_part);
  *vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);


  for (int ipart = 0; ipart < n_part; ipart++) {

    int _n_face;
    int _nEdge;
    int _nEdgePartBound;
    int _n_vtx;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           ipart,
                           &_n_face,
                           &_nEdge,
                           &_nEdgePartBound,
                           &_n_vtx,
                           &_nProc,
                           &_nTPart,
                           &_sFaceEdge,
                           &_sEdgeVtx,
                           &_sEdgeGroup,
                           &_nEdgeGroup2);

    int          *_faceTag;
    int          *_faceEdgeIdx;
    int          *_faceEdge;
    PDM_g_num_t *_face_ln_to_gn;
    int          *_edgeTag;
    int          *_edgeFace;
    int          *_edgeVtxIdx;
    int          *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int          *_edgePartBoundProcIdx;
    int          *_edgePartBoundPartIdx;
    int          *_edgePartBound;
    int          *_vtxTag;
    double       *_vtx;
    PDM_g_num_t *_vtx_ln_to_gn;
    int          *_edgeGroupIdx;
    int          *_edgeGroup;
    PDM_g_num_t  *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &_faceTag,
                           &_faceEdgeIdx,
                           &_faceEdge,
                           &_face_ln_to_gn,
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
                           &_vtx_ln_to_gn,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    (*n_face)[ipart] = _n_face;
    (*face_vtx_idx) [ipart] = (int *) malloc(sizeof(int) * (_n_face + 1));
    (*face_vtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*face_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_face);

    memcpy ((*face_vtx_idx)[ipart], _faceEdgeIdx, (_n_face + 1) * sizeof(int));
    memcpy ((*face_ln_to_gn)[ipart], _face_ln_to_gn, _n_face * sizeof(PDM_g_num_t));

    (*n_vtx)[ipart] = _n_vtx;
    (*vtx_coord)[ipart] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
    (*vtx_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);

    memcpy ((*vtx_coord)[ipart], _vtx, 3 *_n_vtx * sizeof(double));
    memcpy ((*vtx_ln_to_gn)[ipart], _vtx_ln_to_gn, _n_vtx * sizeof(PDM_g_num_t));

    int *_face_vtx = (*face_vtx)[ipart];

    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_n_vtx + 1));

    for (int i = 0; i < _n_vtx + 1; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i];
      int ivtx2 = _edgeVtx[2*i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1; i < _n_vtx + 1; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_n_vtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _n_vtx);
    for (int i = 0; i < _n_vtx; i++) {
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

    for (int i = 0; i < _n_face; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i+1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _face_vtx + idx;

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


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      n_vtxSeg  Number of arguments
 * \param [in]      length   Lenght of square
 * \param [in]      n_part    Number to obtain on this processus
 * \param [in]      post     mesh export status
 * \param [in]      method   Split method
 *
 * \return ppartId  ppart identifier
 *
 */

static void
_create_split_mesh
(
 int               activeRank,
 PDM_MPI_Comm      pdm_mpi_comm,
 double            xmin,
 double            ymin,
 PDM_g_num_t       n_vtxSeg,
 double            length,
 int               n_part,
 PDM_part_split_t  method,
 int               haveRandom,
 int               initRandom,
 PDM_g_num_t      *nGFace,
 PDM_g_num_t      *nGVtx,
 int             **n_face,
 int            ***face_vtx_idx,
 int            ***face_vtx,
 PDM_g_num_t    ***face_ln_to_gn,
 int             **n_vtx,
 double         ***vtx_coord,
 PDM_g_num_t    ***vtx_ln_to_gn
)
{
  struct timeval t_elaps_debut;

  int i_rank;
  int numProcs;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

    double        xmax = xmin + length;
    double        ymax = ymin + length;
    PDM_g_num_t  nx = n_vtxSeg;
    PDM_g_num_t  ny = n_vtxSeg;

    int           dn_face;
    int           dn_vtx;
    int           dNEdge;
    int          *dface_vtx_idx;
    PDM_g_num_t *dface_vtx;
    double       *dvtx_coord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int           nEdgeGroup;
    int          *dEdgeGroupIdx;
    PDM_g_num_t   *dEdgeGroup;


    /*
     *  Create mesh i
     */

    gettimeofday(&t_elaps_debut, NULL);

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
                       &dn_vtx,
                       &dvtx_coord,
                       &dn_face,
                       &dface_vtx_idx,
                       &dface_vtx,
                       &dFaceEdge,
                       &dNEdge,
                       &dEdgeVtx,
                       &dEdgeFace,
                       &nEdgeGroup,
                       &dEdgeGroupIdx,
                       &dEdgeGroup);

    struct timeval t_elaps_fin;

    gettimeofday (&t_elaps_fin, NULL);

    long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
      (t_elaps_debut.tv_usec + 1000000 *
       t_elaps_debut.tv_sec);
    long tranche_elapsed_max = tranche_elapsed;
    double t_elapsed = (double) tranche_elapsed_max/1000000.;
    if (i_rank == 0)
      PDM_printf("[%d] Temps dans creeMaillagePolygone2D : %12.5e\n",
                 i_rank, t_elapsed);

    /*
     *  Create mesh partitions
     */

    int have_dCellPart = 0;

    int *dCellPart = (int *) malloc (dn_face*sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh i
     */

    // int ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    // PDM_part_create (&ppartId,
    //                  pdm_mpi_comm,
    //                  method,
    //                  "PDM_PART_RENUM_CELL_NONE",
    //                  "PDM_PART_RENUM_FACE_NONE",
    //                  nPropertyCell,
    //                  renum_properties_cell,
    //                  nPropertyFace,
    //                  renum_properties_face,
    //                  n_part,
    //                  dn_face,
    //                  dNEdge,
    //                  dn_vtx,
    //                  nEdgeGroup,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  have_dCellPart,
    //                  dCellPart,
    //                  dEdgeFace,
    //                  dEdgeVtxIdx,
    //                  dEdgeVtx,
    //                  NULL,
    //                  dvtx_coord,
    //                  NULL,
    //                  dEdgeGroupIdx,
    //                  dEdgeGroup);

    printf("dn_face = %i | dNEdge = %i | dn_vtx = %i \n", dn_face, dNEdge, dn_vtx);
    PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
                                         method,
                                         "PDM_PART_RENUM_CELL_NONE",
                                         "PDM_PART_RENUM_FACE_NONE",
                                         nPropertyCell,
                                         renum_properties_cell,
                                         nPropertyFace,
                                         renum_properties_face,
                                         n_part,
                                         dn_face,
                                         dNEdge,
                                         dn_vtx,
                                         nEdgeGroup,
                                         dface_vtx_idx,
                                         dFaceEdge,
                                         NULL,
                                         NULL,
                                         have_dCellPart,
                                         dCellPart,
                                         NULL,
                                         dEdgeVtxIdx,
                                         dEdgeVtx,
                                         NULL,
                                         dvtx_coord,
                                         NULL,
                                         dEdgeGroupIdx,
                                         dEdgeGroup);
    free (dCellPart);

    double  *elapsed = NULL;
    double  *cpu = NULL;
    double  *cpu_user = NULL;
    double  *cpu_sys = NULL;

    PDM_part_time_get (ppart,
                       &elapsed,
                       &cpu,
                       &cpu_user,
                       &cpu_sys);

    if (i_rank == 0)
      PDM_printf("[%d] Temps dans ppart : %12.5e\n",
                 i_rank,  elapsed[0]);

    /* Statistiques */

    int    cells_average;
    int    cells_median;
    double cells_std_deviation;
    int    cells_min;
    int    cells_max;
    int    bound_part_faces_average;
    int    bound_part_faces_median;
    double bound_part_faces_std_deviation;
    int    bound_part_faces_min;
    int    bound_part_faces_max;
    int    bound_part_faces_sum;

    PDM_part_stat_get (ppart,
                       &cells_average,
                       &cells_median,
                       &cells_std_deviation,
                       &cells_min,
                       &cells_max,
                       &bound_part_faces_average,
                       &bound_part_faces_median,
                       &bound_part_faces_std_deviation,
                       &bound_part_faces_min,
                       &bound_part_faces_max,
                       &bound_part_faces_sum);

    if (i_rank == 0) {
      PDM_printf ("Statistics :\n");
      PDM_printf ("  - Number of cells :\n");
      PDM_printf ("       * average            : %i\n", cells_average);
      PDM_printf ("       * median             : %i\n", cells_median);
      PDM_printf ("       * standard deviation : %12.5e\n", cells_std_deviation);
      PDM_printf ("       * min                : %i\n", cells_min);
      PDM_printf ("       * max                : %i\n", cells_max);
      PDM_printf ("  - Number of faces exchanging with another partition :\n");
      PDM_printf ("       * average            : %i\n", bound_part_faces_average);
      PDM_printf ("       * median             : %i\n", bound_part_faces_median);
      PDM_printf ("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
      PDM_printf ("       * min                : %i\n", bound_part_faces_min);
      PDM_printf ("       * max                : %i\n", bound_part_faces_max);
      PDM_printf ("       * total              : %i\n", bound_part_faces_sum);
    }

    free (dvtx_coord);
    free (dface_vtx_idx);
    free (dface_vtx);
    free (dFaceEdge);
    free (dEdgeVtxIdx);
    free (dEdgeVtx);
    free (dEdgeFace);
    free (dEdgeGroupIdx);
    free (dEdgeGroup);

    _compute_face_vtx (ppart,
                       n_part,
                       n_face,
                       face_vtx_idx,
                       face_vtx,
                       face_ln_to_gn,
                       n_vtx,
                       vtx_coord,
                       vtx_ln_to_gn);

    PDM_part_free (ppart);

  }
  else {
    *n_face = (int *) malloc(sizeof(int) * n_part);
    *face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
    *face_vtx = (int **) malloc(sizeof(int *) * n_part);
    *face_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    *n_vtx = (int *) malloc(sizeof(int) * n_part);
    *vtx_coord = (double **) malloc(sizeof(double *) * n_part);
    *vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _n_face = 0;
      int _sFaceEdge = 0;
      (*n_face)[ipart] = _n_face;
      (*face_vtx_idx) [ipart] = (int *) malloc(sizeof(int) * (_n_face + 1));
      (*face_vtx_idx) [ipart][0] = 0;
      (*face_vtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*face_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_face);

      int _n_vtx = 0;
      (*n_vtx)[ipart] = _n_vtx;
      (*vtx_coord)[ipart] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
      (*vtx_ln_to_gn)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);
    }
  }

  PDM_MPI_Bcast (nGFace, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (nGVtx, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);

}


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_export_ini_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 const int      n_part,
 int            *n_face,
 int            **face_vtx_idx,
 int            **face_vtx,
 PDM_g_num_t    **face_ln_to_gn,
 int            *n_vtx,
 double         **vtx_coord,
 PDM_g_num_t    **vtx_ln_to_gn
)
{

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  PDM_writer_t *id_cs = PDM_writer_create ("Ensight",
                                           PDM_WRITER_FMT_ASCII,
                                           PDM_WRITER_TOPO_CST,
                                           PDM_WRITER_OFF,
                                           "test_2d_surf_ens",
                                           "mesh1",
                                           pdm_mpi_comm,
                                           PDM_IO_KIND_MPI_SIMPLE,
                                           1.,
                                           NULL);

  /*
   * Creation des variables
   */

  int id_var_num_part;
  int id_var_coo_x;
  int id_var_coo_xyz;
  int id_geom;

  id_var_num_part = PDM_writer_var_create (id_cs,
                                           PDM_WRITER_OFF,
                                           PDM_WRITER_VAR_SCALAR,
                                           PDM_WRITER_VAR_ELEMENTS,
                                           "num_part");

  id_var_coo_x = PDM_writer_var_create (id_cs,
                                        PDM_WRITER_ON,
                                        PDM_WRITER_VAR_SCALAR,
                                        PDM_WRITER_VAR_VERTICES,
                                        "coo_x");

  id_var_coo_xyz = PDM_writer_var_create (id_cs,
                                          PDM_WRITER_ON,
                                          PDM_WRITER_VAR_VECTOR,
                                          PDM_WRITER_VAR_VERTICES,
                                          "coo_xyz");

    /*
     * Creation de la geometrie
     */

    char nom_geom[6];
    strcpy (nom_geom,"mesh1");

    id_geom = PDM_writer_geom_create (id_cs,
                                             nom_geom,
                                             n_part);
    /*
     * Debut des ecritures
     */

    int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

    int *n_part_procs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                       (void *) n_part_procs, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

    free(n_part_procs);

    PDM_writer_step_beg (id_cs, 0.);

    int **_face_nb =  malloc(sizeof(int *) * n_part);
    int **_face_idx =  malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_writer_geom_coord_set (id_cs,
                                 id_geom,
                                 ipart,
                                 n_vtx[ipart],
                                 vtx_coord[ipart],
                                 vtx_ln_to_gn[ipart],
                                 PDM_OWNERSHIP_USER);

      _face_nb[ipart] = malloc(sizeof(int) * n_face[ipart]);
      _face_idx[ipart] = malloc(sizeof(int) * n_face[ipart]);

      for (int j = 0; j < n_face[ipart]; j++) {
        _face_nb[ipart][j]  = face_vtx_idx[ipart][j+1] - face_vtx_idx[ipart][j];
        _face_idx[ipart][j] = face_vtx_idx[ipart][j] + 1;
      }

      PDM_writer_geom_faces_facesom_add (id_cs,
                                         id_geom,
                                         ipart,
                                         n_face[ipart],
                                         _face_idx[ipart],
                                         _face_nb[ipart],
                                         face_vtx[ipart],
                                         face_ln_to_gn[ipart]);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (_face_nb[ipart]);
      free (_face_idx[ipart]);
    }

    free(_face_nb);
    free(_face_idx);

    PDM_writer_geom_write(id_cs,
                          id_geom);

    /* Creation des variables :
       - numero de partition
       - scalaire
       - vecteur
       - tenseur
    */

    PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      val_num_part[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_face[ipart]);
      val_coo_x[ipart]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_vtx[ipart]);
      val_coo_xyz[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * n_vtx[ipart]);
      nsom_part[ipart]    = n_vtx[ipart];

      for (int i = 0; i < n_face[ipart]; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[i_rank];
      }

      for (int i = 0; i < n_vtx[ipart]; i++) {
        val_coo_x[ipart][i]       = vtx_coord[ipart][3*i];
        val_coo_xyz[ipart][3*i  ] = vtx_coord[ipart][3*i  ];
        val_coo_xyz[ipart][3*i+1] = vtx_coord[ipart][3*i+1];
        val_coo_xyz[ipart][3*i+2] = vtx_coord[ipart][3*i+2];
      }

      PDM_writer_var_set (id_cs,
                          id_var_num_part,
                          id_geom,
                          ipart,
                          val_num_part[ipart]);

      PDM_writer_var_set (id_cs,
                          id_var_coo_x,
                          id_geom,
                          ipart,
                          val_coo_x[ipart]);

      PDM_writer_var_set (id_cs,
                          id_var_coo_xyz,
                          id_geom,
                          ipart,
                          val_coo_xyz[ipart]);

    }

    PDM_writer_var_write (id_cs,
                          id_var_num_part);

    PDM_writer_var_free (id_cs,
                         id_var_num_part);

    PDM_writer_var_write (id_cs,
                          id_var_coo_x);

    PDM_writer_var_free (id_cs,
                         id_var_coo_x);

    PDM_writer_var_write (id_cs,
                          id_var_coo_xyz);

    PDM_writer_var_free (id_cs,
                         id_var_coo_xyz);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (val_num_part[ipart]);
      free (val_coo_x[ipart]);
      free (val_coo_xyz[ipart]);
    }

    free (val_num_part);
    free (val_coo_x);
    free (val_coo_xyz);
    free (nsom_part);

    PDM_writer_step_end (id_cs);
    PDM_writer_geom_data_free (id_cs,
                               id_geom);

    PDM_writer_geom_free (id_cs,
                          id_geom);
    PDM_writer_free (id_cs);

    free (debPartProcs);

}


/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{
  PDM_MPI_Init (&argc, &argv);

  /*
   *  Set default values
   */

  PDM_g_num_t      n_vtx_segA = 4;
  double           lengthA  = 1.;
  double           xminA = 0.;
  double           yminA = 0.;
  int              n_partA   = 1;

  int              post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;
  int              haveRandom = 1;
  int              randomTimeInit = 0;

  int              randomMeshAInit = -1;

  int              nProcData = -1;

  int              i_rank;
  int              numProcs;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_segA,
              &lengthA,
              &xminA,
              &yminA,
              &n_partA,
              &post,
              (int *) &method,
              &haveRandom,
              &randomTimeInit,
              &randomMeshAInit,
              &nProcData
              );

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_vtx_segA : %d\n", n_vtx_segA);
    PDM_printf ("  - lengthA : %f\n", lengthA);
    PDM_printf ("  - xminA : %d\n", xminA);
    ;
    PDM_printf ("  - n_partA : %d\n", n_partA);
    PDM_printf ("  - post : %d\n", post);
    PDM_printf ("  - method : %d\n", method);
    PDM_printf ("  - haveRandom : %d\n", haveRandom);
    PDM_printf ("  - randomTimeInit : %d\n", randomTimeInit);
  }

  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t      nGFace;
  PDM_g_num_t      nGVtx;
  int            *n_face;
  int            **face_vtx_idx;
  int            **face_vtx;
  PDM_g_num_t    **face_ln_to_gn;
  int            *n_vtx;
  double         **vtx_coord;
  PDM_g_num_t    **vtx_ln_to_gn;

  int initRandom = 0;

  int n_part = n_partA;

  PDM_MPI_Comm meshComm = PDM_MPI_COMM_WORLD;
  int activeRankMesh = 1;

  if (nProcData > 0 && nProcData < numProcs) {
    int rankInNode = PDM_io_mpi_node_rank (PDM_MPI_COMM_WORLD);

    int nNode = 0;
    int iNode = -1;
    int masterRank = 0;
    if (rankInNode == 0) {
      masterRank = 1;
    }

    int *rankInNodes = malloc(sizeof(int) * numProcs);

    PDM_MPI_Allreduce (&masterRank, &nNode, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
    PDM_MPI_Allgather (&rankInNode, 1, PDM_MPI_INT, rankInNodes, 1, PDM_MPI_INT, PDM_MPI_COMM_WORLD);

    activeRankMesh = 0;

    for (int i = 0; i < i_rank; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && rankInNode == 0) {
        activeRankMesh = 1;
      }
    }

    else {

      if (rankInNode < (nProcData / nNode)) {
        activeRankMesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        activeRankMesh = 1;
      }

    }

    PDM_MPI_Comm_split(PDM_MPI_COMM_WORLD, activeRankMesh, i_rank, &meshComm);

    free (rankInNodes);
  }

  double xmin;
  double ymin;
  double length;
  PDM_g_num_t n_vtx_seg;

  n_vtx_seg = n_vtx_segA;
  length = lengthA;
  xmin = xminA;
  ymin = yminA;
  n_part = n_partA;

  if (randomTimeInit) {
    initRandom =  time( NULL );
  }

  if (randomMeshAInit != -1) {
    initRandom = randomMeshAInit;
  }

  _create_split_mesh (activeRankMesh,
                      meshComm,
                      xmin,
                      ymin,
                      n_vtx_seg,
                      length,
                      n_part,
                      method,
                      haveRandom,
                      initRandom,
                      &nGFace,
                      &nGVtx,
                      &n_face,
                      &face_vtx_idx,
                      &face_vtx,
                      &face_ln_to_gn,
                      &n_vtx,
                      &vtx_coord,
                      &vtx_ln_to_gn);

  if (activeRankMesh) {
    _export_ini_mesh (meshComm,
                      n_part,
                      n_face,
                      face_vtx_idx,
                      face_vtx,
                      face_ln_to_gn,
                      n_vtx,
                      vtx_coord,
                      vtx_ln_to_gn);
  }

  if (nProcData > 0 && nProcData < numProcs) {
    PDM_MPI_Comm_free(&meshComm);
  }

  if (1 == 0) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      printf("vtx_coord :\n");
      for (int j = 0; j < n_vtx[ipart]; j++) {
        printf(""PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", vtx_ln_to_gn[ipart][j],
               vtx_coord[ipart][3*j],
               vtx_coord[ipart][3*j+1],
               vtx_coord[ipart][3*j+2]);
      }
      printf("face_vtx :\n");
      for (int j = 0; j < n_face[ipart]; j++) {
        printf(""PDM_FMT_G_NUM" :", face_ln_to_gn[ipart][j]);
        for (int k = face_vtx_idx[ipart][j]; k < face_vtx_idx[ipart][j+1]; k++) {
          printf(" %d", face_vtx[ipart][k]);
        }
          printf("\n");
      }
    }
  }


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(face_vtx_idx[i_part]);
    free(face_vtx   [i_part]);
    free(face_ln_to_gn[i_part]);
    free(vtx_coord  [i_part]);
    free(vtx_ln_to_gn [i_part]);
  }

  free(n_face);
  free(n_vtx);
  free(face_vtx_idx);
  free(face_vtx   );
  free(face_ln_to_gn);
  free(vtx_coord  );
  free(vtx_ln_to_gn );

  PDM_MPI_Finalize ();

  return 0;

}
