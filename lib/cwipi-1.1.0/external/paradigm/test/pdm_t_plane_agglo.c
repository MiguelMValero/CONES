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
#include "pdm_part_coarse_mesh.h"

#include "pdm_writer.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"


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
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS for split.\n\n"
     "  -pt-scocth       Call PT-Scotch for split.\n\n"
     "  -agglo_metis     Call METIS for agglo.\n\n"
     "  -agglo_scocth    Call Scotch for agglo.\n\n"
     "  -cr     <level>  Coarse rate\n\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
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
 PDM_g_num_t  *n_vtx_seg,
 double        *length,
 int           *n_part,
 double         *cr,
 int           *post,
 int           *method,
 char         **method_agglo,
 int           *haveRandom
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol (argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp (argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *length = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_part = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-cr") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *cr = atof (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-no_random") == 0) {
      *haveRandom = 0;
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
    else if (strcmp (argv[i], "-agglo_scotch") == 0) {
      *method_agglo = (char *) malloc (sizeof (char) * (strlen("PDM_COARSE_MESH_SCOTCH") + 1));
      strcpy(*method_agglo, "PDM_COARSE_MESH_SCOTCH");
    }
    else if (strcmp (argv[i], "-agglo_metis") == 0) {
      *method_agglo = (char *) malloc (sizeof (char) * (strlen("PDM_COARSE_MESH_METIS") + 1));
      strcpy(*method_agglo, "PDM_COARSE_MESH_METIS");
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      n_vtx_seg  Number of arguments
 * \param [in]      length   Lenght of square
 * \param [in]      n_part    Number to obtain on this processus
 * \param [in]      post     mesh export status
 * \param [in]      method   Split method
 *
 * \return ppart_id  ppart identifier
 *
 */

static PDM_part_t *
_create_split_mesh
(
 int               imesh,
 PDM_MPI_Comm      pdm_mpi_comm,
 PDM_g_num_t       n_vtx_seg,
 double            length,
 int               n_part,
 PDM_part_split_t  method,
 int               haveRandom,
 PDM_g_num_t      *nGFace,
 PDM_g_num_t      *nGVtx,
 PDM_g_num_t      *nGEdge,
 int              *n_total_part,
 int              *nEdgeGroup
)
{
  struct timeval t_elaps_debut;

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

  double       xmin = 0.;
  double       xmax = length;
  double       ymin = 0.;
  double       ymax = length;
  PDM_g_num_t  nx   = n_vtx_seg;
  PDM_g_num_t  ny   = n_vtx_seg;
  int          dn_face;
  int          dn_vtx;
  int          dNEdge;
  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  double      *dvtx_coord;
  PDM_g_num_t *dFaceEdge;
  PDM_g_num_t *dEdgeVtx;
  PDM_g_num_t *dEdgeFace;
  int         *dEdgeGroupIdx;
  PDM_g_num_t *dEdgeGroup;

  int          initRandom = 0;

  /*
   *  Create mesh i
   */

  if (imesh == 1) {
    nx *= 2;
    ny *= 2;
  }

  ++initRandom;

  gettimeofday(&t_elaps_debut, NULL);

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
                     nGEdge,
                     &dn_vtx,
                     &dvtx_coord,
                     &dn_face,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dFaceEdge,
                     &dNEdge,
                     &dEdgeVtx,
                     &dEdgeFace,
                     nEdgeGroup,
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
    PDM_printf("[%d] Temps dans creeMaillagePolygone2D %d : %12.5e\n",
           i_rank, imesh, t_elapsed);

  if (0 == 1) {

    PDM_printf ("edgegroup : ");
    for (int i = 0; i < *nEdgeGroup; i++) {
      for (int j = dEdgeGroupIdx[i]; j <  dEdgeGroupIdx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dEdgeGroup[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dface_vtx : ");
    for (int i = 0; i < dn_face; i++) {
      for (int j = dface_vtx_idx[i]; j <  dface_vtx_idx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dface_vtx[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dfaceedge : ");
    for (int i = 0; i < dn_face; i++) {
      for (int j = dface_vtx_idx[i]; j <  dface_vtx_idx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dFaceEdge[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedgevtx : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i+1]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedgeface : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dEdgeVtx[2*i+1]);
      PDM_printf ("\n");
    }
  }

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc (dn_face*sizeof(int));
  int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

  dEdgeVtxIdx[0] = 0;
  for (int i = 0; i < dNEdge; i++) {
    dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
  }

  /*
   *  Split mesh i
   */

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
                                       method,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       "PDM_PART_RENUM_FACE_NONE",
                                       n_property_cell,
                                       renum_properties_cell,
                                       n_property_face,
                                       renum_properties_face,
                                       n_part,
                                       dn_face,
                                       dNEdge,
                                       dn_vtx,
                                       *nEdgeGroup,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dEdgeFace,
                                       dEdgeVtxIdx,
                                       dEdgeVtx,
                                       NULL,
                                       dvtx_coord,
                                       NULL,
                                       dEdgeGroupIdx,
                                       dEdgeGroup);

  free (dcell_part);

  double *elapsed  = NULL;
  double *cpu      = NULL;
  double *cpu_user = NULL;
  double *cpu_sys  = NULL;

  PDM_part_time_get (ppart,
                     &elapsed,
                     &cpu,
                     &cpu_user,
                     &cpu_sys);

  if (i_rank == 0)
    PDM_printf("[%d] Temps dans ppart %d : %12.5e\n",
           i_rank, imesh, elapsed[0]);

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

  /* if (i_rank == 0) { */
  /*   PDM_printf ("Statistics :\n"); */
  /*   PDM_printf ("  - Number of cells :\n"); */
  /*   PDM_printf ("       * average            : %i\n", cells_average);    */
  /*   PDM_printf ("       * median             : %i\n", cells_median);    */
  /*   PDM_printf ("       * standard deviation : %12.5e\n", cells_std_deviation);    */
  /*   PDM_printf ("       * min                : %i\n", cells_min);    */
  /*   PDM_printf ("       * max                : %i\n", cells_max);    */
  /*   PDM_printf ("  - Number of faces exchanging with another partition :\n"); */
  /*   PDM_printf ("       * average            : %i\n", bound_part_faces_average);    */
  /*   PDM_printf ("       * median             : %i\n", bound_part_faces_median);    */
  /*   PDM_printf ("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);    */
  /*   PDM_printf ("       * min                : %i\n", bound_part_faces_min);    */
  /*   PDM_printf ("       * max                : %i\n", bound_part_faces_max);    */
  /*   PDM_printf ("       * total              : %i\n", bound_part_faces_sum);    */
  /* } */

  free (dvtx_coord);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dFaceEdge);
  free (dEdgeVtxIdx);
  free (dEdgeVtx);
  free (dEdgeFace);
  free (dEdgeGroupIdx);
  free (dEdgeGroup);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_face;
    int nEdge;
    int nEdgePartBound;
    int n_vtx;
    int n_proc;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_face,
                           &nEdge,
                           &nEdgePartBound,
                           &n_vtx,
                           &n_proc,
                           n_total_part,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);

  }

  return ppart;
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
 const PDM_MPI_Comm  pdm_mpi_comm,
 PDM_part_t         *ppart,
 const int           n_part
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
                                           "pdm_t_plane_agglo_ens",
                                           "fine_mesh",
                                           pdm_mpi_comm,
                                           PDM_IO_KIND_MPI_SIMPLE,
                                           1.,
                                           NULL );

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

  char nom_geom[] = "mesh1";

  id_geom = PDM_writer_geom_create (id_cs,
                                    nom_geom,
                                    n_part);
  /*
   * Debut des ecritures
   */

  int       **edgeVtxIdx1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
  int       **edgeVtxNB1   = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
  int       **faceEdgeIdx1 = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
  int       **faceEdgeNB1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);

  int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

  int *n_partProcs = (int *) malloc(sizeof(int) * numProcs);

  PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                 (void *) n_partProcs, 1, PDM_MPI_INT,
                 PDM_MPI_COMM_WORLD);

  int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

  debPartProcs[0] = 0;
  for (int i = 0; i < numProcs; i++) {
    debPartProcs[i+1] = debPartProcs[i] + n_partProcs[i];
  }

  free(n_partProcs);

  PDM_writer_step_beg (id_cs, 0.);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_face;
    int nEdge;
    int nEdgePartBound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_face,
                           &nEdge,
                           &nEdgePartBound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);

    int          *face_tag;
    int          *faceEdgeIdx;
    int          *faceEdge;
    PDM_g_num_t *face_ln_to_gn;
    int          *edgeTag;
    int          *edgeFace;
    int          *edgeVtxIdx;
    int          *edgeVtx;
    PDM_g_num_t *edgeLNToGN;
    int          *edgePartBoundProcIdx;
    int          *edgePartBoundPartIdx;
    int          *edgePartBound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *edgeGroupIdx;
    int          *edgeGroup;
    PDM_g_num_t *edgeGroupLNToGN;

    assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

    PDM_part_part_val_get (ppart,
                           i_part,
                           &face_tag,
                           &faceEdgeIdx,
                           &faceEdge,
                           &face_ln_to_gn,
                           &edgeTag,
                           &edgeFace,
                           &edgeVtxIdx,
                           &edgeVtx,
                           &edgeLNToGN,
                           &edgePartBoundProcIdx,
                           &edgePartBoundPartIdx,
                           &edgePartBound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &edgeGroupIdx,
                           &edgeGroup,
                           &edgeGroupLNToGN);

    edgeVtxIdx1[i_part]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
    edgeVtxNB1[i_part]   = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
    faceEdgeIdx1[i_part] = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * n_face);
    faceEdgeNB1[i_part]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * n_face);

    for (int i = 0; i < n_face; i++) {
      faceEdgeNB1[i_part][i] = faceEdgeIdx[i+1] - faceEdgeIdx[i];
      faceEdgeIdx1[i_part][i] = faceEdgeIdx[i] + 1;
    }

    for (int i = 0; i < nEdge; i++) {
      edgeVtxNB1[i_part][i] = edgeVtxIdx[i+1] - edgeVtxIdx[i];
      edgeVtxIdx1[i_part][i] = edgeVtxIdx[i] + 1;
    }

    PDM_writer_geom_coord_set (id_cs,
                               id_geom,
                               i_part,
                               n_vtx,
                               vtx,
                               vtx_ln_to_gn,
                               PDM_OWNERSHIP_USER);

    PDM_writer_geom_cell2d_cellface_add (id_cs,
                                         id_geom,
                                         i_part,
                                         n_face,
                                         nEdge,
                                         edgeVtxIdx1[i_part],
                                         edgeVtxNB1[i_part],
                                         edgeVtx,
                                         faceEdgeIdx1[i_part],
                                         faceEdgeNB1[i_part],
                                         faceEdge,
                                         face_ln_to_gn);
  }

  PDM_writer_geom_write(id_cs,
                        id_geom);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (edgeVtxIdx1[i_part]);
    free (edgeVtxNB1[i_part]);
    free (faceEdgeIdx1[i_part]);
    free (faceEdgeNB1[i_part]);
  }

  free (edgeVtxIdx1);
  free (edgeVtxNB1);
  free (faceEdgeIdx1);
  free (faceEdgeNB1);

  /* Creation des variables :
     - numero de partition
     - scalaire
     - vecteur
     - tenseur
  */

  PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
  PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_face;
    int nEdge;
    int nEdgePartBound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_face,
                           &nEdge,
                           &nEdgePartBound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);

    int          *face_tag;
    int          *faceEdgeIdx;
    int          *faceEdge;
    PDM_g_num_t *face_ln_to_gn;
    int          *edgeTag;
    int          *edgeFace;
    int          *edgeVtxIdx;
    int          *edgeVtx;
    PDM_g_num_t *edgeLNToGN;
    int          *edgePartBoundProcIdx;
    int          *edgePartBoundPartIdx;
    int          *edgePartBound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *edgeGroupIdx;
    int          *edgeGroup;
    PDM_g_num_t *edgeGroupLNToGN;

    assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

    PDM_part_part_val_get (ppart,
                           i_part,
                           &face_tag,
                           &faceEdgeIdx,
                           &faceEdge,
                           &face_ln_to_gn,
                           &edgeTag,
                           &edgeFace,
                           &edgeVtxIdx,
                           &edgeVtx,
                           &edgeLNToGN,
                           &edgePartBoundProcIdx,
                           &edgePartBoundPartIdx,
                           &edgePartBound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &edgeGroupIdx,
                           &edgeGroup,
                           &edgeGroupLNToGN);

    val_num_part[i_part] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_face);
    val_coo_x[i_part]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_vtx);
    val_coo_xyz[i_part]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * n_vtx);
    nsom_part[i_part]    = n_vtx;

    for (int i = 0; i < n_face; i++) {
      val_num_part[i_part][i] = i_part + 1 + debPartProcs[i_rank];
    }

    for (int i = 0; i < n_vtx; i++) {
      val_coo_x[i_part][i]       = vtx[3*i];
      val_coo_xyz[i_part][3*i  ] = vtx[3*i  ];
      val_coo_xyz[i_part][3*i+1] = vtx[3*i+1];
      val_coo_xyz[i_part][3*i+2] = vtx[3*i+2];
    }

    PDM_writer_var_set (id_cs,
                        id_var_num_part,
                        id_geom,
                        i_part,
                        val_num_part[i_part]);

    PDM_writer_var_set (id_cs,
                        id_var_coo_x,
                        id_geom,
                        i_part,
                        val_coo_x[i_part]);

    PDM_writer_var_set (id_cs,
                        id_var_coo_xyz,
                        id_geom,
                        i_part,
                        val_coo_xyz[i_part]);

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

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (val_num_part[i_part]);
    free (val_coo_x[i_part]);
    free (val_coo_xyz[i_part]);
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
 * \brief  Create and split Mesh
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

// static void
// _export_coarse_mesh
// (
//  const PDM_MPI_Comm pdm_mpi_comm,
//  int            cmId,
//  const int      n_part
// )
// {

//   int i_rank;
//   int numProcs;

//   PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
//   PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

//   /*
//    *  Export Mesh to Ensight
//    */

//   int id_cs;


//   id_cs = PDM_writer_create ("Ensight",
//                               PDM_WRITER_FMT_ASCII,
//                               PDM_WRITER_TOPO_CST,
//                               PDM_WRITER_OFF,
//                               "pdm_t_plane_agglo_ens",
//                               "coarse_mesh",
//                               pdm_mpi_comm,
//                               PDM_IO_KIND_MPI_SIMPLE,
//                               1.,
//                               NULL);

//   /*
//    * Creation des variables
//    */

//   int id_var_num_part;
//   int id_var_coo_x;
//   int id_var_coo_xyz;
//   int id_geom;

//   id_var_num_part = PDM_writer_var_create (id_cs,
//                                           PDM_WRITER_OFF,
//                                           PDM_WRITER_VAR_SCALAR,
//                                           PDM_WRITER_VAR_ELEMENTS,
//                                           "num_part");

//   id_var_coo_x = PDM_writer_var_create (id_cs,
//                                        PDM_WRITER_ON,
//                                        PDM_WRITER_VAR_SCALAR,
//                                        PDM_WRITER_VAR_VERTICES,
//                                        "coo_x");

//   id_var_coo_xyz = PDM_writer_var_create (id_cs,
//                                          PDM_WRITER_ON,
//                                          PDM_WRITER_VAR_VECTOR,
//                                          PDM_WRITER_VAR_VERTICES,
//                                          "coo_xyz");

//     /*
//      * Creation de la geometrie
//      */

//   char nom_geom[] = "mesh1";

//   id_geom = PDM_writer_geom_create (id_cs,
//                                  nom_geom,
//                                  n_part);
//   /*
//    * Debut des ecritures
//    */

//   int       **edgeVtxIdx1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
//   int       **edgeVtxNB1   = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
//   int       **faceEdgeIdx1 = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);
//   int       **faceEdgeNB1  = (PDM_l_num_t **) malloc (sizeof(PDM_l_num_t *) * n_part);

//   int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

//   int *n_partProcs = (int *) malloc(sizeof(int) * numProcs);

//   PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
//                  (void *) n_partProcs, 1, PDM_MPI_INT,
//                  PDM_MPI_COMM_WORLD);

//   int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

//   debPartProcs[0] = 0;
//   for (int i = 0; i < numProcs; i++) {
//     debPartProcs[i+1] = debPartProcs[i] + n_partProcs[i];
//   }

//   free(n_partProcs);

//   PDM_writer_step_beg (id_cs, 0.);

//   for (int i_part = 0; i_part < n_part; i_part++) {

//     int n_face;
//     int nEdge;
//     int nEdgeGroup;

//     int nEdgePartBound;
//     int n_vtx;
//     int n_proc;
//     int n_total_part;
//     int sFaceEdge;
//     int sEdgeVtx;
//     int sEdgeGroup;
//     int scoarse_face_to_fine_face;


//     PDM_part_coarse_mesh_part_dim_get(cmId,
//                                       i_part,
//                                       &n_face,
//                                       &nEdge,
//                                       &nEdgePartBound,
//                                       &n_vtx,
//                                       &n_proc,
//                                       &n_total_part,
//                                       &nEdgeGroup,
//                                       &sFaceEdge,
//                                       &sEdgeVtx,
//                                       &sEdgeGroup,
//                                       &scoarse_face_to_fine_face);

//     int          *face_tag;
//     int          *faceEdgeIdx;
//     int          *faceEdge;
//     PDM_g_num_t *face_ln_to_gn;
//     int          *edgeTag;
//     int          *edgeFace;
//     int          *edgeVtxIdx;
//     int          *edgeVtx;
//     PDM_g_num_t *edgeLNToGN;
//     int          *edgePartBoundProcIdx;
//     int          *edgePartBoundPartIdx;
//     int          *edgePartBound;
//     int          *vtx_tag;
//     double       *vtx;
//     PDM_g_num_t *vtx_ln_to_gn;
//     int          *edgeGroupIdx;
//     int          *edgeGroup;
//     PDM_g_num_t *edgeGroupLNToGN;
//     int         *faceInitFaceIdx;
//     int         *faceInitFace;
//     int         *edgeInitEdge;
//     int         *vtxInitVtx;
//     int         *edgeGroupInitEdgeGroup;

//     assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

//     PDM_part_coarse_mesh_part_get (cmId,
//                                    i_part,
//                                   &faceEdgeIdx,
//                                   &faceEdge,
//                                   &face_tag,
//                                   &face_ln_to_gn,
//                                   &faceInitFaceIdx,
//                                   &faceInitFace,
//                                   &edgeFace,
//                                   &edgeVtxIdx,
//                                   &edgeVtx,
//                                   &edgeTag,
//                                   &edgeLNToGN,
//                                   &edgeGroupInitEdgeGroup,
//                                   &edgeInitEdge,
//                                   &vtx,
//                                   &vtx_tag,
//                                   &vtx_ln_to_gn,
//                                   &vtxInitVtx,
//                                   &edgeGroupIdx,
//                                   &edgeGroup,
//                                   &edgeGroupLNToGN,
//                                   &edgePartBoundProcIdx,
//                                   &edgePartBoundPartIdx,
//                                   &edgePartBound);

//     edgeVtxIdx1[i_part]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
//     edgeVtxNB1[i_part]   = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * nEdge);
//     faceEdgeIdx1[i_part] = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * n_face);
//     faceEdgeNB1[i_part]  = (PDM_l_num_t *) malloc (sizeof(PDM_l_num_t) * n_face);

//     for (int i = 0; i < n_face; i++) {
//       faceEdgeNB1[i_part][i] = faceEdgeIdx[i+1] - faceEdgeIdx[i];
//       faceEdgeIdx1[i_part][i] = faceEdgeIdx[i] + 1;
//     }

//     for (int i = 0; i < nEdge; i++) {
//       edgeVtxNB1[i_part][i] = edgeVtxIdx[i+1] - edgeVtxIdx[i];
//       edgeVtxIdx1[i_part][i] = edgeVtxIdx[i] + 1;
//     }

//     PDM_writer_geom_coord_set (id_cs,
//                        id_geom,
//                        i_part,
//                        n_vtx,
//                        vtx,
//                        vtx_ln_to_gn);

//     PDM_writer_geom_cell2d_cellface_add (id_cs,
//                                          id_geom,
//                                          i_part,
//                                          n_face,
//                                          nEdge,
//                                          edgeVtxIdx1[i_part],
//                                          edgeVtxNB1[i_part],
//                                          edgeVtx,
//                                          faceEdgeIdx1[i_part],
//                                          faceEdgeNB1[i_part],
//                                          faceEdge,
//                                          face_ln_to_gn);
//   }

//   PDM_writer_geom_write(id_cs,
//               id_geom);

//   for (int i_part = 0; i_part < n_part; i_part++) {
//     free (edgeVtxIdx1[i_part]);
//     free (edgeVtxNB1[i_part]);
//     free (faceEdgeIdx1[i_part]);
//     free (faceEdgeNB1[i_part]);
//   }

//   free (edgeVtxIdx1);
//   free (edgeVtxNB1);
//   free (faceEdgeIdx1);
//   free (faceEdgeNB1);

//   /* Creation des variables :
//      - numero de partition
//      - scalaire
//      - vecteur
//      - tenseur
//   */

//   PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
//   PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
//   PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

//   for (int i_part = 0; i_part < n_part; i_part++) {

//     int n_face;
//     int nEdge;
//     int nEdgePartBound;
//     int n_vtx;
//     int n_proc;
//     int n_total_part;
//     int sFaceEdge;
//     int sEdgeVtx;
//     int sEdgeGroup;
//     int nEdgeGroup2;

//     PDM_part_part_dim_get (cmId,
//                         i_part,
//                         &n_face,
//                         &nEdge,
//                         &nEdgePartBound,
//                         &n_vtx,
//                         &n_proc,
//                         &n_total_part,
//                         &sFaceEdge,
//                         &sEdgeVtx,
//                         &sEdgeGroup,
//                         &nEdgeGroup2);

//     int          *face_tag;
//     int          *faceEdgeIdx;
//     int          *faceEdge;
//     PDM_g_num_t *face_ln_to_gn;
//     int          *edgeTag;
//     int          *edgeFace;
//     int          *edgeVtxIdx;
//     int          *edgeVtx;
//     PDM_g_num_t *edgeLNToGN;
//     int          *edgePartBoundProcIdx;
//     int          *edgePartBoundPartIdx;
//     int          *edgePartBound;
//     int          *vtx_tag;
//     double       *vtx;
//     PDM_g_num_t *vtx_ln_to_gn;
//     int          *edgeGroupIdx;
//     int          *edgeGroup;
//     PDM_g_num_t *edgeGroupLNToGN;

//     assert (sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

//     PDM_part_part_val_get (cmId,
//                         i_part,
//                         &face_tag,
//                         &faceEdgeIdx,
//                         &faceEdge,
//                         &face_ln_to_gn,
//                         &edgeTag,
//                         &edgeFace,
//                         &edgeVtxIdx,
//                         &edgeVtx,
//                         &edgeLNToGN,
//                         &edgePartBoundProcIdx,
//                         &edgePartBoundPartIdx,
//                         &edgePartBound,
//                         &vtx_tag,
//                         &vtx,
//                         &vtx_ln_to_gn,
//                         &edgeGroupIdx,
//                         &edgeGroup,
//                         &edgeGroupLNToGN);

//     val_num_part[i_part] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_face);
//     val_coo_x[i_part]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * n_vtx);
//     val_coo_xyz[i_part]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * n_vtx);
//     nsom_part[i_part]    = n_vtx;

//     for (int i = 0; i < n_face; i++) {
//       val_num_part[i_part][i] = i_part + 1 + debPartProcs[i_rank];
//     }

//     for (int i = 0; i < n_vtx; i++) {
//       val_coo_x[i_part][i]       = vtx[3*i];
//       val_coo_xyz[i_part][3*i  ] = vtx[3*i  ];
//       val_coo_xyz[i_part][3*i+1] = vtx[3*i+1];
//       val_coo_xyz[i_part][3*i+2] = vtx[3*i+2];
//     }

//     PDM_writer_var_set (id_cs,
//                 id_var_num_part,
//                 id_geom,
//                 i_part,
//                 val_num_part[i_part]);

//     PDM_writer_var_set (id_cs,
//                 id_var_coo_x,
//                 id_geom,
//                 i_part,
//                 val_coo_x[i_part]);

//     PDM_writer_var_set (id_cs,
//                 id_var_coo_xyz,
//                 id_geom,
//                 i_part,
//                 val_coo_xyz[i_part]);

//   }

//   PDM_writer_var_write (id_cs,
//               id_var_num_part);

//   PDM_writer_var_free (id_cs,
//               id_var_num_part);

//   PDM_writer_var_write (id_cs,
//               id_var_coo_x);

//   PDM_writer_var_free (id_cs,
//               id_var_coo_x);

//   PDM_writer_var_write (id_cs,
//               id_var_coo_xyz);

//   PDM_writer_var_free (id_cs,
//               id_var_coo_xyz);

//   for (int i_part = 0; i_part < n_part; i_part++) {
//     free (val_num_part[i_part]);
//     free (val_coo_x[i_part]);
//     free (val_coo_xyz[i_part]);
//   }

//   free (val_num_part);
//   free (val_coo_x);
//   free (val_coo_xyz);
//   free (nsom_part);

//   PDM_writer_step_end (id_cs);
//   PDM_writer_geom_data_free (id_cs,
//                     id_geom);

//   PDM_writer_geom_free (id_cs,
//                id_geom);
//   PDM_writer_free (id_cs);

//   free (debPartProcs);


// }


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

  PDM_g_num_t   n_vtx_seg = 4;
  double        length  = 1.;
  int           n_part   = 1;
  double         cr   = 0.5;
  int           post    = 0;
#ifdef PDM_HAVE_PTSCOTCH
 PDM_part_split_t  method  = PDM_PART_SPLIT_PTSCOTCH;
  const char *agglo_method = "PDM_COARSE_MESH_SCOTCH";
#else
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  =  PDM_PART_SPLIT_PARMETIS;
  const char *agglo_method = "PDM_COARSE_MESH_METIS";
#endif
#endif

  int           haveRandom = 0;
  int           i_rank;
  int           numProcs;

  /*
   *  Read args
   */

  char *_agglo_method = NULL;
  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &n_part,
              &cr,
              &post,
              (int *) &method,
              &_agglo_method,
              &haveRandom);

  if (_agglo_method != NULL) {
    agglo_method = _agglo_method;
  }

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t nGFace;
  PDM_g_num_t nGVtx;
  PDM_g_num_t nGEdge;
  int imesh = 0;
  int n_total_part;
  int nEdgeGroup;

  PDM_part_t *ppart = _create_split_mesh (imesh,
                                          PDM_MPI_COMM_WORLD,
                                          n_vtx_seg,
                                          length,
                                          n_part,
                                          method,
                                          haveRandom,
                                          &nGFace,
                                          &nGVtx,
                                          &nGEdge,
                                          &n_total_part,
                                          &nEdgeGroup);

  _export_ini_mesh (PDM_MPI_COMM_WORLD,
                    ppart,
                    n_part);

  /*
   *  Appel des fonctions d'agglo
   */

  // int cmId;
  const int  have_cell_tag = 0;
  const int  have_face_tag = 0;
  const int  have_vtx_tag = 0;
  const int  have_cell_weight = 0;
  const int  have_face_weight = 0;
  const int  have_face_group = 0;

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_coarse_mesh_t *cm = PDM_part_coarse_mesh_create (PDM_MPI_COMM_WORLD,
                                                       agglo_method,
                                                       "PDM_PART_RENUM_CELL_NONE",
                                                       "PDM_PART_RENUM_FACE_NONE",
                                                       n_property_cell,
                                                       renum_properties_cell,
                                                       n_property_face,
                                                       renum_properties_face,
                                                       n_part,
                                                       n_total_part,
                                                       nEdgeGroup,
                                                       have_cell_tag,
                                                       have_face_tag,
                                                       have_vtx_tag,
                                                       have_cell_weight,
                                                       have_face_weight,
                                                       have_face_group);

  if (_agglo_method != NULL) {
    free (_agglo_method);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_face;
    int nEdge;
    int nEdgePartBound;
    int n_vtx;
    int n_proc;
    int n_total_part1;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_face,
                           &nEdge,
                           &nEdgePartBound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part1,
                           &sFaceEdge,
                           &sEdgeVtx,
                           &sEdgeGroup,
                           &nEdgeGroup2);

    int         *face_tag;
    int         *faceEdgeIdx;
    int         *faceEdge;
    PDM_g_num_t *face_ln_to_gn;
    int         *edgeTag;
    int         *edgeFace;
    int         *edgeVtxIdx;
    int         *edgeVtx;
    PDM_g_num_t *edgeLNToGN;
    int         *edgePartBoundProcIdx;
    int         *edgePartBoundPartIdx;
    int         *edgePartBound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *edgeGroupIdx;
    int         *edgeGroup;
    PDM_g_num_t *edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &face_tag,
                           &faceEdgeIdx,
                           &faceEdge,
                           &face_ln_to_gn,
                           &edgeTag,
                           &edgeFace,
                           &edgeVtxIdx,
                           &edgeVtx,
                           &edgeLNToGN,
                           &edgePartBoundProcIdx,
                           &edgePartBoundPartIdx,
                           &edgePartBound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &edgeGroupIdx,
                           &edgeGroup,
                           &edgeGroupLNToGN);

    int _nc = (int) ((1. - cr) * n_face);
    int nc = PDM_MAX (_nc, 1);

    int _nEdgeGroup = 0;
    PDM_part_coarse_mesh_input (cm,
                                i_part,
                                nc,
                                n_face,
                                nEdge,
                                n_vtx,
                                _nEdgeGroup,
                                nEdgePartBound,
                                faceEdgeIdx,
                                faceEdge,
                                face_tag,
                                NULL,
                                NULL,
                                face_ln_to_gn,
                                edgeFace,
                                edgeVtxIdx,
                                edgeVtx,
                                edgeTag,
                                edgeLNToGN,
                                vtx,
                                vtx_tag,
                                vtx_ln_to_gn,
                                NULL,
                                NULL,
                                NULL,
                                edgePartBoundProcIdx,
                                edgePartBoundPartIdx,
                                edgePartBound);

  }

  PDM_part_coarse_mesh_compute (cm);

  // _export_coarse_mesh (PDM_MPI_COMM_WORLD,
  //                      cmId,
  //                      n_part);

/*
 *  Free meshes
 */
  PDM_part_free(ppart);
  PDM_part_coarse_mesh_free(cm);

 /* for (int imesh = 0; imesh < 2; imesh++) { */
 /*   for (int i_part = 0; i_part < n_part; i_part++) { */
 /*     free (face_vtx[imesh][i_part]); */
 /*   } */
 /*   free (face_vtx[imesh]); */

 /*   PDM_part_free (ppart_id[imesh]); */
 /* } */

 /* free (face_vtx); */
 if (i_rank == 0) {
   PDM_printf ("-- End\n");
 }


 PDM_MPI_Finalize ();


 return 0;

}
