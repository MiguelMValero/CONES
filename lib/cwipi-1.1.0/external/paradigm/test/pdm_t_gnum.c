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
#include "pdm_gnum.h"

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
 double        *length
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
 int               have_random,
 PDM_g_num_t      *n_g_face,
 PDM_g_num_t      *n_g_vtx,
 PDM_g_num_t      *n_g_edge,
 int              *n_total_part,
 int              *n_edge_group
)
{
  struct timeval t_elaps_debut;

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

  double        xmin = 0.;
  double        xmax = length;
  double        ymin = 0.;
  double        ymax = length;
  PDM_g_num_t  nx = n_vtx_seg;
  PDM_g_num_t  ny = n_vtx_seg;
  int           dn_face;
  int           dn_vtx;
  int           dNEdge;
  int          *dface_vtx_idx;
  PDM_g_num_t  *dface_vtx;
  double       *dvtx_coord;
  PDM_g_num_t  *dface_edge;
  PDM_g_num_t  *dedge_vtx;
  PDM_g_num_t  *dedge_face;
  int          *dedge_group_idx;
  PDM_g_num_t  *dedge_group;

  int           init_random = 0;

  /*
   *  Create mesh i
   */

  if (imesh == 1) {
    nx *= 2;
    ny *= 2;
  }

  ++init_random;

  gettimeofday(&t_elaps_debut, NULL);

  PDM_poly_surf_gen (pdm_mpi_comm,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     have_random,
                     init_random,
                     nx,
                     ny,
                     n_g_face,
                     n_g_vtx,
                     n_g_edge,
                     &dn_vtx,
                     &dvtx_coord,
                     &dn_face,
                     &dface_vtx_idx,
                     &dface_vtx,
                     &dface_edge,
                     &dNEdge,
                     &dedge_vtx,
                     &dedge_face,
                     n_edge_group,
                     &dedge_group_idx,
                     &dedge_group);

  // validation
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm, PDM_OWNERSHIP_KEEP);
  // fin validation

  PDM_gnum_set_from_coords (gen_gnum, 0, dn_vtx, dvtx_coord, NULL);

  PDM_gnum_compute (gen_gnum);

  const PDM_g_num_t *_numabs = PDM_gnum_get (gen_gnum, 0);

  for (int j = 0; j < dn_vtx; j++) {
    PDM_printf (PDM_FMT_G_NUM" %12.5e %12.5e %12.5e\n", _numabs[j], dvtx_coord[3*j],
                                                     dvtx_coord[3*j+1],
                                                     dvtx_coord[3*j+2]);
  }

  PDM_gnum_free (gen_gnum);

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
    for (int i = 0; i < *n_edge_group; i++) {
      for (int j = dedge_group_idx[i]; j <  dedge_group_idx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dedge_group[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dface_vtx : ");
    for (int i = 0; i < dn_face; i++) {
      for (int j = dface_vtx_idx[i]; j <  dface_vtx_idx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dface_vtx[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dface_edge : ");
    for (int i = 0; i < dn_face; i++) {
      for (int j = dface_vtx_idx[i]; j <  dface_vtx_idx[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, dface_edge[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedge_vtx : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dedge_vtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dedge_vtx[2*i+1]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedge_face : ");
    for (int i = 0; i < dNEdge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dedge_face[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dedge_face[2*i+1]);
      PDM_printf ("\n");
    }
  }

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc (dn_face*sizeof(int));
  int *dedge_vtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

  dedge_vtxIdx[0] = 0;
  for (int i = 0; i < dNEdge; i++) {
    dedge_vtxIdx[i+1] = 2 + dedge_vtxIdx[i];
  }

  /*
   *  Split mesh i
   */

  // int ppart_id;

  int n_property_cell = 0;
  int *renum_properties_cell = NULL;
  int n_property_face = 0;
  int *renum_properties_face = NULL;

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
                                       *n_edge_group,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dedge_face,
                                       dedge_vtxIdx,
                                       dedge_vtx,
                                       NULL,
                                       dvtx_coord,
                                       NULL,
                                       dedge_group_idx,
                                       dedge_group);

  free (dcell_part);

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
  free (dface_edge);
  free (dedge_vtxIdx);
  free (dedge_vtx);
  free (dedge_face);
  free (dedge_group_idx);
  free (dedge_group);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_face;
    int nEdge;
    int nEdgePartBound;
    int n_vtx;
    int n_proc;
    int sFaceEdge;
    int sEdgeVtx;
    int sEdgeGroup;
    int n_edge_group2;

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
                           &n_edge_group2);

  }

  return ppart;
}

/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[n] : 30,60
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
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;
  int           have_random = 0;

  int           i_rank;
  int           numProcs;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t n_g_face;
  PDM_g_num_t n_g_vtx;
  PDM_g_num_t n_g_edge;
  int imesh = 0;
  int n_total_part;
  int n_edge_group;

  PDM_part_t* ppart = _create_split_mesh (imesh,
                                          PDM_MPI_COMM_WORLD,
                                          n_vtx_seg,
                                          length,
                                          n_part,
                                          method,
                                          have_random,
                                          &n_g_face,
                                          &n_g_vtx,
                                          &n_g_edge,
                                          &n_total_part,
                                          &n_edge_group);


  // PDM_memory_stats_t* ms = PDM_memory_stats_create(3, comm);
  // PDM_memory_stats_add(ms, 0, "Start  : ");

  // int size = 2000000000;
  // double* test = malloc(size * sizeof(double));
  // for(int i = 0; i < size; ++i) {
  //   test[i] = 1.;
  // }
  // PDM_memory_stats_add(ms, 1, "Step 1 : ");


  // free(test);
  // PDM_memory_stats_add(ms, 2, "End    : ");

  // PDM_memory_stats_log(ms);
  // PDM_memory_stats_free(ms);

  PDM_part_free(ppart);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;

}


