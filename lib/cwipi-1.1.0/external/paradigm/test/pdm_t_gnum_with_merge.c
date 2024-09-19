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

#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_timer.h"


/*============================================================================
 * Type definition
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
 int          *n_part,
 double       *length
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

    else if (strcmp (argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        int _n_part = atoi (argv[i]);
        *n_part = (PDM_g_num_t) _n_part;
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
  int num_procs;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &num_procs);

  double        xmin = 0.;
  double        xmax = length;
  double        ymin = 0.;
  double        ymax = length;
  PDM_g_num_t   nx = n_vtx_seg;
  PDM_g_num_t   ny = n_vtx_seg;
  int           dn_face;
  int           dn_vtx;
  int           dn_edge;
  int          *dface_vtx_idx;
  PDM_g_num_t  *dface_vtx;
  double       *dvtx_coord;
  PDM_g_num_t  *dface_edge;
  PDM_g_num_t  *dedge_vtx;
  PDM_g_num_t  *dedge_face;
  int          *dedge_group_idx;
  PDM_g_num_t  *dedge_group;

  int           initRandom = 0;

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
                     have_random,
                     initRandom,
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
                     &dn_edge,
                     &dedge_vtx,
                     &dedge_face,
                     n_edge_group,
                     &dedge_group_idx,
                     &dedge_group);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_coords (gen_gnum, 0, dn_vtx, dvtx_coord, NULL);

  PDM_gnum_compute (gen_gnum);

  const PDM_g_num_t *_numabs2 = PDM_gnum_get (gen_gnum, 0);

  /* const PDM_g_num_t *_numabs2 = NULL; */
  /* int id; */

  /* int nn = 10000; */

  /* for (int i = 0; i < nn; i++) { */

  /*   if (i < nn - 1) { */

  /*     id = PDM_gnum_create (3, 1, PDM_TRUE, 1e-3, pdm_mpi_comm); */
  /*   } */

  /*   else { */
  /*     id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, pdm_mpi_comm); */

  /*   } */

  /*   double *dd = malloc (sizeof(double) * dn_vtx); */

  /*   for (int j = 0; j < dn_vtx; j++) { */
  /*     dd[j] = 1e-5; */
  /*   } */

  /*   if (i < nn - 1) { */
  /*     PDM_gnum_set_from_coords (id, 0, dn_vtx, dvtx_coord, dd); */
  /*   } */
  /*   else { */
  /*     PDM_gnum_set_from_coords (id, 0, dn_vtx, dvtx_coord, NULL); */
  /*   } */

  /*   PDM_gnum_compute (id); */

  /*   _numabs2 = PDM_gnum_get (id, 0); */

  /*   free(dd); */

  /*   if (i < nn - 1) { */
  /*     PDM_gnum_free (id, 0); */
  /*   } */

  /*   FILE *f = fopen("/proc/self/statm", "r"); */

  /*   long int mvirt, mres, mshared, val1, val2, val3; */
  /*   fscanf(f, "%ld %ld %ld %ld %ld %ld", &mvirt, &mres, &mshared, &val1, &val2, &val3); */

  /*   long int m_mvirt, m_mres, m_mshared; */

  /*   PDM_MPI_Allreduce (&mvirt, &m_mvirt, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */
  /*   PDM_MPI_Allreduce (&mres, &m_mres, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */
  /*   PDM_MPI_Allreduce (&mshared, &m_mshared, 1, PDM_MPI_LONG, PDM_MPI_MAX, pdm_mpi_comm); */

  /*   if (i_rank == 0) { */
  /*     printf("mem %d %d : %ld Ko %ld Ko %ld Ko\n", i, dn_vtx, 4*m_mvirt , 4*m_mres, 4*m_mshared); */
  /*     //      printf("mem %d %d : %ld Mo %ld Mo %ld Mo\n", i, dn_vtx, 4*m_mvirt/1024 , 4*m_mres/1024, 4*m_mshared/1024); */
  /*   } */
  /*   fclose(f); */

  /* } */

//  for (int j = 0; j < dn_vtx; j++) {
//    PDM_printf (PDM_FMT_G_NUM" %12.5e %12.5e %12.5e\n", _numabs2[j], dvtx_coord[3*j],
//                                                     dvtx_coord[3*j+1],
//                                                     dvtx_coord[3*j+2]);
//  }

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
    for (int i = 0; i < dn_edge; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, dedge_vtx[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, dedge_vtx[2*i+1]);
      PDM_printf ("\n");
    }

    PDM_printf ("dedge_face : ");
    for (int i = 0; i < dn_edge; i++) {
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
  int *dedge_vtx_idx = (int *) malloc ((dn_edge+1)*sizeof(int));

  dedge_vtx_idx[0] = 0;
  for (int i = 0; i < dn_edge; i++) {
    dedge_vtx_idx[i+1] = 2 + dedge_vtx_idx[i];
  }

  /*
   *  Split mesh i
   */

  // int ppart_id;

  int n_property_cell = 0;
  int *renum_properties_cell = NULL;
  int n_property_face = 0;
  int *renum_properties_face = NULL;

  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((num_procs+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_vtx = (PDM_g_num_t) dn_vtx;

  PDM_MPI_Allgather((void *) &_dn_vtx,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) &(distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    pdm_mpi_comm);

  // mesh->face_distrib[0] = 1;
  distrib[0] = 0;

  for (int i = 1; i < num_procs; i++) {
    distrib[i] +=  distrib[i-1];
  }


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
                                       dn_edge,
                                       dn_vtx,
                                       *n_edge_group,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dedge_face,
                                       dedge_vtx_idx,
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



  free (dvtx_coord);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dface_edge);
  free (dedge_vtx_idx);
  free (dedge_vtx);
  free (dedge_face);
  free (dedge_group_idx);
  free (dedge_group);

  PDM_gen_gnum_t* gen_gnum2 = PDM_gnum_create (3, n_part, PDM_TRUE, 1e-3, pdm_mpi_comm, PDM_OWNERSHIP_KEEP);

  double **char_length = malloc(sizeof(double *) * n_part);

  int *n_vtxs = malloc (sizeof(int) * n_part);
  PDM_g_num_t **vtx_ln_to_gns = malloc (sizeof(PDM_g_num_t *) * n_part);

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

    n_vtxs[i_part] = n_vtx;

    char_length[i_part] = malloc (sizeof(double) * n_vtx);

    for (int j = 0; j < n_vtx; j++) {
      char_length[i_part][j] = HUGE_VAL;
    }

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get(ppart,
                          i_part,
                          &cell_tag,
                          &cell_face_idx,
                          &cell_face,
                          &cell_ln_to_gn,
                          &face_tag,
                          &face_cell,
                          &face_vtx_idx,
                          &face_vtx,
                          &face_ln_to_gn,
                          &face_part_bound_proc_idx,
                          &face_part_bound_part_idx,
                          &face_part_bound,
                          &vtx_tag,
                          &vtx,
                          &vtx_ln_to_gn,
                          &face_group_idx,
                          &face_group,
                          &face_group_ln_to_gn);

    vtx_ln_to_gns[i_part] = vtx_ln_to_gn;

    for (int j = 0; j < nEdge; j++) {
      int i1 = face_vtx_idx[j];
      int n = face_vtx_idx[j+1] - i1;
      for (int k = 0; k < n; k++) {
        int vtx1 = face_vtx[i1+k] - 1;
        int vtx2 = face_vtx[i1+(k+1)%n] - 1;
        double _lEdge = (vtx[3*vtx2  ] - vtx[3*vtx1  ]) * (vtx[3*vtx2  ] - vtx[3*vtx1  ]) +
                        (vtx[3*vtx2+1] - vtx[3*vtx1+1]) * (vtx[3*vtx2+1] - vtx[3*vtx1+1]) +
                        (vtx[3*vtx2+2] - vtx[3*vtx1+2]) * (vtx[3*vtx2+2] - vtx[3*vtx1+2]);
        char_length[i_part][vtx1] = PDM_MIN (char_length[i_part][vtx1], _lEdge);
        char_length[i_part][vtx2] = PDM_MIN (char_length[i_part][vtx2], _lEdge);
      }
    }

    PDM_gnum_set_from_coords (gen_gnum2, i_part, n_vtx, vtx, char_length[i_part]);

  }

  PDM_timer_t *timer = PDM_timer_create();
  PDM_timer_resume(timer);
  PDM_gnum_compute (gen_gnum2);
  PDM_timer_hang_on(timer);
  printf("Compute gnum end %12.5es\n", PDM_timer_elapsed(timer));
  fflush(stdout);
  PDM_timer_free(timer);

  const PDM_g_num_t **_numabs = malloc (sizeof(PDM_g_num_t *) * n_part);

  // Check

  PDM_g_num_t *numabs_init = malloc(sizeof(PDM_g_num_t) * dn_vtx);

  for (int i = 0; i < dn_vtx; i++) {
    numabs_init[i] = distrib[i_rank] + 1 + i;
  }

  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP, 1.,
                                                        &numabs_init,
                                                        NULL,
                                                        &dn_vtx,
                                                        1,
                                                        pdm_mpi_comm);


  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP, 1.,
                                                        vtx_ln_to_gns,
                                                        NULL,
                                                        n_vtxs,
                                                        n_part,
                                                        pdm_mpi_comm);


  for (int i_part = 0; i_part < n_part; i_part++) {

    _numabs[i_part] = PDM_gnum_get (gen_gnum2, i_part);

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

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get(ppart,
                          i_part,
                          &cell_tag,
                          &cell_face_idx,
                          &cell_face,
                          &cell_ln_to_gn,
                          &face_tag,
                          &face_cell,
                          &face_vtx_idx,
                          &face_vtx,
                          &face_ln_to_gn,
                          &face_part_bound_proc_idx,
                          &face_part_bound_part_idx,
                          &face_part_bound,
                          &vtx_tag,
                          &vtx,
                          &vtx_ln_to_gn,
                          &face_group_idx,
                          &face_group,
                          &face_group_ln_to_gn);

  }

  int nElb1 =  PDM_part_to_block_n_elt_block_get (ptb1);

  int nElb2 =  PDM_part_to_block_n_elt_block_get (ptb2);

  assert (nElb1 == nElb2);

  PDM_g_num_t *block_numabs2;
  PDM_g_num_t *block_numabs;
  int *block_stride;

  PDM_part_to_block_exch (ptb1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                          (void **) &_numabs2,
                          &block_stride,
                          (void **) &block_numabs2);

  PDM_part_to_block_exch (ptb2,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                          (void **) _numabs,
                          &block_stride,
                          (void **) &block_numabs);

  for (int i = 0; i < nElb1; i++) {
    if (block_numabs[i] != block_numabs2[i]) {
      PDM_printf("-- diff %d : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" \n",
             i, block_numabs2[i], block_numabs[i]);
      PDM_error (__FILE__, __LINE__, 0, "Error in the generated numbering\n");
    }

  }

//  PDM_g_num_t *n1 = PDM_part_to_block_block_gnum_get (ptb1);
//  PDM_g_num_t *n2 = PDM_part_to_block_block_gnum_get (ptb2);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (char_length[i_part]);
  }
  free (char_length);

  free(_numabs);
  free(block_numabs);
  free(block_numabs2);
  free(n_vtxs);
  free(vtx_ln_to_gns);
  free(numabs_init);
  free(distrib);

  PDM_gnum_free (gen_gnum2);

  PDM_gnum_free (gen_gnum);

  PDM_part_to_block_free (ptb1);
  PDM_part_to_block_free (ptb2);

  return ppart;

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

  PDM_g_num_t   n_vtx_seg = 4;
  double        length  = 1.;
  int           n_part   = 1;
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;
  int           have_random = 0;

  int           i_rank;
  int           num_procs;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &n_part,
              &length);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &num_procs);

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

 PDM_part_free(ppart);

 if (i_rank == 0) {
   PDM_printf ("-- End\n");
 }

 PDM_MPI_Finalize ();


 return 0;

}
