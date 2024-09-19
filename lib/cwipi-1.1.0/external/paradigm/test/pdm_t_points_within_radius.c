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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_para_octree.h"
#include "pdm_timer.h"

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
     "  -t      <level>  Number of Target points (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_faceSeg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   nTgt       Number of Target points
 * \param [inout]   n_part     Number of partitions par process
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_faceSeg,
           double        *length,
           PDM_g_num_t   *nTgt,
           double        *radius,
           int           *sort_close_points,
           int           *n_max_per_leaf,
           int           *randomize,
           int           *repeat_last)
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
        long _n_faceSeg = atol(argv[i]);
        *n_faceSeg = (PDM_g_num_t) _n_faceSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
      }
    }
    else if (strcmp(argv[i], "-r") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *radius = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-mpl") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_max_per_leaf = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-sort") == 0) {
      *sort_close_points = 1;
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-repeat") == 0) {
      *repeat_last = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double _rand01(void) {
  return (double) rand() / (double) RAND_MAX;
}


static void
_gen_cloud
(
 const int           n_pts,
 const double        origin[3],
 const double        length,
 const int           n_rank,
 const int           i_rank,
 double            **pts_coord,
 int                *_n_pts
 )
{
  *_n_pts = (int) (n_pts/n_rank);
  if (i_rank < n_pts%n_rank) {
    (*_n_pts)++;
  }
  *pts_coord = malloc (sizeof(double) * 3 * (*_n_pts));

  for (int i = 0; i < (*_n_pts); i++) {
    for (int j = 0; j < 3; j++) {
      (*pts_coord)[3*i+j] = origin[j] + length * _rand01();
    }
  }
}



static void
_gen_cube_vol
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_faceSeg,
 const double        origin[3],
 const double        length,
 const int           randomize,
 int                *npts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t n_cell = n_faceSeg * n_faceSeg * n_faceSeg;

  PDM_g_num_t *distribCell = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t n_faceFace = n_faceSeg * n_faceSeg;

  // Define distribution
  distribCell[0] = 0;
  PDM_g_num_t stepCell = n_cell / n_rank;
  PDM_g_num_t remainderCell = n_cell % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i] = stepCell;
    const int i1 = i - 1;
    if (i1 < remainderCell)
      distribCell[i] += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i] += distribCell[i-1];
  }

  PDM_g_num_t _dn_cell = distribCell[i_rank+1] - distribCell[i_rank];
  *npts = (int) _dn_cell;

  const double step = length / (double) n_faceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dn_cell);
  *coord = malloc (sizeof(double)      * _dn_cell * 3);

  if (randomize) {
    for (int i = 0; i < *npts; i++) {
      (*g_num)[i] = 1 + i + distribCell[i_rank];
      for (int j = 0; j < 3; j++) {
        (*coord)[3*i+j] = origin[j] + length * _rand01();
      }
    }
  }

  else {
    int _npts = 0;
    for (PDM_g_num_t g = distribCell[i_rank]; g < distribCell[i_rank+1]; g++) {
      PDM_g_num_t i = g % n_faceSeg;
      PDM_g_num_t j = ((g - i) % n_faceFace) / n_faceSeg;
      PDM_g_num_t k = (g - i - n_faceSeg * j) / n_faceFace;

      (*coord)[3 * _npts    ] = (i + 0.5) * step + origin[0];
      (*coord)[3 * _npts + 1] = (j + 0.5) * step + origin[1];
      (*coord)[3 * _npts + 2] = (k + 0.5) * step + origin[2];
      (*g_num)[_npts++] = g + 1;//1 + i + n_faceSeg * j + n_faceFace * k;
    }
  }

  free (distribCell);
}




static void
_points_within_radius
(
 PDM_MPI_Comm         comm,
 const int            n_max_per_leaf,
 const int            sort_close_points,
 const int            n_part_src,
 const int           *n_src,
 const double       **src_coord,
 const PDM_g_num_t  **src_g_num,
 const int            n_part_tgt,
 const int           *n_tgt,
 const double       **tgt_coord,
 const PDM_g_num_t  **tgt_g_num,
 const double       **tgt_radius,
 int               ***close_points_idx,
 PDM_g_num_t       ***close_points_g_num,
 double            ***close_points_dist2
 )
{
  int i_rank;//
  PDM_MPI_Comm_rank (comm, &i_rank);//

  /* Build parallel octree */
  const int depth_max = 31;

  PDM_para_octree_t *octree = PDM_para_octree_create (n_part_src,
                                                      depth_max,
                                                      n_max_per_leaf,
                                                      0,
                                                      comm);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_para_octree_point_cloud_set (octree,
                                     i_part,
                                     n_src[i_part],
                                     src_coord[i_part],
                                     src_g_num[i_part]);
  }

  /* Compute global extents of source and target point clouds */
  double local_min[3] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double local_max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i_part = 0; i_part < n_part_src; i_part++) {
    const double *x = src_coord[i_part];
    for (int i = 0; i < n_src[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    const double *x = tgt_coord[i_part];
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  double global_extents[6];
  PDM_MPI_Allreduce(local_min, global_extents,     3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(local_max, global_extents + 3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_para_octree_build (octree, global_extents);
  //PDM_para_octree_dump (octree);
  PDM_para_octree_dump_times (octree);




  /* Concatenate partitions */
  int _n_tgt = 0;

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    _n_tgt += n_tgt[i_part];
  }

  double      *_tgt_coord   = malloc (sizeof(double)      * _n_tgt * 3);
  PDM_g_num_t *_tgt_g_num   = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  double      *_tgt_radius2 = malloc (sizeof(double)      * _n_tgt);

  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        _tgt_coord[_n_tgt + 3*i + j] = tgt_coord[i_part][3*i + j];
      }
      _tgt_g_num[_n_tgt + i] = tgt_g_num[i_part][i];
      _tgt_radius2[_n_tgt + i] = tgt_radius[i_part][i] * tgt_radius[i_part][i];
    }
    _n_tgt += n_tgt[i_part];
  }

  /* Search source points within radius */
  int         *_close_pts_idx   = NULL;
  PDM_g_num_t *_close_pts_g_num = NULL;
  double      *_close_pts_dist2 = NULL;

  PDM_para_octree_points_within_radius (octree,
                                        sort_close_points,
                                        _n_tgt,
                                        _tgt_coord,
                                        _tgt_g_num,
                                        _tgt_radius2,
                                        &_close_pts_idx,
                                        &_close_pts_g_num,
                                        &_close_pts_dist2);

  if (0) {
    for (int i = 0; i < _n_tgt; i++) {
      printf("tgt point ("PDM_FMT_G_NUM") :", _tgt_g_num[i]);
      for (int j = _close_pts_idx[i]; j < _close_pts_idx[i+1]; j++) {
        printf(" ("PDM_FMT_G_NUM")", _close_pts_g_num[j]);
      }
      printf("\n");
    }
  }
  free (_tgt_coord);
  free (_tgt_g_num);
  free (_tgt_radius2);

  /* Restore partitions */
  *close_points_idx   = (int **)         malloc (sizeof(int *)         * n_part_tgt);
  *close_points_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *close_points_dist2 = (double **)      malloc (sizeof(double *)      * n_part_tgt);
  int idx_part = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {

    int length_close_points = _close_pts_idx[idx_part + n_tgt[i_part]] - _close_pts_idx[idx_part];

    (*close_points_g_num)[i_part] = malloc (sizeof(PDM_g_num_t) * length_close_points);
    (*close_points_dist2)[i_part] = malloc (sizeof(double)      * length_close_points);
    (*close_points_idx)[i_part] = malloc (sizeof(int) * (n_tgt[i_part]+1));
    int *cp_idx = (*close_points_idx)[i_part];
    cp_idx[0] = 0;

    for (int i = 0; i < n_tgt[i_part]; i++) {
      int k = 0;
      for (int j = _close_pts_idx[idx_part+i]; j < _close_pts_idx[idx_part+i+1]; j++) {
        (*close_points_g_num)[i_part][cp_idx[i] + k] = _close_pts_g_num[j];
        (*close_points_dist2)[i_part][cp_idx[i] + k] = _close_pts_dist2[j];
        k++;
      }

      cp_idx[i+1] = cp_idx[i] + k;
    }

    idx_part += n_tgt[i_part];
  }
  free (_close_pts_idx);
  free (_close_pts_g_num);
  free (_close_pts_dist2);

  /* Free parallel octree */
  PDM_para_octree_free (octree);
}



static void
_write_point_cloud
(
 const char        *filename,
 const char        *header,
 const int          n_pts,
 const double       coord[],
 const PDM_g_num_t  g_num[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%12.5e ", coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  fclose(f);
}


static void
_read_point_cloud
(
 const char   *filename,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num
 )
{
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Unable to open %s", filename);
  }

  char line[999];
  while (fgets(line, sizeof(line), f) != NULL) {
    if (strcmp(line,"\n") == 0 && strcmp(line,"\r\n") == 0) {
      continue;
    }

    if (strstr(line, "POINTS") != NULL) {
      int stat = sscanf(line,
                        "%*[^0123456789]%d%*[^0123456789]",
                        n_pts);
      assert (stat);

      *coord = malloc (sizeof(double) * (*n_pts) * 3);
      for (int i = 0; i < *n_pts; i++) {
        fscanf(f, "%lf %lf %lf",
               *coord + 3*i,
               *coord + 3*i + 1,
               *coord + 3*i + 2);
      }
    }

    if (strstr(line, "CELL_DATA") != NULL) {

      *g_num = malloc (sizeof(PDM_g_num_t) * (*n_pts));
      fgets(line, sizeof(line), f);
      fgets(line, sizeof(line), f);
      for (int i = 0; i < *n_pts; i++) {
        fscanf(f, PDM_FMT_G_NUM, *g_num + i);
      }
    }
  }

  fclose(f);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);


  srand (time(NULL) + i_rank);

  /*
   *  Set default values
   */
  PDM_g_num_t n_face_seg        = 10;
  double      length            = 1.;
  PDM_g_num_t n_tgt             = 10;
  double      radius            = 0.1;
  int         sort_close_points = 0;
  int         n_max_per_leaf    = 1;
  int         randomize         = 0;
  int         repeat_last       = 0;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &n_face_seg,
              &length,
              &n_tgt,
              &radius,
              &sort_close_points,
              &n_max_per_leaf,
              &randomize,
              &repeat_last);


  double origin[3] = {0., 0., 0.};
  if (0) {
    if (i_rank == 0) {
      for (int i = 0; i < 3; i++) {
        origin[i] = 2.*_rand01() - 1.;
      }
    }

    PDM_MPI_Bcast (origin, 3, PDM_MPI_DOUBLE, 0, PDM_MPI_COMM_WORLD);
  }

  double aniso[3] = {1., 1., 1.};
  if (0) {
    if (i_rank == 0) {
      for (int i = 0; i < 3; i++) {
        double r = 2.*_rand01() - 1.;
        aniso[i] = pow(2., r);
      }
      printf("aniso = %f %f %f\n", aniso[0], aniso[1], aniso[2]);
    }

    PDM_MPI_Bcast (aniso, 3, PDM_MPI_DOUBLE, 0, PDM_MPI_COMM_WORLD);
  }


  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  double      *tgt_coord = NULL;
  PDM_g_num_t *tgt_g_num = NULL;
  int _n_tgt;

  if (repeat_last) {
    char filename[999];
    sprintf(filename, "tgt_%3.3d.vtk", i_rank);

    _read_point_cloud (filename,
                       &_n_tgt,
                       &tgt_coord,
                       &tgt_g_num);
  }

  else {
    _gen_cloud (n_tgt,
                origin,
                length,
                n_rank,
                i_rank,
                &tgt_coord,
                &_n_tgt);

    for (int i = 0; i < _n_tgt; i++) {
      for (int j = 0; j < 3; j++) {
        tgt_coord[3*i+j] = origin[j] + aniso[j] * (tgt_coord[3*i+j] - origin[j]);
      }
      //printf("[%d] tgt %f %f %f\n", i_rank, tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2]);
    }

    int id_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

    double *tgt_char_length = malloc (sizeof(double) * _n_tgt);

    for (int i = 0; i < _n_tgt; i++) {
      tgt_char_length[i] = length * 1.e-6;
    }

    PDM_gnum_set_from_coords (id_gnum, 0, _n_tgt, tgt_coord, tgt_char_length);

    PDM_gnum_compute (id_gnum);

    tgt_g_num = PDM_gnum_get (id_gnum, 0);

    PDM_gnum_free (id_gnum);
    free (tgt_char_length);
  }


  /*
   *  Search radius
   */
  double *tgt_radius = malloc (sizeof(double) * _n_tgt);
  if (radius < 0) {
    for (int i = 0; i < _n_tgt; i++) {
      tgt_radius[i] = -radius * _rand01();
    }
  }

  else {
    for (int i = 0; i < _n_tgt; i++) {
      tgt_radius[i] = radius;
    }
  }



  /*
   *  Define the source point cloud
   */
  int n_part_src = 1;
  int _n_src;
  double *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;

  if (repeat_last) {
    char filename[999];
    sprintf(filename, "src_%3.3d.vtk", i_rank);

    _read_point_cloud (filename,
                       &_n_src,
                       &src_coord,
                       &src_g_num);
  }

  else {

    _gen_cube_vol (PDM_MPI_COMM_WORLD,
                   n_face_seg,
                   origin,
                   length,
                   randomize,
                   &_n_src,
                   &src_g_num,
                   &src_coord);

    for (int i = 0; i < _n_src; i++) {
      for (int j = 0; j < 3; j++) {
        src_coord[3*i+j] = origin[j] + aniso[j] * (src_coord[3*i+j] - origin[j]);
      }
      //printf("[%d] src %f %f %f\n", i_rank, src_coord[3*i], src_coord[3*i+1], src_coord[3*i+2]);
    }
  }


  if (!repeat_last) {
    char filename[999];

    sprintf(filename, "tgt_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "tgt",
                        _n_tgt,
                        tgt_coord,
                        tgt_g_num);

    sprintf(filename, "src_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "src",
                        _n_src,
                        src_coord,
                        src_g_num);
  }


  PDM_g_num_t n_local[2], n_global[2];
  n_local[0] = (PDM_g_num_t) _n_src;
  n_local[1] = (PDM_g_num_t) _n_tgt;
  PDM_MPI_Reduce (n_local, n_global, 2, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, 0, PDM_MPI_COMM_WORLD);
  /*if (i_rank == 0) {
    //printf("n_procs %d, n_src "PDM_FMT_G_NUM", n_tgt "PDM_FMT_G_NUM", method %d, surf src %d, random src %d\n", n_rank, n_global[0], n_global[1], method, surf_source, randomize);
    printf("n_procs = %d\nn_src_total = "PDM_FMT_G_NUM"\nn_tgt_total = "PDM_FMT_G_NUM"\nrandom src = %d\n", n_rank, n_global[0], n_global[1], randomize);
    }*/


  /*
   *  Compute points within radius
   */
  int         **close_points_idx   = NULL;
  PDM_g_num_t **close_points_g_num = NULL;
  double      **close_points_dist2 = NULL;
  _points_within_radius (PDM_MPI_COMM_WORLD,
                         n_max_per_leaf,
                         sort_close_points,
                         n_part_src,
                         (const int *) &_n_src,
                         (const double **) &src_coord,
                         (const PDM_g_num_t **) &src_g_num,
                         n_part_tgt,
                         (const int *) &_n_tgt,
                         (const double **) &tgt_coord,
                         (const PDM_g_num_t **) &tgt_g_num,
                         (const double **) &tgt_radius,
                         &close_points_idx,
                         &close_points_g_num,
                         &close_points_dist2);

  /*
   *  Check
   */
  if (0) {
    printf("Check:\n");
    for (int i = 0; i < _n_tgt; i++) {
      printf("tgt point ("PDM_FMT_G_NUM") :", tgt_g_num[i]);
      for (int j = close_points_idx[0][i]; j < close_points_idx[0][i+1]; j++) {
        printf(" ("PDM_FMT_G_NUM", %f)", close_points_g_num[0][j], close_points_dist2[0][j]);
      }
      printf("\n");
    }
  }

  if (!randomize) {
    if (i_rank == 0) {
      printf("-- Check\n");
      fflush(stdout);
    }

    int ijk0;
    int ijk_lo[3];
    int ijk_hi[3];
    double cell_side = length / ((double) n_face_seg);
    double coord[3];

    PDM_g_num_t _n_wrong = 0;
    for (int itgt = 0; itgt < _n_tgt; itgt++) {

      double radius2 = tgt_radius[itgt] * tgt_radius[itgt];
      int missed = 0;

      // Find base cell
      for (int idim = 0; idim < 3; idim++) {
        int n_cell = (int) ceil(tgt_radius[itgt] / (cell_side * aniso[idim]));
        ijk0 = (int) floor((tgt_coord[3*itgt+idim] - origin[idim])/ (cell_side * aniso[idim]));
        ijk0 = PDM_MIN (PDM_MAX (ijk0, 0), n_face_seg-1);

        ijk_lo[idim] = PDM_MAX (ijk0 - n_cell, 0);
        ijk_hi[idim] = PDM_MIN (ijk0 + n_cell, n_face_seg-1);
      }

      // inspect search region
      for (int k = ijk_lo[2]; k < ijk_hi[2]; k++) {
        coord[2] = (k + 0.5) * cell_side * aniso[2] + origin[2];
        for (int j = ijk_lo[1]; j < ijk_hi[1]; j++) {
          coord[1] = (j + 0.5) * cell_side * aniso[1] + origin[1];
          for (int i = ijk_lo[0]; i < ijk_hi[0]; i++) {
            coord[0] = (i + 0.5) * cell_side * aniso[0] + origin[0];

            double dist2 = 0.;
            for (int idim = 0; idim < 3; idim++) {
              double delta = tgt_coord[3*itgt+idim] - coord[idim];
              dist2 += delta * delta;
            }

            if (dist2 < radius2) {
              PDM_g_num_t g_num = 1 + i + n_face_seg*(j + n_face_seg*k);
              int found = 0;
              for (int l = close_points_idx[0][itgt]; l < close_points_idx[0][itgt+1]; l++) {
                if (g_num == close_points_g_num[0][l]) {
                  found = 1;
                  break;
                }
              }

              if (!found) {
                missed = 1;
                printf("[%d] ERROR pt ("PDM_FMT_G_NUM") (%f %f %f) : missed point ("PDM_FMT_G_NUM") at dist = %.2g times radius (%f)\n", i_rank, tgt_g_num[itgt], tgt_coord[3*itgt], tgt_coord[3*itgt+1], tgt_coord[3*itgt+2], g_num, sqrt(dist2)/tgt_radius[itgt], tgt_radius[itgt]);
              }
            }

          }
        }
      }

      if (missed) {
        printf("[%d] some missed points for ("PDM_FMT_G_NUM"), only found:", i_rank, tgt_g_num[itgt]);
        for (int l = close_points_idx[0][itgt]; l < close_points_idx[0][itgt+1]; l++) {
          printf(" ("PDM_FMT_G_NUM")", close_points_g_num[0][l]);
        }
        printf("\n");
      }
      _n_wrong += missed;
    }


    PDM_g_num_t n_wrong;
    PDM_MPI_Reduce (&_n_wrong,
                    &n_wrong,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    PDM_MPI_SUM,
                    0,
                    PDM_MPI_COMM_WORLD);
    if (i_rank == 0) {
      PDM_g_num_t wrong_percentage = 100 * n_wrong;
      if (n_tgt > 0) {
        wrong_percentage /= n_tgt;
      }

      printf("\nn_wrong = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM"%%)\n",
             n_wrong,
             n_tgt,
             wrong_percentage);
      fflush(stdout);
      assert (n_wrong == 0);
    }
  }

  /*
   *  Finalize
   */
  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);
  free (tgt_radius);
  for (int i = 0; i < n_part_tgt; i++) {
    free (close_points_idx[i]);
    free (close_points_g_num[i]);
    free (close_points_dist2[i]);
  }
  free (close_points_idx);
  free (close_points_g_num);
  free (close_points_dist2);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
