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
#include "pdm_dcube_gen.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_polygon.h"

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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_r_seg,
           PDM_g_num_t   *n_t_seg,
           PDM_g_num_t   *n_z_seg,
           double        *b,
           int           *n_part,
           int           *post,
           int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nr") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_r_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }

    else if (strcmp(argv[i], "-nt") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_t_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }

    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_z_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }

    else if (strcmp(argv[i], "-b") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *b = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_gen_distributed_mesh
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_r_seg,
 const PDM_g_num_t   n_t_seg,
 const PDM_g_num_t   n_z_seg,
 const double        r_min,
 const double        r_max,
 const double        t_min,
 const double        t_max,
 const double        z_min,
 const double        z_max,
 int                *dn_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *n_face_group,
 PDM_g_num_t       **dface_cell,
 int               **dface_vtx_idx,
 PDM_g_num_t       **dface_vtx,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group,
 double            **dvtx_coord)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);



  PDM_g_num_t n_vtx = n_r_seg * n_t_seg * n_z_seg;
  PDM_g_num_t n_face_rt = (n_r_seg - 1) * (n_t_seg - 1);
  PDM_g_num_t n_face_tz = (n_t_seg - 1) * (n_z_seg - 1);
  PDM_g_num_t n_face_zr = (n_z_seg - 1) * (n_r_seg - 1);
  PDM_g_num_t n_face = n_face_rt * n_z_seg + n_face_tz * n_r_seg + n_face_zr * n_t_seg;
  PDM_g_num_t n_cell = (n_r_seg - 1) * (n_t_seg - 1) * (n_z_seg - 1);
  PDM_g_num_t n_face_lim = 2 * (n_face_rt + n_face_tz + n_face_zr);

  double step_r = (r_max - r_min) / (double) (n_r_seg - 1);
  double step_t = (t_max - t_min) / (double) (n_t_seg - 1);
  double step_z = (z_max - z_min) / (double) (n_z_seg - 1);

  PDM_g_num_t *distrib_vtx      = (PDM_g_num_t *) malloc ((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_face     = (PDM_g_num_t *) malloc ((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_cell     = (PDM_g_num_t *) malloc ((n_rank + 1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *distrib_face_lim = (PDM_g_num_t *) malloc ((n_rank + 1) * sizeof(PDM_g_num_t));

  // Define distribution

  distrib_vtx[0]      = 0;
  distrib_face[0]     = 0;
  distrib_cell[0]     = 0;
  distrib_face_lim[0] = 0;

  PDM_g_num_t step_vtx      = n_vtx / n_rank;
  PDM_g_num_t remainder_vtx = n_vtx % n_rank;

  PDM_g_num_t step_face      = n_face / n_rank;
  PDM_g_num_t remainder_face = n_face % n_rank;

  PDM_g_num_t step_cell      = n_cell / n_rank;
  PDM_g_num_t remainder_cell = n_cell % n_rank;

  PDM_g_num_t step_face_im       = n_face_lim / n_rank;
  PDM_g_num_t remainder_face_lim = n_face_lim % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]     = step_vtx;
    distrib_face[i]    = step_face;
    distrib_cell[i]    = step_cell;
    distrib_face_lim[i] = step_face_im;
    const int i1 = i - 1;
    if (i1 < remainder_vtx)
      distrib_vtx[i]  += 1;
    if (i1 < remainder_face)
      distrib_face[i]  += 1;
    if (i1 < remainder_cell)
      distrib_cell[i]  += 1;
    if (i1 < remainder_face_lim)
      distrib_face_lim[i]  += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib_vtx[i]  += distrib_vtx[i-1];
    distrib_face[i] += distrib_face[i-1];
    distrib_cell[i] += distrib_cell[i-1];
    distrib_face_lim[i] += distrib_face_lim[i-1];
  }

  PDM_g_num_t _dn_cell = distrib_cell[i_rank+1] - distrib_cell[i_rank];
  *dn_cell = (int) _dn_cell;
  PDM_g_num_t _dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];
  *dn_face = (int) _dn_face;
  PDM_g_num_t _dn_vtx = distrib_vtx[i_rank+1]  - distrib_vtx[i_rank];
  *dn_vtx = (int) _dn_vtx;
  PDM_g_num_t _dn_face_lim = distrib_face_lim[i_rank+1] - distrib_face_lim[i_rank];
  int dn_face_lim = (int) _dn_face_lim;

  *n_face_group = 6;

  *dface_cell      = (PDM_g_num_t *) malloc (2*(*dn_face    )    * sizeof(PDM_g_num_t *));
  *dface_vtx_idx   = (int         *) malloc (  (*dn_face + 1)    * sizeof(int         *));
  *dface_vtx       = (PDM_g_num_t *) malloc (4*(*dn_face    )    * sizeof(PDM_g_num_t *));
  *dvtx_coord      = (double      *) malloc (3*(*dn_vtx     )    * sizeof(double      *));
  *dface_group_idx = (int         *) malloc ((*n_face_group + 1) * sizeof(int         *));
  *dface_group     = (PDM_g_num_t *) malloc (dn_face_lim         * sizeof(PDM_g_num_t *));

  PDM_g_num_t  *_dface_cell      = *dface_cell;
  int          *_dface_vtx_idx   = *dface_vtx_idx;
  PDM_g_num_t  *_dface_vtx       = *dface_vtx;
  double       *_dvtx_coord      = *dvtx_coord;
  int          *_dface_group_idx = *dface_group_idx;
  PDM_g_num_t  *_dface_group     = *dface_group;

  _dface_vtx_idx[0] = 0;
  for (int i = 1; i < *dn_face + 1; i++) {
    _dface_vtx_idx[i] = 4 + _dface_vtx_idx[i-1];
  }


  // Coordinates
  PDM_g_num_t n_rt = n_r_seg * n_t_seg;
  const PDM_g_num_t b_vtx_z = distrib_vtx[i_rank] / n_rt;
  const PDM_g_num_t r_vtx_z = distrib_vtx[i_rank] % n_rt;

  const PDM_g_num_t b_vtx_t = r_vtx_z / n_r_seg;
  const PDM_g_num_t b_vtx_r = r_vtx_z % n_r_seg;

  int i_vtx = 0;
  int cpt   = 0;

  for(PDM_g_num_t k = b_vtx_z; k < n_z_seg; k++) {
    PDM_g_num_t _b_vtx_t = 0;
    if (k == b_vtx_z)
      _b_vtx_t = b_vtx_t;

    double z = z_min + k * step_z;

    for(PDM_g_num_t j = _b_vtx_t; j < n_t_seg; j++) {
      PDM_g_num_t _b_vtx_r = 0;

      double t = t_min + j * step_t;
      double ct = cos(t);
      double st = sin(t);
      if ((k == b_vtx_z) && (j == b_vtx_t))
        _b_vtx_r = b_vtx_r;
      for(PDM_g_num_t i = _b_vtx_r; i < n_r_seg; i++) {

        double r = r_min + i * step_r;

        _dvtx_coord[3 * i_vtx    ] = r * ct;
        _dvtx_coord[3 * i_vtx + 1] = r * st;
        _dvtx_coord[3 * i_vtx + 2] = z;
        cpt   += 1;
        i_vtx += 1;
        if (cpt == *dn_vtx)
          break;
      }
      if (cpt == *dn_vtx)
        break;
    }
    if (cpt == *dn_vtx)
      break;
  }

  //
  // face_vtx et face_cell

  cpt = 0;

  PDM_g_num_t distrib_serie[7] = {0,
                                  n_face_rt * n_z_seg,
                                  n_face_rt * n_z_seg + n_face_tz * n_r_seg,
                                  n_face,
                                  0,
                                  0,
                                  0};

  PDM_g_num_t i_serie = 0, r_serie = 0;
  for (int i = 0; i < 3; i++) {
    if (distrib_face[i_rank] >= distrib_serie[i] &&
        distrib_face[i_rank] <  distrib_serie[i+1]) {
      i_serie = i;
      r_serie = distrib_face[i_rank] - distrib_serie[i];
    }
  }


  PDM_g_num_t b1 = 0;
  PDM_g_num_t r1 = 0;

  PDM_g_num_t b2 = 0;
  PDM_g_num_t b3 = 0;

  if (i_serie == 0) {
    b1 = r_serie / n_face_rt;
    r1 = r_serie % n_face_rt;

    b2 = r1 / (n_r_seg - 1);
    b3 = r1 % (n_r_seg - 1);
  }
  else if (i_serie == 1) {
    b1 = r_serie / n_face_tz;
    r1 = r_serie % n_face_tz;

    b2 = r1 / (n_t_seg - 1);
    b3 = r1 % (n_t_seg - 1);
  }
  else {
    b1 = r_serie / n_face_zr;
    r1 = r_serie % n_face_zr;

    b2 = r1 / (n_z_seg - 1);
    b3 = r1 % (n_z_seg - 1);
  }

  if (i_serie == 0) {
    // Faces zmin -> zmax

    for(PDM_g_num_t k = b1; k < n_z_seg; k++) {
      PDM_g_num_t _b2 = 0;
      if (k == b1)
        _b2 = b2;
      for(PDM_g_num_t j = _b2; j < n_t_seg - 1; j++) {
        PDM_g_num_t _b3 = 0;
        if ((k == b1) && (j == b2))
          _b3 = b3;
        for(PDM_g_num_t i = _b3; i < n_r_seg - 1; i++) {

          _dface_vtx[cpt * 4    ] = k * n_rt + (    j * n_r_seg + i + 1);
          _dface_vtx[cpt * 4 + 1] = k * n_rt + ((j+1) * n_r_seg + i + 1);
          _dface_vtx[cpt * 4 + 2] = k * n_rt + ((j+1) * n_r_seg + i + 2);
          _dface_vtx[cpt * 4 + 3] = k * n_rt + (    j * n_r_seg + i + 2);

          if (k == 0) {
            _dface_cell[2*cpt + 0] = j * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = 0;

            _dface_vtx[cpt * 4    ] = k * n_rt + (    j * n_r_seg + i + 2);
            _dface_vtx[cpt * 4 + 1] = k * n_rt + ((j+1) * n_r_seg + i + 2);
            _dface_vtx[cpt * 4 + 2] = k * n_rt + ((j+1) * n_r_seg + i + 1);
            _dface_vtx[cpt * 4 + 3] = k * n_rt + (    j * n_r_seg + i + 1);

          } else if (k == n_z_seg - 1) {
            _dface_cell[2*cpt + 0] = (k-1) * n_face_rt + j * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = 0;

          } else {
            _dface_cell[2*cpt + 0] = (k-1) * n_face_rt + j * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] =     k * n_face_rt + j * (n_r_seg - 1) + i + 1;
          }
          cpt += 1;
          if (cpt == *dn_face)
            break;
        }
        if (cpt == *dn_face)
          break;
      }
      if (cpt == *dn_face)
        break;
    }
    b1 = 0;
    b2 = 0;
    b3 = 0;
  }


  if ((i_serie == 1) || ((i_serie == 0) &&  (cpt != *dn_face))) {
    // Faces xmin -> xmax

    for(PDM_g_num_t i = b1; i < n_r_seg; i++) {
      PDM_g_num_t _b2 = 0;
      if (i == b1)
        _b2 = b2;
      for(PDM_g_num_t k = _b2; k < n_z_seg - 1; k++) {
        PDM_g_num_t _b3 = 0;
        if ((i == b1) && (k == b2))
          _b3 = b3;
        for(PDM_g_num_t j = _b3; j < n_t_seg - 1; j++) {

          _dface_vtx[cpt * 4    ] = (k+1) * n_rt +     j * n_r_seg + i + 1;
          _dface_vtx[cpt * 4 + 1] = (k+1) * n_rt + (j+1) * n_r_seg + i + 1;
          _dface_vtx[cpt * 4 + 2] =     k * n_rt + (j+1) * n_r_seg + i + 1;
          _dface_vtx[cpt * 4 + 3] =     k * n_rt +     j * n_r_seg + i + 1;

          if (i == 0) {
            _dface_cell[2*cpt + 0] = k * n_face_rt + j * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = 0;

            _dface_vtx[cpt * 4    ] =     k * n_rt +     j * n_r_seg + i + 1;
            _dface_vtx[cpt * 4 + 1] =     k * n_rt + (j+1) * n_r_seg + i + 1;
            _dface_vtx[cpt * 4 + 2] = (k+1) * n_rt + (j+1) * n_r_seg + i + 1;
            _dface_vtx[cpt * 4 + 3] = (k+1) * n_rt +     j * n_r_seg + i + 1;

          } else if (i == n_r_seg - 1) {
            _dface_cell[2*cpt + 0] = k * n_face_rt + j * (n_r_seg - 1) + i;
            _dface_cell[2*cpt + 1] = 0;

          } else {
            _dface_cell[2*cpt + 0] = k * n_face_rt + j * (n_r_seg - 1) + i ;
            _dface_cell[2*cpt + 1] = k * n_face_rt + j * (n_r_seg - 1) + i + 1;
          }
          cpt += 1;
          if (cpt == *dn_face)
            break;
        }
        if (cpt == *dn_face)
          break;
      }
      if (cpt == *dn_face)
        break;
    }
    b1 = 0;
    b2 = 0;
    b3 = 0;
  }



  if ((i_serie == 2) || ((i_serie == 1 || i_serie == 0) && (cpt != *dn_face))) {
    // Faces ymin -> ymax

    for(PDM_g_num_t j = b1; j < n_t_seg; j++) {
      PDM_g_num_t _b2 = 0;
      if (j == b1)
        _b2 = b2;
      for(PDM_g_num_t i = _b2; i < n_r_seg - 1; i++) {
        PDM_g_num_t _b3 = 0;
        if ((j == b1) && (i == b2))
          _b3 = b3;
        for(PDM_g_num_t k = _b3; k < n_z_seg - 1; k++) {
          _dface_vtx[cpt * 4    ] =     k * n_rt + j * n_r_seg + i + 1    ;
          _dface_vtx[cpt * 4 + 1] =     k * n_rt + j * n_r_seg + i + 1 + 1;
          _dface_vtx[cpt * 4 + 2] = (k+1) * n_rt + j * n_r_seg + i + 1 + 1;
          _dface_vtx[cpt * 4 + 3] = (k+1) * n_rt + j * n_r_seg + i + 1    ;

          if (j == 0) {
            _dface_cell[2*cpt + 0] = k * n_face_rt + j * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = 0;

            _dface_vtx[cpt * 4    ] = (k+1) * n_rt + j * n_r_seg + i + 1    ;
            _dface_vtx[cpt * 4 + 1] = (k+1) * n_rt + j * n_r_seg + i + 1 + 1;
            _dface_vtx[cpt * 4 + 2] =     k * n_rt + j * n_r_seg + i + 1 + 1;
            _dface_vtx[cpt * 4 + 3] =     k * n_rt + j * n_r_seg + i + 1    ;

          } else if (j == n_t_seg - 1) {
            _dface_cell[2*cpt + 0] =  k * n_face_rt + (j-1) * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = 0;

          } else {
            _dface_cell[2*cpt + 0] = k * n_face_rt + (j-1) * (n_r_seg - 1) + i + 1;
            _dface_cell[2*cpt + 1] = k * n_face_rt +     j * (n_r_seg - 1) + i + 1;
          }
          cpt += 1;
          if (cpt == *dn_face)
            break;
        }
        if (cpt == *dn_face)
          break;
      }
      if (cpt == *dn_face)
        break;
    }
  }


  // Faces limite

  cpt = 0;
  PDM_g_num_t b_face;
  int cpt1 = 0;
  int cpt3 = 0;
  int first_group = 0;

  distrib_serie[1] = n_face_rt;
  distrib_serie[2] = distrib_serie[1] + n_face_rt;
  distrib_serie[3] = distrib_serie[2] + n_face_tz;
  distrib_serie[4] = distrib_serie[3] + n_face_tz;
  distrib_serie[5] = distrib_serie[4] + n_face_zr;
  distrib_serie[6] = n_face_lim;

  for (int i = 0; i < 6; i++) {
    if (distrib_face_lim[i_rank] >= distrib_serie[i] &&
        distrib_face_lim[i_rank] <  distrib_serie[i+1]) {
      i_serie = i;
      r_serie = distrib_face_lim[i_rank] - distrib_serie[i];
    }
  }

  for (int i = 0; i < *n_face_group + 1; i++)
    _dface_group_idx[i] = 0;


  if (i_serie == 0) {
    // Faces zmin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = 0;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_t_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_r_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_r_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[1] = cpt - cpt1;

    /* if (cpt == dn_face_lim) */
    /*   break; */
    first_group = 0;
  }


  if ((i_serie == 1) || ((i_serie == 0) && (cpt != dn_face_lim))) {
    // Faces zmax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_rt * (n_z_seg - 1);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_t_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_r_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_r_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[2] = cpt - cpt1;
     first_group = 0;
  }


  if ((i_serie == 2) || (((i_serie == 0) || (i_serie == 1)) && (cpt != dn_face_lim))) {
    // Faces xmin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_rt * n_z_seg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_z_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_t_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_t_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[3] = cpt - cpt1;
     first_group = 0;
  }


  if ((i_serie == 3) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2)) && (cpt != dn_face_lim))) {
    // Faces xmax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_rt * n_z_seg + n_face_tz * (n_r_seg - 1);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_z_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_t_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_t_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[4] = cpt - cpt1;
    first_group = 0;
  }


  if ((i_serie == 4) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2) || (i_serie == 3)) && (cpt != dn_face_lim))) {
    // Faces ymin

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_rt * n_z_seg + n_face_tz * n_r_seg;

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_r_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_z_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_z_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[5] = cpt - cpt1;
    first_group = 0;
  }


  if ((i_serie == 5) || (((i_serie == 0) || (i_serie == 1)  || (i_serie == 2) || (i_serie == 3) || (i_serie == 4)) && (cpt != dn_face_lim))) {
    // Faces ymax

    if (cpt == 0)
      first_group = 1;

    cpt1 = cpt;

    b_face = n_face_rt * n_z_seg + n_face_tz * n_r_seg + n_face_zr * (n_t_seg - 1);

    cpt3 = 0;
    for(PDM_g_num_t j = 0; j < n_r_seg - 1; j++) {
      for(PDM_g_num_t i = 0; i < n_z_seg - 1; i++) {
        cpt3 += 1;
        if (!first_group || (first_group && ((cpt3 - 1)  >= r_serie))) {
          _dface_group[cpt] = b_face + j * (n_z_seg - 1) + i + 1;
          cpt += 1;
          if (cpt == dn_face_lim)
            break;
        }
      }
      if (cpt == dn_face_lim)
        break;
    }

    _dface_group_idx[6] = cpt - cpt1;
    first_group = 0;

  }

  for (int i = 1; i < *n_face_group + 1; i++)
    _dface_group_idx[i] += _dface_group_idx[i-1];

  free(distrib_vtx);
  free(distrib_face);
  free(distrib_cell);
  free(distrib_face_lim);
}





static void _export_vtk_cloud
(
 const char         *filename,
 int                 nVtx,
 double             *vtxCoord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", nVtx);
  for (int i = 0; i < nVtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%lf ", vtxCoord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", nVtx, 2*nVtx);
  for (int i = 0; i < nVtx; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", nVtx);
  for (int i = 0; i < nVtx; i++) {
    fprintf(f, "1\n");
  }

  fclose(f);
}


static void _export_vtk_surface
(
 const char         *filename,
 int                 nFace,
 int                *faceVtxIdx,
 int                *faceVtx,
 int                *face_data,
 int                 nVtx,
 double             *vtxCoord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", nVtx);
  for (int i = 0; i < nVtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%lf ", vtxCoord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", nFace, nFace + faceVtxIdx[nFace]);
  for (int i = 0; i < nFace; i++) {
    fprintf(f, "%d ", faceVtxIdx[i+1] - faceVtxIdx[i]);
    for (int j = faceVtxIdx[i]; j < faceVtxIdx[i+1]; j++) {
      fprintf(f, "%d ", faceVtx[j] - 1);
    }
    fprintf(f, "\n");
  }

  if (face_data != NULL) {
    fprintf(f, "CELL_DATA %d\n", nFace);
    fprintf(f, "SCALARS scalarf int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < nFace; i++) {
      fprintf(f, "%d\n", face_data[i]);
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

  /*
   *  Set default values
   */

  PDM_g_num_t n_r_seg  = 10;
  PDM_g_num_t n_t_seg  = 5;
  PDM_g_num_t n_z_seg  = 20;
  double      r_min    = 0.3;
  double      r_max    = 0.7;
  double      t_min    = -1.0;
  double      t_max    =  1.0;
  double      z_min    = 0.;
  double      z_max    = 1.0;
  double      b        = 0.5;
  int         n_part   = 1;
  int         post     = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_r_seg,
             &n_t_seg,
             &n_z_seg,
             &b,
             &n_part,
             &post,
             (int *) &method);

  double deg2rad = PDM_PI / 180.;
  t_min *= deg2rad;
  t_max *= deg2rad;
  b = PDM_MIN (1., PDM_MAX (0., b));

  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  // if (i_rank == 0) printf("nr = "PDM_FMT_G_NUM", nt = "PDM_FMT_G_NUM", nz = "PDM_FMT_G_NUM"\n", n_r_seg, n_t_seg, n_z_seg);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell      = NULL;
  int          *dface_vtx_idx   = NULL;
  PDM_g_num_t  *dface_vtx       = NULL;
  double       *dvtx_coord      = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group     = NULL;

  _gen_distributed_mesh (PDM_MPI_COMM_WORLD,
                         n_r_seg,
                         n_t_seg,
                         n_z_seg,
                         r_min,
                         r_max,
                         t_min,
                         t_max,
                         z_min,
                         z_max,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &n_face_group,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dface_group_idx,
                         &dface_group,
                         &dvtx_coord);

  if (1) {
    /* Cell size gradation */
    double dr = r_max - r_min;
    double a = 4.*(1. - b);
    for (int i = 0; i < dn_vtx; i++) {
      double x = dvtx_coord[3*i];
      double y = dvtx_coord[3*i+1];
      double z = dvtx_coord[3*i+2] - 0.5;
      double r = sqrt(x*x + y*y);
      double t = (r - r_min) / dr;
      double r2 = r_min + t*t*t*t*dr;
      dvtx_coord[3*i]   = r2 * x / r;
      dvtx_coord[3*i+1] = r2 * y / r;
      dvtx_coord[3*i+2] = z*(b + a*z*z) + 0.5;
      //dvtx_coord[3*i+2] = 0.5 * (1. + pow(2.*z, 3));
    }
  }

  if (1) {
    /* Bump */
    for (int i = 0; i < dn_vtx; i++) {
      double x = dvtx_coord[3*i];
      double y = dvtx_coord[3*i+1];
      double z = dvtx_coord[3*i+2] - 0.5;
      double r = sqrt(x*x + y*y);
      double dr = 0;
      z *= 8;
      if (fabs(z) < 1) {
        dr = exp(1. - 1. / (1. - z*z));
      }
      //dr = exp(-64.*z*z);
      double r2 = r + 0.05 * dr * (r_max - r);
      dvtx_coord[3*i]   = r2 * x / r;
      dvtx_coord[3*i+1] = r2 * y / r;
    }
  }

  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (i_rank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_t *ppart = PDM_part_create (PDM_MPI_COMM_WORLD,
                                       method,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       "PDM_PART_RENUM_FACE_NONE",
                                       n_property_cell,
                                       renum_properties_cell,
                                       n_property_face,
                                       renum_properties_face,
                                       n_part,
                                       dn_cell,
                                       dn_face,
                                       dn_vtx,
                                       n_face_group,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dface_cell,
                                       dface_vtx_idx,
                                       dface_vtx,
                                       NULL,
                                       dvtx_coord,
                                       NULL,
                                       dface_group_idx,
                                       dface_group);

  free (dcell_part);

  free (dface_cell);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dvtx_coord);
  free (dface_group_idx);
  free (dface_group);

  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t *id_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                               n_point_cloud,
                                                               PDM_MPI_COMM_WORLD,
                                                               PDM_OWNERSHIP_KEEP);

  int **select_face   = malloc (sizeof(int *) * n_part);
  int  *n_select_face = malloc (sizeof(int  ) * n_part);
  int **select_vtx    = malloc (sizeof(int *) * n_part);
  int  *n_select_vtx  = malloc (sizeof(int  ) * n_part);

  int **surface_face_vtx_idx =  malloc (sizeof(int *) * n_part);
  int **surface_face_vtx =  malloc (sizeof(int *) * n_part);
  double **surface_coords = malloc (sizeof(double *) * n_part);

  PDM_g_num_t **surface_face_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surface_vtx_parent_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  const PDM_g_num_t **surface_face_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  const PDM_g_num_t **surface_vtx_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);

  PDM_gen_gnum_t *id_gnum_face = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);
  PDM_gen_gnum_t *id_gnum_vtx  = PDM_gnum_create (3, n_part, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_KEEP);

  double **cell_volume = malloc (sizeof(double *) * n_part);
  double **cell_center = malloc (sizeof(double *) * n_part);

  if (i_rank == 0) {
    printf("-- mesh dist set\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

    n_select_face[i_part] = 0;
    n_select_vtx[i_part] = 0;

    select_face[i_part] = malloc (sizeof(int) * n_face);

    for (int i = 0; i < n_face; i++) {
      select_face[i_part][i] = 0;
    }

    select_vtx[i_part] = malloc (sizeof(int) * n_vtx);

    for (int i = 0; i < n_vtx; i++) {
      select_vtx[i_part][i] = 0;
    }

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t  *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t  *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t  *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t  *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
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

    if (post) {
      int ifile = n_part * i_rank + i_part;
      const char *pref = "";///stck/bandrieu/workspace/paradigma-dev/test/axitransbump/";
      char filename[999];
      sprintf(filename, "%scheck_face_cell_%3.3d.vtk", pref, ifile);

      int *face_data = malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++) {
        //int icel2 = face_cell[2*i+1];
        face_data[i] = (int) face_ln_to_gn[i];//icel2;
      }

      _export_vtk_surface (filename,
                           n_face,
                           face_vtx_idx,
                           face_vtx,
                           face_data,
                           n_vtx,
                           vtx);
      free (face_data);

      int *_face_vtx_idx = malloc (sizeof(int) * (face_group_idx[n_face_group] + 1));
      int *_face_vtx = malloc (sizeof(int) * 4 * face_group_idx[n_face_group]);
      face_data = malloc (sizeof(int) * face_group_idx[n_face_group]);
      _face_vtx_idx[0] = 0;
      for (int i = 0; i < n_face_group; i++) {
        for (int j = face_group_idx[i]; j < face_group_idx[i+1]; j++) {
          _face_vtx_idx[j+1] = _face_vtx_idx[j] + 4;
          int k = face_group[j] - 1;
          for (int l = 0; l < 4; l++) {
            _face_vtx[4*j+l] = face_vtx[4*k+l];
          }
          //face_data[j] = (int) face_ln_to_gn[k];
          face_data[j] = i;
        }
      }

      sprintf(filename, "%scheck_face_lim_%3.3d.vtk", pref, ifile);
      _export_vtk_surface (filename,
                           face_group_idx[n_face_group],
                           _face_vtx_idx,
                           _face_vtx,
                           face_data,
                           n_vtx,
                           vtx);
      free (face_data);
      free (_face_vtx_idx);
      free (_face_vtx);
    }

    int iii = 0;
    for (int i = 0; i < n_face; i++) {
      int icel2 = face_cell[2*i+1];
      if (icel2 == 0) {
        iii++;
        select_face[i_part][i] = 1;
        //select_face[i_part][i] = 0;
      }
    }

    /*int i_group = 2;
    for(int idx_face = face_group_idx[i_group]; idx_face < face_group_idx[i_group+1]; ++idx_face){
      int i_face = face_group[idx_face] - 1;
      select_face[i_part][i_face] = 1;
      }*/
    for(int idx_face = face_group_idx[4]; idx_face < face_group_idx[6]; ++idx_face){
      int i_face = face_group[idx_face] - 1;
      select_face[i_part][i_face] = 0;
    }

    for (int i = 0; i < face_part_bound_proc_idx[n_rank]; i++) {
      select_face[i_part][face_part_bound[4*i]-1] = 0;
    }

    int idx = 1;
    int s_face_vtx = 0;
    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] == 1) {
        select_face[i_part][i] = idx;
        s_face_vtx += (face_vtx_idx[i+1] - face_vtx_idx[i]);
        idx += 1;
      }
    }
    n_select_face[i_part] = idx - 1;

    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] != 0) {
        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          select_vtx[i_part][face_vtx[j]-1] = 1;
        }
      }
    }

    idx = 1;
    for (int i = 0; i < n_vtx; i++) {
      if (select_vtx[i_part][i] == 1) {
        select_vtx[i_part][i] = idx;
        idx += 1;
      }
    }
    n_select_vtx[i_part] = idx - 1;

    surface_face_vtx_idx[i_part] = malloc (sizeof(int) * (n_select_face[i_part] + 1));
    surface_face_vtx_idx[i_part][0] = 0;
    surface_face_vtx[i_part] = malloc (sizeof(int) * s_face_vtx);

    surface_coords[i_part] = malloc (sizeof(double) * 3 * n_select_vtx[i_part]);

    surface_face_parent_gnum[i_part] =
      malloc (sizeof(PDM_g_num_t) * n_select_face[i_part]);
    surface_vtx_parent_gnum[i_part] =
      malloc (sizeof(PDM_g_num_t) * n_select_vtx[i_part]);

    surface_face_gnum[i_part] = NULL;
    surface_vtx_gnum[i_part] = NULL;

    idx = 0;
    int idx1 = 0;
    for (int i = 0; i < n_face; i++) {
      if (select_face[i_part][i] > 0) {
        surface_face_vtx_idx[i_part][idx+1] =
          surface_face_vtx_idx[i_part][idx] + (face_vtx_idx[i+1] - face_vtx_idx[i]);

        surface_face_parent_gnum[i_part][idx] = face_ln_to_gn[i];

        idx += 1;

        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          surface_face_vtx[i_part][idx1++] = select_vtx[i_part][face_vtx[j]-1];
        }
      }
    }

    idx = 0;
    for (int i = 0; i < n_vtx; i++) {
      if (select_vtx[i_part][i] > 0) {

        surface_vtx_parent_gnum[i_part][idx] = vtx_ln_to_gn[i];
        surface_coords[i_part][3*idx  ] = vtx[3*i];
        surface_coords[i_part][3*idx+1] = vtx[3*i+1];
        surface_coords[i_part][3*idx+2] = vtx[3*i+2];

        idx += 1;

      }
    }

    PDM_gnum_set_from_parents (id_gnum_face,
                               i_part,
                               n_select_face[i_part],
                               surface_face_parent_gnum[i_part]);


    PDM_gnum_set_from_parents (id_gnum_vtx,
                               i_part,
                               n_select_vtx[i_part],
                               surface_vtx_parent_gnum[i_part]);

    for (int i = 0; i <  n_select_vtx[i_part]; i++) {

    }

  }

  PDM_gnum_compute (id_gnum_face);

  PDM_gnum_compute (id_gnum_vtx);

  PDM_g_num_t n_g_face_loc = 0;
  PDM_g_num_t n_g_vtx_loc = 0;

  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_vtx = 0;

  for (int i_part = 0; i_part < n_part; i_part++) {
    surface_face_gnum[i_part] = PDM_gnum_get (id_gnum_face, i_part);
    surface_vtx_gnum[i_part] = PDM_gnum_get (id_gnum_vtx, i_part);

    for (int i = 0; i <  n_select_face[i_part]; i++) {
      n_g_face_loc = PDM_MAX(n_g_face_loc, surface_face_gnum[i_part][i]);
    }

    for (int i = 0; i <  n_select_vtx[i_part]; i++) {
      n_g_vtx_loc = PDM_MAX(n_g_vtx_loc, surface_vtx_gnum[i_part][i]);
    }

  }

  PDM_MPI_Allreduce (&n_g_face_loc, &n_g_face, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_MPI_Allreduce (&n_g_vtx_loc, &n_g_vtx, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (id_dist,
                                           n_part);

  PDM_dist_cloud_surf_n_part_cloud_set (id_dist, 0, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    PDM_dist_cloud_surf_surf_mesh_part_set (id_dist,
                                            i_part,
                                            n_select_face[i_part],
                                            surface_face_vtx_idx[i_part],
                                            surface_face_vtx[i_part],
                                            surface_face_gnum[i_part],
                                            n_select_vtx[i_part],
                                            surface_coords[i_part],
                                            surface_vtx_gnum[i_part]);

    int ifile = n_part * i_rank + i_part;
    char filename[999];
    sprintf(filename, "dist_axi_surf_%3.3d.vtk", ifile);
    if (post) {
      //int *face_data = PDM_array_const_int (n_select_face[i_part], i_rank);
      int *face_data = malloc (sizeof(int) * n_select_face[i_part]);
      for (int i = 0; i < n_select_face[i_part]; i++) {
        face_data[i] = (int) surface_face_gnum[i_part][i];
      }
      _export_vtk_surface (filename,
                           n_select_face[i_part],
                           surface_face_vtx_idx[i_part],
                           surface_face_vtx[i_part],
                           face_data,//NULL,
                           n_select_vtx[i_part],
                           surface_coords[i_part]);
      free (face_data);
    }

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

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

    PDM_part_part_val_get (ppart,
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

    const int     isOriented = 0;
    cell_volume[i_part] = malloc(sizeof(double) * n_cell);
    cell_center[i_part] = malloc(sizeof(double) * 3 * n_cell);

    PDM_geom_elem_polyhedra_properties (isOriented,
                                        n_cell,
                                        n_face,
                                       face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume[i_part],
                                        cell_center[i_part],
                                        NULL,
                                        NULL);

    PDM_dist_cloud_surf_cloud_set (id_dist,
                                   0,
                                   i_part,
                                   n_cell,
                                   cell_center[i_part],
                                   cell_ln_to_gn);

    sprintf(filename, "dist_axi_cloud_%3.3d.vtk", ifile);

    if (post) {
      _export_vtk_cloud (filename,
                         n_cell,
                       cell_center[i_part]);
    }
  }

  if (i_rank == 0) {
    printf("-- Dist compute\n");
    fflush(stdout);
  }

  PDM_dist_cloud_surf_compute (id_dist);

  if (i_rank == 0) {
    printf("-- Dist check\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    double      *distance;
    double      *projected;
    PDM_g_num_t *closest_elt_gnum;

    PDM_dist_cloud_surf_get (id_dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_total_part,
                           &scell_face,
                           &sface_vtx,
                           &sface_group,
                           &nEdgeGroup2);

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

    PDM_part_part_val_get (ppart,
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




  /*
   *  Brute force check
   */
  //-->>
  if (1) {
    const double eps_distance = 1.e-6;
    double _vtx_coord[12];

    /* Allgather surface mesh */
    int tn_face = 0;
    int tn_vtx = 0;
    int tl_face_vtx = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      tn_face += n_select_face[i_part];
      tn_vtx += n_select_vtx[i_part];
      tl_face_vtx += surface_face_vtx_idx[i_part][n_select_face[i_part]];
    }

    int *tface_vtx_idx = malloc (sizeof(int) * (tn_face+1));
    int *tface_vtx = malloc (sizeof(int) * tl_face_vtx);
    PDM_g_num_t *tface_g_num = malloc (sizeof(PDM_g_num_t) * tn_face);
    double *tvtx_coord = malloc (sizeof(double) * tn_vtx * 3);
    tn_face = 0;
    tn_vtx = 0;
    tl_face_vtx = 0;
    tface_vtx_idx[0] = 0;
    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_select_face[i_part]; i++) {
        tface_vtx_idx[tn_face+1] = tface_vtx_idx[tn_face];
        for (int j = surface_face_vtx_idx[i_part][i];
             j < surface_face_vtx_idx[i_part][i+1]; j++) {
          tface_vtx[tface_vtx_idx[tn_face+1]++] = tn_vtx + surface_face_vtx[i_part][j];
        }
        tface_g_num[tn_face] = surface_face_gnum[i_part][i];
        tn_face++;
      }
      for (int i = 0; i < n_select_vtx[i_part]; i++) {
        for (int j = 0; j < 3; j++) {
          tvtx_coord[3*tn_vtx + j] = surface_coords[i_part][3*i + j];
        }
        tn_vtx++;
      }
    }

    if (0) {
      int *tface_data = malloc (sizeof(int) * tn_face);
      for (int i = 0; i < tn_face; i++) {
        tface_data[i] = (int) tface_g_num[i];
      }
      /* int idx = 0;
      for (int i_part = 0; i_part < n_part; i_part++) {
        for (int i = 0; i < n_select_face[i_part]; i++) {
          tface_data[idx++] = i_part;
        }
        }*/

      char filename[999];
      sprintf(filename, "tsurf_mesh_%2.2d.vtk", i_rank);

      _export_vtk_surface (filename,
                           tn_face,
                           tface_vtx_idx,
                           tface_vtx,
                           tface_data,
                           tn_vtx/3,
                           tvtx_coord);
      free (tface_data);
    }

    int *all_tn_face = malloc (sizeof(int) * n_rank);
    PDM_MPI_Allgather (&tn_face,    1, PDM_MPI_INT,
                       all_tn_face, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    tn_vtx *= 3;
    int *all_tn_vtx = malloc (sizeof(int) * n_rank);
    PDM_MPI_Allgather (&tn_vtx,    1, PDM_MPI_INT,
                       all_tn_vtx, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *all_tl_face_vtx = malloc (sizeof(int) * n_rank);
    PDM_MPI_Allgather (&tl_face_vtx,    1, PDM_MPI_INT,
                       all_tl_face_vtx, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *shift_face = PDM_array_new_idx_from_sizes_int (all_tn_face, n_rank);
    int *shift_vtx  = PDM_array_new_idx_from_sizes_int (all_tn_vtx, n_rank);
    int *shift_face_vtx = PDM_array_new_idx_from_sizes_int (all_tl_face_vtx, n_rank);

    int gn_face = shift_face[n_rank];
    int gn_vtx = shift_vtx[n_rank];

    PDM_g_num_t *gface_g_num = malloc (sizeof(PDM_g_num_t) * gn_face);
    PDM_MPI_Allgatherv (tface_g_num, tn_face, PDM__PDM_MPI_G_NUM,
                        gface_g_num, all_tn_face, shift_face, PDM__PDM_MPI_G_NUM,
                        PDM_MPI_COMM_WORLD);
    free (tface_g_num);

    double *gvtx_coord = malloc (sizeof(double) * gn_vtx);
    PDM_MPI_Allgatherv (tvtx_coord, tn_vtx, PDM_MPI_DOUBLE,
                        gvtx_coord, all_tn_vtx, shift_vtx, PDM_MPI_DOUBLE,
                        PDM_MPI_COMM_WORLD);
    free (tvtx_coord);
    gn_vtx /= 3;

    int *gface_vtx = malloc (sizeof(int) * gn_face * 4);
    for (int i = 0; i < n_rank; i++) {
      all_tn_face[i] *= 4;
      shift_face[i+1] *= 4;
      shift_vtx[i+1] /= 3;
    }
    PDM_MPI_Allgatherv (tface_vtx, 4*tn_face, PDM_MPI_INT,
                        gface_vtx, all_tn_face, shift_face, PDM_MPI_INT,
                        PDM_MPI_COMM_WORLD);
    free (tface_vtx);
    free (tface_vtx_idx);

    for (int i = 0; i < n_rank; i++) {
      for (int j = shift_face[i]; j < shift_face[i+1]; j++) {
        gface_vtx[j] += shift_vtx[i];
      }
    }

    for (int i = 0; i < n_rank; i++) {
      shift_face[i+1] /= 4;
    }

    int *gface_vtx_idx = malloc (sizeof(int) * (gn_face + 1));
    gface_vtx_idx[0] = 0;
    for (int i = 0; i < gn_face; i++) {
      gface_vtx_idx[i+1] = gface_vtx_idx[i] + 4;
    }

    if (post) {
      int *gface_data = malloc (sizeof(int) * gn_face);
      for (int i = 0; i < gn_face; i++) {
        gface_data[i] = (int) gface_g_num[i];
      }
      char filename[999];
      sprintf(filename, "gsurf_mesh_%2.2d.vtk", i_rank);

      _export_vtk_surface (filename,
                           gn_face,
                           gface_vtx_idx,
                           gface_vtx,
                           gface_data,
                           gn_vtx,
                           gvtx_coord);
      free (gface_data);
    }

    if (1) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        double      *distance;
        double      *projected;
        PDM_g_num_t *closest_elt_gnum;
        PDM_dist_cloud_surf_get (id_dist,
                                 0,
                                 i_part,
                                 &distance,
                                 &projected,
                                 &closest_elt_gnum);

        int n_cell;
        int n_face;
        int n_face_part_bound;
        int n_vtx;
        int n_proc;
        int n_total_part;
        int scell_face;
        int sface_vtx;
        int sface_group;
        int nEdgeGroup2;

        PDM_part_part_dim_get (ppart,
                               i_part,
                               &n_cell,
                               &n_face,
                               &n_face_part_bound,
                               &n_vtx,
                               &n_proc,
                               &n_total_part,
                               &scell_face,
                               &sface_vtx,
                               &sface_group,
                               &nEdgeGroup2);

        int          *cell_tag;
        int          *cell_face_idx;
        int          *cell_face;
        PDM_g_num_t  *cell_ln_to_gn;
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

        PDM_part_part_val_get (ppart,
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

        for (int itgt = 0; itgt < n_cell; itgt++) {
          double dist_min = HUGE_VAL;
          PDM_g_num_t arg_min = 0;
          double proj_min[3] = {0};

          for (int rank = 0; rank < n_rank; rank++) {
            for (int i = shift_face[rank]; i < shift_face[rank+1]; i++) {
              int elt_vtx_n = gface_vtx_idx[i+1] - gface_vtx_idx[i];
              for (int j = 0; j < elt_vtx_n; j++) {
                int ivtx = gface_vtx[gface_vtx_idx[i] + j] - 1;
                for (int k = 0; k < 3; k++) {
                  _vtx_coord[3*j + k] = gvtx_coord[3*ivtx + k];
                }
              }

              double dist;
              double proj[3];
              PDM_polygon_evaluate_position (cell_center[i_part] + 3*itgt,
                                             elt_vtx_n,
                                             _vtx_coord,
                                             proj,
                                             &dist);
              if (dist_min > dist ||
                  (PDM_ABS(dist_min-dist) < 1e-12 && gface_g_num[i] < arg_min)) {
                dist_min = dist;
                arg_min = gface_g_num[i];
                for (int k = 0; k < 3; k++) {
                  proj_min[k] = proj[k];
                }
              }
            }
          }

          // Check
          double d_res = sqrt(distance[itgt]);
          double d_exp = sqrt(dist_min);
          if (closest_elt_gnum[itgt] != arg_min ||
              PDM_ABS(d_res - d_exp) > eps_distance) {
            printf("[%d] part %d, pt "PDM_FMT_G_NUM" : expected ("PDM_FMT_G_NUM", %f), result ("PDM_FMT_G_NUM", %f), relative error = %g\n", i_rank, i_part, cell_ln_to_gn[itgt], arg_min, d_exp, closest_elt_gnum[itgt], d_res, (d_res - d_exp)/d_res);
          }

          else {
            double e_proj = 0.;
            for (int k = 0; k < 3; k++) {
              double delta = proj_min[k] - projected[3*itgt + k];
              e_proj += delta * delta;
            }
            e_proj = sqrt(e_proj);
            if (e_proj > eps_distance) {
              printf("[%d] part %d, pt "PDM_FMT_G_NUM", "PDM_FMT_G_NUM", proj : expected (%f %f %f), result (%f %f %f), error = %f\n", i_rank, i_part, cell_ln_to_gn[itgt], closest_elt_gnum[itgt], proj_min[0], proj_min[1], proj_min[2], projected[3*itgt], projected[3*itgt + 1], projected[3*itgt + 2], e_proj);
            }
          }
        }
      }
    }

    free(gface_vtx);
    free(gface_vtx_idx);

    free (gface_g_num);
    free (gvtx_coord);
    free (shift_face);
    free (shift_vtx);
    free (shift_face_vtx);
    free (all_tn_face);
    free (all_tn_vtx);
    free (all_tl_face_vtx);
  }
  //<<--

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (cell_center[i_part]);
    free (cell_volume[i_part]);
  }
  free (cell_center);
  free (cell_volume);

  PDM_part_free(ppart);

  PDM_dist_cloud_surf_dump_times(id_dist);
  PDM_dist_cloud_surf_free (id_dist);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free (select_face[i_part]);
    free (select_vtx[i_part]);

    free (surface_face_vtx_idx[i_part]);
    free (surface_face_vtx[i_part]);
    free (surface_coords[i_part]);

    free (surface_face_parent_gnum[i_part]);
    free (surface_vtx_parent_gnum[i_part]);

  }

  free (select_face);
  free (select_vtx);

  free (n_select_face);
  free (n_select_vtx);

  free (surface_face_vtx_idx);
  free (surface_face_vtx);
  free (surface_coords);

  free (surface_face_parent_gnum);
  free (surface_vtx_parent_gnum);

  free (surface_face_gnum);
  free (surface_vtx_gnum);

  PDM_gnum_free(id_gnum_face);
  PDM_gnum_free(id_gnum_vtx);

  PDM_MPI_Finalize();

   if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
