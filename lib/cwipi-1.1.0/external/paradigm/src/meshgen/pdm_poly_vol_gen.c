/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_poly_vol_gen.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"


/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static double _rand (void)
{
  return 2*((double) rand() / (double) RAND_MAX) - 1;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Generate a distributed polyhedral mesh
 *
 * \param [in]   pdm_comm         MPI communicator
 * \param [in]   xmin             Minimum x-coordinate
 * \param [in]   ymin             Minimum y-coordinate
 * \param [in]   zmin             Minimum z-coordinate
 * \param [in]   lengthx          Length in the x-direction
 * \param [in]   lengthy          Length in the y-direction
 * \param [in]   lengthz          Length in the z-direction
 * \param [in]   nx               Number of vertices in the x-direction
 * \param [in]   ny               Number of vertices in the y-direction
 * \param [in]   nz               Number of vertices in the z-direction
 * \param [in]   randomize        Enable/disable randomization
 * \param [in]   random_seed      Random seed
 * \param [out]  ng_cell          Global number of cells
 * \param [out]  ng_face          Global number of faces
 * \param [out]  ng_vtx           Global number of vertices
 * \param [out]  n_face_group     Number of face groups
 * \param [out]  dn_cell          Local number of cells
 * \param [out]  dn_face          Local number of faces
 * \param [out]  dn_vtx           Local number of vertices
 * \param [out]  dcell_face_idx   Index of cell-face connectivity (size = \ref dn_cell + 1)
 * \param [out]  dcell_face       Distributed cell-face connectivity (size = \ref dcell_face_idx[\ref dn_cell])
 * \param [out]  dface_cell       Distributed face-cell connectivity (size = 2 * \ref dn_face)
 * \param [out]  dface_vtx_idx    Index of face-vertex connectivity (size = \ref dn_face + 1)
 * \param [out]  dface_vtx        Distributed face-vertex connectivity (size = \ref dface_vtx_idx[\ref dn_face])
 * \param [out]  dvtx_coord       Coordinates of local vertices (size = 3 * \ref dn_vtx)
 * \param [out]  dface_group_idx  Index of dface_group (size = \ref n_face_group + 1)
 * \param [out]  dface_group      Distributed lists of faces in each group (size = \ref dface_group_idx[\ref n_face_group])
 *
 */

void
PDM_poly_vol_gen
(
 PDM_MPI_Comm  comm,
 double        xmin,
 double        ymin,
 double        zmin,
 double        lengthx,
 double        lengthy,
 double        lengthz,
 PDM_g_num_t   nx,
 PDM_g_num_t   ny,
 PDM_g_num_t   nz,
 int           randomize,
 int           random_seed,
 PDM_g_num_t  *ng_cell,
 PDM_g_num_t  *ng_face,
 PDM_g_num_t  *ng_vtx,
 int          *n_face_group,
 int          *dn_cell,
 int          *dn_face,
 int          *dn_vtx,
 int         **dcell_face_idx,
 PDM_g_num_t **dcell_face,
 PDM_g_num_t **dface_cell,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 double      **dvtx_coord,
 int         **dface_group_idx,
 PDM_g_num_t **dface_group
 )
{
  srand (random_seed);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t n_vtx1 = 2*nx;
  PDM_g_num_t n_vtx2 = nx + 1;
  PDM_g_num_t n_vtx3 = n_vtx1 + 2*n_vtx2;
  PDM_g_num_t n_vtx_z_cst = (ny+1)*2*nx + (nx+1)*2*ny + 4;


  PDM_g_num_t n_octo_z_cst = nx * ny;
  PDM_g_num_t n_quadH_z_cst = (nx - 1) * (ny - 1);
  PDM_g_num_t n_tria_z_cst = 2*(nx - 1) + 2*(ny - 1) + 4;
  PDM_g_num_t n_quadV_z_cst = nx*(ny+1) + ny*(nx+1) + 4*nx*ny + 2*(nx-1) + 2*(ny-1) + 8;
  PDM_g_num_t n_faceH_z_cst = n_octo_z_cst + n_quadH_z_cst + n_tria_z_cst;
  PDM_g_num_t n_face_z_cst = n_faceH_z_cst + n_quadV_z_cst;
  PDM_g_num_t ng_face_lim = 2*(n_faceH_z_cst + nz*(2*nx + 1 + 2*ny + 1));

  *ng_vtx = n_vtx_z_cst * (nz + 1);
  *ng_face = nz*n_face_z_cst + n_faceH_z_cst;
  *ng_cell = n_faceH_z_cst * nz;
  *n_face_group = 6;

  /* Define distributions */
  PDM_g_num_t *distrib_vtx  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_cell = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face_lim = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  distrib_vtx[0]      = 0;
  distrib_face[0]     = 0;
  distrib_cell[0]     = 0;
  distrib_face_lim[0] = 0;

  PDM_g_num_t step_vtx       = *ng_vtx / n_rank;
  PDM_g_num_t remainder_vtx  = *ng_vtx % n_rank;

  PDM_g_num_t step_face      = *ng_face / n_rank;
  PDM_g_num_t remainder_face = *ng_face % n_rank;

  PDM_g_num_t step_cell      = *ng_cell / n_rank;
  PDM_g_num_t remainder_cell = *ng_cell % n_rank;

  PDM_g_num_t step_face_lim      = ng_face_lim / n_rank;
  PDM_g_num_t remainder_face_lim = ng_face_lim % n_rank;

  for (int i = 0; i < n_rank; i++) {
    distrib_vtx[i+1] = distrib_vtx[i] + step_vtx;
    if (i < remainder_vtx) {
      distrib_vtx[i+1]++;
    }

    distrib_face[i+1] = distrib_face[i] + step_face;
    if (i < remainder_face) {
      distrib_face[i+1]++;
    }

    distrib_cell[i+1] = distrib_cell[i] + step_cell;
    if (i < remainder_cell) {
      distrib_cell[i+1]++;
    }

    distrib_face_lim[i+1] = distrib_face_lim[i] + step_face_lim;
    if (i < remainder_face_lim) {
      distrib_face_lim[i+1]++;
    }
  }
  *dn_vtx  = (int) (distrib_vtx[i_rank+1]  - distrib_vtx[i_rank]);
  *dn_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);
  *dn_cell = (int) (distrib_cell[i_rank+1] - distrib_cell[i_rank]);
  int dn_face_lim = (int) (distrib_face_lim[i_rank+1] - distrib_face_lim[i_rank]);

  if (0) {
    printf("[%d] dn_cell = %d, dn_face = %d, dn_vtx = %d\n", i_rank, *dn_cell, *dn_face, *dn_vtx);
  }

  /*
   *  Vertices
   */
  *dvtx_coord = malloc (sizeof(double) * (*dn_vtx) * 3);
  double *_dvtx_coord = *dvtx_coord;

  double stepx = lengthx / (double) (3*nx);
  double stepy = lengthy / (double) (3*ny);
  double stepz = 0.;
  if (nz > 0) {
    stepz = lengthz / (double) nz;
  }

  double noise = 0.4 * PDM_MIN (stepx, PDM_MIN (stepy, stepz / 3.));

  for (int ivtx = 0; ivtx < *dn_vtx; ivtx++) {
    PDM_g_num_t g = distrib_vtx[i_rank] + ivtx;
    PDM_g_num_t k = g / n_vtx_z_cst;
    PDM_g_num_t r = g % n_vtx_z_cst;
    PDM_g_num_t i, j;
    double x, y, z;
    z = zmin + k*stepz;

    srand (random_seed + g);

    double rx = noise * _rand();
    double ry = noise * _rand();
    double rz = noise * _rand()/3.;

    if (k == 0 || k == nz) {
      rz = 0.;
    }

    if (r > n_vtx_z_cst - 5) {
      // Corner
      r -= (n_vtx_z_cst - 4);
      j = r / 2;
      i = r % 2;

      x = xmin + i*lengthx;
      y = ymin + j*lengthy;
      rx = 0.;
      ry = 0.;
    }
    else {
      PDM_g_num_t jj = r / n_vtx3;
      PDM_g_num_t s = r % n_vtx3;

      if (s < n_vtx1) {
        // row n_vtx1
        i = 3*(s/2) + 1 + s%2;
        j = 3*jj;
        if (jj == 0 || jj == ny) {
          ry = 0.;
        }
      }
      else {
        // row n_vtx2
        s -= n_vtx1;
        j = 3*jj + s / n_vtx2 + 1;
        i = 3*(s % n_vtx2);
        if (i == 0 || i == 3*nx) {
          rx = 0;
        }
      }
      x = xmin + i*stepx;
      y = ymin + j*stepy;
    }

    _dvtx_coord[3*ivtx    ] = x;
    _dvtx_coord[3*ivtx + 1] = y;
    _dvtx_coord[3*ivtx + 2] = z;

    if (randomize) {
      _dvtx_coord[3*ivtx    ] += rx;
      _dvtx_coord[3*ivtx + 1] += ry;
      _dvtx_coord[3*ivtx + 2] += rz;
    }
  }


  /*
   *  Faces
   */
  PDM_g_num_t idx_quadH = n_octo_z_cst;
  PDM_g_num_t idx_tria1 = idx_quadH + n_quadH_z_cst; // bottom row
  PDM_g_num_t idx_tria2 = idx_tria1 + nx - 1;        // top row
  PDM_g_num_t idx_tria3 = idx_tria2 + nx - 1;        // left column
  PDM_g_num_t idx_tria4 = idx_tria3 + ny - 1;        // right column
  PDM_g_num_t idx_tria5 = idx_tria4 + ny - 1;        // corners
  PDM_g_num_t idx_quadV1 = n_faceH_z_cst;            // _
  PDM_g_num_t idx_quadV2 = idx_quadV1 + nx*(ny + 1); // |
  PDM_g_num_t idx_quadV3 = idx_quadV2 + ny*(nx + 1); // \ (pointing +y)
  PDM_g_num_t idx_quadV4 = idx_quadV3 + nx*ny;       // / (pointing +y)
  PDM_g_num_t idx_quadV5 = idx_quadV4 + nx*ny;       // / (pointing -y)
  PDM_g_num_t idx_quadV6 = idx_quadV5 + nx*ny;       // \ (pointing -y)
  PDM_g_num_t idx_quadV7 = idx_quadV6 + nx*ny;       // bottom row
  PDM_g_num_t idx_quadV8 = idx_quadV7 + nx - 1;      // top row
  PDM_g_num_t idx_quadV9 = idx_quadV8 + nx - 1;      // left column
  PDM_g_num_t idx_quadV10 = idx_quadV9 + ny - 1;     // right column
  PDM_g_num_t idx_quadV11 = idx_quadV10 + ny - 1;    // corners


  *dface_vtx_idx = malloc (sizeof(int) * (*dn_face + 1));
  int *_dface_vtx_idx = *dface_vtx_idx;
  _dface_vtx_idx[0] = 0;

  int s_face_vtx = 8*(*dn_face);
  *dface_vtx = malloc (sizeof(PDM_g_num_t) * s_face_vtx);

  *dface_cell = malloc (sizeof(PDM_g_num_t) * 2 * (*dn_face));

  for (int ifac = 0; ifac < (*dn_face); ifac++) {

    PDM_g_num_t *_dface_vtx = *dface_vtx + _dface_vtx_idx[ifac];
    PDM_g_num_t *_dface_cell = *dface_cell + 2*ifac;

    PDM_g_num_t g = distrib_face[i_rank] + ifac;
    PDM_g_num_t k = g / n_face_z_cst;
    PDM_g_num_t r = g % n_face_z_cst;
    PDM_g_num_t i, j;

    /* Horizontal faces */
    if (r < idx_quadH) {
      // Octagon
      j = r / nx;
      i = r % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + 2*i;
      _dface_vtx[1] = idxv + 2*i + 1;
      _dface_vtx[2] = idxv + n_vtx1 + i + 1;
      _dface_vtx[3] = idxv + n_vtx1 + n_vtx2 + i + 1;
      _dface_vtx[4] = idxv + n_vtx3 + 2*i + 1;
      _dface_vtx[5] = idxv + n_vtx3 + 2*i;
      _dface_vtx[6] = idxv + n_vtx1 + n_vtx2 + i;
      _dface_vtx[7] = idxv + n_vtx1 + i;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 8;

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria1) {
      // Quadrangle
      PDM_g_num_t s = r - idx_quadH;

      j = s / (nx - 1);
      i = s % (nx - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + i + 1;
      _dface_vtx[1] = idxv + n_vtx2 + 2*(i+1);
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx2 + i + 1;
      _dface_vtx[3] = idxv + n_vtx2 + 2*i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_quadH + j*(nx - 1) + i + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria2) {
      // Triangle (bottom row)
      PDM_g_num_t s = r - idx_tria1;

      i = s % (nx - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + 2*i + 1;
      _dface_vtx[1] = idxv + 2*(i+1);
      _dface_vtx[2] = idxv + n_vtx1 + i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria1 + i + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria3) {
      // Triangle (top row)
      r -= idx_tria2;

      i = r % (nx - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + ny*n_vtx3 - n_vtx2 + 1;
      _dface_vtx[0] = idxv + n_vtx2 + 2*i + 1;
      _dface_vtx[1] = idxv + i + 1;
      _dface_vtx[2] = idxv + n_vtx2 + 2*(i+1);
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria2 + i + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria4) {
      // Triangle (left column)
      PDM_g_num_t s = r - idx_tria3;

      j = s % (ny - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx2;
      _dface_vtx[2] = idxv + n_vtx2 + n_vtx1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria3 + j + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria5) {
      // Triangle (right column)
      r -= idx_tria4;

      j = r % (ny - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + nx;
      _dface_vtx[1] = idxv + n_vtx2 + n_vtx1 + nx;
      _dface_vtx[2] = idxv + n_vtx2 + 2*(nx - 1) + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria4 + j + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria5 + 1) {
      // Triangle (bottom-left corner)
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx_z_cst - 4;
      _dface_vtx[1] = idxv;
      _dface_vtx[2] = idxv + n_vtx1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria5 + 1;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria5 + 2) {
      // Triangle (bottom-right corner)
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx1 - 1;
      _dface_vtx[1] = idxv + n_vtx_z_cst - 3;
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx2 - 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria5 + 2;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria5 + 3) {
      // Triangle (top-left corner)
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + (ny-1)*n_vtx3 + n_vtx1 + n_vtx2;
      _dface_vtx[1] = idxv + ny*n_vtx3;
      _dface_vtx[2] = idxv + n_vtx_z_cst - 2;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria5 + 3;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    else if (r < idx_tria5 + 4) {
      // Triangle (top-right corner)
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + (ny-1)*n_vtx3 + n_vtx1 + 2*n_vtx2 - 1;
      _dface_vtx[1] = idxv + n_vtx_z_cst - 1;
      _dface_vtx[2] = idxv + n_vtx_z_cst - 5;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;

      PDM_g_num_t idxc = k*n_faceH_z_cst + idx_tria5 + 4;
      if (k == nz) {
        _dface_cell[0] = idxc - n_faceH_z_cst;
      } else {
        _dface_cell[0] = idxc;
      }
      if (k == 0 || k == nz) {
        _dface_cell[1] = 0;
      } else {
        _dface_cell[1] = idxc - n_faceH_z_cst;
      }
    }

    /* Vertical faces */
    else if (r < idx_quadV2) {
      // Quadrangle _
      PDM_g_num_t s = r - idx_quadV1;
      j = s / nx;
      i = s % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 2*i + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx_z_cst;
      _dface_vtx[2] = idxv + 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      if (j == ny) {
        // flip
        for (int l = 0; l < 2; l++) {
          PDM_g_num_t tmp = _dface_vtx[l];
          _dface_vtx[l] = _dface_vtx[3 - l];
          _dface_vtx[3 - l] = tmp;
        }
      }

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      if (j == 0) {
        _dface_cell[0] = idxc;
        _dface_cell[1] = 0;
      } else if (j == ny) {
        _dface_cell[0] = idxc - nx;
        _dface_cell[1] = 0;
      } else {
        _dface_cell[0] = idxc;
        _dface_cell[1] = idxc - nx;
      }
    }

    else if (r < idx_quadV3) {
      // Quadrangle |
      PDM_g_num_t s = r - idx_quadV2;
      j = s / (nx + 1);
      i = s % (nx + 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + n_vtx1 + i + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx2;
      _dface_vtx[2] = idxv + n_vtx2 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      if (i == nx) {
        // flip
        for (int l = 0; l < 2; l++) {
          PDM_g_num_t tmp = _dface_vtx[l];
          _dface_vtx[l] = _dface_vtx[3 - l];
          _dface_vtx[3 - l] = tmp;
        }
      }

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      if (i == 0) {
        _dface_cell[0] = idxc;
        _dface_cell[1] = 0;
      } else if (i == nx) {
        _dface_cell[0] = idxc - 1;
        _dface_cell[1] = 0;
      } else {
        _dface_cell[0] = idxc;
        _dface_cell[1] = idxc - 1;
      }
    }

    else if (r < idx_quadV4) {
      // Quadrangle \ (pointing +y)
      PDM_g_num_t s = r - idx_quadV3;
      j = s / nx;
      i = s % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + 2*i;
      _dface_vtx[1] = idxv + n_vtx1 + i;
      _dface_vtx[2] = idxv + n_vtx1 + i + n_vtx_z_cst;
      _dface_vtx[3] = idxv + 2*i + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      _dface_cell[0] = idxc;
      if (j == 0) {
        if (i == 0) {
          // bottom-left corner
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria5 + 1;
        } else {
          // bottom row
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria1 + i;
        }
      } else {
        if (i == 0) {
          // left column
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria3 + j;
        } else {
          _dface_cell[1] = k*n_faceH_z_cst + idx_quadH + (j-1)*(nx - 1) + i;
        }
      }
    }

    else if (r < idx_quadV5) {
      // Quadrangle / (pointing +y)
      PDM_g_num_t s = r - idx_quadV4;
      j = s / nx;
      i = s % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + 2*i + 1;
      _dface_vtx[1] = idxv + 2*i + 1 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx1 + i + 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx1 + i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      _dface_cell[0] = idxc;
      if (j == 0) {
        if (i == nx-1) {
          // bottom-right corner
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria5 + 2;
        } else {
          // bottom row
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria1 + i + 1;
        }
      } else {
        if (i == nx-1) {
          // right column
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria4 + j;
        } else {
          _dface_cell[1] = k*n_faceH_z_cst + idx_quadH + (j-1)*(nx - 1) + i + 1;
        }
      }
    }

    else if (r < idx_quadV6) {
      // Quadrangle / (pointing -y)
      PDM_g_num_t s = r - idx_quadV5;
      j = s / nx;
      i = s % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + n_vtx1 + n_vtx2 + i;
      _dface_vtx[1] = idxv + n_vtx3 + 2*i;
      _dface_vtx[2] = idxv + n_vtx3 + 2*i + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx1 + n_vtx2 + i + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      _dface_cell[0] = idxc;
      if (j == ny-1) {
        if (i == 0) {
          // top-left corner
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria5 + 3;
        } else {
          // top row
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria2 + i;
        }
      } else {
        if (i == 0) {
          // left column
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria3 + j + 1;
        } else {
          _dface_cell[1] = k*n_faceH_z_cst + idx_quadH + j*(nx - 1) + i;
        }
      }
    }

    else if (r < idx_quadV7) {
      // Quadrangle \ (pointing -y)
      PDM_g_num_t s = r - idx_quadV6;
      j = s / nx;
      i = s % nx;

      PDM_g_num_t idxv = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv + n_vtx1 + n_vtx2 + i + 1;
      _dface_vtx[1] = idxv + n_vtx1 + n_vtx2 + i + 1 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx3 + 2*i + 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx3 + 2*i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      PDM_g_num_t idxc = k*n_faceH_z_cst + j*nx + i + 1;
      _dface_cell[0] = idxc;
      if (j == ny-1) {
        if (i == nx-1) {
          // top-right corner
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria5 + 4;
        } else {
          // top row
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria2 + i + 1;
        }
      } else {
        if (i == nx-1) {
          // right column
          _dface_cell[1] = k*n_faceH_z_cst + idx_tria4 + j + 1;
        } else {
          _dface_cell[1] = k*n_faceH_z_cst + idx_quadH + j*(nx - 1) + i + 1;
        }
      }
    }

    else if (r < idx_quadV8) {
      // Quadrangle (bottom row)
      PDM_g_num_t s = r - idx_quadV7;
      i = s % (nx - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + 2*i + 1 + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx_z_cst;
      _dface_vtx[2] = idxv + 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria1 + i + 1;
      _dface_cell[1] = 0;
    }

    else if (r < idx_quadV9) {
      // Quadrangle (top row)
      PDM_g_num_t s = r - idx_quadV8;
      i = s % (nx - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + ny*n_vtx3 + 2*i + 1 + 1;
      _dface_vtx[0] = idxv + 1;
      _dface_vtx[1] = idxv + 1 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx_z_cst;
      _dface_vtx[3] = idxv;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria2 + i + 1;
      _dface_cell[1] = 0;
    }

    else if (r < idx_quadV10) {
      // Quadrangle (left column)
      PDM_g_num_t s = r - idx_quadV9;
      j = s % (ny - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx1 + n_vtx2;
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx2+ n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria3 + j + 1;
      _dface_cell[1] = 0;
    }

    else if (r < idx_quadV11) {
      // Quadrangle (right column)
      PDM_g_num_t s = r - idx_quadV10;
      j = s % (ny - 1);

      PDM_g_num_t idxv = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + nx + 1;
      _dface_vtx[0] = idxv;
      _dface_vtx[1] = idxv + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx2+ n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx1 + n_vtx2;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria4 + j + 1;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 7) {
      // Quadrangle bottom-left corner _
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx_z_cst - 4;
      _dface_vtx[1] = idxv + 2*n_vtx_z_cst - 4;
      _dface_vtx[2] = idxv + n_vtx_z_cst;
      _dface_vtx[3] = idxv;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 1;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 6) {
      // Quadrangle bottom-left corner |
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx_z_cst - 4;
      _dface_vtx[1] = idxv + n_vtx1;
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + 2*n_vtx_z_cst - 4;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 1;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 5) {
      // Quadrangle bottom-right corner _
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx1 - 1;
      _dface_vtx[1] = idxv + n_vtx1 - 1 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx_z_cst - 3 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx_z_cst - 3;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 2;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 4) {
      // Quadrangle bottom-right corner |
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx_z_cst - 3;
      _dface_vtx[1] = idxv + n_vtx_z_cst - 3 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx1 + n_vtx2 - 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx1 + n_vtx2 - 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 2;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 3) {
      // Quadrangle top-left corner |
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + (ny-1)*n_vtx3 + n_vtx1 + n_vtx2;
      _dface_vtx[1] = idxv + n_vtx_z_cst - 2;
      _dface_vtx[2] = idxv + n_vtx_z_cst - 2 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + (ny-1)*n_vtx3 + n_vtx1 + n_vtx2 + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 3;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 2) {
      // Quadrangle top-left corner _
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + n_vtx_z_cst - 2;
      _dface_vtx[1] = idxv + ny*n_vtx3;
      _dface_vtx[2] = idxv + ny*n_vtx3 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx_z_cst - 2 + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 3;
      _dface_cell[1] = 0;
    }

    else if (r < n_face_z_cst - 1) {
      // Quadrangle top-right corner |
      PDM_g_num_t idxv = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv + (ny-1)*n_vtx3 + n_vtx1 + 2*n_vtx2 - 1;
      _dface_vtx[1] = idxv + (ny-1)*n_vtx3 + n_vtx1 + 2*n_vtx2 - 1 + n_vtx_z_cst;
      _dface_vtx[2] = idxv + n_vtx_z_cst - 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv + n_vtx_z_cst - 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 4;
      _dface_cell[1] = 0;
    }

    else {
      // Quadrangle top-right corner _
      PDM_g_num_t idxv = (k+1)*n_vtx_z_cst + 1;
      _dface_vtx[0] = idxv - 5;
      _dface_vtx[1] = idxv - 1;
      _dface_vtx[2] = idxv - 1 + n_vtx_z_cst;
      _dface_vtx[3] = idxv - 5 + n_vtx_z_cst;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;

      _dface_cell[0] = k*n_faceH_z_cst + idx_tria5 + 4;
      _dface_cell[1] = 0;
    }

    // Flip horizontal face if z == zmax
    if (r < idx_quadV1 && k == nz) {
      int n_vtx = _dface_vtx_idx[ifac+1] - _dface_vtx_idx[ifac];
      for (int l = 0; l < n_vtx/2; l++) {
        PDM_g_num_t tmp = _dface_vtx[l];
        _dface_vtx[l] = _dface_vtx[n_vtx - l - 1];
        _dface_vtx[n_vtx - l - 1] = tmp;
      }
    }
  }
  *dface_vtx = realloc (*dface_vtx, sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);



  /*
   *  Face groups
   */
  *dface_group_idx = PDM_array_zeros_int (*n_face_group + 1);
  int *_dface_group_idx = *dface_group_idx;

  *dface_group = malloc (sizeof(PDM_g_num_t) * dn_face_lim);
  PDM_g_num_t *_dface_group = *dface_group;

  PDM_g_num_t idx_group1 = n_faceH_z_cst;
  PDM_g_num_t idx_group2 = idx_group1 + n_faceH_z_cst;
  PDM_g_num_t idx_group3 = idx_group2 + nz*(2*nx + 1);
  PDM_g_num_t idx_group4 = idx_group3 + nz*(2*nx + 1);
  PDM_g_num_t idx_group5 = idx_group4 + nz*(2*ny + 1);

  for (int ifac = 0; ifac < dn_face_lim; ifac++) {

    PDM_g_num_t g = distrib_face_lim[i_rank] + ifac;

    if (g < idx_group1) {
      // z = zmin
      _dface_group[_dface_group_idx[1]++] = g + 1;
    }

    else if (g < idx_group2) {
      // z = zmax
      if (_dface_group_idx[2] == 0) {
        _dface_group_idx[2] = _dface_group_idx[1];
      }

      PDM_g_num_t h = g - idx_group1;
      _dface_group[_dface_group_idx[2]++] = nz*n_face_z_cst + h + 1;
    }

    else if (g < idx_group3) {
      // y = ymin
      if (_dface_group_idx[3] == 0) {
        _dface_group_idx[3] = _dface_group_idx[2];
      }

      PDM_g_num_t h = g - idx_group2;
      PDM_g_num_t k = h / (2*nx + 1);
      PDM_g_num_t i = h % (2*nx + 1);
      PDM_g_num_t idxf = k*n_face_z_cst + 1;

      if (i < nx) {
        // _
        _dface_group[_dface_group_idx[3]++] = idxf + idx_quadV1 + i;
      } else if (i < 2*nx - 1) {
        // bottom row
        i -= nx;
        _dface_group[_dface_group_idx[3]++] = idxf + idx_quadV7 + i;
      } else if (i < 2*nx) {
        // bottom-left corner _
        _dface_group[_dface_group_idx[3]++] = idxf + idx_quadV11;
      } else {
        // bottom-right corner _
        _dface_group[_dface_group_idx[3]++] = idxf + idx_quadV11 + 2;
      }
    }

    else if (g < idx_group4) {
      // y = ymax
      if (_dface_group_idx[4] == 0) {
        _dface_group_idx[4] = _dface_group_idx[3];
      }

      PDM_g_num_t h = g - idx_group3;
      PDM_g_num_t k = h / (2*nx + 1);
      PDM_g_num_t i = h % (2*nx + 1);
      PDM_g_num_t idxf = k*n_face_z_cst + 1;

      if (i < nx) {
        // _
        _dface_group[_dface_group_idx[4]++] = idxf + idx_quadV1 + ny*nx + i;
      } else if (i < 2*nx - 1) {
        // top row
        i -= nx;
        _dface_group[_dface_group_idx[4]++] = idxf + idx_quadV8 + i;
      } else if (i < 2*nx) {
        // top-left corner _
        _dface_group[_dface_group_idx[4]++] = idxf + idx_quadV11 + 5;
      } else {
        // top-right corner_
        _dface_group[_dface_group_idx[4]++] = idxf + idx_quadV11 + 7;
      }
    }

    else if (g < idx_group5) {
      //  x = xmin
      if (_dface_group_idx[5] == 0) {
        _dface_group_idx[5] = _dface_group_idx[4];
      }

      PDM_g_num_t h = g - idx_group4;
      PDM_g_num_t k = h / (2*ny + 1);
      PDM_g_num_t j = h % (2*ny + 1);
      PDM_g_num_t idxf = k*n_face_z_cst + 1;

      if (j < ny) {
        // |
        _dface_group[_dface_group_idx[5]++] = idxf + idx_quadV2 + (nx + 1)*j;
      } else if (j < 2*ny - 1) {
        // left column
        j -= ny;
        _dface_group[_dface_group_idx[5]++] = idxf + idx_quadV9 + j;
      } else if (j < 2*ny) {
        // bottom-left corner |
        _dface_group[_dface_group_idx[5]++] = idxf + idx_quadV11 + 1;
      } else {
        // top-left corner |
        _dface_group[_dface_group_idx[5]++] = idxf + idx_quadV11 + 4;
      }
    }

    else {
      //  x = xmax
      if (_dface_group_idx[6] == 0) {
        _dface_group_idx[6] = _dface_group_idx[5];
      }

      PDM_g_num_t h = g - idx_group5;
      PDM_g_num_t k = h / (2*ny + 1);
      PDM_g_num_t j = h % (2*ny + 1);
      PDM_g_num_t idxf = k*n_face_z_cst + 1;

      if (j < ny) {
        // |
        _dface_group[_dface_group_idx[6]++] = idxf + idx_quadV2 + (nx + 1)*j + nx;
      } else if (j < 2*ny - 1) {
        // right column
        j -= ny;
        _dface_group[_dface_group_idx[6]++] = idxf + idx_quadV10 + j;
      } else if (j < 2*ny) {
        // bottom-right corner |
        _dface_group[_dface_group_idx[6]++] = idxf + idx_quadV11 + 3;
      } else {
        // top-right corner |
        _dface_group[_dface_group_idx[6]++] = idxf + idx_quadV11 + 6;
      }
    }
  }

  for (int i = 0; i < *n_face_group; i++) {
    _dface_group_idx[i+1] = PDM_MAX (_dface_group_idx[i+1], _dface_group_idx[i]);
  }


  /*
   *  Cell-face
   */
  *dcell_face_idx = malloc (sizeof(int) * (*dn_cell + 1));
  int *_dcell_face_idx = *dcell_face_idx;
  _dcell_face_idx[0] = 0;

  int s_cell_face = 10 * (*dn_cell);
  *dcell_face = malloc (sizeof(PDM_g_num_t) * s_cell_face);

  for (int icel = 0; icel < (*dn_cell); icel++) {

    PDM_g_num_t *_dcell_face = *dcell_face + _dcell_face_idx[icel];

    PDM_g_num_t g = distrib_cell[i_rank] + icel;
    PDM_g_num_t k = g / n_faceH_z_cst;
    PDM_g_num_t r = g % n_faceH_z_cst;
    PDM_g_num_t i, j;
    PDM_g_num_t idxf = k*n_face_z_cst + 1;

    if (r < idx_quadH) {
      // Octagon
      j = r / nx;
      i = r % nx;

      _dcell_face[0] = idxf + j*nx + i;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV3 + j*nx + i;         // \ +y
      _dcell_face[3] = idxf + idx_quadV1 + j*nx + i;         // _ +y
      _dcell_face[4] = idxf + idx_quadV4 + j*nx + i;         // / +y
      _dcell_face[5] = idxf + idx_quadV2 + j*(nx+1) + i + 1; // | -x
      _dcell_face[6] = idxf + idx_quadV6 + j*nx + i;         // \ -y
      _dcell_face[7] = idxf + idx_quadV1 + (j+1)*nx + i;     // _ -y
      _dcell_face[8] = idxf + idx_quadV5 + j*nx + i;         // / -y
      _dcell_face[9] = idxf + idx_quadV2 + j*(nx+1) + i;     // | +x
      if (i < nx-1) {
        _dcell_face[5] = -_dcell_face[5];
      }
      if (j < ny-1) {
        _dcell_face[7] = -_dcell_face[7];
      }
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }
      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 10;
    }

    else if (r < idx_tria1) {
      // Quadrangle
      PDM_g_num_t s = r - idx_quadH;

      j = s / (nx - 1);
      i = s % (nx - 1);

      _dcell_face[0] = idxf + idx_quadH + j*(nx-1) + i;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = -(idxf + idx_quadV6 + j*nx + i);         // \ -y
      _dcell_face[3] = -(idxf + idx_quadV5 + j*nx + i + 1);     // / -y
      _dcell_face[4] = -(idxf + idx_quadV3 + (j+1)*nx + i + 1); // \ +y
      _dcell_face[5] = -(idxf + idx_quadV4 + (j+1)*nx + i);     // / +y
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }
      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 6;
    }

    else if (r < idx_tria2) {
      // Triangle (bottom row)
      PDM_g_num_t s = r - idx_tria1;

      i = s % (nx - 1);

      _dcell_face[0] = idxf + idx_tria1 + i;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV7 + i;        // bottom row
      _dcell_face[3] = -(idxf + idx_quadV3 + i + 1); // \ +y
      _dcell_face[4] = -(idxf + idx_quadV4 + i);     // / +y
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria3) {
      // Triangle (top row)
      r -= idx_tria2;

      i = r % (nx - 1);

      _dcell_face[0] = idxf + idx_tria2 + i;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV8 + i;                    // top row
      _dcell_face[3] = -(idxf + idx_quadV6 + (ny-1)*nx + i);     // \ -y
      _dcell_face[4] = -(idxf + idx_quadV5 + (ny-1)*nx + i + 1); // / -y
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria4) {
      // Triangle (left column)
      PDM_g_num_t s = r - idx_tria3;

      j = s % (ny - 1);

      _dcell_face[0] = idxf + idx_tria3 + j;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV9 + j;           // left column
      _dcell_face[3] = -(idxf + idx_quadV5 + j*nx);     // / -y
      _dcell_face[4] = -(idxf + idx_quadV3 + (j+1)*nx); // \ +y
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria5) {
      // Triangle (right column)
      r -= idx_tria4;

      j = r % (ny - 1);

      _dcell_face[0] = idxf + idx_tria4 + j;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV10 + j;                   // right column
      _dcell_face[3] = -(idxf + idx_quadV4+ (j+1)*nx + nx - 1 ); // / +y
      _dcell_face[4] = -(idxf + idx_quadV6 + j*nx + nx - 1);     // \ -y
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria5 + 1) {
      // Triangle (bottom-left corner)
      _dcell_face[0] = idxf + idx_tria5;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV11;     // bottom-left corner _
      _dcell_face[3] = -(idxf + idx_quadV3);   // \ +y
      _dcell_face[4] = idxf + idx_quadV11 + 1; // bottom-left corner |
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria5 + 2) {
      // Triangle (bottom-right corner)
      _dcell_face[0] = idxf + idx_tria5 + 1;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV11 + 2;        // bottom-right corner _
      _dcell_face[3] = -(idxf + idx_quadV4 + nx - 1); // / +y
      _dcell_face[4] = idxf + idx_quadV11 + 3;        // bottom-right corner |
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else if (r < idx_tria5 + 3) {
      // Triangle (top-left corner)
      _dcell_face[0] = idxf + idx_tria5 + 2;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV11 + 4;           // top-left corner |
      _dcell_face[3] = -(idxf + idx_quadV5 + (ny-1)*nx); // / -y
      _dcell_face[4] = idxf + idx_quadV11 + 5;           // top-left corner _
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

    else {
      // Triangle (top-right corner)
      _dcell_face[0] = idxf + idx_tria5 + 3;
      _dcell_face[1] = _dcell_face[0] + n_face_z_cst;
      _dcell_face[2] = idxf + idx_quadV11 + 6;   // top-right corner |
      _dcell_face[3] = -(idxf + idx_quadV7 - 1); // \ -y
      _dcell_face[4] = idxf + idx_quadV11 + 7;   // top-right corner _
      if (k < nz-1) {
        _dcell_face[1] = -_dcell_face[1];
      }

      _dcell_face_idx[icel+1] = _dcell_face_idx[icel] + 5;
    }

  }
  *dcell_face = realloc (*dcell_face, sizeof(PDM_g_num_t) * _dcell_face_idx[*dn_cell]);


  /* TOUT DUR: switch orientation of all faces -> normals pointing outwards */
  if (1) {
    PDM_g_num_t tmp[8];
    for (int i = 0; i < *dn_face; i++) {

      PDM_g_num_t *dfv = *dface_vtx + (*dface_vtx_idx)[i];
      int n = (*dface_vtx_idx)[i+1] - (*dface_vtx_idx)[i];
      memcpy(tmp,
             dfv,
             sizeof(PDM_g_num_t) * n);

      for (int j = n-1; j >= 0; j--) {
        dfv[j] = tmp[n-j-1];
      }

    }
  }


  free (distrib_vtx);
  free (distrib_face);
  free (distrib_cell);
  free (distrib_face_lim);
}
