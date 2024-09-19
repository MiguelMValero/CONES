
/*============================================================================
 * Locate points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2005-2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "fvmc_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(FVMC_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_error.h>
#include <bftc_mem.h>
#include <bftc_printf.h>
#include <bftc_timer.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_ho_basis.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"
#include "fvmc_point_location.h"
#include "cwipi_config.h"

#ifdef CWP_HAVE_MKL
#include "mkl.h"
#else
#ifdef CWP_HAVE_BLAS
#include "cblas.h"
#endif
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_locator.h"
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+17
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a non blocking communication
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)
typedef struct  {
  void       **buffer;       /* tmp buffers for nonblocking */
  MPI_Request *MPI_request;  /* MPI requests for nonblocking */
  void        *var;          /* buffers updated in wait */
  const int   *local_list;   /* optional indirection list for
                                tmp_buffer -> buffer */
  int          stride;       /* dimension (1 for scalar,
                                3 for interlaced vector) */
  int          size;         /* size of type*/
  _Bool        reverse;      /* reverse mod for exchange */
} _fvmc_locator_nblocking_t;
#endif

/*----------------------------------------------------------------------------
 * Structure defining a locator
 *----------------------------------------------------------------------------*/

struct _fvmc_locator_t {

  /* Basic information */
  /*-------------------*/
  int  opt_bbox_step;        /* Associated step */
  double    tolerance;         /* Associated tolerance */

  _Bool     locate_on_parents; /* Locate relative to parent element numbers
                                  if true, element_id + 1 in concatenated
                                  sections of same element dimension if false */

  int       dim;               /* Spatial dimension */

#if defined(FVMC_HAVE_MPI)
  MPI_Comm  comm;              /* Associated MPI communicator */
#endif

  int       n_ranks;           /* Number of MPI ranks of distant location */
  int       start_rank;        /* First MPI rank of distant location */

  int       n_intersects;      /* Number of intersecting distant ranks */
  int      *intersect_rank;    /* List of intersecting distant ranks */
#if defined(FVMC_HAVE_MPI)
  _fvmc_locator_nblocking_t *nblockings_send; /*    */
  _fvmc_locator_nblocking_t *nblockings_recv; /*    */
  int                    max_nblockings_send;
  int                    max_nblockings_recv;
#endif
  double   *intersect_extents; /* List of intersecting distant extents */

  fvmc_lnum_t   *local_points_idx;   /* Start index of local points per rank
                                       (size: n_intersects + 1)*/
  fvmc_lnum_t   *local_distribution; /* Start index of distant points per rank
                                       (size: n_rank + 1)*/

  fvmc_lnum_t   *distant_points_idx; /* Start index of distant points per
                                       intersect rank
                                      (size: n_intersects + 1)*/
  fvmc_lnum_t   *distant_distribution; /* Start index of distant points per rank
                                       (size: n_rank + 1)*/

  fvmc_lnum_t   *local_point_ids;        /* Local point index for data received
                                           (with blocs starting at
                                           local_points_idx[] indexes,
                                           0 to n-1 numbering) */

  float         *distant_point_distance; /* Distance of distant points to location
                                            element number */
  fvmc_lnum_t   *distant_point_location; /* Location of distant points by parent
                                           element number (with blocs starting
                                           at distant_points_idx[] indexes) */
  fvmc_coord_t  *distant_point_coords;   /* Coordinates of distant points
                                           (with blocs starting at
                                           distant_points_idx[]*dim indexes) */
  fvmc_coord_t  *distant_point_projected_coords;   /* Coordinates of distant points projected on location
                                                    element (with blocs starting at
                                                    distant_points_idx[]*dim indexes) */

  double  *distant_point_uvw;   /* uvw of the distant points projected on location
                                       element (with blocs starting at
                                       distant_points_idx[]*max_n_node_elt indexes) */

  fvmc_lnum_t    n_interior;         /* Number of local points located */
  fvmc_lnum_t   *interior_list;      /* List (1 to n numbering) of points
                                       located */
  fvmc_lnum_t    n_exterior;         /* Number of local points not located */
  fvmc_lnum_t   *exterior_list;      /* List (1 to n numbering) of points
                                       not located */

  /* Timing information (2 fields/time; 0: total; 1: communication) */

  double  location_wtime[2];       /* Location Wall-clock time */
  double  location_cpu_time[2];    /* Location CPU time */
  double  exchange_wtime[2];       /* Variable exchange Wall-clock time */
  double  exchange_cpu_time[2];    /* Variable exchange CPU time */
  double  issend_wtime[2];       /* Variable exchange Wall-clock time */
  double  issend_cpu_time[2];    /* Variable exchange CPU time */
  double  irecv_wtime[2];       /* Variable exchange Wall-clock time */
  double  irecv_cpu_time[2];    /* Variable exchange CPU time */

};

/*============================================================================
 * Static global variables
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

static fvmc_locator_log_t   *_fvmc_locator_log_func = NULL;

/* global variables associated with communication logging */

static int  _fvmc_locator_log_start_p_comm = 0;
static int  _fvmc_locator_log_end_p_comm = 0;
static int  _fvmc_locator_log_start_g_comm = 0;
static int  _fvmc_locator_log_end_g_comm = 0;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log communication start.
 *
 * parameters:
 *   start_p_comm  <-- event number for the start of the "send/recv" or
 *                     "global send/recv" state
 *   timing        <-> 0: wall-clock total; 1 CPU total;
 *                     2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_start_comm(int      start_p_comm,
                          double   timing[4])
{
  timing[2] = bftc_timer_wtime();
  timing[3] = bftc_timer_cpu_time();

  if(_fvmc_locator_log_func != NULL)
    _fvmc_locator_log_func(start_p_comm, 0, NULL);
}

/*----------------------------------------------------------------------------
 * Log communication end.
 *
 * parameters:
 *   end_p_comm  <-- event number for the end of the "send/recv" or
 *                   "global send/recv" state
 *   timing      <-> 0: wall-clock total; 1 CPU total;
 *                   2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_end_comm(int      end_p_comm,
                        double   timing[4])
{
  if(_fvmc_locator_log_func != NULL)
    _fvmc_locator_log_func(end_p_comm, 0, NULL);

  timing[0] += bftc_timer_wtime() - timing[2];
  timing[1] += bftc_timer_cpu_time() - timing[3];
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_within_extents(int                dim,
                const fvmc_coord_t  coords[],
                const double       extents[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Update element extents with a given vertex
 *
 * parameters:
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_id         <-- vertex index (0 to n-1)
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   elt_extents       <-> extents associated with element:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   elt_initialized   <-> are extents already initialized for this vertex
 *                         (for all element vertices except the first) ?
 *----------------------------------------------------------------------------*/

inline static void
_update_elt_extents(int                 dim,
                    fvmc_lnum_t          vertex_id,
                    const fvmc_lnum_t   *parent_vertex_num,
                    const fvmc_coord_t   vertex_coords[],
                    double              elt_extents[],
                    _Bool              *elt_initialized)
{
  fvmc_lnum_t  i, coord_idx;

  if (parent_vertex_num == NULL)
    coord_idx = vertex_id;
  else
    coord_idx = parent_vertex_num[vertex_id] - 1;

  if (*elt_initialized == false) {
    for (i = 0; i < dim; i++) {
      elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }
    *elt_initialized = true;
  }
  else {
    for (i = 0; i < dim; i++) {
      if (elt_extents[i]       > vertex_coords[(coord_idx * dim) + i])
        elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      if (elt_extents[i + dim] < vertex_coords[(coord_idx * dim) + i])
        elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }

  }

}

/*----------------------------------------------------------------------------
 * Adjust element extents with search tolerance and update global extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   elt_dim     <-- element dimension
 *   tolerance   <-- addition to local extents of each element:
 *                   extent = base_extent * (1 + tolerance)
 *   elt_extents <-> extents associated with element:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

inline static void
_elt_extents_finalize(int               dim,
                      int               elt_dim,
                      double            tolerance,
                      double  *restrict elt_extents)
{
  int i;
  double delta[3];

  for (i = 0; i < dim; i++)
    delta[i] = (elt_extents[i+dim] - elt_extents[i]) * tolerance;

  if (elt_dim < dim) {
    double delta_max = delta[0];  /* for 1d or 2d elements, ensure */
    for (i = 0; i < dim; i++) {   /* search extent "thickness" */
      if (delta[i] > delta_max)
        delta_max = delta[i];
    }
    for (i = 0; i < dim; i++)
      delta[i] = delta_max;
  }

  for (i = 0; i < dim; i++) {
    elt_extents[i]     = elt_extents[i]     - delta[i];
    elt_extents[i+dim] = elt_extents[i+dim] + delta[i];
  }

}

/*----------------------------------------------------------------------------
 * Adjust extents with sub-extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   sub_extents <-> extents associated with element or section:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   extents     <-> optional section or mesh extents, to be updated:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim);
 *                   NULL if unused
 *----------------------------------------------------------------------------*/

inline static void
_update_extents(int               dim,
                double  *restrict sub_extents,
                double  *restrict extents)
{
  int i;

  for (i = 0; i < dim; i++) {
    if (sub_extents[i] < extents[i])
      extents[i] = sub_extents[i];
    if (sub_extents[i+dim] > extents[i+dim])
      extents[i+dim] = sub_extents[i+dim];
  }
}

/*----------------------------------------------------------------------------
 * Compute extents of a nodal mesh representation section
 *
 * parameters:
 *   this_section      <-- pointer to section structure
 *   dim               <-- spatial (coordinates) dimension
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   opt_bbox_step     <-- step to compute high order element extents
 *   tolerance         <-- addition to local extents of each element:
 *                         extent = base_extent * (1 + tolerance)
 *   extents           <-> extents associated with section:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/
static void
_nodal_section_extents(const fvmc_nodal_section_t  *this_section,
                       int                         order,
                       int                         dim,
                       const fvmc_lnum_t           *parent_vertex_num,
                       const fvmc_coord_t           vertex_coords[],
                       int                         opt_bbox_step,
                       double                      tolerance,
                       double                      extents[])
{
  fvmc_lnum_t  i, j, k, face_id, vertex_id;
  double elt_extents[6];

  /* initialize extents in case section is empty */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Extents for polyhedra */

  if (this_section->face_index != NULL) {

    for (i = 0; i < this_section->n_elements; i++) {

      _Bool elt_initialized = false;

      for (j = this_section->face_index[i];
           j < this_section->face_index[i + 1];
           j++) {
        face_id = FVMC_ABS(this_section->face_num[j]) - 1;
        for (k = this_section->vertex_index[face_id];
             k < this_section->vertex_index[face_id + 1];
             k++) {
          vertex_id = this_section->vertex_num[k] - 1;

          _update_elt_extents(dim,
                              vertex_id,
                              parent_vertex_num,
                              vertex_coords,
                              elt_extents,
                              &elt_initialized);

        }
      }

      _elt_extents_finalize(dim, 3, tolerance, elt_extents);
      _update_extents(dim, elt_extents, extents);

    }

  }

  /* Extents for polygons */

  else if (this_section->vertex_index != NULL) {

    fvmc_lnum_t  n_faces = (this_section->n_faces > 0) ?
                          this_section->n_faces : this_section->n_elements;

    for (i = 0; i < n_faces; i++) {

      _Bool elt_initialized = false;

      for (j = this_section->vertex_index[i];
           j < this_section->vertex_index[i + 1];
           j++) {
        vertex_id = this_section->vertex_num[j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(dim, 2, tolerance, elt_extents);
      _update_extents(dim, elt_extents, extents);

    }

  }

  /* Extents for regular elements */

  else {

    for (i = 0; i < this_section->n_elements; i++) {

      _Bool elt_initialized = false;

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      // TODO : if order > 1 add points to compute bounding box

      if (order > 1 && FVMC_ABS (opt_bbox_step) != 1) {

        if ((this_section->type == FVMC_EDGE)){

          const int n_step = opt_bbox_step;
          const double step = 1./(n_step - 1);
          const int n_nodes = this_section->stride;
          int n_vtx = 0;

          n_vtx = n_step + 1;

          double *u      = malloc(sizeof(double) * 1 * n_vtx);
          double *ai     = malloc(sizeof(double) * n_nodes * n_vtx);
          double *xyz    = malloc(sizeof(double) * 3 * n_vtx);
          double *coords = malloc(sizeof(double) * 3 * n_nodes);

          for (int jj = 0; jj < n_step ; jj++) {
            u[jj] = jj * step;
          }
          FVMC_ho_basis (FVMC_EDGE, order, n_nodes, n_vtx, u, ai);

          for (int ielt = 0; ielt < this_section->n_elements; ielt++) {

            for (int jj = 0; jj < n_nodes; jj++) {
              vertex_id = this_section->vertex_num[ielt*n_nodes + jj] - 1;
              int coord_idx;
              if (parent_vertex_num == NULL) {
                coord_idx = vertex_id;
              }
              else {
                coord_idx = parent_vertex_num[vertex_id] - 1;
              }
              for (int kk = 0; kk < 3; kk++) {
                coords[3*jj+kk] = vertex_coords[(coord_idx * 3) + kk];
              }
            }

#ifdef CWP_HAVE_BLAS
            double alpha = 1.;
            double beta = 0.;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        n_vtx, 3, n_nodes,
                        alpha,
                        ai, n_nodes,
                        coords, 3,
                        beta,
                        xyz, 3);
#else
            for (int ii = 0; ii < n_vtx; ii++) {

              for (int kk = 0; kk < 3; kk++) {
                xyz[3*ii + kk] = 0.;
              }

              for (int jj = 0; jj < n_nodes; jj++) {

                for (int kk = 0; kk < 3; kk++) {
                  xyz[3*ii +kk] += ai[ii*n_nodes+jj] * coords[3 * jj + kk];
                }
              }
            }

#endif

            for (int ii = 0; ii < n_vtx; ii++) {

              _update_elt_extents(dim,
                                  0,
                                  NULL,
                                  xyz + 3 *ii,
                                  elt_extents,
                                  &elt_initialized);
            }

          }

          free(u);
          free(ai);
          free(xyz);
          free(coords);

        }
        else if ((this_section->type == FVMC_FACE_TRIA) ||
                 (this_section->type == FVMC_FACE_QUAD)){

          assert (dim == 3);

          const int n_step = opt_bbox_step;
          const double step = 1./(n_step - 1);
          const int n_nodes = this_section->stride;

          int n_vtx = 0;

          if (this_section->type == FVMC_FACE_TRIA) {
            n_vtx = (n_step + 2) * (n_step + 1) /2;
          }
          else if (this_section->type == FVMC_FACE_QUAD) {
            n_vtx = (n_step + 1) * (n_step + 1);
          }

          double *uv = malloc (sizeof(double) * 2 * n_vtx);
          double *ai = malloc (sizeof(double) * n_nodes * n_vtx);
          double *xyz = malloc (sizeof(double) * 3 * n_vtx);
          double *coords =  malloc (sizeof(double) * 3 * n_nodes);

          if (this_section->type == FVMC_FACE_TRIA) {
            int i1 = 0;
            for (int jj = 0; jj < n_step + 1; jj++) {
              double v = jj*step;
              for (int ii = 0; ii < n_step + 1 - jj; ii++) {
                double u = ii*step;
                uv[i1++] = u;
                uv[i1++] = v;
              }
            }
            FVMC_ho_basis (FVMC_FACE_TRIA, order, n_nodes, n_vtx, uv, ai);
          }
          else if (this_section->type == FVMC_FACE_QUAD) {
            int i1 = 0;
            for (int jj = 0; jj < n_step + 1; jj++) {
              double v = jj*step;
              for (int ii = 0; ii < n_step + 1; ii++) {
                double u = ii*step;
                uv[i1++] = u;
                uv[i1++] = v;
              }
            }
            FVMC_ho_basis (FVMC_FACE_QUAD, order, n_nodes, n_vtx, uv, ai);
          }

          for (int ielt = 0; ielt < this_section->n_elements; ielt++) {

            for (int jj = 0; jj < n_nodes; jj++) {
              vertex_id = this_section->vertex_num[ielt*n_nodes + jj] - 1;
              int coord_idx;
              if (parent_vertex_num == NULL) {
                coord_idx = vertex_id;
              }
              else {
                coord_idx = parent_vertex_num[vertex_id] - 1;
              }
              for (int kk = 0; kk < 3; kk++) {
                coords[3*jj+kk] = vertex_coords[(coord_idx * 3) + kk];
              }
            }

#ifdef CWP_HAVE_BLAS
            double alpha = 1.;
            double beta = 0.;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        n_vtx, 3, n_nodes,
                        alpha,
                        ai, n_nodes,
                        coords, 3,
                        beta,
                        xyz, 3);

#else
            for (int ii = 0; ii < n_vtx; ii++) {

              for (int kk = 0; kk < 3; kk++) {
                xyz[3*ii + kk] = 0.;
              }

              for (int jj = 0; jj < n_nodes; jj++) {

                for (int kk = 0; kk < 3; kk++) {
                  xyz[3*ii +kk] += ai[ii*n_nodes+jj] * coords[3 * jj + kk];
                }
              }
            }

#endif

            for (int ii = 0; ii < n_vtx; ii++) {

              _update_elt_extents(dim,
                                  0,
                                  NULL,
                                  xyz + 3 *ii,
                                  elt_extents,
                                  &elt_initialized);
            }

          }


          free (uv);
          free (ai);
          free (xyz);
          free (coords);

        }
        else {
          assert (dim == 3);

          const int n_step = opt_bbox_step;
          const double step = 1./(n_step - 1);
          const int n_nodes = this_section->stride;

          int n_vtx = 0;

          if (this_section->type == FVMC_CELL_TETRA) {
            n_vtx =  (n_step+1)*(n_step+2)*(n_step+3)/6;
          }
          else if (this_section->type == FVMC_CELL_HEXA) {
            n_vtx = (n_step + 1) * (n_step + 1) * (n_step + 1);
          }
          else if (this_section->type == FVMC_CELL_PRISM) {
            n_vtx = (n_step + 1) * (n_step + 2) * (n_step + 1) /2;
          }
          else if (this_section->type == FVMC_CELL_PYRAM) {
            n_vtx = (n_step+1)*(n_step+2)*(2*n_step+3)/6;
          }

          double *uvw = malloc (sizeof(double) * 3 * n_vtx);
          double *ai  = malloc (sizeof(double) * n_nodes * n_vtx);
          double *xyz = malloc (sizeof(double) * 3 * n_vtx);
          double *coords =  malloc (sizeof(double) * 3 * n_nodes);

          if (this_section->type == FVMC_CELL_TETRA) {
            int i1 = 0;
            for (int kk = 0; kk < n_step + 1; kk++){
              double w = kk*step;
              for (int jj = 0; jj < n_step + 1 - kk; jj++) {
                double v = jj*step;
                for (int ii = 0; ii < n_step + 1 - jj - kk; ii++) {
                  double u = ii*step;
                  uvw[i1++] = u;
                  uvw[i1++] = v;
                  uvw[i1++] = w;
                }
              }
            }
            FVMC_ho_basis (FVMC_CELL_TETRA, order, n_nodes, n_vtx, uvw, ai);
          }


          if (this_section->type == FVMC_CELL_HEXA) {
            int i1 = 0;
            for (int kk = 0; kk < n_step + 1; kk++){
              double w = kk*step;
              for (int jj = 0; jj < n_step + 1; jj++) {
                double v = jj*step;
                for (int ii = 0; ii < n_step + 1; ii++) {
                  double u = ii*step;
                  uvw[i1++] = u;
                  uvw[i1++] = v;
                  uvw[i1++] = w;
                }
              }
            }
            FVMC_ho_basis (FVMC_CELL_HEXA, order, n_nodes, n_vtx, uvw, ai);
          }


          if (this_section->type == FVMC_CELL_PRISM) {
            int i1 = 0;
            for (int kk = 0; kk < n_step + 1; kk++){
              double w = kk*step;
              for (int jj = 0; jj < n_step + 1; jj++) {
                double v = jj*step;
                for (int ii = 0; ii < n_step + 1 - jj; ii++) {
                  double u = ii*step;
                  uvw[i1++] = u;
                  uvw[i1++] = v;
                  uvw[i1++] = w;
                }
              }
            }
            FVMC_ho_basis (FVMC_CELL_PRISM, order, n_nodes, n_vtx, uvw, ai);
          }


          if (this_section->type == FVMC_CELL_PYRAM) {
            int i1 = 0;
            for (int kk = 0; kk < n_step + 1; kk++){
              double w = kk*step;
              for (int jj = 0; jj < n_step + 1 - kk; jj++) {
                double v = jj*step;
                for (int ii = 0; ii < n_step + 1 - kk; ii++) {
                  double u = ii*step;
                  uvw[i1++] = u;
                  uvw[i1++] = v;
                  uvw[i1++] = w;
                }
              }
            }
            FVMC_ho_basis (FVMC_CELL_PYRAM, order, n_nodes, n_vtx, uvw, ai);
          }


          for (int ielt = 0; ielt < this_section->n_elements; ielt++) {

            for (int jj = 0; jj < n_nodes; jj++) {
              vertex_id = this_section->vertex_num[ielt*n_nodes + jj] - 1;
              int coord_idx;
              if (parent_vertex_num == NULL) {
                coord_idx = vertex_id;
              }
              else {
                coord_idx = parent_vertex_num[vertex_id] - 1;
              }
              for (int kk = 0; kk < 3; kk++) {
                coords[3*jj+kk] = vertex_coords[(coord_idx * 3) + kk];
              }
            }

#ifdef CWP_HAVE_BLAS
            double alpha = 1.;
            double beta = 0.;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        n_vtx, 3, n_nodes,
                        alpha,
                        ai, n_nodes,
                        coords, 3,
                        beta,
                        xyz, 3);
#else
            for (int ii = 0; ii < n_vtx; ii++) {

              for (int kk = 0; kk < 3; kk++) {
                xyz[3*ii + kk] = 0.;
              }

              for (int jj = 0; jj < n_nodes; jj++) {

                for (int kk = 0; kk < 3; kk++) {
                  xyz[3*ii +kk] += ai[ii*n_nodes+jj] * coords[3 * jj + kk];
                }
              }
            }

#endif

            for (int ii = 0; ii < n_vtx; ii++) {

              _update_elt_extents(dim,
                                  0,
                                  NULL,
                                  xyz + 3 *ii,
                                  elt_extents,
                                  &elt_initialized);
            }

          }


          free (uvw);
          free (ai);
          free (xyz);
          free (coords);


        }
      }

      _elt_extents_finalize(dim,
                            this_section->entity_dim,
                            tolerance,
                            elt_extents);
      _update_extents(dim, elt_extents, extents);





    }
  }

}

/*----------------------------------------------------------------------------
 * Compute extents of a nodal mesh representation
 *
 * parameters:
 *   this_nodal   <-- pointer to mesh representation structure
 *   tolerance    <-- addition to local extents of each element:
 *                    extent = base_extent * (1 + tolerance)
 *   extents      <-> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_nodal_extents(const fvmc_nodal_t  *this_nodal,
               int                  opt_bbox_step,
               double              tolerance,
               double              extents[])
{
  int i, j;
  int dim;
  double section_extents[6];

  if (this_nodal == NULL)
    return;

  dim = this_nodal->dim;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  for (i = 0; i < this_nodal->n_sections; i++) {

    _nodal_section_extents(this_nodal->sections[i],
                           this_nodal->order,
                           this_nodal->dim,
                           this_nodal->parent_vertex_num,
                           this_nodal->vertex_coords,
                           opt_bbox_step,
                           tolerance,
                           section_extents);

    for (j = 0; j < this_nodal->dim; j++) {
      if (section_extents[j] < extents[j])
        extents[j] = section_extents[j];
      if (section_extents[j+dim] > extents[j+dim])
        extents[j+dim] = section_extents[j+dim];
    }

  }

}

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim          <-- space dimension of points to locate
 *   n_points     <-- number of points to locate
 *   point_list   <-- optional indirection array to point_coords
 *                    (1 to n_points numbering)
 *   point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(int                  dim,
               fvmc_lnum_t           n_points,
               const fvmc_lnum_t     point_list[],
               const fvmc_coord_t    point_coords[],
               double                extents[])
{
  int i;
  fvmc_lnum_t j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_list != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_list[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given nodal mesh representation and
 * point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   this_nodal        <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_all_distant(fvmc_locator_t       *this_locator,
                    const fvmc_nodal_t   *this_nodal,
                    int                  dim,
                    fvmc_lnum_t           n_points,
                    const fvmc_lnum_t     point_list[],
                    const fvmc_coord_t    point_coords[])
{
  int i, k, stride;
  int dist_rank, dist_index;
  fvmc_lnum_t j;
  fvmc_lnum_t n_coords_loc, n_coords_dist, n_interior, n_exterior;
  fvmc_lnum_t coord_idx, start_idx;
  fvmc_lnum_t *location_loc, *location_dist;
  fvmc_lnum_t *location, *location_rank_id;
  fvmc_lnum_t *send_index, *send_location;
  fvmc_lnum_t *location_count, *location_shift;
  fvmc_coord_t *coords_dist, *send_coords;
  float *distance, *send_distance;
  float *distance_dist, *distance_loc;
  double extents[6];

  MPI_Status status;

  double comm_timing[4] = {0., 0., 0., 0.};

  int local_max_entity_dim = fvmc_nodal_get_max_entity_dim (this_nodal);
  int distant_max_entity_dim;
  for (i = 0; i < this_locator->n_intersects; i++) {
    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = this_locator->intersect_rank[dist_index];

    MPI_Sendrecv(&local_max_entity_dim, 1, FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 &distant_max_entity_dim, 1, FVMC_MPI_LNUM, dist_rank,
                 FVMC_MPI_TAG, this_locator->comm, &status);
  }


  /* Initialization */

  stride = dim * 2;

  BFTC_MALLOC(send_coords, n_points * dim, fvmc_coord_t);
  BFTC_MALLOC(send_index, n_points, fvmc_lnum_t);

  BFTC_MALLOC(location, n_points, fvmc_lnum_t);
  BFTC_MALLOC(location_rank_id, n_points, fvmc_lnum_t);
  BFTC_MALLOC(distance, n_points, float);

  int local_nodal_order = fvmc_nodal_order_get (this_nodal);
  int distant_nodal_order = -1;

  /* Pas terrible de faire une boucle  */

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = this_locator->intersect_rank[dist_index];

    MPI_Sendrecv(&local_nodal_order, 1, FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 &distant_nodal_order, 1, FVMC_MPI_LNUM, dist_rank,
                 FVMC_MPI_TAG, this_locator->comm, &status);
  }

  fvmc_coord_t *projected_coords = NULL;
  double *uvw = NULL;
  int max_entity_dim = FVMC_MAX(local_max_entity_dim, distant_max_entity_dim);

  if (distant_nodal_order != -1) {
    BFTC_MALLOC(projected_coords, n_points*dim, fvmc_coord_t);
    BFTC_MALLOC(uvw, n_points*max_entity_dim, double); // which one??
  }

  //  int max_entity_dim = fvmc_nodal_3 (this_nodal);

  for (j = 0; j < n_points; j++) {
    location[j] = -1;
    location_rank_id[j] = -1;
    distance[j] = -1.0;
  }

  /* First loop on possibly intersecting distant ranks */
  /*---------------------------------------------------*/

  //  int curr_rank;
  //MPI_Comm_rank (this_locator->comm, &curr_rank);

  /* for (i = 0; i < this_locator->n_intersects; i++) { */

  /*   dist_index = i; /\* Ordering (communication schema) not yet optimized *\/ */
  /*   dist_rank  = this_locator->intersect_rank[dist_index]; */
  /*   printf ("d %d ---> %d\n", curr_rank, dist_rank); */

  /*   fvmc_lnum_t ll1 = 0; */
  /*   fvmc_lnum_t ll2 = 0; */

  /*   MPI_Sendrecv(&ll1, 1, FVMC_MPI_LNUM, dist_rank, 10, */
  /*                &ll2, 1, FVMC_MPI_LNUM, dist_rank, */
  /*                10, this_locator->comm, &status); */

  /*   printf ("f %d ---> %d\n", curr_rank, dist_rank); */
  /* } */
  /* fflush(stdout); */

  /* MPI_Barrier(this_locator->comm); */
  
  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = this_locator->intersect_rank[dist_index];

    /* printf ("%3.3d d_boucle 1 %d ---> %d\n", curr_rank, curr_rank, dist_rank); */
    /* fflush(stdout); */

    /* Prepare and send coords that should fit in each send buffer */
    /* Reset buffers for current intersect rank */

    n_coords_loc = 0;

    for (k = 0; k < stride; k++)
      extents[k] = this_locator->intersect_extents[dist_index*stride + k];

    /* Build partial buffer */

    for (j = 0; j < n_points; j++) {

      if (point_list != NULL)
        coord_idx = point_list[j] - 1;
      else
        coord_idx = j;

      if (_within_extents(dim,
                          &(point_coords[dim*coord_idx]),
                          extents) == true) {

        send_index[n_coords_loc] = j;
        for (k = 0; k < dim; k++)
          send_coords[n_coords_loc*dim + k]
            = point_coords[dim*coord_idx + k];

        n_coords_loc += 1;
      }

    }

    /* Send then receive partial buffer */

    dist_rank = this_locator->intersect_rank[dist_index];

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 &n_coords_dist, 1, FVMC_MPI_LNUM, dist_rank,
                 FVMC_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    BFTC_MALLOC(coords_dist, n_coords_dist*dim, fvmc_coord_t);

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                 coords_dist, (int)(n_coords_dist*dim),
                 FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    /* Now locate received coords on local rank */

    BFTC_MALLOC(location_dist, n_coords_dist, fvmc_lnum_t);
    BFTC_MALLOC(distance_dist, n_coords_dist, float);


    for (j = 0; j < n_coords_dist; j++) {
      location_dist[j] = -1;
      distance_dist[j] = -1.0;
    }
    if (fvmc_nodal_order_get (this_nodal) != -1) {
      for (j = 0; j < n_coords_dist; j++) {
        distance_dist[j] = HUGE_VAL;
      }
    }
    else {
      for (j = 0; j < n_coords_dist; j++) {
        distance_dist[j] = -1;
      }
    }

    fvmc_coord_t *projected_coords_dist = NULL;
    double *uvw_dist = NULL;
    if (local_nodal_order != -1) {
      BFTC_MALLOC(projected_coords_dist, n_coords_dist*dim, fvmc_coord_t);
      BFTC_MALLOC(uvw_dist, n_coords_dist * local_max_entity_dim, double);
    }


    fvmc_point_location_nodal(this_nodal,
                              this_locator->tolerance,
                              this_locator->locate_on_parents,
                              n_coords_dist,
                              coords_dist,
                              projected_coords_dist,
                              uvw_dist,
                              location_dist,
                              distance_dist);

    BFTC_FREE(coords_dist);

    /* Exchange location return information with distant rank */

    BFTC_MALLOC(location_loc, n_coords_loc, fvmc_lnum_t);
    BFTC_MALLOC(distance_loc, n_coords_loc, float);

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(location_dist, (int)n_coords_dist,
                 FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 location_loc, (int)n_coords_loc,
                 FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(distance_dist, (int)n_coords_dist,
                 MPI_FLOAT, dist_rank, FVMC_MPI_TAG,
                 distance_loc, (int)n_coords_loc,
                 MPI_FLOAT, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    fvmc_coord_t *projected_coords_loc = NULL;
    double *uvw_loc = NULL;

    if ((local_nodal_order != -1) && (distant_nodal_order != -1)) {
      BFTC_MALLOC(projected_coords_loc, n_coords_loc * dim, fvmc_coord_t);
      MPI_Sendrecv(projected_coords_dist, (int)n_coords_dist * dim,
                   FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                   projected_coords_loc, (int)n_coords_loc * dim,
                   FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);

      BFTC_MALLOC(uvw_loc, n_coords_loc * distant_max_entity_dim, double);
      printf("sendrecv uvw\n");
      fflush(stdout);
      MPI_Sendrecv(uvw_dist, (int)n_coords_dist * local_max_entity_dim,
                   MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
                   uvw_loc, (int)n_coords_loc * distant_max_entity_dim,
                   MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);
      printf("sendrecv uvw out\n");
      fflush(stdout);
    }

    else if (local_nodal_order != -1) {
      MPI_Send(projected_coords_dist, (int)n_coords_dist * dim,
               FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
               this_locator->comm);

      MPI_Send(uvw_dist, (int)n_coords_dist * local_max_entity_dim,
               MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
               this_locator->comm);
    }

    else if (distant_nodal_order != -1) {
      BFTC_MALLOC(projected_coords_loc, n_coords_loc * dim, fvmc_coord_t);
      MPI_Recv(projected_coords_loc, (int)n_coords_loc * dim,
               FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
               this_locator->comm, &status);

      BFTC_MALLOC(uvw_loc, n_coords_loc * distant_max_entity_dim, double);
      MPI_Recv(uvw_loc, (int)n_coords_loc * distant_max_entity_dim,
               MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
               this_locator->comm, &status);
    }


    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    BFTC_FREE(location_dist);
    BFTC_FREE(distance_dist);
    if (projected_coords_dist != NULL) {
      BFTC_FREE(projected_coords_dist);
    }
    if (uvw_dist != NULL) {
      BFTC_FREE(uvw_dist);
    }

    /* Now update location information */

    for (j = 0; j < n_coords_loc; j++) {

      fvmc_lnum_t l = send_index[j];

      if ((location_loc[j] > -1) &&  ((distance_loc[j] > -0.1)
          && (distance_loc[j] < distance[l] || distance[l] < -0.1))) {
        location_rank_id[l] = i;
        location[l] = location_loc[j];
        distance[l] = distance_loc[j];
        if (projected_coords != NULL) {
          for (int l1 = 0; l1 < dim; l1++) {
            projected_coords[l*dim + l1] = projected_coords_loc[j*dim + l1];
          }
        }
        if (uvw != NULL) {
          for (int l1 = 0; l1 < distant_max_entity_dim; l1++) { // ??
            uvw[l*distant_max_entity_dim + l1] = uvw_loc[j*distant_max_entity_dim + l1];
          }
        }
      }

    }

    BFTC_FREE(location_loc);
    BFTC_FREE(distance_loc);
    if (projected_coords_loc != NULL) {
      BFTC_FREE(projected_coords_loc);
    }
    if (uvw_loc != NULL) {
      BFTC_FREE(uvw_loc);
    }

  }

  /* Reorganize localization information */
  /*-------------------------------------*/

  /* Now that localization is done, the location[] array contains
     either -1 if a point was not localized, or a local index
     (associated with the corresponding rank); the distance[] array
     is not needed anymore now that all comparisons have been done */


  BFTC_MALLOC(location_shift, this_locator->n_intersects, fvmc_lnum_t);
  BFTC_MALLOC(location_count, this_locator->n_intersects, fvmc_lnum_t);

  for (i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    if (location_rank_id[j] > -1)
      location_count[location_rank_id[j]] += 1;
    else
      n_exterior += 1;
  }

  this_locator->n_interior = n_points - n_exterior;
  BFTC_MALLOC(this_locator->interior_list, this_locator->n_interior, fvmc_lnum_t);

  this_locator->n_exterior = n_exterior;
  BFTC_MALLOC(this_locator->exterior_list, this_locator->n_exterior, fvmc_lnum_t);

  if (this_locator->n_intersects > 0)
    location_shift[0] = 0;
  for (i = 1; i < this_locator->n_intersects; i++)
    location_shift[i] = location_shift[i-1] + location_count[i-1];

  for (i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  /* send_index[] will now contain information for all blocks */
  for (j = 0; j < n_points; j++)
    send_index[j] = -1;

  BFTC_MALLOC(send_location, n_points, fvmc_lnum_t);
  BFTC_MALLOC(send_distance, n_points, float);

  fvmc_coord_t *send_projected_coords = NULL;
  if (distant_nodal_order != -1) {
    BFTC_MALLOC(send_projected_coords, n_points * dim, fvmc_coord_t);
  }

  fvmc_coord_t *send_uvw = NULL;
  if (distant_nodal_order != -1) {
    BFTC_MALLOC(send_uvw, n_points * distant_max_entity_dim , double);
  }

  n_interior = 0;
  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    const int l_rank = location_rank_id[j];
    if (l_rank > -1) {
      send_index[location_shift[l_rank] + location_count[l_rank]] = j;
      location_count[l_rank] += 1;
      this_locator->interior_list[n_interior] = j + 1;
      n_interior += 1;
    }
    else {
      this_locator->exterior_list[n_exterior] = j + 1;
      n_exterior += 1;
    }
  }


  /* Second loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  /* Count and organize total number of local and distant points */

  BFTC_MALLOC(this_locator->local_points_idx,
             this_locator->n_intersects + 1,
             fvmc_lnum_t);

  BFTC_MALLOC(this_locator->distant_points_idx,
             this_locator->n_intersects + 1,
             fvmc_lnum_t);

  BFTC_MALLOC(this_locator->distant_distribution,
             this_locator->n_ranks + 1,
             fvmc_lnum_t);

  BFTC_MALLOC(this_locator->local_distribution,
             this_locator->n_ranks + 1,
             fvmc_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->distant_points_idx[0] = 0;

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_coords_loc = location_count[i];

    this_locator->local_points_idx[i+1]
      = this_locator->local_points_idx[i] + n_coords_loc;

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 &n_coords_dist, 1, FVMC_MPI_LNUM, dist_rank,
                 FVMC_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    this_locator->distant_points_idx[i+1]
      = this_locator->distant_points_idx[i] + n_coords_dist;

  }

  if (this_locator->n_intersects > 0) {
    k = 0;
    for (i = 0; i < this_locator->n_ranks; i++) {
      if ((i < this_locator->intersect_rank[k]) &&
          (i > this_locator->intersect_rank[k])) {
        this_locator->distant_distribution[i+1] = this_locator->distant_distribution[i];
        this_locator->local_distribution[i+1]   = this_locator->local_distribution[i];
      }
      else if (i == this_locator->intersect_rank[k]) {
        this_locator->distant_distribution[i+1] = this_locator->distant_points_idx[k+1];
        this_locator->local_distribution[i+1]   = this_locator->local_points_idx[k+1];
        if (k < this_locator->n_intersects - 1)
          k += 1;
      }
    }
  }


  /* Third loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  BFTC_MALLOC(this_locator->local_point_ids,
             this_locator->local_points_idx[this_locator->n_intersects],
             fvmc_lnum_t);

  BFTC_MALLOC(this_locator->distant_point_location,
             this_locator->distant_points_idx[this_locator->n_intersects],
             fvmc_lnum_t);

  BFTC_MALLOC(this_locator->distant_point_distance,
             this_locator->distant_points_idx[this_locator->n_intersects],
             float);

  BFTC_MALLOC(this_locator->distant_point_coords,
             this_locator->distant_points_idx[this_locator->n_intersects] * dim,
             fvmc_coord_t);

  this_locator->distant_point_projected_coords = NULL;

  if (local_nodal_order != -1) {

    BFTC_MALLOC(this_locator->distant_point_projected_coords,
                this_locator->distant_points_idx[this_locator->n_intersects] * dim,
                fvmc_coord_t);
  }

  this_locator->distant_point_uvw = NULL;

  if (local_nodal_order != -1) {

    BFTC_MALLOC(this_locator->distant_point_uvw,
                this_locator->distant_points_idx[this_locator->n_intersects] * local_max_entity_dim,
                double);

  }

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_coords_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_coords_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    start_idx = this_locator->local_points_idx[i];

    for (j = 0; j < n_coords_loc; j++) {
      coord_idx = send_index[location_shift[i] + j];
      this_locator->local_point_ids[start_idx + j] = coord_idx;
      send_location[j] = location[coord_idx];
      send_distance[j] = distance[coord_idx];
      if (send_projected_coords != NULL) {
        for (k = 0; k < dim; k++) {
          send_projected_coords[j*dim + k] = projected_coords[dim*coord_idx + k];
        }
      }
      if (send_uvw != NULL) {
        for (k = 0; k < distant_max_entity_dim; k++) {
          send_uvw[j*distant_max_entity_dim + k] = uvw[distant_max_entity_dim*coord_idx + k]; // ??
        }
      }
      if (point_list != NULL) {
        for (k = 0; k < dim; k++)
          send_coords[j*dim + k]
            = point_coords[dim*(point_list[coord_idx] - 1) + k];
      }
      else {
        for (k = 0; k < dim; k++)
          send_coords[j*dim + k]
            = point_coords[dim*coord_idx + k];
      }
    }

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_location, (int)n_coords_loc,
                 FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 (this_locator->distant_point_location
                  + this_locator->distant_points_idx[i]), (int)n_coords_dist,
                 FVMC_MPI_LNUM, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(send_distance, (int)n_coords_loc,
                 MPI_FLOAT, dist_rank, FVMC_MPI_TAG,
                 (this_locator->distant_point_distance
                  + this_locator->distant_points_idx[i]), (int)n_coords_dist,
                 MPI_FLOAT, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                 (this_locator->distant_point_coords
                  + (this_locator->distant_points_idx[i]*dim)),
                 (int)(n_coords_dist*dim),
                 FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    if ((local_nodal_order != -1) && (distant_nodal_order != -1)) {
      MPI_Sendrecv(send_projected_coords, (int)(n_coords_loc*dim),
                   FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                   (this_locator->distant_point_projected_coords
                    + (this_locator->distant_points_idx[i]*dim)),
                   (int)(n_coords_dist*dim),
                   FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);

      MPI_Sendrecv(send_uvw, (int)(n_coords_loc*distant_max_entity_dim), // ??
                   MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
                   (this_locator->distant_point_uvw
                    + (this_locator->distant_points_idx[i]*local_max_entity_dim)), // ??
                   (int)(n_coords_dist*max_entity_dim),
                   MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);

    }

    else if (distant_nodal_order != -1) {
      MPI_Send(send_projected_coords, (int)(n_coords_loc*dim),
               FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
               this_locator->comm);

      MPI_Send(send_uvw, (int)(n_coords_loc*distant_max_entity_dim),
               MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
               this_locator->comm);
    }

    else if (local_nodal_order != -1) {
      MPI_Recv((this_locator->distant_point_projected_coords
                + (this_locator->distant_points_idx[i]*dim)),
               (int)(n_coords_dist*dim),
               FVMC_MPI_COORD, dist_rank, FVMC_MPI_TAG,
               this_locator->comm, &status);

      MPI_Recv((this_locator->distant_point_uvw
                + (this_locator->distant_points_idx[i]*local_max_entity_dim)), // ??
               (int)(n_coords_dist*local_max_entity_dim),
               MPI_DOUBLE, dist_rank, FVMC_MPI_TAG,
               this_locator->comm, &status);
    }

    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

  }

  BFTC_FREE(location_count);
  BFTC_FREE(location_shift);

  BFTC_FREE(send_index);
  BFTC_FREE(send_location);
  BFTC_FREE(send_distance);
  BFTC_FREE(send_coords);
  if (send_projected_coords != NULL) {
    BFTC_FREE(send_projected_coords);
  }
  if (send_uvw != NULL) {
    BFTC_FREE(send_uvw);
  }

  BFTC_FREE(location_rank_id);

  BFTC_FREE(location);
  BFTC_FREE(distance);
  if (projected_coords != NULL) {
    BFTC_FREE(projected_coords);
  }
  if (uvw != NULL) {
    BFTC_FREE(uvw);
  }

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];

  fflush(stdout);


}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given nodal mesh representation and
 * point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   this_nodal        <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_all_local(fvmc_locator_t       *this_locator,
                  const fvmc_nodal_t   *this_nodal,
                  int                  dim,
                  fvmc_lnum_t           n_points,
                  const fvmc_lnum_t     point_list[],
                  const fvmc_coord_t    point_coords[])
{
  int l, stride;
  fvmc_lnum_t j, k;
  fvmc_lnum_t n_coords, n_interior, n_exterior, coord_idx;
  fvmc_lnum_t *location;
  fvmc_lnum_t location_count;
  fvmc_coord_t *coords;
  float *distance;
  double extents[6];

  /* Initialization */

  stride = this_nodal->dim * 2;

  BFTC_MALLOC(coords, n_points * this_nodal->dim, fvmc_coord_t);

  /* Initialize localization information */
  /*-------------------------------------*/

  n_coords = 0;

  for (k = 0; k < stride; k++)
    extents[k] = this_locator->intersect_extents[k];

  /* Build partial buffer */

  for (j = 0; j < n_points; j++) {

    if (point_list != NULL)
      coord_idx = point_list[j] - 1;
    else
      coord_idx = j;

    if (_within_extents(dim,
                        &(point_coords[dim*coord_idx]),
                        extents) == true) {

      for (k = 0; k < dim; k++)
        coords[n_coords*dim + k]
          = point_coords[dim*coord_idx + k];

      n_coords += 1;
    }

  }

  BFTC_REALLOC(coords, n_coords * this_nodal->dim, fvmc_coord_t);

 /*  Now locate coords */

  BFTC_MALLOC(location, n_coords, fvmc_lnum_t);
  BFTC_MALLOC(distance, n_coords, float);

  for (j = 0; j < n_coords; j++) {
    location[j] = -1;
    distance[j] = -1.0;
  }

  fvmc_point_location_nodal(this_nodal,
                            this_locator->tolerance,
                            this_locator->locate_on_parents,
                            n_coords,
                            coords,
                            NULL,
                            NULL,
                            location,
                            distance);


  /* Reorganize localization information */
  /*-------------------------------------*/

  /* Now that localization is done, the location[] array contains
     either -1 if a point was not localized, or a local index;
     the distance[] array is not needed anymore now that all comparisons have
     been done */


  location_count = 0;

  n_exterior = 0;
  for (j = 0; j < n_coords; j++) {
    if (location[j] > -1)
      location_count += 1;
    else
      n_exterior += 1;
  }

  this_locator->n_interior = n_coords - n_exterior;
  BFTC_MALLOC(this_locator->interior_list, this_locator->n_interior, fvmc_lnum_t);

  this_locator->n_exterior = (n_points - n_coords) + n_exterior;
  BFTC_MALLOC(this_locator->exterior_list, this_locator->n_exterior, fvmc_lnum_t);

  /* Organize total number of "local" and "distant" points */

  BFTC_MALLOC(this_locator->local_points_idx, 2, fvmc_lnum_t);
  BFTC_MALLOC(this_locator->distant_points_idx, 2, fvmc_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->local_points_idx[1] = location_count;

  this_locator->distant_points_idx[0] = 0;
  this_locator->distant_points_idx[1] = location_count;

  this_locator->local_point_ids = NULL; /* Not needed for single-process */

  BFTC_MALLOC(this_locator->distant_point_location, location_count, fvmc_lnum_t);
  BFTC_MALLOC(this_locator->distant_point_distance, location_count, float);
  BFTC_MALLOC(this_locator->distant_point_coords, n_coords * dim, fvmc_coord_t);

  location_count = 0;
  n_interior = 0;
  n_exterior = 0;

  for (j = 0, k = 0; j < n_points; j++) {

    if (point_list != NULL)
      coord_idx = point_list[j] - 1;
    else
      coord_idx = j;

    if (_within_extents(dim,
                        &(point_coords[dim*coord_idx]),
                        extents) == true) {

      if (location[k] > -1) {
        this_locator->distant_point_location[location_count] = location[k];
        this_locator->distant_point_distance[location_count] = distance[k];
        for (l = 0; l < dim; l++) {
          this_locator->distant_point_coords[location_count*dim + l]
            = point_coords[coord_idx*dim + l];
        }
        location_count += 1;
        this_locator->interior_list[n_interior] = j + 1;
        n_interior += 1;

      }
      else {
        this_locator->exterior_list[n_exterior] = j + 1;
        n_exterior += 1;
      }

      k += 1;

    }
    else {
      this_locator->exterior_list[n_exterior] = j + 1;
      n_exterior += 1;
    }

  }

  BFTC_FREE(location);
  BFTC_FREE(distance);
  BFTC_FREE(coords);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get a free request
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   request       --> free request
 *----------------------------------------------------------------------------*/

static int
_get_free_request(_fvmc_locator_nblocking_t **nblockings,
                  int *max_nblockings)
{

  int request = 0;
  if (*nblockings == NULL) {

    *max_nblockings = 4;
    BFTC_MALLOC(*nblockings,
               *max_nblockings,
               _fvmc_locator_nblocking_t);
    for (int i = 0; i < *max_nblockings; i++) {
      (*nblockings)[i].buffer = NULL;
      (*nblockings)[i].MPI_request = NULL;
      (*nblockings)[i].var = NULL;
      (*nblockings)[i].local_list = NULL;
      (*nblockings)[i].reverse = false;
      (*nblockings)[i].stride = 0;
      (*nblockings)[i].size = 0;
    }
  }

  else {
    for (request = 0; request < *max_nblockings; request++) {
      if ((*nblockings)[request].MPI_request == NULL) {
        break;
      }
    }
    if (request == *max_nblockings) {
      *max_nblockings = 2 * (*max_nblockings);
      BFTC_REALLOC(*nblockings,
                  *max_nblockings,
                  _fvmc_locator_nblocking_t);
      for (int i = request; i < *max_nblockings; i++) {
        (*nblockings)[i].buffer = NULL;
        (*nblockings)[i].MPI_request = NULL;
        (*nblockings)[i].var = NULL;
        (*nblockings)[i].local_list = NULL;
        (*nblockings)[i].reverse = false;
        (*nblockings)[i].stride = 0;
        (*nblockings)[i].size = 0;
      }
    }
  }

  return request;
}


/*----------------------------------------------------------------------------
 *
 * Non blocking send
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           <-- variable defined on distant points (distant_var)
 *                     or on local points (local_var) if reverse
 *   local_var     <-- variable defined on local points
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, send is reversed
 *                     (send variable defined on local points)
 *   tag           <-- tag for MPI_issend
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

static void
_issend_point_var_distant(fvmc_locator_t     *this_locator,
                          void              *var,
                          const fvmc_lnum_t  *local_list,
                          MPI_Datatype       datatype,
                          size_t             stride,
                          _Bool              reverse,
                          int                tag,
                          int               *request)
{
  int i, dist_v_count, loc_v_count, size;
  int dist_rank, dist_index;
  fvmc_lnum_t n_points_loc,  n_points_dist;
  size_t dist_v_idx;
  void *dist_v_ptr;

  double comm_timing[4] = {0., 0., 0., 0.};

  MPI_Aint extent;

  /* Get a free  request */

  *request = _get_free_request(&(this_locator->nblockings_send),
                               &(this_locator->max_nblockings_send));

  _fvmc_locator_nblocking_t *nblocking_send = this_locator->nblockings_send + *request;

  BFTC_MALLOC(nblocking_send->MPI_request, this_locator->n_intersects, MPI_Request);
  if (reverse)
    BFTC_MALLOC(nblocking_send->buffer, this_locator->n_intersects, void *);
  nblocking_send->reverse = reverse;

  /* Check extent of datatype */
  MPI_Aint lb;
  MPI_Type_get_extent(datatype, &lb, &extent);
  MPI_Type_size(datatype, &size);

  if (extent != size)
    bftc_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Loop on possibly intersecting distant ranks */
  /*---------------------------------------------*/

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    /* Exchange information */

    if (reverse == false) {

      void *distant_var = var;

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
      dist_v_count = n_points_dist * stride;

      if (distant_var != NULL)
        dist_v_ptr = (void *)(((char *)distant_var) + dist_v_idx);
      else
        dist_v_ptr = NULL;

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Issend(dist_v_ptr, dist_v_count, datatype, dist_rank, tag,
                this_locator->comm, nblocking_send->MPI_request + i);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    }

    else { /* if (reverse == true) */

      /* Initialization */

      void *local_var = var;

      const fvmc_lnum_t *_local_point_ids
        = this_locator->local_point_ids + this_locator->local_points_idx[i];

      n_points_loc =    this_locator->local_points_idx[i+1]
                       - this_locator->local_points_idx[i];

      BFTC_MALLOC(nblocking_send->buffer[i], n_points_loc*size*stride, char);

      loc_v_count = n_points_loc*stride;

      if (local_list == NULL) {
        int k;
        size_t l;
        const size_t nbytes = stride*size;
        for (k = 0; k < n_points_loc; k++) {
          const char *local_v_p
            = (const char *)local_var + _local_point_ids[k]*nbytes;
          char *loc_v_buf_p = (char *)nblocking_send->buffer[i] + k*nbytes;
          for (l = 0; l < nbytes; l++)
            loc_v_buf_p[l] = local_v_p[l];
        }
      }
      else {
        int k;
        size_t l;
        const size_t nbytes = stride*size;
        for (k = 0; k < n_points_loc; k++) {
          const char *local_v_p
            = (const char *)local_var
            + (local_list[_local_point_ids[k]] - 1)*nbytes;
          char *loc_v_buf_p = (char *)nblocking_send->buffer[i] + k*nbytes;
          for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
        }
      }

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Issend(nblocking_send->buffer[i], loc_v_count, datatype, dist_rank, tag,
                this_locator->comm, nblocking_send->MPI_request + i);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    }

  } /* End of loop on possibly intersecting ranks */

  this_locator->issend_wtime[1] += comm_timing[0];
  this_locator->issend_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 *
 * Non blocking receive
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           --> variable defined on distant points (distant_var)
 *                     or on local points (local_var) if reverse
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *   tag           <-- tag for MPI_irecv
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

static void
_irecv_point_var_distant(fvmc_locator_t     *this_locator,
                         void              *var,
                         const fvmc_lnum_t  *local_list,
                         MPI_Datatype       datatype,
                         size_t             stride,
                         _Bool              reverse,
                         int                tag,
                         int                *request)
{
  int i, dist_v_count, loc_v_count, size;
  int dist_rank, dist_index;
  fvmc_lnum_t n_points_loc, n_points_dist;
  size_t dist_v_idx;
  void *dist_v_ptr;

  double comm_timing[4] = {0., 0., 0., 0.};

  MPI_Aint extent;

  /* Get a free  request */

  *request = _get_free_request(&(this_locator->nblockings_recv),
                               &(this_locator->max_nblockings_recv));

  _fvmc_locator_nblocking_t *nblocking_recv = this_locator->nblockings_recv + *request;

  BFTC_MALLOC(nblocking_recv->MPI_request, this_locator->n_intersects, MPI_Request);
  BFTC_MALLOC(nblocking_recv->buffer, this_locator->n_intersects, void *);
  nblocking_recv->var = var;
  nblocking_recv->local_list = local_list;
  nblocking_recv->reverse = reverse;
  nblocking_recv->stride = stride;

  /* Check extent of datatype */

	MPI_Aint lb;
  MPI_Type_get_extent(datatype, &lb, &extent);
  MPI_Type_size(datatype, &size);

  nblocking_recv->size = size;

  if (extent != size)
    bftc_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Loop on possibly intersecting distant ranks */
  /*---------------------------------------------*/

  for (i = 0; i < this_locator->n_intersects; i++) {

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    if (reverse == false) {

      /* Initialization */


      n_points_loc =    this_locator->local_points_idx[i+1]
                      - this_locator->local_points_idx[i];

      loc_v_count = n_points_loc*stride;

      BFTC_MALLOC(nblocking_recv->buffer[i], n_points_loc*size*stride, char);

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Irecv(nblocking_recv->buffer[i], loc_v_count, datatype, dist_rank, tag,
                this_locator->comm, nblocking_recv->MPI_request + i);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    }

    else { /* if (reverse == true) */

      void *distant_var = var;

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
      dist_v_count = n_points_dist * stride;

      loc_v_count = n_points_loc*stride;

      /* Exchange information */

      if (distant_var != NULL)
        dist_v_ptr = (void *)(((char *)distant_var) + dist_v_idx);
      else
        dist_v_ptr = NULL;

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Irecv(dist_v_ptr, dist_v_count, datatype, dist_rank, tag,
                this_locator->comm, nblocking_recv->MPI_request + i);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    }

  } /* End of loop on possibly intersecting ranks */

  this_locator->irecv_wtime[1] += comm_timing[0];
  this_locator->irecv_cpu_time[1] += comm_timing[1];
}


/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_distant(fvmc_locator_t     *this_locator,
                            void              *distant_var,
                            void              *local_var,
                            const fvmc_lnum_t  *local_list,
                            MPI_Datatype       datatype,
                            size_t             stride,
                            _Bool              reverse)
{
  int i, dist_v_count, loc_v_count, size;
  int dist_rank, dist_index;
  int dist_v_flag, loc_v_flag;
  fvmc_lnum_t n_points_loc, n_points_loc_max, n_points_dist;
  size_t dist_v_idx;
  void *dist_v_ptr;
  void *loc_v_buf;

  double comm_timing[4] = {0., 0., 0., 0.};

  MPI_Aint extent;
  MPI_Status status;

  /* Check extent of datatype */

	MPI_Aint lb;
  MPI_Type_get_extent(datatype, &lb, &extent);
  MPI_Type_size(datatype, &size);

  if (extent != size)
    bftc_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Initialization */

  n_points_loc_max = 0;

  for (i = 0; i < this_locator->n_intersects; i++) {
    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];
    if (n_points_loc > n_points_loc_max)
      n_points_loc_max = n_points_loc;
  }

  BFTC_MALLOC(loc_v_buf, n_points_loc_max*size*stride, char);

  /* Loop on possibly intersecting distant ranks */
  /*---------------------------------------------*/

  for (i = 0; i < this_locator->n_intersects; i++) {

    const fvmc_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank = this_locator->intersect_rank[dist_index];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (distant_var != NULL && n_points_dist > 0)
      dist_v_flag = 1;
    else
      dist_v_flag = 0;

    _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&dist_v_flag, 1, MPI_INT, dist_rank, FVMC_MPI_TAG,
                 &loc_v_flag, 1, MPI_INT, dist_rank, FVMC_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    if (loc_v_flag == 1 && (local_var == NULL || n_points_loc == 0))
      bftc_error(__FILE__, __LINE__, 0,
                _("Incoherent arguments to different instances in "
                  "_exchange_point_var().\n"
                  "Send and receive operations do not match "
                  "(dist_rank = %d\n)\n"), dist_rank);

    dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
    dist_v_count = n_points_dist * stride * dist_v_flag;

    if (loc_v_flag > 0)
      loc_v_count = n_points_loc*stride;
    else
      loc_v_count = 0;

    /* Exchange information */

    if (distant_var != NULL)
      dist_v_ptr = (void *)(((char *)distant_var) + dist_v_idx);
    else
      dist_v_ptr = NULL;

    if (reverse == false) {

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(dist_v_ptr, dist_v_count, datatype, dist_rank, FVMC_MPI_TAG,
                   loc_v_buf, loc_v_count, datatype, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p = (char *)local_var + _local_point_ids[k]*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            char *local_v_p =   (char *)local_var
                              + (local_list[_local_point_ids[k]] - 1)*nbytes;
            const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              local_v_p[l] = loc_v_buf_p[l];
          }
        }
      }

    }
    else { /* if (reverse == true) */

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var + _local_point_ids[k]*nbytes;
            char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
        else {
          int k;
          size_t l;
          const size_t nbytes = stride*size;
          for (k = 0; k < n_points_loc; k++) {
            const char *local_v_p
              = (const char *)local_var
                + (local_list[_local_point_ids[k]] - 1)*nbytes;
            char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
            for (l = 0; l < nbytes; l++)
              loc_v_buf_p[l] = local_v_p[l];
          }
        }
      }

      _locator_trace_start_comm(_fvmc_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(loc_v_buf, loc_v_count, datatype, dist_rank, FVMC_MPI_TAG,
                   dist_v_ptr, dist_v_count, datatype, dist_rank, FVMC_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_fvmc_locator_log_end_p_comm, comm_timing);

    }

  } /* End of loop on possibly intersecting ranks */

  BFTC_FREE(loc_v_buf);

  this_locator->exchange_wtime[1] += comm_timing[0];
  this_locator->exchange_cpu_time[1] += comm_timing[1];
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Distribute variable defined on "distant points" to the original ("local")
 * points.
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_local(fvmc_locator_t     *this_locator,
                          void              *distant_var,
                          void              *local_var,
                          const fvmc_lnum_t  *local_list,
                          size_t             type_size,
                          size_t             stride,
                          _Bool              reverse)
{
  fvmc_lnum_t i;
  size_t j;
  fvmc_lnum_t n_points_loc;

  const size_t nbytes = stride*type_size;

  /* Initialization */

  if (this_locator->n_interior == 0)
    return;

  n_points_loc =   this_locator->local_points_idx[1]
                 - this_locator->local_points_idx[0];

  assert(n_points_loc == (  this_locator->distant_points_idx[1]
                          - this_locator->distant_points_idx[0]));

  /* Exchange information */

  if (reverse == false) {

    if (local_list == NULL)
      memcpy((void *) local_var, (const void *) distant_var, (size_t) n_points_loc*nbytes);

    else {
      for (i = 0; i < n_points_loc; i++) {
        char *local_var_p = (char *)local_var + (local_list[i] - 1)*nbytes;
        const char *distant_var_p = (const char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          local_var_p[j] = distant_var_p[j];
      }
    }

  }
  else { /* if (reverse == true) */

    if (local_list == NULL)
      memcpy((void *) distant_var, (const void *) local_var, (size_t) n_points_loc*nbytes);

    else {
      for (i = 0; i < n_points_loc; i++) {
        const char *local_var_p
          = (const char *)local_var + (local_list[i] - 1)*nbytes;
        char *distant_var_p = (char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          distant_var_p[j] = local_var_p[j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   time_type         <-- 0 for total times, 1 for communication times
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *   issend_wtime    --> Variable issend Wall-clock time (or NULL)
 *   issend_cpu_time --> Variable issend CPU time (or NULL)
 *   irecv_wtime    --> Variable irecv Wall-clock time (or NULL)
 *   irecv_cpu_time --> Variable irecv CPU time (or NULL)
 *----------------------------------------------------------------------------*/

static void
_get_times(const fvmc_locator_t  *this_locator,
           int                   time_type,
           double               *location_wtime,
           double               *location_cpu_time,
           double               *exchange_wtime,
           double               *exchange_cpu_time,
           double               *issend_wtime,
           double               *issend_cpu_time,
           double               *irecv_wtime,
           double               *irecv_cpu_time)

{
  const fvmc_locator_t  *_locator = this_locator;

  if (this_locator != NULL) {

    if (location_wtime != NULL)
      *location_wtime = _locator->location_wtime[time_type];
    if (location_cpu_time != NULL)
      *location_cpu_time = _locator->location_cpu_time[time_type];

    if (exchange_wtime != NULL)
      *exchange_wtime = _locator->exchange_wtime[time_type];
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = _locator->exchange_cpu_time[time_type];

    if (issend_wtime != NULL)
      *issend_wtime = _locator->issend_wtime[time_type];
    if (issend_cpu_time != NULL)
      *issend_cpu_time = _locator->issend_cpu_time[time_type];

    if (irecv_wtime != NULL)
      *irecv_wtime = _locator->irecv_wtime[time_type];
    if (irecv_cpu_time != NULL)
      *irecv_cpu_time = _locator->irecv_cpu_time[time_type];

  }
  else {

    if (location_wtime != NULL)
      *location_wtime = 0.;
    if (location_cpu_time != NULL)
      *location_cpu_time = 0.;

    if (exchange_wtime != NULL)
      *exchange_wtime = 0.;
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = 0.;

    if (issend_wtime != NULL)
      *issend_wtime = 0.;
    if (issend_cpu_time != NULL)
      *issend_cpu_time = 0.;

    if (irecv_wtime != NULL)
      *irecv_wtime = 0.;
    if (irecv_cpu_time != NULL)
      *irecv_cpu_time = 0.;

  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/
/*----------------------------------------------------------------------------
 * Creation of a locator structure.
 *
 * Note that depending on the choice of ranks of the associated communicator,
 * distant ranks may in fact be truly distant or not. If n_ranks = 1 and
 * start_rank is equal to the current rank in the communicator, the locator
 * will work only locally.
 *
 * parameters:
 *   opt_bbox_step <-- Discretization for the computation of the
 *                     ho element extents
 *                  extent = base_extent * (1 + tolerance)
 *   tolerance  <-- addition to local extents of each element:
 *                  extent = base_extent * (1 + tolerance)
 *   comm       <-- associated MPI communicator
 *   n_ranks    <-- number of MPI ranks associated with distant location
 *   start_rank <-- first MPI rank associated with distant location
 *
 * returns:
 *   pointer to locator
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)

fvmc_locator_t *
fvmc_locator_create(int       opt_bbox_step,
                    double    tolerance,
                    MPI_Comm  comm,
                    int       n_ranks,
                    int       start_rank)
#else
fvmc_locator_t *
fvmc_locator_create(double  tolerance)
#endif
{
  int  i;
  fvmc_locator_t  *this_locator;

  BFTC_MALLOC(this_locator, 1, fvmc_locator_t);

  this_locator->tolerance = tolerance;
  this_locator->dim = 0;

#if defined(FVMC_HAVE_MPI)
  this_locator->comm = comm;
  this_locator->opt_bbox_step = opt_bbox_step;
  this_locator->n_ranks = n_ranks;
  this_locator->start_rank = start_rank;
  this_locator->nblockings_send = NULL;
  this_locator->nblockings_recv = NULL;
  this_locator->max_nblockings_send = 0;
  this_locator->max_nblockings_recv = 0;
#else
  this_locator->n_ranks = 1;
  this_locator->start_rank = 0;
#endif

  this_locator->n_intersects = 0;
  this_locator->intersect_rank = NULL;
  this_locator->intersect_extents = NULL;

  this_locator->local_points_idx = NULL;
  this_locator->distant_points_idx = NULL;

  this_locator->local_point_ids = NULL;

  this_locator->distant_point_location = NULL;
  this_locator->distant_point_coords = NULL;
  this_locator->distant_point_projected_coords = NULL;
  this_locator->distant_point_uvw = NULL;
  this_locator->distant_point_distance = NULL;

  this_locator->n_interior = 0;
  this_locator->interior_list = NULL;

  this_locator->n_exterior = 0;
  this_locator->exterior_list = NULL;
  this_locator->distant_distribution = NULL;
  this_locator->local_distribution = NULL;
  this_locator->distant_point_distance = NULL;

  for (i = 0; i < 2; i++) {
    this_locator->location_wtime[i] = 0.;
    this_locator->location_cpu_time[i] = 0.;
    this_locator->exchange_wtime[i] = 0.;
    this_locator->exchange_cpu_time[i] = 0.;
    this_locator->issend_wtime[i] = 0.;
    this_locator->issend_cpu_time[i] = 0.;
    this_locator->irecv_wtime[i] = 0.;
    this_locator->irecv_cpu_time[i] = 0.;
  }

  return this_locator;
}

/*----------------------------------------------------------------------------
 * Destruction of a locator structure.
 *
 * parameters:
 *   this_locator <-> locator to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_locator_t *
fvmc_locator_destroy(fvmc_locator_t  * this_locator)
{
  if (this_locator != NULL) {

    BFTC_FREE(this_locator->local_points_idx);
    BFTC_FREE(this_locator->distant_points_idx);

    if (this_locator->local_point_ids != NULL)
      BFTC_FREE(this_locator->local_point_ids);

    BFTC_FREE(this_locator->distant_point_location);
    BFTC_FREE(this_locator->distant_point_distance);
    BFTC_FREE(this_locator->distant_point_coords);

    BFTC_FREE(this_locator->intersect_rank);
    BFTC_FREE(this_locator->intersect_extents);

    BFTC_FREE(this_locator->interior_list);
    BFTC_FREE(this_locator->exterior_list);

    if (this_locator->nblockings_send != NULL)
      BFTC_FREE(this_locator->nblockings_send);

    if (this_locator->nblockings_recv != NULL)
      BFTC_FREE(this_locator->nblockings_recv);

    if (this_locator->distant_distribution != NULL)
      BFTC_FREE(this_locator->distant_distribution);
    if (this_locator->local_distribution != NULL)
      BFTC_FREE(this_locator->local_distribution);
    if (this_locator->distant_point_distance != NULL)
      BFTC_FREE(this_locator->distant_point_distance);

    if (this_locator->distant_point_projected_coords != NULL)
      BFTC_FREE(this_locator->distant_point_projected_coords);

    if (this_locator->distant_point_uvw != NULL)
      BFTC_FREE(this_locator->distant_point_uvw);

    BFTC_FREE(this_locator);
  }

  return NULL;
}
/*----------------------------------------------------------------------------
 * locator size
 *
 * parameters:
 *   this_locator <-> locator to get size
 *
 *----------------------------------------------------------------------------*/

size_t
fvmc_locator_size(const fvmc_locator_t  * this_locator)
{
  size_t il_size = 0;
  if (this_locator != NULL) {
     il_size += sizeof(double);
     il_size += sizeof(_Bool);
     il_size += 4 * sizeof(int);
     il_size += this_locator->n_intersects*sizeof(int);
#if defined(FVMC_HAVE_MPI)
     il_size += 2 * sizeof(int);
#endif
     il_size += this_locator->n_intersects * this_locator->dim * 2 * sizeof(double);
     il_size += (this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t);
     il_size += (this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t);
     il_size += (this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t);
     il_size += (this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t);
     il_size +=  this_locator->local_points_idx[this_locator->n_intersects] * sizeof(fvmc_lnum_t);
     il_size +=  this_locator->distant_points_idx[this_locator->n_intersects] * sizeof(float);
     il_size +=  this_locator->distant_points_idx[this_locator->n_intersects] * sizeof(fvmc_lnum_t);
     il_size +=  this_locator->distant_points_idx[this_locator->n_intersects] *this_locator->dim  * sizeof(fvmc_coord_t);
     il_size +=  sizeof(fvmc_lnum_t);
     il_size +=  this_locator->n_interior *  sizeof(fvmc_lnum_t);
     il_size +=  sizeof(fvmc_lnum_t);
     il_size +=  this_locator->n_exterior *  sizeof(fvmc_lnum_t);
  }
  return il_size;
}

/*----------------------------------------------------------------------------
 * save a locator
 *
 * parameters:
 *   this_locator <-> locator to save
 *
 *----------------------------------------------------------------------------*/


void *
fvmc_locator_pack(void *p, const fvmc_locator_t  * this_locator)
{
  size_t s_pack;
  if (this_locator != NULL) {

    s_pack = sizeof(double);
    memcpy(p,(const void *)&this_locator->tolerance, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(_Bool);
    memcpy(p,(const void *)&this_locator->locate_on_parents, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->dim,  s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->n_ranks, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->start_rank,  s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->n_intersects,  s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->n_intersects*sizeof(int);
    memcpy(p,(const void *)this_locator->intersect_rank, s_pack);
    p = (void *) ((char *) p + s_pack);
#if defined(FVMC_HAVE_MPI)

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->max_nblockings_send, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(int);
    memcpy(p,(const void *)&this_locator->max_nblockings_recv, s_pack);
    p = (void *) ((char *) p + s_pack);
#endif

    s_pack = this_locator->n_intersects * this_locator->dim * 2 * sizeof(double);
    memcpy(p,(void *)this_locator->intersect_extents, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = (this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->local_points_idx, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = (this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->local_distribution, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = (this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->distant_points_idx, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = (this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->distant_distribution, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->local_points_idx[this_locator->n_intersects] * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->local_point_ids, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->distant_points_idx[this_locator->n_intersects] * sizeof(float);
    memcpy(p,(void *)this_locator->distant_point_distance, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->distant_points_idx[this_locator->n_intersects] * sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->distant_point_location, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->distant_points_idx[this_locator->n_intersects] *this_locator->dim  * sizeof(fvmc_coord_t);
    memcpy(p,(void *)this_locator->distant_point_coords, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(fvmc_lnum_t);
    memcpy(p,(const void *)&this_locator->n_interior, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = this_locator->n_interior *  sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->interior_list, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack = sizeof(fvmc_lnum_t);
    memcpy(p,(const void *)&this_locator->n_exterior, s_pack);
    p = (void *) ((char *) p + s_pack);

    s_pack =  this_locator->n_exterior *  sizeof(fvmc_lnum_t);
    memcpy(p,(void *)this_locator->exterior_list, s_pack);
    p = (void *) ((char *) p + s_pack);
  }
  return p;
}

/*----------------------------------------------------------------------------
 * unpack a locator
 *
 * parameters:
 *   this_locator <-> locator to read
 *
 *----------------------------------------------------------------------------*/
/* fonction de base aussi appele dans cwipi */
size_t fvmc_locator_unpack_elem(const void * buffer, void *data,  const size_t data_size)
{
  memcpy(data, buffer, data_size);
  return data_size;
}

size_t
fvmc_locator_unpack(unsigned char *buff, fvmc_locator_t  * this_locator)
{
  size_t cur_pos;
  cur_pos = 0;

  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->tolerance, sizeof(double));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->locate_on_parents, sizeof(_Bool));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->dim, sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->n_ranks, sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->start_rank, sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->n_intersects, sizeof(int));
  BFTC_MALLOC(this_locator->intersect_rank,
	      this_locator->n_intersects,
	      int);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->intersect_rank,this_locator->n_intersects*sizeof(int));
#if defined(FVMC_HAVE_MPI)
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->max_nblockings_send,sizeof(int));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->max_nblockings_recv,sizeof(int));
#endif
  BFTC_MALLOC(this_locator->intersect_extents,
	      this_locator->n_intersects * this_locator->dim * 2 ,
	      double);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->intersect_extents,this_locator->n_intersects * this_locator->dim * 2 * sizeof(double));

  BFTC_MALLOC(this_locator->local_points_idx,
	      this_locator->n_intersects + 1,
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->local_points_idx,(this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->local_distribution,
	      this_locator->n_ranks + 1,
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->local_distribution,(this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->distant_points_idx,
	      this_locator->n_intersects + 1,
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->distant_points_idx,(this_locator->n_intersects + 1) * sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->distant_distribution,
	      this_locator->n_ranks + 1,
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->distant_distribution,(this_locator->n_ranks + 1) * sizeof(fvmc_lnum_t));

  BFTC_MALLOC(this_locator->local_point_ids,
	      this_locator->local_points_idx[this_locator->n_intersects],
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->local_point_ids, this_locator->local_points_idx[this_locator->n_intersects] * sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->distant_point_distance,
	      this_locator->distant_points_idx[this_locator->n_intersects],
	      float);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->distant_point_distance, this_locator->distant_points_idx[this_locator->n_intersects]* sizeof(float));
  BFTC_MALLOC(this_locator->distant_point_location,
	      this_locator->distant_points_idx[this_locator->n_intersects],
	      fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->distant_point_location, this_locator->distant_points_idx[this_locator->n_intersects]* sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->distant_point_coords,
	      this_locator->distant_points_idx[this_locator->n_intersects] *this_locator->dim,
	      fvmc_coord_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->distant_point_coords,this_locator->distant_points_idx[this_locator->n_intersects] *this_locator->dim * sizeof(fvmc_coord_t));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->n_interior, sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->interior_list, this_locator->n_interior, fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->interior_list, this_locator->n_interior *  sizeof(fvmc_lnum_t));
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)&this_locator->n_exterior, sizeof(fvmc_lnum_t));
  BFTC_MALLOC(this_locator->exterior_list, this_locator->n_exterior, fvmc_lnum_t);
  cur_pos += fvmc_locator_unpack_elem((void *)&buff[cur_pos], (void *)this_locator->exterior_list, this_locator->n_exterior *  sizeof(fvmc_lnum_t));
  return cur_pos;
}



/*----------------------------------------------------------------------------
 * Prepare locator for use with a given nodal mesh representation.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   this_nodal        <-- pointer to mesh representation structure
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         1, id of element + 1 in concatenated sections
 *                         of same element dimension if 0
 *   dim               <-- space dimension of points to locate
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_locator_set_nodal(fvmc_locator_t       *this_locator,
                      const fvmc_nodal_t   *this_nodal,
                      int                  locate_on_parents,
                      int                  dim,
                      fvmc_lnum_t           n_points,
                      const fvmc_lnum_t     point_list[],
                      const fvmc_coord_t    point_coords[])
{
  int i;
  int stride2;
  double tolerance;
  double w_start, w_end, cpu_start, cpu_end;
  double extents[12];

#if defined(FVMC_HAVE_MPI)
  int j;
  int stride4;
  int comm_rank, comm_size;
  int locdim[2], globdim[2];
  int n_intersects;
  int  *intersect_rank;
  double *recvbuf;
#endif

  double comm_timing[4] = {0., 0., 0., 0.};
  int mpi_flag = 0;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  /* Other intializations */

  this_locator->locate_on_parents = (_Bool)locate_on_parents;

  this_locator->dim = dim;

  tolerance = FVMC_MAX(this_locator->tolerance, 0.1);

  _nodal_extents(this_nodal,
                 this_locator->opt_bbox_step,
                 tolerance,
                 extents);

  _point_extents(dim,
                 n_points,
                 point_list,
                 point_coords,
                 extents + 2*dim);


  /* Release information if previously present */

  if (this_locator->intersect_rank != NULL)
    BFTC_FREE(this_locator->intersect_rank);
  if (this_locator->intersect_extents != NULL)
    BFTC_FREE(this_locator->intersect_extents);

  if (this_locator->local_points_idx != NULL)
    BFTC_FREE(this_locator->local_points_idx);
  if (this_locator->distant_points_idx != NULL)
    BFTC_FREE(this_locator->distant_points_idx);

  if (this_locator->local_point_ids != NULL)
    BFTC_FREE(this_locator->local_point_ids);

  if (this_locator->distant_point_location != NULL)
    BFTC_FREE(this_locator->distant_point_location);
  if (this_locator->distant_point_coords != NULL)
    BFTC_FREE(this_locator->distant_point_coords);
  if (this_locator->distant_point_projected_coords != NULL)
    BFTC_FREE(this_locator->distant_point_projected_coords);
  if (this_locator->distant_point_uvw != NULL)
    BFTC_FREE(this_locator->distant_point_uvw);

  if (this_locator->interior_list != NULL)
    BFTC_FREE(this_locator->interior_list);
  if (this_locator->exterior_list != NULL)
    BFTC_FREE(this_locator->exterior_list);

  if (this_locator->distant_distribution != NULL)
    BFTC_FREE(this_locator->distant_distribution);
  if (this_locator->local_distribution != NULL)
    BFTC_FREE(this_locator->local_distribution);
  if (this_locator->distant_point_distance != NULL)
    BFTC_FREE(this_locator->distant_point_distance);


  /* Prepare locator (MPI version) */
  /*-------------------------------*/

#if defined(FVMC_HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    /* Check that at least one of the local or distant nodal meshes
       is non-NULL, and at least one of the local or distant
       point sets is non null */

    if (this_nodal != NULL)
      locdim[0] = this_nodal->dim;
    else
      locdim[0] = -1;
    if (n_points > 0)
      locdim[1] = dim;
    else
      locdim[1] = -1;

    _locator_trace_start_comm(_fvmc_locator_log_start_g_comm, comm_timing);

    MPI_Allreduce(locdim, globdim, 2, MPI_INT, MPI_MAX,
                  this_locator->comm);

    _locator_trace_end_comm(_fvmc_locator_log_end_g_comm, comm_timing);


    if (globdim[0] < 0 || globdim[1] < 0)
      return;
    else if (this_nodal != NULL && globdim[1] != this_nodal->dim)
      bftc_error(__FILE__, __LINE__, 0,
                _("Locator trying to use distant space dimension %d\n"
                  "with local mesh dimension %d\n"),
                globdim[1], this_nodal->dim);
    else if (this_nodal == NULL && globdim[0] != dim)
      bftc_error(__FILE__, __LINE__, 0,
                _("Locator trying to use local space dimension %d\n"
                  "with distant mesh dimension %d\n"),
                dim, globdim[0]);

    /* Exchange extent information */

    MPI_Comm_rank(this_locator->comm, &comm_rank);
    MPI_Comm_size(this_locator->comm, &comm_size);

    stride2 = dim * 2; /* Stride for one type of extent */
    stride4 = dim * 4; /* Stride for element and vertex
                          extents, end-to-end */

    BFTC_MALLOC(recvbuf, stride4*comm_size, double);

    _locator_trace_start_comm(_fvmc_locator_log_start_g_comm, comm_timing);

    MPI_Allgather(extents, stride4, MPI_DOUBLE, recvbuf, stride4, MPI_DOUBLE,
                  this_locator->comm);

    _locator_trace_end_comm(_fvmc_locator_log_end_g_comm, comm_timing);

    /* Count and mark possible overlaps */

    n_intersects = 0;
    BFTC_MALLOC(intersect_rank, this_locator->n_ranks, int);

    for (i = 0; i < this_locator->n_ranks; i++) {
      j = this_locator->start_rank + i;
      if (  (_intersect_extents(dim,
                                 extents + (2*dim),
                                 recvbuf + (j*stride4)) == true)
          || (_intersect_extents(dim,
                                 extents,
                                 recvbuf + (j*stride4) + (2*dim)) == true)) {
        intersect_rank[n_intersects] = j;
        n_intersects += 1;
      }

    }

    this_locator->n_intersects = n_intersects;
    BFTC_MALLOC(this_locator->intersect_rank,
               this_locator->n_intersects,
               int);
    BFTC_MALLOC(this_locator->intersect_extents,
               this_locator->n_intersects * stride2,
               double);

    for (i = 0; i < this_locator->n_intersects; i++) {

      this_locator->intersect_rank[i] = intersect_rank[i];

      /* Copy only distant element (and not point) extents */

      for (j = 0; j < stride2; j++)
        this_locator->intersect_extents[i*stride2 + j]
          = recvbuf[intersect_rank[i]*stride4 + j];

    }

    /* Free temporary memory */

    BFTC_FREE(intersect_rank);
    BFTC_FREE(recvbuf);

    _locate_all_distant(this_locator,
                        this_nodal,
                        dim,
                        n_points,
                        point_list,
                        point_coords);

  }

#endif

  /* Prepare locator (local version) */
  /*---------------------------------*/

  if (!mpi_flag) {

    if (this_nodal == NULL || n_points == 0)
      return;

    stride2 = this_nodal->dim * 2;

    /* Count and mark possible overlaps */

    if (_intersect_extents(this_nodal->dim,
                           extents,
                           extents + (2*dim)) == true) {

      this_locator->n_intersects = 1;

      BFTC_MALLOC(this_locator->intersect_rank, 1, int);
      BFTC_MALLOC(this_locator->intersect_extents, stride2, double);

      this_locator->intersect_rank[0] = 0;

      for (i = 0; i < stride2; i++)
        this_locator->intersect_extents[i] = extents[i];

      _locate_all_local(this_locator,
                        this_nodal,
                        dim,
                        n_points,
                        point_list,
                        point_coords);
    }

    else {

      this_locator->n_exterior = n_points;
      BFTC_MALLOC(this_locator->exterior_list,
                 this_locator->n_exterior,
                 fvmc_lnum_t);
      for (i = 0; i < this_locator->n_exterior; i++)
        this_locator->exterior_list[i] = i + 1;

    }

  }

  /* Update local_point_ids values */
  /*-------------------------------*/

  if (   this_locator->n_interior > 0
      && this_locator->local_point_ids != NULL) {

    fvmc_lnum_t  *reduced_index;

    BFTC_MALLOC(reduced_index, n_points, fvmc_lnum_t);

    for (i = 0; i < n_points; i++)
      reduced_index[i] = -1;

    assert(  this_locator->local_points_idx[this_locator->n_intersects]
           == this_locator->n_interior);

    for (i = 0; i < this_locator->n_interior; i++)
      reduced_index[this_locator->interior_list[i] - 1] = i;

    /* Update this_locator->local_point_ids[] so that it refers
       to an index in a dense [0, this_locator->n_interior] subset
       of the local points */

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->local_point_ids[i]
        = reduced_index[this_locator->local_point_ids[i]];

    BFTC_FREE(reduced_index);

  }

  /* If an initial point list was given, update
     this_locator->interior_list and this_locator->exterior_list
     so that they refer to the same point set as that initial
     list (and not to an index within the selected point set) */

  if (point_list != NULL) {

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->interior_list[i]
        = point_list[this_locator->interior_list[i] - 1];

    for (i = 0; i < this_locator->n_exterior; i++)
      this_locator->exterior_list[i]
        = point_list[this_locator->exterior_list[i] - 1];

  }

  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->location_wtime[0] += (w_end - w_start);
  this_locator->location_cpu_time[0] += (cpu_end - cpu_start);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}


/*----------------------------------------------------------------------------
 * Return distribution graph of distant points per distant rank
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   distant points index
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_locator_get_dist_distrib(const fvmc_locator_t  *this_locator)
{
  fvmc_lnum_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_distribution;
  }

  return retval;
}


/*----------------------------------------------------------------------------
 * Return distribution graph of local points per distant rank
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   distant points index
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_locator_get_loc_distrib(const fvmc_locator_t  *this_locator)
{
  fvmc_lnum_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->local_distribution;
  }

  return retval;
}


/*----------------------------------------------------------------------------
 * Return number of distant points after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of distant points.
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_locator_get_n_dist_points(const fvmc_locator_t  *this_locator)
{
  fvmc_lnum_t retval = 0;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_points_idx[this_locator->n_intersects];
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return an array of local element numbers containing (or nearest to)
 * each distant point after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   local element numbers associated with distant points (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_locator_get_dist_locations(const fvmc_locator_t  *this_locator)
{
  const fvmc_lnum_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_ranks != 0)
      retval = this_locator->distant_point_location;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return an array of distance to local element numbers containing
 * each distant point after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   distance to local element numbers associated with distant points
 *----------------------------------------------------------------------------*/

const float *
fvmc_locator_get_dist_distances(const fvmc_locator_t  *this_locator)
{
  const float * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_ranks != 0)
      retval = (float *) this_locator->distant_point_distance;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return an array of coordinates of each distant point after
 * locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   coordinate array associated with distant points (interlaced).
 *----------------------------------------------------------------------------*/

const fvmc_coord_t *
fvmc_locator_get_dist_coords(const fvmc_locator_t  *this_locator)
{
  const fvmc_coord_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_coords;
  }

  return retval;
}


/*----------------------------------------------------------------------------
 * Return an array of coordinates of each distant point projected on the closest element.
 * (or NULL), available for high order nodal
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   coordinate array associated with distant points (interlaced).
 *----------------------------------------------------------------------------*/

const fvmc_coord_t *
fvmc_locator_get_dist_projected_coords(const fvmc_locator_t  *this_locator)
{
  const fvmc_coord_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_projected_coords;
  }

  return retval;

}


/*----------------------------------------------------------------------------
 * Return an array of uvw of each distant point in the closest element.
 * (or NULL), available for high order nodal
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   uvw (size = 3 * n_dist_point, interlaced)
 *----------------------------------------------------------------------------*/

const double *
fvmc_locator_get_dist_uvw(const fvmc_locator_t  *this_locator)
{
  const double * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_uvw;
  }

  return retval;

}

/*----------------------------------------------------------------------------
 * Return number of points located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points located.
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_locator_get_n_interior(const fvmc_locator_t  *this_locator)
{
  fvmc_lnum_t retval = 0;

  if (this_locator != NULL)
    retval = this_locator->n_interior;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return list of points located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points located (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_locator_get_interior_list(const fvmc_locator_t  *this_locator)
{
  return this_locator->interior_list;
}

/*----------------------------------------------------------------------------
 * Return number of points not located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points not located.
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_locator_get_n_exterior(const fvmc_locator_t  *this_locator)
{
  return this_locator->n_exterior;
}

/*----------------------------------------------------------------------------
 * Return list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points not located (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_locator_get_exterior_list(const fvmc_locator_t  *this_locator)
{
  return this_locator->exterior_list;
}

/*----------------------------------------------------------------------------
 * Discard list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points not located (1 to n numbering).
 *----------------------------------------------------------------------------*/

void
fvmc_locator_discard_exterior(fvmc_locator_t  *this_locator)
{
  this_locator->n_exterior = 0;
  BFTC_FREE(this_locator->exterior_list);
}


#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Non blocking send dist_var or (local_var if reverse == true)
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           <-- variable defined on distant points (distant_var)
 *                     or variable defined on local points (local_var)
 *                     size: n_dist_points*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *   tag           <-- tag for MPI_issend
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvmc_locator_issend_point_var(fvmc_locator_t     *this_locator,
                             void              *var,
                             const fvmc_lnum_t  *local_list,
                             size_t             type_size,
                             size_t             stride,
                             int                reverse,
                             int                tag,
                             int               *request)
{
  double w_start, w_end, cpu_start, cpu_end;

  _Bool _reverse = reverse;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  MPI_Datatype datatype = MPI_DATATYPE_NULL;

  if (type_size == sizeof(double))
    datatype = MPI_DOUBLE;
  else if (type_size == sizeof(float))
    datatype = MPI_FLOAT;
  else
    bftc_error(__FILE__, __LINE__, 0,
              _("type_size passed to fvmc_locator_issend_point_var() does\n"
                "not correspond to double or float."));

  assert (datatype != MPI_DATATYPE_NULL);

  _issend_point_var_distant(this_locator,
                            var,
                            local_list,
                            datatype,
                            stride,
                            _reverse,
                            tag,
                            request);

  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->issend_wtime[0] += (w_end - w_start);
  this_locator->issend_cpu_time[0] += (cpu_end - cpu_start);
}


/*----------------------------------------------------------------------------
 * Wait for fvmc_locator_iisend
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvmc_locator_issend_wait(fvmc_locator_t     *this_locator,
                        int               request)
{

  MPI_Status status;

  double w_start, w_end, cpu_start, cpu_end;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  _fvmc_locator_nblocking_t *nblocking_send = this_locator->nblockings_send + request;

  /* MPI_Wait for MPI_issend requests */

  for (int i = 0; i < this_locator->n_intersects; i++) {
    MPI_Wait(nblocking_send->MPI_request + i, &status);
    if (nblocking_send->reverse)
      BFTC_FREE(nblocking_send->buffer[i]);
  }

  /* Free request */

  BFTC_FREE(nblocking_send->MPI_request);

  if (nblocking_send->reverse)
    BFTC_FREE(nblocking_send->buffer);

  nblocking_send->MPI_request = NULL;
  nblocking_send->buffer = NULL;
  nblocking_send->local_list = NULL;
  nblocking_send->var = NULL;
  nblocking_send->reverse = false;
  nblocking_send->size = 0;
  nblocking_send->stride = 0;

  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->issend_wtime[0] += (w_end - w_start);
  this_locator->issend_cpu_time[0] += (cpu_end - cpu_start);

}


/*----------------------------------------------------------------------------
 * Non blocking receive dist_var or (local_var if reverse == true)
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           --> variable defined on distant points (distant_var)
 *                     or variable defined on local points (local_var)
 *                     size: n_dist_points*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *   tag           <-- tag for MPI_irecv
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvmc_locator_irecv_point_var(fvmc_locator_t     *this_locator,
                            void              *var,
                            const fvmc_lnum_t  *local_list,
                            size_t             type_size,
                            size_t             stride,
                            int                reverse,
                            int                tag,
                            int*               request)
{
  double w_start, w_end, cpu_start, cpu_end;

  _Bool _reverse = reverse;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  MPI_Datatype datatype = MPI_DATATYPE_NULL;

  if (type_size == sizeof(double))
    datatype = MPI_DOUBLE;
  else if (type_size == sizeof(float))
    datatype = MPI_FLOAT;
  else
    bftc_error(__FILE__, __LINE__, 0,
              _("type_size passed to fvmc_locator_issend_point_var() does\n"
                "not correspond to double or float."));

  assert (datatype != MPI_DATATYPE_NULL);

  _irecv_point_var_distant(this_locator,
                           var,
                           local_list,
                           datatype,
                           stride,
                           _reverse,
                           tag,
                           request);



  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->irecv_wtime[0] += (w_end - w_start);
  this_locator->irecv_cpu_time[0] += (cpu_end - cpu_start);
  }


/*----------------------------------------------------------------------------
 * Wait for fvmc_locator_irecv
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvmc_locator_irecv_wait(fvmc_locator_t     *this_locator,
                        int               request)
{

  MPI_Status status;

  double w_start, w_end, cpu_start, cpu_end;
  fvmc_lnum_t n_points_loc;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

  _fvmc_locator_nblocking_t *nblocking_recv = this_locator->nblockings_recv + request;


  /* MPI_Wait for MPI_irecv requests */

  for (int i = 0; i < this_locator->n_intersects; i++) {
    MPI_Wait(nblocking_recv->MPI_request + i, &status);

    const fvmc_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
                     - this_locator->local_points_idx[i];

    int size = nblocking_recv->size;
    int stride = nblocking_recv->stride;

    if (!nblocking_recv->reverse) {
      if (nblocking_recv->local_list == NULL) {
        int k;
        size_t l;
        const size_t nbytes = stride*size;
        for (k = 0; k < n_points_loc; k++) {
          char *local_v_p = (char *)nblocking_recv->var + _local_point_ids[k]*nbytes;
          const char *loc_v_buf_p = (const char *)nblocking_recv->buffer[i] + k*nbytes;
          for (l = 0; l < nbytes; l++)
            local_v_p[l] = loc_v_buf_p[l];
        }
      }
      else {
        int k;
        size_t l;
        const size_t nbytes = stride*size;
        for (k = 0; k < n_points_loc; k++) {
          char *local_v_p =   (char *)nblocking_recv->var
            + (nblocking_recv->local_list[_local_point_ids[k]] - 1)*nbytes;
          const char *loc_v_buf_p = (const char *)nblocking_recv->buffer[i] + k*nbytes;
          for (l = 0; l < nbytes; l++)
            local_v_p[l] = loc_v_buf_p[l];
        }
      }
      BFTC_FREE(nblocking_recv->buffer[i]);
    }
  }

  /* Free request */

  BFTC_FREE(nblocking_recv->MPI_request);
  BFTC_FREE(nblocking_recv->buffer);

  nblocking_recv->buffer = NULL;
  nblocking_recv->MPI_request = NULL;
  nblocking_recv->var = NULL;
  nblocking_recv->local_list = NULL;
  nblocking_recv->reverse = false;
  nblocking_recv->size = 0;
  nblocking_recv->stride = 0;

  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->irecv_wtime[0] += (w_end - w_start);
  this_locator->irecv_cpu_time[0] += (cpu_end - cpu_start);

}

#endif

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * The caller should have defined the values of distant_var[] for the
 * distant points, whose coordinates are given by
 * fvmc_locator_get_dist_coords(), and which are located in the elements
 * whose numbers are given by fvmc_locator_get_dist_locations().
 *
 * The local_var[] is defined at the located points (those whose
 * numbers are returned by fvmc_locator_get_interior_list().
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *                     size: n_dist_points*stride
 *   local_var     <-> variable defined on located local points (received)
 *                     size: n_interior*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

void
fvmc_locator_exchange_point_var(fvmc_locator_t     *this_locator,
                               void              *distant_var,
                               void              *local_var,
                               const fvmc_lnum_t  *local_list,
                               size_t             type_size,
                               size_t             stride,
                               int                reverse)
{
  double w_start, w_end, cpu_start, cpu_end;

  int mpi_flag = 0;
  _Bool _reverse = reverse;

  /* Initialize timing */

  w_start = bftc_timer_wtime();
  cpu_start = bftc_timer_cpu_time();

#if defined(FVMC_HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    MPI_Datatype datatype = MPI_DATATYPE_NULL;

    if (type_size == sizeof(double))
      datatype = MPI_DOUBLE;
    else if (type_size == sizeof(float))
      datatype = MPI_FLOAT;
    else
      bftc_error(__FILE__, __LINE__, 0,
                _("type_size passed to fvmc_locator_exchange_point_var() does\n"
                  "not correspond to double or float."));

    assert (datatype != MPI_DATATYPE_NULL);

    _exchange_point_var_distant(this_locator,
                                distant_var,
                                local_var,
                                local_list,
                                datatype,
                                stride,
                                _reverse);

  }

#endif /* defined(FVMC_HAVE_MPI) */

  if (!mpi_flag)
    _exchange_point_var_local(this_locator,
                              distant_var,
                              local_var,
                              local_list,
                              type_size,
                              stride,
                              _reverse);

  /* Finalize timing */

  w_end = bftc_timer_wtime();
  cpu_end = bftc_timer_cpu_time();

  this_locator->exchange_wtime[0] += (w_end - w_start);
  this_locator->exchange_cpu_time[0] += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * In parallel mode, this includes communication time.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *   issend_wtime      --> Variable exchange Wall-clock time (or NULL)
 *   issend_cpu_time   --> Variable exchange CPU time (or NULL)
 *   irecv_wtime       --> Variable exchange Wall-clock time (or NULL)
 *   irecv_cpu_time    --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvmc_locator_get_times(const fvmc_locator_t  *this_locator,
                      double               *location_wtime,
                      double               *location_cpu_time,
                      double               *exchange_wtime,
                      double               *exchange_cpu_time,
                      double               *issend_wtime,
                      double               *issend_cpu_time,
                      double               *irecv_wtime,
                      double               *irecv_cpu_time)
{
  _get_times(this_locator,
             0,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time,
             issend_wtime, issend_cpu_time,
             irecv_wtime, irecv_cpu_time);
}

/*----------------------------------------------------------------------------
 * Return communication timing information.
 *
 * In serial mode, returned times are always zero..
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *   issend_wtime      --> Variable exchange Wall-clock time (or NULL)
 *   issend_cpu_time   --> Variable exchange CPU time (or NULL)
 *   irecv_wtime       --> Variable exchange Wall-clock time (or NULL)
 *   irecv_cpu_time    --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvmc_locator_get_comm_times(const fvmc_locator_t  *this_locator,
                           double               *location_wtime,
                           double               *location_cpu_time,
                           double               *exchange_wtime,
                           double               *exchange_cpu_time,
                           double               *issend_wtime,
                           double               *issend_cpu_time,
                           double               *irecv_wtime,
                           double               *irecv_cpu_time)
{
  _get_times(this_locator,
             1,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time,
             issend_wtime, issend_cpu_time,
             irecv_wtime, irecv_cpu_time);
}

/*----------------------------------------------------------------------------
 * Dump printout of a locator structure.
 *
 * parameters:
 *   this_locator  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvmc_locator_dump(const fvmc_locator_t  *this_locator)
{
  int  i;
  fvmc_lnum_t  j, k;
  const fvmc_lnum_t  *idx, *index, *loc;
  const fvmc_coord_t  *coords;

  const fvmc_locator_t  *_locator = this_locator;

  if (this_locator == NULL)
    return;

  /* Basic information */
  /*-------------------*/

  bftc_printf("\n"
             "Locator:\n\n"
             "Tolerance:                             %f\n"
             "Locate on parents:                     %d\n"
             "Spatial dimension:                     %d\n"
             "Number of ranks of distant location:   %d\n"
             "First rank of distant location:        %d\n"
             "Number of intersecting distant ranks:  %d\n",
             _locator->tolerance,
             (int)_locator->locate_on_parents, _locator->dim,
             _locator->n_ranks, _locator->start_rank,
             _locator->n_intersects);

#if defined(FVMC_HAVE_MPI)
  if (_locator->comm != MPI_COMM_NULL)
    bftc_printf("\n"
               "Associated MPI communicator:           %ld\n",
               (long)(_locator->comm));
#endif

  /* Arrays indexed by rank */
  /*------------------------*/

  for (i = 0; i < _locator->n_intersects; i++) {

    bftc_printf("\n"
               "  Intersection %d with distant rank %d\n\n",
               i+1, _locator->intersect_rank[i]);

    bftc_printf("    Distant rank extents:\n");

    k = i * (_locator->dim) * 2;
    for (j = 0; j < _locator->dim; j++)
      bftc_printf("    [%12.5e, %12.5e]\n",
                 _locator->intersect_extents[k + j],
                 _locator->intersect_extents[k + _locator->dim + j]);

  }

  if (_locator->n_interior > 0) {

    if (_locator->local_point_ids != NULL) {

      bftc_printf("\n  Local point ids (for receiving):\n\n");
      idx = _locator->local_points_idx;
      index = _locator->local_point_ids;
      for (i = 0; i < _locator->n_intersects; i++) {
        if (idx[i+1] > idx[i]) {
          bftc_printf("%6d (idx = %10d) %10d\n",
                     i + 1, idx[i], index[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            bftc_printf("                          %10d\n", index[j]);
        }
        else {
          bftc_printf("%6d (idx = %10d)\n", i + 1, idx[i]);
        }
        bftc_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
      }

    }

  }

  if (_locator->distant_points_idx != NULL) {

    idx = _locator->distant_points_idx;
    loc = _locator->distant_point_location;
    coords = _locator->distant_point_coords;

    if (idx[_locator->n_intersects] > 0)
      bftc_printf("\n  Distant point location:\n\n");

    for (i = 0; i < _locator->n_intersects; i++) {

      if (idx[i+1] > idx[i]) {

        if (_locator->dim == 1) {
          bftc_printf("%6d (idx = %10d) %10d [%12.5e]\n",
                     i + 1, _locator->intersect_rank[i], idx[i],
                     loc[idx[i]], coords[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            bftc_printf("                          %10d [%12.5e]\n",
                       loc[j], coords[j]);
        }
        else if (_locator->dim == 2) {
          bftc_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e]\n",
                     i + 1, idx[i], loc[idx[i]],
                     coords[2*idx[i]], coords[2*idx[i]+1]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            bftc_printf("                          %10d [%12.5e, %12.5e]\n",
                       loc[j], coords[2*j], coords[2*j+1]);
        }
        else if (_locator->dim == 3) {
          bftc_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e, %12.5e]\n",
                     i + 1, idx[i], loc[idx[i]],
                     coords[3*idx[i]], coords[3*idx[i]+1], coords[3*idx[i]+2]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            bftc_printf("                          "
                       "%10d [%12.5e, %12.5e, %12.5e]\n",
                       loc[j], coords[3*j], coords[3*j+1], coords[3*j+2]);
        }

      } /* if (idx[i+1] > idx[i]) */

    }

    if (idx[_locator->n_intersects] > 0)
      bftc_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
  }

  /* Local arrays */
  /*--------------*/

  bftc_printf("\n"
             "  Number of local points successfully located:  %d\n\n",
             _locator->n_interior);

  for (j = 0; j < _locator->n_interior; j++)
    bftc_printf("    %10d\n", _locator->interior_list[j]);

  if  (_locator->n_interior > 0)
    bftc_printf("\n");

  bftc_printf("  Number of local points not located:  %d\n",
             _locator->n_exterior);

  for (j = 0; j < _locator->n_exterior; j++)
    bftc_printf("    %10d\n", _locator->exterior_list[j]);

  if  (_locator->n_exterior > 0)
    bftc_printf("\n");

  /* Timing information */
  /*--------------------*/

  bftc_printf("  Location Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->location_wtime[0], _locator->location_wtime[0]);

  bftc_printf("  Location CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->location_cpu_time[0], _locator->location_cpu_time[0]);

  bftc_printf("  Exchange Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->exchange_wtime[0], _locator->exchange_wtime[0]);

  bftc_printf("  Exchange CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->exchange_cpu_time[0], _locator->exchange_cpu_time[0]);

  bftc_printf("  Issend Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->issend_wtime[0], _locator->issend_wtime[0]);

  bftc_printf("  Issend CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->issend_cpu_time[0], _locator->issend_cpu_time[0]);

  bftc_printf("  Irecv Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->irecv_wtime[0], _locator->irecv_wtime[0]);

  bftc_printf("  Irecv CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->irecv_cpu_time[0], _locator->irecv_cpu_time[0]);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Register communication logging functions for locator instrumentation.
 *
 * By default, locators are not instrumented.
 *
 * Functions using MPE may be defined and used, but other similar systems
 * may be used.
 *
 * parameters:
 *   fct           <-- pointer to logging function
 *   start_p_comm  <-- point to point communication start event number
 *   end_p_comm    <-- point to point communication end event number
 *   start_g_comm  <-- global communication start event number
 *   end_g_comm    <-- global communication end event number
 *----------------------------------------------------------------------------*/

void
fvmc_locator_set_comm_log(fvmc_locator_log_t  *log_function,
                         int                 start_p_comm,
                         int                 end_p_comm,
                         int                 start_g_comm,
                         int                 end_g_comm)
{
  _fvmc_locator_log_func = log_function;

  _fvmc_locator_log_start_p_comm = start_p_comm;
  _fvmc_locator_log_end_p_comm = end_p_comm;
  _fvmc_locator_log_start_g_comm = start_g_comm;
  _fvmc_locator_log_end_g_comm = end_g_comm;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
