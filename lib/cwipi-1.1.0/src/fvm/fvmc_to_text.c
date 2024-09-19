/*============================================================================
 * Write a nodal representation associated with a mesh to file
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2006  EDF

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_file.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_gather.h"
#include "fvmc_io_num.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"
#include "fvmc_writer_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_to_text.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Text writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  bftc_file_t  *file;      /* Output file */

  fvmc_writer_time_dep_t   time_dependency; /* Mesh time dependency */

  int          rank;      /* Rank of current process in communicator */
  int          n_ranks;   /* Number of processes in communicator */

#if defined(FVMC_HAVE_MPI)
  MPI_Comm     comm;      /* Associated MPI communicator */
#endif

} fvmc_to_text_writer_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write slice of a vector of doubles to a text file
 *
 * parameters:
 *   stride         <-- number of values per element
 *   num_start      <-- global number of first element for this slice
 *   num_end        <-- global number of past the last element for this slice
 *   values         <-- pointer to values slice array
 *   f              <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_vector(const size_t             stride,
                    const fvmc_gnum_t         num_start,
                    const fvmc_gnum_t         num_end,
                    const double             values[],
                    bftc_file_t        *const f)
{
  size_t  i, k;
  fvmc_gnum_t  j;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  switch(stride) {

  case 1:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12.5f\n",
                      (unsigned long)j,
                      values[i]);
    break;

  case 2:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12.5f %12.5f\n",
                      (unsigned long)j,
                      values[2*i], values[2*i+1]);
    break;

  case 3:
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12.5f %12.5f %12.5f\n",
                      (unsigned long)j,
                      values[3*i], values[3*i+1], values[3*i+2]);
    break;

  default: /* Fallback, requiring more calls */
    for (i = 0, j = num_start ; j < num_end ; i++, j++) {
      bftc_file_printf(f, "%12lu :", (unsigned long)j);
      for (k = 0 ; k < stride ; k++)
        bftc_file_printf(f, " %12.5f",
                        values[stride*i+k]);
      bftc_file_printf(f, "\n");
    }
    break;
  }

}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write strided global connectivity slice to a text file
 *
 * parameters:
 *   stride           <-- number of vertices per element type
 *   num_start        <-- global number of first element for this slice
 *   num_end          <-- global number of last element for this slice
 *   global_connect_s <-- global connectivity slice array
 *   f                <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_slice_connect_g(const int                stride,
                       const fvmc_gnum_t         num_start,
                       const fvmc_gnum_t         num_end,
                       const fvmc_gnum_t         global_connect_s[],
                       bftc_file_t        *const f)
{
  fvmc_gnum_t  i, j;

  /* If called by non I/O rank, return */
  if (f == NULL)
    return;

  switch(stride) {

  case 2: /* edge */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*2],
                      (unsigned long)global_connect_s[i*2+1]);
    break;

  case 3: /* FVMC_FACE_TRIA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*3],
                      (unsigned long)global_connect_s[i*3+1],
                      (unsigned long)global_connect_s[i*3+2]);
    break;

  case 4: /* FVMC_FACE_QUAD or FVMC_CELL_TETRA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*4],
                      (unsigned long)global_connect_s[i*4+1],
                      (unsigned long)global_connect_s[i*4+2],
                      (unsigned long)global_connect_s[i*4+3]);
    break;

  case 5: /* FVMC_CELL_PYRAM */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*5],
                      (unsigned long)global_connect_s[i*5+1],
                      (unsigned long)global_connect_s[i*5+2],
                      (unsigned long)global_connect_s[i*5+3],
                      (unsigned long)global_connect_s[i*5+4]);
    break;

  case 6: /* FVMC_CELL_PRISM */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*6],
                      (unsigned long)global_connect_s[i*6+1],
                      (unsigned long)global_connect_s[i*6+2],
                      (unsigned long)global_connect_s[i*6+3],
                      (unsigned long)global_connect_s[i*6+4],
                      (unsigned long)global_connect_s[i*6+5]);
    break;

  case 8: /* FVMC_CELL_HEXA */
    for (i = 0, j = num_start ; j < num_end ; i++, j++)
      bftc_file_printf(f,
                      "%12lu : "
                      "%12lu %12lu %12lu %12lu\n"
                      "               "
                      "%12lu %12lu %12lu %12lu\n",
                      (unsigned long)j,
                      (unsigned long)global_connect_s[i*8],
                      (unsigned long)global_connect_s[i*8+1],
                      (unsigned long)global_connect_s[i*8+2],
                      (unsigned long)global_connect_s[i*8+3],
                      (unsigned long)global_connect_s[i*8+4],
                      (unsigned long)global_connect_s[i*8+5],
                      (unsigned long)global_connect_s[i*8+6],
                      (unsigned long)global_connect_s[i*8+7]);
    break;

  default:
    assert(0);
  }

}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write strided local connectivity to a text file
 *
 * parameters:
 *   stride  <-- number of vertices per element type
 *   n_elems <-- number of elements
 *   connect <-- connectivity array
 *   f       <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_connect_l(const int                stride,
                 const fvmc_lnum_t         n_elems,
                 const fvmc_lnum_t         connect[],
                 bftc_file_t        *const f)
{
  fvmc_lnum_t  i;

  switch(stride) {

  case 2: /* edge */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f, "%12lu : %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*2],
                      (unsigned long)connect[i*2+1]);
    break;

  case 3: /* FVMC_FACE_TRIA */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*3],
                      (unsigned long)connect[i*3+1],
                      (unsigned long)connect[i*3+2]);
    break;

  case 4: /* FVMC_FACE_QUAD or FVMC_CELL_TETRA */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*4],
                      (unsigned long)connect[i*4+1],
                      (unsigned long)connect[i*4+2],
                      (unsigned long)connect[i*4+3]);
    break;

  case 5: /* FVMC_CELL_PYRAM */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*5],
                      (unsigned long)connect[i*5+1],
                      (unsigned long)connect[i*5+2],
                      (unsigned long)connect[i*5+3],
                      (unsigned long)connect[i*5+4]);
    break;

  case 6: /* FVMC_CELL_PRISM */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f, "%12lu : %12lu %12lu %12lu %12lu %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*6],
                      (unsigned long)connect[i*6+1],
                      (unsigned long)connect[i*6+2],
                      (unsigned long)connect[i*6+3],
                      (unsigned long)connect[i*6+4],
                      (unsigned long)connect[i*6+5]);
    break;

  case 8: /* FVMC_CELL_HEXA */
    for (i = 0 ; i < n_elems ; i++)
      bftc_file_printf(f,
                      "%12lu : "
                      "%12lu %12lu %12lu %12lu\n"
                      "               "
                      "%12lu %12lu %12lu %12lu\n",
                      (unsigned long)(i+1),
                      (unsigned long)connect[i*8],
                      (unsigned long)connect[i*8+1],
                      (unsigned long)connect[i*8+2],
                      (unsigned long)connect[i*8+3],
                      (unsigned long)connect[i*8+4],
                      (unsigned long)connect[i*8+5],
                      (unsigned long)connect[i*8+6],
                      (unsigned long)connect[i*8+7]);
    break;

  default:
    assert(0);
  }

}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to a text file in parallel mode
 *
 * parameters:
 *   section                      <-- pointer to nodal mesh section structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   rank                         <-- rank in communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller       -- global connectivity slice provided
 *                                    by caller
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra_g(const fvmc_nodal_section_t  *const section,
                          const fvmc_io_num_t         *const global_vertex_num,
                          MPI_Comm                  comm,
                          const int                 rank,
                          const fvmc_gnum_t          global_s_size,
                          const fvmc_gnum_t          global_connect_s_size_caller,
                          fvmc_gnum_t                global_connect_s_caller[],
                          bftc_file_t         *const f)

{
  fvmc_lnum_t  i, j, k, l;
  fvmc_gnum_t  i_s, j_s, k_s;

  fvmc_lnum_t  cell_length, face_length;
  fvmc_lnum_t  face_id;
  fvmc_gnum_t  global_num_start = 1;
  fvmc_gnum_t  global_num_end = 0;
  fvmc_gnum_t  global_cell_face_idx_shift = 0 ;
  fvmc_gnum_t  global_cell_vtx_idx_shift = 0 ;

  fvmc_lnum_t  *face_lengths = NULL;
  fvmc_lnum_t  *cell_vtx_idx = NULL;
  fvmc_lnum_t  *cell_connect = NULL;
  fvmc_gnum_t  *global_cell_face_idx_s = NULL;
  fvmc_gnum_t  *global_face_lengths_s = NULL;
  fvmc_gnum_t  *global_cell_vtx_idx_s = NULL;

  fvmc_gnum_t   global_face_lengths_s_size = 0;
  fvmc_gnum_t   global_face_lengths_s_size_prev = 0;

  fvmc_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvmc_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvmc_gnum_t  *global_connect_s = global_connect_s_caller;

  fvmc_gather_slice_t   *polyhedra_slice = NULL;

  const fvmc_gnum_t   n_g_polyhedra = fvmc_io_num_get_global_count
                                       (section->global_element_num);

  /* Allocate memory for additionnal indexes */

  BFTC_MALLOC(global_cell_face_idx_s, global_s_size + 1, fvmc_gnum_t);
  BFTC_MALLOC(global_cell_vtx_idx_s, global_s_size + 1, fvmc_gnum_t);

  /* Every face should have at least 3 vertices, so cell->vertex connectivity
     should be at least 3 times the size of the cell->face connectivity;
     So we choose 1/3 (rounded up) of the size of the cell->vertex slice
     buffer as the initial size of the face_lengths slice buffer */

  global_face_lengths_s_size = (global_connect_s_size / 3) + 1;
  global_face_lengths_s_size_prev = global_face_lengths_s_size;
  BFTC_MALLOC(global_face_lengths_s, global_face_lengths_s_size, fvmc_gnum_t);

  /* Build local polyhedron indexes and connectivity information */

  BFTC_MALLOC(cell_vtx_idx,
             section->n_elements + 1,
             fvmc_lnum_t);

  BFTC_MALLOC(face_lengths,
             section->face_index[section->n_elements],
             fvmc_lnum_t);

  j_s = 0;
  l = 0;

  cell_vtx_idx[0] = 0;
  for (i = 0 ; i < section->n_elements ; i++) {
    cell_length = 0;
    for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
      face_id = FVMC_ABS(section->face_num[j]) - 1;
      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);
      face_lengths[l++] = face_length;
      cell_length += face_length;
    }
    cell_vtx_idx[i+1] = cell_vtx_idx[i] + cell_length;
  }

  BFTC_MALLOC(cell_connect,
             cell_vtx_idx[section->n_elements],
             fvmc_lnum_t);

  l = 0;

  for (i = 0 ; i < section->n_elements ; i++) {
    for (j = section->face_index[i] ; j < section->face_index[i+1] ; j++) {
      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        for (k = section->vertex_index[face_id] ;
             k < section->vertex_index[face_id+1] ;
             k++)
          cell_connect[l++] = section->vertex_num[k];
      }
      else {
        face_id = -section->face_num[j] - 1;
        k = section->vertex_index[face_id] ;
        cell_connect[l++] = section->vertex_num[k];
          for (k = section->vertex_index[face_id+1] - 1 ;
               k > section->vertex_index[face_id] ;
               k--)
            cell_connect[l++] = section->vertex_num[k];
      }
    }
  }

  /* Export by slices */

  polyhedra_slice = fvmc_gather_slice_create(section->global_element_num,
                                            global_s_size,
                                            comm);

  while (fvmc_gather_slice_advance(polyhedra_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather cell->vertices index */

    fvmc_gather_slice_index(cell_vtx_idx,
                           global_cell_vtx_idx_s,
                           section->global_element_num,
                           comm,
                           polyhedra_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvmc_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_cell_vtx_idx_s,
                                    polyhedra_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFTC_REALLOC(global_connect_s, global_connect_s_size, fvmc_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now gather cell->vertices connectivity */

    fvmc_gather_indexed_numbers(cell_vtx_idx,
                               cell_connect,
                               global_connect_s,
                               global_vertex_num,
                               section->global_element_num,
                               comm,
                               global_cell_vtx_idx_s,
                               polyhedra_slice);

    /* Now build the slice index for number of vertices per face */

    fvmc_gather_slice_index(section->face_index,
                           global_cell_face_idx_s,
                           section->global_element_num,
                           comm,
                           polyhedra_slice);

    /* If the face_lengths slice buffer is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (rank == 0) {
      global_face_lengths_s_size =
        FVMC_MAX(global_face_lengths_s_size,
                global_cell_face_idx_s[global_num_end - global_num_start]);
    }
    MPI_Bcast(&global_face_lengths_s_size, 1, FVMC_MPI_GNUM, 0, comm);
    if (global_face_lengths_s_size_prev < global_face_lengths_s_size) {
      BFTC_REALLOC(global_face_lengths_s,
                  global_face_lengths_s_size, fvmc_gnum_t);
      global_face_lengths_s_size_prev = global_face_lengths_s_size;
    }

    /* Now gather the number of vertices per face */

    fvmc_gather_indexed_numbers(section->face_index,
                               face_lengths,
                               global_face_lengths_s,
                               NULL,
                               section->global_element_num,
                               comm,
                               global_cell_face_idx_s,
                               polyhedra_slice);

    /* Do all printing for cells on rank 0 */

    if (rank == 0) {

      int  line_values;
      char str_num_cell[16];
      char str_num_face[16];
      char str_idx_cell[32];
      char str_idx_face[32];

      /* Print cell connectivity */

      k_s = 0;

      /* Loop on polyhedral cells in slice */

      for (i = 0, i_s = global_num_start ; i_s < global_num_end ; i++, i_s++) {

        /* Loop on cell faces */

        for (j = 0, j_s = global_cell_face_idx_s[i] ;
             j_s < global_cell_face_idx_s[i+1] ;
             j++, j_s++) {

          /* Print cell and face numbers and indexes */

          if (j_s == global_cell_face_idx_s[i]) {
            sprintf(str_num_cell, "%12lu", (unsigned long)i_s);
            sprintf(str_idx_cell, "[%lu] :",
                    (unsigned long)(j_s + global_cell_face_idx_shift + 1));
          }
          else {
            str_num_cell[0] = '\0';
            str_idx_cell[0] = '\0';
          }
          sprintf(str_num_face, "%5u", (unsigned)(j+1));
          sprintf(str_idx_face, "[%lu] :",
                    (unsigned long)(k_s + global_cell_vtx_idx_shift + 1));

          bftc_file_printf(f, "%12s %14s %5s %14s",
                          str_num_cell, str_idx_cell,
                          str_num_face, str_idx_face);

          /* Print face vertex numbers */

          line_values = 0;
          for (k = 0 ; k < (fvmc_lnum_t)global_face_lengths_s[j_s] ; k++) {
            if (line_values > 2) {
              line_values = 0;
              bftc_file_printf(f,"\n%48s", "");
            }
            bftc_file_printf(f, " %12lu",
                            (unsigned long)global_connect_s[k_s++]);
            line_values++;
          }
          bftc_file_printf(f, "\n");

        } /* End of loop on cell faces */

        assert(k_s == global_cell_vtx_idx_s[i+1]);

      } /* End of loop on polyhedral cells in slice */

      if (global_num_end > n_g_polyhedra) {
        str_num_cell[0] = '\0';
        sprintf(str_idx_cell, "[%lu] :",
                (unsigned long)(j_s + global_cell_face_idx_shift + 1));
        str_num_face[0] = '\0';
        sprintf(str_idx_face, "[%lu] :",
                (unsigned long)(k_s + global_cell_vtx_idx_shift + 1));
        bftc_file_printf(f, "%12s %14s %5s %14s",
                        str_num_cell, str_idx_cell,
                        str_num_face, str_idx_face);
      }

      /* Update shift for conversion from slice index to true index */

      global_cell_vtx_idx_shift
        += global_cell_vtx_idx_s[global_num_end - global_num_start];
      global_cell_face_idx_shift
        += global_cell_face_idx_s[global_num_end - global_num_start];

    } /* End of printing for rank 0 for this slice */

  }

  fvmc_gather_slice_destroy(polyhedra_slice);

  /* Free memory */

  BFTC_FREE(global_cell_face_idx_s);
  BFTC_FREE(global_cell_vtx_idx_s);
  BFTC_FREE(global_face_lengths_s);
  BFTC_FREE(cell_vtx_idx);
  BFTC_FREE(face_lengths);
  BFTC_FREE(cell_connect);

  if (global_connect_s != global_connect_s_caller)
    BFTC_FREE(global_connect_s);
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to a text file in serial mode
 *
 * parameters:
 *   section      <-- pointer to nodal mesh section structure
 *   f            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra_l(const fvmc_nodal_section_t  *const section,
                          bftc_file_t                 *const f)

{
  fvmc_lnum_t  i, j, k, l;

  fvmc_lnum_t  connect_length, face_length;
  fvmc_lnum_t  face_id;

  int  face_sgn, line_values;
  char str_num_cell[16];
  char str_num_face[16];
  char str_idx_cell[32];
  char str_idx_face[32];

  /* Print cell connectivity directly, without using extra buffers */

  connect_length = 0;
  j = 0;

  /* Loop on all polyhedral cells */

  for (i = 0 ; i < section->n_elements ; i++) {

    /* Loop on cell faces */

    for (j = section->face_index[i] ;
         j < section->face_index[i+1] ;
         j++) {

      /* Print cell and face numbers and indexes */

      if (j == section->face_index[i]) {
        sprintf(str_num_cell, "%12lu", (unsigned long)(i+1));
        sprintf(str_idx_cell, "[%lu] :",
                (unsigned long)(section->face_index[i] + 1));
      }
      else {
        str_num_cell[0] = '\0';
        str_idx_cell[0] = '\0';
      }
      sprintf(str_num_face, "%5u", (unsigned)(j-section->face_index[i]+1));
      sprintf(str_idx_face, "[%lu] :", (unsigned long)(connect_length+1));

      bftc_file_printf(f, "%12s %14s %5s %14s",
                      str_num_cell, str_idx_cell,
                      str_num_face, str_idx_face);

      /* Print face vertex numbers */

      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        face_sgn = 1;
      }
      else {
        face_id = -section->face_num[j] - 1;
        face_sgn = -1;
      }

      line_values = 0;
      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);
      connect_length += face_length;

      for (k = 0 ; k < face_length ; k++) {
        l =   section->vertex_index[face_id]
            + (face_length + (k*face_sgn))%face_length;
        if (line_values > 2) {
          line_values = 0;
          bftc_file_printf(f,"\n%48s", "");
        }
        bftc_file_printf(f, " %12lu",
                        (unsigned long)section->vertex_num[l]);
        line_values++;
      }
      bftc_file_printf(f, "\n");

    } /* End of loop on cell faces */

  } /* End of loop on polyhedral cells */

  str_num_cell[0] = '\0';
  sprintf(str_idx_cell, "[%lu] :", (unsigned long)(j + 1));
  str_num_face[0] = '\0';
  sprintf(str_idx_face, "[%lu] :", (unsigned long)(connect_length+1));
  bftc_file_printf(f, "%12s %14s %5s %14s",
                  str_num_cell, str_idx_cell,
                  str_num_face, str_idx_face);
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text file in parallel mode
 *
 * parameters:
 *   section                      <-- pointer to nodal mesh section structure
 *   global_vertex_num            <-- pointer to vertex global numbering
 *   comm                         <-- associated MPI communicator
 *   rank                         <-- rank in communicator
 *   n_ranks                      <-- number of processes in communicator
 *   global_s_size                <-- global slice size
 *   global_connect_s_size_caller <-- global connectivity slice size
 *                                    defined by caller
 *   global_connect_s_caller       -- global connectivity slice provided
 *                                    by caller
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polygons_g(const fvmc_nodal_section_t  *const section,
                         const fvmc_io_num_t         *const global_vertex_num,
                         MPI_Comm                  comm,
                         const int                 rank,
                         const fvmc_gnum_t          global_s_size,
                         const fvmc_gnum_t          global_connect_s_size_caller,
                         fvmc_gnum_t                global_connect_s_caller[],
                         bftc_file_t         *const f)

{
  fvmc_lnum_t  i, j;
  fvmc_gnum_t  i_s, j_s;

  fvmc_gnum_t  global_num_start;
  fvmc_gnum_t  global_num_end;
  fvmc_gnum_t  global_idx_shift = 0 ;

  fvmc_gnum_t  *global_idx_s = NULL;

  fvmc_gnum_t   global_connect_s_size = global_connect_s_size_caller;
  fvmc_gnum_t   global_connect_s_size_prev = global_connect_s_size_caller;
  fvmc_gnum_t  *global_connect_s = global_connect_s_caller;

  fvmc_gather_slice_t   *polygons_slice = NULL;

  const fvmc_gnum_t   n_g_polygons = fvmc_io_num_get_global_count
                                       (section->global_element_num);

  /* Allocate memory for additionnal indexes */

  BFTC_MALLOC(global_idx_s, global_s_size + 1, fvmc_gnum_t);

  /* Export by slices */

  polygons_slice = fvmc_gather_slice_create(section->global_element_num,
                                           global_s_size,
                                           comm);

  while (fvmc_gather_slice_advance(polygons_slice,
                                  &global_num_start,
                                  &global_num_end) == 0) {

    /* Gather face->vertices index */

    fvmc_gather_slice_index(section->vertex_index,
                           global_idx_s,
                           section->global_element_num,
                           comm,
                           polygons_slice);

    /* Recompute maximum value of global_num_end for this slice */

    fvmc_gather_resize_indexed_slice(10,
                                    &global_num_end,
                                    &global_connect_s_size,
                                    comm,
                                    global_idx_s,
                                    polygons_slice);

    /* If the buffer passed to this function is too small, allocate a
       larger one; in this case, we may as well keep it for all slices */

    if (global_connect_s_size_prev < global_connect_s_size) {
      if (global_connect_s == global_connect_s_caller)
        global_connect_s = NULL;
      BFTC_REALLOC(global_connect_s, global_connect_s_size, fvmc_gnum_t);
      global_connect_s_size_prev = global_connect_s_size;
    }

    /* Now gather face->vertices connectivity */

    fvmc_gather_indexed_numbers(section->vertex_index,
                               section->vertex_num,
                               global_connect_s,
                               global_vertex_num,
                               section->global_element_num,
                               comm,
                               global_idx_s,
                               polygons_slice);

    /* Do all printing for faces on rank 0 */

    if (rank == 0) {

      int  line_values;
      char str_num_face[16];
      char str_idx_face[32];

      /* Print face connectivity */

      /* Loop on polygonal faces in slice */

      for (i = 0, i_s = global_num_start ; i_s < global_num_end ; i++, i_s++) {

        /* Print cell and face numbers and indexes */

        sprintf(str_num_face, "%12lu", (unsigned long)i_s);
        sprintf(str_idx_face, "[%lu] :",
                (unsigned long)(global_idx_s[i] + global_idx_shift + 1));

        bftc_file_printf(f, "%12s %14s",
                        str_num_face, str_idx_face);

        /* Print face vertex numbers */

        line_values = 0;
        for (j = 0, j_s = global_idx_s[i] ;
             j_s < global_idx_s[i+1] ;
             j++, j_s++) {
          if (line_values > 2) {
            line_values = 0;
            bftc_file_printf(f,"\n%27s", "");
          }
          bftc_file_printf(f, " %12lu",
                          (unsigned long)global_connect_s[j_s]);
          line_values++;
        }
        bftc_file_printf(f, "\n");

      } /* End of loop on polyhedral cells in slice */

      if (global_num_end > n_g_polygons) {
        str_num_face[0] = '\0';
        sprintf(str_idx_face, "[%lu] :",
                (unsigned long)(global_idx_s[i] + global_idx_shift + 1));
        bftc_file_printf(f, "%12s %14s",
                        str_num_face, str_idx_face);
      }

      /* Update shift for conversion from slice index to true index */

      global_idx_shift
        += global_idx_s[global_num_end - global_num_start + 1];

    } /* End of printing for rank 0 for this slice */

  }

  fvmc_gather_slice_destroy(polygons_slice);

  /* Free memory */

  BFTC_FREE(global_idx_s);

  if (global_connect_s != global_connect_s_caller)
    BFTC_FREE(global_connect_s);
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to a text file in serial mode
 *
 * parameters:
 *   section      <-- pointer to nodal mesh section structure
 *   f                            <-- pointer to associated file
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polygons_l(const fvmc_nodal_section_t  *const section,
                         bftc_file_t                 *const f)

{
  fvmc_lnum_t  i, j;

  int  line_values;
  char str_num_face[16];
  char str_idx_face[32];

  /* Print face connectivity directly, without using extra buffers */

  j = 0; /* Initialize here in case section->n_elements = 0 */

  /* Loop on all polygonal faces */

  for (i = 0 ; i < section->n_elements ; i++) {

    /* Print face numbers and indexes */

    sprintf(str_num_face, "%12lu", (unsigned long)(i+1));
    sprintf(str_idx_face, "[%lu] :",
            (unsigned long)(section->vertex_index[i] + 1));

    bftc_file_printf(f, "%12s %14s",
                    str_num_face, str_idx_face);

    /* Print face vertex numbers */

    line_values = 0;

    for (j = section->vertex_index[i] ;
         j < section->vertex_index[i+1] ;
         j++) {
      if (line_values > 2) {
        line_values = 0;
        bftc_file_printf(f,"\n%27s", "");
      }
      bftc_file_printf(f, " %12lu",
                      (unsigned long)section->vertex_num[j]);
      line_values++;
    }
    bftc_file_printf(f, "\n");

  } /* End of loop on polygonal faces */

  str_num_face[0] = '\0';
  sprintf(str_idx_face, "[%lu] :", (unsigned long)(j + 1));
  bftc_file_printf(f, "%12s %14s",
                  str_num_face, str_idx_face);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to text file writer.
 *
 * parameters:
 *   name    <-- base output case name.
 *   options <-- whitespace separated, lowercase options list
 *   comm    <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque text file writer structure.
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)
void *
fvmc_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvmc_writer_time_dep_t         time_dependency,
                        const MPI_Comm                      comm)
#else
void *
fvmc_to_text_init_writer(const char                   *const name,
                        const char                   *const path,
                        const char                   *const options,
                        const fvmc_writer_time_dep_t         time_dependency)
#endif

{
  fvmc_to_text_writer_t  *this_writer = NULL;
  int  rank = 0;

  /* Initialize writer */

  BFTC_MALLOC(this_writer, 1, fvmc_to_text_writer_t);

  this_writer->time_dependency = time_dependency;

  this_writer->rank = 0;
  this_writer->n_ranks = 1;

#if defined(FVMC_HAVE_MPI)
  {
    int mpi_flag, n_ranks;
    MPI_Initialized(&mpi_flag);

    if (mpi_flag) {
      this_writer->comm = comm;
      MPI_Comm_rank(this_writer->comm, &rank);
      MPI_Comm_size(this_writer->comm, &n_ranks);
      this_writer->rank = rank;
      this_writer->n_ranks = n_ranks;
    }
    else
      this_writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(FVMC_HAVE_MPI) */

  if (rank == 0) {

    char * file_name;
    int path_len = 0;
    const char extension[] = ".txt";

    if (path != NULL)
      path_len = strlen(path);

    BFTC_MALLOC(file_name,
               path_len + strlen(name) + strlen(extension) + 1,
               char);

    if (path != NULL)
      strcpy(file_name, path);
    else
      file_name[0] = '\0';

    strcat(file_name, name);
    strcat(file_name, extension);
    this_writer->file = bftc_file_open(file_name,
                                      BFTC_FILE_MODE_WRITE,
                                      BFTC_FILE_TYPE_TEXT);
    BFTC_FREE(file_name);

  }
  else

    this_writer->file = NULL;

  /* Parse options */

  if (rank == 0 && options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1 ; i2 < l_tot && options[i2] != ' ' ; i2++);
      l_opt = i2 - i1 + 1;

      bftc_file_printf(this_writer->file,
                      _("Option: %*s\n"), l_opt, options + i1);

      i1 = i2;

    }

  }

  /* Return writer */

  return this_writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to text file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque text file writer structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvmc_to_text_finalize_writer(void  *this_writer_p)
{
  fvmc_to_text_writer_t *this_writer = (fvmc_to_text_writer_t *)this_writer_p;

  if (this_writer->file != NULL)
    this_writer->file = bftc_file_free(this_writer->file);

  BFTC_FREE(this_writer);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a text file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvmc_to_text_export_nodal(void               *const this_writer_p,
                         const fvmc_nodal_t  *const mesh)
{
  int         section_id;
  fvmc_lnum_t  i, j;

  fvmc_to_text_writer_t  *this_writer = (fvmc_to_text_writer_t *)this_writer_p;
  bftc_file_t  *f = NULL;

  fvmc_gnum_t   global_connect_s_size, global_s_size;
  fvmc_gnum_t   n_g_vertices = 0;
  fvmc_gnum_t   n_g_edges = 0;
  fvmc_gnum_t   n_g_faces = 0;
  fvmc_gnum_t   n_g_cells = 0;
  fvmc_gnum_t  * n_g_elements_section = NULL;

#if defined(FVMC_HAVE_MPI)

  fvmc_gnum_t  *global_connect_s = NULL;
  MPI_Comm    comm = this_writer->comm;

#endif

  const int  rank = this_writer->rank;
  const int  n_ranks = this_writer->n_ranks;

  /* Initialization */
  /*----------------*/

  f = this_writer->file;

  /* Buffer sizes required in parallel mode, global sizes always required */
  /*----------------------------------------------------------------------*/

  BFTC_MALLOC(n_g_elements_section, mesh->n_sections, fvmc_gnum_t);

  fvmc_writer_def_nodal_buf_size(mesh,
                                n_ranks,
                                12,
                                5,
                                &n_g_vertices,
                                n_g_elements_section,
                                &global_s_size,
                                &global_connect_s_size);

  for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

    const fvmc_nodal_section_t  *const  section = mesh->sections[section_id];

    switch(section->entity_dim) {
    case 1:
      n_g_edges += n_g_elements_section[section_id];
      break;
    case 2:
      n_g_faces += n_g_elements_section[section_id];
      break;
    case 3:
      n_g_cells += n_g_elements_section[section_id];
      break;
    default:
      assert(0);
    }

  }

  /* Global indicators */
  /*--------------------*/

  if (rank == 0) {

    if (mesh->name != NULL)
      bftc_file_printf(f, _("\n"
                           "Mesh name: %s\n"),
                      mesh->name);
    else
      bftc_file_printf(f, _("\n"
                           "Unnamed mesh\n"));

    bftc_file_printf(f, _("\n"
                         "Mesh dimension:     %d\n"
                         "Number of domains:  %d\n"
                         "Number of sections:  %d\n"),
                    mesh->dim, mesh->n_doms, mesh->n_sections);

    bftc_file_printf(f, _("\n"
                         "Number of cells:               %d\n"
                         "Number of faces:               %d\n"
                         "Number of edges:               %d\n"
                         "Number of vertices:            %d\n"),
                    n_g_cells,
                    n_g_faces,
                    n_g_edges,
                    n_g_vertices);

  }

  /* Vertex coordinates */
  /*--------------------*/

  {
    const int      stride = mesh->dim;
    const double  *local_coords;
    double        *coords_tmp = NULL;

    if (mesh->parent_vertex_num != NULL) {
      BFTC_MALLOC(coords_tmp, stride * mesh->n_vertices, double);
      for (i = 0 ; i < mesh->n_vertices ; i++) {
        for (j = 0 ; j < stride ; j++)
          coords_tmp[i*stride + j]
            = mesh->vertex_coords[(mesh->parent_vertex_num[i]-1)*stride + j];
      }
      local_coords = coords_tmp;
    }
    else
      local_coords = mesh->vertex_coords;

    if (rank == 0)
      bftc_file_printf(f, _("\nVertex coordinates:\n\n"));

    /* loop on slices in parallel mode, use whole array in serial mode */

#if defined(FVMC_HAVE_MPI)

    if (n_ranks > 1) { /* start of output in parallel mode */

      fvmc_gnum_t   global_num_start = 1;
      fvmc_gnum_t   global_num_end = 0;

      fvmc_gather_slice_t   *vertices_slice = NULL;
      double               *global_coords_s = NULL;

      BFTC_MALLOC(global_coords_s, global_s_size * stride, double);

      vertices_slice = fvmc_gather_slice_create(mesh->global_vertex_num,
                                               global_s_size,
                                               comm);

      while (fvmc_gather_slice_advance(vertices_slice,
                                      &global_num_start,
                                      &global_num_end) == 0) {

        fvmc_gather_array(local_coords,
                         global_coords_s,
                         MPI_DOUBLE,
                         (size_t)stride,
                         mesh->global_vertex_num,
                         comm,
                         vertices_slice);

        if (rank == 0)
          _write_slice_vector(stride,
                              global_num_start,
                              global_num_end,
                              global_coords_s,
                              f);

      }

      fvmc_gather_slice_destroy(vertices_slice);

      BFTC_FREE(global_coords_s);

    } /* end of output in parallel mode */

#endif /* defined(FVMC_HAVE_MPI) */

    if (n_ranks == 1) { /* start of output in serial mode */

      _write_slice_vector(stride,
                          1,
                          (fvmc_gnum_t)(mesh->n_vertices + 1),
                          local_coords,
                          f);

    } /* end of output in serial mode */

    if (coords_tmp != NULL)
      BFTC_FREE(coords_tmp);
  }

  /* Allocate connectivity buffer for use wih all types of elements */

#if defined(FVMC_HAVE_MPI)

  if (n_ranks > 1)
    BFTC_MALLOC(global_connect_s, global_connect_s_size, fvmc_gnum_t);

#endif

  /* Section connectivity */
  /*----------------------*/

  for (section_id = 0 ; section_id < mesh->n_sections ; section_id++) {

    const fvmc_nodal_section_t  *const  section = mesh->sections[section_id];

    if (rank == 0)
      bftc_file_printf(f, _("\nSection: %s\n"
                           "  Number of elements: %lu\n\n"),
                      _(fvmc_elements_type_name[section->type]),
                      (unsigned long)(n_g_elements_section[section_id]));

    /* Output for strided (regular) element types */
    /*--------------------------------------------*/

    if (section->stride > 0) {

#if defined(FVMC_HAVE_MPI)

      if (n_ranks > 1) { /* start of output in parallel mode */

        fvmc_gnum_t   global_num_start;
        fvmc_gnum_t   global_num_end;

        fvmc_gather_slice_t   *elements_slice = NULL;

        elements_slice
          = fvmc_gather_slice_create(section->global_element_num,
                                    global_s_size,
                                    comm);

        while (fvmc_gather_slice_advance(elements_slice,
                                        &global_num_start,
                                        &global_num_end) == 0) {

          fvmc_gather_strided_connect(section->vertex_num,
                                     global_connect_s,
                                     section->stride,
                                     mesh->global_vertex_num,
                                     section->global_element_num,
                                     comm,
                                     elements_slice);

          if (rank == 0)
            _write_slice_connect_g(section->stride,
                                   global_num_start,
                                   global_num_end,
                                   global_connect_s,
                                   f);

        }

        fvmc_gather_slice_destroy(elements_slice);

      } /* end of output in parallel mode */

#endif /* defined(FVMC_HAVE_MPI) */

      if (n_ranks == 1) { /* start of output in serial mode */

        _write_connect_l(section->stride,
                         section->n_elements,
                         section->vertex_num,
                         f);

      } /* end of output in serial mode */

    } /* end of output for strided element types */

    /* Output for polygons */
    /*---------------------*/

    else if (section->type == FVMC_FACE_POLY) {

#if defined(FVMC_HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        _export_nodal_polygons_g(section,
                                 mesh->global_vertex_num,
                                 comm,
                                 rank,
                                 global_s_size,
                                 global_connect_s_size,
                                 global_connect_s,
                                 f);

      }

#endif /* defined(FVMC_HAVE_MPI) */

      if (n_ranks == 1)

        _export_nodal_polygons_l(section, f);

    }

    /* Output for polyhedra */
    /*----------------------*/

    else if (section->type == FVMC_CELL_POLY) {

#if defined(FVMC_HAVE_MPI)

      /* output in parallel mode */

      if (n_ranks > 1) {

        _export_nodal_polyhedra_g(section,
                                  mesh->global_vertex_num,
                                  comm,
                                  rank,
                                  global_s_size,
                                  global_connect_s_size,
                                  global_connect_s,
                                  f);

      }

#endif /* defined(FVMC_HAVE_MPI) */

      if (n_ranks == 1)

        _export_nodal_polyhedra_l(section, f);

    }

  } /* End of loop on sections */

  /* Free buffers */
  /*--------------*/

  BFTC_FREE(n_g_elements_section);

  /* Free buffers required in parallel mode */

#if defined(FVMC_HAVE_MPI)

  if (n_ranks > 1)
    BFTC_FREE(global_connect_s);

#endif /* defined(FVMC_HAVE_MPI) */

  /* Close dump file */
  /*-----------------*/

  if (rank == 0)
    bftc_file_flush(f);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
