/*============================================================================
 * Initialization of a nodal connectivity definition based upon
 * a (possibly partial) descending connectivity
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

#include <bftc_error.h>
#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 * Local headers associated with the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal_from_desc.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Enumeration definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Nodal connectivity reconstruction return code.
 *----------------------------------------------------------------------------*/

typedef enum {

  FVMC_NODAL_FROM_DESC_SUCCESS,  /* Successful reconstruction */
  FVMC_NODAL_FROM_DESC_FAILURE,  /* Reconstruction failure (connectivity pb.) */
  FVMC_NODAL_FROM_DESC_WORIENT   /* Reconstruction with orientation warning
                                   (for possible orientation problem) */
} fvmc_nodal_from_desc_t;

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Global static variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine if a given cell is a prism or a polyhedron
 *
 * This function should only be called on cells having 2 triangular and
 * 3 quadrangular faces; if the triangles are adjacent, the cell
 * is not a prism, and is considered to be a polyhedron. Otherwise, it
 * seems to be truly a prism (supposing it is closed and well oriented).
 *
 * parameters:
 *   cell_id,        <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

inline static fvmc_element_t
_is_prism_or_poly(const fvmc_lnum_t   cell_id,
                  const int          n_face_lists,
                  const fvmc_lnum_t   face_list_shift[],
                  const fvmc_lnum_t  *face_vertex_idx[],
                  const fvmc_lnum_t  *face_vertex_num[],
                  const fvmc_lnum_t   cell_face_idx[],
                  const fvmc_lnum_t   cell_face_num[])
{
  int         vtx_id_1, vtx_id_2;
  fvmc_lnum_t  face_id, fl;
  fvmc_lnum_t  idx, idx_start, idx_end;
  fvmc_lnum_t  n_face_vertices;
  fvmc_lnum_t  vtx, vertex_id_start, vertex_id_end;

  fvmc_lnum_t  vtx_tria[6];
  int         n_trias = 0;

  /* Extract 2 triangles */

  idx_start = cell_face_idx[cell_id]     - 1;
  idx_end   = cell_face_idx[cell_id + 1] - 1;

  for (idx = idx_start ; idx < idx_end ; idx++) {

    face_id = FVMC_ABS(cell_face_num[idx]) - 1;

    for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
    assert(fl > -1);
    face_id -= face_list_shift[fl];

    vertex_id_start = face_vertex_idx[fl][face_id] - 1;
    vertex_id_end   = face_vertex_idx[fl][face_id + 1] - 1;
    n_face_vertices = vertex_id_end - vertex_id_start;

    if (n_face_vertices == 3) {

      if (cell_face_num[idx] > 0) {
        for (vtx = 0 ; vtx < 3 ; vtx++)
          vtx_tria[n_trias*3 + vtx]
            = face_vertex_num[fl][vertex_id_start + vtx];
      }
      else {
        for (vtx = 0 ; vtx < 3 ; vtx++)
          vtx_tria[n_trias*3 + vtx]
            = face_vertex_num[fl][vertex_id_end - 1 - vtx];
      }

      n_trias += 1;

      if (n_trias == 2) /* Once we have found the two triangles, */
        break;          /*  we do not need to look at other faces */
    }

  }

  assert(n_trias == 2);

  /* Are the triangles adjacent ? */

  for (vtx_id_1 = 0; vtx_id_1 < 3; vtx_id_1++) {
    for (vtx_id_2 = 0; vtx_id_2 < 3; vtx_id_2++) {
      if (vtx_tria[3 + vtx_id_2] == vtx_tria[vtx_id_1]) {
        return FVMC_CELL_POLY;
      }
    }
  }

  /* If no adjacent triangles were found, we have a prism */

  return FVMC_CELL_PRISM;
}

/*----------------------------------------------------------------------------
 * Determination of a given cell's type.
 *
 * If the optional cell_vtx_tria[3*4] and cell_vtx_quad[4*6] arrays are given,
 * they are filled with the vertex indexes of the cell's triangle and
 * quadrangle type faces (as long as no simple polyhedral face is encountered
 * and we have no more than 4 triangles or 6 quadrangles), for easier
 * nodal cell reconstuction. Vertex indexes are ordered as usual for
 * an outward pointing normal.
 *
 * parameters:
 *   cell_id,        <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_vtx_tria   --> local triangle definitions (optional, 4 max.)
 *   cell_vtx_quad   --> local quadrangle definitions (optional, 6 max.)
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

inline static fvmc_element_t
_nodal_cell_from_desc(const fvmc_lnum_t   cell_id,
                      const int          n_face_lists,
                      const fvmc_lnum_t   face_list_shift[],
                      const fvmc_lnum_t  *face_vertex_idx[],
                      const fvmc_lnum_t  *face_vertex_num[],
                      const fvmc_lnum_t   cell_face_idx[],
                      const fvmc_lnum_t   cell_face_num[],
                      fvmc_lnum_t        *cell_vtx_tria,
                      fvmc_lnum_t        *cell_vtx_quad)
{
  fvmc_lnum_t  n_cell_faces;
  fvmc_lnum_t  face_id, fl;
  fvmc_lnum_t  idx, idx_start, idx_end;
  fvmc_lnum_t  n_face_vertices;
  fvmc_lnum_t  vtx, vertex_id_start, vertex_id_end;

  fvmc_element_t cell_type;

  fvmc_lnum_t  n_trias = 0;
  fvmc_lnum_t  n_quads = 0;
  fvmc_lnum_t  n_ngons = 0;


  /* Guess connectivity types */
  /*--------------------------*/

  n_cell_faces = cell_face_idx[cell_id + 1] - cell_face_idx[cell_id];

  /* If we have more than 6 faces, we have a general simple polyhedron */

  if (n_cell_faces > 6)
    return FVMC_CELL_POLY;

  /*
    Otherwise, we probably have a "classical" element; for example,
    with 6 faces, we probably have a hexahedron, though we could
    also have a tetrahedron with one of its edges and thus two of
    its faces split. We count vertices per face to check.
  */

  idx_start = cell_face_idx[cell_id]     - 1;
  idx_end   = cell_face_idx[cell_id + 1] - 1;

  for (idx = idx_start ; idx < idx_end ; idx++) {

    face_id = FVMC_ABS(cell_face_num[idx]) - 1;

    for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
    assert(fl > -1);
    face_id -= face_list_shift[fl];

    vertex_id_start = face_vertex_idx[fl][face_id] - 1;
    vertex_id_end   = face_vertex_idx[fl][face_id + 1] - 1;
    n_face_vertices = vertex_id_end - vertex_id_start;

    if (n_face_vertices == 3) {

      if (cell_vtx_tria != NULL && n_trias < 4) {
        if (cell_face_num[idx] > 0) {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_tria[n_trias*3 + vtx]
              = face_vertex_num[fl][vertex_id_start + vtx];
        }
        else {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_tria[n_trias*3 + vtx]
              = face_vertex_num[fl][vertex_id_end - 1 - vtx];
        }
      }

      n_trias += 1;

    }
    else if (n_face_vertices == 4) {

      if (cell_vtx_quad != NULL && n_quads < 6) {
        if (cell_face_num[idx] > 0) {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_quad[n_quads*4 + vtx]
              = face_vertex_num[fl][vertex_id_start + vtx];
        }
        else {
          for (vtx = 0 ; vtx < n_face_vertices ; vtx++)
            cell_vtx_quad[n_quads*4 + vtx]
              = face_vertex_num[fl][vertex_id_end - 1 - vtx];
        }
      }

      n_quads += 1;
    }
    else

      n_ngons += 1;

  }

  /* Return element type */

  if (n_ngons > 0)
    cell_type = FVMC_CELL_POLY;
  else {
    if (n_trias == 0 && n_quads == 6)
      cell_type = FVMC_CELL_HEXA;
    else if (n_trias == 2 && n_quads == 3)
      cell_type = _is_prism_or_poly(cell_id,
                                    n_face_lists,
                                    face_list_shift,
                                    face_vertex_idx,
                                    face_vertex_num,
                                    cell_face_idx,
                                    cell_face_num);
    else if (n_trias == 4) {
      if (n_quads == 0)
        cell_type = FVMC_CELL_TETRA;
      else if (n_quads == 1)
        cell_type = FVMC_CELL_PYRAM;
      else
        cell_type = FVMC_CELL_POLY;
    }
    else
      cell_type = FVMC_CELL_POLY;
  }

  return cell_type;

}

/*----------------------------------------------------------------------------
 * Tetrahedron construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_tetra --> tetrahedron connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvmc_nodal_from_desc_t
_nodal_from_desc_cnv_cel_tetra(const fvmc_lnum_t  cell_vtx_tria[],
                               fvmc_lnum_t  cell_vtx_tetra[])
{
  fvmc_lnum_t  vertex_id, face_id;
  fvmc_lnum_t  direction;
  fvmc_lnum_t  vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("face 1 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bftc_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bftc_printf("face 3 : %d %d %d\n",
             cell_vtx_tria[6], cell_vtx_tria[7], cell_vtx_tria[8]);
  bftc_printf("face 4 : %d %d %d\n",
             cell_vtx_tria[9], cell_vtx_tria[10], cell_vtx_tria[11]);
#endif

  /*
    Base : vertices of the first face; we take the vertices in opposite
    order ("bottom" face numbering of the tetrahedron with outward
    pointing normal).

    *        x 4
    *       /|\
    *      / | \
    *     /  |  \
    *  1 x- -|- -x 3
    *     \  |  /
    *      \ | /
    *       \|/
    *        x 2
  */


  cell_vtx_tetra[0] = cell_vtx_tria[2];
  cell_vtx_tetra[1] = cell_vtx_tria[1];
  cell_vtx_tetra[2] = cell_vtx_tria[0];

  /*
    We have found 3 of 4 vertices; all other triangles should share
    vertex 4, and one of those should share vertices 1 and 2 of the
    base triangle.
  */

  vtx_num_1 = cell_vtx_tetra[0];
  vtx_num_2 = cell_vtx_tetra[1];

  direction = 0;

  for (face_id = 1 ; face_id < 4 ; face_id++) {

    for (vertex_id = 0 ; vertex_id < 3 ; vertex_id++) {

      vtx_num = cell_vtx_tria[face_id*3 + vertex_id];

      if (vtx_num == vtx_num_1) {
        if (cell_vtx_tria[face_id*3 + ((vertex_id+1) % 3)] == vtx_num_2) {
          direction = 1;
          break;
        }
        else if (cell_vtx_tria[face_id*3 + ((vertex_id-1+3) % 3)] == vtx_num_2) {
          direction = -1;
          break;
        }
      }

    }

    if (direction != 0)
      break;

  }

  if (direction == -1)
    warn_orient = true;
  else if (direction == 0)
    return FVMC_NODAL_FROM_DESC_FAILURE;

  cell_vtx_tetra[3]
    = cell_vtx_tria[face_id*3 + ((vertex_id + (3 + (2 * direction))) % 3)];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("tetra : %d %d %d %d\n",
             cell_vtx_tetra[0], cell_vtx_tetra[1],
             cell_vtx_tetra[2], cell_vtx_tetra[3]);
#endif

  if (warn_orient == true)
    return FVMC_NODAL_FROM_DESC_WORIENT;
  else
    return FVMC_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Pyramid construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_quad  <-- quadrangle faces connectivity
 *   cell_vtx_pyram --> pyramid connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvmc_nodal_from_desc_t
_nodal_from_desc_cnv_cel_pyram(const fvmc_lnum_t  cell_vtx_tria[],
                               const fvmc_lnum_t  cell_vtx_quad[],
                               fvmc_lnum_t  cell_vtx_pyram[])
{
  fvmc_lnum_t  vertex_id, face_id;
  fvmc_lnum_t  direction;
  fvmc_lnum_t  vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("face 1 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bftc_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bftc_printf("face 3 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bftc_printf("face 4 : %d %d %d\n",
             cell_vtx_tria[6], cell_vtx_tria[7], cell_vtx_tria[8]);
  bftc_printf("face 5 : %d %d %d\n",
             cell_vtx_tria[9], cell_vtx_tria[10], cell_vtx_tria[11]);
#endif

  /*
    Base : vertices of the quadrangle; we take the vertices in opposite
    order ("bottom" face numbering of the pyramid with outward
    pointing normal).

    *         5 x
    *          /|\
    *         //| \
    *        // |  \
    *     4 x/--|---x 3
    *      //   |  /
    *     //    | /
    *  1 x-------x 2
  */


  cell_vtx_pyram[0] = cell_vtx_quad[3];
  cell_vtx_pyram[1] = cell_vtx_quad[2];
  cell_vtx_pyram[2] = cell_vtx_quad[1];
  cell_vtx_pyram[3] = cell_vtx_quad[0];

  /*
    We have found 4 out of 5 vertices; all 4 triangles should share
    vertex 5, and one of those should share vertices 1 and 2 of the
    base quadrangle.
  */

  vtx_num_1 = cell_vtx_pyram[0];
  vtx_num_2 = cell_vtx_pyram[1];

  direction = 0;

  for (face_id = 0 ; face_id < 4 ; face_id++) {

    for (vertex_id = 0 ; vertex_id < 3 ; vertex_id++) {

      vtx_num = cell_vtx_tria[face_id*3 + vertex_id];

      if (vtx_num == vtx_num_1) {
        if (cell_vtx_tria[face_id*3 + ((vertex_id+1) % 3)] == vtx_num_2) {
          direction = 1;
          break;
        }
        else if (cell_vtx_tria[face_id*3 + ((vertex_id-1+3) % 3)] == vtx_num_2) {
          direction = -1;
          break;
        }
      }

    }

    if (direction != 0)
      break;

  }

  if (direction == -1)
    warn_orient = true;
  else if (direction == 0)
    return FVMC_NODAL_FROM_DESC_FAILURE;

  cell_vtx_pyram[4]
    = cell_vtx_tria[face_id*3 + ((vertex_id + (3 + (2 * direction))) % 3)];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("pyram : %d %d %d %d %d\n",
             cell_vtx_pyram[0], cell_vtx_pyram[1], cell_vtx_pyram[2],
             cell_vtx_pyram[3], cell_vtx_pyram[4]);
#endif

  if (warn_orient == true)
    return FVMC_NODAL_FROM_DESC_WORIENT;
  else
    return FVMC_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Prism (pentahedron) construction
 *
 * parameters:
 *   cell_vtx_tria  <-- triangular faces connectivity
 *   cell_vtx_quad  <-- quadrangle faces connectivity
 *   cell_vtx_prism --> prism connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvmc_nodal_from_desc_t
_nodal_from_desc_cnv_cel_prism(const fvmc_lnum_t  cell_vtx_tria[],
                               const fvmc_lnum_t  cell_vtx_quad[],
                               fvmc_lnum_t  cell_vtx_prism[])
{
  fvmc_lnum_t  vertex_id, face_id;
  fvmc_lnum_t  ipass, direction;
  fvmc_lnum_t  vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("face 1 : %d %d %d\n",
             cell_vtx_tria[0], cell_vtx_tria[1], cell_vtx_tria[2]);
  bftc_printf("face 2 : %d %d %d\n",
             cell_vtx_tria[3], cell_vtx_tria[4], cell_vtx_tria[5]);
  bftc_printf("face 3 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bftc_printf("face 4 : %d %d %d %d\n",
             cell_vtx_quad[4], cell_vtx_quad[5],
             cell_vtx_quad[6], cell_vtx_quad[7]);
  bftc_printf("face 5 : %d %d %d %d\n",
             cell_vtx_quad[8], cell_vtx_quad[9],
             cell_vtx_quad[10], cell_vtx_quad[11]);
#endif

  /*
    Base : vertices of the first triangle; we take the vertices in opposite
    order ("bottom" face numbering of the prism with outward
    pointing normal).

    *  4 x-------x 6
    *    |\     /|
    *    | \   / |
    *  1 x- \-/ -x 3
    *     \ 5x  /
    *      \ | /
    *       \|/
    *        x 2
  */


  cell_vtx_prism[0] = cell_vtx_tria[2];
  cell_vtx_prism[1] = cell_vtx_tria[1];
  cell_vtx_prism[2] = cell_vtx_tria[0];

  /*
    We have found 3 out of 6 vertices; we first seek the quadrangle sharing
    vertices 1 and 2, so as to determine vertices 4 and 5,
    then we seek the quadrangle sharing vertices 2 and 3, so as to
    determine vertices 5 and 6.
  */

  for (ipass = 0 ; ipass < 2 ; ipass++) {

    vtx_num_1 = cell_vtx_prism[    ipass];
    vtx_num_2 = cell_vtx_prism[1 + ipass];

    direction = 0;

    for (face_id = 0 ; face_id < 4 ; face_id++) {

      for (vertex_id = 0 ; vertex_id < 4 ; vertex_id++) {

        vtx_num = cell_vtx_quad[face_id*4 + vertex_id];

        if (vtx_num == vtx_num_1) {
          if (cell_vtx_quad[face_id*4 + ((vertex_id+1) % 4)] == vtx_num_2) {
            direction = 1;
            break;
          }
          else if (   cell_vtx_quad[face_id*4 + ((vertex_id-1+4) % 4)]
                   == vtx_num_2) {
            direction = -1;
            break;
          }
        }

      }

      if (direction != 0)
        break;

    }

    if (direction == -1)
      warn_orient = true;
    else if (direction == 0)
      return FVMC_NODAL_FROM_DESC_FAILURE;

    cell_vtx_prism[3 + ipass]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (3 * direction))) % 4)];

    cell_vtx_prism[4 + ipass]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (2 * direction))) % 4)];

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("prism : %d %d %d %d %d %d\n",
             cell_vtx_prism[0], cell_vtx_prism[1], cell_vtx_prism[2],
             cell_vtx_prism[3], cell_vtx_prism[4], cell_vtx_prism[5]);
#endif

  if (warn_orient == true)
    return FVMC_NODAL_FROM_DESC_WORIENT;
  else
    return FVMC_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Hexahedron construction
 *
 * parameters:
 *   cell_vtx_quad <-- quadrangle faces connectivity
 *   cell_vtx_hexa --> hexahedron connectivity (pre-allocated)
 *----------------------------------------------------------------------------*/

inline static fvmc_nodal_from_desc_t
_nodal_from_desc_cnv_cel_hexa(const fvmc_lnum_t  cell_vtx_quad[],
                              fvmc_lnum_t  cell_vtx_hexa[])
{
  fvmc_lnum_t  vertex_id, face_id;
  fvmc_lnum_t  ipass, direction;
  fvmc_lnum_t  vtx_num, vtx_num_1, vtx_num_2;

  _Bool  warn_orient = false;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("face 1 : %d %d %d %d\n",
             cell_vtx_quad[0], cell_vtx_quad[1],
             cell_vtx_quad[2], cell_vtx_quad[3]);
  bftc_printf("face 2 : %d %d %d %d\n",
             cell_vtx_quad[4], cell_vtx_quad[5],
             cell_vtx_quad[6], cell_vtx_quad[7]);
  bftc_printf("face 3 : %d %d %d %d\n",
             cell_vtx_quad[8], cell_vtx_quad[9],
             cell_vtx_quad[10], cell_vtx_quad[11]);
  bftc_printf("face 4 : %d %d %d %d\n",
             cell_vtx_quad[12], cell_vtx_quad[13],
             cell_vtx_quad[14], cell_vtx_quad[15]);
  bftc_printf("face 5 : %d %d %d %d\n",
             cell_vtx_quad[16], cell_vtx_quad[17],
             cell_vtx_quad[18], cell_vtx_quad[19]);
  bftc_printf("face 6 : %d %d %d %d\n",
             cell_vtx_quad[20], cell_vtx_quad[21],
             cell_vtx_quad[22], cell_vtx_quad[23]);
#endif

  /*
    Base : vertices of the fisrt face; we take the vertices in opposite
    order ("bottom" face numbering of the pyramid with outward
    pointing normal).

    *     8 x-------x 7
    *      /|      /|
    *     / |     / |
    *  5 x-------x6 |
    *    | 4x----|--x 3
    *    | /     | /
    *    |/      |/
    *  1 x-------x 2
  */


  cell_vtx_hexa[0] = cell_vtx_quad[3];
  cell_vtx_hexa[1] = cell_vtx_quad[2];
  cell_vtx_hexa[2] = cell_vtx_quad[1];
  cell_vtx_hexa[3] = cell_vtx_quad[0];

  /*
    We have found 4 out of 8 vertices; we first seek the quadrangle  sharing
    vertices 1 and 2, so as to determine vertices 5 and 6,
    then we seek the quadrangle sharing vertices 3 and 4, so as to
    determine vertices 7 and 8.
  */

  for (ipass = 0 ; ipass < 2 ; ipass++) {

    vtx_num_1 = cell_vtx_hexa[     ipass * 2 ];
    vtx_num_2 = cell_vtx_hexa[1 + (ipass * 2)];

    direction = 0;

    for (face_id = 1 ; face_id < 6 ; face_id++) {

      for (vertex_id = 0 ; vertex_id < 4 ; vertex_id++) {

        vtx_num = cell_vtx_quad[face_id*4 + vertex_id];

        if (vtx_num == vtx_num_1) {
          if (cell_vtx_quad[face_id*4 + ((vertex_id+1) % 4)] == vtx_num_2) {
            direction = 1;
            break;
          }
          else if (   cell_vtx_quad[face_id*4 + ((vertex_id-1+4) % 4)]
                   == vtx_num_2) {
            direction = -1;
            break;
          }
        }

      }

      if (direction != 0)
        break;

    }

    if (direction == -1)
      warn_orient = true;
    else if (direction == 0)
      return FVMC_NODAL_FROM_DESC_FAILURE;

    cell_vtx_hexa[4 + (ipass * 2)]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (3 * direction))) % 4)];

    cell_vtx_hexa[5 + (ipass * 2)]
      = cell_vtx_quad[face_id*4 + ((vertex_id + (4 + (2 * direction))) % 4)];

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("hexa : %d %d %d %d %d %d %d %d\n",
             cell_vtx_hexa[0], cell_vtx_hexa[1], cell_vtx_hexa[2],
             cell_vtx_hexa[3], cell_vtx_hexa[4], cell_vtx_hexa[5],
             cell_vtx_hexa[6], cell_vtx_hexa[7]);
#endif

  if (warn_orient == true)
    return FVMC_NODAL_FROM_DESC_WORIENT;
  else
    return FVMC_NODAL_FROM_DESC_SUCCESS;

}

/*----------------------------------------------------------------------------
 * Determination of the number of vertices defining a given face.
 *
 * parameters:
 *   face_id,        <-- face id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *
 * returns:
 *   number of vertices of face defined by face_id
 *----------------------------------------------------------------------------*/

inline static fvmc_lnum_t
_nodal_face_from_desc_size(const fvmc_lnum_t   face_id,
                           const int          n_face_lists,
                           const fvmc_lnum_t   face_list_shift[],
                           const fvmc_lnum_t  *face_vertex_idx[])
{
  fvmc_lnum_t  fl, _face_id;
  fvmc_lnum_t  n_face_vertices;
  fvmc_lnum_t  vertex_id_start, vertex_id_end;

  /* Compute number of vertices */
  /*----------------------------*/

  _face_id = face_id;
  for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
  assert(fl > -1);
  _face_id -= face_list_shift[fl];

  vertex_id_start = face_vertex_idx[fl][_face_id] - 1;
  vertex_id_end   = face_vertex_idx[fl][_face_id + 1] - 1;
  n_face_vertices = vertex_id_end - vertex_id_start;

  return n_face_vertices;

}

/*----------------------------------------------------------------------------
 * Copy the vertex numbers defining a given face.
 *
 * The face_vtx[] array is filled with the vertex numbers of the face,
 * and should be large enough to receive them.
 *
 * parameters:
 *   face_id,        <-- face id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   face_vtx        --> local face definition
 *----------------------------------------------------------------------------*/

inline static void
_nodal_face_from_desc_copy(const fvmc_lnum_t   face_id,
                           const int          n_face_lists,
                           const fvmc_lnum_t   face_list_shift[],
                           const fvmc_lnum_t  *face_vertex_idx[],
                           const fvmc_lnum_t  *face_vertex_num[],
                           fvmc_lnum_t        *face_vtx)
{
  fvmc_lnum_t  fl, _face_id;
  fvmc_lnum_t  vtx, vertex_id, vertex_id_start, vertex_id_end;

  /* Copy vertex numbers */
  /*---------------------*/

  _face_id = face_id;
  for (fl = n_face_lists - 1 ; _face_id < face_list_shift[fl] ; fl--);
  assert(fl > -1);
  _face_id -= face_list_shift[fl];

  vertex_id_start = face_vertex_idx[fl][_face_id] - 1;
  vertex_id_end   = face_vertex_idx[fl][_face_id + 1] - 1;

  for (vtx = 0, vertex_id = vertex_id_start ;
       vertex_id < vertex_id_end ;
       vtx++, vertex_id++)
    face_vtx[vtx] = face_vertex_num[fl][vertex_id];
}

/*----------------------------------------------------------------------------
 * Polyhedral cells connectivity extraction.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_polys         <-- size of list_poly[]
 *   list_poly       <-- list of polyhedral cells (cell index, 1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_face_list  --> numbers of faces defining polyhedra
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_extract_polyhedra(fvmc_nodal_section_t  *this_section,
                             const fvmc_lnum_t      n_polys,
                             const fvmc_lnum_t      list_poly[],
                             const int             n_face_lists,
                             const fvmc_lnum_t      face_list_shift[],
                             const fvmc_lnum_t     *face_vertex_idx[],
                             const fvmc_lnum_t     *face_vertex_num[],
                             const fvmc_lnum_t      cell_face_idx[],
                             const fvmc_lnum_t      cell_face_num[],
                             fvmc_lnum_t           *cell_face_list[])
{
  fvmc_lnum_t   n_faces, n_cell_faces;
  fvmc_lnum_t   c_cell_face_vertex_idxs, c_cell_face_vertex_nums;
  fvmc_lnum_t   n_cell_face_vertex_nums;
  fvmc_lnum_t   face_counter, face_id, poly_id, cell_id;
  fvmc_lnum_t   fl, i, idx, idx_start, idx_end, num_count;

  fvmc_lnum_t  *local_face_num = NULL;

  int sgn;


  /* Indicators initialization */

  n_faces = face_list_shift[n_face_lists] - face_list_shift[0];

  BFTC_MALLOC(local_face_num, n_faces, fvmc_lnum_t);

  for (face_id = 0 ; face_id < n_faces ; face_id++)
    local_face_num[face_id] = 0;

  /* Flagging of referenced faces and Cells -> Faces indexes */
  /*---------------------------------------------------------*/

  n_cell_faces = 0;

  BFTC_MALLOC(this_section->_face_index, n_polys + 1, fvmc_lnum_t);
  this_section->face_index = this_section->_face_index;
  this_section->_face_index[0] = 0;

  for (poly_id = 0 ; poly_id < n_polys ; poly_id++) {

    cell_id = list_poly[poly_id] - 1;

    idx_start = cell_face_idx[cell_id]     - 1;
    idx_end   = cell_face_idx[cell_id  +1] - 1;

    this_section->_face_index[poly_id + 1]
      = this_section->_face_index[poly_id] + (idx_end - idx_start);

    for (idx = idx_start ; idx < idx_end ; idx++) {

      face_id = FVMC_ABS(cell_face_num[idx]) - 1;

      /* Mark only used values for now, local_face_num[] values
         will be computed later by looping on faces so that
         for nonzero values, local_face_num[i] > local_face_num[j]
         if i > j ; this is important for later parts of this algorithm */

      if (local_face_num[face_id] == 0)
        local_face_num[face_id] = 1;

    }

  }

  /* Counting for faces -> vertices connectivity and local face numbering */

  n_cell_face_vertex_nums = 0;
  for (face_counter = 0 ; face_counter < n_faces ; face_counter++) {

    if (local_face_num[face_counter] != 0) {

      /* Transform local_face_num[] from a marker to a renumbering array */

      n_cell_faces += 1;
      local_face_num[face_counter] = n_cell_faces;

      /* Counting for faces -> vertices connectivity */

      face_id = face_counter;
      for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
      assert(fl > -1);
      face_id -= face_list_shift[fl];

      n_cell_face_vertex_nums += (  face_vertex_idx[fl][face_id+1]
                                  - face_vertex_idx[fl][face_id]);

    }
  }

  /* Cells -> Faces Connectivity (face numbers) */
  /*--------------------------------------------*/

  BFTC_MALLOC(this_section->_face_num,
             this_section->_face_index[n_polys],
             fvmc_lnum_t);
  this_section->face_num = this_section->_face_num;

  num_count = 0;

  for (poly_id = 0 ; poly_id < n_polys ; poly_id++) {

    cell_id = list_poly[poly_id] - 1;

    idx_start = cell_face_idx[cell_id]     - 1;
    idx_end   = cell_face_idx[cell_id  +1] - 1;

    for (idx = idx_start ; idx < idx_end ; idx++) {

      face_id = FVMC_ABS(cell_face_num[idx]) - 1;
      sgn = (cell_face_num[idx] > 0) ? 1 : -1;

      this_section->_face_num[num_count++]
        = sgn * local_face_num[face_id] ;

    }

  }

  assert(num_count == this_section->_face_index[n_polys]);

  /* Faces -> Vertices Connectivity */
  /*--------------------------------*/

  if (cell_face_list != NULL)
    BFTC_MALLOC(*cell_face_list, n_cell_faces, fvmc_lnum_t);

  BFTC_MALLOC(this_section->_vertex_index,
             n_cell_faces + 1,
             fvmc_lnum_t);
  this_section->vertex_index = this_section->_vertex_index;

  BFTC_MALLOC(this_section->_vertex_num,
             n_cell_face_vertex_nums,
             fvmc_lnum_t);
  this_section->vertex_num = this_section->_vertex_num;

  /* Definition */

  c_cell_face_vertex_idxs = 0;
  c_cell_face_vertex_nums = 0;
  this_section->_vertex_index[0] = 0;

  for (face_counter = 0 ; face_counter < n_faces ; face_counter++) {

    if (local_face_num[face_counter] != 0) {

      if (cell_face_list != NULL)
        (*cell_face_list)[local_face_num[face_counter] -1 ] = face_counter + 1;

      face_id = face_counter;
      for (fl = n_face_lists - 1 ; face_id < face_list_shift[fl] ; fl--);
      assert(fl < n_face_lists);
      face_id -= face_list_shift[fl];

      for (i = face_vertex_idx[fl][face_id] - 1 ;
           i < face_vertex_idx[fl][face_id + 1] - 1 ;
           i++)
        this_section->_vertex_num[c_cell_face_vertex_nums++]
          = face_vertex_num[fl][i];

      c_cell_face_vertex_idxs += 1;
      this_section->_vertex_index[c_cell_face_vertex_idxs]
        = c_cell_face_vertex_nums;

    }
  }

  /* No past-the-end index value counted,
     so counter = this_nodal->n_cell_faces */

  assert(c_cell_face_vertex_idxs == n_cell_faces);
  assert(c_cell_face_vertex_nums == n_cell_face_vertex_nums);

  /* Set pointers and connectivity size */

  this_section->connectivity_size = n_cell_face_vertex_nums;

  this_section->n_faces      = n_cell_faces;
  this_section->face_index   = this_section->_face_index;
  this_section->face_num     = this_section->_face_num;
  this_section->vertex_index = this_section->_vertex_index;
  this_section->vertex_num   = this_section->_vertex_num;

  /* Free memory */

  BFTC_FREE(local_face_num);

}

/*----------------------------------------------------------------------------
 * Raise parent element numbering in given sections by one level, as
 * defined by the parent_element_num[] arrays.
 *
 * This is useful when a nodal mesh is extracted from a temporary subset of
 * a parent mesh (corresponding to the parent_element_num[] list, using
 * 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   n_sections         <-- size of sections array
 *   sections           <-- array of sections to reduce
 *   parent_element_num <-- element -> parent element number (1 to n) if
 *                          non-trivial (i.e. if element definitions
 *                          correspond to a subset of the parent mesh),
 *                          NULL otherwise.
 *----------------------------------------------------------------------------*/

static void
_raise_sections_parent_num(const int             n_sections,
                           fvmc_nodal_section_t  *sections[],
                           const fvmc_lnum_t      parent_element_num[])
{
  int  section_id;
  fvmc_lnum_t  element_counter;

  fvmc_nodal_section_t  *section;

  if (parent_element_num == NULL)
    return;

  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL) {
      if (section->_parent_element_num == NULL) {
        BFTC_MALLOC(section->_parent_element_num,
                   section->n_elements,
                   fvmc_lnum_t);
        section->parent_element_num = section->_parent_element_num;
      }
      for (element_counter = 0 ;
           element_counter < section->n_elements ;
           element_counter++)
        section->_parent_element_num[element_counter]
          = parent_element_num[section->parent_element_num[element_counter]
                               - 1];
    }
  }

}

/*----------------------------------------------------------------------------
 * Free parent element numbering arrays when not needed.
 *
 * If the parent mesh is already complete and ordered for a given section
 * (i.e. not based on a partial extraction or then based on the n first
 * elements in the same order), the parent element number is not needed,
 * and may be freed.
 *
 * parameters:
 *   sections        <-- array of sections to reduce
 *   n_sections      <-- size of sections array
 *----------------------------------------------------------------------------*/

static void
_optimize_sections_parent_num(const int             n_sections,
                              fvmc_nodal_section_t  *sections[])
{
  int  section_id;
  fvmc_lnum_t  element_counter;

  fvmc_nodal_section_t  *section;

  /* If the parent mesh is already complete for a given element type
     (i.e. not based on a partial extraction or then based on the
     n first elements), the parent cell number is not needed */

  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL) {
      for (element_counter = 0 ;
           element_counter < section->n_elements ;
           element_counter++) {
        if (section->parent_element_num[element_counter] != element_counter + 1)
          break;
      }
      if (element_counter == section->n_elements) {
        section->parent_element_num = NULL;
        BFTC_FREE(section->_parent_element_num);
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Add nodal mesh structure sections to a nodal mesh.
 *
 * Sections to add are defined by an array of section pointers,
 * which may contain NULL entries. Only the non-NULL entries are
 * added to the nodal mesh structure, and they belong to that structure
 * from this point on (i.e. their lifecycle is based upon that of
 * the nodal mesh structure).
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_sections      <-- size of sections array
 *   sections        <-- array of sections to add
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_add_sections(fvmc_nodal_t          *this_nodal,
                        const int             n_sections,
                        fvmc_nodal_section_t  *sections[])
{
  int  section_id, section_count;

  fvmc_nodal_section_t  *section;

  /* Add sections to nodal mesh structure */

  section_count = 0;
  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL)
      section_count++;
  }

  BFTC_REALLOC(this_nodal->sections,
              this_nodal->n_sections + section_count,
              fvmc_nodal_section_t *);

  section_count = 0;
  for (section_id = 0 ; section_id < n_sections ; section_id++) {
    section = sections[section_id];
    if (section != NULL)
      this_nodal->sections[this_nodal->n_sections + section_count++]
        = section;
  }
  this_nodal->n_sections += section_count;

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert and add cells from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_cells[] argument is non-NULL, cells
 * {extr_cells[0], extr_cells[1], extr_cells[n_extr_cells - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, cells
 * {1, 2, ..., n_extr_cells} are considered.
 *
 * In addition, an optional parent_cell_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_cell_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_cells    <-- count of cells to add
 *   extr_cells      <-- optional filter list of cells to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists + 1
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   parent_cell_num <-- cell -> parent cell number (1 to n) if non-trivial
 *                       (i.e. if cell definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *   cell_face_list  --> numbers of faces defining polyhedra
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_from_desc_add_cells(fvmc_nodal_t        *this_nodal,
                              const fvmc_lnum_t    n_extr_cells,
                              const fvmc_lnum_t    extr_cells[],
                              const int           n_face_lists,
                              const fvmc_lnum_t    face_list_shift[],
                              const fvmc_lnum_t   *face_vertex_idx[],
                              const fvmc_lnum_t   *face_vertex_num[],
                              const fvmc_lnum_t    cell_face_idx[],
                              const fvmc_lnum_t    cell_face_num[],
                              const fvmc_lnum_t    parent_cell_num[],
                              fvmc_lnum_t         *cell_face_list[])
{
  int  type_id;
  fvmc_lnum_t  cell_counter, cell_id;

  fvmc_element_t  cell_type;

  fvmc_lnum_t  cell_vtx_tria[3*4]; /* We will only seek to fill these arrays */
  fvmc_lnum_t  cell_vtx_quad[4*6]; /* for local faces 1-4 and 1-6 at most    */
  fvmc_lnum_t  *p_cell_vertex;

  fvmc_lnum_t  n_elements_type[FVMC_N_ELEMENT_TYPES];
  fvmc_gnum_t  n_g_elements_type[FVMC_N_ELEMENT_TYPES];

  fvmc_nodal_section_t  *section;
  fvmc_nodal_section_t  *sections[FVMC_N_ELEMENT_TYPES];

  fvmc_gnum_t  n_orient_pbs = 0; /* Number of cells with potential (non-
                                   fatal) orientation problem */

  fvmc_nodal_from_desc_t  retcode;

  /* Initialization */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    n_elements_type[type_id] = 0;
    sections[type_id] = NULL;
  }

  /* Guess connectivity types */
  /*--------------------------*/

  for (cell_counter = 0 ; cell_counter < n_extr_cells ; cell_counter++) {

    if (extr_cells != NULL)
      cell_id = extr_cells[cell_counter] - 1;
    else
      cell_id = cell_counter;

    cell_type = _nodal_cell_from_desc(cell_id,
                                      n_face_lists,
                                      face_list_shift,
                                      face_vertex_idx,
                                      face_vertex_num,
                                      cell_face_idx,
                                      cell_face_num,
                                      NULL,
                                      NULL);

    n_elements_type[cell_type] += 1;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("dbg : fvmc_nodal_cells_from_desc\n"
             "n_tetra = %d\n"
             "n_pyram = %d\n"
             "n_prism = %d\n"
             "n_hexa  = %d\n"
             "n_poly  = %d\n",
             n_elements_type[FVMC_CELL_TETRA],
             n_elements_type[FVMC_CELL_PYRAM],
             n_elements_type[FVMC_CELL_PRISM],
             n_elements_type[FVMC_CELL_HEXA],
             n_elements_type[FVMC_CELL_POLY]);
#endif

  /* Set dimensions (and reset local counters) */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++)
    n_g_elements_type[type_id] = n_elements_type[type_id];

  fvmc_parall_counter(n_g_elements_type, FVMC_N_ELEMENT_TYPES);

  int cell_type_i_d = (int) FVMC_CELL_TETRA;
  int cell_type_i_f = (int) FVMC_CELL_POLY;
  int cell_type_i;

  for (cell_type_i = cell_type_i_d;
       cell_type_i <= cell_type_i_f;
       cell_type_i++) {
    if (n_g_elements_type[cell_type_i] > 0) {
      int order = 1;
      sections[cell_type_i] = fvmc_nodal_section_create((fvmc_element_t) cell_type_i, order);
      sections[cell_type_i]->n_elements = n_elements_type[cell_type_i];
      this_nodal->n_cells += n_elements_type[cell_type_i];
    }
    n_elements_type[cell_type_i] = 0;
  }

  /* Main memory allocations */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      if (   section->type != FVMC_FACE_POLY
          && section->type != FVMC_CELL_POLY) {
        section->stride = fvmc_nodal_n_vertices_element((fvmc_element_t) type_id, section->order);
        section->connectivity_size = section->stride * section->n_elements;
        BFTC_MALLOC(section->_vertex_num, section->connectivity_size, fvmc_lnum_t);
        section->vertex_num = section->_vertex_num;
      }
    }
  }

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      BFTC_MALLOC(section->_parent_element_num, section->n_elements, fvmc_lnum_t);
      section->parent_element_num = section->_parent_element_num;
    }
  }

  /* Construction of nodal connectivities */
  /*---------------------------------------*/

  for (cell_counter = 0 ; cell_counter < n_extr_cells ; cell_counter++) {

    if (extr_cells != NULL)
      cell_id = extr_cells[cell_counter] - 1;
    else
      cell_id = cell_counter;

    cell_type = _nodal_cell_from_desc(cell_id,
                                      n_face_lists,
                                      face_list_shift,
                                      face_vertex_idx,
                                      face_vertex_num,
                                      cell_face_idx,
                                      cell_face_num,
                                      cell_vtx_tria,
                                      cell_vtx_quad);

    section = sections[cell_type];

    p_cell_vertex =   section->_vertex_num
                    + (  n_elements_type[cell_type]
                       * fvmc_nodal_n_vertices_element(cell_type, section->order));

    switch (cell_type) {
    case FVMC_CELL_TETRA:
      retcode = _nodal_from_desc_cnv_cel_tetra(cell_vtx_tria,
                                               p_cell_vertex);
      break;
    case FVMC_CELL_PYRAM:
      retcode = _nodal_from_desc_cnv_cel_pyram(cell_vtx_tria,
                                               cell_vtx_quad,
                                               p_cell_vertex);
      break;
    case FVMC_CELL_PRISM:
      retcode = _nodal_from_desc_cnv_cel_prism(cell_vtx_tria,
                                               cell_vtx_quad,
                                               p_cell_vertex);
      break;
    case FVMC_CELL_HEXA:
      retcode = _nodal_from_desc_cnv_cel_hexa(cell_vtx_quad,
                                              p_cell_vertex);
      break;
    default:
      retcode = FVMC_NODAL_FROM_DESC_SUCCESS;
      break;
    }

    /* Temporary value of parent cell num based on local cell id
       (so that the list of polyhedra given as an argument to
       _fvmc_nodal_extract_polyhedra() below is based on the local
       numbering, like the cell->face connectivity */

    section->_parent_element_num[n_elements_type[cell_type]] = cell_id + 1;

    n_elements_type[cell_type] += 1;

    if (retcode == FVMC_NODAL_FROM_DESC_WORIENT)
      n_orient_pbs += 1;
    else
      if (retcode == FVMC_NODAL_FROM_DESC_FAILURE)
        bftc_error(__FILE__, __LINE__, 0,
                  _("Incoherent connectivity for cell %d\n"),
                  cell_id + 1);

  }

  fvmc_parall_counter(&n_orient_pbs, 1);

  if (n_orient_pbs > 0)
    bftc_printf("Warning: Possible nodal connectivity orientation\n"
               "         problems for at least %d cells\n",
               n_orient_pbs);

  /* Extraction of remaining polyhedra */
  /*-----------------------------------*/

  if (sections[FVMC_CELL_POLY] != NULL)
    _fvmc_nodal_extract_polyhedra
      (sections[FVMC_CELL_POLY],
       n_elements_type[FVMC_CELL_POLY],
       sections[FVMC_CELL_POLY]->parent_element_num,
       n_face_lists,
       face_list_shift,
       face_vertex_idx,
       face_vertex_num,
       cell_face_idx,
       cell_face_num,
       cell_face_list);

  /* We can now base the final value of the parent cell number on
     the parent (and not local) numbering */

  _raise_sections_parent_num(FVMC_N_ELEMENT_TYPES, sections, parent_cell_num);

  /* Add sections to nodal mesh structure */

  _optimize_sections_parent_num(FVMC_N_ELEMENT_TYPES, sections);

  _fvmc_nodal_add_sections(this_nodal, FVMC_N_ELEMENT_TYPES, sections);

}

/*----------------------------------------------------------------------------
 * Convert and add faces from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_faces[] argument is non-NULL, faces
 * {extr_faces[0], extr_faces[1], extr_faces[n_extr_faces - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, faces
 * {1, 2, ..., n_extr_faces} are considered.
 *
 * In addition, an optional parent_face_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_face_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_faces    <-- count of faces to add
 *   extr_faces      <-- optional filter list of faces to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex_num <-- face -> vertex numbers (per face list)
 *   parent_face_num <-- face -> parent face number (1 to n) if non-trivial
 *                       (i.e. if face definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_from_desc_add_faces(fvmc_nodal_t        *this_nodal,
                              const fvmc_lnum_t    n_extr_faces,
                              const fvmc_lnum_t    extr_faces[],
                              const int           n_face_lists,
                              const fvmc_lnum_t    face_list_shift[],
                              const fvmc_lnum_t   *face_vertex_idx[],
                              const fvmc_lnum_t   *face_vertex_num[],
                              const fvmc_lnum_t    parent_face_num[])
{
  int  type_id;
  fvmc_lnum_t  face_counter, face_id;

  fvmc_element_t  face_type;

  fvmc_lnum_t  *p_vertex_idx;
  fvmc_lnum_t  *p_vertex_num;

  fvmc_lnum_t  poly_connect_size;
  fvmc_lnum_t  n_face_vertices;
  fvmc_lnum_t  n_elements_type[FVMC_N_ELEMENT_TYPES];
  fvmc_gnum_t  n_g_elements_type[FVMC_N_ELEMENT_TYPES];

  fvmc_nodal_section_t  *section;
  fvmc_nodal_section_t  *sections[FVMC_N_ELEMENT_TYPES];

  /* Initialization */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    n_elements_type[type_id] = 0;
    sections[type_id] = NULL;
  }
  poly_connect_size = 0;

  /* Compute connectivity type */
  /*---------------------------*/

  for (face_counter = 0 ; face_counter < n_extr_faces ; face_counter++) {

    if (extr_faces != NULL)
      face_id = extr_faces[face_counter] - 1;
    else
      face_id = face_counter;

    n_face_vertices = _nodal_face_from_desc_size(face_id,
                                                 n_face_lists,
                                                 face_list_shift,
                                                 face_vertex_idx);

    switch (n_face_vertices) {
    case 3:
      face_type = FVMC_FACE_TRIA;
      break;
    case 4:
      face_type = FVMC_FACE_QUAD;
      break;
    default:
      face_type = FVMC_FACE_POLY;
      poly_connect_size += n_face_vertices;
      break;
    }

    n_elements_type[face_type] += 1;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bftc_printf("dbg : fvmc_nodal_faces_from_desc\n"
             "n_tria = %d\n"
             "n_quad = %d\n"
             "n_poly = %d\n",
             n_elements_type[FVMC_FACE_TRIA],
             n_elements_type[FVMC_FACE_QUAD],
             n_elements_type[FVMC_FACE_POLY]);
#endif

  /* Set dimensions (and reset local counters) */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++)
    n_g_elements_type[type_id] = n_elements_type[type_id];

  fvmc_parall_counter(n_g_elements_type, FVMC_N_ELEMENT_TYPES);

  int face_type_i_d = (int) FVMC_FACE_TRIA;
  int face_type_i_f = (int) FVMC_FACE_POLY;
  int face_type_i;

  for (face_type_i = face_type_i_d;
       face_type_i <= face_type_i_f;
       face_type_i++) {
    if (n_g_elements_type[face_type_i] > 0) {
      int order = 1;
      sections[face_type_i] = fvmc_nodal_section_create((fvmc_element_t) face_type_i, order);
      sections[face_type_i]->n_elements = n_elements_type[face_type_i];
      this_nodal->n_faces += n_elements_type[face_type_i];
    }
    n_elements_type[face_type_i] = 0;
  }

  /* Main memory allocations */

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      if (section->type != FVMC_FACE_POLY) {
        section->stride = fvmc_nodal_n_vertices_element((fvmc_element_t)type_id, section->order);
        section->connectivity_size = section->stride * section->n_elements;
        BFTC_MALLOC(section->_vertex_num, section->connectivity_size, fvmc_lnum_t);
        section->vertex_num = section->_vertex_num;
      }
      else {
        section->stride = fvmc_nodal_n_vertices_element((fvmc_element_t)type_id, section->order);
        section->connectivity_size = poly_connect_size;
        BFTC_MALLOC(section->_vertex_index, section->n_elements + 1, fvmc_lnum_t);
        BFTC_MALLOC(section->_vertex_num, section->connectivity_size, fvmc_lnum_t);
        section->vertex_index = section->_vertex_index;
        section->vertex_num = section->_vertex_num;
        section->_vertex_index[0] = 0;
      }
    }
  }

  for (type_id = 0 ; type_id < FVMC_N_ELEMENT_TYPES ; type_id++) {
    section = sections[type_id];
    if (section != NULL) {
      BFTC_MALLOC(section->_parent_element_num, section->n_elements, fvmc_lnum_t);
      section->parent_element_num = section->_parent_element_num;
    }
  }

  /* Construction of nodal connectivities */
  /*---------------------------------------*/

  for (face_counter = 0 ; face_counter < n_extr_faces ; face_counter++) {

    if (extr_faces != NULL)
      face_id = extr_faces[face_counter] - 1;
    else
      face_id = face_counter;

    n_face_vertices = _nodal_face_from_desc_size(face_id,
                                                 n_face_lists,
                                                 face_list_shift,
                                                 face_vertex_idx);

    switch (n_face_vertices) {
    case 3:
      face_type = FVMC_FACE_TRIA;
      section = sections[face_type];
      p_vertex_num = section->_vertex_num + (n_elements_type[face_type] * 3);
      break;
    case 4:
      face_type = FVMC_FACE_QUAD;
      section = sections[face_type];
      p_vertex_num = section->_vertex_num + (n_elements_type[face_type] * 4);
      break;
    default:
      face_type = FVMC_FACE_POLY;
      section = sections[face_type];
      p_vertex_idx = section->_vertex_index + n_elements_type[face_type];
      *(p_vertex_idx + 1) = *p_vertex_idx + n_face_vertices;
      p_vertex_num = section->_vertex_num + (*p_vertex_idx);
      break;
    }

    _nodal_face_from_desc_copy(face_id,
                               n_face_lists,
                               face_list_shift,
                               face_vertex_idx,
                               face_vertex_num,
                               p_vertex_num);

    section->_parent_element_num[n_elements_type[face_type]] = face_id + 1;

    n_elements_type[face_type] += 1;

  }

  /* We can now base the final value of the parent face number on
     the parent (and not local) numbering */

  _raise_sections_parent_num(FVMC_N_ELEMENT_TYPES, sections, parent_face_num);

  /* Add sections to nodal mesh structure */

  _optimize_sections_parent_num(FVMC_N_ELEMENT_TYPES, sections);

  _fvmc_nodal_add_sections(this_nodal, FVMC_N_ELEMENT_TYPES, sections);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
