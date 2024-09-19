/*============================================================================
 * Ordering of nodal mesh entity lists and connectivity
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2009  EDF

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
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_order.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal_order.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create ordered parent entity list or order existing list.
 *
 * parameters:
 *   _list  <-> pointer to optional list (1 to n numbering) of selected
 *              entities. An existing list is ordered, otherwise one is
 *              created.
 *   list   <-> pointer tooptional list (1 to n numbering) of selected
 *              entities. A shared list is copied to _list and ordered,
 *              a private list (pointed to by _list) is simply ordered.
 *   order  <-- ordering of entities (0 to n-1).
 *   nb_ent <-- number of entities considered.
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_order_parent_list(fvmc_lnum_t        * _list[],
                             const fvmc_lnum_t  * list[],
                             const fvmc_lnum_t    order[],
                             const size_t        nb_ent)
{
  size_t  i;

  fvmc_lnum_t  *ordered_list = NULL;

  BFTC_MALLOC(ordered_list, nb_ent, fvmc_lnum_t);

  if (*list != NULL) {

    for (i = 0 ; i < nb_ent ; i++)
      ordered_list[i] = (*list)[order[i]];

    if (*_list != NULL) {
      for (i = 0 ; i < nb_ent ; i++)
        (*_list)[i] = ordered_list[i];
      BFTC_FREE(ordered_list);
    }
    else
      *_list = ordered_list;

  }
  else {

    assert(*list == NULL);

    for (i = 0 ; i < nb_ent ; i++)
      ordered_list[i] = order[i] + 1;
    *_list = ordered_list;

  }

  *list = *_list;
}

/*----------------------------------------------------------------------------
 * Reorder strided connectivity array.
 *
 * parameters:
 *   connect <-> connectivity array (nb_ent * stride) to be ordered.
 *   order   <-- ordering of entities (0 to n-1).
 *   nb_ent  <-- number of entities considered.
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_order_strided_connect(fvmc_lnum_t          connect[],
                                 const fvmc_lnum_t    order[],
                                 const size_t        stride,
                                 const size_t        nb_ent)
{
  size_t  i, j;
  fvmc_lnum_t  *p1, *p2;

  fvmc_lnum_t  *tmp_connect = NULL;

  BFTC_MALLOC(tmp_connect, nb_ent * stride, fvmc_lnum_t);

  /* Temporary ordered copy */

  for (i = 0 ; i < nb_ent ; i++) {
    p1 = tmp_connect + i*stride;
    p2 = connect + (order[i] * stride);
    for (j = 0 ; j < stride ; j++)
      *p1++ = *p2++;
  }

  /* Now put back in initial location */

  memcpy(connect, tmp_connect, stride * nb_ent * sizeof(fvmc_lnum_t));

  BFTC_FREE(tmp_connect);
}

/*----------------------------------------------------------------------------
 * Reorder indexed connectivity array.
 *
 * parameters:
 *   connect_idx <-> connectivity index array (0 to n -1) to be ordered.
 *   connect_num <-> connectivity numbers array to be ordered.
 *   order       <-- ordering of entities (0 to n-1).
 *   nb_ent      <-- number of entities considered.
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_order_indexed_connect(fvmc_lnum_t          connect_idx[],
                                 fvmc_lnum_t          connect_num[],
                                 const fvmc_lnum_t    order[],
                                 const size_t        nb_ent)
{
  size_t  i, j, nb_ent_max, nb_loc;
  fvmc_lnum_t  *p1, *p2;

  fvmc_lnum_t  *tmp_connect = NULL;

  nb_ent_max = connect_idx[nb_ent] ; /* size of connect_num */
  if (nb_ent > nb_ent_max) /* only if some entities have no connectivity */
    nb_ent_max = nb_ent;

  BFTC_MALLOC(tmp_connect, nb_ent_max, fvmc_lnum_t);

  /* Temporary ordered copy of values */

  p1 = tmp_connect;
  for (i = 0 ; i < nb_ent ; i++) {
    nb_loc = connect_idx[order[i]+1] - connect_idx[order[i]];
    p2 = connect_num + connect_idx[order[i]];
    for (j = 0 ; j < nb_loc ; j++)
      *p1++ = *p2++;
  }

  /* Now put back in initial location */

  memcpy(connect_num, tmp_connect,
         (size_t)(connect_idx[nb_ent]) * sizeof(fvmc_lnum_t));

  /* Index to size : size associated with entity i in position i+1 */

  for (i = nb_ent ; i > 0 ; i--)
    connect_idx[i] = connect_idx[i] - connect_idx[i-1];

  /* Temporary ordered copy of transformed index */

  p1 = tmp_connect;
  *p1++ = 0;
  for (i = 0 ; i < nb_ent ; i++)
    *p1++ = connect_idx[order[i] + 1];

  /* Put back in initial location and re-convert to index*/

  memcpy(connect_idx, tmp_connect, (size_t)(nb_ent + 1) * sizeof(fvmc_lnum_t));

  for (i = 0 ; i < nb_ent ; i++)
    connect_idx[i+1] = connect_idx[i+1] + connect_idx[i];

  BFTC_FREE(tmp_connect);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Locally order cells and associated connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent cells (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_order_cells(fvmc_nodal_t       *this_nodal,
                      const fvmc_gnum_t   parent_global_number[])
{
  fvmc_lnum_t  i;
  fvmc_lnum_t  *order = NULL;
  fvmc_nodal_section_t  *section = NULL;

  if (this_nodal == NULL)
    return;

  /* Order locally if necessary */

  for (i = 0 ; i < this_nodal->n_sections ; i++) {

    section = this_nodal->sections[i];

    if (section->entity_dim == 3) {

      assert(section->global_element_num == NULL);

      if (fvmc_order_local_test(section->parent_element_num,
                               parent_global_number,
                               section->n_elements) == false) {

        order = fvmc_order_local(section->parent_element_num,
                                parent_global_number,
                                section->n_elements);

        _fvmc_nodal_order_parent_list(&(section->_parent_element_num),
                                     &(section->parent_element_num),
                                     order,
                                     section->n_elements);

        if (section->type != FVMC_CELL_POLY) {
          fvmc_nodal_section_copy_on_write(section, false, false, false, true);
          _fvmc_nodal_order_strided_connect(section->_vertex_num,
                                           order,
                                           (size_t)(section->stride),
                                           section->n_elements);
        }
        else {
          fvmc_nodal_section_copy_on_write(section, true, true, false, false);
          _fvmc_nodal_order_indexed_connect(section->_face_index,
                                           section->_face_num,
                                           order,
                                           section->n_elements);
        }

        BFTC_FREE(order);

      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Locally order faces and associated connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent faces (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_order_faces(fvmc_nodal_t       *this_nodal,
                      const fvmc_gnum_t   parent_global_number[])
{
  fvmc_lnum_t  i;
  fvmc_lnum_t  *order = NULL;
  fvmc_nodal_section_t  *section = NULL;

  if (this_nodal == NULL)
    return;

  /* Order locally if necessary */

  for (i = 0 ; i < this_nodal->n_sections ; i++) {

    section = this_nodal->sections[i];

    if (section->entity_dim == 2) {

      assert(section->global_element_num == NULL);

      if (fvmc_order_local_test(section->parent_element_num,
                               parent_global_number,
                               section->n_elements) == false) {

        order = fvmc_order_local(section->parent_element_num,
                                parent_global_number,
                                section->n_elements);

        _fvmc_nodal_order_parent_list(&(section->_parent_element_num),
                                     &(section->parent_element_num),
                                     order,
                                     section->n_elements);

        if (section->type != FVMC_FACE_POLY) {
          fvmc_nodal_section_copy_on_write(section, false, false, false, true);
          _fvmc_nodal_order_strided_connect(section->_vertex_num,
                                           order,
                                           (size_t)(section->stride),
                                           section->n_elements);
        }
        else {
          fvmc_nodal_section_copy_on_write(section, false, false, true, true);
          _fvmc_nodal_order_indexed_connect(section->_vertex_index,
                                           section->_vertex_num,
                                           order,
                                           section->n_elements);
        }

        BFTC_FREE(order);

      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Locally order vertices and update connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent vertices (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_order_vertices(fvmc_nodal_t       *this_nodal,
                         const fvmc_gnum_t   parent_global_number[])
{
  int         i;
  size_t      j;

  fvmc_lnum_t  *order = NULL;
  fvmc_lnum_t  *renumber = NULL;
  fvmc_nodal_section_t  *section = NULL;

  /* Do nothing for trivial cases */

  if (this_nodal == NULL)
    return;

  else if (this_nodal->n_vertices < 2)
    return;

  /* Return if already ordered */

  if (fvmc_order_local_test(this_nodal->parent_vertex_num,
                           parent_global_number,
                           this_nodal->n_vertices) == true)
    return;

  /* Else, we must re-order vertices and update connectivity */

  order = fvmc_order_local(this_nodal->parent_vertex_num,
                          parent_global_number,
                          this_nodal->n_vertices);

  /* Re-order parent list */

  _fvmc_nodal_order_parent_list(&(this_nodal->_parent_vertex_num),
                               &(this_nodal->parent_vertex_num),
                               order,
                               this_nodal->n_vertices);

  /* Calculate renumbering table for associated connectivity
     and free ordering array, no longer needed after that */

  renumber = fvmc_order_local_renumbering(order,
                                         this_nodal->n_vertices);

  BFTC_FREE(order);

  /* Update element connectivities */

  for (i = 0 ; i < this_nodal->n_sections ; i++) {

    section = this_nodal->sections[i];

    fvmc_nodal_section_copy_on_write(section, false, false, false, true);

    for (j = 0 ; j < section->connectivity_size ; j++)
      section->_vertex_num[j] = renumber[section->_vertex_num[j] - 1] + 1;

  }

  /* Free renumbering table */

  BFTC_FREE(renumber);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
