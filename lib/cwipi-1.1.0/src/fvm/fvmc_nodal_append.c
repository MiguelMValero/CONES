/*============================================================================
 * Append sections to a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2008  EDF

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
#include <math.h>

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
#include "fvmc_io_num.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal_append.h"

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
 * Create new section, transferring ownership of the given connectivity
 * and optional parent number arrays to that section.
 *
 * parameters:
 *   n_elements         <-- number of elements in section
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

static fvmc_nodal_section_t *
_transfer_to_section(fvmc_lnum_t      n_elements,
                     fvmc_element_t   t_elt,
                     int              order,
                     int             *ho_uvw_to_local_ordering, 
                     int             *ho_local_to_uvw, 
                     fvmc_lnum_t      face_index[],
                     fvmc_lnum_t      face_num[],
                     fvmc_lnum_t      vertex_index[],
                     fvmc_lnum_t      vertex_num[],
                     fvmc_lnum_t      parent_element_num[])
{
  fvmc_nodal_section_t  *this_section = NULL;

  this_section = fvmc_nodal_section_create(t_elt, order);

  this_section->n_elements = n_elements;


  /* Connectivity */

  if (t_elt == FVMC_CELL_POLY) {
    this_section->_face_index = face_index;
    this_section->_face_num = face_num;
  }

  if (t_elt == FVMC_FACE_POLY || t_elt == FVMC_CELL_POLY)
    this_section->_vertex_index = vertex_index;

  this_section->_vertex_num = vertex_num;

  this_section->_parent_element_num = parent_element_num;

  /* Shared arrays */

  this_section->face_index = this_section->_face_index;
  this_section->face_num = this_section->_face_num;
  this_section->vertex_index = this_section->_vertex_index;
  this_section->vertex_num = this_section->_vertex_num;
  this_section->parent_element_num = this_section->_parent_element_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == FVMC_FACE_POLY) {
    if (this_section->vertex_index != NULL)
      this_section->connectivity_size
        = this_section->vertex_index[this_section->n_elements];
    else
      this_section->connectivity_size = 0;
  }
  else if (this_section->type == FVMC_CELL_POLY) {
    fvmc_lnum_t i, _face_num;
    if (this_section->face_index != NULL) {
      for (i = 0;
           i < this_section->face_index[this_section->n_elements];
           i++) {
        _face_num = FVMC_ABS(this_section->face_num[i]);
        if (_face_num > this_section->n_faces)
          this_section->n_faces = _face_num;
      }
      this_section->connectivity_size
        = this_section->vertex_index[this_section->n_faces];
    }
    else
      this_section->connectivity_size = 0;      
  }

  if (ho_uvw_to_local_ordering != NULL) {

    const int n_nodes = fvmc_nodal_n_vertices_element (t_elt, order);

    this_section->_ho_vertex_num = malloc (sizeof(int) * n_nodes * this_section->n_elements);

    this_section->ho_local_to_user_ordering = malloc (sizeof(int) * n_nodes);
      
    int stride = this_section->entity_dim;

    for (int k = 0; k < n_nodes; k++) {
      const int *_uvw = ho_local_to_uvw + k * stride;
      int idx = 0; 
      for (int l = 0; l < stride; l++) {
        idx += (int) pow((order+1),l) * _uvw[l];
      }
      int local_num = ho_uvw_to_local_ordering[idx];
      this_section->ho_local_to_user_ordering[local_num] = k;
    }

    for (int j = 0; j < this_section->n_elements; j++) {
      int *_ho_vertex_num = this_section->_ho_vertex_num + j * n_nodes;
      const int *_vertex_num = this_section->vertex_num + j * n_nodes;
      for (int k = 0; k < n_nodes; k++) {
        _ho_vertex_num[k] = _vertex_num[this_section->ho_local_to_user_ordering[k]];
      }
    }
  }
  
  return this_section;
}

/*----------------------------------------------------------------------------
 * Create new section, mapping the given connectivity and optional
 * parent number arrays to that section.
 *
 * parameters:
 *   n_elements         <-- number of elements in section
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

static fvmc_nodal_section_t *
_map_to_section(fvmc_lnum_t      n_elements,
                fvmc_element_t   type,
                int              order, 
                int             *ho_uvw_to_local_ordering, 
                int             *ho_local_to_uvw, 
                fvmc_lnum_t      face_index[],
                fvmc_lnum_t      face_num[],
                fvmc_lnum_t      vertex_index[],
                fvmc_lnum_t      vertex_num[],
                fvmc_lnum_t      parent_element_num[])
{
  fvmc_nodal_section_t  *this_section = NULL;

  this_section = fvmc_nodal_section_create(type, order);

  this_section->n_elements = n_elements;

  /* Connectivity */

  if (type == FVMC_CELL_POLY) {
    this_section->face_index = face_index;
    this_section->face_num = face_num;
  }

  if (type == FVMC_FACE_POLY || type == FVMC_CELL_POLY)
    this_section->vertex_index = vertex_index;

  this_section->vertex_num = vertex_num;

  this_section->parent_element_num = parent_element_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == FVMC_FACE_POLY) {
    if (this_section->vertex_index != NULL)
      this_section->connectivity_size
        = this_section->vertex_index[this_section->n_elements];
    else
      this_section->connectivity_size = 0;
  }
  else if (this_section->type == FVMC_CELL_POLY) {
    fvmc_lnum_t i, _face_num; 
    if (this_section->face_index != NULL) {
      for (i = 0;
           i < this_section->face_index[this_section->n_elements];
           i++) {
        _face_num = FVMC_ABS(this_section->face_num[i]);
        if (_face_num > this_section->n_faces)
          this_section->n_faces = _face_num;
      }
      this_section->connectivity_size
        = this_section->vertex_index[this_section->n_faces];
    }
    else
      this_section->connectivity_size = 0;
  }

  if (ho_uvw_to_local_ordering != NULL) {

    const int n_nodes = fvmc_nodal_n_vertices_element (type, order);

    this_section->_ho_vertex_num = malloc (sizeof(int) * n_nodes * this_section->n_elements);

    this_section->ho_local_to_user_ordering = malloc (sizeof(int) * n_nodes);
      
    int stride = this_section->entity_dim;

    for (int k = 0; k < n_nodes; k++) {
      const int *_uvw = ho_local_to_uvw + k * stride;
      int idx = 0; 
      for (int l = 0; l < stride; l++) {
        idx += (int) pow((order+1),l) * _uvw[l];
      }
      int local_num = ho_uvw_to_local_ordering[idx];
      this_section->ho_local_to_user_ordering[local_num] = k;
    }

    for (int j = 0; j < this_section->n_elements; j++) {
      int *_ho_vertex_num = this_section->_ho_vertex_num + j * n_nodes;
      const int *_vertex_num = this_section->vertex_num + j * n_nodes;
      for (int k = 0; k < n_nodes; k++) {
        _ho_vertex_num[k] = _vertex_num[this_section->ho_local_to_user_ordering[k]];
      }
    }
  }

  return this_section;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Append a new section to an existing fvmc_nodal mesh, and transfer
 * ownership of the given connectivity and optional parent number arrays to
 * that section.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_append_by_transfer(fvmc_nodal_t    *this_nodal,
                             fvmc_lnum_t      n_elements,
                             fvmc_element_t   type,
                             fvmc_lnum_t      face_index[],
                             fvmc_lnum_t      face_num[],
                             fvmc_lnum_t      vertex_index[],
                             fvmc_lnum_t      vertex_num[],
                             fvmc_lnum_t      parent_element_num[])
{
  fvmc_nodal_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(this_nodal != NULL);

  n_sections = this_nodal->n_sections;

  /* Create new section */

  BFTC_REALLOC(this_nodal->sections, n_sections + 1, fvmc_nodal_section_t *);
  BFTC_REALLOC(this_nodal->sections_idx, n_sections + 2, int);
  if (n_sections == 0) {
    this_nodal->sections_idx[0] = 0;
  }
  this_nodal->sections_idx[n_sections+1] = this_nodal->sections_idx[n_sections] + n_elements; 

  int  *ho_uvw_to_local_ordering = NULL;
  
  int  *ho_user_to_uvw = NULL;
  
  if (this_nodal->ho_uvw_to_local_ordering != NULL) {
    ho_uvw_to_local_ordering = this_nodal->ho_uvw_to_local_ordering[type];
    ho_user_to_uvw = this_nodal->ho_user_to_uvw[type];
  }

  new_section = _transfer_to_section(n_elements,
                                     type,
                                     this_nodal->order,
                                     ho_uvw_to_local_ordering,
                                     ho_user_to_uvw,
                                     face_index,
                                     face_num,
                                     vertex_index,
                                     vertex_num,
                                     parent_element_num);

  this_nodal->sections[n_sections] = new_section;
  this_nodal->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    this_nodal->n_cells += n_elements;
    break;
  case 2:
    this_nodal->n_faces += n_elements;
    break;
  case 1:
    this_nodal->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------
 * Append a new section to an existing fvmc_nodal mesh, sharing the given
 * given connectivity and optional parent number arrays with the caller.
 *
 * The caller should not destroy or modify the arrays passed to this
 * function until the nodal mesh is destroyed.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_append_shared(fvmc_nodal_t    *this_nodal,
                        fvmc_lnum_t      n_elements,
                        fvmc_element_t   type,
                        fvmc_lnum_t      face_index[],
                        fvmc_lnum_t      face_num[],
                        fvmc_lnum_t      vertex_index[],
                        fvmc_lnum_t      vertex_num[],
                        fvmc_lnum_t      parent_element_num[])
{
  fvmc_nodal_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(this_nodal != NULL);

  n_sections = this_nodal->n_sections;

  /* Create new section */

  BFTC_REALLOC(this_nodal->sections, n_sections + 1, fvmc_nodal_section_t *);
  BFTC_REALLOC(this_nodal->sections_idx, n_sections + 2, int);
  if (n_sections == 0) {
    this_nodal->sections_idx[0] = 0;
  }
  this_nodal->sections_idx[n_sections+1] = this_nodal->sections_idx[n_sections] + n_elements; 

  int  *ho_uvw_to_local_ordering = NULL;
  
  int  *ho_user_to_uvw = NULL;
  
  if (this_nodal->ho_uvw_to_local_ordering != NULL) {
    ho_uvw_to_local_ordering = this_nodal->ho_uvw_to_local_ordering[type];
    ho_user_to_uvw = this_nodal->ho_user_to_uvw[type];
  }
  
  new_section = _map_to_section(n_elements,
                                type,
                                this_nodal->order,
                                ho_uvw_to_local_ordering,
                                ho_user_to_uvw, 
                                face_index,
                                face_num,
                                vertex_index,
                                vertex_num,
                                parent_element_num);

  this_nodal->sections[n_sections] = new_section;
  this_nodal->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    this_nodal->n_cells += n_elements;
    break;
  case 2:
    this_nodal->n_faces += n_elements;
    break;
  case 1:
    this_nodal->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
