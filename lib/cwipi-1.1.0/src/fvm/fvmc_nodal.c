/*============================================================================
 * Main structure for a nodal representation associated with a mesh
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
#include "fvmc_parall.h"
#include "fvmc_tesselation.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"

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

/* Number of vertices associated with each "nodal" element type */

//const int  fvmc_nodal_n_vertices_element[] = {2,   /* Edge */
//                                             3,   /* Triangle */
//                                             4,   /* Quadrangle */
//                                             0,   /* Simple polygon */
//                                             4,   /* Tetrahedron */
//                                             5,   /* Pyramid */
//                                             6,   /* Prism */
//                                             8,   /* Hexahedron */
//                                             0};  /* Simple polyhedron */

/* Number of vertices associated with each "nodal" element type */

static int  fvmc_nodal_n_edges_element[] = {1,   /* Edge */
                                           3,   /* Triangle */
                                           4,   /* Quadrangle */
                                           0,   /* Simple polygon */
                                           6,   /* Tetrahedron */
                                           8,   /* Pyramid */
                                           9,   /* Prism */
                                           12,  /* Hexahedron */
                                           0};  /* Simple polyhedron */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * compute the idirection (u, v, w) to local ordering
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   t_elt                <-- element type
 *
 * returns:
 *   indirection
 *----------------------------------------------------------------------------*/

static int* _uvw_to_local_ordering (const fvmc_nodal_t *this_nodal,
                                    const fvmc_element_t t_elt)
{
  int s_ordering = 0;

  int order = this_nodal->order;
  
  int n_vtx = fvmc_nodal_n_vertices_element (t_elt, this_nodal->order);
  
  switch(t_elt) {
  case FVMC_EDGE:               /* Edge */
    s_ordering = order + 1;
    break;
  case FVMC_FACE_TRIA:          /* Triangle */
  case FVMC_FACE_QUAD:          /* Quadrangle */
    s_ordering = (order + 1) * (order + 1);
    break;
  case FVMC_CELL_TETRA:         /* Tetrahedron */
  case FVMC_CELL_PYRAM:         /* Pyramid */
  case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
  case FVMC_CELL_HEXA:         /* Hexahedron (brick) */
    s_ordering = (order + 1) * (order + 1) * (order + 1);
    break;
  default:  
    bftc_error(__FILE__, __LINE__, 0,
                  _("_uvw_to_local_ordering : high order unavailable "));
  }

  int *_uvw_to_local = (int *) malloc (sizeof(int) * s_ordering); 

  for (int i = 0; i < s_ordering; i++) {
    _uvw_to_local[i] = -1;
  }

  switch(t_elt) {
  case FVMC_EDGE: {
    int _cpt_vtx = 0;
    for (int u = 0; u < order + 1; u++) {
      _uvw_to_local[u] = _cpt_vtx++;
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_FACE_TRIA: {
    int _cpt_vtx = 0;
    for (int v = 0; v < order + 1; v++) {
      for (int u = 0; u < order + 1 - v; u++) {
        int idx = (order + 1) * v + u;

        _uvw_to_local[idx] = _cpt_vtx++;
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_FACE_QUAD: {
    int _cpt_vtx = 0;
    for (int v = 0; v < order + 1; v++) {
      for (int u = 0; u < order + 1; u++) {
        int idx = (order + 1) * v + u;

        _uvw_to_local[idx] = _cpt_vtx++;
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_TETRA: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1 - w; v++) {
        for (int u = 0; u < order + 1 - v - w; u++) {
          int idx =   (order + 1) * (order + 1) * w 
                    + (order + 1) * v
                    + u;

          _uvw_to_local[idx] = _cpt_vtx++;
        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_PYRAM: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1 - w; v++) {
        for (int u = 0; u < order + 1 - w; u++) {
          int idx =   (order + 1) * (order + 1) * w 
                    + (order + 1) * v
                    + u;

          _uvw_to_local[idx] = _cpt_vtx++;
        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_PRISM: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1; v++) {
        for (int u = 0; u < order + 1 - v ; u++) {
          int idx =   (order + 1) * (order + 1) * w 
                    + (order + 1) * v
                    + u;

          _uvw_to_local[idx] = _cpt_vtx++;
        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_HEXA: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1; v++) {
        for (int u = 0; u < order + 1; u++) {
          int idx =   (order + 1) * (order + 1) * w 
                    + (order + 1) * v
                    + u;

          _uvw_to_local[idx] = _cpt_vtx++;
        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  default:  
    bftc_error(__FILE__, __LINE__, 0,
                  _("_uvw_to_local_ordering : high order unavailable "));
  }

  return _uvw_to_local;

}


/*----------------------------------------------------------------------------
 * compute the coordinates of the nodes of the reference element
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   t_elt                <-- element type
 *
 * returns:
 *   coordinates
 *----------------------------------------------------------------------------*/

static double* _local_ref_nodes (const fvmc_nodal_t *this_nodal,
                                 const fvmc_element_t t_elt)
{
  int order = this_nodal->order;
  
  int n_vtx = fvmc_nodal_n_vertices_element (t_elt, this_nodal->order);
  
  double *_coords = (double *) malloc (sizeof(double) * 3 * n_vtx); 

  double step = 1. / (double) order;
  
  switch(t_elt) {
  case FVMC_EDGE: {
    int _cpt_vtx = 0;
    for (int u = 0; u < order + 1; u++) {
      _coords[3 * _cpt_vtx    ] = u * step;
      _coords[3 * _cpt_vtx + 1] = 0;
      _coords[3 * _cpt_vtx + 2] = 0;
    }
    break;
  }
  case FVMC_FACE_TRIA: {
    int _cpt_vtx = 0;
    for (int v = 0; v < order + 1; v++) {
      for (int u = 0; u < order + 1 - v; u++) {
        _coords[3 * _cpt_vtx    ] = u * step;
        _coords[3 * _cpt_vtx + 1] = v * step;
        _coords[3 * _cpt_vtx + 2] = 0;
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_FACE_QUAD: {
    int _cpt_vtx = 0;
    for (int v = 0; v < order + 1; v++) {
      for (int u = 0; u < order + 1; u++) {
        _coords[3 * _cpt_vtx    ] = u * step;
        _coords[3 * _cpt_vtx + 1] = v * step;
        _coords[3 * _cpt_vtx + 2] = 0;
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_TETRA: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1 - w; v++) {
        for (int u = 0; u < order + 1 - v - w; u++) {
          _coords[3 * _cpt_vtx    ] = u * step;
          _coords[3 * _cpt_vtx + 1] = v * step;
          _coords[3 * _cpt_vtx + 2] = w * step;

        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_PYRAM: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1 - w; v++) {
        for (int u = 0; u < order + 1 - w; u++) {
          _coords[3 * _cpt_vtx    ] = u * step;
          _coords[3 * _cpt_vtx + 1] = v * step;
          _coords[3 * _cpt_vtx + 2] = w * step;

        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_PRISM: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1; v++) {
        for (int u = 0; u < order + 1 - v ; u++) {
          _coords[3 * _cpt_vtx    ] = u * step;
          _coords[3 * _cpt_vtx + 1] = v * step;
          _coords[3 * _cpt_vtx + 2] = w * step;

        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  case FVMC_CELL_HEXA: {
    int _cpt_vtx = 0;
    for (int w = 0; w < order + 1; w++) {
      for (int v = 0; v < order + 1; v++) {
        for (int u = 0; u < order + 1; u++) {
          _coords[3 * _cpt_vtx    ] = u * step;
          _coords[3 * _cpt_vtx + 1] = v * step;
          _coords[3 * _cpt_vtx + 2] = w * step;

        }
      }
    }
    assert (_cpt_vtx == n_vtx);
    break;
  }
  default:  
    bftc_error(__FILE__, __LINE__, 0,
                  _("_local_ref_nodes : high order unavailable "));
  }

  return _coords;

}


/*----------------------------------------------------------------------------
 * Compare edges (qsort function).
 *
 * parameters:
 *   x <-> pointer to first edge definition
 *   y <-> pointer to second edge definition
 *
 * returns:
 *   result of strcmp() on group names
 *----------------------------------------------------------------------------*/

static int _compare_edges(const void *x, const void *y)
{
  int retval = 1;

  const fvmc_lnum_t *e0 = (const fvmc_lnum_t *) x;
  const fvmc_lnum_t *e1 = (const fvmc_lnum_t *) y;

  if (e0[0] < e1[0])
    retval = -1;

  else if (e0[0] == e1[0]) {
    if (e0[1] < e1[1])
      retval = -1;
    else if (e0[1] == e1[1])
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Copy a nodal mesh section representation structure, sharing arrays
 * with the original structure.
 *
 * parameters:
 *   this_section <-> pointer to structure that should be copied
 *
 * returns:
 *   pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static fvmc_nodal_section_t *
_fvmc_nodal_section_copy(const fvmc_nodal_section_t *this_section)
{
  fvmc_nodal_section_t  *new_section;

  BFTC_MALLOC(new_section, 1, fvmc_nodal_section_t);

  /* Global information */

  new_section->entity_dim = this_section->entity_dim;

  new_section->n_elements = this_section->n_elements;
  new_section->type = this_section->type;

  /* Connectivity */

  new_section->connectivity_size = this_section->connectivity_size;
  new_section->stride = this_section->stride;

  new_section->n_faces = this_section->n_faces;

  new_section->face_index = this_section->face_index;
  new_section->face_num = this_section->face_num;
  new_section->vertex_index = this_section->vertex_index;
  new_section->vertex_num = this_section->vertex_num;

  new_section->_face_index = NULL;
  new_section->_face_num   = NULL;
  new_section->_vertex_index = NULL;
  new_section->_vertex_num = NULL;

  new_section->tesselation = NULL;  /* TODO: copy tesselation */

  /* Numbering */
  /*-----------*/

  new_section->parent_element_num = this_section->parent_element_num;
  new_section->_parent_element_num = NULL;

  if (this_section->global_element_num != NULL) {
    fvmc_lnum_t n_ent
      = fvmc_io_num_get_local_count(this_section->global_element_num);
    fvmc_gnum_t global_count
      = fvmc_io_num_get_global_count(this_section->global_element_num);
    const fvmc_gnum_t *global_num
      = fvmc_io_num_get_global_num(this_section->global_element_num);

    new_section->global_element_num
      = fvmc_io_num_create_shared(global_num, global_count, n_ent);
  }
  else
    new_section->global_element_num = NULL;

  return (new_section);
}

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh representation section: only the associations
 * (numberings) necessary to redistribution of fields for output are
 * conserved, the full connectivity being no longer useful once it has been
 * output.
 *
 * parameters:
 *   this_section      <-> pointer to structure that should be reduced
 *
 * returns:
 *   true if connectivity has been reduced
 *----------------------------------------------------------------------------*/

static _Bool
_fvmc_nodal_section_reduce(fvmc_nodal_section_t  * this_section)
{
  _Bool retval = false;

  /* If we have a tesselation of polyhedra (face index != NULL),
     we may need to keep the connectivity information, to
     interpolate nodal values to added vertices */

  if (   this_section->tesselation == NULL
      || this_section->_face_index == NULL) {

      /* Connectivity */

    this_section->connectivity_size = 0;

    if (this_section->_face_index != NULL)
      BFTC_FREE(this_section->_face_index);
    this_section->face_index = NULL;

    if (this_section->_face_num != NULL)
      BFTC_FREE(this_section->_face_num);
    this_section->face_num = NULL;

    if (this_section->_vertex_index != NULL)
      BFTC_FREE(this_section->_vertex_index);
    this_section->vertex_index = NULL;

    if (this_section->_vertex_num != NULL)
      BFTC_FREE(this_section->_vertex_num);
    this_section->vertex_num = NULL;

    retval = true;
  }

  if (this_section->tesselation != NULL)
    fvmc_tesselation_reduce(this_section->tesselation);

  return retval;
}

/*----------------------------------------------------------------------------
 * Change entity parent numbering; this is useful when entities of the
 * parent mesh have been renumbered after a nodal mesh representation
 * structure's creation. As the parent_num[] array is defined only when
 * non trivial (i.e. not 1, 2, ..., n), it may be allocated or freed
 * by this function. The return argument corresponds to the new
 * pointer which should replace the parent_num input argument.
 *
 * parameters:
 *   parent_num_size     <-- size of local parent numbering array
 *   new_parent_num      <-- pointer to local parent renumbering array
 *                           ({1, ..., n} <-- {1, ..., n})
 *   parent_num          <-> pointer to local parent numbering array
 *   _parent_num         <-> pointer to local parent numbering array if
 *                           owner, NULL otherwise
 *
 * returns:
 *   pointer to resulting parent_num[] array
 *----------------------------------------------------------------------------*/

static fvmc_lnum_t *
_renumber_parent_num(fvmc_lnum_t         parent_num_size,
                     const fvmc_lnum_t   new_parent_num[],
                     const fvmc_lnum_t   parent_num[],
                     fvmc_lnum_t         _parent_num[])
{
  int  i;
  fvmc_lnum_t  old_num_id;
  fvmc_lnum_t *parent_num_p = _parent_num;
  _Bool trivial = true;

  if (parent_num_size > 0 && new_parent_num != NULL) {

    if (parent_num_p != NULL) {
      for (i = 0; i < parent_num_size; i++) {
        old_num_id = parent_num_p[i] - 1;
        parent_num_p[i] = new_parent_num[old_num_id];
        if (parent_num_p[i] != i+1)
          trivial = false;
      }
    }
    else {
      BFTC_MALLOC(parent_num_p, parent_num_size, fvmc_lnum_t);
      if (parent_num != NULL) {
        for (i = 0; i < parent_num_size; i++) {
          old_num_id = parent_num[i] - 1;
          parent_num_p[i] = new_parent_num[old_num_id];
          if (parent_num_p[i] != i+1)
            trivial = false;
        }
      }
      else {
        for (i = 0; i < parent_num_size; i++) {
          parent_num_p[i] = new_parent_num[i];
          if (parent_num_p[i] != i+1)
            trivial = false;
        }
      }
    }
  }

  if (trivial == true)
    BFTC_FREE(parent_num_p);

  return parent_num_p;
}

/*----------------------------------------------------------------------------
 * Renumber vertices based on those actually referenced, and update
 * connectivity arrays and parent numbering in accordance.
 *
 * The number of vertices assigned to the nodal mesh (this_nodal->n_vertices)
 * is computed and set by this function. If this number was previously
 * non-zero (i.e. vertices have already been assigned to the structure),
 * those vertices are considered as referenced. This is useful if we want
 * to avoid discarding a given set of vertices, such as when building a
 * nodal mesh representation containing only vertices.
 *
 * parameters:
 *   this_nodal <-> nodal mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_vertices(fvmc_nodal_t  *this_nodal)
{
  size_t      i;
  int         section_id;
  fvmc_lnum_t  j;
  fvmc_lnum_t  vertex_id;
  fvmc_lnum_t  n_vertices;
  fvmc_nodal_section_t  *section;

  fvmc_lnum_t  *loc_vertex_num = NULL;
  fvmc_lnum_t   max_vertex_num = 0;

  /* Find maximum vertex reference */
  /*-------------------------------*/

  /* The mesh may already contain direct vertex references
     (as in the case of a "mesh" only containing vertices) */

  if (this_nodal->n_vertices > 0) {
    if (this_nodal->parent_vertex_num != NULL) {
      for (j = 0; j < this_nodal->n_vertices; j++) {
        if (this_nodal->parent_vertex_num[j] > max_vertex_num)
          max_vertex_num = this_nodal->parent_vertex_num[j];
      }
    }
    else
      max_vertex_num = this_nodal->n_vertices;
  }

  /* In most cases, the mesh will reference vertices through elements */

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {
    section = this_nodal->sections[section_id];
    if (this_nodal->parent_vertex_num != NULL) {
      for (i = 0; i < section->connectivity_size; i++) {
        fvmc_lnum_t vertex_num
          = this_nodal->parent_vertex_num[section->vertex_num[i] - 1];
        if (vertex_num > max_vertex_num)
          max_vertex_num = vertex_num;
      }
    }
    else {
      for (i = 0; i < section->connectivity_size; i++) {
        if (section->vertex_num[i] > max_vertex_num)
          max_vertex_num = section->vertex_num[i];
      }
    }
  }

  /* Flag referenced vertices and compute size */
  /*-------------------------------------------*/

  BFTC_MALLOC(loc_vertex_num, max_vertex_num, fvmc_lnum_t);

  for (vertex_id = 0; vertex_id < max_vertex_num; vertex_id++)
    loc_vertex_num[vertex_id] = 0;

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {
    section = this_nodal->sections[section_id];
    if (this_nodal->parent_vertex_num != NULL) {
      for (i = 0; i < section->connectivity_size; i++) {
        vertex_id
          = this_nodal->parent_vertex_num[section->vertex_num[i] - 1] - 1;
        if (loc_vertex_num[vertex_id] == 0)
          loc_vertex_num[vertex_id] = 1;
      }
    }
    else {
      for (i = 0; i < section->connectivity_size; i++) {
        vertex_id = section->vertex_num[i] - 1;
        if (loc_vertex_num[vertex_id] == 0)
          loc_vertex_num[vertex_id] = 1;
      }
    }
  }

  /* Build vertices renumbering */
  /*----------------------------*/

  n_vertices = 0;

  for (vertex_id = 0; vertex_id < max_vertex_num; vertex_id++) {
    if (loc_vertex_num[vertex_id] == 1) {
      n_vertices += 1;
      loc_vertex_num[vertex_id] = n_vertices;
    }
  }
  this_nodal->n_vertices = n_vertices;

  /* Update connectivity and vertex parent numbering */
  /*-------------------------------------------------*/

  /* If all vertices are flagged, no need to renumber */

  if (n_vertices == max_vertex_num)
    BFTC_FREE(loc_vertex_num);

  else {

    /* Update connectivity */

    for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {
      section = this_nodal->sections[section_id];
      if (section->_vertex_num == NULL)
        fvmc_nodal_section_copy_on_write(section, false, false, false, true);
      if (this_nodal->parent_vertex_num != NULL) {
        for (i = 0; i < section->connectivity_size; i++) {
          vertex_id
            = this_nodal->parent_vertex_num[section->vertex_num[i] - 1] - 1;
          section->_vertex_num[i] = loc_vertex_num[vertex_id];
        }
      }
      else {
        for (i = 0; i < section->connectivity_size; i++) {
          vertex_id = section->vertex_num[i] - 1;
          section->_vertex_num[i] = loc_vertex_num[vertex_id];
        }
      }
    }

    /* Build or update vertex parent numbering */

    this_nodal->parent_vertex_num = NULL;
    if (this_nodal->_parent_vertex_num != NULL)
      BFTC_FREE(this_nodal->_parent_vertex_num);

    if (loc_vertex_num != NULL) {
      BFTC_MALLOC(this_nodal->_parent_vertex_num, n_vertices, fvmc_lnum_t);
      for (vertex_id = 0; vertex_id < max_vertex_num; vertex_id++) {
        if (loc_vertex_num[vertex_id] > 0)
          this_nodal->_parent_vertex_num[loc_vertex_num[vertex_id] - 1]
            = vertex_id + 1;
      }
      this_nodal->parent_vertex_num = this_nodal->_parent_vertex_num;
    }
  }

  /* Free renumbering array */

  BFTC_FREE(loc_vertex_num);
}

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure section.
 *
 * parameters:
 *   this_section <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

static void
_fvmc_nodal_section_dump(const fvmc_nodal_section_t  *this_section)
{
  fvmc_lnum_t  n_elements, i, j;
  const fvmc_lnum_t  *idx, *num;

  /* Global indicators */
  /*--------------------*/

  bftc_printf("\n"
             "Entity dimension:     %d\n"
             "Number of elements:   %ld\n"
             "Element type:         %s\n",
             this_section->entity_dim, (long)this_section->n_elements,
             fvmc_elements_type_name[this_section->type]);

  bftc_printf("\n"
             "Connectivity_size:     %lu\n"
             "Stride:                %d\n"
             "Number of faces:       %d\n",
             (unsigned long)(this_section->connectivity_size),
             this_section->stride,
             (long)(this_section->n_faces));

  bftc_printf("\n"
             "Pointers to shareable arrays:\n"
             "  face_index:           %p\n"
             "  face_num:             %p\n"
             "  vertex_index:         %p\n"
             "  vertex_num:           %p\n"
             "  parent_element_num:   %p\n",
             this_section->face_index, this_section->face_num,
             this_section->vertex_index, this_section->vertex_num,
             this_section->parent_element_num);

  bftc_printf("\n"
             "Pointers to local arrays:\n"
             "  _face_index:          %p\n"
             "  _face_num:            %p\n"
             "  _vertex_index:        %p\n"
             "  _vertex_num:          %p\n"
             "  _parent_element_num:  %p\n",
             this_section->_face_index, this_section->_face_num,
             this_section->_vertex_index, this_section->_vertex_num,
             this_section->_parent_element_num);

  if (this_section->face_index != NULL) {
    bftc_printf("\nPolyhedra -> faces (polygons) connectivity:\n\n");
    n_elements = this_section->n_elements;
    idx = this_section->face_index;
    num = this_section->face_num;
    for (i = 0; i < n_elements; i++) {
      bftc_printf("%10d (idx = %10d) %10d\n",
                 i+1, idx[i], num[idx[i]]);
      for (j = idx[i] + 1; j < idx[i + 1]; j++)
        bftc_printf("                              %10d\n", num[j]);
    }
    bftc_printf("      end  (idx = %10d)\n", idx[n_elements]);
  }

  if (this_section->vertex_index != NULL) {
    fvmc_lnum_t  n_faces = (this_section->n_faces > 0) ?
                          this_section->n_faces : this_section->n_elements;
    bftc_printf("\nPolygons -> vertices connectivity:\n\n");
    idx = this_section->vertex_index;
    num = this_section->vertex_num;
    for (i = 0; i < n_faces; i++) {
      bftc_printf("%10d (idx = %10d) %10d\n",
                i + 1, idx[i], num[idx[i]]);
      for (j = idx[i] + 1; j < idx[i + 1]; j++)
        bftc_printf("                              %10d\n", num[j]);
    }
    bftc_printf("      end  (idx = %10d)\n", idx[n_faces]);
  }

  else {
    bftc_printf("\nElements -> vertices connectivity:\n\n");
    n_elements = this_section->n_elements;
    num = this_section->vertex_num;
    switch(this_section->stride) {
    case 2:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d\n",
                   i+1, num[i*2], num[i*2+1]);
      break;
    case 3:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d %10d\n",
                   i+1, num[i*3], num[i*3+1], num[i*3+2]);
      break;
    case 4:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d %10d %10d\n",
                   i+1, num[i*4], num[i*4+1], num[i*4+2],
                   num[i*4+3]);
      break;
    case 5:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d %10d %10d %10d\n",
                   i+1, num[i*5], num[i*5+1], num[i*5+2],
                   num[i*5+3], num[i*5+4]);
      break;
    case 6:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d %10d %10d %10d %10d\n",
                   i+1, num[i*6], num[i*6+1], num[i*6+2],
                   num[i*6+3], num[i*6+4], num[i*6+5]);
      break;
    case 8:
      for (i = 0; i < n_elements; i++)
        bftc_printf("%10d : %10d %10d %10d %10d %10d %10d %10d %10d\n",
                   i+1, num[i*8], num[i*8+1], num[i*8+2], num[i*8+3],
                   num[i*8+4], num[i*8+5], num[i*8+6], num[i*8+7]);
      break;
    default:
      for (i = 0; i < n_elements; i++) {
        bftc_printf("%10d :", i+1);
        for (j = 0; j < this_section->stride; j++)
          bftc_printf(" %10d", num[i*this_section->stride + j]);
        bftc_printf("\n");
      }
    }
  }

  /* Faces tesselation */

  if (this_section->tesselation != NULL)
    fvmc_tesselation_dump(this_section->tesselation);

  /* Numbers of associated elements in the parent mesh */

  bftc_printf("\nLocal element numbers in parent mesh:\n");
  if (this_section->parent_element_num == NULL)
    bftc_printf("\n  Nil\n\n");
  else {
    for (i = 0; i < this_section->n_elements; i++)
      bftc_printf("  %10d %10d\n", i+1, this_section->parent_element_num[i]);
  }

  /* Global element numbers (only for parallel execution) */

  if (this_section->global_element_num != NULL) {
    bftc_printf("\nGlobal element numbers:\n");
    fvmc_io_num_dump(this_section->global_element_num);
  }
}

/*============================================================================
 * Semi-private function definitions (prototypes in fvmc_nodal_priv.h)
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Get number of vertices.
 *
 * returns:
 *   Number of vertices
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_n_vertices_element (fvmc_element_t type, int order)
{
 int n_vtx = 0;
 int _order = order;
 if (order == -1) {
   _order = 1;
 }
 
 switch(type) {
 case FVMC_EDGE:               /* Edge */
   n_vtx = (_order+1);
   break;
 case FVMC_FACE_TRIA:          /* Triangle */
   n_vtx = (_order+1)*(_order+2)/2; 
   break;
 case FVMC_FACE_QUAD:          /* Quadrangle */
   n_vtx = (_order+1)*(_order+1); 
   break;
 case FVMC_FACE_POLY:          /* Simple Polygon */
   n_vtx = -1;
   break;
 case FVMC_CELL_TETRA:         /* Tetrahedron */
   n_vtx = (_order+1)*(_order+2)*(_order+3)/6; 
   break;
 case FVMC_CELL_PYRAM:         /* Pyramid */
   n_vtx = (_order+1)*(_order+2)*(2*_order+3)/6;
   break;
 case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
   n_vtx = (_order+1)*(_order+1)*(_order+2)/2; 
   break;
 case FVMC_CELL_HEXA:         /* Hexahedron (brick) */
   n_vtx = (_order+1)*(_order+1)*(_order+1); 
   break;
 case FVMC_CELL_POLY:          /* Simple Polyhedron (convex or quasi-convex) */
   n_vtx = -1;
   break;
 default:  
   n_vtx = -1;
 }

 return n_vtx;
}

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh section representation structure.
 *
 * parameters:
 *   type <-- type of element defined by this section
 *
 * returns:
 *   pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_section_t *
fvmc_nodal_section_create(const fvmc_element_t  type, int order)
{
  fvmc_nodal_section_t  *this_section;

  BFTC_MALLOC(this_section, 1, fvmc_nodal_section_t);

  /* Global information */

  if (type == FVMC_EDGE)
    this_section->entity_dim = 1;
  else if (type >= FVMC_FACE_TRIA && type <= FVMC_FACE_POLY)
    this_section->entity_dim = 2;
  else
    this_section->entity_dim = 3;

  this_section->n_elements = 0;
  this_section->type = type;
  this_section->order = order;

  this_section->ho_local_to_user_ordering = NULL;
  this_section->_ho_vertex_num = NULL;
  
  /* Connectivity */

  this_section->connectivity_size = 0;

  if (type != FVMC_FACE_POLY && type != FVMC_CELL_POLY)
    this_section->stride = fvmc_nodal_n_vertices_element(type, order);
  else
    this_section->stride = 0;

  this_section->n_faces = 0;
  this_section->face_index = NULL;
  this_section->face_num   = NULL;
  this_section->vertex_index = NULL;
  this_section->vertex_num = NULL;

  this_section->_face_index = NULL;
  this_section->_face_num   = NULL;
  this_section->_vertex_index = NULL;
  this_section->_vertex_num = NULL;

  this_section->tesselation = NULL;

  /* Numbering */
  /*-----------*/

  this_section->parent_element_num = NULL;
  this_section->_parent_element_num = NULL;

  this_section->global_element_num = NULL;

  return (this_section);

}

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh section representation structure.
 *
 * parameters:
 *   this_section <-> pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_nodal_section_t *
fvmc_nodal_section_destroy(fvmc_nodal_section_t  * this_section)
{

  /* Connectivity */

  if (this_section->_face_index != NULL)
    BFTC_FREE(this_section->_face_index);
  if (this_section->_face_num != NULL)
    BFTC_FREE(this_section->_face_num);

  if (this_section->_vertex_index != NULL)
    BFTC_FREE(this_section->_vertex_index);
  if (this_section->_vertex_num != NULL)
    BFTC_FREE(this_section->_vertex_num);

  if (this_section->tesselation != NULL)
    fvmc_tesselation_destroy(this_section->tesselation);

  /* Numbering */
  /*-----------*/

  if (this_section->parent_element_num != NULL) {
    this_section->parent_element_num = NULL;
    BFTC_FREE(this_section->_parent_element_num);
  }

  if (this_section->global_element_num != NULL)
    fvmc_io_num_destroy(this_section->global_element_num);

  
  if (this_section->ho_local_to_user_ordering != NULL)
    free (this_section->ho_local_to_user_ordering);
  
  if (this_section->_ho_vertex_num != NULL)
    free (this_section->_ho_vertex_num);


  /* Main structure destroyed and NULL returned */

  BFTC_FREE(this_section);

  return (this_section);

}

/*----------------------------------------------------------------------------
 * Copy selected shared connectivity information to private connectivity
 * for a nodal mesh section.
 *
 * parameters:
 *   this_section      <-> pointer to section structure
 *   copy_face_index   <-- copy face index (polyhedra only) ?
 *   copy_face_num     <-- copy face numbers (polyhedra only) ?
 *   copy_vertex_index <-- copy vertex index (polyhedra/polygons only) ?
 *   copy_vertex_num   <-- copy vertex numbers ?
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_section_copy_on_write(fvmc_nodal_section_t  *this_section,
                                _Bool                 copy_face_index,
                                _Bool                 copy_face_num,
                                _Bool                 copy_vertex_index,
                                _Bool                 copy_vertex_num)
{
  fvmc_lnum_t  n_faces;
  size_t  i;

  if (copy_face_index == true
      && this_section->face_index != NULL && this_section->_face_index == NULL) {
    BFTC_MALLOC(this_section->_face_index, this_section->n_elements + 1, fvmc_lnum_t);
    for (i = 0; i < (size_t)(this_section->n_elements + 1); i++) {
      this_section->_face_index[i] = this_section->face_index[i];
    }
    this_section->face_index = this_section->_face_index;
  }

  if (copy_face_num == true
      && this_section->face_num != NULL && this_section->_face_num == NULL) {
    n_faces = this_section->face_index[this_section->n_elements];
    BFTC_MALLOC(this_section->_face_num, n_faces, fvmc_lnum_t);
    for (i = 0; i < (size_t)n_faces; i++) {
      this_section->_face_num[i] = this_section->face_num[i];
    }
    this_section->face_num = this_section->_face_num;
  }

  if (   copy_vertex_index == true
      && this_section->vertex_index != NULL
      && this_section->_vertex_index == NULL) {
    if (this_section->n_faces != 0)
      n_faces = this_section->n_faces;
    else
      n_faces = this_section->n_elements;
    BFTC_MALLOC(this_section->_vertex_index, n_faces + 1, fvmc_lnum_t);
    for (i = 0; i < (size_t)n_faces + 1; i++) {
      this_section->_vertex_index[i] = this_section->vertex_index[i];
    }
    this_section->vertex_index = this_section->_vertex_index;
  }

  if (copy_vertex_num == true && this_section->_vertex_num == NULL) {
    BFTC_MALLOC(this_section->_vertex_num,
               this_section->connectivity_size, fvmc_lnum_t);
    for (i = 0; i < this_section->connectivity_size; i++) {
      this_section->_vertex_num[i] = this_section->vertex_num[i];
    }
    this_section->vertex_num = this_section->_vertex_num;
  }

}

/*----------------------------------------------------------------------------
 * Return global number of elements associated with section.
 *
 * parameters:
 *   this_section      <-- pointer to section structure
 *
 * returns:
 *   global number of elements associated with section
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_nodal_section_n_g_elements(const fvmc_nodal_section_t  *this_section)
{
  if (this_section->global_element_num != NULL)
    return fvmc_io_num_get_global_count(this_section->global_element_num);
  else
    return this_section->n_elements;
}

/*----------------------------------------------------------------------------
 * Return global number of vertices associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   global number of vertices associated with nodal mesh
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_nodal_n_g_vertices(const fvmc_nodal_t  *this_nodal)
{
  fvmc_gnum_t  n_g_vertices;

  if (this_nodal->global_vertex_num != NULL)
    n_g_vertices = fvmc_io_num_get_global_count(this_nodal->global_vertex_num);
  else
    n_g_vertices = this_nodal->n_vertices;

  return n_g_vertices;
}

/*----------------------------------------------------------------------------
 * Define cell->face connectivity for strided cell types.
 *
 * parameters:
 *   element_type     <-- type of strided element
 *   n_faces          --> number of element faces
 *   n_face_vertices  --> number of vertices of each face
 *   face_vertices    --> face -> vertex base connectivity (0 to n-1)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_cell_face_connect(fvmc_element_t   element_type,
                            int            *n_faces,
                            int             n_face_vertices[6],
                            int             face_vertices[6][4])
{
  int i, j;

  /* Initialization */

  *n_faces = 0;

  for (i = 0; i < 6; i++) {
    n_face_vertices[i] = 0;
    for (j = 0; j < 4; j++)
      face_vertices[i][j] = 0;
  }

  /* Define connectivity based on element type */

  switch(element_type) {

  case FVMC_CELL_TETRA:
    {
      fvmc_lnum_t _face_vertices[4][3] = {{1, 3, 2},     /*       x 4     */
                                         {1, 2, 4},     /*      /|\      */
                                         {1, 4, 3},     /*     / | \     */
                                         {2, 3, 4}};    /*    /  |  \    */
                                                        /* 1 x- -|- -x 3 */
      for (i = 0; i < 4; i++) {                         /*    \  |  /    */
        n_face_vertices[i] = 3;                         /*     \ | /     */
        for (j = 0; j < 3; j++)                         /*      \|/      */
          face_vertices[i][j] = _face_vertices[i][j];   /*       x 2     */
      }
      *n_faces = 4;
    }
    break;

  case FVMC_CELL_PYRAM:
    {
      fvmc_lnum_t _n_face_vertices[5] = {3, 3, 3, 3, 4};
      fvmc_lnum_t _face_vertices[5][4] = {{1, 2, 5, 0},  /*        5 x       */
                                         {1, 5, 4, 0},  /*         /|\      */
                                         {2, 3, 5, 0},  /*        //| \     */
                                         {3, 4, 5, 0},  /*       // |  \    */
                                         {1, 4, 3, 2}}; /*    4 x/--|---x 3 */
                                                        /*     //   |  /    */
      for (i = 0; i < 5; i++) {                         /*    //    | /     */
        n_face_vertices[i] = _n_face_vertices[i];       /* 1 x-------x 2    */
        for (j = 0; j < 4; j++)
          face_vertices[i][j] = _face_vertices[i][j];
      }
      *n_faces = 5;
    }
    break;

  case FVMC_CELL_PRISM:
    {
      fvmc_lnum_t _n_face_vertices[5] = {3, 3, 4, 4, 4};
      fvmc_lnum_t _face_vertices[5][4] = {{1, 3, 2, 0},  /* 4 x-------x 6 */
                                         {4, 5, 6, 0},  /*   |\     /|   */
                                         {1, 2, 5, 4},  /*   | \   / |   */
                                         {1, 4, 6, 3},  /* 1 x- \-/ -x 3 */
                                         {2, 3, 6, 5}}; /*    \ 5x  /    */
                                                        /*     \ | /     */
      for (i = 0; i < 5; i++) {                         /*      \|/      */
        n_face_vertices[i] = _n_face_vertices[i];       /*       x 2     */
        for (j = 0; j < 4; j++)
          face_vertices[i][j] = _face_vertices[i][j];
      }
      *n_faces = 5;
    }
    break;

  case FVMC_CELL_HEXA:
    {
      fvmc_lnum_t _n_face_vertices[6] = {4, 4, 4, 4, 4, 4};
      fvmc_lnum_t _face_vertices[6][4] = {{1, 4, 3, 2},  /*    8 x-------x 7 */
                                         {1, 2, 6, 5},  /*     /|      /|   */
                                         {1, 5, 8, 4},  /*    / |     / |   */
                                         {2, 3, 7, 6},  /* 5 x-------x6 |   */
                                         {3, 4, 8, 7},  /*   | 4x----|--x 3 */
                                         {5, 6, 7, 8}}; /*   | /     | /    */
      for (i = 0; i < 6; i++) {                         /*   |/      |/     */
        n_face_vertices[i] = _n_face_vertices[i];       /* 1 x-------x 2    */
        for (j = 0; j < 4; j++)
          face_vertices[i][j] = _face_vertices[i][j];
      }
      *n_faces = 6;
    }
    break;

  default:
    assert(0);
  }

  /* Switch from (1, n) to (0, n-1) numbering */

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 4; j++)
      face_vertices[i][j] -= 1;
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh representation structure.
 *
 * parameters:
 *   name   <-- name that should be assigned to the nodal mesh
 *   dim    <-- spatial dimension
 *   order  <-- order
 *
 * returns:
 *  pointer to created nodal mesh representation structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_t *
fvmc_nodal_create(const char  *name,
                  int          dim,
                  int          order)
{
  fvmc_nodal_t  *this_nodal;

  BFTC_MALLOC(this_nodal, 1, fvmc_nodal_t);

  /* Global indicators */

  if (name != NULL) {
    BFTC_MALLOC(this_nodal->name, strlen(name) + 1, char);
    strcpy(this_nodal->name, name);
  }
  else
    this_nodal->name = NULL;

  this_nodal->dim     = dim;
  this_nodal->order   = order;
  this_nodal->num_dom = fvmc_parall_get_rank() + 1;
  this_nodal->n_doms  = fvmc_parall_get_size();
  this_nodal->n_sections = 0;
  this_nodal->sections_idx = NULL;

  this_nodal->ho_uvw_to_local_ordering = NULL;
  this_nodal->ho_user_to_uvw = NULL;
  
  /* Local dimensions */

  this_nodal->n_cells = 0;
  this_nodal->n_faces = 0;
  this_nodal->n_edges = 0;
  this_nodal->n_vertices = 0;

  /* Local structures */

  this_nodal->vertex_coords = NULL;
  this_nodal->_vertex_coords = NULL;

  this_nodal->parent_vertex_num = NULL;
  this_nodal->_parent_vertex_num = NULL;

  this_nodal->global_vertex_num = NULL;

  this_nodal->sections = NULL;

  return (this_nodal);

}

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh representation structure.
 *
 * parameters:
 *   this_nodal  <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_nodal_t *
fvmc_nodal_destroy(fvmc_nodal_t  * this_nodal)
{

  /* Local structures */

  if (this_nodal->name != NULL)
    BFTC_FREE(this_nodal->name);

  if (this_nodal->_vertex_coords != NULL)
    BFTC_FREE(this_nodal->_vertex_coords);

  if (this_nodal->parent_vertex_num != NULL) {
    this_nodal->parent_vertex_num = NULL;
    BFTC_FREE(this_nodal->_parent_vertex_num);
  }

  if (this_nodal->global_vertex_num != NULL)
    fvmc_io_num_destroy(this_nodal->global_vertex_num);

  for (int i = 0; i < this_nodal->n_sections; i++)
    fvmc_nodal_section_destroy(this_nodal->sections[i]);

  if (this_nodal->sections != NULL)
    BFTC_FREE(this_nodal->sections);

  if (this_nodal->sections_idx != NULL)
    BFTC_FREE(this_nodal->sections_idx);
    
  /* Main structure destroyed and NULL returned */

  if (this_nodal->ho_uvw_to_local_ordering != NULL) {
    for (int i = 0; i <  FVMC_N_ELEMENT_TYPES; i++) {
      if (this_nodal->ho_uvw_to_local_ordering[i] != NULL) {
        free (this_nodal->ho_uvw_to_local_ordering[i]);
      }
    }
    free (this_nodal->ho_uvw_to_local_ordering);
  }
  
  if (this_nodal->ho_user_to_uvw != NULL) {
    for (int i = 0; i <  FVMC_N_ELEMENT_TYPES; i++) {
      if (this_nodal->ho_user_to_uvw[i] != NULL) {
        free (this_nodal->ho_user_to_uvw[i]);
      }
    }
    free (this_nodal->ho_user_to_uvw);
  }

  BFTC_FREE(this_nodal);

  return (this_nodal);
}

/*----------------------------------------------------------------------------
 * Copy a nodal mesh representation structure, sharing arrays with the
 * original structure.
 *
 * parameters:
 *   this_nodal  <-> pointer to structure that should be copied
 *
 * returns:
 *   pointer to created nodal mesh representation structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_t *
fvmc_nodal_copy(const fvmc_nodal_t *this_nodal)
{
  int i;
  fvmc_nodal_t  *new_nodal;

  BFTC_MALLOC(new_nodal, 1, fvmc_nodal_t);

  /* Global indicators */

  if (this_nodal->name != NULL) {
    BFTC_MALLOC(new_nodal->name, strlen(this_nodal->name) + 1, char);
    strcpy(new_nodal->name, this_nodal->name);
  }
  else
    new_nodal->name = NULL;

  new_nodal->dim     = this_nodal->dim;
  new_nodal->num_dom = this_nodal->num_dom;
  new_nodal->n_doms  = this_nodal->n_doms;
  new_nodal->n_sections = this_nodal->n_sections;

  /* Local dimensions */

  new_nodal->n_cells = this_nodal->n_cells;
  new_nodal->n_faces = this_nodal->n_faces;
  new_nodal->n_edges = this_nodal->n_edges;
  new_nodal->n_vertices = this_nodal->n_vertices;

  /* Local structures */

  new_nodal->vertex_coords = this_nodal->vertex_coords;
  new_nodal->_vertex_coords = NULL;

  new_nodal->parent_vertex_num = this_nodal->parent_vertex_num;
  new_nodal->_parent_vertex_num = NULL;

  if (this_nodal->global_vertex_num != NULL) {
    fvmc_lnum_t n_ent
      = fvmc_io_num_get_local_count(this_nodal->global_vertex_num);
    fvmc_gnum_t global_count
      = fvmc_io_num_get_global_count(this_nodal->global_vertex_num);
    const fvmc_gnum_t *global_num
      = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);

    new_nodal->global_vertex_num
      = fvmc_io_num_create_shared(global_num, global_count, n_ent);
  }
  else
    new_nodal->global_vertex_num = NULL;


  BFTC_MALLOC(new_nodal->sections,
             new_nodal->n_sections,
             fvmc_nodal_section_t *);
  for (i = 0; i < new_nodal->n_sections; i++)
    new_nodal->sections[i] = _fvmc_nodal_section_copy(this_nodal->sections[i]);

  BFTC_MALLOC(new_nodal->sections_idx,
             new_nodal->n_sections + 1,
             int);
  for (i = 0; i < new_nodal->n_sections+1; i++)
    new_nodal->sections_idx[i] = this_nodal->sections_idx[i];

  return (new_nodal);
}

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh representation structure: only the associations
 * (numberings) necessary to redistribution of fields for output are
 * conserved, the full connectivity being in many cases no longer useful
 * once it has been output. If the del_vertex_num value is set
 * to true, vertex-based values may not be output in parallel mode
 * after this function is called.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be reduced
 *   del_vertex_num    <-- indicates if vertex parent indirection and
 *                         I/O numbering are destroyed (1) or not (0)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_reduce(fvmc_nodal_t  *this_nodal,
                 int           del_vertex_num)
{
  int  i;
  _Bool reduce_vertices = true;

  /* Connectivity */

  for (i = 0; i < this_nodal->n_sections; i++) {
    if (_fvmc_nodal_section_reduce(this_nodal->sections[i]) == false)
      reduce_vertices = false;
  }

  /* Vertices */

  if (reduce_vertices == true) {

    if (this_nodal->_vertex_coords != NULL)
      BFTC_FREE(this_nodal->_vertex_coords);
    this_nodal->vertex_coords = NULL;

  }

  /* Depending on this option, output on vertices may not remain possible */

  if (del_vertex_num > 0) {

    if (this_nodal->parent_vertex_num != NULL) {
      this_nodal->parent_vertex_num = NULL;
      BFTC_FREE(this_nodal->_parent_vertex_num);
    }

    if (this_nodal->global_vertex_num != NULL)
      this_nodal->global_vertex_num
        = fvmc_io_num_destroy(this_nodal->global_vertex_num);

  }

}

/*----------------------------------------------------------------------------
 * Change entity parent numbering; this is useful when entities of the
 * parent mesh have been renumbered after a nodal mesh representation
 * structure's creation.
 *
 * parameters:
 *   this_nodal          <-- nodal mesh structure
 *   new_parent_num      <-- pointer to local parent renumbering array
 *                           ({1, ..., n} <-- {1, ..., n})
 *   entity_dim          <-- 3 for cells, 2 for faces, 1 for edges,
 *                           and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_change_parent_num(fvmc_nodal_t       *this_nodal,
                            const fvmc_lnum_t   new_parent_num[],
                            int                entity_dim)
{
  /* Vertices */

  if (entity_dim == 0) {

    this_nodal->_parent_vertex_num
      = _renumber_parent_num(this_nodal->n_vertices,
                             new_parent_num,
                             this_nodal->parent_vertex_num,
                             this_nodal->_parent_vertex_num);
    this_nodal->parent_vertex_num = this_nodal->_parent_vertex_num;

  }

  /* Other elements */

  else {

    int  i = 0;
    fvmc_nodal_section_t  *section = NULL;

    for (i = 0; i < this_nodal->n_sections; i++) {
      section = this_nodal->sections[i];
      if (section->entity_dim == entity_dim) {
        section->_parent_element_num
          = _renumber_parent_num(section->n_elements,
                                 new_parent_num,
                                 section->parent_element_num,
                                 section->_parent_element_num);
        section->parent_element_num = section->_parent_element_num;
      }
    }

  }

}

/*----------------------------------------------------------------------------
 * Remove entity parent numbering; this is useful for example when we
 * want to assign coordinates or fields to an extracted mesh using
 * arrays relative to the mesh, and not to its parent.
 *
 * This is equivalent to calling fvmc_nodal_change_parent_num(), with
 * 'trivial' (1 o n) new_parent_num[] values.
 *
 * parameters:
 *   this_nodal          <-- nodal mesh structure
 *   entity_dim          <-- 3 for cells, 2 for faces, 1 for edges,
 *                           and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_remove_parent_num(fvmc_nodal_t  *this_nodal,
                            int           entity_dim)
{
  /* Vertices */

  if (entity_dim == 0) {
    this_nodal->parent_vertex_num = NULL;
    if (this_nodal->_parent_vertex_num != NULL)
      BFTC_FREE(this_nodal->_parent_vertex_num);
  }

  /* Other elements */

  else {

    int  i = 0;
    fvmc_nodal_section_t  *section = NULL;

    for (i = 0; i < this_nodal->n_sections; i++) {
      section = this_nodal->sections[i];
      if (section->entity_dim == entity_dim) {
        section->parent_element_num = NULL;
        if (section->_parent_element_num != NULL)
          BFTC_FREE(section->_parent_element_num);
      }
    }

  }

}

/*----------------------------------------------------------------------------
 * Build external numbering for entities based on global numbers.
 *
 * parameters:
 *   this_nodal           <-- nodal mesh structure
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   entity_dim           <-- 3 for cells, 2 for faces, 1 for edges,
 *                            and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_init_io_num(fvmc_nodal_t       *this_nodal,
                      const fvmc_gnum_t   parent_global_numbers[],
                      int                entity_dim)
{
  int  i;
  fvmc_nodal_section_t  *section;

  if (entity_dim == 0)
    this_nodal->global_vertex_num
      = fvmc_io_num_create(this_nodal->parent_vertex_num,
                          parent_global_numbers,
                          this_nodal->n_vertices,
                          0);

  else {
    for (i = 0; i < this_nodal->n_sections; i++) {
      section = this_nodal->sections[i];
      if (section->entity_dim == entity_dim) {
        section->global_element_num
          = fvmc_io_num_create(section->parent_element_num,
                              parent_global_numbers,
                              section->n_elements,
                              0);
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Preset number and list of vertices to assign to a nodal mesh.
 *
 * If the parent_vertex_num argument is NULL, the list is assumed to
 * be {1, 2, ..., n}. If parent_vertex_num is given, it specifies a
 * list of n vertices from a larger set (1 to n numbering).
 *
 * Ownership of the given parent vertex numbering array is
 * transferred to the nodal mesh representation structure.
 *
 * This function should be called before fvmc_nodal_set_shared_vertices()
 * or fvmc_nodal_transfer_vertices() if we want to force certain
 * vertices to appear in the mesh (especially if we want to define
 * a mesh containing only vertices).
 *
 * parameters:
 *   this_nodal        <-> nodal mesh structure
 *   n_vertices        <-- number of vertices to assign
 *   parent_vertex_num <-- parent numbers of vertices to assign
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_define_vertex_list(fvmc_nodal_t  *this_nodal,
                             fvmc_lnum_t    n_vertices,
                             fvmc_lnum_t    parent_vertex_num[])
{
  assert(this_nodal != NULL);

  this_nodal->n_vertices = n_vertices;

  this_nodal->parent_vertex_num = NULL;
  if (this_nodal->_parent_vertex_num != NULL)
    BFTC_FREE(this_nodal->_parent_vertex_num);

  if (parent_vertex_num != NULL) {
    this_nodal->_parent_vertex_num = parent_vertex_num;
    this_nodal->parent_vertex_num = parent_vertex_num;
  }

}

/*----------------------------------------------------------------------------
 * Assign shared vertex coordinates to an extracted nodal mesh,
 * renumbering vertex numbers based on those really referenced,
 * and updating connectivity arrays in accordance.
 *
 * This function should only be called once all element sections
 * have been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   vertex_coords   <-- coordinates of parent vertices (interlaced)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_set_shared_vertices(fvmc_nodal_t        *this_nodal,
                              const fvmc_coord_t   vertex_coords[])
{
  assert(this_nodal != NULL);

  /* Map vertex coordinates to array passed as argument
     (this_nodal->_vertex_coords remains NULL, so only
     the const pointer may be used for a shared array) */

  this_nodal->vertex_coords = vertex_coords;

  /* If the mesh contains only vertices, its n_vertices and
     parent_vertex_num must already have been set, and do not
     require updating */

  if (this_nodal->n_sections == 0)
    return;

  /* Renumber vertices based on those really referenced */

  _renumber_vertices(this_nodal);

}

/*----------------------------------------------------------------------------
 * Assign private vertex coordinates to a nodal mesh,
 * renumbering vertex numbers based on those really referenced,
 * and updating connectivity arrays in accordance.
 *
 * Ownership of the given coordinates array is transferred to
 * the nodal mesh representation structure.
 *
 * This function should only be called once all element sections
 * have been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   vertex_coords   <-- coordinates of parent vertices (interlaced)
 *
 * returns:
 *   updated pointer to vertex_coords (may be different from initial
 *   argument if vertices were renumbered).
 *----------------------------------------------------------------------------*/

fvmc_coord_t *
fvmc_nodal_transfer_vertices(fvmc_nodal_t  *this_nodal,
                            fvmc_coord_t   vertex_coords[])
{
  fvmc_lnum_t  i;
  int         j;

  fvmc_coord_t  *_vertex_coords = vertex_coords;

  assert(this_nodal != NULL);

  /* Renumber vertices based on those really referenced, and
     update connectivity arrays in accordance. */

  _renumber_vertices(this_nodal);

  /* If renumbering is necessary, update connectivity */

  if (this_nodal->parent_vertex_num != NULL) {

    int dim = this_nodal->dim;
    const fvmc_lnum_t *parent_vertex_num = this_nodal->parent_vertex_num;

    BFTC_MALLOC(_vertex_coords, this_nodal->n_vertices * dim, fvmc_coord_t);

    for (i = 0; i < this_nodal->n_vertices; i++) {
      for (j = 0; j < dim; j++)
        _vertex_coords[i*dim + j]
          = vertex_coords[(parent_vertex_num[i]-1)*dim + j];
    }

    BFTC_FREE(vertex_coords);

    this_nodal->parent_vertex_num = NULL;
    if (this_nodal->_parent_vertex_num != NULL)
      BFTC_FREE(this_nodal->_parent_vertex_num);
  }

  this_nodal->_vertex_coords = _vertex_coords;
  this_nodal->vertex_coords = _vertex_coords;

  return _vertex_coords;
}


/*----------------------------------------------------------------------------
 * return type of an element
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element              <-- element (1 to n numbering).
 *
 * returns:
 *   type
 *----------------------------------------------------------------------------*/

fvmc_element_t
fvmc_nodal_get_type_elt(const fvmc_nodal_t  *this_nodal, const int elt)
{
  assert(this_nodal != NULL);

  int _elt = elt-1;

  int elt_section = 0;
  
  while (_elt >= this_nodal->sections_idx[++elt_section]);

  elt_section -= 1 ;

  return this_nodal->sections[elt_section]->type;
}

/*----------------------------------------------------------------------------
 * return local to user numbering
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element              <-- element (1 to n numbering).
 *
 * returns:
 *   local to user numbering
 *----------------------------------------------------------------------------*/

const int*
fvmc_nodal_get_local_to_user_numbering_elt (const fvmc_nodal_t  *this_nodal, const int elt)
{
  assert(this_nodal != NULL);

  int _elt = elt-1;

  int elt_section = 0;
  
  while (_elt >= this_nodal->sections_idx[++elt_section]);

  elt_section -= 1 ;

  return this_nodal->sections[elt_section]->ho_local_to_user_ordering;
}


/*----------------------------------------------------------------------------
 * return internal connectivity
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element              <-- element (1 to n numbering).
 *
 * returns:
 *   type
 *----------------------------------------------------------------------------*/

const int *
fvmc_nodal_get_internal_connec_elt(const fvmc_nodal_t  *this_nodal, const int elt)
{
  assert(this_nodal != NULL);

  int _elt = elt-1;

  int elt_section = 0;

  int *_vertex_num = NULL;

  //  printf("this_nodal->sections_idx %d\n", this_nodal->sections_idx[_elt]);
  
  while (_elt >= this_nodal->sections_idx[++elt_section]) {
    printf("this_nodal->sections_idx %d %d\n", _elt, this_nodal->sections_idx[elt_section]);
  }

  elt_section -= 1 ;
  int elt_loc = _elt - this_nodal->sections_idx[elt_section];

  if (this_nodal->sections[elt_section]->_ho_vertex_num != NULL) {
    _vertex_num = this_nodal->sections[elt_section]->_ho_vertex_num +
      elt_loc * this_nodal->sections[elt_section]->stride;
  }
  else {
    _vertex_num = this_nodal->sections[elt_section]->_vertex_num +
      elt_loc * this_nodal->sections[elt_section]->stride;
  }
  
  return _vertex_num;
}

/*----------------------------------------------------------------------------
 * return connectivity
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element              <-- element (1 to n numbering).
 *
 * returns:
 *   type
 *----------------------------------------------------------------------------*/

const int *
fvmc_nodal_get_connec_elt(const fvmc_nodal_t  *this_nodal, const int elt)
{
  assert(this_nodal != NULL);

  int _elt = elt-1;

  int elt_section = 0;
  
  while (_elt >= this_nodal->sections_idx[++elt_section]);

  elt_section -= 1 ;
  int elt_loc = _elt - this_nodal->sections_idx[elt_section];

  return this_nodal->sections[elt_section]->_vertex_num +
    elt_loc * this_nodal->sections[elt_section]->stride;

}

/*----------------------------------------------------------------------------
 * return the number of nodes of an element
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element              <-- element (1 to n numbering).
 *
 * returns:
 *   number of nodes
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_get_n_node_elt(const fvmc_nodal_t  *this_nodal, const int elt)
{
  assert(this_nodal != NULL);

  int _elt = elt-1;

  int elt_section = 0;
  
  while (_elt >= this_nodal->sections_idx[++elt_section]);

  elt_section -= 1 ;

  return this_nodal->sections[elt_section]->stride;
}

/*----------------------------------------------------------------------------
 * Obtain the name of a nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   pointer to constant string containing the mesh name
 *----------------------------------------------------------------------------*/

const char *
fvmc_nodal_get_name(const fvmc_nodal_t  *this_nodal)
{
  assert(this_nodal != NULL);

  return this_nodal->name;
}

/*----------------------------------------------------------------------------
 * Return spatial dimension of the nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *
 * returns:
 *  spatial dimension.
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_get_dim(const fvmc_nodal_t  *this_nodal)
{
  return this_nodal->dim;
}

/*----------------------------------------------------------------------------
 * Return maximum dimension of entities in a nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *
 * returns:
 *  maximum dimension of entities in mesh (0 to 3)
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_get_max_entity_dim(const fvmc_nodal_t  *this_nodal)
{
  int  section_id;
  int  max_entity_dim = 0;

  assert(this_nodal != NULL);

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {
    const fvmc_nodal_section_t  *section = this_nodal->sections[section_id];
    if (section->entity_dim > max_entity_dim)
      max_entity_dim = section->entity_dim;
  }

  return max_entity_dim;
}

/*----------------------------------------------------------------------------
 * Return number of entities of a given dimension in a nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   entity_dim <-- dimension of entities we want to count (0 to 3)
 *
 * returns:
 *  number of entities of given dimension in mesh
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_nodal_get_n_entities(const fvmc_nodal_t  *this_nodal,
                         int                 entity_dim)
{
  fvmc_lnum_t n_entities;

  assert(this_nodal != NULL);

  switch(entity_dim) {
  case 0:
    n_entities = this_nodal->n_vertices;
    break;
  case 1:
    n_entities = this_nodal->n_edges;
    break;
  case 2:
    n_entities = this_nodal->n_faces;
    break;
  case 3:
    n_entities = this_nodal->n_cells;
    break;
  default:
    n_entities = 0;
  }

  return n_entities;
}

/*----------------------------------------------------------------------------
 * Return global number of vertices associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   global number of vertices associated with nodal mesh
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_nodal_get_n_g_vertices(const fvmc_nodal_t  *this_nodal)
{
  return fvmc_nodal_n_g_vertices(this_nodal);
}

/*----------------------------------------------------------------------------
 * Return global number of elements of a given type associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element_type         <-- type of elements for query
 *
 * returns:
 *   global number of elements of the given type associated with nodal mesh
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_nodal_get_n_g_elements(const fvmc_nodal_t  *this_nodal,
                           fvmc_element_t       element_type)
{
  int  i;
  fvmc_gnum_t  n_g_elements = 0;

  assert(this_nodal != NULL);

  for (i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t  *section = this_nodal->sections[i];
    if (section->type == element_type)
      n_g_elements += fvmc_nodal_section_n_g_elements(section);
  }

  return n_g_elements;
}

/*----------------------------------------------------------------------------
 * Return local number of elements of a given type associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element_type         <-- type of elements for query
 *
 * returns:
 *   local number of elements of the given type associated with nodal mesh
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_nodal_get_n_elements(const fvmc_nodal_t  *this_nodal,
                         fvmc_element_t       element_type)
{
  int  i;
  fvmc_lnum_t  n_elements = 0;

  assert(this_nodal != NULL);

  for (i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t  *section = this_nodal->sections[i];
    if (section->type == element_type)
      n_elements += section->n_elements;
  }

  return n_elements;
}

/*----------------------------------------------------------------------------
 * Return local parent numbering array for all entities of a given
 * dimension in a nodal mesh.
 *
 * The number of entities of the given dimension may be obtained
 * through fvmc_nodal_get_n_entities(), the parent_num[] array is populated
 * with the parent entity numbers of those entities, in order (i.e. in
 * local section order, section by section).
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   entity_dim <-- dimension of entities we are interested in (0 to 3)
 *   parent_num --> entity parent numbering (array must be pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_get_parent_num(const fvmc_nodal_t  *this_nodal,
                         int                 entity_dim,
                         fvmc_lnum_t          parent_num[])
{
  int section_id;
  fvmc_lnum_t i;

  fvmc_lnum_t entity_count = 0;

  assert(this_nodal != NULL);

  /* Entity dimension 0: vertices */

  if (entity_dim == 0) {
    if (this_nodal->parent_vertex_num != NULL) {
      for (i = 0; i < this_nodal->n_vertices; i++)
        parent_num[entity_count++] = this_nodal->parent_vertex_num[i];
    }
    else {
      for (i = 0; i < this_nodal->n_vertices; i++)
        parent_num[entity_count++] = i + 1;
    }
  }

  /* Entity dimension > 0: edges, faces, or cells */

  else {

    for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

      const fvmc_nodal_section_t  *section = this_nodal->sections[section_id];

      if (section->entity_dim == entity_dim) {
        if (section->parent_element_num != NULL) {
          for (i = 0; i < section->n_elements; i++)
            parent_num[entity_count++] = section->parent_element_num[i];
        }
        else {
          for (i = 0; i < section->n_elements; i++)
            parent_num[entity_count++] = i + 1;
        }
      }

    } /* end loop on sections */

  }
}

/*----------------------------------------------------------------------------
 * Compute tesselation a a nodal mesh's sections of a given type, and add the
 * corresponding structure to the mesh representation.
 *
 * If global element numbers are used (i.e. in parallel mode), this function
 * should be only be used after calling fvmc_nodal_init_io_num().
 *
 * If some mesh sections have already been tesselated, their tesselation
 * is unchanged.
 *
 * parameters:
 *   this_nodal  <-> pointer to nodal mesh structure
 *   type        <-> element type that should be tesselated
 *   error_count --> number of elements with a tesselation error
 *                   counter (optional)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_tesselate(fvmc_nodal_t    *this_nodal,
                    fvmc_element_t   type,
                    fvmc_lnum_t     *error_count)
{
  int section_id;
  fvmc_lnum_t section_error_count;

  assert(this_nodal != NULL);

  if (error_count != NULL)
    *error_count = 0;

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    fvmc_nodal_section_t  *section = this_nodal->sections[section_id];

    if (section->type == type && section->tesselation == NULL) {

      section->tesselation = fvmc_tesselation_create(type,
                                                    section->n_elements,
                                                    section->face_index,
                                                    section->face_num,
                                                    section->vertex_index,
                                                    section->vertex_num,
                                                    section->global_element_num);

      fvmc_tesselation_init(section->tesselation,
                           this_nodal->dim,
                           this_nodal->vertex_coords,
                           this_nodal->parent_vertex_num,
                           &section_error_count);

      if (error_count != NULL)
        *error_count += section_error_count;
    }

  }
}

/*----------------------------------------------------------------------------
 * Build a nodal representation structure based on extraction of a
 * mesh's edges.
 *
 * parameters:
 *   name        <-- name to assign to extracted mesh
 *   this_nodal  <-> pointer to nodal mesh structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_t *
fvmc_nodal_copy_edges(const char         *name,
                     const fvmc_nodal_t  *this_nodal)
{
  int i;
  fvmc_lnum_t j, k;
  fvmc_lnum_t n_edges = 0, n_max_edges = 0;
  fvmc_nodal_t *new_nodal = NULL;
  fvmc_nodal_section_t *new_section = NULL;

  BFTC_MALLOC(new_nodal, 1, fvmc_nodal_t);

  /* Global indicators */

  if (name != NULL) {
    BFTC_MALLOC(new_nodal->name, strlen(name) + 1, char);
    strcpy(new_nodal->name, name);
  }
  else
    new_nodal->name = NULL;

  new_nodal->dim     = this_nodal->dim;
  new_nodal->num_dom = this_nodal->num_dom;
  new_nodal->n_doms  = this_nodal->n_doms;
  new_nodal->n_sections = 1;

  /* Local dimensions */

  new_nodal->n_cells = 0;
  new_nodal->n_faces = 0;
  new_nodal->n_edges = 0;
  new_nodal->n_vertices = this_nodal->n_vertices;

  /* Local structures */

  new_nodal->vertex_coords = this_nodal->vertex_coords;
  new_nodal->_vertex_coords = NULL;

  new_nodal->parent_vertex_num = this_nodal->parent_vertex_num;
  new_nodal->_parent_vertex_num = NULL;

  if (this_nodal->global_vertex_num != NULL) {
    fvmc_lnum_t n_ent
      = fvmc_io_num_get_local_count(this_nodal->global_vertex_num);
    fvmc_gnum_t global_count
      = fvmc_io_num_get_global_count(this_nodal->global_vertex_num);
    const fvmc_gnum_t *global_num
      = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);

    new_nodal->global_vertex_num
      = fvmc_io_num_create_shared(global_num, global_count, n_ent);
  }
  else
    new_nodal->global_vertex_num = NULL;

  /* Counting step */

  for (i = 0; i < this_nodal->n_sections; i++) {
    const fvmc_nodal_section_t *this_section = this_nodal->sections[i];
    if (this_section->vertex_index == NULL)
      n_max_edges += (  fvmc_nodal_n_edges_element[this_section->type]
                      * this_section->n_elements);
    else if (this_section->type == FVMC_FACE_POLY)
      n_max_edges += this_section->vertex_index[this_section->n_elements];
    else if (this_section->type == FVMC_CELL_POLY)
      n_max_edges += this_section->vertex_index[this_section->n_faces];
  }

  BFTC_MALLOC(new_nodal->sections, 1, fvmc_nodal_section_t *);

  int order = -1;
  
  new_section = fvmc_nodal_section_create(FVMC_EDGE,order);
  new_nodal->sections[0] = new_section;

  BFTC_MALLOC(new_section->_vertex_num, n_max_edges*2, fvmc_lnum_t);

  /* Add edges */

  for (i = 0; i < this_nodal->n_sections; i++) {

    const fvmc_nodal_section_t *this_section = this_nodal->sections[i];

    if (this_section->order != -1) {
        bftc_error(__FILE__, __LINE__, 0,
                  _("fvmc_nodal_copy_edges : element order > 1 is not taking into account"));
    }

    
    if (   this_section->type == FVMC_FACE_POLY
        || this_section->type == FVMC_CELL_POLY) {

      fvmc_lnum_t n_faces = this_section->type == FVMC_FACE_POLY ?
        this_section->n_elements : this_section->n_faces;

      for (j = 0; j < n_faces; j++) {
        const fvmc_lnum_t face_start_id = this_section->vertex_index[j];
        const fvmc_lnum_t n_face_edges
          = this_section->vertex_index[j+1] - this_section->vertex_index[j];
        for (k = 0; k < n_face_edges; k++) {
          new_section->_vertex_num[n_edges*2]
            = this_section->vertex_num[face_start_id + k];
          new_section->_vertex_num[n_edges*2 + 1]
            = this_section->vertex_num[face_start_id + (k + 1)%n_face_edges];
          n_edges += 1;
        }
      }

    }
    else {

      fvmc_lnum_t edges[2][12];

      fvmc_lnum_t n_elt_edges = fvmc_nodal_n_edges_element[this_section->type];
      fvmc_lnum_t n_elts = this_section->n_elements;
      fvmc_lnum_t stride = this_section->stride;

      switch (this_section->type) {

      case FVMC_EDGE:
      case FVMC_FACE_TRIA:
      case FVMC_FACE_QUAD:
        for (j = 0; j < n_elt_edges; j++) {
          edges[0][j] = j;
          edges[1][j] = (j+1)%n_elt_edges;
        }
        break;

      case FVMC_CELL_TETRA:
        edges[0][0] = 0; edges[1][0] = 1;
        edges[0][1] = 1; edges[1][1] = 2;
        edges[0][2] = 2; edges[1][2] = 0;
        edges[0][3] = 0; edges[1][3] = 3;
        edges[0][4] = 1; edges[1][4] = 3;
        edges[0][5] = 2; edges[1][5] = 3;
        break;

      case FVMC_CELL_PYRAM:
        edges[0][0] = 0; edges[1][0] = 1;
        edges[0][1] = 1; edges[1][1] = 2;
        edges[0][2] = 2; edges[1][2] = 3;
        edges[0][3] = 3; edges[1][3] = 0;
        edges[0][4] = 0; edges[1][4] = 4;
        edges[0][5] = 1; edges[1][5] = 4;
        edges[0][6] = 2; edges[1][6] = 4;
        edges[0][7] = 3; edges[1][7] = 4;
        break;

      case FVMC_CELL_PRISM:
        edges[0][0] = 0; edges[1][0] = 1;
        edges[0][1] = 1; edges[1][1] = 2;
        edges[0][2] = 2; edges[1][2] = 0;
        edges[0][3] = 0; edges[1][3] = 3;
        edges[0][4] = 1; edges[1][4] = 4;
        edges[0][5] = 2; edges[1][5] = 5;
        edges[0][6] = 3; edges[1][6] = 4;
        edges[0][7] = 4; edges[1][7] = 5;
        edges[0][8] = 5; edges[1][8] = 3;
        break;

      case FVMC_CELL_HEXA:
        edges[0][0] = 0; edges[1][0] = 1;
        edges[0][1] = 1; edges[1][1] = 2;
        edges[0][2] = 2; edges[1][2] = 3;
        edges[0][3] = 3; edges[1][3] = 0;
        edges[0][4] = 0; edges[1][4] = 4;
        edges[0][5] = 1; edges[1][5] = 5;
        edges[0][6] = 2; edges[1][6] = 6;
        edges[0][7] = 3; edges[1][7] = 7;
        edges[0][8] = 4; edges[1][8] = 5;
        edges[0][9] = 5; edges[1][9] = 6;
        edges[0][10] = 6; edges[1][10] = 7;
        edges[0][11] = 7; edges[1][11] = 4;
        break;

      default:
        assert(0);
        edges[0][0] = -1; /* For nonempty default clause */
      }

      for (j = 0; j < n_elts; j++) {
        const fvmc_lnum_t *_vertex_num = this_section->vertex_num + (j*stride);
        for (k = 0; k < n_elt_edges; k++) {
          new_section->_vertex_num[n_edges*2] = _vertex_num[edges[0][k]];
          new_section->_vertex_num[n_edges*2 + 1] = _vertex_num[edges[1][k]];
          n_edges += 1;
        }
      }
    }
  } /* End of loop on sections */

  assert(n_edges == n_max_edges);

  /* Ensure edges are oriented in the same direction */

  if (this_nodal->global_vertex_num != NULL) {

    const fvmc_gnum_t *v_num_g
      = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);

    for (j = 0; j < n_max_edges; j++) {
      fvmc_lnum_t vnum_1 = new_section->_vertex_num[j*2];
      fvmc_lnum_t vnum_2 = new_section->_vertex_num[j*2 + 1];
      if (v_num_g[vnum_1 - 1] > v_num_g[vnum_2 - 1]) {
        new_section->_vertex_num[j*2] = vnum_2;
        new_section->_vertex_num[j*2 + 1] = vnum_1;
      }

    }

  }
  else {

    for (j = 0; j < n_max_edges; j++) {
      fvmc_lnum_t vnum_1 = new_section->_vertex_num[j*2];
      fvmc_lnum_t vnum_2 = new_section->_vertex_num[j*2 + 1];
      if (vnum_1 > vnum_2) {
        new_section->_vertex_num[j*2] = vnum_2;
        new_section->_vertex_num[j*2 + 1] = vnum_1;
      }
    }

  }

  /* Sort and remove duplicates
     (use qsort rather than fvmc_order_local_s() so as to sort in place) */

  qsort(new_section->_vertex_num,
        n_max_edges,
        sizeof(fvmc_lnum_t) * 2,
        &_compare_edges);

  {
    fvmc_lnum_t vn_1_p = -1;
    fvmc_lnum_t vn_2_p = -1;

    n_edges = 0;

    for (j = 0; j < n_max_edges; j++) {

      fvmc_lnum_t vn_1 = new_section->_vertex_num[j*2];
      fvmc_lnum_t vn_2 = new_section->_vertex_num[j*2 + 1];

      if (vn_1 != vn_1_p || vn_2 != vn_2_p) {
        new_section->_vertex_num[n_edges*2]     = vn_1;
        new_section->_vertex_num[n_edges*2 + 1] = vn_2;
        vn_1_p = vn_1;
        vn_2_p = vn_2;
        n_edges += 1;
      }
    }
  }

  /* Resize edge connectivity to adjust to final size */

  BFTC_REALLOC(new_section->_vertex_num, n_edges*2, fvmc_lnum_t);
  new_section->vertex_num = new_section->_vertex_num;

  new_section->n_elements = n_edges;
  new_nodal->n_edges = n_edges;

  /* Build  global edge numbering if necessary */

  if (new_nodal->n_doms > 1) {

    fvmc_gnum_t *edge_vertices_g; /* edges -> global vertices */

    BFTC_MALLOC(edge_vertices_g, n_edges*2, fvmc_gnum_t);

    if (this_nodal->global_vertex_num != NULL) {
      const fvmc_gnum_t *v_num_g
        = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);
      for (j = 0; j < n_edges; j++) {
        edge_vertices_g[j*2]   = v_num_g[new_section->_vertex_num[j*2] - 1];
        edge_vertices_g[j*2+1] = v_num_g[new_section->_vertex_num[j*2+1] - 1];
      }
    }
    else {
      for (j = 0; j < n_edges; j++) {
        edge_vertices_g[j*2]     = new_section->_vertex_num[j*2];
        edge_vertices_g[j*2 + 1] = new_section->_vertex_num[j*2 + 1];
      }
    }

    new_section->global_element_num
      = fvmc_io_num_create_from_adj_s(NULL, edge_vertices_g, n_edges, 2);


    BFTC_FREE(edge_vertices_g);
  };

  return (new_nodal);
}

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure.
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_dump(const fvmc_nodal_t  *this_nodal)
{
  fvmc_lnum_t  i;
  fvmc_lnum_t  num_vertex = 1;
  const fvmc_coord_t  *coord = this_nodal->vertex_coords;

  /* Global indicators */
  /*--------------------*/

  bftc_printf("\n"
             "Mesh name:\"%s\"\n",
             this_nodal->name);

  bftc_printf("\n"
             "Mesh dimension:               %d\n"
             "Domain number:                %d\n"
             "Number of domains:            %d\n"
             "Number of sections:           %d\n",
             this_nodal->dim, this_nodal->num_dom, this_nodal->n_doms,
             this_nodal->n_sections);

  bftc_printf("\n"
             "Number of cells:               %d\n"
             "Number of faces:               %d\n"
             "Number of edges:               %d\n"
             "Number of vertices:            %d\n",
            this_nodal->n_cells,
            this_nodal->n_faces,
            this_nodal->n_edges,
            this_nodal->n_vertices);

  if (this_nodal->n_vertices > 0) {

    bftc_printf("\n"
               "Pointers to shareable arrays:\n"
               "  vertex_coords:        %p\n"
               "  parent_vertex_num:    %p\n",
               this_nodal->vertex_coords,
               this_nodal->parent_vertex_num);

    bftc_printf("\n"
               "Pointers to local arrays:\n"
               "  _vertex_coords:       %p\n"
               "  _parent_vertex_num:   %p\n",
               this_nodal->_vertex_coords,
               this_nodal->_parent_vertex_num);

    /* Output coordinates depending on parent numbering */

    if (this_nodal->parent_vertex_num == NULL) {

      bftc_printf("\nVertex coordinates:\n\n");
      switch(this_nodal->dim) {
      case 1:
        for (i = 0; i < this_nodal->n_vertices; i++)
          bftc_printf("%10d : %12.5f\n",
                     num_vertex++, (double)(coord[i]));
        break;
      case 2:
        for (i = 0; i < this_nodal->n_vertices; i++)
          bftc_printf("%10d : %12.5f %12.5f\n",
                     num_vertex++, (double)(coord[i*2]),
                     (double)(coord[i*2+1]));
        break;
      case 3:
        for (i = 0; i < this_nodal->n_vertices; i++)
          bftc_printf("%10d : %12.5f %12.5f %12.5f\n",
                     num_vertex++, (double)(coord[i*3]),
                     (double)(coord[i*3+1]), (double)(coord[i*3+2]));
        break;
      default:
        bftc_printf("coordinates not output\n"
                   "dimension = %d unsupported\n", this_nodal->dim);
      }

    }
    else { /* if (this_nodal->parent_vertex_num != NULL) */

      bftc_printf("\nVertex parent and coordinates:\n\n");

      switch(this_nodal->dim) {
      case 1:
        for (i = 0; i < this_nodal->n_vertices; i++) {
          coord =   this_nodal->vertex_coords
                  + (this_nodal->parent_vertex_num[i]-1);
          bftc_printf("%10d : %12.5f\n",
                     num_vertex++, (double)(coord[0]));
        }
        break;
      case 2:
        for (i = 0; i < this_nodal->n_vertices; i++) {
          coord =   this_nodal->vertex_coords
                  + ((this_nodal->parent_vertex_num[i]-1)*2);
          bftc_printf("%10d : %12.5f %12.5f\n",
                     num_vertex++, (double)(coord[0]), (double)(coord[1]));
        }
        break;
      case 3:
        for (i = 0; i < this_nodal->n_vertices; i++) {
          coord =   this_nodal->vertex_coords
                  + ((this_nodal->parent_vertex_num[i]-1)*3);
          bftc_printf("%10d : %12.5f %12.5f %12.5f\n",
                     num_vertex++, (double)(coord[0]), (double)(coord[1]),
                     (double)(coord[2]));
        }
        break;
      default:
        bftc_printf("coordinates not output\n"
                   "dimension = %d unsupported\n", this_nodal->dim);
      }

    }

  }

  /* Global vertex numbers (only for parallel execution) */
  if (this_nodal->global_vertex_num != NULL) {
    bftc_printf("\nGlobal vertex numbers:\n\n");
    fvmc_io_num_dump(this_nodal->global_vertex_num);
  }

  /* Dump element sections */
  /*-----------------------*/
  
  for (i = 0; i < this_nodal->n_sections; i++)
    _fvmc_nodal_section_dump(this_nodal->sections[i]);
  
}


/* Copie dans des tableaux de toutes les connectivites sauf polyedres pour couplage cwipi */

void
fvmc_nodal_get_vertex(const fvmc_nodal_t  *this_nodal,
                     fvmc_lnum_t *n_elts,
                     fvmc_lnum_t **vertices_index, 
                     fvmc_lnum_t **vertices)
{
  int ind_poly = 0;
  int s_vertices = 0;
  
  *n_elts = 0;
  for (int i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t *section = this_nodal->sections[i];
    if (section->type != FVMC_CELL_POLY) {
      if (ind_poly == 1)
        bftc_error(__FILE__, __LINE__, 0,
                  _("Fonction temporaire pour couplage saturne/cedre :"
                    " Les polyedres sont attendus a la fin"));
      *n_elts += section->n_elements;
      if (section->stride == 0)
        s_vertices += section->vertex_index[section->n_elements];
      else
        s_vertices += section->stride * section->n_elements;
    }
    if (section->type == FVMC_CELL_POLY)
      bftc_error(__FILE__, __LINE__, 0,
                _("This function is not implemented yet for polyedra"));
  }
  
  BFTC_MALLOC(*vertices_index, *n_elts+1, fvmc_lnum_t);
  BFTC_MALLOC(*vertices, s_vertices, fvmc_lnum_t);
  
  (*vertices_index)[0] = 0;
  *n_elts = 0;
  
  for (int i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t *section = this_nodal->sections[i];
    if (section->type != FVMC_CELL_POLY) {
      int section_index = (*vertices_index)[*n_elts];
      int n_previous_elts = *n_elts;
      if (section->stride != 0) {
        for (int j = 0; j < section->n_elements; j++) {  
          (*vertices_index)[n_previous_elts + j + 1] = (*vertices_index)[n_previous_elts + j] + section->stride;
        }
        for (int k = 0; k <  section->n_elements * section->stride; k++)
          (*vertices)[section_index + k] = section->vertex_num[k];
      }
      else {
        for (int j = 0; j < section->n_elements; j++) 
          (*vertices_index)[n_previous_elts + j + 1] =  (*vertices_index)[n_previous_elts + j] + 
            section->vertex_index[j+1] - section->vertex_index[j];
        for (int k = 0; k < section->vertex_index[section->n_elements]; k++)
          (*vertices)[section_index + k] = section->vertex_num[k];
      }
      *n_elts += section->n_elements;
    }
  }
}


void
fvmc_nodal_get_coords(const fvmc_nodal_t  *this_nodal,
                     fvmc_lnum_t *n_vertex,
                     double **coords)
{
  *n_vertex = this_nodal->n_vertices;
  
  BFTC_MALLOC(*coords,this_nodal->n_vertices * this_nodal->dim , double);
  
  for (int i = 0; i < this_nodal->dim * this_nodal->n_vertices; i++) 
    (*coords)[i] = this_nodal->vertex_coords[i];

} 


//void
//fvmc_nodal_get_poly_vertex(const fvmc_nodal_t  *this_nodal,
//                          fvmc_lnum_t *n_elts,
//                          fvmc_lnum_t **faces,  
//                          fvmc_lnum_t **faces_index, 
//                          fvmc_lnum_t **vertices, 
//                          fvmc_lnum_t **vertices_index)
//{
//  bftc_error(__FILE__, __LINE__, 0,
//            _("fvmc_nodal_get_poly_verte : This function is not implemented yet"));
//
  /*   int s_faces = 0; */
/*   int s_vertices = 0; */

/*   *n_elts = 0; */
/*   for (int i = 0; i < this_nodal->n_sections; i++) { */
/*     fvmc_nodal_section_t *section = this_nodal->sections[i]; */
/*     if (section->type == FVMC_CELL_POLY) { */
/*       *n_elts += section->n_elements; */
/*       s_faces += section->face_index[section->n_elements]; */
/*     } */
/*   } */

/*   BFTC_MALLOC(*faces_index, *n_elts, fvmc_lnum_t); */
/*   BFTC_MALLOC(*faces, s_faces, fvmc_lnum_t); */

/*   faces_index[0] = 0; */
/*   for (int i = 0; i < this_nodal->n_sections; i++) { */
/*     fvmc_nodal_section_t *section = this_nodal->sections[i]; */
/*     if (section->type == FVMC_CELL_POLY) { */
/*       for (int j = 0; j < section->n_elements; j++)  */
/*         faces_index[j+1] =  faces_index[j] +  */
/*           section->faces_index[j+1] - section->faces_index[j]; */
/*       for (int k = 0; k < section->n_elements * section->vertex_index[section->n_elements]; k++) */
/*         vertices[vertices_index[j] + k] = section->vertex_num[k]; */
/*     } */
/*   } */
//}

/*----------------------------------------------------------------------------
 * Set high order ordering
 *
 * TODO: doc about orderin for each element
 *
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *   t_elt      <-- type of element
 *   n_nodes    <-- number of nodes
 *   uvw_grid   <-- uvw grid 
 *
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_ho_ordering_set (fvmc_nodal_t  *this_nodal,
                            const fvmc_element_t t_elt,
                            const int n_nodes,
                            const int *uvw_grid)
{

  int expected_n_nodes = fvmc_nodal_n_vertices_element (t_elt, this_nodal->order);
  int order = this_nodal->order;
  
  if (n_nodes != expected_n_nodes) {
    bftc_error(__FILE__, __LINE__, 0,
               _("The number of nodes is not that expected\n"));
  }
  
  if (this_nodal->ho_uvw_to_local_ordering == NULL) {
    this_nodal->ho_uvw_to_local_ordering = (int **) malloc (sizeof(int *) * FVMC_N_ELEMENT_TYPES);
    for (int i = 0; i < FVMC_N_ELEMENT_TYPES; i++) {
      this_nodal->ho_uvw_to_local_ordering[i] = NULL;
    }
  }
  
  if (this_nodal->ho_user_to_uvw == NULL) {
    this_nodal->ho_user_to_uvw = (int **) malloc (sizeof(int *) * FVMC_N_ELEMENT_TYPES);
    for (int i = 0; i < FVMC_N_ELEMENT_TYPES; i++) {
      this_nodal->ho_user_to_uvw[i] = NULL;
    }
  }

  if (this_nodal->ho_uvw_to_local_ordering[t_elt] == NULL) {
    this_nodal->ho_uvw_to_local_ordering[t_elt] = _uvw_to_local_ordering (this_nodal, t_elt);
  }

  if (this_nodal->ho_user_to_uvw[t_elt] == NULL) {
    
    int stride = 0;

    switch(t_elt) {
    case FVMC_EDGE:               /* Edge */
      stride = 1;
      break;
    case FVMC_FACE_TRIA:          /* Triangle */
    case FVMC_FACE_QUAD:          /* Quadrangle */
      stride = 2;
      break;
    case FVMC_CELL_TETRA:         /* Tetrahedron */
    case FVMC_CELL_PYRAM:         /* Pyramid */
    case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
    case FVMC_CELL_HEXA:         /* Hexahedron (brick) */
      stride = 3;
      break;
    default:  
      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_nodal_ho_ordering_set : high order unavailable for this element type\n"));
    }

    this_nodal->ho_user_to_uvw[t_elt] = malloc (sizeof(int) * n_nodes * stride);
    memcpy(this_nodal->ho_user_to_uvw[t_elt], uvw_grid, n_nodes * stride);
    
  }

  int *_ho_uvw_to_local_ordering = this_nodal->ho_uvw_to_local_ordering[t_elt];
  for (int i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t  *_section = this_nodal->sections[i];

    if (_section->type == t_elt) {
      if (_section->_ho_vertex_num == NULL) {
        _section->_ho_vertex_num = malloc (sizeof(int) * n_nodes * _section->n_elements);
      }

      if (_section->ho_local_to_user_ordering == NULL) {
        _section->ho_local_to_user_ordering = malloc (sizeof(int) * n_nodes);
      }
      
      int stride = _section->entity_dim;

      for (int k = 0; k < n_nodes; k++) {
        const int *_uvw = uvw_grid + k * stride;
        int idx = 0; 
        for (int l = 0; l < stride; l++) {
          idx += (int) pow((order+1),l) * _uvw[l];
        }
        int local_num = _ho_uvw_to_local_ordering[idx];
        _section->ho_local_to_user_ordering[local_num] = k;
      }
        
      for (int j = 0; j < _section->n_elements; j++) {
        int *_ho_vertex_num = _section->_ho_vertex_num + j * n_nodes;
        const int *_vertex_num = _section->vertex_num + j * n_nodes;
        for (int k = 0; k < n_nodes; k++) {
          _ho_vertex_num[k] = _vertex_num[_section->ho_local_to_user_ordering[k]];
        }
      }
      
    }
  }
  
}

/*----------------------------------------------------------------------------
 * Set high order ordering from the coordinates of the nodes of the reference element
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *   t_elt      <-- type of element
 *   n_nodes    <-- number of nodes
 *   coords     <-- coordinates of the nodes of the reference element
 *
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_ho_ordering_from_ref_elt_set (fvmc_nodal_t  *this_nodal,
                                         const fvmc_element_t t_elt,
                                         const int n_nodes,
                                         const double *coords)
{

  int stride = 0;
  int order = this_nodal->order;
  
  switch(t_elt) {
  case FVMC_EDGE:               /* Edge */
    stride = 1;
    break;
  case FVMC_FACE_TRIA:          /* Triangle */
  case FVMC_FACE_QUAD:          /* Quadrangle */
    stride = 2;
    break;
  case FVMC_CELL_TETRA:         /* Tetrahedron */
  case FVMC_CELL_PYRAM:         /* Pyramid */
  case FVMC_CELL_PRISM:         /* Prism (pentahedron) */
  case FVMC_CELL_HEXA:         /* Hexahedron (brick) */
    stride = 3;
    break;
  default:  
    bftc_error(__FILE__, __LINE__, 0,
                  _("_uvw_to_local_ordering : high order unavailable "));
  }

  int *_uvw_grid = malloc(sizeof(int) * stride * n_nodes);

  double *_local_coords = _local_ref_nodes (this_nodal, t_elt);
  int* _uvw_to_local = _uvw_to_local_ordering (this_nodal, t_elt);
  int* _local_to_uvw = malloc(sizeof(int) * stride * n_nodes);

  for (int i = 0; i < stride * n_nodes; i++) {
    if (_uvw_to_local[i] > -1) {
      int w = 0;
      int v = 0;
      int u = 0;
      int restw = _uvw_to_local[i];
      int restv = _uvw_to_local[i];

      if (stride > 2) {
        w = _uvw_to_local[i] / ((order+1) * (order+1));
        restw = _uvw_to_local[i] % ((order+1) * (order+1));
      }
      
      if (stride > 1) {
        v = restw / (order+1);
        restv = restw % (order+1);
      }

      u = restv % (order+1);
      
      _local_to_uvw [stride * _uvw_to_local[i]]     = u;
      if (stride > 1) {
        _local_to_uvw [stride * _uvw_to_local[i] + 1] = v;
        if (stride > 2) {
          _local_to_uvw [stride * _uvw_to_local[i] + 2] = w;
        }
      }
    }  
  }
  
  //TODO: n^2 complexity call pdm_merge_points for a n*log(n) complexity

  for (int i = 0; i < n_nodes; i++) {
    double x1 = coords[3*i];
    double y1 = coords[3*i+1];
    double z1 = coords[3*i+2];
    int is_find = 0;
    
    for (int j = 0; j < n_nodes; j++) {
      double x2 = _local_coords[3*i];
      double y2 = _local_coords[3*i+1];
      double z2 = _local_coords[3*i+2];

      double dist2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);

      if ( dist2 <= 1e-8) {
        for (int k = 0; k < stride; k++) {
          _uvw_grid[i * stride + k] = _local_to_uvw[j * stride + k];
        }
        is_find = 1;
        break;
      }
    }
    if (!is_find) {
      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_nodal_ho_ordering_from_ref_elt_set : the node (%12.5e, %12.5e, %12.5e) is not a node of the reference element\n"), x1, y1, z1);

    }
  }
  
  fvmc_nodal_ho_ordering_set (this_nodal,
                              t_elt,
                              n_nodes,
                              _uvw_grid);

  //
  // fvmc_nodal_ho_ref_elt_coords (t_elt, coords_ref);

  free (_uvw_grid);
  free (_local_coords);
  free (_uvw_to_local);
  free (_local_to_uvw);
  
}

/*----------------------------------------------------------------------------
 * Return order
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *
 * return:
 *   order
 *
 *----------------------------------------------------------------------------*/

int 
fvmc_nodal_order_get (const fvmc_nodal_t  *this_nodal)
{
  return this_nodal->order;
}


/*----------------------------------------------------------------------------
 * Return maximum number of nodes in an element
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *
 * return:
 *   max_stride
 *
 *----------------------------------------------------------------------------*/

int 
fvmc_nodal_max_n_node_elt (const fvmc_nodal_t  *this_nodal)
{
  int max_stride = 0;

  for (int i = 0; i < this_nodal->n_sections; i++) {
    fvmc_nodal_section_t  *_section = this_nodal->sections[i];
    max_stride = FVMC_MAX (max_stride, _section->stride);
  }
    
  return max_stride;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
