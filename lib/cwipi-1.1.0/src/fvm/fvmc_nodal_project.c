/*============================================================================
 * Triangulation of nodal mesh sections
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
#include <math.h>
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
#include "fvmc_io_num.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_nodal_project.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal_triangulate.h"

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
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new mesh section from projection of a given extruded section
 * to its base plane.
 *
 * The algorithm used here always leads to one edge per face.
 * The global numbering of the new section is not generated here, as it may
 * simply be transferred from the initial section.
 *
 * parameters:
 *   dim               <-- spatial dimension
 *   vertex_coords     <-- associated vertex coordinates array
 *   parent_vertex_num <-- optional indirection to vertex coordinates
 *   base_section      <-- pointer to structure that should be triangulated
 *   error_count       --> number of triangulation errors counter (optional)
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static fvmc_nodal_section_t *
_faces_to_edges(int                         dim,
                int                         chosen_axis,
                const fvmc_coord_t           vertex_coords[],
                const fvmc_lnum_t            parent_vertex_num[],
                const fvmc_nodal_section_t  *base_section,
                _Bool                      *selected_vertices,
                fvmc_lnum_t                 *error_count)
{
  fvmc_lnum_t n_vertices, n_elements;
  fvmc_lnum_t i, j, k, vertex_id;

  fvmc_nodal_section_t *ret_section = NULL;

  /* Initialization */

  if (error_count != NULL)
    *error_count = 0;

  n_elements = base_section->n_elements;

  if (base_section->order > 1) {
    bftc_error(__FILE__, __LINE__, 0, "element order > 1 is not tanking into account\n");
  }
  
  /* Create new section */
  /*--------------------*/

  int order = 1;
  
  ret_section = fvmc_nodal_section_create(FVMC_EDGE, order);

  ret_section->n_elements = base_section->n_elements;
  ret_section->stride = 2;
  ret_section->connectivity_size =   ret_section->stride
                                   * ret_section->n_elements;

  BFTC_MALLOC(ret_section->_vertex_num,
             ret_section->connectivity_size,
             fvmc_lnum_t);
  ret_section->vertex_num = ret_section->_vertex_num;

  if (base_section->parent_element_num != NULL) {

    BFTC_MALLOC(ret_section->_parent_element_num,
               ret_section->n_elements,
               fvmc_lnum_t);
    ret_section->parent_element_num = ret_section->_parent_element_num;

  }

  /* Main loop on section face elements */
  /*------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    fvmc_lnum_t tmp_id[2], tmp_selected_edge[2];
    fvmc_coord_t tmp_edge_center;

    fvmc_lnum_t selected_edge[2] = {0, 0};
    fvmc_coord_t selected_edge_center = 0;

    /* Compute number of vertices on each face */

    if (base_section->vertex_index != NULL) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      vertex_id = base_section->vertex_index[i];
    }
    else {
      n_vertices = base_section->stride;
      vertex_id = base_section->stride * i;
    }

    /* Initialize the edge selection */

    selected_edge[0] = base_section->vertex_num[vertex_id + n_vertices - 1];
    selected_edge[1] = base_section->vertex_num[vertex_id];

    if (parent_vertex_num == NULL)
      selected_edge_center = 0.5 *
        ( vertex_coords[(selected_edge[0]-1)*dim + chosen_axis]
        + vertex_coords[(selected_edge[1]-1)*dim + chosen_axis] );

    else {

      tmp_id[0] = parent_vertex_num[selected_edge[0] - 1];
      tmp_id[1] = parent_vertex_num[selected_edge[1] - 1];

      selected_edge_center = 0.5 *
        ( vertex_coords[(tmp_id[0]-1) * dim + chosen_axis]
        + vertex_coords[(tmp_id[1]-1) * dim + chosen_axis] );

    }

    /* Loop on face's vertices */

    for (k = 1; k < n_vertices ; k++) {

      tmp_selected_edge[0] = base_section->vertex_num[vertex_id + k - 1];
      tmp_selected_edge[1] = base_section->vertex_num[vertex_id + k];

      if (parent_vertex_num == NULL)
        tmp_edge_center = 0.5 *
          (  vertex_coords[(tmp_selected_edge[0]-1)*dim + chosen_axis]
           + vertex_coords[(tmp_selected_edge[1]-1)*dim + chosen_axis]);

      else {

        tmp_id[0] = parent_vertex_num[tmp_selected_edge[0] - 1];
        tmp_id[1] = parent_vertex_num[tmp_selected_edge[1] - 1];

        tmp_edge_center
          = 0.5 * (  vertex_coords[(tmp_id[0]-1)*dim + chosen_axis]
                   + vertex_coords[(tmp_id[1]-1)*dim + chosen_axis]);

      }

      if (tmp_edge_center < selected_edge_center) {

        selected_edge_center = tmp_edge_center;

        selected_edge[0] = tmp_selected_edge[0];
        selected_edge[1] = tmp_selected_edge[1];

      }

    } /* End of loop on element's vertices */

    selected_vertices[selected_edge[0]-1] = true;
    selected_vertices[selected_edge[1]-1] = true;

    /* Define vertex_num */

    for (j = 0; j < 2; j++)
      ret_section->_vertex_num[i*2 + j] = selected_edge[j];

    if (base_section->parent_element_num != NULL)
      ret_section->_parent_element_num[i]
        = base_section->parent_element_num[i];

  } /* End of loop on elements of the section */

  /* Return new section */

  return ret_section;
}

/*----------------------------------------------------------------------------
 * Delete unused vertices after extracting edges from faces.
 *
 * parameters:
 *   this_nodal        <-> nodal mesh structure
 *   selected_vertices <-- number of triangulation errors counter (optional)
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static void
_compact_mesh(fvmc_nodal_t   *this_nodal,
              _Bool         *selected_vertices)
{
  int i, j, section_id;

  int idx = 0;

  fvmc_lnum_t   *old_to_new = NULL;
  fvmc_lnum_t   *new_to_old = NULL;
  fvmc_lnum_t   *new_parent_vtx_num = NULL;
  fvmc_coord_t  *new_coords = NULL;
  fvmc_io_num_t *new_vtx_io_num = NULL;

  fvmc_lnum_t n_vertices_new = 0;
  fvmc_lnum_t n_vertices_old = this_nodal->n_vertices;

  const int dim = this_nodal->dim;
  const fvmc_gnum_t *old_global_vtx_num = NULL;

  /* Count the new number of vertices */

  for (i = 0; i < n_vertices_old; i++) {
    if (selected_vertices[i] == true)
      idx++;
  }

  n_vertices_new = idx;


  BFTC_MALLOC(new_to_old, n_vertices_new, fvmc_lnum_t);
  BFTC_MALLOC(old_to_new, n_vertices_old, fvmc_lnum_t);

  for (i = 0, idx = 0; i < n_vertices_old; i++) {

    old_to_new[i] = -1;

    if (selected_vertices[i] == true) {
      new_to_old[idx] = i + 1;
      old_to_new[i] = idx + 1;
      idx++;
    }

  }

  if (n_vertices_new != n_vertices_old) {

    /* Adapt nodal structure if required */
    /* --------------------------------- */

    if (this_nodal->_vertex_coords != NULL) {

      /* generate a new coordinate array */

      BFTC_MALLOC(new_coords, dim * n_vertices_new, fvmc_coord_t);

      if (this_nodal->_parent_vertex_num != NULL) {
        BFTC_FREE(this_nodal->_parent_vertex_num);
        this_nodal->parent_vertex_num = NULL;
      }

      for (i = 0, idx = 0; i < n_vertices_old; i++) {

        if (selected_vertices[i] == true) {
          for (j = 0; j < dim; j++)
            new_coords[idx*dim + j] = this_nodal->vertex_coords[i*dim + j];
          idx++;
        }

      }

    } /* End if _coords != NULL */

    else { /* this_nodal->_coords == NULL */

      if (this_nodal->parent_vertex_num != NULL) {

        /* generate a new parent_vertex_num */

        BFTC_MALLOC(new_parent_vtx_num, n_vertices_new, fvmc_lnum_t);

        for (i = 0, idx = 0; i < n_vertices_old; i++) {

          if (selected_vertices[i] == true)
            new_parent_vtx_num[idx++] = this_nodal->parent_vertex_num[i];

        }

        if (this_nodal->_parent_vertex_num != NULL)
          BFTC_FREE(this_nodal->_parent_vertex_num);

        this_nodal->_parent_vertex_num = new_parent_vtx_num;
        this_nodal->parent_vertex_num = this_nodal->_parent_vertex_num;


      } /* End if parent_vertex_num != NULL */

    } /* End if this_nodal->_coords == NULL */

    /* Adapt _vertex_num in each section */
    /* --------------------------------- */

    for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

      fvmc_nodal_section_t *section = this_nodal->sections[section_id];
      const fvmc_lnum_t n_connect = section->stride * section->n_elements;

      if (section->type == FVMC_EDGE) {

        if (section->_vertex_num == NULL)
          BFTC_MALLOC(section->_vertex_num, n_connect, fvmc_lnum_t);

        for (j = 0; j < n_connect; j++)
          section->_vertex_num[j] = old_to_new[section->vertex_num[j] - 1];

        section->vertex_num = section->_vertex_num;
      }

    }

  } /* End if n_vertices_old != n_vertices_new */

  /* Adapt global_vertex_num if required */
  /* ----------------------------------- */

  if (this_nodal->global_vertex_num != NULL) {

    old_global_vtx_num = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);

    new_vtx_io_num = fvmc_io_num_create(new_to_old,
                                       old_global_vtx_num,
                                       n_vertices_new,
                                       0);

    fvmc_io_num_destroy(this_nodal->global_vertex_num);

  }

  this_nodal->global_vertex_num = new_vtx_io_num;
  this_nodal->n_vertices = n_vertices_new;

  BFTC_FREE(old_to_new);
  BFTC_FREE(new_to_old);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Project an extruded mesh to its base plane.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be cut in edges
 *   chosen_axis       <-- indicate which axis is selected to extract edges
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_project(fvmc_nodal_t  *this_nodal,
                  int           chosen_axis,
                  fvmc_lnum_t   *error_count)
{
  int i;

  fvmc_lnum_t n_vertices = 0;
  fvmc_lnum_t n_edges = 0;
  fvmc_lnum_t section_error_count = 0;

  _Bool *selected_vertices = NULL;

  assert(this_nodal != NULL);

  n_vertices = this_nodal->n_vertices;
  BFTC_MALLOC(selected_vertices, n_vertices, _Bool);

  for (i = 0; i < n_vertices; i++)
    selected_vertices[i] = false;

  /* Now extract edges and update new section list */

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvmc_nodal_section_t *e_section = NULL;
    fvmc_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->entity_dim == 2) {

      e_section = _faces_to_edges(this_nodal->dim,
                                  chosen_axis,
                                  this_nodal->vertex_coords,
                                  this_nodal->parent_vertex_num,
                                  _section,
                                  selected_vertices,
                                  &section_error_count);

      if (error_count != NULL)
        *error_count += section_error_count;

      /* Transfer global numbering if present */

      if (_section->global_element_num != NULL) {
        e_section->global_element_num = _section->global_element_num;
        _section->global_element_num = NULL;
      }

      fvmc_nodal_section_destroy(_section);
      this_nodal->sections[i] = e_section;

      n_edges += e_section->n_elements;

    }

  }

  /* Adapt parent_vertex_num, coords and global_vertex_num */

  _compact_mesh(this_nodal,
                selected_vertices);

  this_nodal->n_faces = 0;
  this_nodal->n_edges = n_edges;

  BFTC_FREE(selected_vertices);
}

/*----------------------------------------------------------------------------
 * Reduce the spatial dimension of a mesh, discarding the last coordinate
 * component.
 *
 * The mesh's spatial dimension is reduced by 1.
 *
 * parameters:
 *   this_nodal <-> pointer to structure that projected
 *   matrix     <-- projection matrix
 *                  3D -> 2D: (a11, a12, a13, a21, a22, a23)
 *                  2D -> 1D: (a11, a12)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_project_coords(fvmc_nodal_t  *this_nodal,
                         double       matrix[])
{
  int i;

  int new_dim = 0, old_dim, entity_dim = 0;
  fvmc_lnum_t n_vertices = 0;
  fvmc_coord_t *new_coords = NULL;

  assert(this_nodal != NULL);

  n_vertices = this_nodal->n_vertices;
  new_dim = this_nodal->dim - 1;
  old_dim = this_nodal->dim;

  entity_dim = fvmc_nodal_get_max_entity_dim(this_nodal);

  if (entity_dim > new_dim)
    bftc_error(__FILE__, __LINE__, 0,
              _("Projecting coordinates is not allowed for a mesh\n"
                "containing entities of dimension %d, as its\n"
                "spatial dimension would be reduced to %d"),
              entity_dim, new_dim);

  BFTC_MALLOC(new_coords, n_vertices*new_dim, fvmc_coord_t);

  if (old_dim == 3) {

    if (this_nodal->parent_vertex_num  == NULL) {
      for (i = 0; i < n_vertices; i++) {
        const double *c = this_nodal->vertex_coords + i*3;
        new_coords[i*2]     = matrix[0]*c[0] + matrix[1]*c[1] + matrix[2]*c[2];
        new_coords[i*2 + 1] = matrix[3]*c[0] + matrix[4]*c[1] + matrix[5]*c[2];
      }
    }
    else {
      for (i = 0; i < n_vertices; i++) {
        const double *c =   this_nodal->vertex_coords
                          + (this_nodal->parent_vertex_num[i]-1)*3;
        new_coords[i*2]     = matrix[0]*c[0] + matrix[1]*c[1] + matrix[2]*c[2];
        new_coords[i*2 + 1] = matrix[3]*c[0] + matrix[4]*c[1] + matrix[5]*c[2];
      }
    }

  }

  else if (old_dim == 2) {

    if (this_nodal->parent_vertex_num  == NULL) {
      for (i = 0; i < n_vertices; i++) {
        const double *c = this_nodal->vertex_coords + i*2;
        new_coords[i] = matrix[0]*c[0] + matrix[1]*c[1];
      }
    }
    else {
      for (i = 0; i < n_vertices; i++) {
        const double *c =   this_nodal->vertex_coords
                          + (this_nodal->parent_vertex_num[i]-1)*2;
        new_coords[i] = matrix[0]*c[0] + matrix[1]*c[1];
      }
    }

  }

  else
    bftc_error(__FILE__, __LINE__, 0,
              _("Projecting coordinates is only allowed for a mesh\n"
                "of initial spatial dimension %d"),
              old_dim);

  /* Update structure */

  this_nodal->dim = new_dim;

  if (this_nodal->_vertex_coords != NULL)
    BFTC_FREE(this_nodal->_vertex_coords);

  this_nodal->parent_vertex_num = NULL;
  if (this_nodal->_parent_vertex_num != NULL)
    BFTC_FREE(this_nodal->_parent_vertex_num);

  this_nodal->vertex_coords = new_coords;
  this_nodal->_vertex_coords = new_coords;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
