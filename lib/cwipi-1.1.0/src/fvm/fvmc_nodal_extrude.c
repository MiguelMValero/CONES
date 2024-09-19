/*============================================================================
 * Extrusion of a nodal representation associated with a mesh
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
#include "fvmc_io_num.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_nodal_extrude.h"

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
 * Extrude strided section.
 *
 * Entity parent numbering is removed.
 *
 * parameters:
 *   this_section      <-> pointer to structure that should be extruded
 *----------------------------------------------------------------------------*/

static void
_extrude_strided_section(fvmc_nodal_section_t  * this_section,
                         const fvmc_lnum_t       n_layers)
{
  int stride;
  size_t  connectivity_size;
  int k;
  fvmc_lnum_t  i, j;
  fvmc_lnum_t  base_vertex_id;
  fvmc_lnum_t  element_shift, layer_shift, bottom_shift, top_shift;
  fvmc_lnum_t *vertex_num;

  fvmc_lnum_t n_elements = this_section->n_elements;
  const fvmc_lnum_t  n_planes = n_layers + 1;

  /* Build new connectivity */

  stride = this_section->stride * 2;
  connectivity_size = this_section->n_elements * stride * n_layers;

  BFTC_MALLOC(vertex_num, connectivity_size, fvmc_lnum_t);
  this_section->connectivity_size = 0;

  for (i = 0; i < this_section->n_elements; i++) {
    element_shift = n_layers * stride * i;
    for (j = 0; j < n_layers; j++) {
      layer_shift = j * stride;
      bottom_shift = element_shift + layer_shift;
      top_shift = element_shift + layer_shift + this_section->stride;
      for (k = 0; k < this_section->stride; k++) {
        base_vertex_id
          = this_section->vertex_num[this_section->stride*i + k] - 1;
        vertex_num[bottom_shift + k]
          = n_planes*base_vertex_id + j + 1;
        vertex_num[top_shift + k]
          = n_planes*base_vertex_id + j + 2;
      }
    }
  }

  this_section->connectivity_size = connectivity_size;
  /* Replace old connectivity */

  if (this_section->_vertex_num != NULL)
    BFTC_FREE(this_section->_vertex_num);

  this_section->_vertex_num = vertex_num;
  this_section->vertex_num = this_section->_vertex_num;

  this_section->connectivity_size = connectivity_size;

  /* Remove old parent numbering */

  this_section->parent_element_num = NULL;
  if (this_section->_parent_element_num != NULL)
    BFTC_FREE(this_section->_parent_element_num);

  /* Update global_numbering */

  if (this_section->global_element_num != NULL) {

    /* Create new global numbering */

    fvmc_gnum_t  *global_element_num = NULL;

    const fvmc_gnum_t *old_global_element_num
      = fvmc_io_num_get_global_num(this_section->global_element_num);

    BFTC_MALLOC(global_element_num, n_elements*n_layers, fvmc_gnum_t);

    for (i = 0; i < n_elements; i++) {
      fvmc_gnum_t  base_num = (  (old_global_element_num[i]-1)
                              * (fvmc_gnum_t)n_layers) + 1;
      for (j = 0; j < n_layers; j++)
        global_element_num[i*n_layers + j] = base_num + (fvmc_gnum_t)j;
    }

    /* Remplace old global number with new */

    fvmc_io_num_destroy(this_section->global_element_num);
    this_section->global_element_num = fvmc_io_num_create(NULL,
                                                         global_element_num,
                                                         n_elements * n_layers,
                                                         0);

  }

  /* Update section info */

  this_section->n_elements *= n_layers;

  switch(this_section->type) {
  case FVMC_EDGE:
    this_section->type = FVMC_FACE_QUAD;
    break;
  case FVMC_FACE_TRIA:
    this_section->type = FVMC_CELL_PRISM;
    break;
  case FVMC_FACE_QUAD:
    this_section->type = FVMC_CELL_HEXA;
    break;
  default:
    assert(0);
  }

  this_section->entity_dim += 1;
  this_section->stride *=2;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extrude nodal mesh.
 *
 * Vertex and element parent numbering is removed if present.
 *
 * Note: layout of new elements in memory is such that the definitions
 *       of all elements extruded from a same ancestor are contiguous.
 *       that is, {e_1, e_2, ..., e_n} leads to
 *       {e_1_layer_1, ..., e_1_layer_m, e_2_layer_1, ... e_n_layer_m}
 *
 * parameters:
 *   this_section      <-> pointer to structure that should be extruded
 *   n_layers          <-> number of extruded layers
 *   extrusion_vectors <-> length and direction of extrusion for each vertex;
 *                         size: mesh_spatial_dim . n_vertices
 *   distribution      <-> optional distribution of resulting vertices
 *                         along each extrusion vector (size: n_layers + 1)
 *                         with values ranging from 0 to 1, or NULL for
 *                         a regular distribution.
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_extrude(fvmc_nodal_t        *const this_nodal,
                  const fvmc_lnum_t          n_layers,
                  const fvmc_coord_t         extrusion_vectors[],
                  const fvmc_coord_t         distribution[])
{
  int dim, k;
  fvmc_lnum_t  i, j;
  fvmc_lnum_t  n_vertices;
  fvmc_lnum_t  vertex_shift;
  fvmc_coord_t  *_distrib = NULL;
  fvmc_coord_t  *new_coords = NULL;

  const fvmc_lnum_t  n_planes = n_layers + 1;
  const fvmc_coord_t  *distrib = distribution;
  const fvmc_coord_t  *old_coords = NULL;

  assert(this_nodal != NULL);
  assert(extrusion_vectors != NULL || this_nodal->n_vertices == 0);

  dim = this_nodal->dim;

  /* Check that no section is of too high dimension */

  for (i = 0; i < this_nodal->n_sections; i++) {
    const fvmc_nodal_section_t *_section = this_nodal->sections[i];
    if (_section->entity_dim >= dim)
      bftc_error(__FILE__, __LINE__, 0,
                _("Dimension of mesh \"%s\" section %d equals %d\n"
                  "with mesh spatial dimension %d prior to extrusion\n"
                  "when it should be smaller."),
                this_nodal->name, i+1, _section->entity_dim, dim);
  }

  /* Set distribution if necessary */

  if (distribution == NULL) {
    BFTC_MALLOC(_distrib, n_planes, fvmc_coord_t);
    for (i = 0; i < n_planes; i++)
      _distrib[i] = ((double)i) / ((double)n_layers);
    distrib = _distrib;
  }

  /* Compute new coordinates */

  n_vertices = this_nodal->n_vertices;
  old_coords = this_nodal->vertex_coords;

  BFTC_MALLOC(new_coords, n_planes*n_vertices*dim, fvmc_coord_t);

  if (this_nodal->_parent_vertex_num != NULL) {

    for (i = 0; i < n_vertices; i++) {
      const double *_old_coords
        = old_coords + ((this_nodal->parent_vertex_num[i]-1) * dim);
      vertex_shift = n_planes * dim * i;
      for (j = 0; j < n_planes; j++) {
        for (k = 0; k < dim; k++) {
          new_coords[vertex_shift + (j*dim) + k]
            =   _old_coords[k]
              + (extrusion_vectors[i*dim + k] * distrib[j]);
        }
      }
    }

  }
  else {

    for (i = 0; i < n_vertices; i++) {
      vertex_shift = n_planes * dim * i;
      for (j = 0; j < n_planes; j++) {
        for (k = 0; k < dim; k++) {
          new_coords[vertex_shift + (j*dim) + k]
            =   old_coords[i*dim + k]
              + (extrusion_vectors[i*dim + k] * distrib[j]);
        }
      }
    }

  }

  /* Replace old coords with new */

  if (this_nodal->_vertex_coords != NULL)
    BFTC_FREE(this_nodal->_vertex_coords);

  this_nodal->_vertex_coords = new_coords;
  this_nodal->vertex_coords = this_nodal->_vertex_coords;

  this_nodal->parent_vertex_num = NULL;
  if (this_nodal->_parent_vertex_num != NULL)
    BFTC_FREE(this_nodal->_parent_vertex_num);

  /* Update global numbering */

  if (this_nodal->global_vertex_num != NULL) {

    /* Create new global numbering */

    fvmc_gnum_t  *global_vertex_num = NULL;

    const fvmc_gnum_t *old_global_vertex_num
      = fvmc_io_num_get_global_num(this_nodal->global_vertex_num);

    BFTC_MALLOC(global_vertex_num, n_planes*n_vertices, fvmc_gnum_t);

    for (i = 0; i < n_vertices; i++) {
      fvmc_gnum_t  base_num = (  (old_global_vertex_num[i]-1)
                              * (fvmc_gnum_t)n_planes) + 1;
      for (j = 0; j < n_planes; j++)
        global_vertex_num[i*n_planes + j] = base_num + (fvmc_gnum_t)j;
    }

    /* Remplace old global number with new */

    fvmc_io_num_destroy(this_nodal->global_vertex_num);
    this_nodal->global_vertex_num = fvmc_io_num_create(NULL,
                                                      global_vertex_num,
                                                      n_vertices * n_planes,
                                                      0);

  }

  /* We may now update the number of vertices */

  this_nodal->n_vertices = n_vertices * n_planes;

  /* Extrude element definitions */

  this_nodal->n_cells = 0;
  this_nodal->n_faces = 0;
  this_nodal->n_edges = 0;

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvmc_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->stride > 0)
      _extrude_strided_section(_section,
                               n_layers);

    else
      bftc_error(__FILE__, __LINE__, 0,
                _("Extrusion of non strided sections not implemented yet."));

    switch(_section->entity_dim) {
    case 3:
      this_nodal->n_cells += _section->n_elements;
      break;
    case 2:
      this_nodal->n_faces += _section->n_elements;
      break;
    default:
      assert(0);
    }

  }

  /* If the mesh contains only vertices and no elements, a section
     of edges linking "extruded" vertices should be created */

  if (this_nodal->n_vertices != 0 && this_nodal->n_sections == 0)
    bftc_error(__FILE__, __LINE__, 0,
              _("Extrusion of vertices only to edges not implemented yet."));

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
