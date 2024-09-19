#ifndef __FVMC_NODAL_PRIV_H__
#define __FVMC_NODAL_PRIV_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_tesselation.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a mesh section
 *----------------------------------------------------------------------------*/

typedef struct _fvmc_nodal_section_t {

  /* Basic information */
  /*-------------------*/

  int         entity_dim;          /* Entity dimension */

  fvmc_lnum_t  n_elements;          /* Number of elements */

  fvmc_element_t  type;             /* Element types */
  int             order;            /* Element order */
  int             *ho_local_to_user_ordering;     /* Local element ordering (local to user ordering) */
  
  /* Connectivity */
  /*--------------*/

  size_t      connectivity_size;   /* Size of vertex_num array;
                                      for strided elements:
                                       (n_elements * stride)
                                      for polygons:
                                       (vertex_index[n_elements])
                                      for polyhedra:
                                       (vertex_index[n_faces]) */

  int         stride;              /* Element size for regular elements
                                      (0 for polygons and polyhedra) */

  fvmc_lnum_t  n_faces;             /* Number of faces defining polyhedra */

  /* Pointers to connectivity arrays, which may be shared */

  const fvmc_lnum_t  *face_index;   /* polyhedron -> faces index (O to n-1);
                                      size: n_elements + 1 */
  const fvmc_lnum_t  *face_num;     /* polyhedron -> face numbers (1 to n, signed,
                                      > 0 for outwards pointing face normal
                                      < 0 for inwards pointing face normal);
                                      size: face_index[n_elements] */

  const fvmc_lnum_t  *vertex_index; /* polygon face -> vertices index (O to n-1);
                                      size: n_faces + 1 */

  const fvmc_lnum_t  *vertex_num;   /* vertex numbers (1 to n);
                                      size: connectivity_size */

  /* Pointers to local connectivity arrays, if owner */

  fvmc_lnum_t  *_face_index;        /* face_index if owner, NULL if shared */
  fvmc_lnum_t  *_face_num;          /* face_num if owner, NULL if shared */
  fvmc_lnum_t  *_vertex_index;      /* vertex_index if owner, NULL if shared */
  fvmc_lnum_t  *_vertex_num;        /* vertex numbers if owner, NULL if shared */

  fvmc_lnum_t  *_ho_vertex_num;     /* vertex numbers odered from _ho_ordering */

  /* Auxiliary structure used to define subdivision of elements into
     simpler element types (usually polygons to triangles and
     polyhedra to tetrahedra and pyramids) */

  fvmc_tesselation_t  *tesselation;

  /* Numbering */
  /*-----------*/

  const fvmc_lnum_t  *parent_element_num; /* Local numbers (1 to n) of local
                                            elements in the parent mesh,
                                            associated with the section's
                                            elements.

                                            This array is necessary to redis-
                                            tribute output fields when the
                                            section has been either associated
                                            with an unsorted mixed mesh,
                                            renumbered, or is associated with a
                                            subset of a more complete mesh,
                                            such as a clip plane. When used for
                                            a subset, it also defines the lists
                                            of elements of the parent mesh
                                            belonging to that subset.

                                            This array is present only when non
                                            "trivial" (i.e. not 1, 2, ..., n). */

  fvmc_lnum_t    *_parent_element_num;    /* pointer to parent_element_num if
                                            owner, NULL otherwise */

  fvmc_io_num_t  *global_element_num;     /* Global element numbers */

} fvmc_nodal_section_t;

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

struct _fvmc_nodal_t {

  /* Global indicators */
  /*-------------------*/

  char  *name;                 /* Mesh name */

  int    dim;                  /* Spatial dimension */
  int    num_dom;              /* Local domain number */
  int    n_doms;               /* Global number of domains */
  int    n_sections;           /* Number of sections */

  int    order;                /* order */
  int   **ho_uvw_to_local_ordering;  /* (U, V, W) ordering to local element ordering 
                                        for each reference element type */
  int   **ho_user_to_uvw;            /* user ordering to (U, V, W) */

  /* Local dimensions */
  /*------------------*/

  /* Total number of cells, faces, edges, and vertices */
  
  fvmc_lnum_t  n_cells;
  fvmc_lnum_t  n_faces;
  fvmc_lnum_t  n_edges;
  fvmc_lnum_t  n_vertices;

  fvmc_lnum_t *sections_idx;
  
  /* Vertex definitions; */
  /*---------------------*/

  const fvmc_coord_t  *vertex_coords;    /* pointer to  vertex coordinates
                                           (always interlaced:
                                           x1, y1, z1, x2, y2, z2, ...) */
  fvmc_coord_t        *_vertex_coords;   /* pointer to vertex coordinates if
                                           owner (for use with own algorithms) */

  const fvmc_lnum_t  *parent_vertex_num; /* Local numbers (1 to n) of local
                                           vertices in the parent mesh.

                                           This array is necessary to redis-
                                           tribute output fields when a nodal
                                           mesh has been renumbered or is
                                           associated with a subset of a more
                                           complete mesh, such as a clip plane
                                           (in which case it also defines the
                                           lists of vertices of the parent
                                           mesh in that subset).

                                           This array is present only when non
                                           "trivial" (i.e. not 1, 2, ..., n). */

  fvmc_lnum_t    *_parent_vertex_num;    /* pointer to parent_vertex_num if
                                           owner, NULL otherwise */

  fvmc_io_num_t  *global_vertex_num;     /* Global vertex numbering */

  /* Mesh connectivity */
  /*-------------------*/

  fvmc_nodal_section_t  **sections;  /* Array of section descriptions */

};

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh section representation structure.
 *
 * parameters:
 *   type <-- type of element defined by this section
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_section_t *
fvmc_nodal_section_create(const fvmc_element_t  type, int order);

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh section representation structure.
 *
 * parameters:
 *   this_section <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_nodal_section_t *
fvmc_nodal_section_destroy(fvmc_nodal_section_t  * this_section);

/*----------------------------------------------------------------------------
 * Copy selected shared connectivity information to private connectivity
 * for a nodal mesh section .
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
                                _Bool                 copy_vertex_num);

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
fvmc_nodal_section_n_g_elements(const fvmc_nodal_section_t  *this_section);

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
fvmc_nodal_n_g_vertices(const fvmc_nodal_t  *this_nodal);

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
                            int             face_vertices[6][4]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_NODAL_PRIV_H__ */
