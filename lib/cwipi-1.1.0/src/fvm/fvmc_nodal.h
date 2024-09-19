#ifndef __FVMC_NODAL_H__
#define __FVMC_NODAL_H__

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_io_num.h"

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
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

typedef struct _fvmc_nodal_t fvmc_nodal_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/
int
fvmc_nodal_n_vertices_element (fvmc_element_t type, int order);

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh representation structure.
 *
 * parameters:
 *   name  <-- name that should be assigned to the nodal mesh
 *   dim   <-- spatial dimension
 *   order <-- order
 *
 * returns:
 *  pointer to created nodal mesh representation structure
 *----------------------------------------------------------------------------*/

fvmc_nodal_t *
fvmc_nodal_create(const char  *name,
                  int          dim,
                  int          order);

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
fvmc_nodal_destroy(fvmc_nodal_t  *this_nodal);

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
fvmc_nodal_copy(const fvmc_nodal_t *this_nodal);

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
                 int           del_vertex_num);

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
fvmc_nodal_change_parent_num(fvmc_nodal_t        *this_nodal,
                            const fvmc_lnum_t    new_parent_num[],
                            int                 entity_dim);

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
                            int           entity_dim);

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
fvmc_nodal_init_io_num(fvmc_nodal_t        *this_nodal,
                      const fvmc_gnum_t    parent_global_numbers[],
                      int                 entity_dim);

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
                             fvmc_lnum_t    parent_vertex_num[]);

/*----------------------------------------------------------------------------
 * Assign shared vertex coordinates to an extracted nodal mesh,
 * renumbering vertex numbers based on those really referenced,
 * and updating connectivity arrays in accordance.
 *
 * This function should be called once all element sections have
 * been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   vertex_coords   <-- coordinates of parent vertices (interlaced)
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_set_shared_vertices(fvmc_nodal_t        *this_nodal,
                              const fvmc_coord_t   vertex_coords[]);

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
                            fvmc_coord_t   vertex_coords[]);

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
fvmc_nodal_get_name(const fvmc_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return spatial dimension of the nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *
 * returns:
 *  spatial dimension
 *----------------------------------------------------------------------------*/

int
fvmc_nodal_get_dim(const fvmc_nodal_t  *this_nodal);

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
fvmc_nodal_get_max_entity_dim(const fvmc_nodal_t  *this_nodal);

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
                         int                 entity_dim);

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
fvmc_nodal_get_n_g_vertices(const fvmc_nodal_t  *this_nodal);

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
                           fvmc_element_t       element_type);

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
                         fvmc_element_t       element_type);

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
                         fvmc_lnum_t          parent_num[]);

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
                    fvmc_lnum_t     *error_count);

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
                     const fvmc_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure.
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_dump(const fvmc_nodal_t  *this_nodal);


/* Ajout temporaire....*/

void
fvmc_nodal_get_vertex(const fvmc_nodal_t  *this_nodal,
                     fvmc_lnum_t *n_elts,
                     fvmc_lnum_t **vertices_index, 
                     fvmc_lnum_t **vertices);

void
fvmc_nodal_get_coords(const fvmc_nodal_t  *this_nodal,
                     fvmc_lnum_t *n_vertex,
                     double **coords); 
void
fvmc_nodal_get_poly_vertex(const fvmc_nodal_t  *this_nodal,
                          fvmc_lnum_t *n_elts,
                          fvmc_lnum_t **faces,  
                          fvmc_lnum_t **faces_index, 
                          fvmc_lnum_t **vertices, 
                          fvmc_lnum_t **vertices_index);

/*----------------------------------------------------------------------------
 * Set high order ordering
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *   t_elt      <-- type of element
 *   n_nodes    <-- number of nodes
 *   ordering   <-- ordering
 *
 *----------------------------------------------------------------------------*/

void
fvmc_nodal_ho_ordering_set (fvmc_nodal_t  *this_nodal,
                            const fvmc_element_t t_elt,
                            const int n_nodes,
                            const int *uvw_grid);

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
                                         const double *coords);

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
fvmc_nodal_order_get (const fvmc_nodal_t  *this_nodal);

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
fvmc_nodal_get_n_node_elt(const fvmc_nodal_t  *this_nodal, const int elt);

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
fvmc_nodal_get_type_elt(const fvmc_nodal_t  *this_nodal, const int elt);

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
fvmc_nodal_get_local_to_user_numbering_elt (const fvmc_nodal_t  *this_nodal, const int elt);

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
fvmc_nodal_get_internal_connec_elt(const fvmc_nodal_t  *this_nodal, const int elt);

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
fvmc_nodal_get_connec_elt(const fvmc_nodal_t  *this_nodal, const int elt);



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
fvmc_nodal_max_n_node_elt (const fvmc_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_NODAL_H__ */
