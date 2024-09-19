#ifndef __FVMC_TESSELATION_H__
#define __FVMC_TESSELATION_H__

/*============================================================================
 * Structure describing a mesh section's subdivision into simpler elements
 *
 * This is mostly useful to replace polygons or polyhedra by simpler
 * elements such as triangles, tetrahedra, and prisms upon data export.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2006-2007  EDF

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

#if defined(FVMC_HAVE_MPI)
#include <mpi.h>
#endif

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

#define FVMC_TESSELATION_N_SUB_TYPES_MAX 2

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a tesselation of a mesh section.
 *----------------------------------------------------------------------------*/

/*
  Pointer to tesselation structure. The structure
  itself is private, and is defined in fvmc_tesselation.c
*/

typedef struct _fvmc_tesselation_t fvmc_tesselation_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a mesh section's subdivision into simpler elements.
 *
 * The structure contains pointers to the mesh section's connectivity,
 * (passed as arguments), which is not copied. This structure should thus
 * always be destroyed before the mesh section to which it relates.
 *
 * Unused connectivity array arguments should be set to NULL (such as
 * face_index[] and face_num[] for 2D or regular (strided) elements,
 * and vertex_index[] for strided elements.
 *
 * At this stage, the structure does not yet contain tesselation information.
 *
 * parameters:
 *   element_type       <-- type of elements considered
 *   n_elements         <-- number of elements
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          dimension [n_elements + 1]
 *   face_num           <-- element -> face numbers (1 to n, signed,
 *                           > 0 for outwards pointing face normal
 *                           < 0 for inwards pointing face normal);
 *                          dimension: [face_index[n_elements]], or NULL
 *   vertex_index       <-- element face -> vertices index (O to n-1);
 *                          dimension: [n_cell_faces + 1], [n_elements + 1],
 *                          or NULL depending on face_index and vertex_index
 *   vertex_num         <-- element -> vertex connectivity (1 to n)
 *   global_element_num <-- global element numbers (NULL in serial mode)
 *
 * returns:
 *  pointer to created mesh section tesselation structure
 *----------------------------------------------------------------------------*/

fvmc_tesselation_t *
fvmc_tesselation_create(fvmc_element_t        element_type,
                       fvmc_lnum_t           n_elements,
                       const fvmc_lnum_t     face_index[],
                       const fvmc_lnum_t     face_num[],
                       const fvmc_lnum_t     vertex_index[],
                       const fvmc_lnum_t     vertex_num[],
                       const fvmc_io_num_t  *global_element_num);

/*----------------------------------------------------------------------------
 * Destruction of a mesh section tesselation structure.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_tesselation_t *
fvmc_tesselation_destroy(fvmc_tesselation_t  * this_tesselation);

/*----------------------------------------------------------------------------
 * Tesselate a mesh section referred to by an fvmc_tesselation_t structure.
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   dim                <-- spatial dimension
 *   vertex_coords      <-- associated vertex coordinates array
 *   parent_vertex_num  <-- optional indirection to vertex coordinates
 *   error_count        --> number of elements with a tesselation error
 *                          counter (optional)
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_init(fvmc_tesselation_t  *this_tesselation,
                     int                 dim,
                     const fvmc_coord_t   vertex_coords[],
                     const fvmc_lnum_t    parent_vertex_num[],
                     fvmc_lnum_t         *error_count);

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh polygon splitting representation structure;
 * only the associations (numberings) necessary to redistribution of fields
 * for output are conserved, the full connectivity being no longer useful
 * once it has been output.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be reduced
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_reduce(fvmc_tesselation_t  * this_tesselation);

/*----------------------------------------------------------------------------
 * Return number of parent elements of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of parent elements
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_n_elements(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return global number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

fvmc_gnum_t
fvmc_tesselation_n_g_vertices_add(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return (local) number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_n_vertices_add(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return number of resulting sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of resulting sub-types of the tesselation
 *----------------------------------------------------------------------------*/

int
fvmc_tesselation_n_sub_types(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return given sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

fvmc_element_t
fvmc_tesselation_sub_type(const fvmc_tesselation_t  *this_tesselation,
                         int                       sub_type_id);

/*----------------------------------------------------------------------------
 * Return number of elements of a given sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_n_sub_elements(const fvmc_tesselation_t  *this_tesselation,
                               fvmc_element_t             sub_type);

/*----------------------------------------------------------------------------
 * Obtain the global and maximum number of elements of a given sub-type
 * of a tesselation.
 *
 * parameters:
 *   this_tesselation    <-- tesselation structure
 *   sub_type_id         <-- index of sub-type in tesselation (0 to n-1)
 *   n_sub_elements_glob --> global number of sub-elements of the given type
 *   n_sub_elements_max  --> maximum number of sub-elements per element
 *                           of the given type (for all ranks)
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_get_global_size(const fvmc_tesselation_t  *this_tesselation,
                                fvmc_element_t             sub_type,
                                fvmc_gnum_t               *n_sub_elements_glob,
                                fvmc_lnum_t               *n_sub_elements_max);

/*----------------------------------------------------------------------------
 * Return global numbering of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   pointer to global numbering of added vertices for this tesselation,
 *   or NULL if no added vertices are present.
 *----------------------------------------------------------------------------*/

const fvmc_io_num_t *
fvmc_tesselation_global_vertex_num(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Compute coordinates of added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   vertex_coords      --> coordinates of added vertices
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_vertex_coords(const fvmc_tesselation_t  *this_tesselation,
                              fvmc_coord_t               vertex_coords[]);

/*----------------------------------------------------------------------------
 * Return index of sub-elements associated with each element of a given
 * sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   index of sub-elements associated with each element (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

const fvmc_lnum_t *
fvmc_tesselation_sub_elt_index(const fvmc_tesselation_t  *this_tesselation,
                              fvmc_element_t             sub_type);

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute index values corresponding to given range of indices,
 * for an element -> sub-element value distribution.
 *
 * This index is used mainly to gather a decoded tesselation connectivity
 * or element -> sub-element data distribution, expanding the corresponding
 * data only on the given range.
 * Only the index values in the start_id to end_id range are set by
 * this function, starting with index[start_id] = 0.
 *
 * parameters:
 *   this_tesselation      <-- tesselation structure
 *   connect_type          <-- destination element type
 *   stride                <-- number of associated values per sub-element
 *   start_id              <-- start index of polyhedra subset in parent section
 *   buffer_limit          <-- maximum number of sub-elements of destination
 *                             element type allowable for vertex_num[] buffer
 *   global_num_end        <-> past the end (maximum + 1) parent element
 *                             global number (reduced on return if required
 *                             by buffer_size limits)
 *   index                 --> sub-element index
 *   comm                  <-- associated MPI communicator
 *
 * returns:
 *   polyhedron index end corresponding to decoded range
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_range_index_g(const fvmc_tesselation_t  *this_tesselation,
                              fvmc_element_t             connect_type,
                              int                       stride,
                              fvmc_lnum_t                start_id,
                              fvmc_lnum_t                buffer_limit,
                              fvmc_gnum_t               *global_num_end,
                              fvmc_lnum_t                index[],
                              MPI_Comm                  comm);

/*----------------------------------------------------------------------------
 * Decode tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers and computing unneeded element
 * connectivities when exporting data in slices, this function may decode
 * a partial connectivity range, starting at polygon index start_id and ending
 * either when the indicated buffer size is attained, or the global element
 * number corresponding to a given polygon exceeds a given value.
 * It returns the effective polygon index end.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of polygons subset in parent section
 *   buffer_limit       <-- maximum number of sub-elements of destination
 *                          element type allowable for vertex_num[] buffer
 *   global_num_end     <-> past the end (maximum + 1) parent element
 *                          global number (reduced on return if required
 *                          by buffer_limit)
 *   extra_vertex_base  <-- starting number for added vertices
 *   global_vertex_num  <-- global vertex numbering
 *   vertex_num         --> sub-element (global) vertex connectivity
 *   comm               <-- associated MPI communicator
 *
 * returns:
 *   polygon index corresponding to end of decoded range
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_decode_g(const fvmc_tesselation_t  *this_tesselation,
                         fvmc_element_t             connect_type,
                         fvmc_lnum_t                start_id,
                         fvmc_lnum_t                buffer_limit,
                         fvmc_gnum_t               *global_num_end,
                         const fvmc_io_num_t       *global_vertex_num,
                         fvmc_gnum_t                extra_vertex_base,
                         fvmc_gnum_t                vertex_num[],
                         MPI_Comm                  comm);

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Decode tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers and computing unneeded element
 * connectivities, this function may decode a partial connectivity range,
 * starting at polygon index start_id and ending either when the indicated
 * buffer size or the last polygon is attained.
 * It returns the effective polygon index end.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of polygons subset in parent section
 *   buffer_limit       <-- maximum number of sub-elements of destination
 *                          element type allowable for vertex_num[] buffer
 *   extra_vertex_base  <-- starting number for added vertices
 *   vertex_num         --> sub-element (global) vertex connectivity
 *
 * returns:
 *   polygon index corresponding to end of decoded range
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_tesselation_decode(const fvmc_tesselation_t  *this_tesselation,
                       fvmc_element_t             connect_type,
                       fvmc_lnum_t                start_id,
                       fvmc_lnum_t                buffer_limit,
                       fvmc_lnum_t                extra_vertex_base,
                       fvmc_lnum_t                vertex_num[]);

/*----------------------------------------------------------------------------
 * Distribute "per element" data from the base elements to their tesselation.
 *
 * The same data array is used for input and output, so as to avoid requiring
 * excess allocation in typical use cases (extracting data from a parent mesh
 * to a buffer and distributing it as per its tesselation).
 * The data array should be at least of size:
 * [sub_elt_index[end_id] - sub_elt_index[start_id] * size
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of elements subset in parent section
 *   end_id             <-- end index of elements subset in parent section
 *   size               <-- data size for each element (sizeof(type)*stride)
 *   data               <-> undistributed data in, distributed data out
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_distribute(const fvmc_tesselation_t  *this_tesselation,
                           fvmc_element_t             connect_type,
                           fvmc_lnum_t                start_id,
                           fvmc_lnum_t                end_id,
                           size_t                    size,
                           void                     *data);

/*----------------------------------------------------------------------------
 * Compute field values at added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   vertex_coords    <-- coordinates of added vertices
 *   src_dim          <-- dimension of source data
 *   src_dim_shift    <-- source data dimension shift (start index)
 *   dest_dim         <-- destination data dimension (1 if non interlaced)
 *   start_id         <-- added vertices start index
 *   end_id           <-- added vertices past the end index
 *   src_interlace    <-- indicates if source data is interlaced
 *   src_datatype     <-- source data type (float, double, or int)
 *   dest_datatype    <-- destination data type (float, double, or int)
 *   n_parent_lists   <-- number of parent lists (if parent_num != NULL)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   parent_num       <-- if n_parent_lists > 0, parent entity numbers
 *   src_data         <-- array of source arrays (at least one, with one per
 *                        source dimension if non interlaced, times one per
 *                        parent list if multiple parent lists, with
 *                        x_parent_1, y_parent_1, ..., x_parent_2, ...) order
 *   dest_data        --> destination buffer
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_vertex_values(const fvmc_tesselation_t  *this_tesselation,
                              int                       src_dim,
                              int                       src_dim_shift,
                              int                       dest_dim,
                              fvmc_lnum_t                start_id,
                              fvmc_lnum_t                end_id,
                              fvmc_interlace_t           src_interlace,
                              fvmc_datatype_t            src_datatype,
                              fvmc_datatype_t            dest_datatype,
                              int                       n_parent_lists,
                              const fvmc_lnum_t          parent_num_shift[],
                              const fvmc_lnum_t          parent_num[],
                              const void         *const src_data[],
                              void               *const dest_data);

/*----------------------------------------------------------------------------
 * Dump printout of a mesh section tesselation structure.
 *
 * parameters:
 *   this_tesselation  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvmc_tesselation_dump(const fvmc_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_TESSELATION_H__ */
