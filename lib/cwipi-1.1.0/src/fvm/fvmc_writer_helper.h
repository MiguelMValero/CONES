#ifndef __FVMC_WRITER_HELPER_H__
#define __FVMC_WRITER_HELPER_H__

/*============================================================================
 * Helper types and functions for mesh and field writers
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2006  EDF

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

#if defined(FVMC_HAVE_MPI )
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_gather.h"
#include "fvmc_nodal.h"
#include "fvmc_writer.h"

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
 * FVM nodal to writer section translation list
 *----------------------------------------------------------------------------*/

typedef struct _fvmc_writer_section_t {

  struct _fvmc_writer_section_t  *next;  /* Pointer to next element
                                           in list (NULL at end) */

  const fvmc_nodal_section_t  *section;  /* Corresponding section in mesh */

  fvmc_gnum_t  extra_vertex_base;        /* Start global number of added
                                           vertices (for tesselation) */

  fvmc_lnum_t  num_shift;                /* Element number shift when no
                                           parent lists are used */
  fvmc_element_t  type;                  /* Corresponding element type (may
                                           differ from  section->type when
                                           using tesselations) */
  int order;

  _Bool   continues_previous;           /* Indicates if the corresponding FVM
                                           nodal section should be appended
                                           to the previous one on output */

} fvmc_writer_section_t;

/*----------------------------------------------------------------------------
 * FVM nodal to writer field output helper
 *----------------------------------------------------------------------------*/

/*
  Pointer to structure keeping track of the status of a writer's field
  output state. The structure itself is private, and is defined in fvmc_writer.c
*/

typedef struct _fvmc_writer_field_helper_t fvmc_writer_field_helper_t;

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build list of sections to output
 *
 * Depending on whether multiple sections of a given element type
 * are allowed or not, sections may be ordered in different ways.
 * Discarded element types are not added to the list.
 *
 * parameters:
 *   mesh                 <-- pointer to nodal mesh structure
 *   group_same_type      <-- group sections of the same type
 *   discard_polygons     <-- ignore polygonal sections
 *   discard_polyhedra    <-- ignore polyhedral sections
 *   divide_polygons      <-- tesselate polygonal sections
 *   divide_polyhedra     <-- tesselate polyhedral sections
 *
 * returns:
 *   array of section translations (must be freed by caller),
 *   or NULL if section list is completely empty
 *----------------------------------------------------------------------------*/

fvmc_writer_section_t *
fvmc_writer_export_list(const fvmc_nodal_t  *mesh,
                       _Bool               group_same_type,
                       _Bool               discard_polygons,
                       _Bool               discard_polyhedra,
                       _Bool               divide_polygons,
                       _Bool               divide_polyhedra);

/*----------------------------------------------------------------------------
 * Create field writer helper structure.
 *
 * Local values are initialized, ang lobal values are set to zero
 * (they may be initialized by calling fvmc_writer_field_helper_init_g()).
 *
 * parameters:
 *   mesh                 <-- pointer to nodal mesh structure
 *   section_list         <-- point to export section list helper structure
 *   field_dim            <-- indicates output field dimension
 *   interlace            <-- indicates if output is interlaced
 *   location             <-- indicates if field is cell or node based
 *   datatype             <-- datatype of destination buffers
 *
 * returns:
 *   pointer to allocated and initialized field writer helper
 *----------------------------------------------------------------------------*/

fvmc_writer_field_helper_t *
fvmc_writer_field_helper_create(const fvmc_nodal_t          *mesh,
                               const fvmc_writer_section_t *section_list,
                               int                         field_dim,
                               fvmc_interlace_t             interlace,
                               fvmc_datatype_t              datatype,
                               fvmc_writer_var_loc_t        location);

/*----------------------------------------------------------------------------
 * Destroy FVM field writer helper.
 *
 * parameters:
 *   helper <-> pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_writer_field_helper_t *
fvmc_writer_field_helper_destroy(fvmc_writer_field_helper_t *helper);

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize global values for an fvmc_writer_field_helper structure.
 *
 * Internal buffers for gathering of data to rank 0 of the given
 * communicator are also allocated.
 *
 * parameters:
 *   helper        <-> pointer to structure that should be initialized
 *   section_list  <-- point to export section list helper structure
 *   mesh          <-- pointer to nodal mesh structure
 *   comm          <-- associated MPI communicator
 *
 * returns:
 *   pointer to allocated and initialized field writer helper
 *----------------------------------------------------------------------------*/

void
fvmc_writer_field_helper_init_g(fvmc_writer_field_helper_t   *helper,
                               const fvmc_writer_section_t  *section_list,
                               const fvmc_nodal_t           *mesh,
                               MPI_Comm                     comm);

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return sizes associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *   input_size               --> Total field locations in input (or NULL)
 *   output_size              --> Total field locations in output (or NULL)
 *   max_grouped_elements_out --> Max. field locations in a single group
 *                                (elements of a given type if sections are
 *                                grouped, elements of a given section
 *                                otherwise; NULL if unused)
 *   min_output_buffer_size   --> Minimum required buffer size (or NULL)
 *----------------------------------------------------------------------------*/

void
fvmc_writer_field_helper_get_size(const fvmc_writer_field_helper_t  *helper,
                                 size_t  *input_size,
                                 size_t  *output_size,
                                 size_t  *max_grouped_elements_out,
                                 size_t  *min_output_buffer_size);

/*----------------------------------------------------------------------------
 * Return the output dimension associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *
 * returns:
 *   field dimension associated with helper
 *----------------------------------------------------------------------------*/

int
fvmc_writer_field_helper_field_dim(const fvmc_writer_field_helper_t  *helper);

/*----------------------------------------------------------------------------
 * Return the output datatype associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *
 * returns:
 *   output datatype associated with helper
 *----------------------------------------------------------------------------*/

fvmc_datatype_t
fvmc_writer_field_helper_datatype(const fvmc_writer_field_helper_t  *helper);

/*----------------------------------------------------------------------------
 * Partially distribute field values to an output buffer.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   export_section     <-- pointer to section helper structure
 *   src_dim            <-- dimension of source data
 *   src_dim_shift      <-- source data dimension shift (start index)
 *   src_interlace      <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   datatype           <-- input data type
 *   field_values       <-- pointer to input array
 *   output_buffer      <-- pointer to output buffer
 *                          (working array only for ranks > 0)
 *   output_buffer_size <-- size of output buffer (in datatype units)
 *   output_size        --> size of output upon return (in datatype units)
 *
 * returns:
 *   0 if values were distributed to the output buffer, 1 if the end of the
 *   section has already been reached and no values were left to distribute.
 *----------------------------------------------------------------------------*/

int
fvmc_writer_field_helper_step_e(fvmc_writer_field_helper_t   *helper,
                               const fvmc_writer_section_t  *export_section,
                               int                          src_dim,
                               int                          src_dim_shift,
                               fvmc_interlace_t              src_interlace,
                               int                          n_parent_lists,
                               const fvmc_lnum_t             parent_num_shift[],
                               fvmc_datatype_t               datatype,
                               const void            *const field_values[],
                               void                        *output_buffer,
                               size_t                       output_buffer_size,
                               size_t                      *output_size);

/*----------------------------------------------------------------------------
 * Partially distribute per node field values to an output buffer.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   mesh               <-- pointer to associated mesh
 *   src_dim            <-- dimension of source data
 *   src_dim_shift      <-- source data dimension shift (start index)
 *   src_interlace      <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   datatype           <-- input data type
 *   field_values       <-- pointer to input array
 *   output_buffer      <-- pointer to output buffer
 *                          (working array only for ranks > 0)
 *   output_buffer_size <-- size of output buffer (in datatype units)
 *   output_size        --> size of output upon return (in datatype units)
 *
 * returns:
 *   0 if values were distributed to the output buffer, 1 if the end of the
 *   section has already been reached and no values were left to distribute.
 *----------------------------------------------------------------------------*/

int
fvmc_writer_field_helper_step_n(fvmc_writer_field_helper_t   *helper,
                               const fvmc_nodal_t           *mesh,
                               int                          src_dim,
                               int                          src_dim_shift,
                               fvmc_interlace_t              src_interlace,
                               int                          n_parent_lists,
                               const fvmc_lnum_t             parent_num_shift[],
                               fvmc_datatype_t               datatype,
                               const void            *const field_values[],
                               void                        *output_buffer,
                               size_t                       output_buffer_size,
                               size_t                      *output_size);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_WRITER_HELPER_H__ */
