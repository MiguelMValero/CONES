#ifndef __FVMC_TO_CGNS_H__
#define __FVMC_TO_CGNS_H__

#if defined(HAVE_CGNS)

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CGNS files
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2005-2007  EDF

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to CGNS file writer.
 *
 * Options are:
 *   discard_polygons    do not output polygons or related values
 *   discard_polyhedra   do not output polyhedra or related values
 *   divide_polygons     tesselate polygons with triangles
 *
 * As CGNS does not handle polyhedral elements, polyhedra are automatically
 * tesselated with tetrahedra and pyramids (adding a vertex near each
 * polyhedron's center) unless discarded.
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque CGNS writer structure.
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)

void *
fvmc_to_cgns_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvmc_writer_time_dep_t   time_dependency,
                        MPI_Comm                comm);

#else

void *
fvmc_to_cgns_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvmc_writer_time_dep_t   time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to CGNS file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque CGNS writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvmc_to_cgns_finalize_writer(void  *this_writer_p);

/*----------------------------------------------------------------------------
 * Associate new time step with a CGNS geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvmc_to_cgns_set_mesh_time(void     *this_writer_p,
                          int       time_step,
                          double    time_value);

/*----------------------------------------------------------------------------
 * Indicate if elements of a given type in a mesh associated with a given
 * CGNS file writer need to be tesselated.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *   element_type  <-- element type we are interested in
 *
 * returns:
 *   1 if tesselation of the given element type is needed, 0 otherwise
 *----------------------------------------------------------------------------*/

int
fvmc_to_cgns_needs_tesselation(fvmc_writer_t       *this_writer_p,
                              const fvmc_nodal_t  *mesh,
                              fvmc_element_t       element_type);

/*----------------------------------------------------------------------------
 * Write nodal mesh to a CGNS file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvmc_to_cgns_export_nodal(void               *this_writer_p,
                         const fvmc_nodal_t  *mesh);

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a CGNS file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvmc_to_cgns_export_field(void                   *this_writer_p,
                         const fvmc_nodal_t      *mesh,
                         const char             *name,
                         fvmc_writer_var_loc_t    location,
                         int                     dimension,
                         fvmc_interlace_t         interlace,
                         int                     n_parent_lists,
                         const fvmc_lnum_t        parent_num_shift[],
                         fvmc_datatype_t          datatype,
                         int                     time_step,
                         double                  time_value,
                         const void       *const field_values[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* HAVE_CGNS */

#endif /* __FVMC_TO_CGNS_H__ */
