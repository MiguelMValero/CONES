#ifndef __FVMC_WRITER_PRIV_H__
#define __FVMC_WRITER_PRIV_H__

/*============================================================================
 * Private types for mesh and field writers
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

/*----------------------------------------------------------------------------*/

#include "fvmc_config.h"

#if defined(FVMC_HAVE_MPI )
#include <mpi.h>
#endif

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

/*
 * Writer format implementation and functionality info
 */

#define FVMC_WRITER_FORMAT_USE_EXTERNAL    (1 << 0)

#define FVMC_WRITER_FORMAT_HAS_POLYGON     (1 << 1)
#define FVMC_WRITER_FORMAT_HAS_POLYHEDRON  (1 << 2)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int
(fvmc_writer_n_version_strings_t) (int format_index);

typedef const char *
(fvmc_writer_version_string_t)(int format_index,
                              int string_index,
                              int compile_time_version);

#if defined(FVMC_HAVE_MPI)

typedef void *
(fvmc_writer_init_t) (const char             *name,
                     const char             *path,
                     const char             *options,
                     fvmc_writer_time_dep_t   time_dependency,
                     MPI_Comm                comm);

#else

typedef void *
(fvmc_writer_init_t) (const char             *name,
                     const char             *path,
                     const char             *options,
                     fvmc_writer_time_dep_t   time_dependency);

#endif /* defined(FVMC_HAVE_MPI) */

typedef void *
(fvmc_writer_finalize_t) (void  *this_writer);

typedef void
(fvmc_writer_set_mesh_time_t) (void    *this_writer,
                              int      time_step,
                              double   time_value);

typedef int
(fvmc_writer_needs_tesselation_t) (fvmc_writer_t       *this_writer,
                                  const fvmc_nodal_t  *mesh,
                                  fvmc_element_t       element_type);

typedef void
(fvmc_writer_export_nodal_t) (void               *this_writer,
                             const fvmc_nodal_t  *mesh);

typedef void
(fvmc_writer_export_field_t) (void                   *this_writer,
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

typedef void
(fvmc_writer_flush_t) (fvmc_writer_t  *this_writer);

/*----------------------------------------------------------------------------
 * Format information structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                     name[32];     /* Format name */
  char                     version[16];  /* Format version (if defined) */
  int                      info_mask;    /* Additional format info */
  fvmc_writer_time_dep_t    max_time_dep; /* Maximum time dependency level
                                            possible with this format */

  fvmc_writer_n_version_strings_t  *n_version_strings_func;
  fvmc_writer_version_string_t     *version_string_func;
  fvmc_writer_init_t               *init_func;
  fvmc_writer_finalize_t           *finalize_func;
  fvmc_writer_set_mesh_time_t      *set_mesh_time_func;
  fvmc_writer_needs_tesselation_t  *needs_tesselation_func;
  fvmc_writer_export_nodal_t       *export_nodal_func;
  fvmc_writer_export_field_t       *export_field_func;
  fvmc_writer_flush_t              *flush_func;

} fvmc_writer_format_t;

/*----------------------------------------------------------------------------
 * Structure defining a writer definition
 *----------------------------------------------------------------------------*/

struct _fvmc_writer_t {

  char                   *name;           /* Writer name */
  fvmc_writer_format_t    *format;         /* Output format */
  fvmc_writer_time_dep_t   time_dep;       /* Geometry time dependency */
  void                   *format_writer;  /* Format-specific writer */

  double                  mesh_wtime;     /* Meshes output Wall-clock time */
  double                  mesh_cpu_time;  /* Meshes output CPU time */
  double                  field_wtime;    /* Fields output Wall-clock time */
  double                  field_cpu_time; /* Fields output CPU time */

};

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute recommended buffer sizes to input or output a nodal mesh
 * definition by slices. This is especially useful when gathering the mesh for
 * output by slices using standard I/O in parallel mode.
 *
 * The global number of vertices and elements of each slice may also
 * be returned, if the pointers n_g_vertices and n_g_elements_section
 * are non-NULL respectively.
 *
 * The number of slices indicated is a minimum, and only a target;
 * computation is based primarily on cell and face connectivity, and the
 * target should be met for strided connectivities on those types of elements
 * only. Using an "optimistic" (i.e. small) mean number of vertices per
 * polyhedra or polygon will typically lead to requiring more slices, as
 * the connectivity slice size returned will be smaller than that truly
 * required for the corresponding slice size.
 * Slice sizes required for edges connectivity will meet the target only
 * when the global numbers of cells and faces given are zero, so as to
 * avoid generating too large connectivity slice sizes for cells should a mesh
 * contain both (as for example a hexahedral connectivity slice is 8 times
 * larger than the corresponding slice size, while an edges connectivity is
 * only 2 times as large).
 *
 * parameters:
 *   this_nodal                 <-- pointer to nodal mesh structure
 *   n_slices                   <-- target number of slices required
 *   n_polyhedron_vertices_mean <-- estimate of the mean number of vertices
 *                                  per polyhedron
 *   n_polygon_vertices_mean    <-- estimate of the mean number of vertices
 *                                  per polygon
 *   n_g_vertices               --> global number of vertices (or NULL)
 *   n_g_elements_section       --> array for global number of elements per
 *                                  section (or NULL)
 *   global_s_size              --> maximum number of entities defined per slice
 *   global_connect_s_size      --> maximum number of connectivity values
 *                                  per slice
 *----------------------------------------------------------------------------*/

void
fvmc_writer_def_nodal_buf_size(const fvmc_nodal_t  *this_nodal,
                              int                 n_slices,
                              int                 n_polyhedron_vertices_mean,
                              int                 n_polygon_vertices_mean,
                              fvmc_gnum_t         *n_g_vertices,
                              fvmc_gnum_t          n_g_elements_section[],
                              fvmc_gnum_t         *global_s_size,
                              fvmc_gnum_t         *global_connect_s_size);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_WRITER_PRIV_H__ */
