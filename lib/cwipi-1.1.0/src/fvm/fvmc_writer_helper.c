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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_error.h>
#include <bftc_file.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_convert_array.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_writer.h"
#include "fvmc_writer_helper.h"

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

#define FVMC_WRITER_MIN_ELEMENTS     32
#define FVMC_WRITER_MIN_SUB_ELEMENTS 32


/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM nodal to writer field output helper
 *----------------------------------------------------------------------------*/

struct _fvmc_writer_field_helper_t {

  int                         field_dim;       /* Field dimension */
  fvmc_interlace_t             interlace;       /* Field interlaced or not */
  fvmc_datatype_t              datatype;        /* Output datatype */
  fvmc_writer_var_loc_t        location;        /* Variable location */

  fvmc_gnum_t                  input_size;     /* Total input support size
                                                 (field values / dimension) */
  fvmc_gnum_t                  output_size;     /* Total output support size
                                                  (field values / dimension) */

  /* Local dimensions */

  fvmc_lnum_t  n_vertices;                  /* Local number of base vertices */
  fvmc_lnum_t  n_vertices_add;              /* Local number of added vertices */
  fvmc_lnum_t  max_vertices_add;            /* Maximum local number of added
                                              vertices per mesh section */

  fvmc_lnum_t  max_grouped_elements;        /* Maximum local number of elements
                                              per grouped sections */
  fvmc_lnum_t  max_grouped_elements_out;    /* Maximum local number of output
                                              elements per grouped sections */

  fvmc_lnum_t  max_section_elements;        /* Maximum local number of elements
                                              per section */
  fvmc_lnum_t  max_section_elements_out;    /* Maximum local number of
                                              output elements per section
                                              (= max_section_elements if
                                              no tesselation is used) */
  fvmc_lnum_t  n_sub_elements_max;          /* Maximum number of sub-elements
                                              per element (1 if no tesselation
                                              is used, global maximum in
                                              parallel mode) */

  int         n_added_vertex_sections;     /* Number of polyhedral sections
                                              adding vertices */
  int        *added_vertex_section;        /* List of polyhedral sections
                                              adding vertices */


  /* State related fields */

  fvmc_lnum_t                  start_id;      /* Local section start */
  const fvmc_writer_section_t *last_section;  /* Current section pointer */

  int         last_added_vertex_section;     /* Index of current added vertex
                                                section in list */

#if defined(FVMC_HAVE_MPI)

  /* Global dimensions */

  fvmc_gnum_t  n_vertices_g;                /* Global number of base vertices */
  fvmc_gnum_t  n_vertices_add_g;            /* Global number of added vertices */
  fvmc_gnum_t  max_vertices_add_g;          /* Maximum global number of added
                                              vertices per mesh section */

  fvmc_gnum_t  max_grouped_elements_g;      /* Maximum global number of elements
                                              per grouped sections */
  fvmc_gnum_t  max_grouped_elements_out_g;  /* Maximum global number of output
                                              elements per grouped sections */

  fvmc_gnum_t  max_section_elements_g;      /* Maximum global number of elements
                                              per section */
  fvmc_gnum_t  max_section_elements_out_g;  /* Maximum global number of
                                              output elements per section
                                              (= max_section_elements_g if
                                              no tesselation is used) */

  /* Arrays for parallel data preparation and distribution */

  size_t          local_buffer_size;        /* Max. values in local buffer */

  size_t          local_index_size;         /* Max. values in local index */
  size_t          global_index_size_g;      /* Max. values in global index */

  fvmc_lnum_t  *local_idx;                   /* Local tesselation index */
  fvmc_gnum_t  *global_idx;                  /* Global tesselation index */

  void *local_buffer;                       /* Local buffer (parallel mode) */

  /* Additionnal parallel state */

  MPI_Comm                    comm;         /* Associated MPI communicator */
  int                         rank;         /* Rank in communicator */
  fvmc_gather_slice_t         *slice;        /* Slice structure for
                                               partial gathers */

  fvmc_gnum_t                  start_num_g;  /* Global section start */

#endif
};

/*============================================================================
 * Static and constant variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare sections by associated type
 *
 * parameters:
 *   s1 <-- pointer to first section
 *   s2 <-- pointer to second section
 *
 * returns:
 *   difference of section types
 *----------------------------------------------------------------------------*/

static int
_compare_sections(const void *s1_p,
                  const void *s2_p)
{
  const fvmc_writer_section_t *s1 = (const fvmc_writer_section_t *) s1_p;
  const fvmc_writer_section_t *s2 = (const fvmc_writer_section_t *) s2_p;

  return ((int)(s1->type) - (int)(s2->type));
}

#if defined(FVMC_HAVE_MPI )

/*----------------------------------------------------------------------------
 * If necessary, reduce a local end index to adjust it to a given global
 * past the end global number.
 *
 * parameters:
 *   this_section          <-- fvmc_nodal section structure
 *   end_id                <-> end_id that may require reduction
 *   global_num_end        <-- past the end (maximum + 1) parent element
 *                             global number
 *----------------------------------------------------------------------------*/

static inline void
_limit_end_id_g(const fvmc_nodal_section_t  *this_section,
                fvmc_lnum_t                 *end_id,
                fvmc_gnum_t                  global_num_end)
{
  int last_id;

  const fvmc_gnum_t *global_element_num
    = fvmc_io_num_get_global_num(this_section->global_element_num);

  for (last_id = *end_id - 1;
       (   last_id > 0
        && global_element_num[last_id] >= global_num_end);
       last_id --);

  if (last_id > -1)
    *end_id = last_id + 1;
}

/*----------------------------------------------------------------------------
 * Partially distribute per element field values to an output buffer
 * in parallel mode.
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

static int
_field_helper_step_eg(fvmc_writer_field_helper_t   *helper,
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
                      size_t                      *output_size)
{
  fvmc_writer_field_helper_t *h = helper;

  int  retval = 0;

  fvmc_gnum_t  slice_output_size = 0;
  fvmc_gnum_t  end_num_g = 0;

  int  stride = 1;
  fvmc_lnum_t  num_shift = 0;

  const fvmc_nodal_section_t  *section = export_section->section;
  const fvmc_lnum_t  *parent_entity_num = section->parent_element_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE)
    stride = h->field_dim;

  if (n_parent_lists == 0)
    num_shift = export_section->num_shift;

  /* When at the beginning of a slice, setup state */

  if (h->start_num_g == 1) {

    /* prepare slice; when FVM sections are appended, we must change slice
       ranges for each section; this implies destroying and re-creating
       slices (otherwise we may simply re_initialize the slice) */

    if (h->slice != NULL && export_section != h->last_section) {
      h->slice = fvmc_gather_slice_destroy(h->slice);
      if (h->global_idx != NULL)
        BFTC_FREE(h->global_idx);
    }

    if (h->slice == NULL) {
      size_t output_buffer_base_size = output_buffer_size / stride;
      h->slice = fvmc_gather_slice_create(section->global_element_num,
                                         output_buffer_base_size,
                                         h->comm);
      if (export_section->type != section->type)
        BFTC_MALLOC(h->global_idx, output_buffer_base_size + 1, fvmc_gnum_t);
    }

    else /* if we are here, we have changed dimension index, not section */
      fvmc_gather_slice_reinitialize(h->slice);

    h->start_id = 0;

  }

  /* If slice end has not been attained yet, gather to buffer,
     update, and exit */

  if (fvmc_gather_slice_advance(h->slice,
                               &(h->start_num_g),
                               &end_num_g) == 0) {

    /* Standard section */
    /*------------------*/

    if (export_section->type == section->type) {

      /* Convert array for whole section on first pass
         (reuse for other slices) */

      if (h->start_num_g == 1)
        fvmc_convert_array(src_dim,
                          src_dim_shift,
                          stride,
                          num_shift,
                          section->n_elements + num_shift,
                          src_interlace,
                          datatype,
                          h->datatype,
                          n_parent_lists,
                          parent_num_shift,
                          parent_entity_num,
                          field_values,
                          h->local_buffer);

      /* Gather strided values for current slice */

      fvmc_gather_array(h->local_buffer,
                       output_buffer,
                       fvmc_datatype_to_mpi[h->datatype],
                       stride,
                       section->global_element_num,
                       h->comm,
                       h->slice);

      if (h->rank == 0)
        slice_output_size = end_num_g - h->start_num_g;

      h->start_num_g = end_num_g;

    }

    /* Tesselated section */
    /*--------------------*/

    else if (export_section->type != section->type) {

      fvmc_lnum_t  end_id;
      fvmc_gnum_t  local_buffer_base_size = h->local_buffer_size / stride;
      fvmc_gnum_t  global_s_size_loc = output_buffer_size;

      const fvmc_tesselation_t  *tesselation = section->tesselation;

      /* Build element->sub_elements index */

      end_id
        = fvmc_tesselation_range_index_g(tesselation,
                                        export_section->type,
                                        stride,
                                        h->start_id,
                                        local_buffer_base_size,
                                        &end_num_g,
                                        h->local_idx,
                                        h->comm);

      /* Check if the maximum id returned on some ranks leads to a
         lower end_num_g than initially required (due to the
         local buffer being too small) and adjust slice if necessary */

      fvmc_gather_slice_limit(h->slice, &end_num_g);

      /* Gather element->sub-elements index */

      fvmc_gather_slice_index(h->local_idx,
                             h->global_idx,
                             section->global_element_num,
                             h->comm,
                             h->slice);

      /* Recompute maximum value of end_num_g for this slice */

      fvmc_gather_resize_indexed_slice(FVMC_WRITER_MIN_ELEMENTS,
                                      &end_num_g,
                                      &global_s_size_loc,
                                      h->comm,
                                      h->global_idx,
                                      h->slice);

      _limit_end_id_g(section, &end_id, end_num_g);

      /* The buffer passed to this function should not be too small,
         as its size should have been computed in a consistent manner,
         but we check here in case it was not defined correctly by the
         calling writer. */

      if (global_s_size_loc > output_buffer_size) {
        bftc_error(__FILE__, __LINE__, 0,
                  _("Output buffer too small:\n"
                    "Current size = %lu, minimum size required = %lu."),
                  (unsigned long)output_buffer_size,
                  (unsigned long)(global_s_size_loc * stride));
      }

      /* Extract data from parent; note that the following call to
         fvmc_tesselation_distribute() will modify the var_buffer
         array (which avoids requiring an additional buffer) */

      fvmc_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        h->start_id + num_shift,
                        end_id + num_shift,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        h->local_buffer);

      /* Distribute data first then gather indexed values */

      fvmc_tesselation_distribute(section->tesselation,
                                 export_section->type,
                                 h->start_id,
                                 end_id,
                                 fvmc_datatype_size[h->datatype] * stride,
                                 h->local_buffer);

      /* Now gather distributed values */

      fvmc_gather_indexed(h->local_buffer,
                         output_buffer,
                         fvmc_datatype_to_mpi[h->datatype],
                         h->local_idx,
                         section->global_element_num,
                         h->comm,
                         h->global_idx,
                         h->slice);

      if (h->rank == 0)
        slice_output_size = h->global_idx[end_num_g - h->start_num_g] / stride;

      h->start_id = end_id;
      h->start_num_g = end_num_g;

    } /* End for tesselated section */

  } /* End for this slice */

  /* If we have reached the end of the section, reset indexes */

  else {

    h->last_section = export_section;
    h->start_num_g  = 1;
    h->start_id     = 0;

    retval = 1;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  return retval;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Partially distribute per element field values to an output buffer
 * in serial mode.
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

static int
_field_helper_step_el(fvmc_writer_field_helper_t   *helper,
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
                      size_t                      *output_size)
{
  fvmc_writer_field_helper_t *h = helper;

  int  retval = 0;

  fvmc_gnum_t  slice_output_size = 0;

  int  stride = 1;
  fvmc_lnum_t  end_id = 0;
  fvmc_lnum_t  num_shift = 0;

  size_t output_buffer_base_size = output_buffer_size;

  const fvmc_nodal_section_t  *section = export_section->section;
  const fvmc_lnum_t  *parent_entity_num = section->parent_element_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE) {
    stride = h->field_dim;
    output_buffer_base_size /= h->field_dim;
  }

  if (n_parent_lists == 0)
    num_shift = export_section->num_shift;

  /* If slice end has not been attained yet, extract values,
     update, and exit */

  if (h->start_id < section->n_elements) {

    /* Standard section */
    /*------------------*/

    if (export_section->type == section->type) {

      end_id = FVMC_MIN(h->start_id + (fvmc_lnum_t)output_buffer_base_size,
                       section->n_elements);

      fvmc_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        h->start_id + num_shift,
                        end_id + num_shift,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        output_buffer);

      slice_output_size = end_id - h->start_id;

    }

    /* Tesselated section */
    /*--------------------*/

    else if (export_section->type != section->type) {

      fvmc_lnum_t  n_sub_elements_max = 0;

      const fvmc_tesselation_t  *tesselation = section->tesselation;
      const fvmc_lnum_t *sub_element_idx
        = fvmc_tesselation_sub_elt_index(section->tesselation,
                                        export_section->type);

      fvmc_lnum_t  output_buffer_size_min
        = fvmc_tesselation_n_sub_elements(section->tesselation,
                                         export_section->type);

      fvmc_tesselation_get_global_size(section->tesselation,
                                      export_section->type,
                                      NULL,
                                      &n_sub_elements_max);

      output_buffer_size_min = FVMC_MIN(output_buffer_size_min,
                                       (  n_sub_elements_max
                                        * FVMC_WRITER_MIN_SUB_ELEMENTS));

      /* The buffer passed to this function should not be too small,
         as its size should have been computed in a consistent manner,
         but we check here in case it was not defined correctly by the
         calling writer. */

      if ((size_t)output_buffer_size_min > output_buffer_base_size) {
        bftc_error(__FILE__, __LINE__, 0,
                  _("Output buffer too small:\n"
                    "Current size = %lu, minimum size required = %lu."),
                  (unsigned long)output_buffer_size,
                  (unsigned long)(output_buffer_size_min * stride));
      }

      /* compute index end to fill output_buffer after distribution */

      for (end_id = h->start_id;
           (   end_id < section->n_elements
            &&  (sub_element_idx[end_id] <   (fvmc_lnum_t)output_buffer_base_size
                                           + sub_element_idx[h->start_id]));
           end_id++);
      if (  sub_element_idx[end_id] - sub_element_idx[h->start_id]
          > (fvmc_lnum_t)output_buffer_base_size)
        end_id--;

      /* Extract data from parent; note that the following call to
         fvmc_tesselation_distribute() will modify the var_buffer
         array (which avoids requiring an additional buffer) */

      fvmc_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        h->start_id + num_shift,
                        end_id + num_shift,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        output_buffer);

      /* distribute data to tesselation and write values */

      fvmc_tesselation_distribute(tesselation,
                                 export_section->type,
                                 h->start_id,
                                 end_id,
                                 fvmc_datatype_size[h->datatype] * stride,
                                 output_buffer);

      slice_output_size = (  sub_element_idx[end_id]
                           - sub_element_idx[h->start_id]);

    } /* End for tesselated section */

    h->start_id = end_id;

  } /* End for this slice */

  /* If we have reached the end of the section, reset indexes */

  else {

    h->last_section = export_section;
    h->start_id     = 0;

    retval = 1;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  return retval;
}

#if defined(FVMC_HAVE_MPI )

/*----------------------------------------------------------------------------
 * Partially distribute per node field values to an output buffer
 * in parallel mode.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   mesh               <-- pointer to nodal mesh
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

static int
_field_helper_step_ng(fvmc_writer_field_helper_t   *helper,
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
                      size_t                      *output_size)
{
  fvmc_writer_field_helper_t *h = helper;

  int  retval = 0;

  fvmc_gnum_t  slice_output_size = 0;
  fvmc_gnum_t  end_num_g = 0;

  int  stride = 1;

  const fvmc_lnum_t  *parent_entity_num = mesh->parent_vertex_num;
  const fvmc_io_num_t  *global_entity_num = mesh->global_vertex_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE)
    stride = h->field_dim;

  /* When at the beginning of a slice, setup state */

  if (h->start_num_g == 1) {

    /* prepare slice; when polyhedral FVM sections are appended (for
       extra vertices), we must change slice ranges for each section;
       this implies destroying and re-creating slices (otherwise we
       may simply re_initialize the slice) */

    if (h->slice != NULL && h->n_added_vertex_sections > 0) {
      h->slice = fvmc_gather_slice_destroy(h->slice);
      if (h->global_idx != NULL)
        BFTC_FREE(h->global_idx);
    }

    if (h->slice == NULL) {
      size_t output_buffer_base_size = output_buffer_size / stride;
      h->slice = fvmc_gather_slice_create(global_entity_num,
                                         output_buffer_base_size,
                                         h->comm);
    }

    else /* if we are here, we have changed dimension index, not section */
      fvmc_gather_slice_reinitialize(h->slice);

    h->start_id = 0;

  }

  /* Main vertices */
  /*---------------*/

  /* If slice end has not been attained yet, gather to buffer,
     update, and exit */

  if (   (h->start_num_g < h->n_vertices_g + 1)
      && fvmc_gather_slice_advance(h->slice,
                                  &(h->start_num_g),
                                  &end_num_g) == 0) {

    /* Convert array for whole section on first pass
       (reuse for other slices) */

    if (h->start_num_g == 1)
      fvmc_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        0,
                        mesh->n_vertices,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        h->local_buffer);

    /* Gather strided values for current slice */

    fvmc_gather_array(h->local_buffer,
                     output_buffer,
                     fvmc_datatype_to_mpi[h->datatype],
                     stride,
                     global_entity_num,
                     h->comm,
                     h->slice);

    if (h->rank == 0)
      slice_output_size = end_num_g - h->start_num_g;

    h->start_num_g = end_num_g;

  } /* End for this slice */

  /* Additional vertices in case of tesselation */
  /*--------------------------------------------*/

  else if (h->last_added_vertex_section < h->n_added_vertex_sections) {

    fvmc_gnum_t n_g_vertices_section;

    const fvmc_nodal_section_t  *section
      = mesh->sections[h->added_vertex_section[h->last_added_vertex_section]];

    /* Temporarily shift h->start_num_g to section value */

    h->start_num_g -= h->n_vertices_g ;

    n_g_vertices_section = fvmc_tesselation_n_g_vertices_add(section->tesselation);

    /* Interpolate values if both source and destination are floating point */

    if (   (h->datatype == FVMC_DOUBLE || h->datatype == FVMC_FLOAT)
        && (datatype == FVMC_DOUBLE || datatype == FVMC_FLOAT)) {

      /* Build a new slice if changing section */

      if (h->start_num_g == 1) {
        size_t output_buffer_base_size = output_buffer_size / stride;
        global_entity_num = section->global_element_num ;
        h->slice = fvmc_gather_slice_destroy(h->slice);
        h->slice = fvmc_gather_slice_create(global_entity_num,
                                           output_buffer_base_size,
                                           h->comm);
      }

      /* If slice end has not been attained yet, gather to buffer,
         update, and exit */

      if (   (h->start_num_g < h->n_vertices_g + 1)
           && fvmc_gather_slice_advance(h->slice,
                                       &(h->start_num_g),
                                       &end_num_g) == 0) {

        /* Convert array for whole section on first pass
           (reuse for other slices) */

        if (h->start_num_g == 1) {

          fvmc_lnum_t n_vertices_section
            = fvmc_tesselation_n_vertices_add(section->tesselation);

          fvmc_tesselation_vertex_values(section->tesselation,
                                        src_dim,
                                        src_dim_shift,
                                        stride,
                                        0,
                                        n_vertices_section,
                                        src_interlace,
                                        datatype,
                                        h->datatype,
                                        n_parent_lists,
                                        parent_num_shift,
                                        mesh->parent_vertex_num,
                                        field_values,
                                        h->local_buffer);

        }

        /* Gather strided values for current slice */

        fvmc_gather_array(h->local_buffer,
                         output_buffer,
                         fvmc_datatype_to_mpi[h->datatype],
                         stride,
                         global_entity_num,
                         h->comm,
                         h->slice);

      }

    }   /* End for floating point values */

    /* Set to zero if source or destination is not floating point */

    else {

      unsigned char *output_buffer_v = (unsigned char *) output_buffer;

      end_num_g = h->start_num_g + output_buffer_size / stride;
      if (end_num_g > n_g_vertices_section + 1)
        end_num_g = n_g_vertices_section + 1;

      slice_output_size = end_num_g - h->start_num_g;

      if (h->rank == 0) {

        size_t ii;
        size_t slice_output_size_c = slice_output_size * stride;

        slice_output_size_c *= fvmc_datatype_size[datatype];
        for (ii = 0; ii < slice_output_size_c; ii++)
          output_buffer_v[ii] = 0.;

      }

    }  /* End for non-floating point values */

    /* Advance to next section if necessary */

    if (h->rank == 0)
      slice_output_size = end_num_g - h->start_num_g;

    h->start_num_g = end_num_g;

    if (end_num_g > n_g_vertices_section) {
      h->last_added_vertex_section += 1;
      h->start_num_g = 1;
    }

    /* Shift h->start_num_g back to a value greater than h->n_vertices_g */

    h->start_num_g += h->n_vertices_g ;

  }

  /* If we have reached the end of this vertex set, reset indexes */

  else {

    h->start_num_g  = 1;
    h->start_id     = 0;

    h->last_added_vertex_section = 0;

    retval = 1;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  return retval;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Partially distribute per node field values to an output buffer
 * in serial mode.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   mesh               <-- pointer to nodal mesh
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

static int
_field_helper_step_nl(fvmc_writer_field_helper_t   *helper,
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
                      size_t                      *output_size)
{
  fvmc_writer_field_helper_t *h = helper;

  int  retval = 0;

  fvmc_gnum_t  slice_output_size = 0;
  fvmc_lnum_t  end_id = 0;

  int  stride = 1;

  const fvmc_lnum_t  *parent_entity_num = mesh->parent_vertex_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE)
    stride = h->field_dim;

  /* Main vertices */
  /*---------------*/

  if (h->start_id < mesh->n_vertices) {

    fvmc_lnum_t n_entities = mesh->n_vertices;

    end_id = h->start_id + (output_buffer_size / stride);
    end_id = FVMC_MIN(end_id, n_entities);

    fvmc_convert_array(src_dim,
                      src_dim_shift,
                      stride,
                      h->start_id,
                      end_id,
                      src_interlace,
                      datatype,
                      h->datatype,
                      n_parent_lists,
                      parent_num_shift,
                      parent_entity_num,
                      field_values,
                      output_buffer);

    slice_output_size = end_id - h->start_id;

    h->start_id = end_id;

  }

  /* Additional vertices in case of tesselation */
  /*--------------------------------------------*/

  else if (h->last_added_vertex_section < h->n_added_vertex_sections) {

    size_t ii;
    size_t slice_output_size_c;
    fvmc_lnum_t n_vertices_section;

    const fvmc_nodal_section_t  *section
      = mesh->sections[h->added_vertex_section[h->last_added_vertex_section]];

    /* Temporarily shift h->start_id to section value */

    h->start_id -= h->n_vertices ;

    n_vertices_section = fvmc_tesselation_n_vertices_add(section->tesselation);

    end_id = h->start_id + output_buffer_size / stride;
    if (end_id > n_vertices_section)
      end_id = n_vertices_section;

    slice_output_size = end_id - h->start_id;
    slice_output_size_c = slice_output_size * stride;

    if (   (h->datatype == FVMC_DOUBLE || h->datatype == FVMC_FLOAT)
        && (datatype == FVMC_DOUBLE || datatype == FVMC_FLOAT)) {

      fvmc_tesselation_vertex_values(section->tesselation,
                                    src_dim,
                                    src_dim_shift,
                                    stride,
                                    h->start_id,
                                    end_id,
                                    src_interlace,
                                    datatype,
                                    h->datatype,
                                    n_parent_lists,
                                    parent_num_shift,
                                    mesh->parent_vertex_num,
                                    field_values,
                                    output_buffer);

    }
    else {
      unsigned char *output_buffer_v = (unsigned char *) output_buffer;
      slice_output_size_c *= fvmc_datatype_size[datatype];
      for (ii = 0; ii < slice_output_size_c; ii++)
        output_buffer_v[ii] = 0.;
    }

    /* Advance to next section if necessary */

    h->start_id = end_id;

    if (end_id >= n_vertices_section) {
      h->last_added_vertex_section += 1;
      h->start_id = 0;
    }

    /* Shift h->start_id back to a value greater than h->n_vertices */

    h->start_id += h->n_vertices ;

  }

  /* If we have reached the end of this vertex set, reset indexes */

  else {

    h->start_id     = 0;

    h->last_added_vertex_section = 0;

    retval = 1;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  return retval;
}

/*============================================================================
 * Semi-private function definitions (prototypes in fvmc_writer_helper.h)
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
                       _Bool               divide_polyhedra)
{
  int  i, j;
  int  n_sections = 0;
  int  n_sub_types = 0;
  fvmc_element_t sub_type[FVMC_TESSELATION_N_SUB_TYPES_MAX];
  fvmc_writer_section_t *export_list = NULL;

  fvmc_lnum_t  num_shift = 0;
  fvmc_gnum_t  extra_vertex_base = fvmc_nodal_n_g_vertices(mesh) + 1;
  const int  export_dim = fvmc_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  n_sections = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvmc_nodal_section_t  *const  section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (section->entity_dim == export_dim) {

      /* Optionally discard polygons or polyhedra */
      if (   (section->type == FVMC_FACE_POLY && discard_polygons == true)
          || (section->type == FVMC_CELL_POLY && discard_polyhedra == true))
        continue;

      /* Possibly add multiple sections when dividing polygons or polyhedra
         (ignore those sections if no tesselation present) */
      else if (   (section->type == FVMC_FACE_POLY && divide_polygons == true)
               || (section->type == FVMC_CELL_POLY && divide_polyhedra == true)) {
        if (section->tesselation != NULL)
          n_sections += fvmc_tesselation_n_sub_types(section->tesselation);
      }

      else
        n_sections += 1;

    }

  }

  /* If no sections are present no list is returned */

  if (n_sections == 0)
    return NULL;

  BFTC_MALLOC(export_list, n_sections, fvmc_writer_section_t);

  for (i = 0 ; i < n_sections - 1 ; i++)
    (export_list[i]).next = export_list + i + 1;
  (export_list[n_sections - 1]).next = NULL;

  /* Build unsorted list */

  n_sections = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvmc_tesselation_t    *tesselation = NULL;
    const fvmc_nodal_section_t  *const  section = mesh->sections[i];

    /* Ignore sections with entity dimension other than the highest */
    if (section->entity_dim != export_dim)
      continue;

    /* Ignore polygonal or polyhedra sections if they should be discarded */
    if (   (section->type == FVMC_FACE_POLY && discard_polygons == true)
        || (section->type == FVMC_CELL_POLY && discard_polyhedra == true))
      continue;

    /* Depending on section type and tesselation,
       we may have multiple sub-types */

    n_sub_types = 1;
    sub_type[0] = section->type;

    if (   (section->type == FVMC_FACE_POLY && divide_polygons == true)
        || (section->type == FVMC_CELL_POLY && divide_polyhedra == true)) {

      if (section->tesselation != NULL) {
        tesselation = section->tesselation;
        n_sub_types = fvmc_tesselation_n_sub_types(section->tesselation);
        for (j = 0; j < n_sub_types; j++)
          sub_type[j] = fvmc_tesselation_sub_type(section->tesselation, j);
      }
      else
        continue; /* ignore section if no tesselation present */

    }

    for (j = 0; j < n_sub_types; j++) {

      (export_list[n_sections]).section = section;
      (export_list[n_sections]).type = sub_type[j];
      (export_list[n_sections]).continues_previous = false;
      if (tesselation == NULL)
        (export_list[n_sections]).extra_vertex_base = 0;
      else
        (export_list[n_sections]).extra_vertex_base = extra_vertex_base;
      (export_list[n_sections]).num_shift = num_shift;

      n_sections++;

    }

    if (tesselation != NULL)
      extra_vertex_base += fvmc_tesselation_n_g_vertices_add(tesselation);

    num_shift += section->n_elements;

  }

  /* Now order list so as to fuse sections of similar type */

  if (group_same_type == true && n_sections > 1) {

    qsort(export_list, n_sections, sizeof(fvmc_writer_section_t),
          _compare_sections);

    for (i = 1; i < n_sections; i++) {
      if ((export_list[i-1]).type == (export_list[i]).type)
        (export_list[i]).continues_previous = true;
    }

  }

  for (i = 0; i < n_sections - 1; i++)
    (export_list[i]).next = &(export_list[i+1]);
  export_list[n_sections - 1].next = NULL;

  return export_list;
}

/*----------------------------------------------------------------------------
 * Create field writer helper structure.
 *
 * Local values are initialized, ang lobal values are set to zero
 * (they may be initialized by calling fvmc_writer_field_helper_init_g()).
 *
 * The mesh argument is not used when location is FVMC_WRITER_PER_ELEMENT,
 * so NULL may be given instead in this case. The section_list
 * argument is used even when location is FVMC_WRITER_PER_NODE to determine
 * if polyhedra are are divided (using extra vertices). Thus, NULL
 * may be given instead only if we are sure the corresponding exported mesh
 * does not contain divided polyhedra.
 *
 * parameters:
 *   mesh            <-- pointer to nodal mesh structure
 *   section_list    <-- point to export section list helper structure
 *   field_dim       <-- indicates output field dimension
 *   interlace       <-- indicates if output is interlaced
 *   location        <-- indicates if field is cell or node based
 *   datatype        <-- datatype of destination buffers
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
                               fvmc_writer_var_loc_t        location)
{
  fvmc_writer_field_helper_t *h = NULL;

  /* Initialize structure */

  BFTC_MALLOC(h, 1, fvmc_writer_field_helper_t);

  /* Initialize structure */

  h->field_dim = field_dim;
  h->interlace = interlace;
  h->datatype = datatype;
  h->location = location;

  /* Zero fields that will be computed later */

  h->input_size = 0;
  h->output_size = 0;

  h->n_vertices = 0;
  h->n_vertices_add = 0;
  h->max_vertices_add = 0;

  h->max_grouped_elements = 0;
  h->max_grouped_elements_out = 0;
  h->max_section_elements = 0;
  h->max_section_elements_out = 0;

  h->n_sub_elements_max = 1;

  h->n_added_vertex_sections = 0;
  h->added_vertex_section = NULL;

  /* State values */

  h->start_id = 0;
  h->last_section = NULL;

  h->last_added_vertex_section = 0;

#if defined(FVMC_HAVE_MPI)

  h->n_vertices_g = 0;
  h->n_vertices_add_g = 0;
  h->max_vertices_add_g = 0;
  h->max_grouped_elements_g = 0;
  h->max_grouped_elements_out_g = 0;

  h->max_section_elements_g = 0;
  h->max_section_elements_out_g = 0;

  h->local_buffer_size= 0;
  h->local_index_size = 0;
  h->global_index_size_g = 0;

  h->local_buffer = NULL;
  h->local_idx = NULL;
  h->global_idx = NULL;

  h->comm = MPI_COMM_NULL;
  h->rank = -1;
  h->slice = NULL;

  h->start_num_g = 1;

#endif /* defined(FVMC_HAVE_MPI) */

  /* Compute local dimensions */

  if (location == FVMC_WRITER_PER_ELEMENT) {

    const fvmc_writer_section_t *export_section = section_list;

    fvmc_lnum_t n_grouped_elements = 0;
    fvmc_lnum_t n_grouped_elements_out = 0;

    while (export_section != NULL) {

      const fvmc_nodal_section_t *section = export_section->section;

      fvmc_lnum_t n_elements;
      fvmc_lnum_t n_sub_elements;
      fvmc_lnum_t n_sub_elements_max = 1;

      if (export_section->continues_previous == false) {
        n_grouped_elements = 0;
        n_grouped_elements_out = 0;
      }

      n_elements = section->n_elements;

      if (export_section->type == section->type) { /* Regular section */
        n_sub_elements = n_elements;
      }
      else { /* Tesselated section */
        fvmc_tesselation_get_global_size(section->tesselation,
                                        export_section->type,
                                        NULL,
                                        &n_sub_elements_max);
        n_sub_elements = fvmc_tesselation_n_sub_elements(section->tesselation,
                                                        export_section->type);
      }

      n_grouped_elements += section->n_elements;
      n_grouped_elements_out += n_sub_elements;

      /* Update global helper values */

      h->input_size  += section->n_elements;
      h->output_size += n_sub_elements;

      h->max_grouped_elements     = FVMC_MAX(h->max_grouped_elements,
                                            n_grouped_elements);
      h->max_grouped_elements_out = FVMC_MAX(h->max_grouped_elements_out,
                                            n_grouped_elements_out);

      h->max_section_elements     = FVMC_MAX(h->max_section_elements,
                                            section->n_elements);
      h->max_section_elements_out = FVMC_MAX(h->max_section_elements_out,
                                            n_sub_elements);

      h->n_sub_elements_max = FVMC_MAX(h->n_sub_elements_max,
                                      n_sub_elements_max);

      /* continue with next section */

      export_section = export_section->next;

    } /* End of loop on sections */

  }

  else if (location == FVMC_WRITER_PER_NODE) {

    int n_added_vertex_sections_max = 0;

    const fvmc_writer_section_t *export_section;

    h->n_vertices   = mesh->n_vertices;

    h->input_size  = h->n_vertices;
    h->output_size = h->n_vertices;

    /* Determine if polyhedra are tesselated */

    for (export_section = section_list;
         export_section != NULL;
         export_section = export_section->next) {

      const fvmc_nodal_section_t *section = export_section->section;

      if (   export_section->type != section->type
          && section->type == FVMC_CELL_POLY)
        n_added_vertex_sections_max += 1;

    }

    /* In case of polyhedra tesselation, count for added vertices;
       at this stage, we know that polyhedra are tesselated */

    if (n_added_vertex_sections_max > 0) {

      int ii, jj;

      for (ii = 0, jj = 0; ii < mesh->n_sections; ii++) {
        if ((mesh->sections[ii])->type == FVMC_CELL_POLY)
          h->n_added_vertex_sections += 1;
      }

      BFTC_MALLOC(h->added_vertex_section, h->n_added_vertex_sections, int);

      for (ii = 0, jj = 0; ii < mesh->n_sections; ii++) {

        const fvmc_nodal_section_t *section = mesh->sections[ii];

        if (section->type == FVMC_CELL_POLY) {

          fvmc_lnum_t n_vertices_add
            = fvmc_tesselation_n_vertices_add(section->tesselation);

          h->output_size += n_vertices_add;

          h->n_vertices_add += n_vertices_add;
          h->max_vertices_add = FVMC_MAX(h->max_vertices_add,
                                        n_vertices_add);

          h->added_vertex_section[jj] = ii;
          jj++;

        }

      }
    }

  }

  return h;
}

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
fvmc_writer_field_helper_destroy(fvmc_writer_field_helper_t *helper)
{
  fvmc_writer_field_helper_t *h = helper;

  if (h->added_vertex_section != NULL)
    BFTC_FREE(h->added_vertex_section);

#if defined(FVMC_HAVE_MPI)

  if (h->slice != NULL)
    h->slice = fvmc_gather_slice_destroy(h->slice);

  BFTC_FREE(h->global_idx);
  BFTC_FREE(h->local_idx);
  BFTC_FREE(h->local_buffer);

#endif

  BFTC_FREE(h);

  return h;
}

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
                               MPI_Comm                     comm)
{
  int comm_size = 1;
  fvmc_writer_field_helper_t *h = helper;

  h->input_size = 0; /* Previously computed, so reset */
  h->output_size = 0;

  /* Compute global dimensions */

  if (h->location == FVMC_WRITER_PER_ELEMENT) {

    const fvmc_writer_section_t *export_section = section_list;

    fvmc_gnum_t n_grouped_elements_g = 0;
    fvmc_gnum_t n_grouped_elements_out_g = 0;

    while (export_section != NULL) {

      const fvmc_nodal_section_t *section = export_section->section;

      fvmc_gnum_t n_g_elements;
      fvmc_gnum_t n_g_sub_elements;
      fvmc_lnum_t n_sub_elements_max = 1;

      size_t local_buffer_size = section->n_elements;
      size_t local_index_size = 0;
      size_t global_index_size_g = 0;

      if (export_section->continues_previous == false) {
        n_grouped_elements_g = 0;
        n_grouped_elements_out_g = 0;
      }

      if (section->global_element_num != NULL) {
        n_g_elements
          = fvmc_io_num_get_global_count(section->global_element_num);
      }
      else {
        n_g_elements = section->n_elements;
      }

      if (export_section->type == section->type) { /* Regular section */
        n_g_sub_elements = n_g_elements;
      }
      else { /* Tesselated section */
        fvmc_tesselation_get_global_size(section->tesselation,
                                        export_section->type,
                                        &n_g_sub_elements,
                                        &n_sub_elements_max);
      }

      n_grouped_elements_g += n_g_elements;
      n_grouped_elements_out_g += n_g_sub_elements;

      /* Update global dimensions */

      h->input_size += n_g_elements;
      h->output_size += n_g_sub_elements;

      h->max_grouped_elements_g     = FVMC_MAX(h->max_grouped_elements_g,
                                              n_grouped_elements_g);
      h->max_grouped_elements_out_g = FVMC_MAX(h->max_grouped_elements_out_g,
                                              n_grouped_elements_out_g);

      h->max_section_elements_g     = FVMC_MAX(h->max_section_elements_g,
                                              n_g_elements);
      h->max_section_elements_out_g = FVMC_MAX(h->max_section_elements_out_g,
                                              n_g_sub_elements);

      /* Update buffer sizes */

      h->local_buffer_size = FVMC_MAX(h->local_buffer_size,
                                     (size_t)section->n_elements);

      if (export_section->type != section->type) { /* Tesselated section */
        local_buffer_size = FVMC_MAX(section->n_elements,
                                    (  n_sub_elements_max
                                     * FVMC_WRITER_MIN_SUB_ELEMENTS));
        local_index_size = section->n_elements + 1;
        global_index_size_g = n_g_elements + 1;
      }
      h->local_buffer_size = FVMC_MAX(h->local_buffer_size,
                                     local_buffer_size);
      h->local_index_size = FVMC_MAX(h->local_index_size,
                                    local_index_size);
      h->global_index_size_g = FVMC_MAX(h->global_index_size_g,
                                       global_index_size_g);

      /* continue with next section */

      export_section = export_section->next;

    } /* End of loop on sections */

  }

  else if (h->location == FVMC_WRITER_PER_NODE) {

    h->n_vertices_g = fvmc_nodal_n_g_vertices(mesh);

    h->input_size = h->n_vertices_g;
    h->output_size = h->n_vertices_g;

    h->local_buffer_size = mesh->n_vertices;

    /* If polyhedra are tesselated */

    if (h->n_added_vertex_sections > 0) {

      int ii;

      for (ii = 0; ii < h->n_added_vertex_sections; ii++) {

        fvmc_gnum_t n_g_vertices_add;

        const fvmc_nodal_section_t  *section
          = mesh->sections[h->added_vertex_section
                           [h->last_added_vertex_section]];

        n_g_vertices_add
          = fvmc_tesselation_n_g_vertices_add(section->tesselation);

        h->output_size += n_g_vertices_add;

        h->n_vertices_add_g += n_g_vertices_add;
        h->max_vertices_add_g = FVMC_MAX(h->max_vertices_add_g,
                                        n_g_vertices_add);

        /* Update buffer sizes */

        h->local_buffer_size = FVMC_MAX(h->local_buffer_size,
                                       (size_t)section->n_elements);

      }
    }

  }

  /* With interlaced multidimensional arrays,
     multiply local buffer size by the stride */

  if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE) {
    h->local_buffer_size *= h->field_dim;
  }

  /* Get info on the current MPI communicator */

  if (comm != MPI_COMM_NULL) {
    h->comm = comm;
    MPI_Comm_rank(comm, &(h->rank));
    MPI_Comm_size(comm, &comm_size);
  }

  if (comm_size < 2)
    h->rank = -1;

  /* Allocate private arrays if this has not been done yet;
     Note that the global index buffer for tesselations can only
     be allocated once the output buffer's size is known, so it
     is deferred to the helper stepping function */

  if (h->local_buffer != NULL) {
    BFTC_FREE(h->global_idx);
    BFTC_FREE(h->local_idx);
    BFTC_FREE(h->local_buffer);
  }

  BFTC_MALLOC(h->local_buffer,
             h->local_buffer_size*fvmc_datatype_size[h->datatype],
             char);

  if (h->n_sub_elements_max > 1)
    BFTC_MALLOC(h->local_idx, h->local_index_size, fvmc_lnum_t);

}

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
                                 size_t  *min_output_buffer_size)
{
  const fvmc_writer_field_helper_t *h = helper;

  assert(h != NULL);

  if (input_size != NULL)
    *input_size = h->input_size;
  if (output_size != NULL)
    *output_size = h->output_size;

  if (max_grouped_elements_out != NULL)
    *max_grouped_elements_out = h->max_grouped_elements_out;

#if defined(FVMC_HAVE_MPI)

  if (h->rank > -1) {

    if (max_grouped_elements_out != NULL)
      *max_grouped_elements_out = h->max_grouped_elements_out_g;

  }

#endif

  if (min_output_buffer_size != NULL) {

    size_t min_size = 0;

    if (h->n_sub_elements_max > 1) {
      min_size = h->n_sub_elements_max * FVMC_WRITER_MIN_SUB_ELEMENTS;
      min_size = FVMC_MIN(min_size, h->output_size);
    }

    if (h->output_size > 0)
      min_size = FVMC_MAX(FVMC_WRITER_MIN_ELEMENTS, min_size);

    if (min_size > h->output_size)
      min_size = h->output_size;

    /* With interlaced multidimensional arrays,
       multiply buffer sizes by the stride */

    if (h->field_dim > 1 && h->interlace == FVMC_INTERLACE)
      min_size *= h->field_dim;

    *min_output_buffer_size = min_size;

  }

}

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
fvmc_writer_field_helper_field_dim(const fvmc_writer_field_helper_t  *helper)
{
  return helper->field_dim;
}

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
fvmc_writer_field_helper_datatype(const fvmc_writer_field_helper_t  *helper)
{
  return helper->datatype;
}

/*----------------------------------------------------------------------------
 * Partially distribute per element field values to an output buffer.
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
                               size_t                      *output_size)
{
  assert(helper != NULL);

#if defined(FVMC_HAVE_MPI)

  if (helper->rank > -1) {

    return _field_helper_step_eg(helper,
                                 export_section,
                                 src_dim,
                                 src_dim_shift,
                                 src_interlace,
                                 n_parent_lists,
                                 parent_num_shift,
                                 datatype,
                                 field_values,
                                 output_buffer,
                                 output_buffer_size,
                                 output_size);

  }

#endif /* defined(FVMC_HAVE_MPI) */

  return _field_helper_step_el(helper,
                               export_section,
                               src_dim,
                               src_dim_shift,
                               src_interlace,
                               n_parent_lists,
                               parent_num_shift,
                               datatype,
                               field_values,
                               output_buffer,
                               output_buffer_size,
                               output_size);
}

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
                               size_t                      *output_size)
{
  assert(helper != NULL);

#if defined(FVMC_HAVE_MPI)

  if (helper->rank > -1) {

    return _field_helper_step_ng(helper,
                                 mesh,
                                 src_dim,
                                 src_dim_shift,
                                 src_interlace,
                                 n_parent_lists,
                                 parent_num_shift,
                                 datatype,
                                 field_values,
                                 output_buffer,
                                 output_buffer_size,
                                 output_size);

  }

#endif /* defined(FVMC_HAVE_MPI) */

  return _field_helper_step_nl(helper,
                               mesh,
                               src_dim,
                               src_dim_shift,
                               src_interlace,
                               n_parent_lists,
                               parent_num_shift,
                               datatype,
                               field_values,
                               output_buffer,
                               output_buffer_size,
                               output_size);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
